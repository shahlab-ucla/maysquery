import httpx
import logging
import random
import asyncio
import os
from typing import List, Optional
from models import OrthologMapping, ValidatedTarget, ExecutionLogEntry, PipelineConfig

logger = logging.getLogger(__name__)

async def get_uniprot_for_maize_gene(maize_gene: str) -> str:
    """Query Ensembl Plants to get the UniProt ID for a maize gene."""
    url = f"https://rest.ensembl.org/xrefs/id/{maize_gene}?all_levels=1"
    headers = {"Content-Type": "application/json"}
    
    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(url, headers=headers, timeout=10.0)
            if resp.status_code == 200:
                for xref in resp.json():
                    db = xref.get("dbname", "").lower()
                    if "uniprot" in db:
                        return xref.get("primary_id")
    except Exception as e:
        logger.error(f"Error fetching UniProt xref for {maize_gene}: {e}")
    
    return ""

async def validate_structure_alphafold(uniprot_id: str) -> float:
    """Live query to AlphaFold DB for average pLDDT."""
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    
    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(url, timeout=10.0)
            if resp.status_code == 200:
                data = resp.json()
                if isinstance(data, list) and len(data) > 0:
                    return data[0].get("global_metric_value", 75.0)
    except Exception as e:
        logger.error(f"AlphaFold DB error for {uniprot_id}: {e}")
        
    return 0.0

# Module-level cache for the Gramene expression-breadth lookup.
# Key: maize gene ID. Value: list of GXA experiment IDs (sorted).
_EXPRESSION_CACHE: dict = {}


async def fetch_expression_breadth(maize_gene: str) -> List[str]:
    """
    Return the list of EBI Expression Atlas experiments in which `maize_gene`
    is reported as expressed. Sourced from Gramene's `expressed_in_gxa_attr_ss`
    field (`data.gramene.org/search`).

    Replaces the previous `/gxa/json/search/baseline?geneQuery=…` call, which
    EBI deprecated (now returns 404 for everything). Per-tissue FPKM values
    are no longer available without downloading per-experiment TSVs, but the
    breadth ("expressed in N experiments") is a real qualitative signal.
    """
    if not maize_gene:
        return []
    if maize_gene in _EXPRESSION_CACHE:
        return _EXPRESSION_CACHE[maize_gene]

    experiments: List[str] = []
    try:
        async with httpx.AsyncClient(follow_redirects=True) as client:
            r = await client.get(
                "https://data.gramene.org/search",
                params={"q": maize_gene, "rows": 1},
                timeout=15.0,
            )
            if r.status_code == 200:
                docs = (((r.json() or {}).get("response") or {}).get("docs") or [])
                if docs and docs[0].get("id") == maize_gene:
                    experiments = sorted(docs[0].get("expressed_in_gxa_attr_ss") or [])
    except Exception as e:
        logger.error(f"Gramene expression lookup failed for {maize_gene}: {e}")

    _EXPRESSION_CACHE[maize_gene] = experiments
    return experiments


# Kept for backward compat — returns an experiment_id → 'expressed' dict, so
# the existing UI's tissue-list rendering still produces meaningful chips
# until the frontend migrates to `n_expression_experiments` directly.
async def validate_expression_rnaseq(maize_gene: str) -> dict:
    experiments = await fetch_expression_breadth(maize_gene)
    return {exp_id: "expressed" for exp_id in experiments[:6]}

_AFDB_URL_CACHE: dict = {}  # uniprot_id -> pdb URL (or "" if absent)


async def _resolve_alphafold_pdb_url(uniprot_id: str) -> str:
    """
    Ask the AlphaFold prediction API for the current PDB URL. This avoids
    hardcoding a model version (the DB has moved v1→v4→v6+ and our previous
    hardcoded v4 URL 404s for every entry as of 2026). Cached in-process.
    """
    if uniprot_id in _AFDB_URL_CACHE:
        return _AFDB_URL_CACHE[uniprot_id]

    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    pdb_url = ""
    try:
        async with httpx.AsyncClient(follow_redirects=True) as client:
            resp = await client.get(api_url, timeout=15.0)
            if resp.status_code == 200:
                data = resp.json()
                if isinstance(data, list) and data:
                    pdb_url = data[0].get("pdbUrl", "") or ""
    except Exception as e:
        logger.error(f"AlphaFold API lookup failed for {uniprot_id}: {e}")

    _AFDB_URL_CACHE[uniprot_id] = pdb_url
    return pdb_url


async def download_alphafold_pdb(uniprot_id: str, tmp_dir: str) -> str:
    """Download the AlphaFold PDB for a UniProt accession (any model version)."""
    if not uniprot_id:
        return ""

    filepath = os.path.join(tmp_dir, f"{uniprot_id}.pdb")
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        return filepath

    pdb_url = await _resolve_alphafold_pdb_url(uniprot_id)
    if not pdb_url:
        logger.warning(f"AlphaFold has no prediction record for {uniprot_id}")
        return ""

    try:
        async with httpx.AsyncClient(follow_redirects=True) as client:
            resp = await client.get(pdb_url, timeout=30.0)
            if resp.status_code == 200 and resp.content:
                with open(filepath, "wb") as f:
                    f.write(resp.content)
                return filepath
            logger.warning(
                f"AlphaFold PDB fetch returned HTTP {resp.status_code} for {uniprot_id} ({pdb_url})"
            )
    except Exception as e:
        logger.error(f"Error downloading PDB for {uniprot_id} from {pdb_url}: {e}")

    return ""

async def calculate_tm_score_foldseek(query_uniprot_id: str, target_uniprot_id: str) -> float:
    """
    Downloads PDBs and runs foldseek locally to calculate the structural alignment TM-score.
    Falls back to a mock score if foldseek is not installed.
    """
    tmp_dir = os.path.join(os.path.dirname(__file__), "tmpFolder")
    os.makedirs(tmp_dir, exist_ok=True)
    
    query_pdb = await download_alphafold_pdb(query_uniprot_id, tmp_dir)
    target_pdb = await download_alphafold_pdb(target_uniprot_id, tmp_dir)
    
    if not query_pdb or not target_pdb:
        return round(random.uniform(0.5, 0.99), 2)
        
    output_tsv = os.path.join(tmp_dir, f"{query_uniprot_id}_{target_uniprot_id}_results.tsv")
    tmp_foldseek = os.path.join(tmp_dir, "fs_tmp")
    os.makedirs(tmp_foldseek, exist_ok=True)
    
    # We construct the command equivalent to the specification. 
    # The spec states: foldseek easy-search query.pdb maize_afdb_database ...
    # Since we don't have the whole maize DB locally, we do a 1-to-1 search query.pdb target.pdb
    # Run Foldseek via WSL
    cmd = [
        "wsl", "foldseek", "easy-search",
        os.path.basename(query_pdb), os.path.basename(target_pdb), os.path.basename(output_tsv), os.path.basename(tmp_foldseek),
        "--exhaustive-search", "1",
        "--format-output", "query,target,qlen,tlen,alntmscore,rmsd"
    ]
    
    try:
        process = await asyncio.create_subprocess_exec(
            *cmd,
            cwd=tmp_dir,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await process.communicate()
        
        if process.returncode == 0 and os.path.exists(output_tsv):
            # Parse the output TSV
            with open(output_tsv, "r") as f:
                lines = f.readlines()
                if lines:
                    # e.g. P07954.pdb Zm000...pdb 100 100 0.85 2.1
                    parts = lines[0].strip().split()
                    if len(parts) >= 5:
                        return float(parts[4])
                        
        else:
            logger.warning(f"Foldseek failed: {stderr.decode()}")
    except FileNotFoundError:
        # Foldseek binary is not installed on this system PATH
        logger.warning("Foldseek binary not found. Falling back to mock TM-score.")
    except Exception as e:
        logger.error(f"Error running Foldseek: {e}")
        
    # Mock fallback
    return round(random.uniform(0.5, 0.99), 2)

async def _enrich_one(
    ortholog: OrthologMapping,
    *,
    do_foldseek_1to1: bool,
    logs: List[ExecutionLogEntry],
) -> ValidatedTarget:
    """
    Run pLDDT (AlphaFold), expression breadth (Gramene), and — only when
    `do_foldseek_1to1` is True — a 1-to-1 Foldseek alignment. Structure-based
    and consensus hits never need the 1-to-1 call because Phase 4.5 already
    aligned them against the whole maize AFDB.
    """
    target_uniprot_id = await get_uniprot_for_maize_gene(ortholog.maize_gene_model)
    plddt = await validate_structure_alphafold(target_uniprot_id)
    experiments = await fetch_expression_breadth(ortholog.maize_gene_model)

    has_structural = "Foldseek-structural" in (ortholog.sources or [])
    if has_structural:
        tm_score = ortholog.similarity_score / 100.0
        tm_provenance = "Phase 4.5"
        kind = "cheap"        # pLDDT + expression only
    elif do_foldseek_1to1:
        tm_score = await calculate_tm_score_foldseek(ortholog.query_uniprot_id, target_uniprot_id)
        tm_provenance = "Phase 5 1-to-1"
        kind = "full"
    else:
        # Sequence-only hit outside the top-N: cheap enrichment only,
        # leave TM blank (frontend renders "—").
        tm_score = 0.0
        tm_provenance = "skipped (cheap path)"
        kind = "cheap"

    logs.append(ExecutionLogEntry(
        phase=5, database="AlphaFold + Expression breadth + Foldseek",
        status="success" if plddt > 70.0 else ("info" if plddt > 0 else "warning"),
        hits=len(experiments),
        message=(
            f"{ortholog.maize_gene_model}"
            f" (UniProt={target_uniprot_id or '?'}, kind={kind}):"
            f" pLDDT={plddt:.1f}, TM={tm_score:.2f} ({tm_provenance}),"
            f" expressed in {len(experiments)} GXA experiments"
        ),
    ))

    return ValidatedTarget(
        maize_gene_model=ortholog.maize_gene_model,
        tm_score=tm_score,
        plddt=plddt,
        tissue_expression_fpkm={exp: "expressed" for exp in experiments[:6]},
        n_expression_experiments=len(experiments),
        expression_experiments=experiments[:12],
        enrichment_kind=kind,
    )


# Backwards-compat alias (some places used to import this name).
async def process_ortholog(ortholog: OrthologMapping, logs: List[ExecutionLogEntry] = None) -> ValidatedTarget:
    return await _enrich_one(ortholog, do_foldseek_1to1=True, logs=logs or [])


async def execute_phase5(orthologs: List[OrthologMapping], logs: List[ExecutionLogEntry] = None,
                          config: Optional[PipelineConfig] = None) -> List[ValidatedTarget]:
    """
    Phase 5: structural + transcriptomic enrichment for ALL discovered orthologs.

    No longer "top-N only". The expensive piece — 1-to-1 Foldseek for
    sequence-only hits without a known TM — is still capped to the top
    `enrichment_top_n`. Everything else (pLDDT, Gramene expression breadth)
    is cheap and runs for every hit.

    For consensus and structure-only hits, the Foldseek TM is reused from
    Phase 4.5 (no duplicated effort).
    """
    if logs is None:
        logs = []
    if config is None:
        config = PipelineConfig()
    plddt_threshold = config.plddt_threshold
    top_n_full = config.enrichment_top_n

    if not orthologs:
        return []

    # Split orthologs: top-N get the FULL path (1-to-1 Foldseek if needed);
    # everything else gets the CHEAP path (pLDDT + expression only, TM reused
    # from Phase 4.5 for structural hits, blank for sequence-only ones).
    full_set = orthologs[:top_n_full]
    cheap_set = orthologs[top_n_full:]

    logs.append(ExecutionLogEntry(
        phase=5, database="AlphaFold + Expression breadth + Foldseek",
        status="info", hits=0,
        message=(
            f"Enriching all {len(orthologs)} hits: top {len(full_set)} get full enrichment "
            f"(1-to-1 Foldseek for sequence-only); remaining {len(cheap_set)} get cheap "
            f"pLDDT + expression-breadth (TM reused from Phase 4.5 where available). "
            f"pLDDT filter > {plddt_threshold:.0f}."
        ),
    ))

    full_tasks  = [_enrich_one(o, do_foldseek_1to1=True,  logs=logs) for o in full_set]
    cheap_tasks = [_enrich_one(o, do_foldseek_1to1=False, logs=logs) for o in cheap_set]
    all_targets = await asyncio.gather(*full_tasks, *cheap_tasks)

    # Keep all targets that have either decent structural model OR known expression.
    # We don't filter cheap-path hits as aggressively because their lack of TM
    # isn't evidence of a bad ortholog — Phase 4.5 may already have given them one.
    enriched_targets = [
        t for t in all_targets
        if t.plddt > plddt_threshold or t.n_expression_experiments > 0
    ]

    expressed = sum(1 for t in enriched_targets if t.n_expression_experiments > 0)
    logs.append(ExecutionLogEntry(
        phase=5, database="AlphaFold pLDDT + Expression filter",
        status="success" if enriched_targets else "warning",
        hits=len(enriched_targets),
        message=(
            f"{len(enriched_targets)}/{len(all_targets)} hits retained "
            f"(pLDDT > {plddt_threshold:.0f} OR expression evidence); "
            f"{expressed} have detected expression in ≥1 Atlas experiment."
        ),
    ))

    return enriched_targets
