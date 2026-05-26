"""
phase4_5.py — Structure-guided ortholog discovery.

Where Phase 4 finds maize orthologs by sequence homology (Ensembl/PLAZA/BioMart
+ optional HMMER fallback), Phase 4.5 runs Foldseek `easy-search` against the
indexed maize AlphaFold proteome. This rescues "hidden" orthologs whose primary
sequence has diverged below standard detection thresholds but whose 3D fold is
still recognizable.

Requires the maize AF + foldseek index built by `install_maize_afdb.py`.
If the index is missing, this phase emits a warning log and returns no
mappings — the rest of the pipeline still runs.
"""

from __future__ import annotations

import asyncio
import logging
import os
import re
import subprocess
from typing import List, Optional, Dict

import httpx

from models import ProteinCandidate, OrthologMapping, ExecutionLogEntry, PipelineConfig
from install_maize_afdb import (
    is_db_ready, DB_PREFIX, to_wsl_path,
)
from phase5 import download_alphafold_pdb

logger = logging.getLogger(__name__)

# AlphaFold PDB filenames look like AF-<UniProt>-F<frag>-model_v<ver>.pdb
_AF_FILE_RE = re.compile(r"AF-([A-Z0-9]+)-F\d+-model_v\d+")

# Default semaphore (overridden per-call by PipelineConfig.foldseek_concurrency)
_FOLDSEEK_SEMAPHORE = asyncio.Semaphore(2)
_FOLDSEEK_SEMAPHORE_SIZE = 2


def _get_semaphore(desired_size: int) -> asyncio.Semaphore:
    """Lazily rebuild the module-global semaphore when concurrency changes."""
    global _FOLDSEEK_SEMAPHORE, _FOLDSEEK_SEMAPHORE_SIZE
    desired_size = max(1, int(desired_size))
    if desired_size != _FOLDSEEK_SEMAPHORE_SIZE:
        _FOLDSEEK_SEMAPHORE = asyncio.Semaphore(desired_size)
        _FOLDSEEK_SEMAPHORE_SIZE = desired_size
    return _FOLDSEEK_SEMAPHORE


async def _foldseek_search_against_maize(
    query_pdb_path: str,
    query_id: str,
    tm_threshold: float,
    max_hits: int,
) -> List[Dict]:
    """
    Run `wsl foldseek easy-search` against the indexed maize AFDB.

    Search conditioning (whole-chain orthology discovery):
      * `--alignment-type 2`  → TM-align mode. Required for *full-protein*
        TM-scores (`qtmscore` normalised by query length, `ttmscore` normalised
        by target length). The default mode (3Di Smith-Waterman) only gives
        `alntmscore` which is normalised by aligned-region length and tends to
        over-score short partial alignments.
      * `-e 10`               → loose E-value upstream filter; we post-filter
        by TM-score downstream.
      * `--max-seqs 300`      → Foldseek default; passes top-300 prefilter
        candidates to the alignment stage.
      * Post-filter: `max(qtmscore, ttmscore) >= tm_threshold`. This catches
        both true full-length orthologs (both TMs high) AND the case where
        the query is a sub-domain of a larger maize protein (qtmscore high,
        ttmscore low) — at the cost of also catching the inverse (query has
        extra domains; qtmscore low, ttmscore high). Both are useful leads.

    What this CANNOT find: a query enzyme whose single conserved catalytic
    domain is fused with a totally different maize partner, where the
    catalytic-domain TM is high but the chain-level TM is below threshold
    because the rest of the proteins disagree. That requires domain-sliced
    re-search (Phase 6's territory) and isn't done at discovery time yet.

    Returns deduplicated, TM-filtered list of:
      {target_uniprot, qtm, ttm, tm_score (=max of the two), rmsd, evalue,
       prob, lddt, fident, qlen, tlen, alnlen}
    """
    if not is_db_ready():
        return []

    tmp_dir = os.path.join(os.path.dirname(__file__), "tmpFolder")
    os.makedirs(tmp_dir, exist_ok=True)
    out_tsv = os.path.join(tmp_dir, f"{query_id}_struct_hits.tsv")
    fs_tmp = os.path.join(tmp_dir, f"fs_tmp_struct_{query_id}")
    os.makedirs(fs_tmp, exist_ok=True)

    # Column order MUST match parser indexing below.
    FMT_COLS = "query,target,qlen,tlen,alnlen,qtmscore,ttmscore,alntmscore,rmsd,evalue,prob,lddt,fident"

    cmd = [
        "wsl", "foldseek", "easy-search",
        to_wsl_path(query_pdb_path),
        to_wsl_path(DB_PREFIX),
        to_wsl_path(out_tsv),
        to_wsl_path(fs_tmp),
        "--alignment-type", "2",          # TM-align: gives full-protein qtmscore/ttmscore
        "--format-output", FMT_COLS,
        "-e", "10",
        "--max-seqs", "300",
    ]

    async with _FOLDSEEK_SEMAPHORE:  # populated by execute_phase4_5 via _get_semaphore
        try:
            res = await asyncio.to_thread(
                subprocess.run, cmd,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600,
            )
        except subprocess.TimeoutExpired:
            logger.error(f"foldseek easy-search timed out for {query_id}")
            return []
        except FileNotFoundError as e:
            logger.error(f"foldseek launch failed for {query_id}: {e!r}")
            return []
        except Exception as e:
            logger.error(f"foldseek error for {query_id}: {type(e).__name__}: {e}")
            return []

    if res.returncode != 0:
        logger.error(
            f"foldseek easy-search failed for {query_id} "
            f"(rc={res.returncode}): {res.stderr.decode(errors='replace')[:400]}"
        )
        return []

    hits: List[Dict] = []
    seen: set[str] = set()
    if not os.path.exists(out_tsv):
        return hits

    # Column indices into the TSV (must match FMT_COLS above)
    QUERY, TARGET, QLEN, TLEN, ALNLEN, QTMSCORE, TTMSCORE, ALNTMSCORE, RMSD, EVALUE, PROB, LDDT, FIDENT = range(13)

    with open(out_tsv, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 13:
                continue
            target_name = parts[TARGET]
            m = _AF_FILE_RE.search(target_name)
            if not m:
                continue
            uniprot = m.group(1)
            if uniprot in seen:
                continue
            try:
                qtm = float(parts[QTMSCORE])
                ttm = float(parts[TTMSCORE])
            except ValueError:
                continue
            tm_score = max(qtm, ttm)
            if tm_score < tm_threshold:
                continue
            seen.add(uniprot)

            def _f(idx):
                try: return float(parts[idx])
                except ValueError: return None
            def _i(idx):
                try: return int(parts[idx])
                except ValueError: return None

            hits.append({
                "target_uniprot": uniprot,
                "tm_score": tm_score,
                "qtm": qtm,
                "ttm": ttm,
                "rmsd": _f(RMSD),
                "evalue": _f(EVALUE),
                "prob": _f(PROB),
                "lddt": _f(LDDT),
                "fident": _f(FIDENT),
                "qlen": _i(QLEN),
                "tlen": _i(TLEN),
                "alnlen": _i(ALNLEN),
            })
            if len(hits) >= max_hits:
                break

    # Sort by qtm desc so the strongest query-coverage hits rank first
    hits.sort(key=lambda h: h["qtm"], reverse=True)

    # Cleanup scratch
    try:
        os.remove(out_tsv)
    except OSError:
        pass
    try:
        import shutil
        shutil.rmtree(fs_tmp, ignore_errors=True)
    except Exception:
        pass

    return hits


# In-process cache for UniProt → maize gene-model lookups (cheap, can stale across restarts)
_UNIPROT_TO_GENE_CACHE: Dict[str, Optional[str]] = {}


async def _uniprot_to_maize_gene(client: httpx.AsyncClient, uniprot_id: str) -> Optional[str]:
    """
    Map a maize UniProt accession to its Zm00001eb* gene model.

    The previous implementation used `rest.ensembl.org/xrefs/id/{uniprot}`,
    which does NOT accept UniProt accessions (400 'ID not found') — same bug
    as Phase 4's old Ensembl path. We now use UniProt's own REST entry and
    look for the `EnsemblPlants`, `Gramene`, or `MaizeGDB` cross-references
    whose `properties.GeneId` is the Zm gene model.
    """
    if uniprot_id in _UNIPROT_TO_GENE_CACHE:
        return _UNIPROT_TO_GENE_CACHE[uniprot_id]

    gene_id: Optional[str] = None
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        resp = await client.get(url, timeout=15.0)
        if resp.status_code == 200:
            d = resp.json()
            for xref in d.get("uniProtKBCrossReferences", []) or []:
                db = xref.get("database") or ""
                if db not in ("EnsemblPlants", "Gramene", "MaizeGDB"):
                    continue
                # Properties usually carry GeneId; for MaizeGDB the id itself is the gene
                for p in xref.get("properties", []) or []:
                    if p.get("key") == "GeneId":
                        v = (p.get("value") or "").split(".")[0]
                        if v.startswith("Zm0"):
                            gene_id = v
                            break
                if gene_id:
                    break
                # MaizeGDB xref puts the gene model in `id` directly
                rid = (xref.get("id") or "").split(".")[0]
                if db == "MaizeGDB" and rid.startswith("Zm0"):
                    gene_id = rid
                    break
    except Exception as e:
        logger.error(f"UniProt→gene lookup failed for {uniprot_id}: {e}")

    _UNIPROT_TO_GENE_CACHE[uniprot_id] = gene_id
    return gene_id


async def execute_phase4_5(
    proteins: List[ProteinCandidate],
    tm_threshold: float = 0.5,
    max_hits_per_query: int = 10,
    logs: Optional[List[ExecutionLogEntry]] = None,
    config: Optional[PipelineConfig] = None,
) -> List[OrthologMapping]:
    """
    For each pan-life query protein with an AlphaFold model, search the
    indexed maize AlphaFold proteome with foldseek and return TM-filtered
    `OrthologMapping`s.

    Returns [] (logging a warning) when the maize DB hasn't been built yet,
    so the rest of the pipeline keeps working.
    """
    if logs is None:
        logs = []
    if config is not None:
        # Config overrides any direct kwargs the caller passed
        tm_threshold = config.foldseek_tm_threshold
        max_hits_per_query = config.foldseek_max_hits_per_query
        _get_semaphore(config.foldseek_concurrency)

    if not proteins:
        return []

    if not is_db_ready():
        logs.append(ExecutionLogEntry(
            phase=4, database="Foldseek structural discovery",
            status="warning", hits=0,
            message=(
                "Maize AlphaFold index not built — Phase 4.5 disabled. "
                "Run `python install_maize_afdb.py` or POST /api/maize_afdb/install to enable."
            ),
        ))
        return []

    logs.append(ExecutionLogEntry(
        phase=4, database="Foldseek structural discovery",
        status="info", hits=0,
        message=(
            f"Whole-chain TM-align (alignment-type 2) vs indexed maize AlphaFold proteome: "
            f"{len(proteins)} query structures, keep hits with max(qTM,tTM)≥{tm_threshold}, "
            f"top {max_hits_per_query}/query. "
            "Note: single-domain matches in multi-domain queries may rank low — "
            "Phase 6 handles domain re-search for already-enriched hits."
        ),
    ))

    tmp_dir = os.path.join(os.path.dirname(__file__), "tmpFolder")
    os.makedirs(tmp_dir, exist_ok=True)

    # Step 1: download query PDBs + run foldseek searches in parallel (bounded)
    async def search_one(p: ProteinCandidate):
        query_pdb = await download_alphafold_pdb(p.uniprot_accession, tmp_dir)
        if not query_pdb:
            logs.append(ExecutionLogEntry(
                phase=4, database="Foldseek structural discovery",
                status="warning", hits=0,
                message=f"No AlphaFold model for {p.uniprot_accession} — skipping structural search",
            ))
            return p.uniprot_accession, []
        hits = await _foldseek_search_against_maize(
            query_pdb, p.uniprot_accession, tm_threshold, max_hits_per_query,
        )
        msg = f"{p.uniprot_accession} → {len(hits)} structural maize hits"
        if hits:
            h = hits[0]
            msg += (
                f" (top qTM={h['qtm']:.2f} tTM={h['ttm']:.2f}"
                + (f" prob={h['prob']:.2f}" if h.get('prob') is not None else "")
                + (f" lDDT={h['lddt']:.2f}" if h.get('lddt') is not None else "")
                + ")"
            )
        logs.append(ExecutionLogEntry(
            phase=4, database="Foldseek structural discovery",
            status="success" if hits else "warning",
            hits=len(hits),
            message=msg,
        ))
        return p.uniprot_accession, hits

    per_protein = await asyncio.gather(*(search_one(p) for p in proteins))

    # Collect every unique target UniProt to resolve in one batch
    all_targets = {h["target_uniprot"] for _, hs in per_protein for h in hs}

    async with httpx.AsyncClient() as client:
        gene_map = await asyncio.gather(
            *(_uniprot_to_maize_gene(client, u) for u in all_targets)
        )
    uniprot_to_gene = dict(zip(all_targets, gene_map))

    # Step 2: convert hits → OrthologMapping
    # similarity_score is stored as max(qTM, tTM) * 100 to fit the existing
    # 0..100 schema; the qTM/tTM breakdown + LDDT/prob are encoded into
    # plaza_orthogroup so the UI/CSV/HTML can surface the full provenance.
    mappings: List[OrthologMapping] = []
    for query_uniprot, hits in per_protein:
        for h in hits:
            t_uniprot = h["target_uniprot"]
            gene = uniprot_to_gene.get(t_uniprot) or f"UNIPROT:{t_uniprot}"
            qtm = h.get("qtm")
            ttm = h.get("ttm")
            prob = h.get("prob")
            lddt = h.get("lddt")
            tag_bits = [f"Foldseek qTM={qtm:.2f}"]
            if ttm is not None:  tag_bits.append(f"tTM={ttm:.2f}")
            if prob is not None: tag_bits.append(f"p={prob:.2f}")
            if lddt is not None: tag_bits.append(f"lDDT={lddt:.2f}")
            mappings.append(OrthologMapping(
                query_uniprot_id=query_uniprot,
                maize_gene_model=gene,
                plaza_orthogroup=" ".join(tag_bits),
                similarity_score=round(h["tm_score"] * 100, 1),
                sources=["Foldseek-structural"],
                consensus_score=1,
                source_evidence={"Foldseek-structural": {
                    "method": "Foldseek 3D structural alignment vs maize AlphaFold proteome",
                    "query_uniprot": query_uniprot,
                    "target_uniprot": t_uniprot,
                    "qtm": qtm, "ttm": ttm,
                    "tm_score": h["tm_score"],
                    "prob": prob, "lddt": lddt,
                    "rmsd": h.get("rmsd"),
                    "evalue": h.get("evalue"),
                    "alnlen": h.get("alnlen"),
                }},
            ))

    logs.append(ExecutionLogEntry(
        phase=4, database="Foldseek structural discovery",
        status="success" if mappings else "warning",
        hits=len(mappings),
        message=f"Phase 4.5 complete: {len(mappings)} structural ortholog candidates "
                f"across {len(all_targets)} unique maize targets",
    ))

    return mappings
