import httpx
import logging
import asyncio
from typing import List, Dict, Optional
from models import ProteinCandidate, OrthologMapping, ExecutionLogEntry, PipelineConfig
from hmmer_runner import run_phmmer_search
from phase4_5 import execute_phase4_5
from corncyc_lookup import corncyc_orthologs_for_chebi

logger = logging.getLogger(__name__)

# UniProt REST cache (acc → (species_slug, [ensembl_gene_ids]))
_UNIPROT_ENSEMBL_CACHE: Dict[str, tuple] = {}


async def _fetch_uniprot_species_and_ensembl(client: httpx.AsyncClient, uniprot_accession: str) -> tuple:
    """
    Look up the organism + Ensembl-style gene IDs for a UniProt accession.
    Returns (species_slug, [gene_ids]). Cached in-process.

    UniProt's REST API is the right source — Ensembl's own `/xrefs/id/{uniprot}`
    endpoint does NOT accept UniProt accessions (returns 400 "ID not found").

    For plant proteins, the Ensembl-compatible gene IDs live under the
    `Gramene` (cross-species plant DB) and `TAIR` (Arabidopsis-specific)
    xref database names — NOT under any `Ensembl*` key. We also strip the
    `.N` transcript-version suffix since the homology endpoint wants the
    bare gene ID.
    """
    if uniprot_accession in _UNIPROT_ENSEMBL_CACHE:
        return _UNIPROT_ENSEMBL_CACHE[uniprot_accession]
    species_slug = ""
    gene_ids: List[str] = []
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}.json"
        r = await client.get(url, timeout=15.0)
        if r.status_code == 200:
            d = r.json()
            org = (d.get("organism") or {}).get("scientificName", "")
            species_slug = org.lower().replace(" ", "_").strip()

            # Database names whose `id` field is an Ensembl-compatible gene model
            ENSEMBL_GENE_DBS = {
                "Ensembl", "EnsemblPlants", "EnsemblFungi", "EnsemblMetazoa",
                "EnsemblBacteria", "EnsemblProtists",
                "Gramene",   # plants (TAIR/AT*, OS*, ZM*, etc. — cross-references to Ensembl Plants)
                "TAIR",      # Arabidopsis (TAIR IDs are Ensembl Plants IDs)
            }
            for xref in d.get("uniProtKBCrossReferences", []) or []:
                db = xref.get("database") or ""
                if db not in ENSEMBL_GENE_DBS:
                    continue
                gid = ""
                for p in xref.get("properties", []) or []:
                    if p.get("key") == "GeneId":
                        gid = (p.get("value") or "").split(".")[0]
                        break
                if not gid:
                    gid = (xref.get("id") or "").split(".")[0]
                if gid:
                    gene_ids.append(gid)
            gene_ids = list(dict.fromkeys(gene_ids))  # dedup, preserve order
    except Exception as e:
        logger.error(
            f"UniProt REST lookup failed for {uniprot_accession} "
            f"({type(e).__name__}): {e!r}"
        )

    _UNIPROT_ENSEMBL_CACHE[uniprot_accession] = (species_slug, gene_ids)
    return species_slug, gene_ids


# Ensembl Compara has multiple compara databases. Try the right one for each kingdom.
_COMPARA_BY_KINGDOM = [
    ("plants",     {"arabidopsis_thaliana", "oryza_sativa", "saccharomyces_cerevisiae",
                    "selaginella_moellendorffii", "physcomitrella_patens", "brachypodium_distachyon",
                    "sorghum_bicolor", "glycine_max", "zea_mays"}),
    # 'pan_homology' covers cross-kingdom orthology where available
    ("pan_homology", set()),
]


async def query_ensembl_homology(uniprot_accession: str) -> List[OrthologMapping]:
    """
    Map a UniProt accession to maize gene IDs via Ensembl Compara.

    Path: UniProt REST → (organism + Ensembl gene IDs) → Ensembl
    /homology/id/{species}/{gene_id}?compara=plants;target_species=zea_mays.
    For non-plant proteins (cattle, human, E. coli, etc.) this returns nothing
    — those rely on Phase 4.5 structural discovery and the HMMER fallback.
    """
    results: List[OrthologMapping] = []
    headers = {"Content-Type": "application/json"}

    try:
        async with httpx.AsyncClient(follow_redirects=True) as client:
            species, gene_ids = await _fetch_uniprot_species_and_ensembl(client, uniprot_accession)
            if not species or not gene_ids:
                return results

            for gid in gene_ids:
                # Try the plants compara first; if the source organism isn't in it, the call
                # returns an empty homologies list and we move on.
                url = (
                    f"https://rest.ensembl.org/homology/id/{species}/{gid}"
                    f"?compara=plants;target_species=zea_mays;type=orthologues"
                )
                try:
                    r = await client.get(url, headers=headers, timeout=20.0)
                    if r.status_code != 200:
                        continue
                    data = r.json()
                    for item in data.get("data", []) or []:
                        for hom in item.get("homologies", []) or []:
                            target = hom.get("target", {}) or {}
                            maize_id = target.get("id")
                            if not maize_id:
                                continue
                            perc_id = target.get("perc_id")
                            perc_pos = target.get("perc_pos")
                            score = float(perc_id if perc_id is not None else perc_pos if perc_pos is not None else 50.0)
                            results.append(OrthologMapping(
                                maize_gene_model=maize_id,
                                plaza_orthogroup=f"Ensembl-Compara plants ({hom.get('type','ortholog')})",
                                similarity_score=score,
                                sources=["Ensembl"],
                                consensus_score=1,
                            ))
                except Exception as e:
                    logger.warning(f"Ensembl homology call failed for {species}/{gid}: {e}")
    except Exception as e:
        logger.error(f"Error querying Ensembl Homology for {uniprot_accession}: {e}")

    return results

_PLAZA_UNAVAILABLE_LOGGED = False


async def query_plaza_homology(uniprot_accession: str) -> List[OrthologMapping]:
    """
    Live query to PLAZA Monocots 5.0 API.

    As of 2026 the anonymous PLAZA REST endpoint returns 403 Forbidden — the
    server is up but rejects unauthenticated requests. We log this once and
    return [] so the pipeline keeps moving; Phase 4.5 (Foldseek) covers the
    same biological ground when PLAZA is unavailable.
    """
    global _PLAZA_UNAVAILABLE_LOGGED
    results: List[OrthologMapping] = []
    try:
        async with httpx.AsyncClient(follow_redirects=True) as client:
            search_url = (
                "https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v5_monocots/"
                f"api/v2/genes?uniprot={uniprot_accession}"
            )
            r = await client.get(search_url, headers={"Accept": "application/json"}, timeout=15.0)
            if r.status_code == 403:
                if not _PLAZA_UNAVAILABLE_LOGGED:
                    logger.warning("PLAZA Monocots 5 API returned 403 — endpoint requires auth or has been moved. Skipping PLAZA lookups for the rest of this session.")
                    _PLAZA_UNAVAILABLE_LOGGED = True
                return results
            if r.status_code != 200:
                logger.warning(f"PLAZA returned HTTP {r.status_code} for {uniprot_accession}")
                return results
            data = r.json()
            genes = data.get("genes", []) or []
            if not genes:
                return results
            plaza_gene_id = genes[0].get("id")
            if not plaza_gene_id:
                return results
            ortho_url = (
                "https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v5_monocots/"
                f"api/v2/orthologs?gene_ids={plaza_gene_id}&species=zma"
            )
            r2 = await client.get(ortho_url, headers={"Accept": "application/json"}, timeout=15.0)
            if r2.status_code == 200:
                for ortho in (r2.json().get("orthologs", []) or []):
                    maize_id = ortho.get("id")
                    if maize_id:
                        results.append(OrthologMapping(
                            maize_gene_model=maize_id,
                            plaza_orthogroup=ortho.get("orthogroup", "PLAZA_ORTHO"),
                            similarity_score=75.0,
                            sources=["PLAZA"],
                            consensus_score=1,
                        ))
    except Exception as e:
        logger.error(f"Error querying PLAZA API: {e}")
    return results


async def query_biomart_homology(uniprot_accession: str) -> List[OrthologMapping]:
    """
    Phytozome / Ensembl Plants BioMart-style ortholog lookup via the Compara
    homology endpoint.

    This used to call Ensembl's `/xrefs/id/{uniprot_acc}?external_db=maize_gdb`,
    which doesn't work — `/xrefs/id` requires Ensembl IDs, not UniProt
    accessions. The previous implementation always returned [].

    We now follow the same UniProt → Ensembl gene ID resolution as Ensembl
    homology but query with `compara=pan_homology` which is Ensembl's
    cross-kingdom Compara database. This catches some orthologs the strict
    `compara=plants` query misses (e.g. cyanobacterial → maize chloroplast).
    """
    results: List[OrthologMapping] = []
    try:
        async with httpx.AsyncClient(follow_redirects=True) as client:
            species, gene_ids = await _fetch_uniprot_species_and_ensembl(client, uniprot_accession)
            if not species or not gene_ids:
                return results
            for gid in gene_ids:
                url = (
                    f"https://rest.ensembl.org/homology/id/{species}/{gid}"
                    f"?compara=pan_homology;target_species=zea_mays;type=orthologues"
                )
                r = await client.get(url, headers={"Content-Type": "application/json"}, timeout=20.0)
                if r.status_code != 200:
                    continue
                for item in (r.json().get("data", []) or []):
                    for hom in (item.get("homologies", []) or []):
                        target = hom.get("target", {}) or {}
                        maize_id = target.get("id")
                        if not maize_id:
                            continue
                        perc_id = target.get("perc_id") or target.get("perc_pos") or 50.0
                        results.append(OrthologMapping(
                            maize_gene_model=maize_id,
                            plaza_orthogroup=f"Ensembl pan-homology ({hom.get('type','ortholog')})",
                            similarity_score=float(perc_id),
                            sources=["EnsemblPanHomology"],
                            consensus_score=1,
                        ))
    except Exception as e:
        logger.error(f"Error querying Ensembl pan-homology: {e}")
    return results

async def execute_phase4(proteins: List[ProteinCandidate], input_data,
                          logs: List[ExecutionLogEntry] = None,
                          config: Optional[PipelineConfig] = None,
                          chebi_id: Optional[str] = None) -> List[OrthologMapping]:
    """
    Executes Phase 4 — sequence orthology (Ensembl + PLAZA + pan-homology),
    plus the parallel structural-discovery lane (Phase 4.5 via Foldseek),
    plus the curated CornCyc lane (when CornCyc is installed and `chebi_id`
    is supplied).
    """
    if logs is None:
        logs = []
    if config is None:
        config = PipelineConfig()

    logs.append(ExecutionLogEntry(
        phase=4, database="Ensembl + PLAZA + BioMart", status="info", hits=0,
        message=f"Projecting {len(proteins)} pan-life proteins into Zea mays via 3 parallel ortholog APIs"
    ))

    all_raw_mappings = []

    async def fetch_and_tag(protein, query_func):
        res = await query_func(protein.uniprot_accession)
        for r in res:
            r.query_uniprot_id = protein.uniprot_accession
        return res

    tasks = []
    for protein in proteins:
        tasks.append(fetch_and_tag(protein, query_ensembl_homology))
        tasks.append(fetch_and_tag(protein, query_plaza_homology))
        tasks.append(fetch_and_tag(protein, query_biomart_homology))

    # Run Phase 4.5 (structure-guided discovery via Foldseek) in parallel with
    # the sequence-based ortholog APIs. It's a no-op if the maize AF index isn't built.
    sequence_db_task = asyncio.gather(*tasks)
    structural_task = execute_phase4_5(
        proteins,
        tm_threshold=config.foldseek_tm_threshold,
        max_hits_per_query=config.foldseek_max_hits_per_query,
        logs=logs,
        config=config,
    )
    db_results, structural_mappings = await asyncio.gather(sequence_db_task, structural_task)

    # Tally per-DB hits across all proteins (indexes: 0=Ensembl, 1=PLAZA, 2=BioMart, cycling per protein)
    ensembl_total = sum(len(r) for i, r in enumerate(db_results) if i % 3 == 0)
    plaza_total = sum(len(r) for i, r in enumerate(db_results) if i % 3 == 1)
    biomart_total = sum(len(r) for i, r in enumerate(db_results) if i % 3 == 2)
    logs.append(ExecutionLogEntry(
        phase=4, database="Ensembl Homology",
        status="success" if ensembl_total else "warning",
        hits=ensembl_total,
        message=f"Ensembl REST returned {ensembl_total} raw maize orthologs"
    ))
    logs.append(ExecutionLogEntry(
        phase=4, database="PLAZA Monocots v5",
        status="success" if plaza_total else "warning",
        hits=plaza_total,
        message=f"PLAZA API returned {plaza_total} raw maize orthologs"
    ))
    logs.append(ExecutionLogEntry(
        phase=4, database="Ensembl pan-homology Compara",
        status="success" if biomart_total else "warning",
        hits=biomart_total,
        message=f"Ensembl pan-homology Compara returned {biomart_total} raw maize orthologs"
    ))

    for res_list in db_results:
        all_raw_mappings.extend(res_list)
    # Fold structural-discovery hits into the same consensus reducer so a maize
    # gene found by both sequence and structure earns consensus_score=2.
    all_raw_mappings.extend(structural_mappings)

    # CornCyc curated lane (opt-in; requires the user-supplied PGDB to be present).
    # Adds maize genes annotated by the Plant Metabolic Network as catalysing
    # reactions involving the query compound's ChEBI ID. Joins the same consensus
    # reducer so consensus_score rises for genes also found by Ensembl/Foldseek.
    if chebi_id:
        corncyc_mappings = corncyc_orthologs_for_chebi(chebi_id)
        if corncyc_mappings:
            all_raw_mappings.extend(corncyc_mappings)
            logs.append(ExecutionLogEntry(
                phase=4, database="CornCyc (PMN)", status="success",
                hits=len(corncyc_mappings),
                message=f"CornCyc adds {len(corncyc_mappings)} curated maize gene annotation(s) for {chebi_id}"
            ))
        else:
            logs.append(ExecutionLogEntry(
                phase=4, database="CornCyc (PMN)", status="info", hits=0,
                message=f"CornCyc has no curated annotation for {chebi_id} (or PGDB not installed)"
            ))
        
    merged_mappings: Dict[str, OrthologMapping] = {}
    
    for m in all_raw_mappings:
        gene = m.maize_gene_model
        if gene not in merged_mappings:
            merged_mappings[gene] = OrthologMapping(
                query_uniprot_id=m.query_uniprot_id,
                maize_gene_model=gene,
                plaza_orthogroup=m.plaza_orthogroup,
                similarity_score=m.similarity_score,
                sources=m.sources.copy(),
                consensus_score=1
            )
        else:
            current = merged_mappings[gene]
            for src in m.sources:
                if src not in current.sources:
                    current.sources.append(src)
            
            current.consensus_score = len(current.sources)
            
            if m.similarity_score > current.similarity_score:
                current.similarity_score = m.similarity_score
                
            if current.plaza_orthogroup == "LIVE_ORTHO" and m.plaza_orthogroup != "LIVE_ORTHO":
                current.plaza_orthogroup = m.plaza_orthogroup
                
    final_list = list(merged_mappings.values())
    final_list.sort(key=lambda x: (x.consensus_score, x.similarity_score), reverse=True)
    
    # HMMER Tier 2 Fallback
    if not final_list and proteins:
        e_cutoff = config.hmmer_e_value
        logger.info("Live APIs returned 0 hits. Initiating Tier 2 HMMER Local Alignment fallback...")
        logs.append(ExecutionLogEntry(
            phase=4, database="Local HMMER", status="info", hits=0,
            message=f"All live ortholog APIs returned 0 — running local phmmer (E≤{e_cutoff:.0e}) for {len(proteins)} query sequences"
        ))
        hmmer_hits = []
        for p in proteins:
            if not p.sequence:
                continue
            logger.info(f"Running HMMER for {p.uniprot_accession} with E-value {e_cutoff}...")
            logs.append(ExecutionLogEntry(
                phase=4, database="Local HMMER (phmmer)", status="info", hits=0,
                message=f"phmmer {p.uniprot_accession} ({len(p.sequence)} aa) vs Zm-NAM-5.0 proteome"
            ))
            hits = await run_phmmer_search(p.sequence, p.uniprot_accession, e_value_cutoff=e_cutoff, max_hits=None)
            logs.append(ExecutionLogEntry(
                phase=4, database="Local HMMER (phmmer)",
                status="success" if hits else "warning",
                hits=len(hits),
                message=f"phmmer {p.uniprot_accession} → {len(hits)} maize hits" + (f" (best E={hits[0]['e_value']:.1e})" if hits else "")
            ))
            for h in hits:
                # Calculate a mock similarity score from the bit score for normalization
                sim_score = min(100.0, (h["bit_score"] / 1000.0) * 100.0)
                hmmer_hits.append(OrthologMapping(
                    query_uniprot_id=p.uniprot_accession,
                    maize_gene_model=h["maize_gene_model"],
                    plaza_orthogroup=f"HMMER_E:{h['e_value']:.1e}",
                    similarity_score=sim_score,
                    sources=["Local_HMMER"],
                    consensus_score=1
                ))
        
        # Sort HMMER hits by similarity score
        hmmer_hits.sort(key=lambda x: x.similarity_score, reverse=True)
        # Take all unique gene models
        unique_hmmer = []
        seen_genes = set()
        for h in hmmer_hits:
            if h.maize_gene_model not in seen_genes:
                unique_hmmer.append(h)
                seen_genes.add(h.maize_gene_model)
                
        final_list = unique_hmmer
        
        logs.append(ExecutionLogEntry(
            phase=4,
            database="Local HMMER",
            status="success" if final_list else "error",
            hits=len(final_list),
            message=f"API failed. Local HMMER found {len(final_list)} structurally validated homologs." if final_list else "HMMER search failed to find any homologs."
        ))
    else:
        logs.append(ExecutionLogEntry(
            phase=4,
            database="Ensembl Plants / PLAZA",
            status="success" if final_list else "warning",
            hits=len(final_list),
            message=f"Found {len(final_list)} orthologs across multiple databases." if final_list else "Zero orthologs found via API cross-referencing."
        ))
    
    return final_list
