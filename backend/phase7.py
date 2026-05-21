import httpx
import logging
from typing import List, Optional
from models import ValidatedTarget, AdvancedHomologyTarget, EnsemblOrtholog, ExecutionLogEntry, PipelineConfig

logger = logging.getLogger(__name__)

async def fetch_ensembl_homologs(maize_gene_id: str, max_results: int = 10) -> List[EnsemblOrtholog]:
    url = f"http://rest.ensembl.org/homology/id/zea_mays/{maize_gene_id}?compara=plants"
    orthologs = []
    
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(url, headers={"Accept": "application/json"}, timeout=15.0)
            print(f"[Phase 7] HTTP {response.status_code} for {url}")
            if response.status_code == 200:
                data = response.json()
                if data and "data" in data and len(data["data"]) > 0:
                    homologies = data["data"][0].get("homologies", [])
                    for h in homologies:
                        # Only take orthologs (ignore within-species paralogs if needed, or keep top ones)
                        target = h.get("target", {})
                        if target:
                            orthologs.append(EnsemblOrtholog(
                                species=target.get("species", "unknown"),
                                gene_id=target.get("id", ""),
                                protein_id=target.get("protein_id", ""),
                                percent_identity=target.get("perc_id", 0.0)
                            ))
                    print(f"[Phase 7] {maize_gene_id} fetched {len(orthologs)} orthologs")
            else:
                print(f"[Phase 7] Error {response.status_code}: {response.text[:200]}")
                logger.warning(f"Ensembl API returned {response.status_code} for {maize_gene_id}")
    except Exception as e:
        print(f"[Phase 7] Exception: {e}")
        logger.error(f"Error querying Ensembl API for {maize_gene_id}: {e}")
        
    # Sort by percent identity descending
    orthologs.sort(key=lambda x: x.percent_identity, reverse=True)
    return orthologs[:max_results]

async def execute_phase7(targets: List[ValidatedTarget], logs: List[ExecutionLogEntry] = None,
                          config: Optional[PipelineConfig] = None) -> List[AdvancedHomologyTarget]:
    if logs is None:
        logs = []
    if config is None:
        config = PipelineConfig()
    max_orthologs = config.compara_max_orthologs_per_target

    logs.append(ExecutionLogEntry(
        phase=7, database="Ensembl Compara (plants)", status="info", hits=0,
        message=f"Cross-checking {len(targets)} validated targets against pan-plant Compara orthologs (top {max_orthologs}/target)"
    ))

    advanced_targets = []
    print(f"[Phase 7] Executing on {len(targets)} targets")

    total_ensembl_hits = 0
    for t in targets:
        print(f"[Phase 7] Processing {t.maize_gene_model}")
        # Ensembl Compara check
        ensembl_orthologs = await fetch_ensembl_homologs(t.maize_gene_model, max_results=max_orthologs)
        total_ensembl_hits += len(ensembl_orthologs)

        logs.append(ExecutionLogEntry(
            phase=7, database="Ensembl Compara (plants)",
            status="success" if ensembl_orthologs else "warning",
            hits=len(ensembl_orthologs),
            message=f"{t.maize_gene_model} → {len(ensembl_orthologs)} pan-plant orthologs"
            + (f" (top: {ensembl_orthologs[0].species} {ensembl_orthologs[0].percent_identity:.1f}%)" if ensembl_orthologs else "")
        ))

        # PLAZA Synteny check (mock removed)
        plaza_orthogroup = None
        plaza_synteny_blocks = []

        advanced_targets.append(AdvancedHomologyTarget(
            maize_gene_model=t.maize_gene_model,
            ensembl_orthologs=ensembl_orthologs,
            plaza_synteny_blocks=plaza_synteny_blocks,
            plaza_orthogroup=plaza_orthogroup
        ))

    logs.append(ExecutionLogEntry(
        phase=7,
        database="Ensembl Compara (plants)",
        status="success" if advanced_targets else "warning",
        hits=total_ensembl_hits,
        message=f"Phase 7 complete: {len(targets)} targets mapped to {total_ensembl_hits} pan-plant orthologs" if targets else "No targets reached Phase 7 (Phase 5 returned nothing)"
    ))

    return advanced_targets
