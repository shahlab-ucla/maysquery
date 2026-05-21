import httpx
import logging
from typing import List, Tuple, Optional
from models import ChemicalEntity, ReactionNetwork, ProteinCandidate, GOTerm, ExecutionLogEntry, PipelineConfig
from chebi_utils import expand_chebi_conjugate_forms
import asyncio

logger = logging.getLogger(__name__)

# Default per-category cap; overridable per-request via PipelineConfig.uniprot_size_per_category.
UNIPROT_SIZE_PER_CATEGORY = 25


async def query_uniprot_category(query_str: str, category_name: str,
                                  size: int = UNIPROT_SIZE_PER_CATEGORY) -> Tuple[List[ProteinCandidate], int]:
    """
    Query UniProt REST API and return (proteins, total_matching).
    `total_matching` comes from the X-Total-Results header so callers can
    report 'showing N of M' truncation.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"

    params = {
        "query": query_str,
        "format": "json",
        "size": int(size),
        "fields": "accession,sequence,go_id",
    }

    results: List[ProteinCandidate] = []
    total = -1
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(url, params=params, timeout=20.0)
            if response.status_code == 200:
                try:
                    total = int(response.headers.get("X-Total-Results", "-1"))
                except ValueError:
                    total = -1
                data = response.json()
                for entry in data.get("results", []):
                    accession = entry.get("primaryAccession", "Unknown")
                    sequence = entry.get("sequence", {}).get("value", "")

                    go_terms = []
                    for db_ref in entry.get("uniProtKBCrossReferences", []):
                        if db_ref.get("database") == "GO":
                            go_id = db_ref.get("id")
                            go_name = "Unknown"
                            for prop in db_ref.get("properties", []):
                                if prop.get("key") == "GoTerm":
                                    raw_name = prop.get("value", "")
                                    go_name = raw_name.split(":", 1)[1] if ":" in raw_name else raw_name
                                    break
                            go_terms.append({"id": go_id, "name": go_name})

                    results.append(ProteinCandidate(
                        uniprot_accession=accession,
                        sequence=sequence,
                        go_terms=go_terms,
                        category=category_name,
                    ))
            else:
                logger.warning(f"UniProt API returned {response.status_code} for {query_str}")
    except Exception as e:
        logger.error(f"Error querying UniProt for {category_name}: {e}")

    return results, total

def _chebi_or_clause(chebi_ids: List[str]) -> str:
    """Build a UniProt query fragment: (chebi:"CHEBI:..." OR chebi:"CHEBI:...")."""
    inner = " OR ".join(f'chebi:"{c}"' for c in chebi_ids)
    return f"({inner})"


def _trunc_note(returned: int, total: int) -> str:
    if total < 0:                  # header missing
        return ""
    if total > returned:
        return f" (showing top {returned} of {total} matching — bump UNIPROT_SIZE_PER_CATEGORY for more)"
    return f" (all {total} returned)"


async def execute_phase3(chemical_entity: ChemicalEntity, reactions: List[ReactionNetwork],
                          logs: List[ExecutionLogEntry] = None,
                          config: Optional[PipelineConfig] = None) -> List[ProteinCandidate]:
    if logs is None:
        logs = []
    if config is None:
        config = PipelineConfig()
    """
    Phase 3: Reaction → Protein Pool.

    Queries UniProt for Enzymes (via Rhea IDs from Phase 2, or ChEBI fallback),
    Transporters (KW-0813), and Receptors (KW-0675) in parallel. ChEBI queries
    are expanded to include the full conjugate-acid/base family so we don't
    miss proteins annotated against a different protonation state than the
    one PubChem returned.
    """
    chebi_family = await expand_chebi_conjugate_forms(chemical_entity.chebi_id, max_depth=config.chebi_expansion_depth)
    if not chebi_family:
        chebi_family = [chemical_entity.chebi_id]
    chebi_clause = _chebi_or_clause(chebi_family)
    family_note = (
        f" (CHEBI family: {', '.join(chebi_family)})"
        if len(chebi_family) > 1 else ""
    )

    # 1. Enzymes: prefer Rhea IDs from Phase 2 (more specific); else ChEBI+keyword
    rhea_queries = [f'rhea:"{r.rhea_id.lower()}"' for r in reactions]
    if rhea_queries:
        rhea_combined = " OR ".join(rhea_queries)
        enzyme_query = f"({rhea_combined}) AND reviewed:true"
        enzyme_desc = f"{len(rhea_queries)} Rhea IDs"
    else:
        enzyme_query = f"{chebi_clause} AND keyword:KW-0255 AND reviewed:true"
        enzyme_desc = f"ChEBI fallback{family_note}"

    transport_query = f"{chebi_clause} AND keyword:KW-0813 AND reviewed:true"
    receptor_query  = f"{chebi_clause} AND keyword:KW-0675 AND reviewed:true"

    size = config.uniprot_size_per_category
    logs.append(ExecutionLogEntry(
        phase=3, database="UniProtKB Swiss-Prot", status="info", hits=0,
        message=(
            f"Querying UniProt for Enzymes ({enzyme_desc}), Transporters (KW-0813), "
            f"Receptors (KW-0675) in parallel — top {size}/category{family_note}"
        ),
    ))

    results = await asyncio.gather(
        query_uniprot_category(enzyme_query, "Enzyme", size=size),
        query_uniprot_category(transport_query, "Transporter", size=size),
        query_uniprot_category(receptor_query, "Receptor", size=size),
    )
    (enzyme_list, enz_total), (transporter_list, tra_total), (receptor_list, rec_total) = results

    logs.append(ExecutionLogEntry(
        phase=3, database="UniProtKB Swiss-Prot",
        status="success" if enzyme_list else "warning",
        hits=len(enzyme_list),
        message=f"Enzyme query returned {len(enzyme_list)} reviewed entries{_trunc_note(len(enzyme_list), enz_total)}",
    ))
    logs.append(ExecutionLogEntry(
        phase=3, database="UniProtKB Swiss-Prot",
        status="success" if transporter_list else "warning",
        hits=len(transporter_list),
        message=f"Transporter query returned {len(transporter_list)} reviewed entries{_trunc_note(len(transporter_list), tra_total)}",
    ))
    logs.append(ExecutionLogEntry(
        phase=3, database="UniProtKB Swiss-Prot",
        status="success" if receptor_list else "warning",
        hits=len(receptor_list),
        message=f"Receptor query returned {len(receptor_list)} reviewed entries{_trunc_note(len(receptor_list), rec_total)}",
    ))

    proteins = [p for sublist in (enzyme_list, transporter_list, receptor_list) for p in sublist]

    unique_proteins = []
    seen = set()
    for p in proteins:
        if p.uniprot_accession not in seen:
            seen.add(p.uniprot_accession)
            unique_proteins.append(p)

    if not unique_proteins:
        logger.warning(f"No proteins found in UniProt for {chemical_entity.chebi_id}")

    logs.append(ExecutionLogEntry(
        phase=3, database="UniProtKB",
        status="success" if unique_proteins else "warning",
        hits=len(unique_proteins),
        message=f"Deduplicated to {len(unique_proteins)} unique proteins",
    ))
    return unique_proteins


async def execute_phase3_by_ec(ec_number: str, logs: List[ExecutionLogEntry] = None,
                                config: Optional[PipelineConfig] = None) -> List[ProteinCandidate]:
    if logs is None:
        logs = []
    if config is None:
        config = PipelineConfig()
    """Phase 3: Skips Phase 1 & 2 by directly querying UniProt with an EC number."""
    size = config.uniprot_size_per_category
    ec_query = f'ec:{ec_number} AND reviewed:true'
    logs.append(ExecutionLogEntry(
        phase=3, database="UniProtKB Swiss-Prot", status="info", hits=0,
        message=f"Direct EC query: ec:{ec_number} AND reviewed:true (top {size})",
    ))
    proteins, total = await query_uniprot_category(ec_query, f"EC:{ec_number}", size=size)

    if not proteins:
        logger.warning(f"No proteins found in UniProt for EC Number {ec_number}")

    logs.append(ExecutionLogEntry(
        phase=3, database="UniProtKB Swiss-Prot",
        status="success" if proteins else "error",
        hits=len(proteins),
        message=f"Fetched {len(proteins)} pan-taxonomic proteins for EC:{ec_number}{_trunc_note(len(proteins), total)}",
    ))
    return proteins
