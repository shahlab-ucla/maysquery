import httpx
import logging
from typing import List, Tuple, Optional
from models import ChemicalEntity, ReactionNetwork, ExecutionLogEntry, PipelineConfig
from chebi_utils import expand_chebi_conjugate_forms
import urllib.parse

logger = logging.getLogger(__name__)

# Default Phase 2 query cap. Overridable per-request via PipelineConfig.rhea_fetch_limit.
# Bumped from 10 → 100 historically: Rhea reactions are cheap to fetch, and the
# previous LIMIT 10 silently truncated common metabolites (ATP, NAD+, etc.).
RHEA_FETCH_LIMIT = 100

async def query_kegg_pathways(chebi_id: str) -> Tuple[List[str], List[str]]:
    """
    Map ChEBI -> KEGG compound -> KEGG pathways. Returns (pathway_ids, pathway_names).
    Pathway names come from a second KEGG `list/pathway/{id}` call so the UI
    can show 'Citrate cycle (TCA cycle)' instead of just 'path:map00020'.
    """
    numeric_id = chebi_id.split(":")[-1]
    pathway_ids: List[str] = []
    pathway_names: List[str] = []

    try:
        async with httpx.AsyncClient(follow_redirects=True) as client:
            # ChEBI -> KEGG compound ID via the proper cross-reference endpoint.
            # The old `find/compound/{n}` was a text search and could match unrelated
            # compounds whose KEGG ID happened to start with the same digits
            # (e.g. CHEBI:18012 → C18012 Flaviolin instead of C00122 Fumarate).
            r_conv = await client.get(f"https://rest.kegg.jp/conv/compound/chebi:{numeric_id}", timeout=10.0)
            kegg_cpd_id = None
            if r_conv.status_code == 200 and r_conv.text.strip():
                # Each line: "chebi:18012\tcpd:C00122"
                first = r_conv.text.strip().split("\n")[0]
                if "\t" in first:
                    kegg_cpd_id = first.split("\t")[1]  # e.g. "cpd:C00122"

            if not kegg_cpd_id:
                return pathway_ids, pathway_names

            # KEGG compound -> KEGG pathway IDs
            r_link = await client.get(f"https://rest.kegg.jp/link/pathway/{kegg_cpd_id}", timeout=10.0)
            if r_link.status_code == 200:
                for line in r_link.text.strip().split("\n"):
                    if line and "\t" in line:
                        pid = line.split("\t")[1]
                        if pid not in pathway_ids:
                            pathway_ids.append(pid)

            if not pathway_ids:
                return pathway_ids, pathway_names

            # Pathway IDs -> human-readable names. KEGG `list pathway` accepts batched IDs.
            joined = "+".join(p.replace("path:", "") for p in pathway_ids)
            r_names = await client.get(f"https://rest.kegg.jp/list/{joined}", timeout=15.0)
            if r_names.status_code == 200:
                # Each line: "path:mapXXXXX\tHuman readable name"
                id_to_name = {}
                for line in r_names.text.strip().split("\n"):
                    parts = line.split("\t", 1)
                    if len(parts) == 2:
                        key = parts[0] if parts[0].startswith("path:") else f"path:{parts[0]}"
                        id_to_name[key] = parts[1].strip()
                pathway_names = [id_to_name.get(pid, pid) for pid in pathway_ids]
    except Exception as e:
        logger.error(f"Error querying KEGG: {e}")

    return pathway_ids, pathway_names

def _rhea_values_clause(chebi_ids: List[str]) -> str:
    """Build a SPARQL VALUES clause like `VALUES ?target_chebi { chebi:18012 chebi:29806 }`."""
    tokens = " ".join(f"chebi:{cid.split(':')[-1]}" for cid in chebi_ids)
    return f"VALUES ?target_chebi {{ {tokens} }}"


async def _rhea_count(chebi_ids: List[str]) -> int:
    """Cheap COUNT query so we can report 'returned X of Y' truncation."""
    url = "https://sparql.rhea-db.org/sparql"
    values = _rhea_values_clause(chebi_ids)
    query = f"""
    PREFIX rh:<http://rdf.rhea-db.org/>
    PREFIX chebi:<http://purl.obolibrary.org/obo/CHEBI_>
    SELECT (COUNT(DISTINCT ?reaction) AS ?n) WHERE {{
      {values}
      ?reaction rdfs:subClassOf rh:Reaction .
      ?reaction rh:side ?s . ?s rh:contains ?p . ?p rh:compound ?c .
      ?c rh:chebi ?target_chebi .
    }}
    """
    try:
        async with httpx.AsyncClient() as client:
            r = await client.post(
                url,
                data=urllib.parse.urlencode({"query": query}),
                headers={"Accept": "application/sparql-results+json",
                         "Content-Type": "application/x-www-form-urlencoded"},
                timeout=15.0,
            )
            if r.status_code == 200:
                return int(r.json()["results"]["bindings"][0]["n"]["value"])
    except Exception as e:
        logger.warning(f"Rhea COUNT failed: {e}")
    return -1  # unknown


async def query_rhea_sparql(chebi_ids: List[str], limit: int = RHEA_FETCH_LIMIT) -> Tuple[List[ReactionNetwork], int]:
    """
    Query the Rhea SPARQL endpoint for reactions involving ANY of the given
    ChEBI IDs (the seed plus its conjugate-acid/base forms).

    Also fetches `equation` (human-readable text), `label`, and EC numbers
    via GROUP_CONCAT so reactions show meaningful titles instead of bare
    'RHEA:12345' identifiers.
    Returns (reactions, total_available).
    """
    url = "https://sparql.rhea-db.org/sparql"
    values = _rhea_values_clause(chebi_ids)
    query = f"""
    PREFIX rh:<http://rdf.rhea-db.org/>
    PREFIX chebi:<http://purl.obolibrary.org/obo/CHEBI_>
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>

    SELECT DISTINCT ?reaction ?equation ?label ?isTransport ?isBalanced
                    (GROUP_CONCAT(DISTINCT ?ec; separator="|") AS ?ecs)
    WHERE {{
      {values}
      ?reaction rdfs:subClassOf rh:Reaction .
      ?reaction rh:equation ?equation .
      ?reaction rh:isTransport ?isTransport .
      ?reaction rh:isChemicallyBalanced ?isBalanced .
      OPTIONAL {{ ?reaction rdfs:label ?label . }}
      OPTIONAL {{ ?reaction rh:ec ?ec . }}

      ?reaction rh:side ?reactionSide .
      ?reactionSide rh:contains ?participant .
      ?participant rh:compound ?compound .
      ?compound rh:chebi ?target_chebi .
    }}
    GROUP BY ?reaction ?equation ?label ?isTransport ?isBalanced
    LIMIT {int(limit)}
    """

    results: List[ReactionNetwork] = []
    total = -1
    try:
        async with httpx.AsyncClient() as client:
            encoded_query = urllib.parse.urlencode({"query": query})
            headers = {
                "Accept": "application/sparql-results+json",
                "Content-Type": "application/x-www-form-urlencoded",
            }
            response = await client.post(url, data=encoded_query, headers=headers, timeout=30.0)

            if response.status_code == 200:
                bindings = response.json().get("results", {}).get("bindings", [])
                for b in bindings:
                    r_id = b["reaction"]["value"].split("/")[-1]
                    eq = (b.get("equation") or {}).get("value") or None
                    lbl = (b.get("label") or {}).get("value") or None
                    ecs_raw = (b.get("ecs") or {}).get("value") or ""
                    ecs = [e.rsplit("/", 1)[-1] for e in ecs_raw.split("|") if e]
                    results.append(ReactionNetwork(
                        rhea_id=f"RHEA:{r_id}",
                        pathway_names=[],
                        is_transport=b["isTransport"]["value"] == "true",
                        is_balanced=b["isBalanced"]["value"] == "true",
                        equation=eq,
                        label=lbl,
                        ec_numbers=ecs,
                    ))
            else:
                logger.warning(f"Rhea SPARQL returned {response.status_code}: {response.text[:200]}")
    except Exception as e:
        logger.error(f"Error querying Rhea: {e}")

    if len(results) >= limit:
        total = await _rhea_count(chebi_ids)

    return results, total

async def execute_phase2(chemical_entity: ChemicalEntity, logs: List[ExecutionLogEntry] = None,
                         config: Optional[PipelineConfig] = None) -> List[ReactionNetwork]:
    if logs is None:
        logs = []
    if config is None:
        config = PipelineConfig()

    """
    Executes Phase 2: ChEBI -> KEGG Pathways & Rhea Reactions.
    Completely drops all mock fallbacks in favor of live DB traversal.
    """
    logs.append(ExecutionLogEntry(
        phase=2, database="KEGG REST", status="info", hits=0,
        message=f"Looking up KEGG pathways for {chemical_entity.chebi_id}"
    ))
    pathway_ids, pathway_names = await query_kegg_pathways(chemical_entity.chebi_id)
    logs.append(ExecutionLogEntry(
        phase=2, database="KEGG REST",
        status="success" if pathway_ids else "warning",
        hits=len(pathway_ids),
        message=f"KEGG returned {len(pathway_ids)} pathway IDs"
            + (f" (e.g. '{pathway_names[0]}')" if pathway_names else "")
    ))

    # Expand the seed ChEBI ID to include conjugate-acid/base forms.
    # Critical because PubChem returns neutral parents (e.g. fumaric acid
    # CHEBI:18012) while Rhea uses ionised forms (fumarate(2-) CHEBI:29806).
    chebi_family = await expand_chebi_conjugate_forms(chemical_entity.chebi_id, max_depth=config.chebi_expansion_depth)
    if not chebi_family:
        chebi_family = [chemical_entity.chebi_id]
    expansion_note = (
        f" (expanded to {len(chebi_family)} conjugate forms: {', '.join(chebi_family)})"
        if len(chebi_family) > 1 else ""
    )
    logs.append(ExecutionLogEntry(
        phase=2, database="Rhea SPARQL", status="info", hits=0,
        message=f"Executing SPARQL for chemically-balanced reactions involving {chemical_entity.chebi_id}{expansion_note}"
    ))

    reactions, total_available = await query_rhea_sparql(chebi_family, limit=config.rhea_fetch_limit)

    # Pathway info is per-metabolite (KEGG query), but we copy onto each reaction
    # so the existing template code doesn't break.
    for r in reactions:
        r.pathway_names = pathway_names
        r.pathway_ids = pathway_ids

    if not reactions:
        logger.warning(f"No valid chemical reactions found in Rhea for {chemical_entity.chebi_id}.")

    if total_available > len(reactions):
        truncation_msg = f" (showing top {len(reactions)} of {total_available} available — bump RHEA_FETCH_LIMIT to see more)"
    elif total_available == len(reactions):
        truncation_msg = f" (all {total_available} returned)"
    else:
        truncation_msg = ""

    logs.append(ExecutionLogEntry(
        phase=2,
        database="Rhea SPARQL",
        status="success" if reactions else "warning",
        hits=len(reactions),
        message=f"Rhea returned {len(reactions)} reaction networks{truncation_msg}"
            + (f". Top: {reactions[0].rhea_id}" if reactions else "")
    ))

    return reactions
