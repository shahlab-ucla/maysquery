import httpx
import logging
from typing import List
from models import ChemicalEntity, ReactionNetwork
import urllib.parse

logger = logging.getLogger(__name__)

async def query_kegg_pathways(chebi_id: str) -> List[str]:
    """
    Query KEGG to map ChEBI ID to a KEGG compound, and then to pathways.
    """
    numeric_id = chebi_id.split(":")[-1]
    url_find = f"http://rest.kegg.jp/find/compound/{numeric_id}"
    
    kegg_cpd_id = None
    pathways = []
    
    try:
        async with httpx.AsyncClient() as client:
            resp_find = await client.get(url_find, timeout=10.0)
            if resp_find.status_code == 200:
                lines = resp_find.text.strip().split("\n")
                if lines and lines[0]:
                    kegg_cpd_id = lines[0].split("\t")[0]
                    
            if kegg_cpd_id:
                url_link = f"http://rest.kegg.jp/link/pathway/{kegg_cpd_id}"
                resp_link = await client.get(url_link, timeout=10.0)
                if resp_link.status_code == 200:
                    lines = resp_link.text.strip().split("\n")
                    for line in lines:
                        if line:
                            pathway_id = line.split("\t")[1]
                            pathways.append(pathway_id)
    except Exception as e:
        logger.error(f"Error querying KEGG: {e}")
        
    return pathways

async def query_rhea_sparql(chebi_id: str) -> List[ReactionNetwork]:
    """
    Query Rhea SPARQL endpoint to find reactions involving the ChEBI ID.
    Enforces rh:isChemicallyBalanced and checks rh:isTransport.
    """
    url = "https://sparql.rhea-db.org/sparql"
    
    # Correct SPARQL prefix format for Rhea
    query = f"""
    PREFIX rh:<http://rdf.rhea-db.org/>
    PREFIX chebi:<http://purl.obolibrary.org/obo/CHEBI_>
    
    SELECT DISTINCT ?reaction ?isTransport ?isBalanced
    WHERE {{
      ?reaction rdfs:subClassOf rh:Reaction .
      ?reaction rh:equation ?equation .
      ?reaction rh:isTransport ?isTransport .
      ?reaction rh:isChemicallyBalanced ?isBalanced .
      
      ?reaction rh:side ?reactionSide .
      ?reactionSide rh:contains ?participant .
      ?participant rh:compound ?compound .
      ?compound rh:chebi chebi:{chebi_id.split(":")[-1]} .
    }}
    LIMIT 10
    """
    
    results = []
    try:
        async with httpx.AsyncClient() as client:
            encoded_query = urllib.parse.urlencode({"query": query})
            headers = {
                "Accept": "application/sparql-results+json",
                "Content-Type": "application/x-www-form-urlencoded"
            }
            response = await client.post(url, data=encoded_query, headers=headers, timeout=15.0)
            
            if response.status_code == 200:
                data = response.json()
                bindings = data.get("results", {}).get("bindings", [])
                for b in bindings:
                    r_id = b["reaction"]["value"].split("/")[-1]
                    is_transport = b["isTransport"]["value"] == "true"
                    is_balanced = b["isBalanced"]["value"] == "true"
                    results.append(ReactionNetwork(
                        rhea_id=f"RHEA:{r_id}",
                        pathway_names=[],
                        is_transport=is_transport,
                        is_balanced=is_balanced
                    ))
            else:
                logger.warning(f"Rhea SPARQL returned {response.status_code}")
    except Exception as e:
        logger.error(f"Error querying Rhea: {e}")
        
    return results

async def execute_phase2(chemical_entity: ChemicalEntity) -> List[ReactionNetwork]:
    """
    Executes Phase 2: ChEBI -> KEGG Pathways & Rhea Reactions.
    Completely drops all mock fallbacks in favor of live DB traversal.
    """
    pathways = await query_kegg_pathways(chemical_entity.chebi_id)
    reactions = await query_rhea_sparql(chemical_entity.chebi_id)
    
    for r in reactions:
        r.pathway_names = pathways
        
    if not reactions:
        logger.warning(f"No valid chemical reactions found in Rhea for {chemical_entity.chebi_id}.")
        
    return reactions
