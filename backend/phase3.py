import httpx
import logging
from typing import List
from models import ChemicalEntity, ReactionNetwork, ProteinCandidate
import asyncio

logger = logging.getLogger(__name__)

async def query_uniprot_category(query_str: str, category_name: str) -> List[ProteinCandidate]:
    """
    Query UniProt REST API to find proteins based on a query string.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    
    params = {
        "query": query_str,
        "format": "json",
        "size": 5, # limit to 5 per category to prevent massive downstream processing
        "fields": "accession,sequence,go_id"
    }
    
    results = []
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(url, params=params, timeout=15.0)
            if response.status_code == 200:
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
                                    if ":" in raw_name:
                                        go_name = raw_name.split(":", 1)[1]
                                    else:
                                        go_name = raw_name
                                    break
                            go_terms.append({"id": go_id, "name": go_name})
                            
                    results.append(ProteinCandidate(
                        uniprot_accession=accession,
                        sequence=sequence,
                        go_terms=go_terms,
                        category=category_name
                    ))
            else:
                logger.warning(f"UniProt API returned {response.status_code} for {query_str}")
    except Exception as e:
        logger.error(f"Error querying UniProt for {category_name}: {e}")
        
    return results

async def execute_phase3(chemical_entity: ChemicalEntity, reactions: List[ReactionNetwork]) -> List[ProteinCandidate]:
    """
    Executes Phase 3: Reaction -> Protein Pool
    Queries UniProt for Enzymes, Transporters, and Receptors completely dynamically.
    """
    tasks = []
    
    # 1. Enzymes: We use the RHEA IDs found in Phase 2 if available.
    # Otherwise, fallback to a general ChEBI search for enzymes.
    rhea_queries = []
    for r in reactions:
        # e.g., RHEA:10604 -> rhea:10604
        rhea_queries.append(f'rhea:"{r.rhea_id.lower()}"')
        
    if rhea_queries:
        rhea_combined = " OR ".join(rhea_queries)
        enzyme_query = f"({rhea_combined}) AND reviewed:true"
    else:
        # Fallback to ChEBI + Enzyme keyword
        enzyme_query = f'chebi:"{chemical_entity.chebi_id}" AND keyword:KW-0255 AND reviewed:true'
        
    tasks.append(query_uniprot_category(enzyme_query, "Enzyme"))
    
    # 2. Transporters: ChEBI + Transport keyword (KW-0813)
    transport_query = f'chebi:"{chemical_entity.chebi_id}" AND keyword:KW-0813 AND reviewed:true'
    tasks.append(query_uniprot_category(transport_query, "Transporter"))
    
    # 3. Receptors: ChEBI + Receptor keyword (KW-0675)
    receptor_query = f'chebi:"{chemical_entity.chebi_id}" AND keyword:KW-0675 AND reviewed:true'
    tasks.append(query_uniprot_category(receptor_query, "Receptor"))
    
    results = await asyncio.gather(*tasks)
    
    # Flatten list of lists
    proteins = [p for sublist in results for p in sublist]
    
    # Deduplicate by accession
    unique_proteins = []
    seen = set()
    for p in proteins:
        if p.uniprot_accession not in seen:
            seen.add(p.uniprot_accession)
            unique_proteins.append(p)
            
    if not unique_proteins:
        logger.warning(f"No proteins found in UniProt for {chemical_entity.chebi_id}")
        
    return unique_proteins

async def execute_phase3_by_ec(ec_number: str) -> List[ProteinCandidate]:
    """
    Executes Phase 3: Skips Phase 1 & 2 by directly querying UniProt with an EC number.
    """
    ec_query = f'ec:{ec_number} AND reviewed:true'
    proteins = await query_uniprot_category(ec_query, f"EC:{ec_number}")
    
    if not proteins:
        logger.warning(f"No proteins found in UniProt for EC Number {ec_number}")
        
    return proteins
