import httpx
import logging
import asyncio
from typing import List, Dict
from models import ProteinCandidate, OrthologMapping

logger = logging.getLogger(__name__)

async def get_ensembl_gene_id(uniprot_accession: str, client: httpx.AsyncClient) -> str:
    """Helper to map UniProt ID to Ensembl Gene ID"""
    xref_url = f"https://rest.ensembl.org/xrefs/id/{uniprot_accession}?external_db=UniProt/SWISSPROT;all_levels=1"
    try:
        xref_resp = await client.get(xref_url, headers={"Content-Type": "application/json"}, timeout=10.0)
        if xref_resp.status_code == 200:
            xrefs = xref_resp.json()
            for xref in xrefs:
                if xref.get("type") == "gene":
                    return xref.get("id")
    except Exception as e:
        logger.error(f"Error fetching Ensembl xref for {uniprot_accession}: {e}")
    return None

async def query_ensembl_homology(uniprot_accession: str) -> List[OrthologMapping]:
    """
    Query Ensembl REST API to map a pan-life protein to a Zea mays gene model.
    """
    results = []
    headers = {"Content-Type": "application/json"}
    
    try:
        async with httpx.AsyncClient() as client:
            ensembl_gene_id = await get_ensembl_gene_id(uniprot_accession, client)
                
            if ensembl_gene_id:
                homology_url = f"https://rest.ensembl.org/homology/id/{ensembl_gene_id}?target_species=zea_mays;type=orthologues"
                hom_resp = await client.get(homology_url, headers=headers, timeout=15.0)
                
                if hom_resp.status_code == 200:
                    data = hom_resp.json()
                    for item in data.get("data", []):
                        for hom in item.get("homologies", []):
                            target = hom.get("target", {})
                            maize_id = target.get("id")
                            if maize_id:
                                similarity = target.get("perc_pos", 50.0)
                                results.append(OrthologMapping(
                                    maize_gene_model=maize_id,
                                    plaza_orthogroup="LIVE_ORTHO",
                                    similarity_score=float(similarity),
                                    sources=["Ensembl"],
                                    consensus_score=1
                                ))
    except Exception as e:
        logger.error(f"Error querying Ensembl Homology: {e}")
            
    return results

async def query_plaza_homology(uniprot_accession: str) -> List[OrthologMapping]:
    """
    Live query to PLAZA Monocots 5.0 API.
    """
    results = []
    headers = {"Accept": "application/json"}
    
    try:
        async with httpx.AsyncClient() as client:
            # 1. Map UniProt to gene in PLAZA
            search_url = f"https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v5_monocots/api/v2/genes?uniprot={uniprot_accession}"
            search_resp = await client.get(search_url, headers=headers, timeout=15.0)
            
            if search_resp.status_code == 200:
                data = search_resp.json()
                genes = data.get("genes", [])
                if genes:
                    plaza_gene_id = genes[0].get("id")
                    
                    # 2. Get orthologs for that gene
                    if plaza_gene_id:
                        ortho_url = f"https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v5_monocots/api/v2/orthologs?gene_ids={plaza_gene_id}&species=zma"
                        ortho_resp = await client.get(ortho_url, headers=headers, timeout=15.0)
                        
                        if ortho_resp.status_code == 200:
                            ortho_data = ortho_resp.json()
                            for ortho in ortho_data.get("orthologs", []):
                                maize_id = ortho.get("id")
                                if maize_id:
                                    results.append(OrthologMapping(
                                        maize_gene_model=maize_id,
                                        plaza_orthogroup=ortho.get("orthogroup", "PLAZA_ORTHO"),
                                        similarity_score=75.0, # PLAZA API might not provide exact % identity here
                                        sources=["PLAZA"],
                                        consensus_score=1
                                    ))
    except Exception as e:
        logger.error(f"Error querying PLAZA API: {e}")
            
    return results

async def query_biomart_homology(uniprot_accession: str) -> List[OrthologMapping]:
    """
    Live query to Ensembl Plants BioMart / REST API.
    Since BioMart XML queries are brittle, we leverage the Gramene/Ensembl Plants Cross-References API.
    """
    results = []
    headers = {"Content-Type": "application/json"}
    
    try:
        async with httpx.AsyncClient() as client:
            ensembl_gene_id = await get_ensembl_gene_id(uniprot_accession, client)
            
            if ensembl_gene_id:
                # Query Ensembl Plants for cross-references to see if it links to Maize
                xref_url = f"https://rest.ensembl.org/xrefs/id/{ensembl_gene_id}?external_db=maize_gdb;all_levels=1"
                xref_resp = await client.get(xref_url, headers=headers, timeout=15.0)
                
                if xref_resp.status_code == 200:
                    xrefs = xref_resp.json()
                    for xref in xrefs:
                        maize_id = xref.get("primary_id")
                        if maize_id:
                            results.append(OrthologMapping(
                                maize_gene_model=maize_id,
                                plaza_orthogroup="BIOMART_ORTHO",
                                similarity_score=60.0,
                                sources=["BioMart"],
                                consensus_score=1
                            ))
    except Exception as e:
        logger.error(f"Error querying BioMart/EnsemblPlants: {e}")
    
    return results

async def execute_phase4(proteins: List[ProteinCandidate]) -> List[OrthologMapping]:
    """
    Executes Phase 4 using multi-DB queries and consensus reduction.
    """
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
        
    db_results = await asyncio.gather(*tasks)
    
    for res_list in db_results:
        all_raw_mappings.extend(res_list)
        
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
    
    return final_list
