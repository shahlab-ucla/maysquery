import httpx
import logging
from typing import List, Optional
from models import MetaboliteInput, ChemicalEntity

logger = logging.getLogger(__name__)

async def cmm_search(input_data: MetaboliteInput) -> List[dict]:
    """
    Query CEU Mass Mediator API to find putative compounds for a given m/z.
    """
    url = "https://ceumass.eps.uspceu.es/api/v3/batch"
    
    adducts = input_data.adducts
    if not adducts:
        if input_data.mode == "positive":
            adducts = ["M+H", "M+Na"]
        else:
            adducts = ["M-H", "M+Cl"]
            
    payload = {
        "masses": [input_data.mz],
        "tolerance": input_data.tolerance_ppm,
        "tolerance_mode": "ppm",
        "ion_mode": input_data.mode,
        "adducts": adducts,
        "databases": ["HMDB", "KEGG"]
    }
    
    try:
        async with httpx.AsyncClient() as client:
            response = await client.post(url, json=payload, timeout=10.0)
            if response.status_code == 200:
                data = response.json()
                return data.get("results", [])
            else:
                logger.warning(f"CMM API returned {response.status_code}")
                return []
    except Exception as e:
        logger.error(f"Error querying CMM: {e}")
        return []

async def fallback_metabolomics_workbench(mz: float, adducts: List[str], tolerance: float) -> str:
    """
    Fallback to Metabolomics Workbench if CMM fails.
    """
    adduct = adducts[0] if adducts else "M-H"
    safe_adduct = adduct.replace("+", "%2B")
    url = f"https://www.metabolomicsworkbench.org/rest/moverz/MB/{mz}/{safe_adduct}/{tolerance}/json"
    
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(url, timeout=10.0)
            if response.status_code == 200:
                data = response.json()
                if isinstance(data, list) and len(data) > 0:
                    return data[0].get("name", "")
                elif isinstance(data, dict) and len(data) > 0:
                    first_key = list(data.keys())[0]
                    return data[first_key].get("name", "")
    except Exception as e:
        logger.error(f"MW fallback failed: {e}")
    return ""

async def fetch_pubchem_data(compound_name: str) -> Optional[dict]:
    """
    Query PubChem PUG-REST for CID, SMILES, and Monoisotopic Mass.
    """
    # Replace spaces for URL
    safe_name = compound_name.replace(" ", "%20")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{safe_name}/property/IsomericSMILES,MonoisotopicMass/JSON"
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(url, timeout=10.0)
            if response.status_code == 200:
                data = response.json()
                props = data.get("PropertyTable", {}).get("Properties", [])
                if props:
                    return props[0]
            return None
    except Exception as e:
        logger.error(f"Error querying PubChem: {e}")
        return None

async def map_pubchem_to_chebi(cid: str) -> Optional[str]:
    """
    Map PubChem CID to ChEBI ID using PubChem RDF or REST cross-references.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(url, timeout=10.0)
            if response.status_code == 200:
                data = response.json()
                xrefs = data.get("InformationList", {}).get("Information", [])
                for info in xrefs:
                    registry_ids = info.get("RegistryID", [])
                    for reg in registry_ids:
                        if "CHEBI:" in reg:
                            return reg
            return None
    except Exception as e:
        logger.error(f"Error querying ChEBI mapping: {e}")
        return None

async def execute_phase1(input_data: MetaboliteInput) -> ChemicalEntity:
    """
    Executes Phase 1 of the pipeline: mz -> CMM -> PubChem -> ChEBI.
    Now uses STRICT API integration with MW fallback.
    """
    name = ""
    try:
        cmm_results = await cmm_search(input_data)
        if cmm_results and len(cmm_results) > 0:
            compounds = cmm_results[0].get("compounds", [])
            if compounds:
                name = compounds[0].get("name", "")
    except Exception as e:
        logger.error(f"CMM error: {e}")
            
    if not name:
        logger.warning("CMM lookup failed. Triggering Metabolomics Workbench fallback...")
        name = await fallback_metabolomics_workbench(input_data.mz, input_data.adducts, input_data.tolerance_ppm)
        
    if not name:
        logger.error("Could not resolve m/z to a chemical name.")
        raise ValueError(f"CMM and MW lookups failed for m/z {input_data.mz}")
            
    return await fetch_chemical_entity_by_name(name)

async def fetch_chemical_entity_by_name(name: str) -> ChemicalEntity:
    pubchem_data = await fetch_pubchem_data(name)
    if not pubchem_data:
        logger.error(f"Could not resolve compound name '{name}' in PubChem.")
        raise ValueError(f"PubChem lookup failed for '{name}'")
        
    cid = str(pubchem_data["CID"])
    smiles = pubchem_data.get("IsomericSMILES", "")
    mass = pubchem_data.get("MonoisotopicMass", 0.0)
    
    chebi_id = await map_pubchem_to_chebi(cid)
    if not chebi_id:
        logger.warning(f"No CHEBI mapping found for PubChem CID {cid}. Using placeholder.")
        chebi_id = f"PUBCHEM:{cid}"
        
    return ChemicalEntity(
        chebi_id=chebi_id,
        pubchem_cid=cid,
        smiles=smiles,
        monoisotopic_mass=mass
    )
