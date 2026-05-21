import httpx
import logging
from typing import List, Optional
from models import MetaboliteInput, ChemicalEntity, ExecutionLogEntry
from urllib.parse import quote

logger = logging.getLogger(__name__)

async def cmm_search(input_data: MetaboliteInput) -> List[dict]:
    """
    Query CEU Mass Mediator API to find putative compounds for a given m/z.

    The v3 /batch endpoint requires `metabolites_type`, `masses_mode`, and
    `deuterium` in addition to the obvious m/z fields — omitting any of them
    causes the Tomcat backend to throw a NullPointerException (HTTP 500).
    Each entry in `results` is the compound itself (no nested array).
    """
    url = "https://ceumass.eps.uspceu.es/api/v3/batch"

    # CMM rejects (with HTTP 500 + NullPointerException) ANY adduct it doesn't
    # recognize, taking down the whole batch. The set below is empirically
    # verified against the v3 endpoint as of 2026 — be careful adding more.
    adducts = input_data.adducts
    if not adducts:
        if input_data.mode == "positive":
            adducts = ["M+H", "M+Na", "M+K", "M+NH4", "M+H-H2O"]
        else:
            adducts = ["M-H", "M+Cl", "M+Br"]

    payload = {
        "metabolites_type": "all-except-peptides",
        "databases": ["all-except-mine"],
        "masses_mode": "mz",
        "ion_mode": input_data.mode,
        "adducts": adducts,
        "deuterium": False,
        "tolerance": float(input_data.tolerance_ppm),
        "tolerance_mode": "ppm",
        "masses": [float(input_data.mz)],
    }

    try:
        async with httpx.AsyncClient() as client:
            response = await client.post(url, json=payload, timeout=20.0)
            if response.status_code == 200:
                data = response.json()
                return data.get("results", [])
            logger.warning(
                f"CMM API returned HTTP {response.status_code}. "
                f"Body[:200]={response.text[:200]!r}"
            )
            return []
    except Exception as e:
        logger.error(f"Error querying CMM ({type(e).__name__}): {e}")
        return []


def _parse_mw_tsv(body: str) -> List[dict]:
    """
    Parse Metabolomics Workbench moverz TSV output into a list of dicts
    sorted by abs(delta) ascending. The endpoint returns text/html but the
    body is actually tab-separated values with one header row.
    """
    rows = []
    lines = [ln for ln in body.splitlines() if ln.strip()]
    if not lines:
        return rows
    header = [h.strip() for h in lines[0].split("\t")]
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) < len(header):
            continue
        row = dict(zip(header, parts))
        try:
            row["_delta"] = abs(float(row.get("Delta", "9.99")))
        except ValueError:
            row["_delta"] = 9.99
        rows.append(row)
    rows.sort(key=lambda r: r["_delta"])
    return rows


async def fallback_metabolomics_workbench(mz: float, adducts: List[str], tolerance: float, mode: str = "negative") -> str:
    """
    Fallback to Metabolomics Workbench /rest/moverz endpoint.

    Quirks of the MW endpoint:
      * It 302-redirects from /rest/moverz/.../json to the underlying PHP
        script — httpx must be told to follow_redirects.
      * Despite the `/json` URL suffix, it returns tab-separated text
        with a header row; r.json() always fails.
      * Tolerance must be a whole number (ppm).
    """
    if not adducts:
        adducts = ["M+H", "M+Na"] if mode == "positive" else ["M-H", "M+Cl"]

    # MW expects integer ppm tolerance — round up so we don't lose hits
    tol = max(1, int(round(tolerance)))

    async with httpx.AsyncClient(follow_redirects=True) as client:
        for adduct in adducts:
            safe_adduct = adduct.replace("+", "%2B")
            url = f"https://www.metabolomicsworkbench.org/rest/moverz/MB/{mz}/{safe_adduct}/{tol}/json"
            try:
                resp = await client.get(url, timeout=20.0)
                if resp.status_code != 200:
                    logger.warning(f"MW {adduct} returned HTTP {resp.status_code}")
                    continue
                rows = _parse_mw_tsv(resp.text)
                if rows and rows[0].get("Name"):
                    logger.info(f"MW {adduct}: top hit '{rows[0]['Name']}' (Δ={rows[0].get('Delta','?')})")
                    return rows[0]["Name"]
            except Exception as e:
                logger.error(f"MW {adduct} failed ({type(e).__name__}): {e}")
    return ""

async def fetch_pubchem_data(compound_name: str, logs: List[ExecutionLogEntry] = None) -> Optional[dict]:
    """
    Query PubChem PUG-REST for CID, SMILES, and Monoisotopic Mass.
    """
    if logs is None:
        logs = []
    
    encoded_name = quote(compound_name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/MonoisotopicMass,IsomericSMILES/JSON"
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

async def execute_phase1(input_data: MetaboliteInput, logs: List[ExecutionLogEntry] = None) -> ChemicalEntity:
    """
    Executes Phase 1 of the pipeline: mz -> CMM -> PubChem -> ChEBI.
    Now uses STRICT API integration with MW fallback.
    """
    if logs is None:
        logs = []

    logs.append(ExecutionLogEntry(
        phase=1, database="CEU Mass Mediator", status="info", hits=0,
        message=f"Querying CMM for m/z {input_data.mz} ({input_data.mode}, ±{input_data.tolerance_ppm}ppm)"
    ))

    name = ""
    try:
        cmm_results = await cmm_search(input_data)
        # CMM v3 /batch returns a flat list — each entry IS the matched compound.
        # Many candidates can share the same isomeric formula (e.g. at m/z
        # 115.0037 [M-H]: formylpyruvate, fumaric acid, maleic acid). Rank by
        # finalScore desc, then by smallest |error_ppm|.
        if cmm_results:
            scored = sorted(
                cmm_results,
                key=lambda r: (-(r.get("finalScore") or -999), abs(r.get("error_ppm") or 999)),
            )
            top = scored[0]
            name = top.get("name", "") or ""
            top_names = [r.get("name", "?") for r in scored[:5] if r.get("name")]
            adduct = top.get("adduct", "")
            err_ppm = top.get("error_ppm", "?")
        logs.append(ExecutionLogEntry(
            phase=1, database="CEU Mass Mediator",
            status="success" if name else "warning",
            hits=len(cmm_results) if cmm_results else 0,
            message=(
                f"CMM returned {len(cmm_results)} candidate(s); top: '{name}' "
                f"(adduct={adduct}, Δ={err_ppm}ppm)"
                + (f". Other candidates: {', '.join(top_names[1:])}" if len(top_names) > 1 else "")
                if name else
                f"CMM returned no compounds at m/z {input_data.mz} in {input_data.mode} mode"
            )
        ))
    except Exception as e:
        logger.error(f"CMM error: {e}")
        logs.append(ExecutionLogEntry(
            phase=1, database="CEU Mass Mediator", status="error", hits=0,
            message=f"CMM query failed: {type(e).__name__}: {e}"
        ))

    if not name:
        logs.append(ExecutionLogEntry(
            phase=1, database="Metabolomics Workbench", status="info", hits=0,
            message=f"CMM miss — falling back to Metabolomics Workbench (trying multiple adducts)"
        ))
        name = await fallback_metabolomics_workbench(
            input_data.mz, input_data.adducts, input_data.tolerance_ppm, input_data.mode
        )
        logs.append(ExecutionLogEntry(
            phase=1, database="Metabolomics Workbench",
            status="success" if name else "error",
            hits=1 if name else 0,
            message=f"MW resolved top hit: '{name}'" if name else f"MW returned no matches at m/z {input_data.mz} in {input_data.mode} mode"
        ))

    if not name:
        logger.error("Could not resolve m/z to a chemical name.")
        opposite = "negative" if input_data.mode == "positive" else "positive"
        raise ValueError(
            f"No compound matched m/z {input_data.mz} in {input_data.mode} mode "
            f"(tolerance ±{input_data.tolerance_ppm}ppm). "
            f"CMM and Metabolomics Workbench both returned 0 candidates. "
            f"If you expected a hit, try {opposite} mode or widen the tolerance — "
            f"e.g. 115.0037 is fumaric acid as [M-H]⁻ in negative mode."
        )

    return await fetch_chemical_entity_by_name(name, logs, input_data)

async def fetch_chemical_entity_by_name(name: str, logs: List[ExecutionLogEntry] = None, input_data: MetaboliteInput = None) -> ChemicalEntity:
    if logs is None:
        logs = []

    logs.append(ExecutionLogEntry(
        phase=1, database="PubChem PUG-REST", status="info", hits=0,
        message=f"Resolving compound name '{name}' to PubChem CID"
    ))
    pubchem_data = await fetch_pubchem_data(name, logs)
    if not pubchem_data:
        logs.append(ExecutionLogEntry(
            phase=1, database="PubChem PUG-REST", status="error", hits=0,
            message=f"PubChem lookup failed for '{name}'"
        ))
        logger.error(f"Could not resolve compound name '{name}' in PubChem.")
        raise ValueError(f"PubChem lookup failed for '{name}'")

    cid = str(pubchem_data["CID"])
    smiles = pubchem_data.get("IsomericSMILES", "")
    mass = pubchem_data.get("MonoisotopicMass", 0.0)
    logs.append(ExecutionLogEntry(
        phase=1, database="PubChem PUG-REST", status="success", hits=1,
        message=f"PubChem CID {cid} resolved (mass={mass}, SMILES={smiles[:40]}{'...' if len(smiles) > 40 else ''})"
    ))

    logs.append(ExecutionLogEntry(
        phase=1, database="PubChem→ChEBI xref", status="info", hits=0,
        message=f"Mapping CID {cid} to ChEBI identifier"
    ))
    chebi_id = await map_pubchem_to_chebi(cid)
    if not chebi_id:
        logger.warning(f"No CHEBI mapping found for PubChem CID {cid}. Using placeholder.")
        chebi_id = f"PUBCHEM:{cid}"
        logs.append(ExecutionLogEntry(
            phase=1, database="PubChem→ChEBI xref", status="warning", hits=0,
            message=f"No ChEBI mapping — using placeholder {chebi_id}"
        ))
    else:
        logs.append(ExecutionLogEntry(
            phase=1, database="PubChem→ChEBI xref", status="success", hits=1,
            message=f"Mapped to {chebi_id}"
        ))

    return ChemicalEntity(
        chebi_id=chebi_id,
        pubchem_cid=cid,
        smiles=smiles,
        monoisotopic_mass=mass
    )
