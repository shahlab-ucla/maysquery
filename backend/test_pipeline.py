import asyncio
from models import MetaboliteInput, ChemicalEntity
from main import run_single_pipeline
from phase2 import execute_phase2
from phase3 import execute_phase3
from phase4 import execute_phase4
from phase5 import execute_phase5

async def test_chemical(name: str, chebi_id: str, pubchem_cid: str, mass: float):
    print(f"--- Testing {name} (ChEBI: {chebi_id}) ---")
    
    # Simulating Phase 1 output directly to bypass the failing CMM API
    entity = ChemicalEntity(
        name=name,
        monoisotopic_mass=mass,
        chebi_id=chebi_id,
        pubchem_cid=pubchem_cid
    )
    
    try:
        reactions = await execute_phase2(entity)
        print(f"Phase 2: Found {len(reactions)} Reactions")
        
        proteins = await execute_phase3(entity, reactions)
        print(f"Phase 3: Found {len(proteins)} Proteins")
        
        orthologs = await execute_phase4(proteins)
        print(f"Phase 4: Found {len(orthologs)} Orthologs")
        
        targets = await execute_phase5(orthologs)
        print(f"Phase 5: Found {len(targets)} Final Targets")
        
        if targets:
            print(f"Top Target: {targets[0].maize_gene_model} (pLDDT: {targets[0].plddt})")
    except Exception as e:
        print(f"Pipeline failed for {name}: {e}")
    print()

async def main():
    # Bypass CMM (which is currently returning 500) and start from Chemical Entity
    await test_chemical("Fumarate", "CHEBI:18012", "1045", 116.011)
    await test_chemical("Citrate", "CHEBI:16947", "311", 192.027)

if __name__ == "__main__":
    asyncio.run(main())
