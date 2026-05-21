import asyncio
from models import OrthologMapping
from phase5 import execute_phase5
import os

async def main():
    # Pass an ortholog mapping with a known query UniProt ID (Fumarase: P07954)
    # Target maize model Zm00001eb016240 usually maps back to a UniProt in our fallback dict
    mappings = [
        OrthologMapping(
            query_uniprot_id="P07954", 
            maize_gene_model="Zm00001eb016240", 
            plaza_orthogroup="ORTHO_001",
            similarity_score=87.0,
            sources=["Ensembl", "PLAZA"],
            consensus_score=2
        )
    ]
    
    print("Executing Phase 5 structural validation...")
    results = await execute_phase5(mappings)
    
    print(f"Validated {len(results)} targets.")
    for res in results:
        print(f"Gene: {res.maize_gene_model:<18} | pLDDT: {res.plddt:>5} | TM-Score: {res.tm_score} | Expr: {res.tissue_expression_fpkm}")

if __name__ == "__main__":
    asyncio.run(main())
