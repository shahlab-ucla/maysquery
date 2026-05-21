import asyncio
from models import ProteinCandidate
from phase4 import execute_phase4

async def main():
    # Provide the proteins that match our fallback mocks to ensure hits
    proteins = [
        ProteinCandidate(uniprot_accession="P07954", sequence="MOCK", go_terms=[], category="Enzyme"),
        ProteinCandidate(uniprot_accession="Q13183", sequence="MOCK", go_terms=[], category="Enzyme"),
        ProteinCandidate(uniprot_accession="P04424", sequence="MOCK", go_terms=[], category="Enzyme"),
        ProteinCandidate(uniprot_accession="Q768R5", sequence="MOCK", go_terms=[], category="Transporter")
    ]
    
    print("Executing Phase 4 with multi-DB queries...")
    results = await execute_phase4(proteins)
    
    print(f"Found {len(results)} merged ortholog mappings.")
    for res in results:
        print(f"Gene: {res.maize_gene_model:<18} | Consensus: {res.consensus_score} | Similarity: {res.similarity_score:>5} | Sources: {res.sources}")

if __name__ == "__main__":
    asyncio.run(main())
