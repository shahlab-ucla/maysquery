from pydantic import BaseModel, Field
from typing import List, Optional

# Phase 1: Chemical Entity Standardization
class MetaboliteInput(BaseModel):
    query_type: str = Field("mz", description="Search type: mz, chemical, or ec")
    mz: Optional[float] = None
    mode: Optional[str] = Field("negative", description="Ion mode: positive or negative")
    adducts: Optional[List[str]] = []
    tolerance_ppm: Optional[float] = 5.0
    chemical_name: Optional[str] = None
    ec_number: Optional[str] = None
    spatial_fdr: Optional[float] = None

class ChemicalEntity(BaseModel):
    chebi_id: str
    pubchem_cid: Optional[str] = None
    smiles: Optional[str] = None
    monoisotopic_mass: float

# Phase 2: Reaction Networks
class ReactionNetwork(BaseModel):
    rhea_id: str
    pathway_names: List[str]
    is_transport: bool
    is_balanced: bool

class GOTerm(BaseModel):
    id: str
    name: str

# Phase 3: Protein Pool
class ProteinCandidate(BaseModel):
    uniprot_accession: str
    sequence: str
    go_terms: List[GOTerm]
    category: str = Field(..., description="Enzyme, Transporter, or Receptor")

# Phase 4: Zea mays Projection
class OrthologMapping(BaseModel):
    query_uniprot_id: str = ""
    maize_gene_model: str
    plaza_orthogroup: str
    similarity_score: float
    sources: List[str] = []
    consensus_score: int = 1

# Phase 5: Structural & Transcriptomic Validation
class ValidatedTarget(BaseModel):
    maize_gene_model: str
    tm_score: float
    plddt: float
    tissue_expression_fpkm: dict

# Phase 6: Domain Structural Homology
class DomainValidatedTarget(BaseModel):
    maize_gene_model: str
    query_uniprot_id: str
    pfam_domain_id: str
    pfam_domain_name: str
    domain_start: int
    domain_end: int
    domain_tm_score: float
