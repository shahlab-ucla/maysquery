from pydantic import BaseModel, Field
from typing import List, Optional, Any
from datetime import datetime

class ExecutionLogEntry(BaseModel):
    phase: int
    database: str
    status: str  # e.g., 'success', 'error', 'info'
    hits: int
    message: str
    timestamp: str = Field(default_factory=lambda: datetime.utcnow().isoformat())


class PipelineConfig(BaseModel):
    """
    All tunable cutoffs/caps/thresholds in one place. Sent with every pipeline
    request from the Configuration tab; phases read from here instead of
    module-level constants. Defaults match the per-module historical defaults
    so behaviour is unchanged if a client omits the config entirely.
    """
    # ----- Phase 1: Chemical identification -----
    cmm_tolerance_ppm: float = Field(5.0, ge=1.0, le=50.0,
        description="m/z mass tolerance (ppm) for CMM and Metabolomics Workbench lookups. Lower = stricter; higher = more candidates.")

    # ----- Phase 2: Reactions -----
    rhea_fetch_limit: int = Field(100, ge=10, le=1000,
        description="Max Rhea reactions to fetch per metabolite. Most metabolites have <50; common cofactors (ATP, NAD+) can have many hundreds.")
    chebi_expansion_depth: int = Field(2, ge=0, le=4,
        description="Hops through the ChEBI conjugate-acid/base graph. 2 covers diprotic acids (fumarate, succinate). 0 disables expansion (use the seed ChEBI only).")

    # ----- Phase 3: Protein pool -----
    uniprot_size_per_category: int = Field(25, ge=1, le=200,
        description="Top-N reviewed UniProt entries per category (Enzyme / Transporter / Receptor). Each protein fans out into 3 ortholog API calls + 1 Foldseek search downstream.")

    # ----- Phase 4: Sequence orthology -----
    hmmer_e_value: float = Field(1e-5, gt=0, le=1.0,
        description="Fallback HMMER (phmmer) E-value cutoff. Used only when live ortholog APIs (Ensembl/PLAZA/BioMart) return zero. Lower = stricter homology.")

    # ----- Phase 4.5: Structural discovery (Foldseek vs maize AFDB) -----
    foldseek_tm_threshold: float = Field(0.5, ge=0.2, le=0.95,
        description="Minimum max(qTM, tTM) to accept a structural ortholog. 0.5 = conventional 'same fold' threshold; 0.4 catches more distant orthologs at the cost of more noise.")
    foldseek_max_hits_per_query: int = Field(10, ge=1, le=100,
        description="Top-N maize structural hits to keep per pan-life query protein.")
    foldseek_concurrency: int = Field(2, ge=1, le=8,
        description="Parallel Foldseek processes. Foldseek itself is multi-threaded; 2 saturates a typical 8-core workstation.")

    # ----- Phase 5: Enrichment -----
    enrichment_top_n: int = Field(10, ge=1, le=100,
        description="Top-N discovered orthologs to enrich with pLDDT, tissue expression, and 1-to-1 Foldseek TM (for sequence-only hits).")
    plddt_threshold: float = Field(70.0, ge=0.0, le=100.0,
        description="Minimum AlphaFold pLDDT to keep an enriched target. <70 = uncertain backbone; <50 = likely disordered.")

    # ----- Phase 7: Pan-plant Compara -----
    compara_max_orthologs_per_target: int = Field(10, ge=1, le=50,
        description="Top-N Ensembl Compara plant orthologs to return per validated maize target.")


# Phase 1: Chemical Entity Standardization
class MetaboliteInput(BaseModel):
    query_type: str = Field("mz", description="Search type: mz, chemical, or ec")
    mz: Optional[float] = None
    mode: str = "negative"
    adducts: List[str] = []
    tolerance_ppm: float = 5.0          # kept for backward compat; superseded by pipeline_config.cmm_tolerance_ppm
    chemical_name: Optional[str] = None
    ec_number: Optional[str] = None
    hmmer_e_value: float = 1e-5         # kept for backward compat; superseded by pipeline_config.hmmer_e_value
    spatial_fdr: Optional[float] = None
    pipeline_config: PipelineConfig = Field(default_factory=PipelineConfig)

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
    equation: Optional[str] = None       # Human-readable equation, e.g. "(S)-malate = fumarate + H2O"
    label: Optional[str] = None          # Reaction label / common name if available
    ec_numbers: List[str] = Field(default_factory=list)
    pathway_ids: List[str] = Field(default_factory=list)   # KEGG pathway IDs (path:map01100, ...)

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
    tissue_expression_fpkm: dict = Field(default_factory=dict)
    # Replaces the dead EBI /gxa/json/search/baseline endpoint. Sourced from
    # Gramene which lists every Expression Atlas experiment the gene was
    # detected in (qualitative "is expressed?" signal across N experiments).
    n_expression_experiments: int = 0
    expression_experiments: List[str] = Field(default_factory=list)
    # Provenance: which Phase 5 sub-steps actually ran for this target.
    enrichment_kind: str = "full"   # "full" | "cheap" | "tm_only"

# Phase 6: Domain Structural Homology
class DomainValidatedTarget(BaseModel):
    maize_gene_model: str
    query_uniprot_id: str
    pfam_domain_id: str
    pfam_domain_name: str
    domain_start: int
    domain_end: int
    domain_tm_score: float

# Phase 7: Advanced Homology Validation (Tier 1)
class EnsemblOrtholog(BaseModel):
    species: str
    gene_id: str
    protein_id: Optional[str] = None
    percent_identity: float

class AdvancedHomologyTarget(BaseModel):
    maize_gene_model: str
    ensembl_orthologs: List[EnsemblOrtholog]
    plaza_synteny_blocks: Optional[List[str]] = None
    plaza_orthogroup: Optional[str] = None
