# Maysquery User Guide

Maysquery automates the complex bioinformatics workflow required to trace a putative metabolite feature back to a specific candidate gene in *Zea mays*. This guide details the 6-phase pipeline, the data sources queried, and the fallback strategies employed to ensure robustness.

## 1. Input Modes

Maysquery supports three primary entry points, allowing flexibility depending on what data you already have:

*   **m/z Value:** The standard starting point for untargeted metabolomics. Requires an m/z float, ionization mode (positive/negative), and adducts. Triggers the full 6-phase pipeline.
*   **Chemical Name:** Bypasses Phase 1 (m/z resolution) and starts directly by looking up the chemical structure in PubChem.
*   **EC Number:** Bypasses Phases 1 and 2. Directly searches UniProt for proteins annotated with the given Enzyme Commission number (e.g., 4.2.1.2) to find maize orthologs.

## 2. Pipeline Phases & Database Queries

### Phase 1: Chemical Entity Standardization
**Goal:** Resolve an m/z value into a specific chemical entity with a defined mass, structure, and database identifiers.
*   **Primary Database:** CEU Mass Mediator (CMM) API. Used to predict the most likely metabolite name based on m/z, adducts, and tolerance.
*   **Fallback Strategy:** If CMM is down (e.g., returns a 500 error), the pipeline automatically falls back to the **Metabolomics Workbench (MW) REST API** (`moverz` endpoint) to resolve the feature.
*   **Standardization:** Once a name is resolved, the pipeline queries **PubChem (PUG-REST)** to retrieve the CID, SMILES string, and Monoisotopic Mass. It then cross-references PubChem to find the **ChEBI ID**.

### Phase 2: Reaction Networks
**Goal:** Identify the biochemical reactions in which the chemical entity participates.
*   **Primary Database:** **Rhea SPARQL Endpoint**.
*   **Query Logic:** The pipeline uses the resolved ChEBI ID to query Rhea for all metabolic and transport reactions involving that compound. It flags whether the reaction is a transport event and retrieves associated pathway names.

### Phase 3: Protein Pool (Pan-life)
**Goal:** Gather a list of proteins (from any organism) known to catalyze the identified reactions.
*   **Primary Database:** **UniProt REST API**.
*   **Query Logic:**
    *   **Enzymes:** Searches UniProt using the specific Rhea IDs found in Phase 2 (`rhea:"RHEA:XXXX"`).
    *   **Transporters & Receptors:** Falls back to searching UniProt using the ChEBI ID combined with specific UniProt keywords (e.g., `KW-0813` for transport, `KW-0675` for receptors).
    *   **Fallback:** If Phase 1/2 are bypassed via "EC Number" input, UniProt is queried directly using `ec:X.X.X.X`.

### Phase 4: Ortholog Projection
**Goal:** Map the pan-life protein pool to specific *Zea mays* orthologs.
*   **Primary Databases:** **Ensembl Plants REST API** and **PLAZA API**.
*   **Query Logic:** The pipeline queries Ensembl Plants to find homologous genes in maize. It independently queries PLAZA to find orthogroups containing the UniProt ID, filtering for maize members. A consensus score is generated if both databases agree on the ortholog mapping.

### Phase 5: Structural & Transcriptomic Validation
**Goal:** Filter the candidate maize genes using 3D structural homology and evaluate their tissue expression.
*   **Structural Validation:**
    *   Queries the **AlphaFold DB API** to retrieve the PDB structures for both the query UniProt protein and the candidate maize protein.
    *   Evaluates the `pLDDT` (confidence) score.
    *   **Fallback:** If the AlphaFold PDBs are successfully downloaded, it runs a local `foldseek easy-search` to calculate a whole-protein structural TM-Score.
*   **Transcriptomic Validation:**
    *   Queries the **EBI Expression Atlas (RNA-seq)** for baseline tissue expression (FPKM) of the maize gene.

### Phase 6: Domain Structural Homology
**Goal:** Perform highly targeted structural alignment on catalytic domains rather than whole proteins, preventing false rejections caused by domain shuffling or regulatory sequence variations.
*   **Primary Database:** **InterPro REST API**.
*   **Query Logic:**
    1.  Queries InterPro to identify the longest **Pfam** domain on the query protein.
    2.  Extracts the exact start and end residue coordinates of that domain.
    3.  Dynamically slices the downloaded AlphaFold PDB to isolate just the domain structure.
    4.  Runs a localized Foldseek alignment between the *sliced* domain PDB and the *full* maize target PDB.
    5.  Outputs a highly accurate **Domain TM-Score**.

## 3. Batch Processing
For high-throughput analysis, use the **Batch Processing** tab. You can download the CSV template, fill it with your m/z values, modes, and parameters, and upload it. The pipeline will process queries sequentially and generate downloadable HTML and CSV reports summarizing the Top Maize Targets for each feature.
