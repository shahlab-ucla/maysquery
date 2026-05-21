# Maysquery

Maysquery is an automated, high-throughput pipeline designed to map untargeted metabolomics data (m/z values) to candidate enzyme-coding genes in *Zea mays* (maize). By bridging the gap between raw mass spectrometry features and plant genetics, Maysquery accelerates the discovery of metabolic pathways and functional genomics in maize.

## Features
- **Flexible Entry Points:** Initiate queries using m/z values, known Chemical Names, or Enzyme Commission (EC) numbers.
- **Batch Processing:** Upload CSV files to process hundreds of features automatically.
- **Multi-Database Integration:** Automatically queries CEU Mass Mediator, Metabolomics Workbench, PubChem, ChEBI, Rhea, UniProt, Ensembl Plants, PLAZA, AlphaFold, Expression Atlas, and InterPro.
- **Structural Validation:** Utilizes Foldseek for precise, domain-specific 3D structural homology validation to eliminate false-positive orthologs.
- **Interactive Visualization:** Generates a dynamic Cytoscape Knowledge Graph showing the relationships from m/z ➔ Compound ➔ Reaction ➔ Protein ➔ Maize Gene.

## Installation

### Prerequisites
- Python 3.9+
- Foldseek binary installed and available on your system PATH (`foldseek`).

### Installing Foldseek
Foldseek is strictly required for the 3D structural homology validations performed in Phases 5 and 6. 

**Option A: Conda/Mamba (Recommended)**
```bash
conda install -c conda-forge -c bioconda foldseek
```

**Option B: Precompiled Binaries**
```bash
# Download and extract the latest release for your platform from GitHub
# Example for Linux:
wget https://mmseqs.com/foldseek/foldseek-linux-sse2.tar.gz
tar xvzf foldseek-linux-sse2.tar.gz
export PATH=$(pwd)/foldseek/bin/:$PATH
```
*Note: The backend includes an `install_foldseek.py` utility that the server attempts to run on startup if the binary is missing, but manual installation is recommended.*

### Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/shahlab-ucla/maysquery.git
   cd maysquery/backend
   ```
2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the FastAPI server:
   ```bash
   uvicorn main:app --reload
   ```
4. Access the UI: Open `http://localhost:8000` in your web browser.

## Documentation
For detailed information on the pipeline's logic, database queries, and fallback strategies, please refer to the [USER_GUIDE.md](USER_GUIDE.md).
