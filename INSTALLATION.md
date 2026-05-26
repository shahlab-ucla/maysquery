# Maysquery — Installation Guide

Maysquery is a Python (FastAPI) backend with a static-file frontend. The
heavy lifting is delegated to two command-line bioinformatics tools — **HMMER**
(`phmmer`) and **Foldseek** — which do not ship for Windows natively. On
Windows we run those inside a **WSL Ubuntu** sandbox; on macOS and Linux we
install them natively via Homebrew or apt.

This document covers:
- [Windows (with WSL bootstrap for new users)](#windows)
- [macOS](#macos)
- [Linux](#linux)
- [Optional: maize AlphaFold structural-discovery index](#optional-maize-alphafold-structural-discovery-index)
- [Troubleshooting](#troubleshooting)

---

## Prerequisites (all platforms)

- Python **3.10 or newer** on your system PATH (`python --version`)
- Git (for cloning)
- ~3 GB free disk for the base install. The optional structural-discovery
  index needs another ~5 GB of download + ~1–2 GB of persistent index.
- Outbound HTTPS to: `rest.uniprot.org`, `rest.ensembl.org`, `data.gramene.org`,
  `ceumass.eps.uspceu.es`, `pubchem.ncbi.nlm.nih.gov`, `sparql.rhea-db.org`,
  `rest.kegg.jp`, `alphafold.ebi.ac.uk`, `ftp.ebi.ac.uk`, `www.ebi.ac.uk`.

---

## Windows

### Step 0 — install WSL (Windows Subsystem for Linux) if you don't already have it

WSL ships with Windows 10 build 19041+ and all of Windows 11 but is not
enabled by default. To check whether you already have it:

```powershell
wsl --status
```

If you get `WSL is not installed`, open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

That command:
1. Enables the Windows Subsystem for Linux feature
2. Enables Virtual Machine Platform
3. Downloads the latest Linux kernel
4. Sets WSL 2 as the default
5. Installs the Ubuntu distribution

**Restart your computer** when the installer says to. On first boot after the
restart, a terminal window opens asking you to create a UNIX username and
password — pick anything you can remember (you may need it later for
`sudo` inside WSL). Once you see the Ubuntu prompt, type `exit` to return to
Windows. WSL is now ready.

### Step 1 — clone the repo

In regular PowerShell (not the WSL terminal). **Pick a directory you own** —
your user profile, `Documents`, or `C:\src` are all fine. Do **not** clone
into `C:\Windows\System32` even if PowerShell happens to start there; Windows
file-permission and antivirus rules make that path painful.

```powershell
cd $HOME\Documents      # or any directory you own
git clone https://github.com/shahlab-ucla/maysquery.git
cd maysquery
```

### Step 2 — allow PowerShell to run local scripts (one-time)

Windows blocks unsigned local `.ps1` scripts by default. The error looks
like `setup.ps1 cannot be loaded because running scripts is disabled on this
system`. Run this **once per PowerShell window** to allow our scripts:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

`-Scope Process` is the safest choice — the change only affects the current
PowerShell session, no persistent system setting. Confirm with **Y** when
prompted. (If you want it persistent for your account, use `-Scope CurrentUser`
instead.)

### Step 3 — run the automated setup

```powershell
.\setup.ps1
```

This walks through six steps:
1. Confirm WSL is available
2. Install `hmmer` + `foldseek` inside WSL Ubuntu
3. Create a native-Windows Python venv in `backend/venv/`
4. `pip install -r backend/requirements.txt`
5. Pre-download the ~12 MB *Zea mays* reference proteome (used by the HMMER
   fallback path)
6. Optionally build the ~1–2 GB maize AlphaFold structural-discovery index
   used by Phase 4.5 (you'll be prompted — answer `N` and do it later if
   you want to start exploring quickly)

If any step fails the script prints the failing command and continues — you
can run individual steps manually (see the [manual fallback](#manual-fallback-windows)).

### Step 4 — launch the app

```powershell
.\run.ps1
```

Your default browser opens to `http://127.0.0.1:8008/static/index.html`.
Leave the PowerShell window running while you use the app; `Ctrl-C` to stop.

> `setup.ps1` and `run.ps1` must both be invoked from the **repo root** (the
> directory containing `backend/`, `setup.ps1`, `run.ps1`). If you see
> `The term '.\run.ps1' is not recognized` you almost certainly `cd`'d into
> `backend/` first — `cd ..` and try again.

### Manual fallback (Windows)

If `setup.ps1` fails, do the steps by hand:

```powershell
# Inside WSL Ubuntu (`wsl` to enter):
sudo apt-get update -y
sudo apt-get install -y hmmer curl tar
curl -sL https://github.com/steineggerlab/foldseek/releases/download/10-941cd33/foldseek-linux-avx2.tar.gz | tar xz -C /tmp
sudo cp /tmp/foldseek/bin/foldseek /usr/local/bin/
phmmer -h    # should print help
foldseek -h  # should print help
exit         # back to PowerShell

# Back in PowerShell, at the repo root:
cd backend
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install --upgrade pip
pip install -r requirements.txt
.\venv\Scripts\uvicorn.exe main:app --reload --host 127.0.0.1 --port 8008
```

---

## macOS

> The Unix install scripts (`setup.sh`, `run.sh`) are documented procedure
> rather than turn-key installers — they're not currently exercised in CI.

### Step 1 — install Homebrew if you don't already have it

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### Step 2 — clone + setup

```bash
git clone https://github.com/shahlab-ucla/maysquery.git
cd maysquery
./setup.sh
```

The script installs Foldseek (`brewsci/bio/foldseek`) and HMMER, creates the
Python venv, installs deps, and prompts for the optional maize AFDB build.

### Step 3 — launch

```bash
./run.sh
```

The script tries to open your default browser; if not, point it to
`http://127.0.0.1:8008/static/index.html`.

### Manual fallback (macOS)

```bash
brew install brewsci/bio/foldseek hmmer
cd backend
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
uvicorn main:app --reload --host 127.0.0.1 --port 8008
```

---

## Linux

### Step 1 — install the bio-binaries

**Debian / Ubuntu:**
```bash
sudo apt-get update
sudo apt-get install -y hmmer curl tar python3-venv
# Foldseek — choose AVX2 if your CPU supports it
curl -sL https://github.com/steineggerlab/foldseek/releases/download/10-941cd33/foldseek-linux-avx2.tar.gz \
  | tar xz -C /tmp
sudo cp /tmp/foldseek/bin/foldseek /usr/local/bin/
```

**Fedora / RHEL:**
```bash
sudo dnf install -y hmmer curl tar python3
# Foldseek install same as above
```

If your CPU is older than ~2013 (no AVX2), use the `sse2` build instead:
`foldseek-linux-sse2.tar.gz`.

### Step 2 — clone + setup

```bash
git clone https://github.com/shahlab-ucla/maysquery.git
cd maysquery
./setup.sh
```

### Step 3 — launch

```bash
./run.sh
```

Then open `http://127.0.0.1:8008/static/index.html`.

---

## Optional: CornCyc curated maize pathway annotation

CornCyc (the Plant Metabolic Network's *Zea mays* Pathway/Genome Database) is
license-restricted but free for non-commercial use. Once enabled, it powers:

- A **curated discovery lane** (4th lane alongside sequence, structure, and
  consensus) — maize genes annotated by PMN as catalysing reactions
  involving the query compound, even when sequence/structure search misses them.
- A dedicated **"CornCyc maize pathway context"** section in the dashboard
  showing every PlantCyc pathway the compound participates in, with the
  list of maize genes per pathway (linked to MaizeGDB).
- New `corncyc_*` columns in the CSV report and a green-tinted block in the
  HTML report.

To enable:

1. Agree to the PMN license and download the CornCyc PGDB tarball from
   <https://plantcyc.org/database_imported/> (latest tested: CornCyc 13.0.0).
2. Extract the archive somewhere on disk. The flatfiles must end up at
   `<dir>/<version>/data/` (the structure the PMN ships in).
3. Either:
   - Drop the extracted folder at the repo root as `corncyc/<version>/data/`
     (default), or
   - Set the `CORNCYC_DIR` environment variable to the directory that
     contains `<version>/data/`.

The pipeline auto-detects the PGDB on first use and loads it in ~0.3 s. A
banner above the pipeline tracker shows the current status (✓ loaded with
gene/reaction/pathway counts, or ○ "not installed" with the expected path).

`backend/data/corncyc/` is excluded from git so you don't accidentally
redistribute the licensed PGDB.

### License attribution

The PMN CornCyc license requires that any derived data carry attribution
to the source database, its authors, and the modifications made. The
Maysquery source includes a complete attribution file at
[`CORNCYC_ATTRIBUTION.txt`](CORNCYC_ATTRIBUTION.txt) that documents:

1. The source database (CornCyc, PMN, version, URLs)
2. The author list and recommended citation
3. Every modification Maysquery applies to CornCyc data (parsing,
   gene-ID case normalisation, ChEBI cross-reference indexing, etc.)

When CornCyc is enabled and a Maysquery query produces CornCyc-derived
output, the HTML report grows a green attribution footer and the CSV
grows a `#`-comment preamble that point readers to that file. The
on-screen dashboard shows attribution under each pathway-context block.

## Optional: maize AlphaFold structural-discovery index

Phase 4.5 (the "hidden ortholog" structural-discovery lane) needs an
indexed copy of the *Zea mays* AlphaFold proteome on disk. It's optional —
the rest of the pipeline runs fine without it, you just lose structure-based
discovery (Phase 4 sequence search and Phase 5 1-to-1 Foldseek alignment
still work).

You can build it at any time from inside `backend/`:

```bash
# Windows
.\venv\Scripts\python.exe install_maize_afdb.py

# macOS / Linux
./venv/bin/python install_maize_afdb.py
```

Or trigger it from the running app — there's a banner with a **Build Index**
button above the pipeline tracker. The build runs in the background and
streams progress to the terminal panel.

Footprint:
- ~5 GB download of `UP000007305_4577_MAIZE_v6.tar` from EBI's AlphaFold FTP
- ~3 GB of decompressed PDB files (removed after indexing if you accept the
  default cleanup)
- ~1–2 GB persistent Foldseek index in `backend/data/foldseek_maize_index/`

Time: ~20–40 minutes total on a fast SSD + 100 Mb/s link. The download is
resumable, so a killed run picks up where it left off.

---

## Troubleshooting

### `setup.ps1 cannot be loaded because running scripts is disabled`

PowerShell blocks unsigned scripts by default. Run this once in the same
PowerShell window before invoking `setup.ps1` or `run.ps1`:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

### `Form data requires "python-multipart" to be installed`

You're running an older clone whose `requirements.txt` predated the fix.
Either `git pull && .\backend\venv\Scripts\pip install -r backend\requirements.txt`,
or install just the missing dep:

```powershell
.\backend\venv\Scripts\python.exe -m pip install python-multipart
```

### AFDB install died with `tarfile.ReadError: unexpected end of data`

A previous download got interrupted and left a truncated tar on disk. The
current installer detects this and deletes the broken tar automatically — so
just run `install_maize_afdb.py` (or click **Build Index** in the UI) again.
The downloader now resumes from where it left off and retries on transient
network errors up to 5 times.

### "Port 8008 already in use"

Edit `run.ps1` / `run.sh` and change the port. Cache-busted JS/CSS in
`backend/static/index.html` are version-pinned but the port is not.

### "Foldseek binary not found" warnings in the log

The pipeline still runs but falls back to a mocked TM-score. Verify with:
```bash
# Windows
wsl foldseek version

# macOS / Linux
foldseek version
```

### Expression-Atlas / Compara look empty

Some Phase 4–7 features depend on third-party APIs (UniProt, Ensembl
Compara, Gramene, EBI Expression Atlas, AlphaFold). When those go down or
get reorganised, the pipeline logs the failure and degrades gracefully —
check the terminal panel in the UI for `[ERROR]` / `[WARNING]` entries.

### PLAZA API returns 403

Known: the anonymous PLAZA Monocots 5.0 endpoint started rejecting requests
in 2026. The pipeline logs this once per session and continues — Phase 4
still has Ensembl Compara + pan-homology and the structural-discovery lane
covers the same biological ground.

### "Cannot find 'backend' directory"

Both `setup.ps1` / `setup.sh` and `run.ps1` / `run.sh` must be run from
the repository root (the directory containing `backend/`).
