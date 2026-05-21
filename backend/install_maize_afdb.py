"""
install_maize_afdb.py — One-time setup of a local Foldseek structural database
built from the AlphaFold predicted proteome of *Zea mays* (UP000007305, taxid
4577).

The resulting index powers Phase 4.5 (structure-guided ortholog discovery): for
each pan-life query enzyme, we run `foldseek easy-search query.pdb maize_db`
to surface maize proteins whose sequence has diverged below sequence-search
detection limits but whose 3D fold is conserved.

Disk footprint (as of AlphaFold v6, ~40k maize predictions):
  - Tar download:          ~5 GB (resumable)
  - Extracted PDBs:        ~3 GB (intermediates auto-removed if cleanup=True)
  - Foldseek index:        ~1–2 GB
Time on a fast SSD + 100 Mb/s link: ~20–40 minutes total (download dominates).

CLI usage (interactive):
    python install_maize_afdb.py

Programmatic (used by /api/maize_afdb/install):
    from install_maize_afdb import ensure_db_ready
    await ensure_db_ready(interactive=False, progress_callback=cb)
"""

from __future__ import annotations

import asyncio
import gzip
import logging
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from typing import Awaitable, Callable, Optional

import httpx

logger = logging.getLogger(__name__)

# ----- Paths -----

BACKEND_DIR = Path(__file__).resolve().parent
DATA_DIR = BACKEND_DIR / "data"
PDB_DIR = DATA_DIR / "maize_alphafold_pdbs"
INDEX_DIR = DATA_DIR / "foldseek_maize_index"
TAR_PATH = DATA_DIR / "maize_alphafold_proteome.tar"
DB_PREFIX = INDEX_DIR / "maize_db"
DB_MARKER = INDEX_DIR / "maize_db.dbtype"  # foldseek createdb output marker
PDB_DONE_MARKER = PDB_DIR / ".extraction_complete"

# Override via env var so e.g. CI can point to a smaller test tar.
MAIZE_PROTEOME_URL = os.environ.get(
    "MAIZE_AFDB_URL",
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000007305_4577_MAIZE_v6.tar",
)

# Progress callback signature: cb(message: str, stage: str | None, pct: float | None)
ProgressCB = Callable[..., Optional[Awaitable[None]]]


# ----- Path helpers -----

def to_wsl_path(win_path: str | Path) -> str:
    """Convert a Windows path to its /mnt/<drive>/... WSL equivalent."""
    p = Path(win_path).resolve()
    parts = p.parts
    if len(parts) < 1 or ":" not in parts[0]:
        return str(p).replace("\\", "/")
    drive = parts[0][0].lower()
    rest = "/".join(parts[1:])
    return f"/mnt/{drive}/{rest}"


# ----- Public status API -----

def is_db_ready() -> bool:
    """True if the foldseek-indexed maize DB is built and usable."""
    return DB_MARKER.exists() and any(INDEX_DIR.glob("maize_db*"))


def pdb_count() -> int:
    """Number of extracted maize .pdb files on disk (after gunzip)."""
    if not PDB_DIR.exists():
        return 0
    return sum(1 for _ in PDB_DIR.glob("AF-*-model_*.pdb"))


def get_status() -> dict:
    """Snapshot for the /api/maize_afdb/status endpoint."""
    return {
        "ready": is_db_ready(),
        "pdb_dir": str(PDB_DIR),
        "index_dir": str(INDEX_DIR),
        "pdb_count": pdb_count(),
        "index_marker_present": DB_MARKER.exists(),
        "source_url": MAIZE_PROTEOME_URL,
    }


# ----- Progress callback plumbing -----

async def _emit(cb: Optional[ProgressCB], message: str, stage: Optional[str] = None, pct: Optional[float] = None):
    if cb is None:
        print(message, flush=True)
        return
    res = cb(message, stage=stage, pct=pct)
    if asyncio.iscoroutine(res):
        await res


# ----- Step 1: download -----

async def download_proteome(progress: Optional[ProgressCB] = None) -> Path:
    """Stream-download the maize AlphaFold proteome tar."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if TAR_PATH.exists() and TAR_PATH.stat().st_size > 1_000_000_000:
        await _emit(progress, f"Reusing existing tar at {TAR_PATH} ({TAR_PATH.stat().st_size / 1e9:.2f} GB)")
        return TAR_PATH

    await _emit(progress, f"Downloading {MAIZE_PROTEOME_URL} (this can take many minutes)", stage="download", pct=0.0)

    # Append to a partial file so re-runs resume — only if the server supports Range.
    resume_from = TAR_PATH.stat().st_size if TAR_PATH.exists() else 0
    headers = {"Range": f"bytes={resume_from}-"} if resume_from else {}
    mode = "ab" if resume_from else "wb"

    async with httpx.AsyncClient(follow_redirects=True, timeout=None) as client:
        async with client.stream("GET", MAIZE_PROTEOME_URL, headers=headers) as resp:
            if resp.status_code not in (200, 206):
                raise RuntimeError(f"Download failed: HTTP {resp.status_code} {resp.reason_phrase}")
            total = int(resp.headers.get("content-length", 0)) + resume_from
            downloaded = resume_from
            last_emit = 0
            with open(TAR_PATH, mode) as f:
                async for chunk in resp.aiter_bytes(chunk_size=8 * 1024 * 1024):
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total and downloaded - last_emit > 50 * 1024 * 1024:
                        pct = 100 * downloaded / total
                        await _emit(
                            progress,
                            f"Downloaded {downloaded / 1e9:.2f} / {total / 1e9:.2f} GB ({pct:.1f}%)",
                            stage="download", pct=pct,
                        )
                        last_emit = downloaded

    await _emit(progress, f"Download complete: {TAR_PATH.stat().st_size / 1e9:.2f} GB", stage="download", pct=100.0)
    return TAR_PATH


# ----- Step 2: extract + gunzip -----

async def extract_and_decompress(progress: Optional[ProgressCB] = None) -> int:
    """Extract the tar, gunzip *.pdb.gz, drop unneeded metadata files."""
    if PDB_DONE_MARKER.exists():
        n = pdb_count()
        await _emit(progress, f"Reusing existing extraction ({n} PDBs).", stage="extract", pct=100.0)
        return n

    PDB_DIR.mkdir(parents=True, exist_ok=True)
    await _emit(progress, "Extracting tar archive...", stage="extract", pct=0.0)

    def _extract():
        with tarfile.open(TAR_PATH) as tar:
            tar.extractall(PDB_DIR)
    await asyncio.to_thread(_extract)

    await _emit(progress, "Gunzipping PDB models and pruning non-PDB files...", stage="extract", pct=40.0)

    def _post_process() -> int:
        count = 0
        # Decompress PDB models only
        for gz in PDB_DIR.glob("AF-*-model_*.pdb.gz"):
            pdb = gz.with_suffix("")
            if not pdb.exists():
                with gzip.open(gz, "rb") as fin, open(pdb, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
            gz.unlink(missing_ok=True)
            count += 1
        # Drop other AF artefacts to save disk (~2x smaller)
        for pattern in ("AF-*-confidence_*.json.gz", "AF-*-predicted_aligned_error_*.json.gz",
                        "AF-*-model_*.cif.gz", "AF-*-model_*.cif"):
            for f in PDB_DIR.glob(pattern):
                f.unlink(missing_ok=True)
        return count

    n = await asyncio.to_thread(_post_process)
    PDB_DONE_MARKER.write_text(str(n))
    await _emit(progress, f"Extracted and decompressed {n} maize PDB models.", stage="extract", pct=100.0)
    return n


# ----- Step 3: build foldseek index -----

async def build_foldseek_index(progress: Optional[ProgressCB] = None) -> None:
    """Run `foldseek createdb` + `createindex` on the extracted PDBs."""
    INDEX_DIR.mkdir(parents=True, exist_ok=True)

    pdb_dir_wsl = to_wsl_path(PDB_DIR)
    db_prefix_wsl = to_wsl_path(DB_PREFIX)
    tmp_dir = INDEX_DIR / "tmp_createindex"
    tmp_dir.mkdir(exist_ok=True)
    tmp_dir_wsl = to_wsl_path(tmp_dir)

    await _emit(progress, "Running foldseek createdb (this can take 5–10 minutes)...", stage="index", pct=10.0)

    create_cmd = ["wsl", "foldseek", "createdb", pdb_dir_wsl, db_prefix_wsl]
    res = await asyncio.to_thread(subprocess.run, create_cmd, capture_output=True, timeout=3600)
    if res.returncode != 0:
        raise RuntimeError(
            f"foldseek createdb failed (rc={res.returncode}): "
            f"{res.stderr.decode(errors='replace')[:500]}"
        )

    await _emit(progress, "Running foldseek createindex (this can take 5–10 more minutes)...", stage="index", pct=60.0)
    index_cmd = ["wsl", "foldseek", "createindex", db_prefix_wsl, tmp_dir_wsl]
    res = await asyncio.to_thread(subprocess.run, index_cmd, capture_output=True, timeout=3600)
    if res.returncode != 0:
        # createindex is an optimization, not strictly required — log but don't abort.
        logger.warning(
            f"foldseek createindex returned rc={res.returncode}; the DB is still queryable. "
            f"stderr={res.stderr.decode(errors='replace')[:300]}"
        )
        await _emit(progress, "createindex skipped (queries will still work, just slower).", stage="index", pct=90.0)

    # Clean up createindex scratch dir
    shutil.rmtree(tmp_dir, ignore_errors=True)

    await _emit(progress, "Foldseek index built successfully.", stage="index", pct=100.0)


# ----- Optional cleanup -----

def cleanup_intermediates(keep_pdbs: bool = False, keep_tar: bool = False) -> None:
    """Remove the tar and (optionally) the unpacked PDB dir after indexing."""
    if not keep_tar and TAR_PATH.exists():
        TAR_PATH.unlink(missing_ok=True)
    if not keep_pdbs and PDB_DIR.exists():
        shutil.rmtree(PDB_DIR, ignore_errors=True)


# ----- High-level orchestration -----

async def ensure_db_ready(
    interactive: bool = True,
    progress_callback: Optional[ProgressCB] = None,
    cleanup: bool = True,
) -> bool:
    """
    Make sure the maize foldseek index is ready. If it isn't, download +
    extract + index. In interactive mode prompts the user first.
    Returns True iff the DB ended up ready.
    """
    if is_db_ready():
        await _emit(progress_callback, f"Maize foldseek index already present at {INDEX_DIR}.")
        return True

    if interactive:
        print()
        print("=" * 64)
        print(" Maize AlphaFold Structural Discovery Database — Setup")
        print("=" * 64)
        print(f" Source:    {MAIZE_PROTEOME_URL}")
        print(f" Disk:      ~5 GB tar download + ~3 GB extracted PDBs + ~1–2 GB foldseek index")
        print(f" Time:      ~20–40 minutes (download dominates)")
        print(f" Location:  {DATA_DIR}")
        print()
        resp = input("Download + index now? (y/n): ").strip().lower()
        if resp not in ("y", "yes"):
            print("Skipped. Structure-guided discovery (Phase 4.5) will be disabled.")
            return False

    if not _wsl_foldseek_available():
        msg = (
            "WSL foldseek binary not found — cannot build the index. "
            "Run setup.ps1 (or follow INSTALLATION.md) first."
        )
        await _emit(progress_callback, msg, stage="error")
        return False

    try:
        await download_proteome(progress_callback)
        await extract_and_decompress(progress_callback)
        await build_foldseek_index(progress_callback)
        if cleanup:
            cleanup_intermediates(keep_pdbs=False, keep_tar=False)
            await _emit(progress_callback, "Removed intermediate tar + raw PDBs to reclaim disk.")
    except Exception as e:
        await _emit(progress_callback, f"Setup failed: {type(e).__name__}: {e}", stage="error")
        logger.exception("Maize AFDB setup failed")
        return False

    return is_db_ready()


def _wsl_foldseek_available() -> bool:
    try:
        r = subprocess.run(["wsl", "foldseek", "version"],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=15)
        return r.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


# ----- CLI entry -----

if __name__ == "__main__":
    sys.exit(0 if asyncio.run(ensure_db_ready(interactive=True)) else 1)
