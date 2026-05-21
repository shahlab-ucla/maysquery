import os
import subprocess
import httpx
import gzip
import logging
import asyncio
from typing import List, Dict

logger = logging.getLogger(__name__)

HMMER_DB_URL = "http://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/zea_mays/pep/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.pep.all.fa.gz"

async def check_and_install_hmmer():
    """Checks if HMMER is installed via WSL, and installs it if missing."""
    try:
        res = subprocess.run(["wsl", "phmmer", "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if res.returncode == 0:
            return True
    except FileNotFoundError:
        logger.error("WSL is not installed or available.")
        return False

    logger.info("Installing HMMER in WSL...")
    try:
        # Run apt-get install as root without prompting for password (works in most default WSL setups)
        subprocess.run(["wsl", "-u", "root", "apt-get", "update"], check=True)
        subprocess.run(["wsl", "-u", "root", "apt-get", "install", "-y", "hmmer"], check=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to install HMMER via WSL: {e}")
        return False

async def download_maize_proteome(data_dir: str):
    fasta_gz_path = os.path.join(data_dir, "zea_mays_proteome.fasta.gz")
    fasta_path = os.path.join(data_dir, "zea_mays_proteome.fasta")

    if os.path.exists(fasta_path):
        return fasta_path

    logger.info("Downloading Maize Proteome for HMMER...")
    try:
        async with httpx.AsyncClient() as client:
            async with client.stream("GET", HMMER_DB_URL) as response:
                response.raise_for_status()
                with open(fasta_gz_path, "wb") as f:
                    async for chunk in response.aiter_bytes():
                        f.write(chunk)
                        
        logger.info("Extracting Proteome FASTA...")
        with gzip.open(fasta_gz_path, "rb") as f_in:
            with open(fasta_path, "wb") as f_out:
                f_out.write(f_in.read())
                
        os.remove(fasta_gz_path)
        return fasta_path
    except Exception as e:
        logger.error(f"Failed to download/extract Maize Proteome: {e}")
        return None

async def run_phmmer_search(query_sequence: str, query_id: str, e_value_cutoff: float = 1e-5, max_hits: int = None) -> List[Dict]:
    """
    Runs phmmer via WSL against the local Maize proteome.
    """
    if not await check_and_install_hmmer():
        return []

    data_dir = os.path.join(os.path.dirname(__file__), "data")
    os.makedirs(data_dir, exist_ok=True)
    
    db_path = await download_maize_proteome(data_dir)
    if not db_path:
        return []

    query_fasta = os.path.join(data_dir, f"{query_id}_query.fasta")
    results_tbl = os.path.join(data_dir, f"{query_id}_phmmer.tbl")
    
    with open(query_fasta, "w") as f:
        f.write(f">{query_id}\n{query_sequence}\n")

    # Run phmmer via WSL
    # Note: Using relative paths or WSL-translated paths if running inside data_dir
    # Since wsl preserves cwd, we can run it safely using basenames if we set cwd=data_dir
    cmd = [
        "wsl", "phmmer",
        "--tblout", os.path.basename(results_tbl),
        "-E", str(e_value_cutoff),
        os.path.basename(query_fasta),
        os.path.basename(db_path)
    ]
    
    logger.info(f"Running HMMER: {' '.join(cmd)}")

    try:
        # Use asyncio.to_thread for cross-platform reliability — Windows asyncio's
        # subprocess_exec is fragile under uvicorn's reload mode (NotImplementedError
        # with the default SelectorEventLoop). Synchronous subprocess.run in a
        # worker thread sidesteps the issue.
        result = await asyncio.to_thread(
            subprocess.run,
            cmd,
            cwd=data_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=300,
        )

        if result.returncode != 0:
            stderr_text = result.stderr.decode(errors="replace").strip()
            stdout_text = result.stdout.decode(errors="replace").strip()
            logger.error(
                f"phmmer failed (rc={result.returncode}). "
                f"stderr={stderr_text!r} stdout={stdout_text!r} cmd={' '.join(cmd)!r}"
            )
            return []

        hits = []
        if os.path.exists(results_tbl):
            with open(results_tbl, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) >= 6:
                        target_name = parts[0]
                        # Maize gene models in Ensembl are usually Zm0000...
                        # Format is like: Zm00001eb016240_T001. We want the gene model, so we split by '_'
                        gene_model = target_name.split('_')[0] if '_' in target_name else target_name

                        e_value = float(parts[4])
                        bit_score = float(parts[5])

                        if e_value <= e_value_cutoff:
                            # Avoid duplicate gene models from alternative transcripts
                            if not any(h['maize_gene_model'] == gene_model for h in hits):
                                hits.append({
                                    "maize_gene_model": gene_model,
                                    "e_value": e_value,
                                    "bit_score": bit_score
                                })

            # Sort by bit score descending
            hits.sort(key=lambda x: x["bit_score"], reverse=True)
            if max_hits is not None:
                hits = hits[:max_hits]

        return hits
    except FileNotFoundError as e:
        logger.error(
            f"phmmer launch failed — 'wsl' or 'phmmer' binary not found on PATH. "
            f"Detail: {e!r}. Cmd was: {' '.join(cmd)!r}"
        )
        return []
    except subprocess.TimeoutExpired:
        logger.error(f"phmmer timed out after 300s. Cmd: {' '.join(cmd)!r}")
        return []
    except Exception as e:
        logger.error(
            f"Error executing phmmer ({type(e).__name__}): {e!r}. Cmd: {' '.join(cmd)!r}"
        )
        return []
    finally:
        # Cleanup
        if os.path.exists(query_fasta):
            os.remove(query_fasta)
        if os.path.exists(results_tbl):
            os.remove(results_tbl)
