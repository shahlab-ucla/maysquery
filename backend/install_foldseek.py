import os
import platform
import subprocess
import urllib.request
import tarfile
import zipfile
import shutil
import logging

logger = logging.getLogger(__name__)

def is_foldseek_installed() -> bool:
    try:
        subprocess.run(["foldseek", "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False

def install_foldseek():
    """Downloads and installs the correct Foldseek binary for the platform."""
    sys_plat = platform.system().lower()
    machine = platform.machine().lower()
    
    url = ""
    archive_ext = ""
    
    if sys_plat == "windows":
        url = "https://mmseqs.com/foldseek/foldseek-windows-avx2.zip"
        archive_ext = ".zip"
    elif sys_plat == "darwin": # macOS
        url = "https://mmseqs.com/foldseek/foldseek-mac-universal.tar.gz"
        archive_ext = ".tar.gz"
    elif sys_plat == "linux":
        url = "https://mmseqs.com/foldseek/foldseek-linux-sse2.tar.gz"
        archive_ext = ".tar.gz"
    else:
        print(f"Unsupported OS for auto-install: {sys_plat}")
        return False
        
    tmp_archive = os.path.join(os.path.dirname(__file__), f"foldseek_tmp{archive_ext}")
    print(f"Downloading Foldseek from {url}...")
    
    try:
        urllib.request.urlretrieve(url, tmp_archive)
        print("Download complete. Extracting...")
        
        extract_dir = os.path.join(os.path.dirname(__file__), "foldseek_extracted")
        os.makedirs(extract_dir, exist_ok=True)
        
        if archive_ext == ".zip":
            with zipfile.ZipFile(tmp_archive, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
        elif archive_ext == ".tar.gz":
            with tarfile.open(tmp_archive, 'r:gz') as tar_ref:
                tar_ref.extractall(extract_dir)
                
        # Move the binary to the active venv
        # The archive usually creates a 'foldseek' folder containing 'bin/foldseek'
        # e.g., foldseek_extracted/foldseek/bin/foldseek(.exe)
        bin_source_dir = os.path.join(extract_dir, "foldseek", "bin")
        if not os.path.exists(bin_source_dir):
            # Try just extract_dir in case it extracted flat
            bin_source_dir = extract_dir
            
        executable_name = "foldseek.exe" if sys_plat == "windows" else "foldseek"
        source_exe = os.path.join(bin_source_dir, executable_name)
        
        if not os.path.exists(source_exe):
            print(f"Could not find {executable_name} in extracted files.")
            return False
            
        # Determine target venv bin path
        # If running inside a venv, sys.prefix is the venv dir
        import sys
        if sys_plat == "windows":
            target_bin_dir = os.path.join(sys.prefix, "Scripts")
        else:
            target_bin_dir = os.path.join(sys.prefix, "bin")
            
        if not os.path.exists(target_bin_dir):
            os.makedirs(target_bin_dir, exist_ok=True)
            
        target_exe = os.path.join(target_bin_dir, executable_name)
        
        print(f"Installing {executable_name} to {target_bin_dir}...")
        shutil.copy2(source_exe, target_exe)
        
        if sys_plat != "windows":
            os.chmod(target_exe, 0o755)
            
        print("Foldseek successfully installed!")
        
    except Exception as e:
        print(f"Failed to install Foldseek: {e}")
        return False
    finally:
        # Cleanup
        if os.path.exists(tmp_archive):
            os.remove(tmp_archive)
        if os.path.exists(extract_dir):
            shutil.rmtree(extract_dir)
            
    return True

def check_and_prompt_foldseek():
    """Checks for Foldseek and prompts user if missing."""
    if not is_foldseek_installed():
        print("\n--- FOLDSEEK REQUIRED ---")
        print("Foldseek structural alignment binary was not found on your system PATH.")
        print("This is required for Phase 5 of the pipeline to perform 3D structural validation.")
        resp = input("Would you like to automatically download and install Foldseek to your virtual environment? (y/n): ")
        if resp.lower().strip() in ['y', 'yes']:
            install_foldseek()
        else:
            print("Skipping installation. Phase 5 alignments will fail.")
            
if __name__ == "__main__":
    check_and_prompt_foldseek()
