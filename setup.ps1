<#
.SYNOPSIS
    Automated Setup Script for DBsearch (Maize-Homology Pipeline)
    This script sets up WSL, installs the required Linux binaries (HMMER & Foldseek),
    and sets up the Windows Python backend.
#>

Write-Host "========================================================" -ForegroundColor Cyan
Write-Host "      DBsearch Automated Setup & Installation           " -ForegroundColor Cyan
Write-Host "========================================================" -ForegroundColor Cyan
Write-Host ""

# 1. Check WSL Installation
Write-Host "[1/7] Checking Windows Subsystem for Linux (WSL)..." -ForegroundColor Yellow
try {
    $null = wsl --status 2>&1
    if ($LASTEXITCODE -eq 0) {
        Write-Host "WSL is already installed and configured." -ForegroundColor Green
    } else {
        Write-Host "WSL is not installed or not set as default." -ForegroundColor Red
        Write-Host "Please run PowerShell as Administrator and execute: wsl --install" -ForegroundColor Yellow
        Write-Host "After your computer restarts, re-run this setup.ps1 script." -ForegroundColor Yellow
        exit
    }
} catch {
    Write-Host "WSL is not installed." -ForegroundColor Red
    Write-Host "Please run PowerShell as Administrator and execute: wsl --install" -ForegroundColor Yellow
    exit
}

# 2. Install Linux Dependencies inside WSL (HMMER & Foldseek)
Write-Host "`n[2/7] Installing HMMER (phmmer) and Foldseek inside WSL Ubuntu..." -ForegroundColor Yellow
try {
    Write-Host "Updating WSL package lists..."
    wsl -u root apt-get update -y
    Write-Host "Installing HMMER..."
    wsl -u root apt-get install -y hmmer curl tar
    
    Write-Host "Installing Foldseek..."
    wsl -u root bash -c "curl -sL https://github.com/steineggerlab/foldseek/releases/download/10-941cd33/foldseek-linux-avx2.tar.gz | tar xz -C /tmp && cp /tmp/foldseek/bin/foldseek /usr/local/bin/ && rm -rf /tmp/foldseek"
    
    Write-Host "Verifying Foldseek installation..."
    $fsCheck = wsl foldseek -h 2>&1
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Linux binaries installed successfully." -ForegroundColor Green
    } else {
        Write-Host "Warning: Foldseek might not have installed properly." -ForegroundColor Yellow
    }
} catch {
    Write-Host "Failed to install Linux binaries inside WSL. Please check WSL permissions." -ForegroundColor Red
}

# 3. Setup Python Virtual Environment (Windows Native)
Write-Host "`n[3/7] Setting up Windows Python Environment..." -ForegroundColor Yellow
$backendDir = Join-Path $PWD "backend"

if (-not (Test-Path $backendDir)) {
    Write-Host "Error: Cannot find 'backend' directory. Make sure you are running this from the root of the DBsearch repo." -ForegroundColor Red
    exit
}

Set-Location $backendDir

if (-not (Test-Path "venv")) {
    Write-Host "Creating Python virtual environment..."
    python -m venv venv
} else {
    Write-Host "Virtual environment 'venv' already exists." -ForegroundColor Green
}

# 4. Install Python Requirements
Write-Host "`n[4/7] Installing Python Dependencies..." -ForegroundColor Yellow
try {
    .\venv\Scripts\python.exe -m pip install --upgrade pip
    .\venv\Scripts\pip.exe install -r requirements.txt
    Write-Host "Python dependencies installed successfully." -ForegroundColor Green
} catch {
    Write-Host "Failed to install Python dependencies." -ForegroundColor Red
}

# 5. Download Zea mays Proteome Cache
Write-Host "`n[5/7] Downloading Zea mays Proteome for HMMER (approx. 12MB)..." -ForegroundColor Yellow
try {
    # Using python script to trigger the download logic inside hmmer_runner.py
    $cacheScript = @"
import asyncio
import os
from hmmer_runner import download_maize_proteome

async def main():
    data_dir = os.path.join(os.getcwd(), 'data')
    os.makedirs(data_dir, exist_ok=True)
    await download_maize_proteome(data_dir)
    print('Proteome successfully downloaded and extracted.')

asyncio.run(main())
"@
    $cacheScript | Out-File -FilePath "download_cache.py" -Encoding UTF8
    .\venv\Scripts\python.exe download_cache.py
    Remove-Item "download_cache.py"
} catch {
    Write-Host "Warning: Failed to pre-download proteome. The app will download it on the first HMMER run." -ForegroundColor Yellow
}

# 6. (Optional) Build the maize AlphaFold + Foldseek structural-discovery index
Write-Host "`n[6/7] Maize AlphaFold structural-discovery index..." -ForegroundColor Yellow
Write-Host "      This powers Phase 4.5 (structure-guided ortholog discovery)." -ForegroundColor Gray
Write-Host "      Cost: ~5 GB tar download + ~1-2 GB final foldseek index + 20-40 min one-time build." -ForegroundColor Gray
Write-Host "      You can skip now and build it later via the UI button or:" -ForegroundColor Gray
Write-Host "          .\venv\Scripts\python.exe install_maize_afdb.py" -ForegroundColor Gray
$buildAfdb = Read-Host "      Build the maize structural index now? (y/N)"
if ($buildAfdb -match '^(y|Y|yes|YES)$') {
    try {
        .\venv\Scripts\python.exe install_maize_afdb.py
    } catch {
        Write-Host "      Maize AFDB build failed. You can retry later via the UI." -ForegroundColor Yellow
    }
} else {
    Write-Host "      Skipped. Phase 4.5 will be disabled until the index is built." -ForegroundColor Gray
}

# 7. (Optional) CornCyc curated maize PGDB
$repoRoot = Split-Path -Parent $PSCommandPath
$cornCycRoot = Join-Path $repoRoot "corncyc"
Write-Host "`n[7/7] CornCyc curated maize pathway annotation..." -ForegroundColor Yellow
Write-Host "      This powers the curated discovery lane (4th lane) and the" -ForegroundColor Gray
Write-Host "      'CornCyc maize pathway context' section in the dashboard." -ForegroundColor Gray
Write-Host "      Cost: ~50-150 MB download + license-acceptance form on plantcyc.org." -ForegroundColor Gray
Write-Host "      Time: ~5 min (license form is the slow part)." -ForegroundColor Gray
Write-Host ""

if (Test-Path (Join-Path $cornCycRoot "default-version")) {
    Write-Host "      CornCyc already detected at $cornCycRoot — skipping." -ForegroundColor Green
} else {
    Write-Host "      You can skip now and install CornCyc later via the UI banner's" -ForegroundColor Gray
    Write-Host "      'Install instructions' button, or by following the manual steps" -ForegroundColor Gray
    Write-Host "      in INSTALLATION.md." -ForegroundColor Gray
    $wantCornCyc = Read-Host "      Open the PMN download page now and walk through CornCyc install? (y/N)"
    if ($wantCornCyc -match '^(y|Y|yes|YES)$') {
        Write-Host ""
        Write-Host "      Opening https://www.plantcyc.org/downloads in your browser..." -ForegroundColor Cyan
        try { Start-Process "https://www.plantcyc.org/downloads" } catch { Write-Host "      (couldn't open browser; visit the URL manually)" -ForegroundColor Yellow }
        Write-Host ""
        Write-Host "      Next steps (in the browser):" -ForegroundColor White
        Write-Host "        1. Sign in / register with PMN (free for non-commercial use)" -ForegroundColor Gray
        Write-Host "        2. Accept the license agreement for CornCyc" -ForegroundColor Gray
        Write-Host "        3. Download the tarball (typically corncyc-13.0.0.tar.gz)" -ForegroundColor Gray
        Write-Host ""
        Write-Host "      Then extract and move it into this repo:" -ForegroundColor White
        Write-Host "        tar.exe -xzf corncyc-13.0.0.tar.gz" -ForegroundColor Yellow
        Write-Host "        Move-Item corncyc $cornCycRoot" -ForegroundColor Yellow
        Write-Host ""
        Write-Host "      Final directory layout should look like:" -ForegroundColor White
        Write-Host "        $cornCycRoot\13.0.0\data\compounds.dat" -ForegroundColor Gray
        Write-Host ""
        Read-Host "      Press Enter once you've finished extracting (or just press Enter to skip)"

        # Verify
        if (Test-Path $cornCycRoot) {
            try {
                $statusJson = .\venv\Scripts\python.exe -c "import sys; sys.path.insert(0,'.'); import json; from corncyc_loader import get_status; print(json.dumps(get_status()))"
                $status = $statusJson | ConvertFrom-Json
                if ($status.available) {
                    Write-Host "      ✓ CornCyc detected at $($status.data_dir)" -ForegroundColor Green
                    if ($status.compounds) {
                        Write-Host "        $($status.compounds) compounds, $($status.reactions) reactions, $($status.pathways) pathways, $($status.maize_genes) maize genes" -ForegroundColor Gray
                    }
                } else {
                    Write-Host "      ⚠ CornCyc not yet detected at $cornCycRoot. You can drop it there later" -ForegroundColor Yellow
                    Write-Host "        and click 'Check again' in the UI banner — no restart needed." -ForegroundColor Yellow
                }
            } catch {
                Write-Host "      (couldn't auto-verify; restart the app and check the banner)" -ForegroundColor Yellow
            }
        } else {
            Write-Host "      No $cornCycRoot found yet. Extract the tarball there later and the app" -ForegroundColor Yellow
            Write-Host "      will auto-detect on next start (or click 'Check again' in the banner)." -ForegroundColor Yellow
        }
    } else {
        Write-Host "      Skipped. The curated lane will stay disabled until CornCyc is installed." -ForegroundColor Gray
    }
}

Write-Host "`n========================================================" -ForegroundColor Cyan
Write-Host "                  SETUP COMPLETE!                       " -ForegroundColor Green
Write-Host "========================================================" -ForegroundColor Cyan
Write-Host "To launch the DBsearch app at any time, simply double-click or run:" -ForegroundColor White
Write-Host "    .\run.ps1" -ForegroundColor Yellow
Write-Host ""
