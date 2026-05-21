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
Write-Host "[1/6] Checking Windows Subsystem for Linux (WSL)..." -ForegroundColor Yellow
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
Write-Host "`n[2/6] Installing HMMER (phmmer) and Foldseek inside WSL Ubuntu..." -ForegroundColor Yellow
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
Write-Host "`n[3/6] Setting up Windows Python Environment..." -ForegroundColor Yellow
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
Write-Host "`n[4/6] Installing Python Dependencies..." -ForegroundColor Yellow
try {
    .\venv\Scripts\python.exe -m pip install --upgrade pip
    .\venv\Scripts\pip.exe install -r requirements.txt
    Write-Host "Python dependencies installed successfully." -ForegroundColor Green
} catch {
    Write-Host "Failed to install Python dependencies." -ForegroundColor Red
}

# 5. Download Zea mays Proteome Cache
Write-Host "`n[5/6] Downloading Zea mays Proteome for HMMER (approx. 12MB)..." -ForegroundColor Yellow
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
Write-Host "`n[6/6] Maize AlphaFold structural-discovery index..." -ForegroundColor Yellow
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

Write-Host "`n========================================================" -ForegroundColor Cyan
Write-Host "                  SETUP COMPLETE!                       " -ForegroundColor Green
Write-Host "========================================================" -ForegroundColor Cyan
Write-Host "To launch the DBsearch app at any time, simply double-click or run:" -ForegroundColor White
Write-Host "    .\run.ps1" -ForegroundColor Yellow
Write-Host ""
