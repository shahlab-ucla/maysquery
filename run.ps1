<#
.SYNOPSIS
    Launcher Script for DBsearch
#>

Write-Host "========================================================" -ForegroundColor Cyan
Write-Host "                 Launching DBsearch                     " -ForegroundColor Cyan
Write-Host "========================================================" -ForegroundColor Cyan

$backendDir = Join-Path $PWD "backend"
if (-not (Test-Path $backendDir)) {
    Write-Host "Error: Cannot find 'backend' directory. Make sure you are running this from the DBsearch root." -ForegroundColor Red
    Pause
    exit
}

Set-Location $backendDir

if (-not (Test-Path "venv\Scripts\uvicorn.exe")) {
    Write-Host "Error: The virtual environment is missing or uvicorn is not installed." -ForegroundColor Red
    Write-Host "Please run .\setup.ps1 first." -ForegroundColor Yellow
    Pause
    exit
}

Write-Host "Starting FastAPI Backend server..." -ForegroundColor Green
Write-Host "The application will open in your default browser." -ForegroundColor White
Write-Host "Press CTRL+C in this terminal to stop the server." -ForegroundColor Yellow
Write-Host ""

# Open the browser to the static index page
Start-Process "http://127.0.0.1:8008/static/index.html"

# Run Uvicorn
.\venv\Scripts\uvicorn.exe main:app --reload --host 127.0.0.1 --port 8008
