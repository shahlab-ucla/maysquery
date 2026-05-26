#!/usr/bin/env bash
# Maysquery automated setup for macOS and Linux.
#
# Installs Foldseek + HMMER natively (no WSL needed on these platforms),
# creates a Python venv, installs Python deps, and optionally builds the
# maize AlphaFold structural-discovery index.
#
# This script is NOT tested in CI; treat it as a documented procedure rather
# than a turn-key installer. If a step fails, see INSTALLATION.md for the
# manual equivalent.

set -e

cyan()   { printf '\033[36m%s\033[0m\n' "$1"; }
yellow() { printf '\033[33m%s\033[0m\n' "$1"; }
green()  { printf '\033[32m%s\033[0m\n' "$1"; }
red()    { printf '\033[31m%s\033[0m\n' "$1"; }

cyan '========================================================'
cyan '   Maysquery Automated Setup & Installation (Unix)     '
cyan '========================================================'

# ---------- 1. Platform detection ----------
yellow '[1/7] Detecting platform...'
OS=$(uname -s)
case "$OS" in
  Darwin) PLATFORM=mac;   green "  -> macOS" ;;
  Linux)  PLATFORM=linux; green "  -> Linux ($(uname -m))" ;;
  *)      red "  -> Unsupported OS: $OS. Aborting."; exit 1 ;;
esac

# ---------- 2. Bio-binary install (Foldseek + HMMER) ----------
yellow '[2/7] Installing Foldseek + HMMER...'
need_install() { ! command -v "$1" >/dev/null 2>&1; }

if [ "$PLATFORM" = "mac" ]; then
  if ! command -v brew >/dev/null 2>&1; then
    red "  Homebrew is required on macOS. Install from https://brew.sh and re-run."
    exit 1
  fi
  if need_install foldseek; then brew install brewsci/bio/foldseek; else green "  foldseek already installed"; fi
  if need_install phmmer;   then brew install hmmer;                else green "  hmmer already installed";    fi
else
  # Linux — prefer apt where available, fall back to manual binary download
  if command -v apt-get >/dev/null 2>&1; then
    if need_install phmmer; then sudo apt-get update -y && sudo apt-get install -y hmmer curl tar; else green "  hmmer already installed"; fi
  elif command -v dnf >/dev/null 2>&1; then
    if need_install phmmer; then sudo dnf install -y hmmer curl tar; else green "  hmmer already installed"; fi
  else
    yellow "  No supported package manager. Install hmmer + curl + tar manually before continuing."
  fi
  if need_install foldseek; then
    yellow "  Downloading Foldseek..."
    TMPDIR=$(mktemp -d)
    # Pick AVX2 if /proc/cpuinfo says so, else SSE2
    if grep -q avx2 /proc/cpuinfo 2>/dev/null; then
      FOLDSEEK_BUILD="avx2"
    else
      FOLDSEEK_BUILD="sse2"
    fi
    curl -sL "https://github.com/steineggerlab/foldseek/releases/download/10-941cd33/foldseek-linux-${FOLDSEEK_BUILD}.tar.gz" \
      | tar xz -C "$TMPDIR"
    sudo cp "$TMPDIR/foldseek/bin/foldseek" /usr/local/bin/
    rm -rf "$TMPDIR"
  else
    green "  foldseek already installed"
  fi
fi

if command -v foldseek >/dev/null 2>&1 && command -v phmmer >/dev/null 2>&1; then
  green "  Bio-binaries installed."
else
  red "  WARNING: foldseek and/or phmmer not on PATH after install. The pipeline will fall back to mocked scores."
fi

# ---------- 3. Python venv ----------
yellow '[3/7] Setting up Python virtual environment...'
cd "$(dirname "$0")/backend" || { red "Cannot cd into ./backend"; exit 1; }

if [ ! -d venv ]; then
  python3 -m venv venv
  green "  Created venv"
else
  green "  venv already exists"
fi

# ---------- 4. Python deps ----------
yellow '[4/7] Installing Python dependencies...'
./venv/bin/python -m pip install --upgrade pip
./venv/bin/pip install -r requirements.txt
green "  Done."

# ---------- 5. Maize proteome cache (small, for HMMER fallback) ----------
yellow '[5/7] Pre-downloading the Zea mays reference proteome (~12 MB) for HMMER fallback...'
./venv/bin/python - <<'PY'
import asyncio, os
from hmmer_runner import download_maize_proteome
async def main():
    data_dir = os.path.join(os.getcwd(), 'data')
    os.makedirs(data_dir, exist_ok=True)
    await download_maize_proteome(data_dir)
    print('  Proteome downloaded and extracted.')
asyncio.run(main())
PY

# ---------- 6. Optional: maize AlphaFold structural-discovery index ----------
yellow '[6/7] Maize AlphaFold structural-discovery index'
echo  '      This powers Phase 4.5 (structure-guided ortholog discovery).'
echo  '      Cost: ~5 GB tar download + ~1-2 GB final foldseek index + 20-40 min one-time build.'
echo  '      You can skip now and build it later via the UI button or:'
echo  '          ./venv/bin/python install_maize_afdb.py'
read -r -p "      Build the maize structural index now? (y/N) " ans
case "$ans" in
  y|Y|yes|YES) ./venv/bin/python install_maize_afdb.py || yellow "      Build failed — retry later via the UI." ;;
  *)           echo  "      Skipped. Phase 4.5 will be disabled until the index is built." ;;
esac

# ---------- 7. Optional: CornCyc curated maize PGDB ----------
REPO_ROOT_ABS=$(cd .. && pwd)
CORNCYC_ROOT="$REPO_ROOT_ABS/corncyc"
yellow '[7/7] CornCyc curated maize pathway annotation'
echo  '      Powers the curated discovery lane + dashboard pathway-context section.'
echo  '      Cost: ~50-150 MB download + license-acceptance form on plantcyc.org.'

if [ -f "$CORNCYC_ROOT/default-version" ]; then
  green "      CornCyc already detected at $CORNCYC_ROOT — skipping."
else
  echo  '      You can skip now and install CornCyc later via the UI banner'
  echo  "      ('Install instructions' button) or by following INSTALLATION.md."
  read -r -p "      Open the PMN download page now and walk through install? (y/N) " ans
  case "$ans" in
    y|Y|yes|YES)
      cyan  '      Opening https://www.plantcyc.org/downloads in your browser...'
      if   command -v xdg-open >/dev/null 2>&1; then xdg-open "https://www.plantcyc.org/downloads" >/dev/null 2>&1 &
      elif command -v open     >/dev/null 2>&1; then open    "https://www.plantcyc.org/downloads" >/dev/null 2>&1 &
      else yellow "      (couldn't open browser; visit the URL manually)"
      fi
      echo
      echo  '      Next steps (in the browser):'
      echo  '        1. Sign in / register with PMN (free for non-commercial use)'
      echo  '        2. Accept the license agreement for CornCyc'
      echo  '        3. Download the tarball (typically corncyc-13.0.0.tar.gz)'
      echo
      echo  '      Then extract and move it into this repo:'
      yellow '        tar xzf corncyc-13.0.0.tar.gz'
      yellow "        mv corncyc \"$CORNCYC_ROOT\""
      echo
      echo  '      Final directory layout should look like:'
      echo  "        $CORNCYC_ROOT/13.0.0/data/compounds.dat"
      echo
      read -r -p "      Press Enter once you've finished extracting (or just press Enter to skip) "
      if [ -d "$CORNCYC_ROOT" ]; then
        STATUS_JSON=$(./venv/bin/python -c "import sys; sys.path.insert(0,'.'); import json; from corncyc_loader import get_status; print(json.dumps(get_status()))" 2>/dev/null || echo '{}')
        if echo "$STATUS_JSON" | grep -q '"available": true'; then
          green "      ✓ CornCyc detected"
          echo  "      $STATUS_JSON"
        else
          yellow "      ⚠ CornCyc not yet detected at $CORNCYC_ROOT — drop it there later and click 'Check again' in the UI."
        fi
      else
        yellow "      No $CORNCYC_ROOT found yet. Extract the tarball there later; the app will auto-detect."
      fi
      ;;
    *) echo  '      Skipped. The curated lane will stay disabled until CornCyc is installed.' ;;
  esac
fi

cyan '========================================================'
green '                  SETUP COMPLETE!                       '
cyan '========================================================'
echo "To launch:   ./run.sh"
echo "Or manually: cd backend && ./venv/bin/uvicorn main:app --reload --host 127.0.0.1 --port 8008"
