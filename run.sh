#!/usr/bin/env bash
# Maysquery launcher for macOS and Linux.
#
# Boots the FastAPI server, optionally opens the UI in your default browser.
# Mirrors the Windows run.ps1.

set -e

cd "$(dirname "$0")/backend" || { echo "Cannot find backend/ next to this script"; exit 1; }

if [ ! -x venv/bin/uvicorn ]; then
  echo "venv is missing or uvicorn isn't installed."
  echo "Run ./setup.sh first."
  exit 1
fi

URL="http://127.0.0.1:8008/static/index.html"

# Try to open the browser (best-effort; harmless if it fails)
( sleep 1
  if   command -v xdg-open >/dev/null 2>&1; then xdg-open "$URL"
  elif command -v open     >/dev/null 2>&1; then open "$URL"
  fi ) &

echo "Starting Maysquery on $URL ..."
echo "Press Ctrl-C to stop."
exec ./venv/bin/uvicorn main:app --reload --host 127.0.0.1 --port 8008
