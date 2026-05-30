#!/usr/bin/env bash
# Build the Flet static site for local testing.
#
# Usage: ./build_static.sh
# Then open http://localhost:23455 in your browser.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Place dist *outside* src/ so flet publish never packages the previous
# build's app.tar.gz back into the new tarball.
DIST_DIR="${SCRIPT_DIR}/dist"

echo "==> Clearing Python bytecode cache..."
find "${SCRIPT_DIR}/src" -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

echo "==> Building static site..."
flet publish "${SCRIPT_DIR}/src/main.py" --distpath "${DIST_DIR}" "$@"
echo "==> Build complete: ${DIST_DIR}"
echo "    To serve locally:"
echo "    python -m http.server 23455 -d ${DIST_DIR}"
