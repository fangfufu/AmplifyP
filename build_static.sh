#!/usr/bin/env bash
# Build the Flet static site for local testing.
#
# Usage: ./build_static.sh
# Then open http://localhost:23455 in your browser.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIST_DIR="${SCRIPT_DIR}/src/dist"

echo "==> Building static site..."
flet publish "${SCRIPT_DIR}/src/main.py" "$@"

echo "==> Patching Flet static site for FilePicker support..."
python "${SCRIPT_DIR}/patch_flet_web.py" "${DIST_DIR}"

echo "==> Build complete: ${DIST_DIR}"
echo "    To serve locally:"
echo "    python -m http.server 23455 -d ${DIST_DIR}"
