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
"${SCRIPT_DIR}/.venv/bin/flet" publish "${SCRIPT_DIR}/src/main.py" --distpath "${DIST_DIR}" "$@"
echo "==> Injecting browser unload warning script..."
"${SCRIPT_DIR}/.venv/bin/python" -c "
import pathlib
index_path = pathlib.Path('${DIST_DIR}/index.html')
if index_path.exists():
    content = index_path.read_text(encoding='utf-8')
    js_snippet = '<script>\nwindow.addEventListener(\"beforeunload\", (event) => {\n    event.preventDefault();\n    event.returnValue = \"\";\n});\n</script>\n'
    if '</head>' in content and js_snippet not in content:
        content = content.replace('</head>', js_snippet + '</head>')
        index_path.write_text(content, encoding='utf-8')
        print('Successfully injected beforeunload script into index.html')
"
echo "==> Build complete: ${DIST_DIR}"
echo "    To serve locally:"
echo "    python -m http.server 23455 -d ${DIST_DIR}"
