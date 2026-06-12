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

# Temporarily copy build icons to root of assets so Flet publish auto-generates PWA icons
cp "${SCRIPT_DIR}/src/assets/images/icon.png" "${SCRIPT_DIR}/src/assets/icon.png" || true
cp "${SCRIPT_DIR}/src/assets/images/favicon.png" "${SCRIPT_DIR}/src/assets/favicon.png" || true

GIT_SHA=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
GIT_FULL_SHA=$(git rev-parse HEAD 2>/dev/null || echo "unknown")
echo "GIT_SHA = \"${GIT_SHA}\"" > "${SCRIPT_DIR}/src/amplifyp/gui/git_sha.py"
echo "GIT_FULL_SHA = \"${GIT_FULL_SHA}\"" >> "${SCRIPT_DIR}/src/amplifyp/gui/git_sha.py"

echo "==> Building static site..."
flet publish "${SCRIPT_DIR}/src/main.py" \
  --distpath "${DIST_DIR}" \
  --app-name "AmplifyP" \
  --app-short-name "AmplifyP" \
  --app-description "Simulate Polymerase Chain Reaction (PCR) and predict DNA amplification products." \
  --pwa-background-color "#FFFFFF" \
  --pwa-theme-color "#0175C2" \
  "$@"

# Clean up temporary root assets files and git_sha.py
rm -f "${SCRIPT_DIR}/src/assets/icon.png" "${SCRIPT_DIR}/src/assets/favicon.png" || true
rm -f "${SCRIPT_DIR}/src/amplifyp/gui/git_sha.py" || true


echo "==> Injecting browser unload warning script and git SHA..."
GIT_SHA=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
echo "${GIT_SHA}" > "${DIST_DIR}/.git-sha"
python -c "
import pathlib
index_path = pathlib.Path('${DIST_DIR}/index.html')
if index_path.exists():
    content = index_path.read_text(encoding='utf-8')
    js_snippet = '<script>\nwindow.addEventListener(\"beforeunload\", (event) => {\n    event.preventDefault();\n    event.returnValue = \"\";\n});\nwindow.__APP_SHA__ = \"${GIT_SHA}\";\n</script>\n'
    if '</head>' in content and js_snippet not in content:
        content = content.replace('</head>', js_snippet + '</head>')
        index_path.write_text(content, encoding='utf-8')
        print('Successfully injected beforeunload script into index.html')
"
echo "==> Copying custom assets to dist..."
mkdir -p "${DIST_DIR}/images" "${DIST_DIR}/assets/images"
cp -R "${SCRIPT_DIR}/src/assets/images/"* "${DIST_DIR}/images/" || true
cp -R "${SCRIPT_DIR}/src/assets/images/"* "${DIST_DIR}/assets/images/" || true

echo "==> Replacing Flet loading animation with app icon..."
mkdir -p "${DIST_DIR}/icons"
cp "${SCRIPT_DIR}/src/assets/images/icon.png" "${DIST_DIR}/icons/loading-animation.png"

echo "==> Build complete: ${DIST_DIR}"
echo "    To serve locally:"
echo "    python -m http.server 23455 -d ${DIST_DIR}"
