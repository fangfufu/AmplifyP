#!/usr/bin/env bash
# Build the Flet Linux binary and package as tar.gz and AppImage.
#
# Usage: ./build_linux.sh [--install-deps]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Disable rich output and skip Flutter Doctor Android toolchain checks
export FLET_CLI_SKIP_FLUTTER_DOCTOR=1
if [[ "${CI:-}" = "true" ]]; then
  export FLET_CLI_NO_RICH_OUTPUT=1
fi

RESTORE_FLUTTER_ANDROID=false
if command -v flutter >/dev/null 2>&1; then
  FLUTTER_CFG="$(flutter config 2>&1 || true)"
  if [[ "${FLUTTER_CFG}" != *"enable-android: false"* ]]; then
    RESTORE_FLUTTER_ANDROID=true
    flutter config --no-enable-android >/dev/null 2>&1 || true
  fi
fi

cleanup() {
  if [[ "${RESTORE_FLUTTER_ANDROID:-false}" = "true" ]] && command -v flutter >/dev/null 2>&1; then
    flutter config --enable-android >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT INT TERM

# Source the virtual environment if it exists and is not already sourced
if [[ -z "${VIRTUAL_ENV:-}" ]] && [[ -d ".venv" ]]; then
  echo "==> Sourcing virtual environment..."
  # shellcheck disable=SC1091
  source .venv/bin/activate
fi

FLET_BIN="flet"
if [[ -n "${VIRTUAL_ENV:-}" ]]; then
  if [[ ! -x "${VIRTUAL_ENV}/bin/flet" ]]; then
    echo "==> Installing project dependencies into virtual environment..."
    pip install -e ".[dev]"
  fi
  FLET_BIN="${VIRTUAL_ENV}/bin/flet"
elif ! command -v flet >/dev/null 2>&1; then
  echo "Error: 'flet' command not found. Please activate a virtual environment or install dependencies." >&2
  exit 1
fi

INSTALL_DEPS=false
for arg in "$@"; do
  if [[ "$arg" = "--install-deps" ]]; then
    INSTALL_DEPS=true
  fi
done

if [[ "$INSTALL_DEPS" = true ]]; then
  echo "==> Installing system dependencies..."
  sudo apt-get update
  sudo apt-get install -y --no-install-recommends \
    binutils clang cmake llvm lld ninja-build pkg-config \
    libgtk-3-dev libsecret-1-0 libsecret-1-dev libunwind-dev \
    libasound2-dev libgstreamer1.0-dev \
    libgstreamer-plugins-base1.0-dev libgstreamer-plugins-bad1.0-dev \
    libmpv-dev mpv wget
fi

echo "==> Generating Git SHA..."
python scripts/gen_git_sha.py

echo "==> Building Flet Linux binary..."
rm -rf build/linux build/AmplifyP
"${FLET_BIN}" build linux . -o build/linux --project AmplifyP --yes

echo "==> Moving build artefacts..."
mv build/linux build/AmplifyP

echo "==> Creating .desktop file..."
cat << 'EOF' | sed 's/^ *//' > build/AmplifyP/AmplifyP.desktop
[Desktop Entry]
Type=Application
Name=AmplifyP
Comment=Simulate Polymerase Chain Reaction (PCR) and predict DNA amplification products.
Exec=AmplifyP %U
Icon=AmplifyP
Categories=Science;Education;
Terminal=false
EOF

echo "==> Copying icons..."
cp src/assets/images/icon.png build/AmplifyP/AmplifyP.png
cp src/assets/images/icon.png build/AmplifyP/.DirIcon

echo "==> Creating AppRun script..."
cat << 'EOF' | sed 's/^ *//' > build/AmplifyP/AppRun
#!/bin/sh
HERE="$(dirname "$(readlink -f "${0}")")"
export LD_LIBRARY_PATH="${HERE}/lib${LD_LIBRARY_PATH+:${LD_LIBRARY_PATH}}"
exec "${HERE}/AmplifyP" "$@"
EOF
chmod +x build/AmplifyP/AppRun

echo "==> Packaging tarball..."
tar -czf amplifyp-linux.tar.gz -C build AmplifyP

echo "==> Packaging as AppImage..."
if [[ ! -f appimagetool-x86_64.AppImage ]]; then
  wget -q https://github.com/AppImage/appimagetool/releases/download/continuous/appimagetool-x86_64.AppImage -O appimagetool-x86_64.AppImage.tmp
  mv appimagetool-x86_64.AppImage.tmp appimagetool-x86_64.AppImage
  chmod +x appimagetool-x86_64.AppImage
fi

# Ensure squashfs-root is clean before extraction
rm -rf squashfs-root
./appimagetool-x86_64.AppImage --appimage-extract

export ARCH=x86_64
APPIMAGE_ARGS=()
if [[ -n "${UPDATE_INFORMATION:-}" ]]; then
  APPIMAGE_ARGS+=("-u" "${UPDATE_INFORMATION}")
fi
APPIMAGE_NAME="${APPIMAGE_NAME:-amplifyp${VERSION:+-${VERSION}}-x86_64.AppImage}"
./squashfs-root/AppRun "${APPIMAGE_ARGS[@]}" build/AmplifyP "${APPIMAGE_NAME}"

# Clean up extracted appimagetool dir and temporary flet build directories
rm -rf squashfs-root
rm -rf src/build src/dist
rm -f src/amplifyp/gui/git_sha.py

echo "==> Build complete: amplifyp-linux.tar.gz and ${APPIMAGE_NAME}"
