#!/usr/bin/env bash
# Build the Flet Linux binary and package as tar.gz and AppImage.
#
# Usage: ./build_linux.sh [--install-deps]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Source the virtual environment if it exists and is not already sourced
if [ -z "${VIRTUAL_ENV:-}" ] && [ -d ".venv" ]; then
  echo "==> Sourcing virtual environment..."
  # shellcheck disable=SC1091
  source .venv/bin/activate
fi

INSTALL_DEPS=false
for arg in "$@"; do
  if [ "$arg" = "--install-deps" ]; then
    INSTALL_DEPS=true
  fi
done

if [ "$INSTALL_DEPS" = true ]; then
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
flet build linux src -o build/linux --project AmplifyP --yes

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
export LD_LIBRARY_PATH="${HERE}/lib:${LD_LIBRARY_PATH}"
exec "${HERE}/AmplifyP" "$@"
EOF
chmod +x build/AmplifyP/AppRun

echo "==> Packaging tarball..."
tar -czf amplifyp-linux.tar.gz -C build AmplifyP

echo "==> Packaging as AppImage..."
if [ ! -f appimagetool-x86_64.AppImage ]; then
  wget -q https://github.com/AppImage/appimagetool/releases/download/continuous/appimagetool-x86_64.AppImage
  chmod +x appimagetool-x86_64.AppImage
fi

# Ensure squashfs-root is clean before extraction
rm -rf squashfs-root
./appimagetool-x86_64.AppImage --appimage-extract

export ARCH=x86_64
./squashfs-root/AppRun build/AmplifyP amplifyp-linux.AppImage

# Clean up extracted appimagetool dir and temporary flet build directories
rm -rf squashfs-root
rm -rf src/build src/dist

echo "==> Build complete: amplifyp-linux.tar.gz and amplifyp-linux.AppImage"
