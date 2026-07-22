#!/usr/bin/env bash
# Set up Linux development environment for AmplifyP.
#
# Usage: ./scripts/setup_linux.sh [--system-deps]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
cd "${PROJECT_DIR}"

HEADER="==> "
BLUE="\033[0;34m"
GREEN="\033[0;32m"
NC="\033[0m"

log() { echo -e "${HEADER}${BLUE}$*${NC}"; }
ok()  { echo -e "${HEADER}${GREEN}✓$*${NC}"; }

# Parse flags
INSTALL_SYSTEM_DEPS=false
SKIP_VENV=false
for arg in "$@"; do
  case "$arg" in
    --system-deps)  INSTALL_SYSTEM_DEPS=true ;;
    --skip-venv)    SKIP_VENV=true ;;
    *)              echo "Unknown option: $arg"; exit 1 ;;
  esac
done

# 1. System dependencies (optional, explicit opt-in)
if [[ "$INSTALL_SYSTEM_DEPS" = true ]]; then
  log "Installing system dependencies (apt)..."

  if ! command -v sudo &>/dev/null; then
    SUDO=""
  else
    SUDO="sudo"
  fi

  $SUDO apt-get update -qq
  $SUDO apt-get install -y --no-install-recommends \
    binutils clang cmake llvm lld ninja-build pkg-config \
    libgtk-3-dev libsecret-1-0 libsecret-1-dev libunwind-dev \
    libasound2-dev libgstreamer1.0-dev \
    libgstreamer-plugins-base1.0-dev libgstreamer-plugins-bad1.0-dev \
    libmpv-dev mpv wget \
    tesseract-ocr

  # Install Playwright system deps
  log "Installing Playwright system dependencies..."
  if [[ -d ".venv" ]]; then
    source .venv/bin/activate 2>/dev/null || true
    if command -v playwright &>/dev/null; then
      playwright install --with-deps chromium 2>/dev/null || true
    fi
  fi

  ok "System dependencies installed"
fi

# 2. Python version check
log "Checking Python version..."
PYTHON_CMD="${PYTHON_CMD:-python3}"

if ! "$PYTHON_CMD" -c "import sys; exit(0 if sys.version_info >= (3,12) else 1)" 2>/dev/null; then
  echo "Error: Python 3.12+ is required (found: $($PYTHON_CMD --version 2>&1))"
  exit 1
fi

ok "Python $($PYTHON_CMD --version 2>&1)"

# 3. Virtual environment
if [[ "$SKIP_VENV" = false ]]; then
  if [[ -d ".venv" ]]; then
    log "Virtual environment already exists at .venv"
  else
    log "Creating virtual environment at .venv..."
    "$PYTHON_CMD" -m venv .venv
    ok "Virtual environment created"
  fi

  log "Activating virtual environment..."
  # shellcheck disable=SC1091
  source .venv/bin/activate
  ok "Virtual environment active"
else
  if [[ -z "${VIRTUAL_ENV:-}" ]]; then
    echo "Error: --skip-venv specified but no VIRTUAL_ENV set. Activate a venv first."
    exit 1
  fi
  log "Using existing virtual environment: ${VIRTUAL_ENV}"
fi

# 4. Upgrade pip
log "Upgrading pip..."
pip install --upgrade pip -q
ok "pip $(pip --version)"

# 5. Install package in editable mode
log "Installing AmplifyP in editable mode with dev and e2e dependencies..."
pip install -e ".[dev,e2e]" -q
ok "Package installed"

# 6. Playwright browsers
log "Installing Playwright browser binaries..."
if command -v playwright &>/dev/null; then
  playwright install chromium
  playwright install-deps chromium 2>/dev/null || true
  ok "Playwright browsers installed"
else
  echo "  Skipping Playwright browsers (playwright not installed)"
fi

# 7. Git hooks with prek
if command -v prek &>/dev/null; then
  log "Setting up prek git hook..."
  if [[ ! -d ".git/hooks" ]]; then
    echo "  Skipping: .git/hooks not found (not a git repo?)"
  else
    # prek can manage its own hooks; just verify it works
    if prek run --all-files --dry-run &>/dev/null || true; then
      ok "prek available for hook checks"
    else
      echo "  prek installed but dry-run failed (may need system deps)"
    fi
  fi
else
  echo "  prek not installed (skip git hooks setup)"
fi

# 8. Verify installation
log "Verifying installation..."
verify_ok=true

for cmd in pytest pyright ruff prek playwright; do
  if command -v "$cmd" &>/dev/null; then
    ok "$cmd available"
  else
    echo "  ✗ $cmd not found"
    verify_ok=false
  fi
done

# 9. Summary
echo ""
echo "========================================="
echo "  Development environment ready"
echo "========================================="
echo ""
echo "  Quick commands:"
echo "    Run tests:        pytest"
echo "    Run e2e tests:    playwright install && pytest --run-slow -m e2e"
echo "    Type check:       pyright"
echo "    Lint/format:      prek run --all-files"
echo "    Run GUI (desktop): flet run -r src/main.py"
echo "    Run GUI (web):     flet run -w -r src/main.py"
echo ""
echo "  To activate:        source .venv/bin/activate"
echo ""

if [[ "$verify_ok" = false ]]; then
  echo "  Warning: Some tools missing. Check output above."
  echo ""
fi
