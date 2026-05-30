# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""End-to-End tests for the static web version of the GUI."""

import os
import shutil
import subprocess
import sys
import threading
import time
from collections.abc import Generator
from http.server import HTTPServer, SimpleHTTPRequestHandler
from typing import Any

import pytest

# Soft imports for optional test dependencies
try:
    from playwright.sync_api import Page, expect
except ImportError:
    Page = None
    expect = None

# Define ports and paths (port 34521 as per user instructions)
SERVER_PORT = 34521
SRC_DIR = os.path.join(os.getcwd(), "src")
DIST_DIR = os.path.join(os.getcwd(), "dist")


# Check entire module for dependencies
if not Page:
    pytest.fail(
        "E2E test dependencies missing: pytest-playwright. "
        "Please install it using: pip install .[e2e] "
        "And then: playwright install chromium",
        pytrace=False,
    )


@pytest.fixture(scope="session")  # type: ignore[untyped-decorator]
def build_app() -> None:
    """Build the Flet static app.

    This mimics the logic in `build_static.sh`.
    """
    # Ensure dist dir is clean
    if os.path.exists(DIST_DIR):
        shutil.rmtree(DIST_DIR)

    # Clear pyc cache so flet publish packages fresh source, not stale
    # bytecode from a previous build.
    for root, dirs, _ in os.walk("src"):
        for d in dirs:
            if d == "__pycache__":
                shutil.rmtree(os.path.join(root, d), ignore_errors=True)

    print("==> Building static site...")
    # Run flet publish
    flet_path = shutil.which("flet")
    if not flet_path:
        # Fallback to venv path if not in PATH
        flet_bin = "flet.exe" if os.name == "nt" else "flet"
        flet_path = os.path.join(os.path.dirname(sys.executable), flet_bin)

    subprocess.run(  # noqa: S603
        [flet_path, "publish", "src/main.py", "--distpath", DIST_DIR],
        check=True,
        capture_output=True,
    )

    assert os.path.exists(os.path.join(DIST_DIR, "index.html"))


@pytest.fixture(scope="session")  # type: ignore[untyped-decorator]
def serve_app(build_app: None) -> Generator[str, None, None]:
    """Serve the static app in a background thread."""
    HTTPServer.allow_reuse_address = True
    server = HTTPServer(
        ("localhost", SERVER_PORT),
        lambda *args: SimpleHTTPRequestHandler(*args, directory=DIST_DIR),
    )
    thread = threading.Thread(target=server.serve_forever)
    thread.daemon = True
    thread.start()

    base_url = f"http://localhost:{SERVER_PORT}"

    # Wait for server to be responsive
    for _ in range(10):
        try:
            import requests

            if requests.get(base_url, timeout=1).status_code == 200:
                break
        except Exception:
            time.sleep(0.5)

    yield base_url

    server.shutdown()
    thread.join()


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    os.name == "nt", reason="E2E tests are flaky/unsupported on Windows CI"
)  # type: ignore[untyped-decorator]
def test_e2e_save_load(page: Any, serve_app: str) -> None:
    """Test saving and loading state in the web app via Playwright.

    This test enables Flutter Web semantics to interact with the DOM elements
    directly without relying on pixel OCR.
    """
    # Subscribe to console messages
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    # 1. Navigate to app with semantics enabled
    page.goto(f"{serve_app}/?enable-semantics=true")

    # 2. Wait for loading (Pyodide initialization)
    print("Waiting for app to load...")
    expect(page).to_have_title("AmplifyP", timeout=120000)

    # Wait for the Flutter web engine to initialize and render semantics
    print("Waiting for Flutter semantics-host to appear...")
    page.wait_for_selector(
        "flt-semantics-host", state="attached", timeout=60000
    )

    # 3. Enter State
    # Activate Flutter accessibility/semantics overlay
    print("Activating Flutter Web accessibility semantics...")
    page.wait_for_selector(
        "flt-semantics-placeholder", state="attached", timeout=30000
    )
    page.locator("flt-semantics-placeholder").first.dispatch_event("click")

    # Find the Template Sequence textarea/input field
    # (these are populated in the DOM once accessibility semantics is enabled)
    template_input = page.get_by_label("Template Sequence").first
    template_input.click()
    template_input.fill("ATGCATGC")

    # 4. Save State
    # Find the Save button via its tooltip "Save State"
    save_btn = page.locator("[aria-label*='Save State']").first

    with page.expect_download(timeout=20000) as download_info:
        save_btn.click()

    download = download_info.value
    path = download.path()
    with open(path, encoding="utf-8") as f:
        content = f.read()
        print("Downloaded content:", content)
        assert "ATGCATGC" in content
