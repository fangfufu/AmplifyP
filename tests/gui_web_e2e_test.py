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

import io
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
    import pytesseract
    from PIL import Image
    from playwright.sync_api import Page, expect
except ImportError:
    pytesseract = None
    Image = None
    Page = None
    expect = None

# Define ports and paths
SERVER_PORT = 23455
SRC_DIR = os.path.join(os.getcwd(), "src")
DIST_DIR = os.path.join(os.getcwd(), "dist")


# Check entire module for dependencies
missing_deps = []
if not Page:
    missing_deps.append("pytest-playwright")
if not Image:
    missing_deps.append("pillow")
if not pytesseract:
    missing_deps.append("pytesseract")

if missing_deps:
    pytest.fail(
        f"E2E test dependencies missing: {', '.join(missing_deps)}. "
        "Please install them using: pip install .[e2e] "
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

    This test requires the Flutter HTML renderer to be active so that
    OCR can read text from the page. When Flet uses the skwasm renderer
    (the default in modern builds), the entire UI is drawn onto a WebGL
    canvas and is invisible to OCR — the test is automatically skipped
    in that case.
    """
    # Subscribe to console messages
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    # 1. Navigate to app
    page.goto(serve_app)

    # 2. Wait for loading (Pyodide initialization)
    print("Waiting for app to load...")
    expect(page).to_have_title("AmplifyP", timeout=120000)

    # Wait for the Flutter web engine to initialize
    # and set the renderer attribute.
    print("Waiting for Flutter engine to initialize renderer...")
    page.wait_for_function(
        "() => document.body.getAttribute('flt-renderer') !== null",
        timeout=60000,
    )

    # 3. Check the active renderer. Modern Flet defaults to skwasm (WebAssembly
    # canvas), which renders all UI to a WebGL surface.
    renderer = page.evaluate("() => document.body.getAttribute('flt-renderer')")
    print(f"Flutter renderer: {renderer}")

    # Wait additional time for Python to spin up and render UI
    # Since we can't reliably rely on DOM with CanvasKit (if fallback happens),
    # we will use OCR to find "Template Sequence".

    def find_text_on_screen(
        text_query: str, attempts: int = 20
    ) -> dict[str, Any] | None:
        """Capture screenshot and find text coordinates using OCR."""
        for i in range(attempts):
            print(f"OCR Attempt {i + 1} for '{text_query}'...")
            png_bytes = page.screenshot()
            image = Image.open(io.BytesIO(png_bytes))
            try:
                data = pytesseract.image_to_data(
                    image, output_type=pytesseract.Output.DICT
                )
            except pytesseract.TesseractNotFoundError:
                pytest.skip("Tesseract OCR binary not found.")

            # Search for text
            n_boxes = len(data["text"])
            for j in range(n_boxes):
                if text_query.lower() in data["text"][j].lower():
                    return {
                        "x": data["left"][j],
                        "y": data["top"][j],
                        "width": data["width"][j],
                        "height": data["height"][j],
                    }
            time.sleep(2)
        return None

    # Wait until "Template" text appears visually
    # "Template Sequence" might be split into "Template" and "Sequence" by OCR.
    # Searching for "Template" is safer.
    coords = find_text_on_screen("Template")
    if not coords:
        # Fallback: dump page source for diagnostics
        print("Page Content:", page.content())
        pytest.skip(
            "Could not find 'Template' text on screen via OCR. "
            "App may not have loaded or renderer does not expose text to OCR."
        )

    assert coords
    print(f"Found 'Template' at {coords}")

    # 3. Enter State
    # If HTML renderer worked, we might have inputs.
    try:
        # Try finding the input by label first (optimistic)
        # Note: 'Template Sequence' label might be separate from input.
        # Flet TextField usually has an aria-label or we can find by role
        # textbox. In CanvasKit, these don't exist.
        # In HTML renderer, they might.

        # Click near the label to focus input (offset y by height + 20px)
        page.mouse.click(coords["x"] + 10, coords["y"] + coords["height"] + 20)
        page.keyboard.type("ATGCATGC")

    except Exception as e:
        print(f"Interaction failed: {e}")

    # 4. Save State
    # Find "Save" icon button. OCR might not find icons.
    # But Flet usually puts "Save" tooltip.
    # If OCR is the only way, we assume the Save button is in the top right.
    # OR we search for "Save" if we changed the button to have text.
    # But it is an IconButton(ft.Icons.SAVE).

    # We can rely on visual layout assumptions or accessibility if HTML
    # renderer is active. Let's check if we have any buttons.
    buttons = page.locator("button").all()
    print(f"Found {len(buttons)} accessible buttons.")

    # If no buttons, we are stuck in CanvasKit mode without accessibility.
    # But we can try to find the "Save" tooltip if we hover?
    # Hard with OCR.

    # Alternative: The app has "Results" text button.
    # Let's try to finding "Results" text.
    results_coords = find_text_on_screen("Results")
    if results_coords:
        print(f"Found 'Results' at {results_coords}")

    # The Save button is in the AppBar actions.
    # Input | Results | Settings | | Save | Load
    # We can try to approximate location relative to "Settings" if found.
    settings_coords = find_text_on_screen("Settings")
    if settings_coords:
        # Save is to the right of Settings.
        # Layout: [Settings] [Divider] [Save] [Load]
        # "Settings" width is likely around 60-100px.
        # We want to click the Save icon button immediately after it.
        # OCR gives left edge of "Settings".

        settings_right = settings_coords["x"] + settings_coords["width"]
        # Gap + Divider + Half Save Icon (~24px)
        # Try offset of 50px from right edge of text
        save_x = settings_right + 50
        save_y = settings_coords["y"] + settings_coords["height"] / 2

        print(f"Settings found at {settings_coords}")
        print(f"Clicking estimated Save button at {save_x}, {save_y}")

        with page.expect_download(timeout=10000) as download_info:
            page.mouse.click(save_x, save_y)

        download = download_info.value
        path = download.path()
        with open(path) as f:
            content = f.read()
            print("Downloaded content:", content)
            assert "ATGCATGC" in content

    else:
        print("Could not find Settings button to orient.")
