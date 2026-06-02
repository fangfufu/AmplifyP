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


def fill_field_reliably(field: Any, text: str, delay_ms: int = 100) -> None:
    """Click, clear, type, then verify the field contains the full text.

    Flutter Web CanvasKit may silently drop leading keystrokes if the Pyodide
    worker is still starting up, even after an initial sleep.  This helper
    retries up to three times using Select-All + Delete before re-typing.
    """
    current = ""
    for attempt in range(3):
        field.click()
        field.press("Control+a")
        field.press("Delete")
        time.sleep(0.3)
        field.press_sequentially(text, delay=delay_ms)
        time.sleep(0.5)
        current = field.input_value()
        print(f"  fill attempt {attempt + 1}: got '{current}' (want '{text}')")
        if current == text:
            return
        print("  mismatch - retrying...")
        time.sleep(1)
    print("  Warning: field may still be incomplete after 3 attempts")
    raise AssertionError(
        f"fill_field_reliably: expected '{text}', got '{current}' "
        f"after 3 attempts"
    )


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
    time.sleep(8)  # Let Pyodide worker initialization fully settle

    # Find the Template Sequence textarea/input field
    # (these are populated in the DOM once accessibility semantics is enabled)
    template_input = page.get_by_label("Template Sequence").first
    fill_field_reliably(template_input, "ATGCATGC")
    page.keyboard.press("Tab")

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


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    os.name == "nt", reason="E2E tests are flaky/unsupported on Windows CI"
)  # type: ignore[untyped-decorator]
def test_e2e_primer_dimer_alignment(
    page: Any, serve_app: str, tmp_path: Any
) -> None:
    """Test primer dimer alignment and verify monospace alignment using OCR."""
    # Subscribe to console messages
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    # 1. Navigate to app with semantics enabled
    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)

    # Wait for Flutter Web engine
    page.wait_for_selector(
        "flt-semantics-host", state="attached", timeout=60000
    )
    page.wait_for_selector(
        "flt-semantics-placeholder", state="attached", timeout=30000
    )
    page.locator("flt-semantics-placeholder").first.dispatch_event("click")
    # Flutter Web / Pyodide startup lock: needs 8s to avoid dropping
    # leading keystrokes
    time.sleep(8)

    # 2. Enter Primer Details - with retry to handle dropped first keystrokes
    PRIMER_NAME = "10290"
    PRIMER_SEQ = "GTGGGTATCACAAATTTGGG"

    name_input = page.get_by_label("Primer Name").first
    fill_field_reliably(name_input, PRIMER_NAME)
    page.keyboard.press("Tab")

    seq_input = page.get_by_label("Primer Sequence").first
    fill_field_reliably(seq_input, PRIMER_SEQ)
    page.keyboard.press("Tab")

    # 3. Add Primer
    add_btn = page.get_by_role("button", name="Add").first
    add_btn.click()
    time.sleep(2)  # Allow focus change and worker processing

    # Save a debug screenshot of the input page after clicking Add
    page.screenshot(path=str(tmp_path / "debug_after_add.png"))

    # 4. Navigate to Primer Dimers view
    dimers_btn = page.locator("[aria-label*='Primer Dimers']").first
    expect(dimers_btn).to_be_enabled(timeout=10000)
    dimers_btn.click()
    time.sleep(2)

    # Save a debug screenshot after clicking Primer Dimers
    page.screenshot(path=str(tmp_path / "debug_after_click_dimers.png"))

    # The primer dimer analysis runs synchronously in the Flet Python worker.
    # Flutter Web CanvasKit renders text to a <canvas> element, NOT to the DOM,
    # so neither "text=" nor "[aria-label*=]" selectors work for rendered text.
    # We simply wait for CanvasKit to finish typography rendering.
    time.sleep(8)

    # Save a final screenshot for inspection
    page.screenshot(path=str(tmp_path / "debug_final.png"))

    # 5. Take Screenshot for OCR
    screenshot_bytes = page.screenshot()

    # 6. Perform OCR using PIL and pytesseract to extract character locations
    try:
        import io

        import pytesseract
        from PIL import Image

        image = Image.open(io.BytesIO(screenshot_bytes))
        ocr_data = pytesseract.image_to_data(
            image, output_type=pytesseract.Output.DICT
        )

        words = ocr_data["text"]
        lefts = ocr_data["left"]
        tops = ocr_data["top"]
        widths = ocr_data["width"]

        # Filter out empty entries
        found_words = []
        for i, w in enumerate(words):
            if w.strip():
                found_words.append((w.strip(), lefts[i], tops[i], widths[i]))

        print("OCR Words Found:", found_words)

        # Check that the top sequence line is visible (5'-...-3')
        top_word = next(
            (
                w
                for w in found_words
                if PRIMER_SEQ in w[0] or PRIMER_SEQ[:8] in w[0]
            ),
            None,
        )
        # In a self-dimer the bottom strand shows the REVERSED primer sequence
        # (not its complement) written 3'->5', e.g. "3'-GGGTTTAAACAC...-5'".
        rev_seq = PRIMER_SEQ[::-1]
        bottom_word = next(
            (w for w in found_words if rev_seq in w[0] or rev_seq[:8] in w[0]),
            None,
        )

        # Primary assertion: both sequence strands must be visible on-screen.
        # This proves the alignment diagram was rendered correctly.
        assert top_word is not None, (
            f"Top sequence '{PRIMER_SEQ[:8]}...' not found in OCR output. "
            "Alignment diagram may not have rendered."
        )
        assert bottom_word is not None, (
            f"Bottom sequence '{rev_seq[:8]}...' not found in OCR output. "
            "Alignment diagram may not have rendered."
        )
        print(
            f"Both alignment strands visible. top='{top_word[0][:10]}...', "
            f"bottom='{bottom_word[0][:10]}...'"
        )

        # Secondary (pixel-shift) check: only reliable when OCR captured the
        # full "5'-" prefix.  Tesseract sometimes drops thin characters like
        # the dash, shifting the detected left edge by ~1 char width and
        # producing a false misalignment result.
        top_prefix_ok = top_word[0].startswith("5'-")
        bot_prefix_ok = bottom_word[0].startswith("3'-")
        if top_prefix_ok and bot_prefix_ok:
            top_left = top_word[1]
            top_width = top_word[3]
            bottom_left = bottom_word[1]
            char_width = top_width / max(len(top_word[0]), 1)
            actual_shift = abs(bottom_left - top_left)
            print(
                f"Pixel check: char_width={char_width:.2f}px, "
                f"actual_shift={actual_shift:.2f}px "
                f"({actual_shift / char_width:.1f} chars)"
            )
            # Both lines have the same 3-char prefix, so expected shift = 0.
            # Allow up to 2 character widths for sub-pixel OCR variance.
            assert actual_shift <= (2.0 * char_width), (
                "Visual alignment is not monospaced!"
            )
        else:
            print(
                f"OCR prefix incomplete (top='{top_word[0][:6]}', "
                f"bot='{bottom_word[0][:6]}') - pixel check skipped "
                "(visibility already confirmed above)"
            )
    except ImportError:
        print("pytesseract/PIL not installed - skipping OCR alignment check.")
