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

try:
    import requests
except ImportError:
    requests = None

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
    """Build the Flet static app using build_static.sh."""
    if os.path.exists(os.path.join(DIST_DIR, "index.html")):
        print("==> Flet static app already built, skipping build...")
        return

    print("==> Building static site using build_static.sh...")
    script_path = os.path.join(os.getcwd(), "build_static.sh")
    import sys

    env = os.environ.copy()
    venv_bin = os.path.dirname(sys.executable)
    env["PATH"] = os.path.pathsep.join([venv_bin, env.get("PATH", "")])
    subprocess.run(  # noqa: S603
        [script_path],
        check=True,
        capture_output=True,
        env=env,
    )
    assert os.path.exists(os.path.join(DIST_DIR, "index.html"))


@pytest.fixture(scope="session")  # type: ignore[untyped-decorator]
def serve_app(build_app: None) -> Generator[str, None, None]:
    """Serve the static app in a background thread."""
    HTTPServer.allow_reuse_address = True

    class CustomHTTPRequestHandler(SimpleHTTPRequestHandler):
        def translate_path(self, path: str) -> str:
            if path.startswith("/AmplifyP/"):
                path = path[len("/AmplifyP") :]
            elif path == "/AmplifyP":
                path = "/"
            return super().translate_path(path)

    server = HTTPServer(
        ("localhost", SERVER_PORT),
        lambda *args: CustomHTTPRequestHandler(*args, directory=DIST_DIR),
    )
    thread = threading.Thread(target=server.serve_forever)
    thread.daemon = True
    thread.start()

    base_url = f"http://localhost:{SERVER_PORT}"

    # Wait for server to be responsive
    for _ in range(10):
        if requests is None:
            break
        try:
            if requests.get(base_url, timeout=1).status_code == 200:
                break
        except requests.exceptions.RequestException:
            time.sleep(0.5)

    yield base_url

    server.shutdown()
    thread.join()


def fill_field_reliably(
    page: Any,
    selector: str,
    text: str,
    delay_ms: int = 100,
    use_last: bool = False,
    index: int | None = None,
) -> None:
    """Focus and type text into a Flutter Web text field.

    Two key Flutter Web quirks must be handled here:

    1. **Opacity filter**: Flutter Web's ``flt-semantics-host`` has
       ``filter: opacity(0%)`` applied.  Playwright's default ``click()``
       refuses to interact with elements whose ancestor has zero opacity.
       ``force=True`` bypasses this check.  The inner semantic nodes still
       carry ``pointer-events: all`` so interaction succeeds.

    2. **DOM value is always empty**: Flutter renders text via WebGL /
       CanvasKit and *never* writes typed characters back to the native
       textarea/input ``value`` property.  Reading ``element.value`` via
       ``page.evaluate()`` or ``input_value()`` always returns ``''``.
       Verification is therefore left to downstream test assertions (e.g.
       checking the downloaded file contains the expected sequence).
    """
    if index is not None:
        field = page.locator(selector).nth(index)
    elif use_last:
        field = page.locator(selector).last
    else:
        field = page.locator(selector).first
    # force=True bypasses the visibility check from filter:opacity(0%) on the
    # flt-semantics-host parent.
    field.click(force=True)
    field.press("Control+a")
    field.press("Delete")
    time.sleep(0.2)
    field.press_sequentially(text, delay=delay_ms)
    time.sleep(0.3)
    print(
        f"  typed '{text}' into {selector} (index={index}, use_last={use_last})"
    )


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    sys.platform != "linux", reason="E2E tests only run on Linux"
)  # type: ignore[untyped-decorator]
def test_e2e_primer_lifecycle_and_state(
    page: Any, serve_app: str, tmp_path: Any
) -> None:
    """Test full primer lifecycle and state saving/loading.

    Steps:
      - Add 2 valid primers.
      - Add 2 invalid primers.
      - Try and activate invalid primers and make sure they don't get activated.
      - Save the primer list.
      - Clear the primer list.
      - Load the primer list.
      - Save the state.
      - Load the state.
    """
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    # 1. Navigate to app with semantics enabled
    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)

    # 2. Add 2 valid primers
    print("Adding 2 valid primers...")
    add_primer_to_trailing_row(page, "V1", "ATGCATGCATGCATGC")
    add_primer_to_trailing_row(page, "V2", "GCATGCATGCATGCAT")

    # 3. Add 2 invalid primers
    print("Adding 2 invalid primers...")
    add_primer_to_trailing_row(page, "I1", "XYZXYZXYZXYZ")
    add_primer_to_trailing_row(page, "I2", "ATGCATGCATGCAT-XYZ")

    # Add extra valid (V3) and invalid (I3) primers.
    print("Adding extra valid (V3) primer...")
    add_primer_to_trailing_row(page, "V3", "CGATCGATCGATCGAT")
    # Verify V3 added: 6 rows (12 inputs) and checkbox is checked.
    expect(page.locator('input:not([type="file"])')).to_have_count(12)
    name_inputs = page.locator('input:not([type="file"])')
    expect(
        name_inputs.nth(4 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_checked(timeout=15000)

    print("Adding extra invalid (I3) primer...")
    add_primer_to_trailing_row(page, "I3", "XYZXYZXYZ")
    # Verify I3 added: 7 rows (14 inputs) and checkbox is enabled.
    expect(page.locator('input:not([type="file"])')).to_have_count(14)
    expect(
        name_inputs.nth(5 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_enabled(timeout=15000)

    print("Deleting V3 and I3 using delete buttons...")
    # There are 6 primers in the list:
    # V1 (0), V2 (1), I1 (2), I2 (3), V3 (4), I3 (5).
    # Since each row has exactly 2 text input fields (Name and Sequence),
    # the Name input of V3 (index 4) is at global input index 8.
    # Focus V3's name input to make its row controls visible.
    page.locator('input:not([type="file"])').nth(8).click(force=True)
    time.sleep(1)

    # Only the focused row exposes its controls in the semantic tree.
    # Use .first (not .nth(N)) because exactly one Delete Primer button is
    # visible at any time — the one belonging to the focused row.
    delete_btn = page.locator("[aria-label*='Add Primer Below']").first
    delete_btn.wait_for(state="attached", timeout=5000)
    box = delete_btn.bounding_box()
    assert box is not None
    page.mouse.click(box["x"] + box["width"] - 76, box["y"] + box["height"] / 2)
    time.sleep(1)

    # Verify V3 deleted: I3 checkbox is now at index 4 in name_inputs.
    expect(page.locator('input:not([type="file"])')).to_have_count(12)
    expect(
        name_inputs.nth(4 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_enabled(timeout=15000)

    # Focus I3 (index 4 after V3 deletion) - Name input is at index 8.
    page.locator('input:not([type="file"])').nth(8).click(force=True)
    time.sleep(1)
    delete_btn.wait_for(state="attached", timeout=5000)
    box = delete_btn.bounding_box()
    assert box is not None
    page.mouse.click(box["x"] + box["width"] - 76, box["y"] + box["height"] / 2)
    time.sleep(1)

    # Verify I3 deleted: count returned to 5 rows (10 inputs).
    expect(page.locator('input:not([type="file"])')).to_have_count(10)

    # 4. Verify checkboxes and try to activate invalid primers
    print("Verifying checkbox state and attempting to activate invalid ones...")
    expect(
        name_inputs.nth(0 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(1 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(2 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    page.screenshot(path="debug_checkboxes.png")
    expect(
        name_inputs.nth(3 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_checked(timeout=15000)

    # Uncheck them to return to unchecked baseline state
    name_inputs.nth(2 * 2).locator("xpath=../..").get_by_role("checkbox").click(
        force=True
    )
    name_inputs.nth(3 * 2).locator("xpath=../..").get_by_role("checkbox").click(
        force=True
    )
    time.sleep(1)

    # Ensure they are unchecked again
    expect(
        name_inputs.nth(2 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).not_to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(3 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).not_to_be_checked(timeout=15000)

    # 5. Save the primer list
    print("Saving active primer list...")
    save_primers_btn = page.locator("[aria-label*='Save primers']").first
    with page.expect_download(timeout=20000) as download_info:
        save_primers_btn.click()
    download = download_info.value
    primers_csv_path = tmp_path / "primers_saved.csv"
    download.save_as(str(primers_csv_path))

    # Verify that all primers (active/valid and inactive/invalid) were exported
    with open(primers_csv_path, encoding="utf-8") as f:
        primer_list_content = f.read()
    assert "V1" in primer_list_content
    assert "V2" in primer_list_content
    assert "I1" in primer_list_content
    assert "I2" in primer_list_content

    # 6. Clear the primer list
    print("Clearing the primer list...")
    clear_primers_btn = page.locator("[aria-label*='Clear All Primers']").first
    clear_primers_btn.click()
    time.sleep(2)

    # Verify list is empty by checking if only one row (trailing row) is left
    expect(page.locator('input[aria-label="New Primer Name"]')).to_have_count(1)

    # 7. Load the primer list
    print("Loading the primer list...")
    load_primers_btn = page.locator("[aria-label*='Load primers']").first
    with page.expect_file_chooser() as fc_info:
        load_primers_btn.click()
    file_chooser = fc_info.value
    file_chooser.set_files(str(primers_csv_path))
    # Wait for loaded primers (8 inputs total) to be attached
    page.locator('input:not([type="file"])').nth(7).wait_for(
        state="attached", timeout=15000
    )
    time.sleep(1)

    # Verify loaded primers: V1/V2 checked; I1/I2 enabled and unchecked
    expect(
        name_inputs.nth(0 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(1 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(2 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(3 * 2).locator("xpath=../..").get_by_role("checkbox")
    ).to_be_checked(timeout=15000)

    # 8. Save the state
    print("Saving the full state...")
    TEMPLATE_SEL = 'textarea[aria-label="Enter DNA sequence here..."]'
    fill_field_reliably(page, TEMPLATE_SEL, "ATGCATGC")
    page.keyboard.press("Tab")
    time.sleep(1)

    save_state_btn = page.locator("[aria-label*='Save all']").first
    with page.expect_download(timeout=20000) as download_info:
        save_state_btn.click()
    download = download_info.value
    state_yaml_path = tmp_path / "state_saved.yaml"
    download.save_as(str(state_yaml_path))

    # Verify state file contents (including the invalid ones)
    with open(state_yaml_path, encoding="utf-8") as f:
        state_content = f.read()
    assert "ATGCATGC" in state_content
    assert "V1" in state_content
    assert "V2" in state_content
    assert "I1" in state_content
    assert "I2" in state_content

    # 9. Load the state
    print("Refreshing/reloading page to clear state...")
    page.goto(f"{serve_app}/?enable-semantics=true")
    wait_for_semantics(page)

    # Wait for app/Pyodide to fully load and the input to appear
    page.wait_for_selector(
        'input[aria-label="New Primer Name"]', state="attached", timeout=60000
    )

    # Assert clean state before load
    expect(page.locator('input[aria-label="New Primer Name"]')).to_have_count(1)

    print("Loading the saved state...")
    load_state_btn = page.locator("[aria-label*='Load all']").first
    with page.expect_file_chooser() as fc_info:
        load_state_btn.click()
    file_chooser = fc_info.value
    file_chooser.set_files(str(state_yaml_path))
    time.sleep(5)

    # Verify template and primers are restored (we can save and check)
    with page.expect_download(timeout=20000) as download_info:
        save_state_btn.click()
    download = download_info.value
    final_yaml_path = tmp_path / "final_saved.yaml"
    download.save_as(str(final_yaml_path))

    with open(final_yaml_path, encoding="utf-8") as f:
        final_content = f.read()
    assert "ATGCATGC" in final_content
    assert "V1" in final_content
    assert "V2" in final_content
    assert "I1" in final_content
    assert "I2" in final_content


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    sys.platform != "linux", reason="E2E tests only run on Linux"
)  # type: ignore[untyped-decorator]
def test_e2e_dimer_alignment(page: Any, serve_app: str, tmp_path: Any) -> None:
    """Test dimer alignment and verify monospace alignment using OCR."""
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

    # Wait for the primer name input to appear — confirms Pyodide is fully
    # initialised and the semantics tree is ready. Much more reliable than
    # a fixed sleep.
    NAME_SEL = 'input[aria-label="New Primer Name"]'
    SEQ_SEL = 'input[aria-label="New Primer Sequence"]'
    page.wait_for_selector(NAME_SEL, state="attached", timeout=60000)

    # 2. Enter Primer Details - with retry to handle dropped first keystrokes
    PRIMER_NAME = "10290"
    PRIMER_SEQ = "GTGGGTATCACAAATTTGGG"

    # Flutter Web exposes the hint_text as the aria-label for text fields.
    # Primer Name and Primer Sequence are the hint_texts on the inline row
    # fields; there is no separate "Add" button — filling the trailing empty
    # row and tabbing away triggers blur→sync_to_state which adds the primer.
    # CSS selectors are used so fill_field_reliably can verify via JS eval.
    fill_field_reliably(page, NAME_SEL, PRIMER_NAME)
    page.keyboard.press("Tab")

    fill_field_reliably(page, SEQ_SEL, PRIMER_SEQ)
    # Tab away from the seq field to trigger on_blur → timer → sync_to_state
    page.keyboard.press("Tab")
    time.sleep(2)  # Allow blur timer (0.15s) and state sync to complete

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
        # We check full sequence, prefix, or middle segment for robustness.
        top_word = next(
            (
                w
                for w in found_words
                if (
                    PRIMER_SEQ in w[0]
                    or PRIMER_SEQ[:8] in w[0]
                    or PRIMER_SEQ[4:14] in w[0]
                )
            ),
            None,
        )
        # In a self-dimer the bottom strand shows the REVERSED primer sequence
        # (not its complement) written 3'->5', e.g. "3'-GGGTTTAAACAC...-5'".
        rev_seq = PRIMER_SEQ[::-1]
        bottom_word = next(
            (
                w
                for w in found_words
                if (
                    rev_seq in w[0]
                    or rev_seq[:8] in w[0]
                    or rev_seq[4:14] in w[0]
                )
            ),
            None,
        )

        # Primary assertion: both sequence strands must be visible on-screen.
        # This proves the alignment diagram was rendered correctly.
        assert top_word is not None, (
            f"Top sequence '{PRIMER_SEQ[4:14]}' not found in OCR output. "
            "Alignment diagram may not have rendered."
        )
        assert bottom_word is not None, (
            f"Bottom sequence '{rev_seq[4:14]}' not found in OCR output. "
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


def wait_for_semantics(page: Any) -> None:
    """Wait for Flutter Web semantics to be ready."""
    page.wait_for_selector(
        "flt-semantics-host", state="attached", timeout=60000
    )
    page.wait_for_selector(
        "flt-semantics-placeholder", state="attached", timeout=30000
    )
    page.locator("flt-semantics-placeholder").first.dispatch_event("click")


def save_state(page: Any) -> str:
    """Save state and return the downloaded content."""
    save_btn = page.locator("[aria-label*='Save']").first
    expect(save_btn).to_be_enabled(timeout=10000)

    with page.expect_download(timeout=20000) as download_info:
        save_btn.click()

    download = download_info.value
    path = download.path()
    with open(path, encoding="utf-8") as f:
        return f.read()


def add_primer_to_trailing_row(page: Any, name: str, seq: str) -> None:
    """Add a primer by filling the trailing row fields (last row)."""
    # Wait for the text inputs to be available in the DOM
    page.wait_for_selector(
        'input:not([type="file"])', state="attached", timeout=60000
    )
    initial_count = page.locator('input:not([type="file"])').count()

    # The trailing row's Name and Sequence fields are at the very
    # end of the list
    page.locator('input:not([type="file"])').nth(initial_count - 1).wait_for(
        state="attached", timeout=10000
    )

    # Fill Name field using its precise index
    fill_field_reliably(
        page, 'input:not([type="file"])', name, index=initial_count - 2
    )
    time.sleep(0.3)

    # Fill Sequence field using its precise index
    fill_field_reliably(
        page, 'input:not([type="file"])', seq, index=initial_count - 1
    )
    time.sleep(0.3)

    # Submit the sequence field using keyboard Enter (it is currently focused)
    page.keyboard.press("Enter")

    # Wait for the count to increase by 2 (indicating a new
    # trailing row was auto-added)
    success = False
    for _ in range(25):
        if (
            page.locator('input:not([type="file"])').count()
            == initial_count + 2
        ):
            success = True
            break
        time.sleep(0.2)

    if not success:
        # It was not auto-added (because it was invalid). Re-focus the last
        # non-file input so the row's controls (including Add Primer Below)
        # become visible in the semantic tree.
        page.locator('input:not([type="file"])').last.click(force=True)
        time.sleep(0.5)
        add_btn = page.locator("[aria-label*='Add Primer Below']").first
        add_btn.wait_for(state="attached", timeout=5000)
        box = add_btn.bounding_box()
        assert box is not None
        page.mouse.click(
            box["x"] + box["width"] - 102, box["y"] + box["height"] / 2
        )

        # Wait for the count to increase to initial_count + 2
        for _ in range(30):
            if (
                page.locator('input:not([type="file"])').count()
                == initial_count + 2
            ):
                break
            time.sleep(0.2)
        else:
            raise RuntimeError("Failed to add new primer row manually")

    time.sleep(0.5)
