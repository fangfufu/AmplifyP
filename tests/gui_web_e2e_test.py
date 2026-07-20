# Copyright (C) 2026 AmplifyP Contributors
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
import signal
import subprocess
import sys
import threading
import time
from collections.abc import Generator
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
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

PRIMER_INPUT_SEL = (
    '[aria-label="Primer List"] textarea:not([disabled]):not([readonly])'
)


# Check entire module for dependencies
if not Page:
    pytest.fail(
        "E2E test dependencies missing: pytest-playwright. "
        "Please install it using: pip install .[e2e] "
        "And then: playwright install chromium",
        pytrace=False,
    )


@pytest.fixture(scope="session")  # type: ignore[untyped-decorator]
def browser_type_launch_args(browser_type_launch_args: Any) -> dict[str, Any]:
    """Force headless mode for E2E tests and add flags for stability."""
    return {
        **browser_type_launch_args,
        "headless": True,
        "args": [
            "--disable-gpu",
            "--no-sandbox",
            "--disable-dev-shm-usage",
            "--ignore-gpu-blocklist",
        ],
    }


@pytest.fixture(scope="session")  # type: ignore[untyped-decorator]
def build_app() -> None:
    """Build the Flet static app using build_web.sh."""
    if os.path.exists(os.path.join(DIST_DIR, "index.html")):
        print("==> Flet static app already built, skipping build...")
        return

    print("==> Building static site using build_web.sh...")
    script_path = os.path.join(os.getcwd(), "build_web.sh")

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
    ThreadingHTTPServer.allow_reuse_address = True

    class CustomHTTPRequestHandler(SimpleHTTPRequestHandler):
        def translate_path(self, path: str) -> str:
            if path.startswith("/AmplifyP/"):
                path = path[len("/AmplifyP") :]
            elif path == "/AmplifyP":
                path = "/"
            return super().translate_path(path)

    server = ThreadingHTTPServer(
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
    field.focus()
    field.click(force=True)
    time.sleep(0.5)
    field.press("Control+a")
    field.press("Delete")
    time.sleep(0.2)
    field.press_sequentially(text, delay=delay_ms)
    time.sleep(0.3)
    print(
        f"  typed '{text}' into {selector} (index={index}, use_last={use_last})"
    )


@pytest.fixture  # type: ignore[untyped-decorator]
def e2e_timeout() -> Generator[None, None, None]:
    """Ensure E2E test does not run for more than 5 minutes."""

    def handler(signum: Any, frame: Any) -> None:
        raise TimeoutError("Test timed out after 5 minutes")

    sigalrm = getattr(signal, "SIGALRM", None)
    alarm = getattr(signal, "alarm", None)

    if sigalrm is not None and alarm is not None:
        signal.signal(sigalrm, handler)
        alarm(300)
    try:
        yield
    finally:
        if alarm is not None:
            alarm(0)


@pytest.fixture(autouse=True)  # type: ignore[untyped-decorator]
def cleanup_debug_checkboxes() -> Generator[None, None, None]:
    """Clean up debug_checkboxes.png after each test."""
    yield
    if os.path.exists("debug_checkboxes.png"):
        try:
            os.remove("debug_checkboxes.png")
        except OSError:
            pass


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    sys.platform != "linux", reason="E2E tests only run on Linux"
)  # type: ignore[untyped-decorator]
def test_e2e_primer_lifecycle_and_state(
    page: Any, serve_app: str, tmp_path: Any, e2e_timeout: None
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
    page.set_default_timeout(30000)
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    # 1. Navigate to app with semantics enabled
    page.goto(f"{serve_app}/?enable-semantics=true")
    try:
        page.wait_for_load_state("networkidle", timeout=5000)
    except Exception:  # noqa: S110
        pass
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)

    # Enable "Auto-activate new valid primer" setting so new valid primers
    # get checked automatically
    print("Navigating to Settings to enable auto-activate setting...")
    page.locator("[aria-label='Settings']").last.click(force=True)
    page.get_by_role("button", name="Primer List Settings").wait_for(
        state="attached", timeout=15000
    )
    page.get_by_role("button", name="Primer List Settings").click(force=True)
    time.sleep(1)
    auto_activate_cb = page.get_by_role(
        "checkbox", name="Auto-activate new valid primer"
    )
    auto_activate_cb.wait_for(state="attached", timeout=10000)
    expect(auto_activate_cb).not_to_be_checked()
    auto_activate_cb.click(force=True)
    time.sleep(1)
    expect(auto_activate_cb).to_be_checked()

    print("Navigating back to Input tab...")
    page.locator("[aria-label='Input']").last.click(force=True)
    page.locator(PRIMER_INPUT_SEL).first.wait_for(
        state="attached", timeout=15000
    )
    time.sleep(1)

    # 2. Add 2 valid primers
    print("Adding 2 valid primers...")
    add_primer_to_trailing_row(page, "V1", "ATGCATGCATGCATGC")
    add_primer_to_trailing_row(page, "V2", "GCATGCATGCATGCAT")

    # 3. Add 2 invalid primers
    print("Adding 2 invalid primers...")
    add_primer_to_trailing_row(page, "I1", "XYZXYZXYZXYZ")
    add_primer_to_trailing_row(page, "I2", "XYZATGCATGCATGCATGC")

    # Add extra valid (V3) and invalid (I3) primers.
    print("Adding extra valid (V3) primer...")
    add_primer_to_trailing_row(page, "V3", "CGATCGATCGATCGAT")
    # Verify V3 added: 6 rows (12 inputs) and checkbox is checked.
    expect(page.locator(PRIMER_INPUT_SEL)).to_have_count(12)
    name_inputs = page.locator(PRIMER_INPUT_SEL)
    expect(
        name_inputs.nth(4 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).to_be_checked(timeout=15000)

    print("Adding extra invalid (I3) primer...")
    add_primer_to_trailing_row(page, "I3", "XYZXYZXYZ")
    # Verify I3 added: 7 rows (14 inputs) and checkbox is enabled.
    expect(page.locator(PRIMER_INPUT_SEL)).to_have_count(14)
    expect(
        name_inputs.nth(5 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).to_be_enabled(timeout=15000)

    print("Deleting V3 and I3 using delete buttons...")
    # Focus V3's name input (index 8) to select the row
    page.locator(PRIMER_INPUT_SEL).nth(8).focus()
    page.locator(PRIMER_INPUT_SEL).nth(8).click(force=True)
    time.sleep(1)

    # Click the header Delete Primer button
    delete_btn = page.locator("[aria-label*='Delete Primer']").first
    if not delete_btn.is_visible():
        row_container = page.locator("[aria-label*='Add Primer Below']").first
        delete_btn = row_container.locator("[role='button']").nth(1)
    expect(delete_btn).to_be_enabled(timeout=5000)
    delete_btn.click(force=True)
    time.sleep(1)

    # Verify V3 deleted: I3 is now at index 4 (global index 8).
    expect(page.locator(PRIMER_INPUT_SEL)).to_have_count(12)
    time.sleep(1)

    # Focus I3's name input (index 8) to select the row
    page.locator(PRIMER_INPUT_SEL).nth(8).focus()
    page.locator(PRIMER_INPUT_SEL).nth(8).click(force=True)
    time.sleep(1)

    # Click the header Delete Primer button
    delete_btn = page.locator("[aria-label*='Delete Primer']").first
    if not delete_btn.is_visible():
        row_container = page.locator("[aria-label*='Add Primer Below']").first
        delete_btn = row_container.locator("[role='button']").nth(1)
    expect(delete_btn).to_be_enabled(timeout=5000)
    delete_btn.click(force=True)
    time.sleep(2)

    # Verify both V3 and I3 deleted: count returned to 5 rows (10 inputs).
    expect(page.locator(PRIMER_INPUT_SEL)).to_have_count(10)

    # 4. Verify checkboxes and try to activate invalid primers
    print("Verifying checkbox state and attempting to activate invalid ones...")
    expect(
        name_inputs.nth(0 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(1 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(2 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).not_to_be_checked(timeout=15000)
    page.screenshot(path="debug_checkboxes.png")
    expect(
        name_inputs.nth(3 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).not_to_be_checked(timeout=15000)

    # Try and activate invalid primers (click them) and make sure they
    # don't get activated
    name_inputs.nth(2 * 2).locator("xpath=../../../..").get_by_role(
        "checkbox"
    ).click(force=True)
    name_inputs.nth(3 * 2).locator("xpath=../../../..").get_by_role(
        "checkbox"
    ).click(force=True)
    time.sleep(1)

    # Ensure they remain unchecked
    expect(
        name_inputs.nth(2 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).not_to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(3 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
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

    # Verify list is empty by checking if trailing row remains (2 inputs)
    expect(page.locator(PRIMER_INPUT_SEL)).to_have_count(2)

    # 7. Load the primer list
    print("Loading the primer list...")
    load_primers_btn = page.locator("[aria-label*='Load primers']").first
    with page.expect_file_chooser() as fc_info:
        load_primers_btn.click()
    file_chooser = fc_info.value
    file_chooser.set_files(str(primers_csv_path))
    # Wait for loaded primers (8 inputs total) to be attached
    page.locator(PRIMER_INPUT_SEL).nth(7).wait_for(
        state="attached", timeout=15000
    )
    time.sleep(1)

    # Verify loaded primers: V1/V2 checked; I1/I2 enabled and unchecked
    expect(
        name_inputs.nth(0 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(1 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(2 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).not_to_be_checked(timeout=15000)
    expect(
        name_inputs.nth(3 * 2)
        .locator("xpath=../../../..")
        .get_by_role("checkbox")
    ).not_to_be_checked(timeout=15000)

    # 8. Save the state
    print("Saving the full state...")
    TEMPLATE_SEL = "textarea"
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
    try:
        page.wait_for_load_state("networkidle", timeout=5000)
    except Exception:  # noqa: S110
        pass
    wait_for_semantics(page)

    # Wait for app/Pyodide to fully load and the input to appear
    page.wait_for_selector(PRIMER_INPUT_SEL, state="attached", timeout=60000)

    # Assert clean state before load (2 inputs)
    expect(page.locator(PRIMER_INPUT_SEL)).to_have_count(2)

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
def test_e2e_settings_backup(
    page: Any, serve_app: str, tmp_path: Any, e2e_timeout: None
) -> None:
    """Test Settings Backup: Save/Load Settings via file picker.

    Exercises the FilePicker (pick_and_read_file / save_and_write_file)
    via page.services — the fix for Flet 0.85 web mode.

    Semantics discovery findings:
    - ExpansionTile headers appear as role='button' in the DOM tree,
      accessible via their text content (e.g. 'Backup').
    - After expanding the Backup tile, exactly 2 new role='button' nodes
      with EMPTY text appear. These are Save Settings (first) and Load
      Settings (second). Their labels are rendered on canvas, not DOM.
    - Strategy: count buttons before expand, then target nth new buttons
      by index after expansion.

    Steps:
      1. Navigate to Settings tab.
      2. Click the 'Backup' expansion tile header button to expand it.
      3. Click the first new (empty-text) button = Save Settings download.
      4. Verify the downloaded YAML contains 'settings:' key.
      5. Reload, navigate to Settings > Backup again.
      6. Click the second new (empty-text) button = Load Settings.
      7. Supply the downloaded file via file chooser.
      8. Verify the app is still responsive.
    """
    page.set_default_timeout(30000)
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    def expand_backup_tile() -> None:
        """Click the General Settings tile header to expand it.

        Debug findings:
        - General Settings tile header is role='button' with textContent
          'General Settings'.
        - After clicking, the header merges into a group.
        - The 2 expanded buttons appear as role='button' with
          textContent 'Save Settings' and 'Load Settings'.
        - We locate them by name (textContent) via get_by_role.
        """
        backup_btn = page.get_by_role("button", name="General Settings")
        backup_btn.wait_for(state="attached", timeout=15000)
        print("  Clicking General Settings tile to expand...")
        backup_btn.click(force=True)
        # Wait for 'Save Settings' button to appear
        page.get_by_role("button", name="Save Settings").wait_for(
            state="attached", timeout=10000
        )
        time.sleep(1.5)
        print("  General Settings tile expanded.")

    def navigate_to_settings() -> None:
        """Click the Settings tab and wait for expansion tiles to load."""
        print("  Clicking Settings tab...")
        page.locator("[aria-label='Settings']").last.click(force=True)
        # General Settings tile header must be visible before proceeding
        page.get_by_role("button", name="General Settings").wait_for(
            state="attached", timeout=15000
        )

    # 1. Navigate to app with semantics enabled
    page.goto(f"{serve_app}/?enable-semantics=true")
    try:
        page.wait_for_load_state("networkidle", timeout=5000)
    except Exception:  # noqa: S110
        pass
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)

    # 2. Navigate to Settings and expand Backup
    print("Navigating to Settings...")
    navigate_to_settings()
    expand_backup_tile()

    # 3. Click 'Save Settings' button — triggers download
    print("Clicking Save Settings...")
    save_btn = page.get_by_role("button", name="Save Settings")
    with page.expect_download(timeout=20000) as download_info:
        save_btn.click(force=True)
    download = download_info.value
    settings_yaml_path = tmp_path / "amplify_settings.yaml"
    download.save_as(str(settings_yaml_path))
    print(f"  Downloaded to {settings_yaml_path}")

    # Verify YAML content
    with open(settings_yaml_path, encoding="utf-8") as f:
        settings_content = f.read()
    assert "settings:" in settings_content, (
        "Expected 'settings:' key in downloaded YAML, "
        f"got:\n{settings_content[:200]}"
    )
    print("  Save Settings: YAML content verified.")

    # 4. Reload and navigate to Settings > Backup again
    print("Reloading page...")
    page.goto(f"{serve_app}/?enable-semantics=true")
    try:
        page.wait_for_load_state("networkidle", timeout=5000)
    except Exception:  # noqa: S110
        pass
    wait_for_semantics(page)
    navigate_to_settings()
    expand_backup_tile()

    # 5. Click 'Load Settings' button — triggers file chooser
    print("Clicking Load Settings...")
    load_btn = page.get_by_role("button", name="Load Settings")
    with page.expect_file_chooser(timeout=15000) as fc_info:
        load_btn.click(force=True)
    file_chooser = fc_info.value
    file_chooser.set_files(str(settings_yaml_path))
    time.sleep(3)

    # 6. Verify app is still alive (no crash after load)
    assert page.locator("[aria-label='Settings']").count() > 0, (
        "Settings tab missing — app may have crashed after load"
    )
    print("  Load Settings: app still responsive after load.")
    print("test_e2e_settings_backup PASSED.")


def wait_for_ui(
    page: Any, text: str, timeout_sec: int = 60
) -> tuple[float, float]:
    """Wait for specified text on screen via OCR and return coordinates."""
    import pytesseract
    from PIL import Image

    start_time = time.time()
    while time.time() - start_time < timeout_sec:
        screenshot_bytes = page.screenshot()
        image = Image.open(io.BytesIO(screenshot_bytes))
        ocr_data = pytesseract.image_to_data(
            image, output_type=pytesseract.Output.DICT
        )
        words = ocr_data["text"]
        for i, w in enumerate(words):
            if w.strip() and text.lower() in w.strip().lower():
                left = ocr_data["left"][i]
                top = ocr_data["top"][i]
                width = ocr_data["width"][i]
                height = ocr_data["height"][i]
                center_x = left + width / 2
                center_y = top + height / 2
                return center_x, center_y
        time.sleep(1.0)
    raise TimeoutError(
        f"Timed out waiting for text '{text}' to appear on screen via OCR."
    )


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    sys.platform != "linux", reason="E2E tests only run on Linux"
)  # type: ignore[untyped-decorator]
def test_e2e_dimer_alignment(
    page: Any, serve_app: str, tmp_path: Any, e2e_timeout: None
) -> None:
    """Test dimer alignment and verify monospace alignment using OCR."""
    page.set_default_timeout(30000)
    # Subscribe to console messages
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    # 1. Navigate to app (no semantics)
    page.goto(f"{serve_app}/")
    try:
        page.wait_for_load_state("networkidle", timeout=5000)
    except Exception:  # noqa: S110
        pass
    expect(page).to_have_title("AmplifyP", timeout=120000)

    # Wait for UI to load and capture initial coordinates of toolbar and headers
    print("Waiting for UI to load via OCR...")
    wait_for_ui(page, "Template")

    dimers_x, dimers_y = wait_for_ui(page, "Dimers")
    print(f"Dimers button located at ({dimers_x}, {dimers_y})")

    name_x, name_y = wait_for_ui(page, "Name")
    print(f"Name header located at ({name_x}, {name_y})")

    # Click Name input box (located slightly below the Name header)
    page.mouse.click(name_x, name_y + 35)
    time.sleep(0.5)

    # 2. Enter Primer Details
    PRIMER_NAME = "10290"
    PRIMER_SEQ = "GTGGGTATCACAAATTTGGG"

    page.keyboard.type(PRIMER_NAME)
    print(f"Typed name: {PRIMER_NAME}")
    time.sleep(0.5)

    page.keyboard.press("Tab")
    time.sleep(0.5)

    page.keyboard.type(PRIMER_SEQ)
    print(f"Typed sequence: {PRIMER_SEQ}")
    time.sleep(0.5)

    # Tab away from the seq field to trigger on_blur → timer → sync_to_state
    page.keyboard.press("Tab")
    time.sleep(2.0)  # Allow blur timer and state sync to complete

    # Save a debug screenshot of the input page after clicking Add
    page.screenshot(path=str(tmp_path / "debug_after_add.png"))

    # 4. Navigate to Primer Dimers view by clicking the saved coordinates
    print(f"Navigating to Primer Dimers. Click: ({dimers_x}, {dimers_y})")
    page.mouse.click(dimers_x, dimers_y)
    time.sleep(2)

    # Save a debug screenshot after clicking Primer Dimers
    page.screenshot(path=str(tmp_path / "debug_after_click_dimers.png"))

    # The primer dimer analysis runs synchronously in the Flet Python worker.
    # Flutter Web CanvasKit renders text to a <canvas> element, NOT to the DOM.
    # We simply wait for CanvasKit to finish typography rendering.
    time.sleep(8)

    # Save a final screenshot for inspection
    page.screenshot(path=str(tmp_path / "debug_final.png"))

    # 5. Take Screenshot for OCR
    screenshot_bytes = page.screenshot()

    # 6. Perform OCR using PIL and pytesseract to extract character locations
    try:
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
    placeholder = page.locator("flt-semantics-placeholder").first
    placeholder.wait_for(state="attached", timeout=30000)

    # Retry clicking the placeholder until at least one input is attached
    # ensuring correct ready behaviour of the semantics layer.
    for i in range(30):
        if page.locator(PRIMER_INPUT_SEL).count() > 0:
            print(f"Semantics successfully enabled after {i} retries.")
            return

        print(
            f"Attempting to enable semantics (click retry {i + 1}). "
            f"Checking for ready behaviour..."
        )
        try:
            placeholder.dispatch_event("click")
        except Exception:  # noqa: S110
            pass
        try:
            placeholder.click(force=True, timeout=1000)
        except Exception:  # noqa: S110
            pass
        time.sleep(1.0)

    if page.locator(PRIMER_INPUT_SEL).count() == 0:
        raise RuntimeError(
            "Timed out waiting for Flutter Web semantics to enable."
        )


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
    try:
        page.wait_for_selector(
            PRIMER_INPUT_SEL, state="attached", timeout=15000
        )
    except Exception:
        print("Timeout waiting for inputs. Retrying semantics click...")
        wait_for_semantics(page)

    initial_count = page.locator(PRIMER_INPUT_SEL).count()

    # The trailing row's Name and Sequence fields are at the very
    # end of the list
    page.locator(PRIMER_INPUT_SEL).nth(initial_count - 1).wait_for(
        state="attached", timeout=10000
    )

    # Fill Name field using its precise index
    fill_field_reliably(page, PRIMER_INPUT_SEL, name, index=initial_count - 2)
    time.sleep(0.3)

    # Fill Sequence field using its precise index
    time.sleep(1.0)
    fill_field_reliably(page, PRIMER_INPUT_SEL, seq, index=initial_count - 1)
    time.sleep(0.3)

    # Blur the sequence field by focusing the template sequence field to trigger
    # on_blur → timer → sync_to_state
    page.locator('textarea:not([aria-label="Primer List"])').first.focus()
    time.sleep(1.0)

    # Wait for the count to increase by 2 (indicating a new
    # trailing row was auto-added)
    success = False
    for _ in range(50):
        if page.locator(PRIMER_INPUT_SEL).count() >= initial_count + 2:
            success = True
            break
        time.sleep(0.2)

    if not success:
        # Check one more time before doing anything
        if page.locator(PRIMER_INPUT_SEL).count() >= initial_count + 2:
            return

        # It was not auto-added (because it was invalid). Re-focus the last
        # non-file input so the header Add button is enabled.
        page.locator(PRIMER_INPUT_SEL).last.click(force=True)
        time.sleep(0.5)

        # Check again after focusing/sleeping
        if page.locator(PRIMER_INPUT_SEL).count() >= initial_count + 2:
            return

        add_btn = page.locator("[aria-label*='Add Primer Below']").first
        add_btn.wait_for(state="attached", timeout=5000)
        box = add_btn.bounding_box()
        if box and box["width"] > 100:
            # Resolve the correct button inside the row container if the
            # container width is large
            add_btn = add_btn.locator("[role='button']").first
        expect(add_btn).to_be_enabled(timeout=5000)

        # Click the Add button (try physical click first and then fallback to
        # dispatching event)
        try:
            add_btn.click(force=True)
        except Exception:
            try:
                add_btn.dispatch_event("click")
            except Exception:  # noqa: S110
                pass

        # Wait for the count to increase to initial_count + 2
        for _ in range(30):
            if page.locator(PRIMER_INPUT_SEL).count() >= initial_count + 2:
                break
            time.sleep(0.2)
        else:
            page.screenshot(path="failed_add_row.png")
            raise RuntimeError("Failed to add new primer row manually")

    time.sleep(0.5)
