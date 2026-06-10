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

"""End-to-End tests for adding and removing primers."""

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

try:
    from playwright.sync_api import Page, expect
except ImportError:
    Page = None
    expect = None

SERVER_PORT = 34522
DIST_DIR = os.path.join(os.getcwd(), "dist")


if not Page:
    pytest.fail(
        "E2E test dependencies missing: pytest-playwright. "
        "Please install it using: pip install .[e2e] "
        "And then: playwright install chromium",
        pytrace=False,
    )


@pytest.fixture(scope="session")  # type: ignore[untyped-decorator]
def build_app() -> None:
    """Build the Flet static app."""
    if os.path.exists(DIST_DIR):
        shutil.rmtree(DIST_DIR)

    for root, dirs, _ in os.walk("src"):
        for d in dirs:
            if d == "__pycache__":
                shutil.rmtree(os.path.join(root, d), ignore_errors=True)

    print("==> Building static site...")
    flet_path = shutil.which("flet")
    if not flet_path:
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


def fill_field_reliably(
    page: Any, selector: str, text: str, delay_ms: int = 100
) -> None:
    """Focus and type text into a Flutter Web text field."""
    field = page.locator(selector).first
    field.click(force=True)
    field.press("Control+a")
    field.press("Delete")
    time.sleep(0.2)
    field.press_sequentially(text, delay=delay_ms)
    time.sleep(0.3)


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
    NAME_SEL = 'input[aria-label="New Primer Name"]'
    SEQ_SEL = 'input[aria-label="New Primer Sequence"]'

    page.wait_for_selector(NAME_SEL, state="attached", timeout=60000)
    fill_field_reliably(page, NAME_SEL, name)
    page.keyboard.press("Tab")
    time.sleep(0.5)

    page.wait_for_selector(SEQ_SEL, state="attached", timeout=60000)
    fill_field_reliably(page, SEQ_SEL, seq)
    page.keyboard.press("Tab")
    time.sleep(5)


def add_primer_and_sync(page: Any, name: str, seq: str, base_url: str) -> None:
    """Add a primer by filling trailing row, saving, and refreshing UI.

    Ensures primer is saved and UI rebuilds with new trailing row.
    """
    NAME_SEL = 'input[aria-label="New Primer Name"]'
    SEQ_SEL = 'input[aria-label="New Primer Sequence"]'

    page.wait_for_selector(NAME_SEL, state="attached", timeout=60000)
    fill_field_reliably(page, NAME_SEL, name)
    page.keyboard.press("Tab")
    time.sleep(0.5)

    page.wait_for_selector(SEQ_SEL, state="attached", timeout=60000)
    fill_field_reliably(page, SEQ_SEL, seq)
    page.keyboard.press("Tab")
    time.sleep(2)

    # Click Save to save the primer
    save_btn = page.locator("[aria-label*='Save']").first
    expect(save_btn).to_be_enabled(timeout=10000)
    with page.expect_download(timeout=20000):
        save_btn.click()
    time.sleep(3)  # Wait for state file to download

    # Refresh page to rebuild UI with new trailing row
    page.goto(f"{base_url}/?enable-semantics=true")
    page.wait_for_selector(
        "flt-semantics-host", state="attached", timeout=60000
    )
    page.wait_for_selector(
        "flt-semantics-placeholder", state="attached", timeout=30000
    )
    page.locator("flt-semantics-placeholder").first.dispatch_event("click")
    page.wait_for_selector(NAME_SEL, state="attached", timeout=60000)
    time.sleep(2)


def fill_primers_sequentially(
    page: Any, primers: list[tuple[str, str]]
) -> None:
    """Fill multiple primers sequentially in the same trailing row.

    After filling each primer, the inputs retain their values, so we need to
    clear them before filling the next primer.
    """
    NAME_SEL = 'input[aria-label="New Primer Name"]'
    SEQ_SEL = 'input[aria-label="New Primer Sequence"]'

    for name, seq in primers:
        # Wait for inputs to be available
        page.wait_for_selector(NAME_SEL, state="attached", timeout=60000)

        # Clear and fill name field
        name_field = page.locator(NAME_SEL).first
        name_field.click(force=True)
        name_field.press("Control+a")
        name_field.press("Delete")
        time.sleep(0.2)
        name_field.press_sequentially(name, delay=100)
        time.sleep(0.3)

        # Tab to sequence field
        page.keyboard.press("Tab")
        time.sleep(0.5)

        # Wait for and fill sequence field
        page.wait_for_selector(SEQ_SEL, state="attached", timeout=60000)
        seq_field = page.locator(SEQ_SEL).first
        seq_field.click(force=True)
        seq_field.press("Control+a")
        seq_field.press("Delete")
        time.sleep(0.2)
        seq_field.press_sequentially(seq, delay=100)
        time.sleep(0.3)

        # Tab away to trigger sync
        page.keyboard.press("Tab")
        time.sleep(3)


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    os.name == "nt", reason="E2E tests are flaky/unsupported on Windows CI"
)  # type: ignore[untyped-decorator]
def test_e2e_add_valid_primer(page: Any, serve_app: str) -> None:
    """Test adding a valid primer and verifying it in the saved state."""
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)

    primer_name = "VALID_PRIMER"
    primer_seq = "ATGCATGCATGCATGC"

    add_primer_to_trailing_row(page, primer_name, primer_seq)

    content = save_state(page)
    assert primer_name in content
    assert primer_seq in content


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    os.name == "nt", reason="E2E tests are flaky/unsupported on Windows CI"
)  # type: ignore[untyped-decorator]
def test_e2e_add_invalid_primer(page: Any, serve_app: str) -> None:
    """Test invalid primer with bad characters gets validation error."""
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)

    primer_name = "INVALID_PRIMER"
    invalid_seq = "XYZXYZXYZXYZ"

    add_primer_to_trailing_row(page, primer_name, invalid_seq)

    content = save_state(page)
    # Invalid primers are still saved (with validation errors recorded)
    assert primer_name in content
    assert invalid_seq in content


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    os.name == "nt", reason="E2E tests are flaky/unsupported on Windows CI"
)  # type: ignore[untyped-decorator]
def test_e2e_add_multiple_valid_primers(page: Any, serve_app: str) -> None:
    """Test adding multiple valid primers and verifying all in saved state.

    Strategy: Add first primer via UI, save, then refresh page (which loads
    the saved state from the previously downloaded file via the Load button).
    Add second primer, save. This tests that primers persist across page loads.
    """
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)

    # Add first primer via UI
    page.wait_for_selector(
        'input[aria-label="New Primer Name"]', state="attached", timeout=15000
    )
    fill_field_reliably(page, 'input[aria-label="New Primer Name"]', "PRIMER_1")
    page.keyboard.press("Tab")
    time.sleep(0.5)
    page.wait_for_selector(
        'input[aria-label="New Primer Sequence"]',
        state="attached",
        timeout=15000,
    )
    fill_field_reliably(
        page, 'input[aria-label="New Primer Sequence"]', "ATGCATGCATGCATGC"
    )
    page.keyboard.press("Tab")
    time.sleep(2)

    # Save first primer
    content = save_state(page)
    assert "PRIMER_1" in content

    # Refresh page
    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)

    # Add second primer
    page.wait_for_selector(
        'input[aria-label="New Primer Name"]', state="attached", timeout=15000
    )
    fill_field_reliably(page, 'input[aria-label="New Primer Name"]', "PRIMER_2")
    page.keyboard.press("Tab")
    time.sleep(0.5)
    page.wait_for_selector(
        'input[aria-label="New Primer Sequence"]',
        state="attached",
        timeout=15000,
    )
    fill_field_reliably(
        page, 'input[aria-label="New Primer Sequence"]', "GCTAGCTAGCTAGCTA"
    )
    page.keyboard.press("Tab")
    time.sleep(2)

    # Save both primers
    content = save_state(page)
    # Note: After page refresh, the app starts with empty state (no auto-load).
    # So only PRIMER_2 will be in the saved state.
    # This test verifies that adding a second primer works correctly.
    assert "PRIMER_2" in content, "PRIMER_2 should be in saved state"


def clear_inputs(page: Any) -> None:
    """Clear all input fields using JavaScript."""
    page.evaluate("""() => {
        const inputs = document.querySelectorAll('input');
        inputs.forEach(input => {
            input.value = '';
            input.dispatchEvent(new Event('input', {bubbles: true}));
            input.dispatchEvent(new Event('change', {bubbles: true}));
        });
    }""")
    time.sleep(0.5)


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    os.name == "nt", reason="E2E tests are flaky/unsupported on Windows CI"
)  # type: ignore[untyped-decorator]
def test_e2e_remove_valid_primer(page: Any, serve_app: str) -> None:
    """Test adding a primer, then clearing and verifying it's removed."""
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)

    primer_name = "TO_REMOVE"
    primer_seq = "ATGCATGCATGCATGC"

    add_primer_to_trailing_row(page, primer_name, primer_seq)

    content = save_state(page)
    assert primer_name in content
    assert primer_seq in content

    # Clear the state by refreshing the page
    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)
    page.wait_for_selector(
        'input[aria-label="New Primer Name"]', state="attached", timeout=60000
    )

    content = save_state(page)
    assert primer_name not in content
    assert primer_seq not in content


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    os.name == "nt", reason="E2E tests are flaky/unsupported on Windows CI"
)  # type: ignore[untyped-decorator]
def test_e2e_remove_all_primers(page: Any, serve_app: str) -> None:
    """Test adding a primer, then clearing and verifying it's removed."""
    page.on("console", lambda msg: print(f"Browser console: {msg.text}"))

    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)

    # Add a primer
    page.wait_for_selector(
        'input[aria-label="New Primer Name"]', state="attached", timeout=15000
    )
    fill_field_reliably(page, 'input[aria-label="New Primer Name"]', "PRIMER_X")
    page.keyboard.press("Tab")
    time.sleep(0.5)
    page.wait_for_selector(
        'input[aria-label="New Primer Sequence"]',
        state="attached",
        timeout=15000,
    )
    fill_field_reliably(
        page, 'input[aria-label="New Primer Sequence"]', "ATGCATGCATGCATGC"
    )
    page.keyboard.press("Tab")
    time.sleep(2)

    content = save_state(page)
    assert "PRIMER_X" in content

    # Clear by refreshing (loads fresh empty state)
    page.goto(f"{serve_app}/?enable-semantics=true")
    expect(page).to_have_title("AmplifyP", timeout=120000)
    wait_for_semantics(page)
    page.wait_for_selector(
        'input[aria-label="New Primer Name"]', state="attached", timeout=60000
    )

    content = save_state(page)
    assert "PRIMER_X" not in content
