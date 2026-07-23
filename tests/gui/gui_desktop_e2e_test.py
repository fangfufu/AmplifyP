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

"""Desktop End-to-End tests for CLI arguments."""

import os
import subprocess
import sys

import pytest


@pytest.mark.e2e  # type: ignore[untyped-decorator]
@pytest.mark.skipif(
    sys.platform != "linux",
    reason="Desktop E2E tests only run on Linux (with X11/Xvfb)",
)  # type: ignore[untyped-decorator]
def test_desktop_e2e_auto_close() -> None:
    """Test launching desktop app with a state file and auto-close."""
    state_file = os.path.join(
        os.path.dirname(__file__),
        "..",
        "examples",
        "save_states",
        "amplify_4_full_example.yaml",
    )
    assert os.path.exists(state_file), f"State file {state_file} does not exist"

    # Use the same Python interpreter that is running the tests
    cmd = [
        sys.executable,
        "src/main.py",
        "--state",
        state_file,
        "--auto-close",
    ]

    # Run the command with streams redirected to DEVNULL to avoid hangs
    result = subprocess.run(  # noqa: S603
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env=os.environ.copy(),
        check=False,
        timeout=30,
    )

    assert result.returncode == 0, f"App exited with code {result.returncode}."
