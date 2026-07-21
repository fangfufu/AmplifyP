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

"""Tests for GUI utility functions."""

import subprocess
import sys
from unittest.mock import MagicMock

import flet as ft

from amplifyp.gui.utils.data_helpers import clean_sequence, format_sequence
from amplifyp.gui.utils.gui_helpers import show_error_dialog


def test_clean_sequence() -> None:
    """Test clean_sequence removes whitespaces and newlines."""
    assert clean_sequence(" A T \n G C ") == "ATGC"
    assert clean_sequence("") == ""


def test_format_sequence() -> None:
    """Test format_sequence wraps long sequences at wrap_length."""
    seq = "ATGC" * 10  # 40 bp
    formatted = format_sequence(seq, wrap_length=10)
    assert formatted == "ATGCATGCAT\nGCATGCATGC\nATGCATGCAT\nGCATGCATGC"


def test_show_error_dialog() -> None:
    """Test show_error_dialog correctly manages overlay and OK click."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.overlay = []

    show_error_dialog(mock_page, "Error Title", "Test Message")

    # 1. Verify dialog was added to overlay and set to open
    assert len(mock_page.overlay) == 1
    dialog = mock_page.overlay[0]
    assert isinstance(dialog, ft.AlertDialog)
    assert dialog.open is True

    # 2. Click the OK button action
    assert len(dialog.actions) == 1
    ok_button = dialog.actions[0]
    ok_button.on_click(MagicMock())

    # Verify clicking OK closes the dialog (open = False)
    assert dialog.open is False

    # 3. Simulate Flet's on_dismiss event triggered on dismissal
    dialog.on_dismiss(MagicMock())

    # Verify dialog was cleaned up and removed from overlay
    assert dialog not in mock_page.overlay


def test_get_version_and_sha() -> None:
    """Test get_git_sha, get_full_sha, and get_version."""
    from unittest.mock import MagicMock, patch

    from amplifyp.gui.utils.system import get_full_sha, get_git_sha, get_version

    # Test under mocked 'amplifyp.gui.git_sha' module
    mock_git_sha = MagicMock()
    mock_git_sha.GIT_SHA = "pkgsha"
    mock_git_sha.GIT_FULL_SHA = "pkgfullsha"
    with patch.dict(sys.modules, {"amplifyp.gui.git_sha": mock_git_sha}):
        assert get_git_sha() == "pkgsha"
        assert get_full_sha() == "pkgfullsha"

    # Test under mocked 'js' module
    mock_js = MagicMock()
    mock_js.window = MagicMock()
    mock_js.window.__APP_SHA__ = "mockedsha"

    # Hide amplifyp.gui.git_sha to test js fallback
    with patch.dict(sys.modules, {"amplifyp.gui.git_sha": None, "js": mock_js}):
        # Trigger module reload if needed, or intercept imports
        assert get_git_sha() == "mockedsha"
        assert get_full_sha() == "mockedsha"

    # Test get_version uses __version__ defined in __init__.py
    from amplifyp import __version__

    version_str = get_version()
    assert __version__ in version_str


def test_git_fallback_to_dot_git() -> None:
    """Test git utility fallback to reading .git directory directly."""
    from typing import Any
    from unittest.mock import mock_open, patch

    # Mock subprocess.run to raise SubprocessError.
    # Mock sys.modules to remove git_sha and js.
    # Mock os.path.exists and open to simulate .git.
    with (
        patch("subprocess.run") as mock_run,
        patch("os.path.exists") as mock_exists,
        patch(
            "builtins.open", mock_open(read_data="ref: refs/heads/main\n")
        ) as mock_file_open,
        patch.dict(sys.modules, {"amplifyp.gui.git_sha": None, "js": None}),
    ):
        mock_run.side_effect = subprocess.SubprocessError("git failed")

        # Make os.path.exists return True for HEAD and ref file.
        def exists_side_effect(path: str) -> bool:
            return ".git" in path or "git-sha" in path

        mock_exists.side_effect = exists_side_effect

        # Let's customize the file contents read
        file_contents = {
            "HEAD": "ref: refs/heads/test-branch\n",
            "test-branch": "a1b2c3d4e5f6g7h8i9j0a1b2c3d4e5f6g7h8i9j0\n",
        }

        # We can implement a side effect for mock_file_open
        def open_side_effect(file: Any, *args: Any, **kwargs: Any) -> Any:
            filename = str(file)
            for key, val in file_contents.items():
                if filename.endswith(key):
                    # Return mock file handle supporting context manager.
                    mo = mock_open(read_data=val)
                    return mo.return_value
            raise FileNotFoundError(f"File not found: {filename}")

        mock_file_open.side_effect = open_side_effect

        from amplifyp.gui.utils.system import get_full_sha, get_git_sha

        # Verify it resolves correctly
        assert get_git_sha() == "a1b2c3d"
        assert get_full_sha() == "a1b2c3d4e5f6g7h8i9j0a1b2c3d4e5f6g7h8i9j0"

        # Detached HEAD state
        file_contents_detached = {
            "HEAD": "9999999999999999999999999999999999999999\n"
        }

        def open_side_effect_detached(
            file: Any, *args: Any, **kwargs: Any
        ) -> Any:
            filename = str(file)
            for key, val in file_contents_detached.items():
                if filename.endswith(key):
                    mo = mock_open(read_data=val)
                    return mo.return_value
            raise FileNotFoundError(f"File not found: {filename}")

        mock_file_open.side_effect = open_side_effect_detached

        assert get_git_sha() == "9999999"
        assert get_full_sha() == "9999999999999999999999999999999999999999"
