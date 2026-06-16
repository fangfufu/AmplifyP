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

"""Tests for GUI utility functions."""

from unittest.mock import MagicMock

import flet as ft

from amplifyp.gui.util import clean_sequence, format_sequence, show_error_dialog


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
    import sys
    from unittest.mock import MagicMock, patch

    from amplifyp.gui.util import get_full_sha, get_git_sha, get_version

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
