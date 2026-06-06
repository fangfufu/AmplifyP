# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

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
