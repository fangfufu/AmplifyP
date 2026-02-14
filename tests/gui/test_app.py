"""Tests for the GUI app."""

from unittest.mock import MagicMock

import flet as ft

from amplifyp.gui.app import main


def test_gui_main() -> None:
    """Test the main function of the GUI app."""
    mock_page = MagicMock(spec=ft.Page)

    main(mock_page)

    assert mock_page.title == "AmplifyP - Hello World"
    assert mock_page.vertical_alignment == ft.MainAxisAlignment.CENTER
    mock_page.add.assert_called()
