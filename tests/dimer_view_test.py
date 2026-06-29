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

"""Tests for Dimer View."""

from unittest.mock import MagicMock

import flet as ft
import pytest

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.views.dimer import DimerView
from amplifyp.settings import PrimerDimerSettings


def test_get_primer_dimer_settings() -> None:
    """Test getting primer dimer settings from GUI state."""
    settings = GUISettings()
    # Check defaults
    pd_settings = settings.get_primer_dimer_settings()
    assert isinstance(pd_settings, PrimerDimerSettings)
    assert pd_settings.min_overlap == 3
    assert pd_settings.threshold == pytest.approx(60.0)

    # Modify settings in state
    settings["pd_min_overlap"] = "5"
    settings["pd_threshold"] = "75.0"
    pd_settings = settings.get_primer_dimer_settings()
    assert pd_settings.min_overlap == 5
    assert pd_settings.threshold == pytest.approx(75.0)


def test_dimer_view_no_primers() -> None:
    """Test DimerView when no primers are active."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()
    view = DimerView(mock_page, input_data, settings)

    view.run_analysis()

    # result_list should contain a single container showing no dimers
    assert len(view.result_list.controls) == 1
    control = view.result_list.controls[0]
    assert not isinstance(control, ft.Text), getattr(control, "value", "")
    text_control = control.content
    assert "No primer dimers detected" in text_control.value


def test_dimer_view_with_dimers() -> None:
    """Test DimerView runs analysis and populates UI with alignments."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    # Add two active primers known to dimerize
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "GCATGCATGC", "active": True},
    ]
    # Set overlap and threshold to pass the dimer check
    settings["pd_min_overlap"] = "3"
    settings["pd_threshold"] = "30.0"

    view = DimerView(mock_page, input_data, settings)
    view.run_analysis()

    # result_list should contain cards representing matches
    assert len(view.result_list.controls) >= 1
    first_card = view.result_list.controls[0]
    assert isinstance(first_card, ft.Card)

    # Let's inspect card content
    card_container = first_card.content
    assert isinstance(card_container, ft.Container)

    col_content = card_container.content
    assert isinstance(col_content, ft.Column)

    # Title column
    title_col = col_content.controls[0]
    assert isinstance(title_col, ft.Column)
    title_row = title_col.controls[0]
    assert isinstance(title_row, ft.Row)
    title_text = title_row.controls[0].value
    assert "P1" in title_text or "P2" in title_text

    # Code block with alignment
    code_container = col_content.controls[1]
    assert isinstance(code_container, ft.Container)
    alignment_row = code_container.content
    assert isinstance(alignment_row, ft.Row)

    diagram_text = alignment_row.controls[0]
    assert isinstance(diagram_text, ft.Text)

    # Check that sequence indicators are present in the text spans
    text_values = "".join(span.text for span in diagram_text.spans)
    assert "5'-" in text_values
    assert "3'-" in text_values
    # Mid line should contain bond characters like '|' or ':'
    assert "|" in text_values or ":" in text_values
    # Check that the primer names are present as labels
    assert "P1: " in text_values
