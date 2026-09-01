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

    # Title row
    title_row = col_content.controls[0]
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
    assert "P1" in text_values


def test_dimer_view_max_limit_warning() -> None:
    """Test DimerView warning when dimer count exceeds max render limit."""
    from typing import Any
    from unittest.mock import patch

    from amplifyp.dimer import PrimerDimerGenerator

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "GCATGCATGC", "active": True},
    ]

    view = DimerView(mock_page, input_data, settings)
    mock_dimer = MagicMock()
    # Create 120 mock dimers
    fake_dimers = [mock_dimer] * 120

    def mock_analyse(self: Any) -> None:
        self.primer_dimers = fake_dimers

    with (
        patch.object(PrimerDimerGenerator, "analyse_primers", mock_analyse),
        patch(
            "amplifyp.gui.views.dimer.dimer_view.DimerCard",
            return_value=ft.Container(),
        ),
    ):
        success = view.run_analysis()
        assert success is True
        # First control should be warning container
        assert len(view.result_list.controls) == 101
        warning_container = view.result_list.controls[0]
        assert isinstance(warning_container, ft.Container)
        assert (
            "Warning: 120 primer dimers detected"
            in warning_container.content.value
        )

        # Run again to hit cached state key
        success2 = view.run_analysis()
        assert success2 is True


def test_dimer_view_error_handling() -> None:
    """Test DimerView exception handling during analysis."""
    from unittest.mock import patch

    from amplifyp.dimer import PrimerDimerGenerator

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "GCATGCATGC", "active": True},
    ]

    view = DimerView(mock_page, input_data, settings)
    with (
        patch.object(
            PrimerDimerGenerator,
            "analyse_primers",
            side_effect=RuntimeError("Calculation failure"),
        ),
        patch(
            "amplifyp.gui.views.dimer.dimer_view.show_error_dialog"
        ) as mock_err_dlg,
    ):
        success = view.run_analysis()
        assert success is False
        mock_err_dlg.assert_called_once()
        assert len(view.result_list.controls) == 1
        assert "Error running analysis" in view.result_list.controls[0].value


def test_header_set_update_available() -> None:
    """Test AppHeader.set_update_available updates text and url callback."""
    from unittest.mock import patch

    from amplifyp.gui.views.header import AppHeader

    mock_page = MagicMock(spec=ft.Page)
    mock_page.launch_url = MagicMock()
    settings = GUISettings()

    header = AppHeader(
        settings=settings,
        on_switch_input=MagicMock(),
        on_switch_settings=MagicMock(),
        on_switch_about=MagicMock(),
        on_pcr_click=MagicMock(),
        on_dimers_click=MagicMock(),
        on_save=MagicMock(),
        on_load=MagicMock(),
        pcr_button_ref=ft.Ref[ft.FilledButton](),
        dimers_button_ref=ft.Ref[ft.FilledButton](),
    )

    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        header.set_update_available("v1.0.0")
        assert "(Update v1.0.0 available!)" in header.version_text.value
        assert header.version_text.on_click is not None
        header.version_text.on_click(MagicMock())

    mock_page.launch_url.assert_called_once_with(
        "https://github.com/fangfufu/AmplifyP/releases"
    )


def test_dimer_card_self_dimer_no_names() -> None:
    """Test DimerCard with show_names=False for self-dimer."""
    from amplifyp.dimer import PrimerDimer
    from amplifyp.dna import Primer
    from amplifyp.gui.views.dimer.dimer_card import DimerCard

    p = Primer(sequence="ATGCATGC", name="P1")
    d = PrimerDimer(
        primer_1=p,
        primer_2=p,
        quality=80.0,
        overlap=8,
        p1_pos=0,
    )
    settings = GUISettings()
    card = DimerCard(d, settings, show_names=False)
    assert card is not None
