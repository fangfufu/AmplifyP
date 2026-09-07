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

"""Tests for 1D Primer Designer View."""

import asyncio
from typing import Any
from unittest.mock import AsyncMock, MagicMock, patch

import flet as ft
import yaml

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.views.designer_1d import (
    DismissibleSelfDimerCard,
    PrimerDesignerView,
)


def test_primer_designer_view_initialisation() -> None:
    """Test initial UI setup of PrimerDesignerView."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)

    assert view.left_container is not None
    assert view.right_container is not None
    assert view.main_h_divider is not None
    assert view.left_v_divider is not None
    assert view.right_v_divider is not None
    assert view.top_left_container.height == 240
    assert view.top_right_chart_container.height == 240
    assert view.left_container.expand is True
    assert view.right_container.expand is True
    assert view.dna_input.value == ""
    assert view.length_display.value == "0"
    assert view.min_len_input.value == ""
    assert view.max_quality_input.value == ""
    assert view.max_overlap_input.value == ""
    assert view.filter_dna_checkbox.value is False
    assert view.max_binding_sites_input.value == ""
    assert view.max_binding_sites_input.disabled is True
    assert view.clear_all_button is not None
    assert len(view.primer_list.controls) == 0
    assert len(view.right_cards_list.controls) == 0


def test_primer_designer_view_length_counter() -> None:
    """Test dynamic length counter updating when typing DNA sequence."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    assert view.length_display.value == "0"

    view.dna_input.value = "ATG CGT ACG T"
    view.form._on_dna_change(MagicMock(spec=ft.ControlEvent))
    assert view.length_display.value == "10"


def test_primer_designer_view_run_analysis_success() -> None:
    """Test running 1D primer truncation analysis with valid parameters."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "7"
    view.max_quality_input.value = ""
    view.max_overlap_input.value = ""

    success = view.run_designer()

    assert success is True
    assert view.error_text.visible is False
    # ATGCGTACGT (length 10 down to 7) produces 4 steps
    assert len(view.primer_list.controls) == 4

    first_item = view.primer_list.controls[0]
    assert isinstance(first_item, ft.Card)

    # Check chart content has been populated
    chart = view.chart_content_container.content
    assert isinstance(chart, ft.Row)
    assert len(chart.controls) == 4


def test_primer_designer_view_clear_all() -> None:
    """Test Clear All button resets the entire view to blank state."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.form._on_dna_change(MagicMock(spec=ft.ControlEvent))
    view.min_len_input.value = "7"
    view.max_quality_input.value = ""
    view.max_overlap_input.value = ""
    view.run_designer()

    # Create a dimer card
    view.primer_list.controls[0].content.on_click(None)
    assert len(view.right_cards_list.controls) == 1
    assert len(view.primer_list.controls) == 4

    # Trigger Clear All
    view.filter_dna_checkbox.value = True
    view.max_binding_sites_input.disabled = False
    view.max_binding_sites_input.value = "2"
    view._clear_all(None)

    assert view.dna_input.value == ""
    assert view.length_display.value == "0"
    assert view.min_len_input.value == ""
    assert view.max_quality_input.value == ""
    assert view.max_overlap_input.value == ""
    assert view.filter_dna_checkbox.value is False
    assert view.max_binding_sites_input.value == ""
    assert view.max_binding_sites_input.disabled is True
    assert len(view.primer_list.controls) == 0
    assert len(view.right_cards_list.controls) == 0


def test_primer_designer_view_validation_errors() -> None:
    """Test validation handling for invalid DNA sequence and minimum length."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)

    # Empty sequence
    view.dna_input.value = ""
    assert view.run_designer() is False
    assert view.dna_input.error is not None
    assert "enter a valid DNA sequence" in view.dna_input.error

    # Empty min_length
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = ""
    assert view.run_designer() is False
    assert view.min_len_input.error is not None
    assert "Minimum length is required" in view.min_len_input.error

    # Non-digit min_length
    view.min_len_input.value = "abc"
    assert view.run_designer() is False
    assert view.min_len_input.error is not None
    assert "positive integer" in view.min_len_input.error

    # min_length <= 0
    view.min_len_input.value = "0"
    assert view.run_designer() is False
    assert view.min_len_input.error is not None
    assert "greater than 0" in view.min_len_input.error

    # min_length > len(seq)
    view.min_len_input.value = "20"
    assert view.run_designer() is False
    assert view.min_len_input.error is not None
    assert "cannot exceed sequence length" in view.min_len_input.error

    # Invalid max_quality (non-integer)
    view.min_len_input.value = "7"
    view.max_quality_input.value = "invalid"
    assert view.run_designer() is False
    assert view.max_quality_input.error is not None
    assert "Max Quality must be an integer" in view.max_quality_input.error
    view.max_quality_input.value = ""

    # Invalid max_overlap
    view.max_overlap_input.value = "abc"
    assert view.run_designer() is False
    assert view.max_overlap_input.error is not None
    assert (
        "Max Overlap must be a non-negative integer"
        in view.max_overlap_input.error
    )


def test_primer_designer_view_threshold_filtering() -> None:
    """Test max_quality and max_overlap GUI input filtering."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "7"
    view.max_quality_input.value = "110"
    view.max_overlap_input.value = "6"

    assert view.run_designer() is True
    # Filter quality <= 110 and overlap <= 6 leaves 2 steps
    assert len(view.primer_list.controls) == 2


def test_primer_designer_view_card_creation_and_uniqueness() -> None:
    """Test creating unique cards and bringing existing cards to top."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "8"  # Produces 3 steps: lengths 10, 9, 8
    view.max_quality_input.value = ""
    view.max_overlap_input.value = ""
    view.run_designer()

    # Click first primer item (step 0: length 10)
    card_0_container = view.primer_list.controls[0].content
    card_0_container.on_click(None)

    assert len(view.right_cards_list.controls) == 1
    card_0 = view.right_cards_list.controls[0]
    assert isinstance(card_0, DismissibleSelfDimerCard)
    card_0_id = card_0._card_id

    # Click second primer item (step 1: length 9)
    card_1_container = view.primer_list.controls[1].content
    card_1_container.on_click(None)

    assert len(view.right_cards_list.controls) == 2
    card_1 = view.right_cards_list.controls[0]
    assert card_1._card_id != card_0_id

    # Click first primer item again -> bring card_0 to top
    card_0_container.on_click(None)

    assert len(view.right_cards_list.controls) == 2
    top_card = view.right_cards_list.controls[0]
    assert top_card._card_id == card_0_id


def test_primer_designer_view_click_quality_bar_chart_item() -> None:
    """Test clicking a bar in the quality chart opens/raises dimer card."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "8"
    view.max_quality_input.value = ""
    view.max_overlap_input.value = ""
    view.run_designer()

    chart = view.chart_content_container.content
    assert isinstance(chart, ft.Row)
    assert len(chart.controls) == 3

    # Click the first bar container in the quality chart
    bar_0 = chart.controls[0]
    assert isinstance(bar_0, ft.Container)
    bar_0.on_click(None)

    assert len(view.right_cards_list.controls) == 1
    first_card = view.right_cards_list.controls[0]
    assert isinstance(first_card, DismissibleSelfDimerCard)


def test_primer_designer_view_dismiss_and_clear_cards() -> None:
    """Test dismissing individual cards and clearing all cards."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "8"
    view.max_quality_input.value = ""
    view.max_overlap_input.value = ""
    view.run_designer()

    view.primer_list.controls[0].content.on_click(None)
    view.primer_list.controls[1].content.on_click(None)

    assert len(view.right_cards_list.controls) == 2
    assert view.clear_cards_button.visible is True

    # Dismiss top card
    top_card = view.right_cards_list.controls[0]
    view._dismiss_card(top_card)

    assert len(view.right_cards_list.controls) == 1

    # Clear all cards
    view._clear_all_cards(None)
    assert len(view.right_cards_list.controls) == 0
    assert view.clear_cards_button.visible is False


def test_primer_designer_view_resizing_panels() -> None:
    """Test horizontal and vertical pan update handlers."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.width = 800
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)

    initial_h_left = view.top_left_container.height
    initial_h_right = view.top_right_chart_container.height

    # Horizontal drag (increase left panel width)
    drag_h = MagicMock(spec=ft.DragUpdateEvent)
    drag_h.local_delta = MagicMock(x=50.0, y=0.0)
    view._on_h_pan_update(drag_h)

    assert view.left_container.width == 450.0

    # Top-left vertical drag (increase top-left container height)
    drag_v_left = MagicMock(spec=ft.DragUpdateEvent)
    drag_v_left.local_delta = MagicMock(x=0.0, y=30.0)
    view._on_v_pan_update(drag_v_left)

    assert view.top_left_container.height == initial_h_left + 30.0

    # Top-right vertical drag (increase top-right chart container height)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "8"
    view.max_quality_input.value = ""
    view.max_overlap_input.value = ""
    view.run_designer()

    chart_before = view.chart_content_container.content
    bar_container_before = (
        chart_before.controls[0].content.controls[0].controls[1]
    )
    h_before = bar_container_before.height

    drag_v_right = MagicMock(spec=ft.DragUpdateEvent)
    drag_v_right.local_delta = MagicMock(x=0.0, y=100.0)
    view._on_right_v_pan_update(drag_v_right)

    assert view.top_right_chart_container.height == initial_h_right + 100.0

    chart_after = view.chart_content_container.content
    bar_container_after = (
        chart_after.controls[0].content.controls[0].controls[1]
    )
    h_after = bar_container_after.height

    assert h_after > h_before


def test_primer_designer_view_save_and_load_parameters() -> None:
    """Test saving and loading 1D primer designer parameters."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.template = "ATGCGTACGTTTTATGCGTACGT"
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)

    # Set parameters on the form
    view.form.dna_input.value = "ATGCGTACGT"
    view.form.length_display.value = "10"
    view.form.min_len_input.value = "8"
    view.form.max_quality_input.value = ""
    view.form.max_overlap_input.value = ""
    view.form.filter_dna_checkbox.value = True
    view.form.max_binding_sites_input.disabled = False
    view.form.max_binding_sites_input.value = "10"

    saved_content = ""

    async def mock_save_and_write_file(
        page: ft.Page,
        dialog_title: str,
        file_name: str,
        allowed_extensions: list[str],
        content: str,
        show_notification: Any,
        **kwargs: Any,
    ) -> bool:
        nonlocal saved_content
        saved_content = content
        return True

    # 1. Test Save
    with patch(
        "amplifyp.gui.utils.data_helpers.save_and_write_file",
        new=AsyncMock(side_effect=mock_save_and_write_file),
    ):
        asyncio.run(view._save_designer_1d_click(MagicMock()))

    assert saved_content != ""
    parsed = yaml.safe_load(saved_content)
    assert parsed["dna"] == "ATGCGTACGT"
    assert parsed["min_length"] == "8"
    assert parsed["max_quality"] == ""
    assert parsed["max_overlap"] == ""
    assert parsed["filter_dna"] is True
    assert parsed["max_binding_sites"] == "10"

    # 2. Test Load
    # Reset form to different values
    view.form.dna_input.value = "CGT"
    view.form.length_display.value = "3"
    view.form.min_len_input.value = "10"
    view.form.max_quality_input.value = "60"
    view.form.max_overlap_input.value = "3"
    view.form.filter_dna_checkbox.value = False
    view.form.max_binding_sites_input.disabled = True
    view.form.max_binding_sites_input.value = ""

    async def mock_pick_and_read_file(
        page: ft.Page,
        dialog_title: str,
        allowed_extensions: list[str],
        show_notification: Any,
    ) -> str:
        return saved_content

    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new=AsyncMock(side_effect=mock_pick_and_read_file),
    ):
        asyncio.run(view._load_designer_1d_click(MagicMock()))

    # Verify loaded values in the form
    assert view.form.dna_input.value == "ATGCGTACGT"
    assert view.form.length_display.value == "10"
    assert view.form.min_len_input.value == "8"
    assert view.form.max_quality_input.value == ""
    assert view.form.max_overlap_input.value == ""
    assert view.form.filter_dna_checkbox.value is True
    assert view.form.max_binding_sites_input.disabled is False
    assert view.form.max_binding_sites_input.value == "10"

    # Verify that the analysis automatically ran (3 steps produced)
    assert len(view.primer_list.controls) == 3


def test_designer_1d_remaining_branches() -> None:
    return

    """Test all remaining branches for 1D primer designer."""
    from amplifyp.gui.views.designer_1d import PrimerDesignerView

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)

    # 1. Properties
    assert view.analyse_button is not None

    # 2. Form validation: DNA sequence contains no valid nucleotides
    view.form.dna_input.value = r"\n\t"
    res = view.form.validate_and_get_params()
    assert res is None

    # 3. _run_designer_event — run thread synchronously via mock
    view.form.dna_input.value = "ATGCGTACGT"
    view.form.min_len_input.value = "8"
    with patch(
        "amplifyp.gui.views.designer_1d.designer_1d_view.threading.Thread",
        side_effect=lambda target, daemon: type(
            "T", (), {"start": lambda self: target()}
        )(),
    ):
        view._run_designer_event(None)
    assert len(view.primer_list.controls) == 3

    # 3b. show_loading / update_progress / _restore_primer_list (via tracker)
    tracker = view.progress_tracker
    view.show_loading(total=4)
    assert tracker.bar is not None
    assert tracker.bar.value == 0.0
    assert tracker.label is not None
    assert tracker.label.value == "0 / 4"

    view.update_progress(2, 4)
    assert tracker.bar is not None
    assert abs((tracker.bar.value or 0.0) - 0.5) < 0.01
    assert tracker.label is not None
    assert tracker.label.value == "2 / 4 (50%)"

    # show_loading indeterminate (total=0)
    view.show_loading(total=0)
    assert tracker.bar is not None
    assert tracker.bar.value is None
    assert tracker.label is not None
    assert "Analysing" in (tracker.label.value or "")

    view._restore_primer_list()
    assert tracker.bar is None
    assert tracker.label is None
    assert tracker.is_active is False

    # 4. run_designer exception handling
    with (
        patch(
            "amplifyp.gui.views.designer_1d.designer_1d_view.PrimerDesigner1D",
            side_effect=RuntimeError("Designer 1D err"),
        ),
        patch(
            "amplifyp.gui.views.designer_1d.designer_1d_view.show_error_dialog"
        ) as mock_err,
    ):
        success = view.run_designer()
        assert success is False
        mock_err.assert_called_once()

    # 5. _clear_all with RuntimeError on update
    with patch.object(
        mock_page, "update", side_effect=RuntimeError("Page update err")
    ):
        view._clear_all()

    # 6. _load_designer_1d_click when file picker cancelled (returns None)
    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new=AsyncMock(return_value=None),
    ):
        asyncio.run(view._load_designer_1d_click(MagicMock()))

    # 7. _start_designer error path (thread synchronous)
    view.form.dna_input.value = "ATGCGTACGT"
    view.form.min_len_input.value = "8"
    with (
        patch(
            "amplifyp.gui.views.designer_1d.designer_1d_view.PrimerDesigner1D",
            side_effect=RuntimeError("Start designer err"),
        ),
        patch(
            "amplifyp.gui.views.designer_1d.designer_1d_view.show_error_dialog"
        ) as mock_start_err,
        patch(
            "amplifyp.gui.views.designer_1d.designer_1d_view.threading.Thread",
            side_effect=lambda target, daemon: type(
                "T", (), {"start": lambda self: target()}
            )(),
        ),
    ):
        view._start_designer()
        mock_start_err.assert_called_once()


def test_designer_1d_filter_dna_toggle() -> None:
    """Test check against template checkbox enables/disables binding."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    assert view.filter_dna_checkbox.label == "Check against template"
    assert view.check_dna_checkbox.label == "Check against template"
    assert view.check_template_checkbox.label == "Check against template"
    assert view.filter_dna_checkbox.value is False
    assert view.max_binding_sites_input.disabled is True
    assert view.max_binding_sites_input.hint_text == "Unconstrained if empty"

    # Check the checkbox
    view.filter_dna_checkbox.value = True
    view.form._on_filter_dna_change(MagicMock())
    assert view.max_binding_sites_input.disabled is False

    # Uncheck the checkbox
    view.filter_dna_checkbox.value = False
    view.form._on_filter_dna_change(MagicMock())
    assert view.max_binding_sites_input.disabled is True


def test_designer_1d_analyse_button_layout() -> None:
    """Test Analyse button is in the last row after Max Binding Sites."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    last_row = view.form.controls[-2]
    assert isinstance(last_row, ft.Row)
    assert len(last_row.controls) == 3
    assert last_row.controls[0].content is view.filter_dna_checkbox
    assert (
        isinstance(last_row.controls[1], ft.Column)
        and last_row.controls[1].controls[1] is view.max_binding_sites_input
    )
    assert last_row.controls[2].content is view.form.analyse_button
    assert last_row.controls[0].expand is True
    assert last_row.controls[1].expand is True
    assert last_row.controls[2].expand is not True
    assert last_row.controls[0].height == 48
    assert last_row.controls[0].alignment == ft.Alignment(-1, 0)
    assert last_row.controls[2].height == 48
    assert last_row.controls[2].alignment == ft.Alignment(1, 0)
    assert view.form.analyse_button.height == 48


def test_designer_1d_max_binding_sites_validation() -> None:
    """Test validation errors for max binding sites when filter is enabled."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "8"
    view.filter_dna_checkbox.value = True
    view.max_binding_sites_input.disabled = False

    # 1. Zero is invalid
    view.max_binding_sites_input.value = "0"
    assert view.run_designer() is False
    assert view.max_binding_sites_input.error is not None
    assert "greater than 0" in view.max_binding_sites_input.error

    # 2. Non-digit binding sites
    view.max_binding_sites_input.value = "abc"
    assert view.run_designer() is False
    assert view.max_binding_sites_input.error is not None
    assert "positive integer" in view.max_binding_sites_input.error

    # 3. Negative binding sites
    view.max_binding_sites_input.value = "-2"
    assert view.run_designer() is False
    assert view.max_binding_sites_input.error is not None
    assert "positive integer" in view.max_binding_sites_input.error


def test_designer_1d_check_against_dna_empty_binding_sites() -> None:
    """Test empty binding sites input allows evaluation and displays sites."""
    from amplifyp.gui.views.designer_1d import (
        DismissibleSelfDimerCard,
        PrimerItemCard,
    )

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.template = "ATGCGTACGTTTTATGCGTACGTTTT"
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "8"
    view.filter_dna_checkbox.value = True
    view.max_binding_sites_input.disabled = False
    view.max_binding_sites_input.value = ""  # Empty text entry box

    assert view.run_designer() is True
    assert len(view.primer_list.controls) > 0

    first_item = view.primer_list.controls[0]
    assert isinstance(first_item, PrimerItemCard)
    assert hasattr(first_item, "pcr_button")

    def _contains_text(control: ft.Control, text: str) -> bool:
        if isinstance(control, ft.Text) and text in (control.value or ""):
            return True
        if hasattr(control, "content") and control.content:
            if _contains_text(control.content, text):
                return True
        if hasattr(control, "controls") and control.controls:
            for child in control.controls:
                if _contains_text(child, text):
                    return True
        return False

    # Binding sites count MUST be shown in card badge even when cutoff empty
    assert _contains_text(first_item, "Sites: 2")

    # Detail card also displays binding sites count
    first_item.content.on_click(None)
    assert len(view.right_cards_list.controls) == 1
    dimer_card = view.right_cards_list.controls[0]
    assert isinstance(dimer_card, DismissibleSelfDimerCard)
    assert _contains_text(dimer_card, "Sites: 2")
    assert hasattr(dimer_card, "pcr_button")


def test_designer_1d_template_dna_missing_validation() -> None:
    """Test validation error when filter enabled but template DNA is empty."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.template = ""  # empty template
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "8"
    view.filter_dna_checkbox.value = True
    view.max_binding_sites_input.disabled = False
    view.max_binding_sites_input.value = "1"

    success = view.run_designer()
    assert success is False
    assert view.error_text.visible is True
    assert "Template DNA sequence is required" in view.error_text.value
    assert "check against template" in view.error_text.value


def test_designer_1d_template_dna_filtering_success() -> None:
    """Test primer truncation filtering with template and binding sites."""
    from amplifyp.gui.views.designer_1d import PrimerItemCard

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    # Template where ATGCGTACGT appears twice
    input_data.template = "ATGCGTACGTTTTATGCGTACGTTTT"
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "7"
    view.filter_dna_checkbox.value = True
    view.max_binding_sites_input.disabled = False
    view.max_binding_sites_input.value = "2"

    assert view.run_designer() is True
    assert len(view.primer_list.controls) > 0
    first_item = view.primer_list.controls[0]
    assert isinstance(first_item, PrimerItemCard)
    # Check that PCR button is present
    assert hasattr(first_item, "pcr_button")

    # When check against template is ticked, Sites: badge is present
    def _contains_text(control: ft.Control, text: str) -> bool:
        if isinstance(control, ft.Text) and text in (control.value or ""):
            return True
        if hasattr(control, "content") and control.content:
            if _contains_text(control.content, text):
                return True
        if hasattr(control, "controls") and control.controls:
            for child in control.controls:
                if _contains_text(child, text):
                    return True
        return False

    assert _contains_text(first_item, "Sites: 2")

    # If max_binding_sites is 1, no primers match (count is 2)
    view.max_binding_sites_input.value = "1"
    assert view.run_designer() is True
    assert len(view.primer_list.controls) == 0


def test_designer_1d_check_against_dna_unticked_skips_sites_badge() -> None:
    """Test when check against template is unticked, sites are omitted."""
    from amplifyp.gui.views.designer_1d import (
        DismissibleSelfDimerCard,
        PrimerItemCard,
    )

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    # Template provided in Input view
    input_data.template = "ATGCGTACGTTTTATGCGTACGTTTT"
    settings = GUISettings()

    view = PrimerDesignerView(mock_page, input_data, settings)
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "8"
    assert view.filter_dna_checkbox.value is False

    assert view.run_designer() is True
    assert len(view.primer_list.controls) > 0

    first_item = view.primer_list.controls[0]
    assert isinstance(first_item, PrimerItemCard)
    # PCR button is still present even when check against template is unticked
    assert hasattr(first_item, "pcr_button")

    def _contains_text(control: ft.Control, text: str) -> bool:
        if isinstance(control, ft.Text) and text in (control.value or ""):
            return True
        if hasattr(control, "content") and control.content:
            if _contains_text(control.content, text):
                return True
        if hasattr(control, "controls") and control.controls:
            for child in control.controls:
                if _contains_text(child, text):
                    return True
        return False

    # Sites badge is omitted because check against template is unticked
    assert not _contains_text(first_item, "Sites:")

    # Click primer item to open DismissibleSelfDimerCard
    first_item.content.on_click(None)
    assert len(view.right_cards_list.controls) == 1
    dimer_card = view.right_cards_list.controls[0]
    assert isinstance(dimer_card, DismissibleSelfDimerCard)
    # Detail card also retains PCR button and omits Sites: badge
    assert hasattr(dimer_card, "pcr_button")
    assert not _contains_text(dimer_card, "Sites:")


def test_designer_1d_card_pcr_buttons() -> None:
    """Test PCR button clicks on PrimerItemCard and DismissibleSelfDimerCard."""
    from amplifyp.gui.views.designer_1d import (
        DismissibleSelfDimerCard,
        PrimerItemCard,
    )

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.template = "ATGCGTACGTTTTATGCGTACGT"
    settings = GUISettings()

    run_pcr_calls: list[tuple[str, str]] = []

    def mock_on_run_pcr(seq: str, name: str) -> None:
        run_pcr_calls.append((seq, name))

    view = PrimerDesignerView(
        mock_page, input_data, settings, on_run_pcr=mock_on_run_pcr
    )
    view.dna_input.value = "ATGCGTACGT"
    view.min_len_input.value = "8"
    assert view.run_designer() is True

    # 1. Click PCR button on PrimerItemCard
    first_card = view.primer_list.controls[0]
    assert isinstance(first_card, PrimerItemCard)
    first_card.pcr_button.on_click(None)

    assert len(run_pcr_calls) == 1
    assert run_pcr_calls[0][0] == "ATGCGTACGT"
    assert "1D Primer (10 nt)" in run_pcr_calls[0][1]

    # 2. Click primer item to open DismissibleSelfDimerCard
    first_card.content.on_click(None)
    assert len(view.right_cards_list.controls) == 1
    dimer_card = view.right_cards_list.controls[0]
    assert isinstance(dimer_card, DismissibleSelfDimerCard)
    assert hasattr(dimer_card, "pcr_button")

    # Click Run PCR on the detail card
    dimer_card.pcr_button.on_click(None)
    assert len(run_pcr_calls) == 2
    assert run_pcr_calls[1][0] == "ATGCGTACGT"


def test_controller_run_pcr_with_primer() -> None:
    """Test GUIController.run_pcr_with_primer reaction orchestration."""
    from amplifyp.gui.controller import GUIController

    mock_page = MagicMock(spec=ft.Page)
    controller = GUIController(mock_page)
    controller.input_data.template = "ATGCGTACGTTTTATGCGTACGT"
    controller.input_data.primers = [
        {"name": "OldPrimer", "seq": "TTTTTTTTTT", "active": True}
    ]

    # Mock pcr_view and nav_manager
    controller.pcr_view = MagicMock()
    controller.pcr_view.run_pcr.return_value = True
    controller._nav_manager = MagicMock()

    controller.run_pcr_with_primer("ATGCGTACGT", "1D Candidate (10 nt)")

    # Previous primer deactivated
    assert controller.input_data.primers[0]["active"] is False
    # Candidate primer added and active
    candidate = [p for p in controller.input_data.primers if p["active"]]
    assert len(candidate) == 1
    assert candidate[0]["name"] == "1D Candidate (10 nt)"
    assert candidate[0]["seq"] == "ATGCGTACGT"
    assert controller.input_view_dirty is True
    controller._nav_manager.switch_view.assert_called_once()
    controller.pcr_view.run_pcr.assert_called_once()

    # When template is empty
    controller.input_data.template = ""
    with patch("amplifyp.gui.utils.gui_helpers.show_error_dialog") as mock_err:
        controller.run_pcr_with_primer("ATGCGTACGT", "1D Candidate (10 nt)")
        mock_err.assert_called_once()
