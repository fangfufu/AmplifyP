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

"""Tests for 2D Primer Designer View GUI components."""

import asyncio
import threading
from typing import Any, cast
from unittest.mock import AsyncMock, MagicMock, patch

import flet as ft
import pytest
import yaml

from amplifyp.dna import DNA, DNADirection
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.views.designer_2d import (
    Designer2DForm,
    Designer2DView,
    Dismissible2DCard,
)
from amplifyp.primer_designer_2d import FilterMetric, PrimerDesigner2D


def _assert_2d_form_defaults(view: Designer2DView) -> None:
    """Assert every 2D designer form field is at its default value."""
    assert view.form.fwd_dna_input.value == ""
    assert view.form.fwd_length_display.value == "0"
    assert view.form.fwd_min_len_input.value == ""
    assert view.form.rev_dna_input.value == ""
    assert view.form.rev_length_display.value == "0"
    assert view.form.rev_min_len_input.value == ""
    assert view.form.max_quality_input.value == ""
    assert view.form.max_overlap_input.value == ""
    assert view.form.filter_dna_checkbox.value is False
    assert view.form.max_amplicons_input.value == ""
    assert view.form.max_amplicons_input.disabled is True


def _sync_run_task(func: Any, *args: Any) -> None:
    """Execute an event-loop task synchronously (unit-test page double).

    The ProgressTracker animation flush loop is skipped: it only stops
    when the analysis finishes, which cannot happen while it blocks the
    test thread.
    """
    if getattr(func, "__name__", "") == "_flush_task":
        return
    asyncio.run(func(*args))


def _run_2d_analysis() -> Designer2DView:
    """Create a 2D view, fill the form, and run the designer synchronously."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.run_task = _sync_run_task
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGATGC"
    view.form.rev_min_len_input.value = "8"
    view.form.max_quality_input.value = ""
    view.form.max_overlap_input.value = ""

    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        view._run_designer_event()
    return view


def test_designer_2d_view_initialisation() -> None:
    """Test initial UI setup of Designer2DView."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)

    assert view.left_container is not None
    assert view.right_container is not None
    assert view.main_h_divider is not None
    assert view.left_v_divider is not None
    assert view.left_container.expand == 1
    assert view.right_container.expand == 1
    _assert_2d_form_defaults(view)
    assert view.form.clear_all_button is not None
    assert len(view.right_cards_list.controls) == 0


def test_designer_2d_view_length_counter() -> None:
    """Test dynamic length counter updating when typing DNA sequences."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)
    assert view.form.fwd_length_display.value == "0"
    assert view.form.rev_length_display.value == "0"

    view.form.fwd_dna_input.value = "ATG CGT ACG T"
    view.form._on_fwd_dna_change(MagicMock(spec=ft.ControlEvent))
    assert view.form.fwd_length_display.value == "10"

    view.form.rev_dna_input.value = "CGT ACG AT GC"
    view.form._on_rev_dna_change(MagicMock(spec=ft.ControlEvent))
    assert view.form.rev_length_display.value == "10"


def test_designer_2d_view_reverse_complement_button() -> None:
    """Test reverse complement button on reverse candidate primer input."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)

    # Empty sequence click does nothing
    view.form.rev_dna_input.value = ""
    view.form._on_reverse_complement_click(None)
    assert view.form.rev_dna_input.value == ""
    assert view.form.rev_length_display.value == "0"

    # Valid sequence is reverse complemented and updates length
    view.form.rev_dna_input.value = "ATGC"
    view.form.rev_length_display.value = "4"
    view.form.rev_dna_input.error = "Error"
    view.form._on_reverse_complement_click(None)
    assert view.form.rev_dna_input.value == "GCAT"
    assert view.form.rev_length_display.value == "4"
    assert view.form.rev_dna_input.error is None

    # Reverse complement again restores original
    view.form._on_reverse_complement_click(None)
    assert view.form.rev_dna_input.value == "ATGC"


def test_designer_2d_form_validation_success() -> None:
    """Test successful input validation in Designer2DForm."""
    settings = GUISettings()
    form = Designer2DForm(settings=settings, on_submit_callback=lambda: None)

    form.fwd_dna_input.value = "ATGCGTACGT"
    form.fwd_min_len_input.value = "8"
    form.rev_dna_input.value = "CGTACGATGC"
    form.rev_min_len_input.value = "8"
    form.max_quality_input.value = "50.0"
    form.max_overlap_input.value = "5"

    (
        fwd_dna,
        fwd_min_len,
        rev_dna,
        rev_min_len,
        threshold,
        max_overlap,
        filter_metric,
        filter_dna_enabled,
        max_amplicons,
    ) = form.validate_and_get_params()

    assert fwd_dna.seq_upper == "ATGCGTACGT"
    assert fwd_min_len == 8
    assert rev_dna.seq_upper == "CGTACGATGC"
    assert rev_min_len == 8
    assert threshold == 50.0
    assert max_overlap == 5
    assert filter_metric == FilterMetric.MAX
    assert filter_dna_enabled is False
    assert max_amplicons is None


def test_designer_2d_form_validation_errors() -> None:
    """Test validation failure cases in Designer2DForm."""
    settings = GUISettings()
    form = Designer2DForm(settings=settings, on_submit_callback=lambda: None)

    # Empty forward sequence and empty reverse sequence simultaneously
    form.fwd_dna_input.value = ""
    form.rev_dna_input.value = ""
    with pytest.raises(ValueError, match="Input validation failed"):
        form.validate_and_get_params()
    assert (
        form.fwd_dna_input.error
        == "Forward candidate primer sequence cannot be empty"
    )
    assert (
        form.rev_dna_input.error
        == "Reverse candidate primer sequence cannot be empty"
    )

    # Empty min lengths
    form.fwd_dna_input.value = "ATGCGTACGT"
    form.rev_dna_input.value = "CGTACGATGC"
    form.fwd_min_len_input.value = ""
    form.rev_min_len_input.value = ""
    with pytest.raises(ValueError, match="Input validation failed"):
        form.validate_and_get_params()
    assert form.fwd_min_len_input.error == "Must be > 0"
    assert form.rev_min_len_input.error == "Must be > 0"

    # Forward min length exceeds sequence length
    form.fwd_dna_input.value = "ATGC"
    form.fwd_min_len_input.value = "10"
    form.rev_dna_input.value = "CGTACGATGC"
    form.rev_min_len_input.value = "8"
    with pytest.raises(ValueError, match="Input validation failed"):
        form.validate_and_get_params()
    assert form.fwd_min_len_input.error == "Exceeds sequence length (4)"

    # Invalid quality filter
    form.fwd_dna_input.value = "ATGCGTACGT"
    form.fwd_min_len_input.value = "8"
    form.rev_dna_input.value = "CGTACGATGC"
    form.rev_min_len_input.value = "8"
    form.max_quality_input.value = "-5.0"
    with pytest.raises(ValueError, match="Input validation failed"):
        form.validate_and_get_params()
    assert form.max_quality_input.error == "Must be >= 0"

    # Invalid max_amplicons when filter_dna is enabled
    form.max_quality_input.value = ""
    form.filter_dna_checkbox.value = True
    form.max_amplicons_input.value = "0"
    with pytest.raises(ValueError, match="Input validation failed"):
        form.validate_and_get_params()
    assert form.max_amplicons_input.error == "Must be > 0"

    form.max_amplicons_input.value = "-1"
    with pytest.raises(ValueError, match="Input validation failed"):
        form.validate_and_get_params()
    assert form.max_amplicons_input.error == "Must be > 0"

    form.max_amplicons_input.value = "invalid"
    with pytest.raises(ValueError, match="Input validation failed"):
        form.validate_and_get_params()
    assert form.max_amplicons_input.error == "Must be > 0"

    # Empty max_amplicons is valid (unconstrained)
    form.max_amplicons_input.value = ""
    res = form.validate_and_get_params()
    assert res[7] is True  # filter_dna_enabled
    assert res[8] is None  # max_amplicons

    # Positive integer max_amplicons is valid
    form.max_amplicons_input.value = "2"
    res = form.validate_and_get_params()
    assert res[7] is True
    assert res[8] == 2


def test_designer_2d_view_run_analysis_and_grid() -> None:
    """Test running 2D analysis populates grid and allows card creation."""
    view = _run_2d_analysis()
    assert view._cached_designer is not None
    # 3 forward lengths (10, 9, 8) x 3 reverse lengths (10, 9, 8) = 9 steps
    assert len(view._cached_designer) == 9

    # Check grid updated
    cell_containers = view.results_grid._cell_containers
    assert len(cell_containers) == 9

    # Click first cell (10, 10)
    step_10_10 = view._cached_designer.get_step(0)
    view.results_grid._on_cell_click(step_10_10, (10, 10))

    # Verify a detail card was created on right panel
    assert len(view.right_cards_list.controls) == 1
    card = view.right_cards_list.controls[0]
    assert isinstance(card, Dismissible2DCard)
    assert (
        "Forward: 10 nt, Reverse: 10 nt" in card._card_id
        or card.step == step_10_10
    )
    # Dimer subcontainers on the active card default to 3 pairs
    assert len(card.dimer_subcontainers.controls) == 3
    labels_default = [
        cast(
            ft.Text,
            cast(
                ft.Row,
                cast(ft.Column, cast(ft.Container, col).content).controls[0],
            ).controls[0],
        ).value
        for col in card.dimer_subcontainers.controls
    ]
    assert labels_default == [
        "Forward Self-Dimer (Fwd-Fwd)",
        "Reverse Self-Dimer (Rev-Rev)",
        "Forward-Reverse Cross-Dimer (Fwd-Rev)",
    ]

    # Enabling setting and calling update_ui updates active card's
    # subcontainers to 4 pairs
    card.settings["designer_2d_show_rev_fwd"] = True
    view.update_ui()
    assert len(card.dimer_subcontainers.controls) == 4
    labels_enabled = [
        cast(
            ft.Text,
            cast(
                ft.Row,
                cast(ft.Column, cast(ft.Container, col).content).controls[0],
            ).controls[0],
        ).value
        for col in card.dimer_subcontainers.controls
    ]
    assert labels_enabled == [
        "Forward Self-Dimer (Fwd-Fwd)",
        "Reverse Self-Dimer (Rev-Rev)",
        "Forward-Reverse Cross-Dimer (Fwd-Rev)",
        "Reverse-Forward Cross-Dimer (Rev-Fwd)",
    ]

    # Disabling setting and calling card.update_ui reverts to 3 pairs
    card.settings["designer_2d_show_rev_fwd"] = False
    card.update_ui()
    assert len(card.dimer_subcontainers.controls) == 3


def test_designer_2d_view_clear_all() -> None:
    """Test Clear All resets all 2D parameters, grid results, and cards."""
    view = _run_2d_analysis()
    assert view._cached_designer is not None
    assert len(view.results_grid._cell_containers) == 9

    # Add a card
    step = view._cached_designer.get_step(0)
    view.results_grid._on_cell_click(step, (10, 10))
    assert len(view.right_cards_list.controls) == 1

    view.form.filter_dna_checkbox.value = True
    view.form.max_amplicons_input.disabled = False
    view.form.max_amplicons_input.value = "3"

    # Clear all
    view._clear_all(None)

    _assert_2d_form_defaults(view)
    assert view._cached_designer is None
    assert len(view.results_grid._cell_containers) == 0
    assert len(view.right_cards_list.controls) == 0


def test_dismissible_2d_card_dismiss_and_clear() -> None:
    """Test card dismiss callback and clear all cards."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)

    fwd_dna = DNA("ATGCGTACGT", direction=DNADirection.FWD)
    rev_dna = DNA("CGTACGATGC", direction=DNADirection.REV)
    designer = PrimerDesigner2D(fwd_dna, 8, rev_dna, 8)

    step = designer.get_step(0)

    # Select step to add card
    view._on_grid_step_selected(step)
    assert len(view.right_cards_list.controls) == 1

    # Dismiss card
    card = view.right_cards_list.controls[0]
    view._dismiss_card(card)
    assert len(view.right_cards_list.controls) == 0

    # Add again and clear all
    view._on_grid_step_selected(step)
    assert len(view.right_cards_list.controls) == 1
    view._clear_all_cards(MagicMock())
    assert len(view.right_cards_list.controls) == 0


def test_grid_cell_click_brings_existing_card_to_top() -> None:
    """Test clicking an existing cell re-orders its detail card to top."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)

    fwd_dna = DNA("ATGCGTACGT", direction=DNADirection.FWD)
    rev_dna = DNA("CGTACGATGC", direction=DNADirection.REV)
    designer = PrimerDesigner2D(fwd_dna, 8, rev_dna, 8)

    step_1 = designer.get_step(0)
    step_2 = designer.get_step(1)

    # Click step_1 then step_2
    view._on_grid_step_selected(step_1)
    view._on_grid_step_selected(step_2)

    assert len(view._active_cards) == 2
    assert view._active_cards[0].step == step_2
    assert view._active_cards[1].step == step_1

    # Click step_1 again -> step_1 should move to top (index 0)
    view._on_grid_step_selected(step_1)
    assert len(view._active_cards) == 2
    assert view._active_cards[0].step == step_1
    assert view._active_cards[1].step == step_2
    assert view.right_cards_list.controls[0].step == step_1


def test_grid_2d_results_view_colour_mapping() -> None:
    """Test Grid2DResultsView applies cell background colours.

    Verifies active scheme applies background colours to grid cells.
    """
    settings = GUISettings()
    settings["designer_2d_colour_scheme"] = "Cool-Warm"

    grid = Designer2DView(
        MagicMock(spec=ft.Page), GUIInput(), settings
    ).results_grid

    fwd_dna = DNA("ATGCGTACGT", direction=DNADirection.FWD)
    rev_dna = DNA("CGTACGATGC", direction=DNADirection.REV)
    designer = PrimerDesigner2D(fwd_dna, 8, rev_dna, 8)

    grid.update_grid(designer)
    assert len(grid._cell_containers) > 0

    first_key = next(iter(grid._cell_containers.keys()))
    first_container = grid._cell_containers[first_key]
    assert first_container.bgcolor is not None
    assert first_container.bgcolor.startswith("#")

    # Select cell -> border changes to PRIMARY but bgcolor remains preserved
    step = designer.get_step(0)
    grid._on_cell_click(step, first_key)
    assert first_container.bgcolor == grid._cell_bg_colours[first_key]


def test_designer_2d_tile_in_settings_view() -> None:
    """Test SettingsView integrates Designer2DTile and updates setting value."""
    from amplifyp.gui.views.settings.settings_view import SettingsView

    mock_page = MagicMock(spec=ft.Page)
    settings = GUISettings()
    settings_view = SettingsView(mock_page, settings)

    assert hasattr(settings_view, "designer_2d_tile")
    assert settings_view.set_designer_2d_colour_scheme is not None
    assert settings_view.set_designer_2d_colour_scheme.value == "Blue-Orange"
    settings_view.set_designer_2d_colour_scheme.value = "Traffic Light"
    mock_event = MagicMock()
    mock_event.control = settings_view.set_designer_2d_colour_scheme
    settings_view._on_change_handler(mock_event)

    assert settings["designer_2d_colour_scheme"] == "Traffic Light"

    assert settings_view.set_designer_2d_show_rev_fwd is not None
    assert settings_view.set_designer_2d_show_rev_fwd.value is False
    settings_view.set_designer_2d_show_rev_fwd.value = True
    mock_event_cb = MagicMock()
    mock_event_cb.control = settings_view.set_designer_2d_show_rev_fwd
    settings_view._on_change_handler(mock_event_cb)

    assert settings["designer_2d_show_rev_fwd"] is True


def test_designer_2d_view_save_and_load_parameters() -> None:
    """Test saving and loading 2D primer designer parameters."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)

    # Set parameters on the form
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGATGC"
    view.form.rev_min_len_input.value = "8"
    view.form.max_quality_input.value = ""
    view.form.max_overlap_input.value = ""

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
        asyncio.run(view._save_designer_2d_click(MagicMock()))

    assert saved_content != ""
    parsed = yaml.safe_load(saved_content)
    assert parsed["fwd_dna"] == "ATGCGTACGT"
    assert parsed["fwd_min_length"] == "8"
    assert parsed["rev_dna"] == "CGTACGATGC"
    assert parsed["rev_min_length"] == "8"
    assert parsed["max_quality"] == ""
    assert parsed["max_overlap"] == ""
    assert parsed["filter_dna"] is False
    assert parsed["max_amplicons"] == ""

    # 2. Test Load
    # Reset form to different values
    view.form.fwd_dna_input.value = "CGT"
    view.form.fwd_min_len_input.value = "10"
    view.form.rev_dna_input.value = "ATG"
    view.form.rev_min_len_input.value = "10"
    view.form.max_quality_input.value = "60.0"
    view.form.max_overlap_input.value = "3"
    view.form.filter_dna_checkbox.value = True
    view.form.max_amplicons_input.disabled = False
    view.form.max_amplicons_input.value = "5"

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
        asyncio.run(view._load_designer_2d_click(MagicMock()))

    # Verify loaded values in the form
    assert view.form.fwd_dna_input.value == "ATGCGTACGT"
    assert view.form.fwd_length_display.value == "10"
    assert view.form.fwd_min_len_input.value == "8"
    assert view.form.rev_dna_input.value == "CGTACGATGC"
    assert view.form.rev_length_display.value == "10"
    assert view.form.rev_min_len_input.value == "8"
    assert view.form.max_quality_input.value == ""
    assert view.form.max_overlap_input.value == ""
    assert view.form.filter_dna_checkbox.value is False
    assert view.form.max_amplicons_input.value == ""
    assert view.form.max_amplicons_input.disabled is True

    # Verify that the analysis automatically ran (populates _cached_designer)
    assert view._cached_designer is not None
    assert len(view._cached_designer) == 9


def test_designer_2d_and_base_remaining_branches() -> None:
    """Test remaining branches in designer 2D, grid, and base classes."""
    from amplifyp.dna import Primer
    from amplifyp.gui.views.designer.designer_card_helpers import (
        build_primer_summary_row,
        format_primer_properties,
    )
    from amplifyp.gui.views.designer.designer_form_base import BaseDesignerForm

    mock_page = MagicMock(spec=ft.Page)
    mock_page.width = 800.0
    mock_page.height = 600.0
    mock_page.run_task = _sync_run_task
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)

    # 1. Properties
    assert view.form.quality_filter_input is not None
    assert view.form.overlap_filter_input is not None

    # 2. Form base methods & validations
    base_form = BaseDesignerForm(
        settings=settings,
        on_submit_callback=MagicMock(),
        on_clear_error_callback=MagicMock(),
    )

    base_form._on_submit_event(None)
    assert cast(MagicMock, base_form.on_submit_callback).called

    # _clear_field_error with page update
    mock_field = ft.TextField(value="abc")
    mock_field.error = "Err"
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        base_form._clear_field_error(MagicMock(control=mock_field))
        assert mock_field.error is None

    # _clear_field_error when page.update raises RuntimeError
    mock_field.error = "Err"
    with (
        patch.object(mock_page, "update", side_effect=RuntimeError("Err")),
        patch.object(ft.Control, "page", new=property(lambda self: mock_page)),
    ):
        base_form._clear_field_error(MagicMock(control=mock_field))

    # show_field_error with page
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        base_form.show_field_error(mock_field, "Test error")
        assert mock_field.error == "Test error"

    # show_error
    base_form.show_error("General error")
    assert base_form.error_text.value == "General error"

    # validate_max_quality
    base_form.max_quality_input.value = "-5"
    _q_val, is_q_v = base_form.validate_max_quality(int_only=True)
    assert is_q_v is False
    assert "non-negative integer" in (base_form.max_quality_input.error or "")

    base_form.max_quality_input.value = "-5"
    _q_val, is_q_v = base_form.validate_max_quality(int_only=False)
    assert is_q_v is False
    assert "non-negative" in (base_form.max_quality_input.error or "")

    base_form.max_quality_input.value = "abc"
    _q_val, is_q_v = base_form.validate_max_quality(int_only=False)
    assert is_q_v is False
    assert "must be an integer" in (base_form.max_quality_input.error or "")

    # validate_max_overlap
    base_form.max_overlap_input.value = "-3"
    _o_val, is_o_v = base_form.validate_max_overlap()
    assert is_o_v is False
    assert "non-negative" in (base_form.max_overlap_input.error or "")

    # 3. 2D Form validations
    # (negative min length and invalid overlap with page update)
    view.form.fwd_dna_input.value = "ATGC"
    view.form.fwd_min_len_input.value = "-1"
    view.form.rev_dna_input.value = "ATGC"
    view.form.rev_min_len_input.value = "0"
    view.form.max_overlap_input.value = "invalid"

    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        with pytest.raises(ValueError):
            view.form.validate_and_get_params()

    # 4. Designer2DView run_designer error handling
    # ValueError from validate_and_get_params
    view._run_designer_event()

    # PrimerDesigner2D exception
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGATGC"
    view.form.rev_min_len_input.value = "8"
    view.form.max_quality_input.value = ""
    view.form.max_overlap_input.value = ""

    with (
        patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.PrimerDesigner2D",
            side_effect=RuntimeError("2D Designer error"),
        ),
        patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.show_error_dialog"
        ) as mock_err_dlg,
        patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
            side_effect=lambda target, daemon: MagicMock(start=target),
        ),
    ):
        view._run_designer_event()
        mock_err_dlg.assert_called_once()

    # 5. Designer2DView UI updates & dismiss with RuntimeError
    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        view._run_designer_event()
    assert view._cached_designer is not None

    with patch.object(mock_page, "update", side_effect=RuntimeError("Err")):
        view.update_ui()
        view._clear_all()

    # 6. _load_designer_2d_click when cancelled (returns None)
    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new=AsyncMock(return_value=None),
    ):
        asyncio.run(view._load_designer_2d_click(MagicMock()))

    # 7. DesignerViewBase: pan resize and YAML error handling
    # Horizontal pan when left_container.width is a float
    view.left_container.width = 400.0
    view._on_h_pan_update(MagicMock(local_delta=MagicMock(x=20.0)))
    assert view.left_container.width >= 400.0

    # Vertical pan when top_left_container.height is None
    view.top_left_container.height = None
    view._on_v_pan_update(MagicMock(local_delta=MagicMock(y=20.0)))
    assert view.top_left_container.height is not None

    # Load YAML with invalid types / malformed YAML
    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new=AsyncMock(return_value="[1, 2, 3]"),  # Not a dict
    ):
        res = asyncio.run(view._load_parameters_yaml("Test"))
        assert res is None

    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new=AsyncMock(return_value=":::invalid yaml"),
    ):
        res = asyncio.run(view._load_parameters_yaml("Test"))
        assert res is None

    # 8. Card helpers: format_primer_properties exception and copy button
    p_bad = Primer("N" * 10, name="Bad")
    with patch.object(
        settings, "calculate_primer_tm", side_effect=RuntimeError("Tm err")
    ):
        tm_t, _pct_t = format_primer_properties(p_bad, settings)
        assert tm_t == "Tm: N/A"

    col_btn = build_primer_summary_row("Test Primer", p_bad, settings)
    copy_btn = col_btn.controls[1].controls[1]

    captured_task = None

    def mock_run_task(task_fn: Any, *args: Any) -> None:
        nonlocal captured_task
        captured_task = task_fn(*args)

    mock_page.run_task = mock_run_task
    # Trigger copy click with page
    copy_btn.on_click(MagicMock(page=mock_page))
    assert captured_task is not None
    with patch.object(ft.Clipboard, "set", new=AsyncMock()) as mock_clip:
        asyncio.run(captured_task)
        mock_clip.assert_called_once_with(p_bad.seq)

    # Trigger with RuntimeError on page
    with patch.object(
        mock_page, "run_task", side_effect=RuntimeError("Task err")
    ):
        copy_btn.on_click(MagicMock(page=mock_page))

    # The copy-button section above rebinds run_task to a capturing mock
    # (mock_run_task), which would create an un-awaited ProgressTracker
    # flush coroutine on any later show_loading call. Reset to a no-op.
    mock_page.run_task = MagicMock()

    # 9. Dismissible2DCard: settings change with page update and string boolean
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGATGC"
    view.form.rev_min_len_input.value = "8"
    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        view._run_designer_event()
    assert view._cached_designer is not None

    step = view._cached_designer.all_steps[0]
    card = Dismissible2DCard(
        card_id="card_1",
        step=step,
        settings=settings,
        dismiss_callback=MagicMock(),
    )
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        card.on_settings_change()
        card.update_ui()

    card.settings = {"designer_2d_show_rev_fwd": "false"}  # type: ignore[assignment]
    sub_col_false = card._build_dimer_subcontainers(14, 12)
    assert len(sub_col_false.controls) == 3

    card.settings = {"designer_2d_show_rev_fwd": "true"}  # type: ignore[assignment]
    sub_col = card._build_dimer_subcontainers(14, 12)
    assert len(sub_col.controls) == 4

    # 10. Grid2DResults: clear_grid with page, empty steps with RuntimeError,
    # missing cell ("N/A"), and on_cell_click
    grid = view.results_grid
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        grid.clear_grid()

        # Empty steps with update exception
        mock_empty_designer = MagicMock(all_steps=[])
        with patch.object(mock_page, "update", side_effect=RuntimeError("Err")):
            grid.update_grid(mock_empty_designer)
        assert (
            "No valid 2D truncation"
            in grid.content_column.controls[1].content.value
        )

        # Missing cell in matrix (steps at (8,8) and (9,9) creating N/A cells)
        step_8_8 = view._cached_designer.all_steps[0]
        step_9_9 = view._cached_designer.all_steps[-1]
        mock_diagonal_designer = MagicMock()
        mock_diagonal_designer.all_steps = [step_8_8, step_9_9]
        grid.update_grid(mock_diagonal_designer)

        # show_loading with a known total shows a ProgressBar at 0%
        tracker = grid.progress_tracker
        grid.show_loading(total=6)
        assert tracker.bar is not None
        assert tracker.bar.value == 0.0
        assert tracker.label is not None
        assert tracker.label.value == "0 / 6"

        # show_loading with total=0 shows indeterminate bar
        grid.show_loading(total=0)
        assert tracker.bar is not None
        assert tracker.bar.value is None  # indeterminate
        assert tracker.label is not None
        assert "Analysing" in (tracker.label.value or "")

        # update_progress advances bar and label
        grid.show_loading(total=6)
        grid.update_progress(3, 6)
        assert tracker.bar is not None
        assert abs((tracker.bar.value or 0.0) - 0.5) < 0.01
        assert tracker.label is not None
        assert tracker.label.value == "3 / 6 (50%)"

        # update_progress is a no-op when hidden
        grid.clear_grid()
        grid.update_progress(1, 6)  # should not raise

        # on_cell_click
        key = (len(step.fwd_fwd.primer_1.seq), len(step.rev_rev.primer_1.seq))
        grid._on_cell_click(step, key)

    # 11. DesignerViewBase: _bring_card_to_top_or_add & _dismiss_card
    view._bring_card_to_top_or_add("card_1", lambda: card)
    view._bring_card_to_top_or_add("card_1", lambda: card)  # Move to top branch

    with patch.object(
        mock_page, "update", side_effect=RuntimeError("Update err")
    ):
        view._bring_card_to_top_or_add(
            "card_new",
            lambda: Dismissible2DCard("card_new", step, settings, MagicMock()),
        )
        view._bring_card_to_top_or_add(
            "card_new",
            lambda: Dismissible2DCard("card_new", step, settings, MagicMock()),
        )
        view._dismiss_card(card)

    # YAML load exception
    with (
        patch(
            "amplifyp.gui.utils.data_helpers.pick_and_read_file",
            new=AsyncMock(return_value="valid: yaml"),
        ),
        patch("yaml.safe_load", side_effect=Exception("YAML load fail")),
    ):
        res = asyncio.run(view._load_parameters_yaml("Test"))
        assert res is None


def test_designer_2d_form_checkbox_toggle() -> None:
    """Test toggling check against template checkbox enables/disables input."""
    settings = GUISettings()
    form = Designer2DForm(settings=settings, on_submit_callback=lambda: None)

    assert form.check_template_checkbox is form.filter_dna_checkbox
    assert form.check_dna_checkbox is form.filter_dna_checkbox
    assert form.max_amplicon_input is form.max_amplicons_input

    # Initially disabled
    assert form.max_amplicons_input.disabled is True

    # Toggle on
    form.filter_dna_checkbox.value = True
    form._on_filter_dna_change(MagicMock())
    assert form.max_amplicons_input.disabled is False

    # Set error, then toggle off
    form.max_amplicons_input.error = "Error"
    form.filter_dna_checkbox.value = False
    form._on_filter_dna_change(MagicMock())
    assert form.max_amplicons_input.disabled is True
    assert form.max_amplicons_input.error is None


def test_designer_2d_view_template_and_amplicons() -> None:
    """Test template DNA evaluation and amplicon filtering in Designer2DView."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.run_task = _sync_run_task
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "ACGTACGCAT"
    view.form.rev_min_len_input.value = "8"
    view.form.max_quality_input.value = ""
    view.form.max_overlap_input.value = ""
    view.form.filter_dna_checkbox.value = True

    # 1. Missing template -> shows form error
    view._run_designer_event()
    assert view.form.error_text.visible is True
    assert "Template DNA sequence is required" in view.form.error_text.value
    assert view._cached_designer is None

    # 2. Provide template in input_data -> unconstrained amplicons
    input_data.template = "ATGCGTACGTTTTATGCGTACGTTTTATGCGTACGT"
    view.form.max_amplicons_input.value = ""
    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        view._run_designer_event()

    assert view._cached_designer is not None
    assert len(view._cached_designer.all_steps) == 9
    for step in view._cached_designer.all_steps:
        assert step.amplicon_count is not None
        assert step.amplicon_count >= 1

    # Select step and verify card displays Amplicons badge
    first_step = view._cached_designer.all_steps[0]
    view._on_grid_step_selected(first_step)
    assert len(view._active_cards) == 1
    card = cast(Dismissible2DCard, view._active_cards[0])
    assert card.amplicon_count == first_step.amplicon_count

    # 3. Constrained max_amplicons filter (e.g. max 1 amplicon)
    view.form.max_amplicons_input.value = "1"
    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        view._run_designer_event()
    assert view._cached_designer is not None
    for step in view._cached_designer.all_steps:
        assert step.amplicon_count is not None
        assert step.amplicon_count <= 1


def test_designer_2d_view_run_pcr_callback() -> None:
    """Test Dismissible2DCard Run PCR button triggers PCR callback."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.run_task = _sync_run_task
    input_data = GUIInput()
    settings = GUISettings()
    pcr_mock = MagicMock()

    view = Designer2DView(mock_page, input_data, settings, on_run_pcr=pcr_mock)
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "ACGTACGCAT"
    view.form.rev_min_len_input.value = "8"

    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        view._run_designer_event()
    assert view._cached_designer is not None

    step = view._cached_designer.all_steps[0]
    view._on_grid_step_selected(step)
    card = cast(Dismissible2DCard, view._active_cards[0])

    assert hasattr(card, "pcr_button")
    card.pcr_button.on_click(MagicMock())
    pcr_mock.assert_called_once()
    args = pcr_mock.call_args[0]
    assert args[0] == step.fwd_fwd.primer_1.seq
    assert args[2] == step.rev_rev.primer_1.seq


def test_controller_run_pcr_with_primer_pair() -> None:
    """Test GUIController run_pcr_with_primer_pair sets primers and runs PCR."""
    from amplifyp.gui.controller import GUIController

    mock_page = MagicMock(spec=ft.Page)
    mock_page.overlay = []
    mock_page.window = MagicMock()

    with (
        patch("amplifyp.gui.controller.UpdateManager"),
        patch("amplifyp.gui.controller.NotificationHelper"),
    ):
        ctrl = GUIController(mock_page)

    # 1. Missing template -> error dialog
    ctrl.input_data.template = ""
    with patch("amplifyp.gui.utils.gui_helpers.show_error_dialog") as mock_err:
        ctrl.run_pcr_with_primer_pair("ATGCGTACGT", "Fwd", "CGTACGATGC", "Rev")
        mock_err.assert_called_once()
        assert "Template Required" in mock_err.call_args[0][1]

    # 2. Valid template -> sets active primers and executes PCR
    ctrl.input_data.template = "ATGCGTACGTTTTATGCGTACGTTTTATGCGTACGT"
    ctrl.input_data.primers = [
        {"name": "Old", "seq": "TTTTTTTTTT", "active": True}
    ]
    ctrl.pcr_view = MagicMock()
    ctrl.pcr_view.run_pcr.return_value = True
    ctrl._nav_manager = MagicMock()

    ctrl.run_pcr_with_primer_pair("ATGCGTACGT", "Fwd 1", "CGTACGATGC", "Rev 1")
    ctrl.pcr_view.run_pcr.assert_called_once()

    active_primers = [p for p in ctrl.input_data.primers if p.get("active")]
    assert len(active_primers) == 2
    assert active_primers[0]["name"] == "Fwd 1"
    assert active_primers[0]["seq"] == "ATGCGTACGT"
    assert active_primers[1]["name"] == "Rev 1"
    assert active_primers[1]["seq"] == "CGTACGATGC"


def test_schedule_on_event_loop_fallbacks() -> None:
    """Test _schedule_on_event_loop running-loop and detached-page paths."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.run_task.side_effect = RuntimeError("detached")
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)

    # Detached page + no running loop: callback discarded, not executed
    # on a temporary worker-thread loop.
    executed = False

    async def _discarded() -> None:
        nonlocal executed
        executed = True

    view._schedule_on_event_loop(_discarded)
    assert executed is False

    # Running loop on the current thread: task created on that loop.
    async def _drive() -> None:
        nonlocal executed
        view._schedule_on_event_loop(_discarded)
        await asyncio.sleep(0)
        assert executed is True

    asyncio.run(_drive())


def test_designer_2d_analyse_button_toggles_to_abort() -> None:
    """Test Analyse button becomes Abort during a run and restores after."""
    mock_page = MagicMock(spec=ft.Page)
    view = Designer2DView(mock_page, GUIInput(), GUISettings())

    # Simulate an in-flight analysis
    view._analysis_running = True
    view._cancel_event = threading.Event()
    view._set_button_abort_mode(True)
    button = view.form.analyse_button
    assert button.content == "Abort"
    assert button.icon == ft.Icons.STOP
    assert button.disabled is False

    # Click while running routes to abort, not a new run
    view._run_designer_event()
    assert view._cancel_event is not None
    assert view._cancel_event.is_set()

    # Finishing the run restores the button
    view._analysis_running = False
    asyncio.run(view._on_analysis_finished())
    assert button.content == "Analyse"
    assert button.icon == ft.Icons.PLAY_ARROW
    assert view._cancel_event is None
    assert view._analysis_running is False


def test_designer_2d_abort_returns_partial_results() -> None:
    """Test aborting mid-run keeps steps analysed so far."""
    from amplifyp.primer_designer_2d import PrimerDesigner2D as RealDesigner

    mock_page = MagicMock(spec=ft.Page)
    view = Designer2DView(mock_page, GUIInput(), GUISettings())
    # 3 fwd lengths (10, 9, 8) x 2 rev lengths (10, 9) = 6 combinations
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGTACG"
    view.form.rev_min_len_input.value = "9"

    def factory(**kwargs: Any) -> Any:
        cancel = kwargs.get("cancel_event")
        orig_progress = kwargs.get("on_progress")

        def progress(done: int, total: int) -> None:
            if cancel is not None:
                cancel.set()
            if orig_progress is not None:
                orig_progress(done, total)

        kwargs["on_progress"] = progress
        return RealDesigner(**kwargs)

    with (
        patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.PrimerDesigner2D",
            side_effect=factory,
        ),
        patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
            side_effect=lambda target, daemon: type(
                "T", (), {"start": lambda self: target()}
            )(),
        ),
        patch.object(
            view,
            "_schedule_on_event_loop",
            side_effect=lambda func, *args: asyncio.run(func(*args)),
        ),
    ):
        view._run_designer_event()

    # Analysis aborted after the first combination; partial step retained.
    assert view._cached_designer is not None
    assert view._cached_designer.aborted is True
    assert len(view._cached_designer.all_steps) == 1

    # Button restored to Analyse mode after the aborted run.
    assert view.form.analyse_button.content == "Analyse"
    assert view.form.analyse_button.icon == ft.Icons.PLAY_ARROW
    assert view._analysis_running is False
    assert view._cancel_event is None


def test_designer_form_base_build_filter_row_with_button() -> None:
    """Test _build_filter_row with include_analyse_button."""
    settings = GUISettings()
    form = Designer2DForm(settings, on_submit_callback=lambda: None)
    row = form._build_filter_row(
        extra_controls=[ft.Text("Extra")], include_analyse_button=True
    )
    assert len(row.controls) == 4


def test_designer_2d_properties_and_branches() -> None:
    """Test 2D designer edge cases and update exception handling."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()
    view = Designer2DView(mock_page, input_data, settings)

    # 1. _run_designer_event with missing template and update RuntimeError
    mock_page.update.side_effect = RuntimeError("update failed")
    view.form.fwd_dna_input.value = "ATGCATGCATGC"
    view.form.fwd_min_len_input.value = "10"
    view.form.rev_dna_input.value = "CGTACGTACGTA"
    view.form.rev_min_len_input.value = "10"
    view.form.filter_dna_checkbox.value = True
    view.input_data.template = ""
    view._run_designer_event()
    assert "Template DNA sequence is required" in (
        view.form.error_text.value or ""
    )

    # 2. _run_designer_event with template and update error on abort setup
    view.input_data.template = "ATGCATGCATGCATGCATGC"
    with (
        patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.PrimerDesigner2D"
        ),
        patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
            side_effect=lambda target, daemon: type(
                "T", (), {"start": lambda self: None}
            )(),
        ),
    ):
        view._run_designer_event()
        assert view._analysis_running is True
    view._analysis_running = False

    # 3. _on_analysis_finished catches RuntimeError
    asyncio.run(view._on_analysis_finished())
    assert view.form.analyse_button.disabled is False

    # 3b. _on_analysis_finished with stale cancel_event returns early
    current_event = threading.Event()
    stale_event = threading.Event()
    view._cancel_event = current_event
    view._set_button_abort_mode(True)
    asyncio.run(view._on_analysis_finished(stale_event))
    assert view._cancel_event is current_event
    assert view.form.analyse_button.content == "Abort"
    # Call with matching cancel_event restores button and clears cancel event
    asyncio.run(view._on_analysis_finished(current_event))
    assert view._cancel_event is None
    assert view.form.analyse_button.content == "Analyse"

    # 4. _handle_run_pcr with no on_run_pcr and missing template
    view.on_run_pcr = None
    view.input_data.template = ""
    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.show_error_dialog"
    ) as mock_err:
        view._handle_run_pcr("ATGC", "Fwd", "CGTA", "Rev")
        mock_err.assert_called_once()

    # 5. Designer2DForm filter_dna_checkbox change updates page
    form = view.form
    with patch.object(type(form), "page", new=property(lambda s: mock_page)):
        mock_page.update.reset_mock()
        mock_page.update.side_effect = None
        form.filter_dna_checkbox.value = False
        form._on_filter_dna_change(MagicMock())
        mock_page.update.assert_called()

        # 6. Designer2DForm _on_reverse_complement_click on invalid DNA
        from amplifyp.errors import InvalidDNASequenceError

        form.rev_dna_input.value = "ATGC"
        with patch(
            "amplifyp.gui.views.designer_2d.designer_2d_form.DNA",
            side_effect=InvalidDNASequenceError(["X"]),
        ):
            mock_page.update.side_effect = RuntimeError("update error")
            form._on_reverse_complement_click(None)
            assert "Invalid DNA sequence" in (form.rev_dna_input.error or "")

        # 7. Designer2DForm _on_reverse_complement_click success updates page
        mock_page.update.side_effect = None
        form.rev_dna_input.value = "ATGC"
        form._on_reverse_complement_click(None)
        assert form.rev_dna_input.value == "GCAT"
        assert form.rev_dna_input.error is None
        mock_page.update.assert_called()


def test_designer_2d_grid_matrix_increasing_order() -> None:
    """Test 2D grid matrix runs increasing length in both dimensions."""
    mock_page = MagicMock(spec=ft.Page)
    settings = GUISettings()
    view = Designer2DView(mock_page, GUIInput(), settings)

    fwd_dna = DNA("ATGCGTACGT", direction=DNADirection.FWD)  # length 10
    rev_dna = DNA("CGTACGATGC", direction=DNADirection.REV)  # length 10
    designer = PrimerDesigner2D(fwd_dna, 8, rev_dna, 8)  # lengths 8, 9, 10

    view.results_grid.update_grid(designer)

    # Traverse to Column with grid_rows:
    # content_column.controls[1] -> Container -> Row -> Container -> Column
    outer_container = view.results_grid.content_column.controls[1]
    assert isinstance(outer_container, ft.Container)
    inner_row = outer_container.content
    assert isinstance(inner_row, ft.Row)
    inner_container = inner_row.controls[0]
    assert isinstance(inner_container, ft.Container)
    grid_rows_column = inner_container.content
    assert isinstance(grid_rows_column, ft.Column)
    grid_rows = grid_rows_column.controls

    # Check header row (columns: Rev \ Fwd, 8 nt, 9 nt, 10 nt)
    header_row = grid_rows[0]
    assert isinstance(header_row, ft.Row)
    header_texts = [
        cell.content.value
        for cell in header_row.controls
        if isinstance(cell, ft.Container) and isinstance(cell.content, ft.Text)
    ]
    assert header_texts == ["Rev \\ Fwd", "8 nt", "9 nt", "10 nt"]

    # Check data rows (rows: 8 nt, 9 nt, 10 nt)
    data_rows = grid_rows[1:]
    row_headers = [
        row.controls[0].content.value
        for row in data_rows
        if isinstance(row, ft.Row)
        and isinstance(row.controls[0], ft.Container)
        and isinstance(row.controls[0].content, ft.Text)
    ]
    assert row_headers == ["8 nt", "9 nt", "10 nt"]


def test_designer_2d_dimer_card_title_and_badges() -> None:
    """Test 2D panel title is renamed and cards omit mean metrics."""
    from amplifyp.gui.views.designer_2d.dismissible_2d_card import (
        Dismissible2DCard,
    )

    mock_page = MagicMock(spec=ft.Page)
    settings = GUISettings()
    view = Designer2DView(mock_page, GUIInput(), settings)

    # Check panel title
    assert view.right_title.value == "2D Primer Pair Dimer Cards"

    fwd_dna = DNA("ATGCGTACGT", direction=DNADirection.FWD)
    rev_dna = DNA("CGTACGATGC", direction=DNADirection.REV)
    designer = PrimerDesigner2D(fwd_dna, 8, rev_dna, 8)
    step = designer.get_step(0)

    card = Dismissible2DCard(
        card_id="test_card",
        step=step,
        settings=settings,
        dismiss_callback=MagicMock(),
    )

    # Check title metric badges in card header row
    card_container = card.content
    assert isinstance(card_container, ft.Container)
    card_col = card_container.content
    assert isinstance(card_col, ft.Column)
    title_row = card_col.controls[0]
    assert isinstance(title_row, ft.Row)

    badge_texts = []
    for ctrl in title_row.controls:
        if isinstance(ctrl, ft.Container) and isinstance(ctrl.content, ft.Text):
            badge_texts.append(ctrl.content.value or "")

    # Max Quality and Max Overlap must be present
    assert any("Max Quality:" in text for text in badge_texts)
    assert any("Max Overlap:" in text for text in badge_texts)

    # Mean Quality and Mean Overlap must NOT be present
    assert not any("Mean Quality:" in text for text in badge_texts)
    assert not any("Mean Overlap:" in text for text in badge_texts)


def test_dismissible_2d_card_separator_setting() -> None:
    """Test Dismissible2DCard formats boundaries using configured separator."""
    from amplifyp.dna import DNA, DNADirection
    from amplifyp.gui.settings import GUISettings
    from amplifyp.gui.views.designer_2d.dismissible_2d_card import (
        Dismissible2DCard,
    )
    from amplifyp.primer_designer_2d import PrimerDesigner2D

    fwd_dna = DNA("ATGCGTACGT", direction=DNADirection.FWD)
    rev_dna = DNA("CGTACGATGC", direction=DNADirection.REV)
    designer = PrimerDesigner2D(fwd_dna, 8, rev_dna, 8)
    step = designer.get_step(0)

    # Test Space separator
    settings_space = GUISettings()
    settings_space["dimer_sequence_separator"] = "Space (' ')"
    card_space = Dismissible2DCard(
        "card_sp", step, settings_space, dismiss_callback=MagicMock()
    )
    first_sub_space = card_space.dimer_subcontainers.controls[0]
    first_diag_space = first_sub_space.content.controls[1].content.controls[0]
    spans_space = "".join(span.text for span in first_diag_space.spans)
    assert "5' " in spans_space
    assert " 3'" in spans_space

    # Test Dash separator
    settings_dash = GUISettings()
    settings_dash["dimer_sequence_separator"] = "Dash ('-')"
    card_dash = Dismissible2DCard(
        "card_ds", step, settings_dash, dismiss_callback=MagicMock()
    )
    first_sub_dash = card_dash.dimer_subcontainers.controls[0]
    first_diag_dash = first_sub_dash.content.controls[1].content.controls[0]
    spans_dash = "".join(span.text for span in first_diag_dash.spans)
    assert "5'-" in spans_dash
    assert "-3'" in spans_dash


def test_designer_2d_view_template_circular_affects_results() -> None:
    """Test that input_data.template_circular changes 2D designer results."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.run_task = _sync_run_task

    fwd_seq = "ATGCATGCAT"
    rev_seq = "CGTACGTACG"
    template_str = "A" * 10 + rev_seq + "C" * 20 + fwd_seq + "T" * 10

    # 1. Linear template (template_circular = False) with max_amplicons = 5
    input_linear = GUIInput()
    input_linear.template = template_str
    input_linear.template_circular = False

    view_linear = Designer2DView(mock_page, input_linear, GUISettings())
    view_linear.form.fwd_dna_input.value = fwd_seq
    view_linear.form.fwd_min_len_input.value = "10"
    view_linear.form.rev_dna_input.value = rev_seq
    view_linear.form.rev_min_len_input.value = "10"
    view_linear.form.filter_dna_checkbox.value = True
    view_linear.form.max_amplicons_input.value = "5"

    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        view_linear._run_designer_event()

    assert view_linear._cached_designer is not None
    assert len(view_linear._cached_designer.all_steps) == 1
    assert view_linear._cached_designer.all_steps[0].amplicon_count == 4

    # 2. Circular template (template_circular = True) with max_amplicons = 5
    input_circ = GUIInput()
    input_circ.template = template_str
    input_circ.template_circular = True

    view_circ = Designer2DView(mock_page, input_circ, GUISettings())
    view_circ.form.fwd_dna_input.value = fwd_seq
    view_circ.form.fwd_min_len_input.value = "10"
    view_circ.form.rev_dna_input.value = rev_seq
    view_circ.form.rev_min_len_input.value = "10"
    view_circ.form.filter_dna_checkbox.value = True
    view_circ.form.max_amplicons_input.value = "5"

    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        view_circ._run_designer_event()

    assert view_circ._cached_designer is not None
    # 16 amplicons on circular template > max_amplicons of 5 -> filtered out
    assert len(view_circ._cached_designer.all_steps) == 0
