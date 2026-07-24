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
from typing import Any
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
    assert view.form.fwd_dna_input.value == ""
    assert view.form.fwd_min_len_input.value == "18"
    assert view.form.rev_dna_input.value == ""
    assert view.form.rev_min_len_input.value == "18"
    assert view.form.filter_metric_dropdown.value == "MAX"
    assert len(view.right_cards_list.controls) == 0


def test_designer_2d_form_validation_success() -> None:
    """Test successful input validation in Designer2DForm."""
    settings = GUISettings()
    form = Designer2DForm(settings=settings, on_submit_callback=lambda: None)

    form.fwd_dna_input.value = "ATGCGTACGT"
    form.fwd_min_len_input.value = "8"
    form.rev_dna_input.value = "CGTACGATGC"
    form.rev_min_len_input.value = "8"
    form.quality_filter_input.value = "50.0"
    form.overlap_filter_input.value = "5"
    form.filter_metric_dropdown.value = "MEAN"

    (
        fwd_dna,
        fwd_min_len,
        rev_dna,
        rev_min_len,
        threshold,
        max_overlap,
        filter_metric,
    ) = form.validate_and_get_params()

    assert fwd_dna.seq_upper == "ATGCGTACGT"
    assert fwd_min_len == 8
    assert rev_dna.seq_upper == "CGTACGATGC"
    assert rev_min_len == 8
    assert threshold == 50.0
    assert max_overlap == 5
    assert filter_metric == FilterMetric.MEAN


def test_designer_2d_form_validation_errors() -> None:
    """Test validation failure cases in Designer2DForm."""
    settings = GUISettings()
    form = Designer2DForm(settings=settings, on_submit_callback=lambda: None)

    # Empty forward sequence and empty reverse sequence simultaneously
    form.fwd_dna_input.value = ""
    form.rev_dna_input.value = ""
    with pytest.raises(ValueError, match="Input validation failed"):
        form.validate_and_get_params()
    assert form.fwd_dna_input.error == "Forward DNA sequence cannot be empty"
    assert form.rev_dna_input.error == "Reverse DNA sequence cannot be empty"

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
    form.quality_filter_input.value = "-5.0"
    with pytest.raises(ValueError, match="Input validation failed"):
        form.validate_and_get_params()
    assert form.quality_filter_input.error == "Must be >= 0"


def test_designer_2d_view_run_analysis_and_grid() -> None:
    """Test running 2D analysis populates grid and allows card creation."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGATGC"
    view.form.rev_min_len_input.value = "8"
    view.form.quality_filter_input.value = ""
    view.form.overlap_filter_input.value = ""
    view.form.filter_metric_dropdown.value = "MAX"

    view._run_designer_event()

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
        "Forward: 10 bp, Reverse: 10 bp" in card._card_id
        or card.step == step_10_10
    )


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
    view.form.quality_filter_input.value = ""
    view.form.overlap_filter_input.value = ""
    view.form.filter_metric_dropdown.value = "MEAN"

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
    assert parsed["quality_filter"] == ""
    assert parsed["overlap_filter"] == ""
    assert parsed["filter_metric"] == "MEAN"

    # 2. Test Load
    # Reset form to different values
    view.form.fwd_dna_input.value = "CGT"
    view.form.fwd_min_len_input.value = "10"
    view.form.rev_dna_input.value = "ATG"
    view.form.rev_min_len_input.value = "10"
    view.form.quality_filter_input.value = "60.0"
    view.form.overlap_filter_input.value = "3"
    view.form.filter_metric_dropdown.value = "MAX"

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
    assert view.form.fwd_min_len_input.value == "8"
    assert view.form.rev_dna_input.value == "CGTACGATGC"
    assert view.form.rev_min_len_input.value == "8"
    assert view.form.quality_filter_input.value == ""
    assert view.form.overlap_filter_input.value == ""
    assert view.form.filter_metric_dropdown.value == "MEAN"

    # Verify that the analysis automatically ran (should populate
    # _cached_designer)
    assert view._cached_designer is not None
    assert len(view._cached_designer) == 9
