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

from unittest.mock import MagicMock

import flet as ft
import pytest

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

    # Empty forward sequence
    form.fwd_dna_input.value = ""
    with pytest.raises(
        ValueError, match="Forward DNA sequence cannot be empty"
    ):
        form.validate_and_get_params()

    # Forward min length exceeds sequence length
    form.fwd_dna_input.value = "ATGC"
    form.fwd_min_len_input.value = "10"
    with pytest.raises(
        ValueError, match="Forward min length exceeds sequence length"
    ):
        form.validate_and_get_params()

    # Invalid quality filter
    form.fwd_dna_input.value = "ATGCGTACGT"
    form.fwd_min_len_input.value = "8"
    form.rev_dna_input.value = "CGTACGATGC"
    form.rev_min_len_input.value = "8"
    form.quality_filter_input.value = "-5.0"
    with pytest.raises(ValueError, match="Quality filter must be non-negative"):
        form.validate_and_get_params()


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
