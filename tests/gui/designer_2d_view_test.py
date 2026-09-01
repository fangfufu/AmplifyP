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
    assert view.form.fwd_min_len_input.value == ""
    assert view.form.rev_dna_input.value == ""
    assert view.form.rev_min_len_input.value == ""
    assert view.form.max_quality_input.value == ""
    assert view.form.max_overlap_input.value == ""
    assert view.form.clear_all_button is not None
    assert len(view.right_cards_list.controls) == 0


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
    ) = form.validate_and_get_params()

    assert fwd_dna.seq_upper == "ATGCGTACGT"
    assert fwd_min_len == 8
    assert rev_dna.seq_upper == "CGTACGATGC"
    assert rev_min_len == 8
    assert threshold == 50.0
    assert max_overlap == 5
    assert filter_metric == FilterMetric.MAX


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
    view.form.max_quality_input.value = ""
    view.form.max_overlap_input.value = ""

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
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    settings = GUISettings()

    view = Designer2DView(mock_page, input_data, settings)
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGATGC"
    view.form.rev_min_len_input.value = "8"
    view.form.max_quality_input.value = ""
    view.form.max_overlap_input.value = ""

    view._run_designer_event()
    assert view._cached_designer is not None
    assert len(view.results_grid._cell_containers) == 9

    # Add a card
    step = view._cached_designer.get_step(0)
    view.results_grid._on_cell_click(step, (10, 10))
    assert len(view.right_cards_list.controls) == 1

    # Clear all
    view._clear_all(None)

    assert view.form.fwd_dna_input.value == ""
    assert view.form.fwd_min_len_input.value == ""
    assert view.form.rev_dna_input.value == ""
    assert view.form.rev_min_len_input.value == ""
    assert view.form.max_quality_input.value == ""
    assert view.form.max_overlap_input.value == ""
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

    # 2. Test Load
    # Reset form to different values
    view.form.fwd_dna_input.value = "CGT"
    view.form.fwd_min_len_input.value = "10"
    view.form.rev_dna_input.value = "ATG"
    view.form.rev_min_len_input.value = "10"
    view.form.max_quality_input.value = "60.0"
    view.form.max_overlap_input.value = "3"

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
    assert view.form.max_quality_input.value == ""
    assert view.form.max_overlap_input.value == ""

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
    ):
        view._run_designer_event()
        mock_err_dlg.assert_called_once()

    # 5. Designer2DView UI updates & dismiss with RuntimeError
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

    # 9. Dismissible2DCard: settings change with page update and string boolean
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGATGC"
    view.form.rev_min_len_input.value = "8"
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
