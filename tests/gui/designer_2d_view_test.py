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
from amplifyp.gui.colours import GUIColours
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

    async def _load_and_drain() -> None:
        loop = asyncio.get_running_loop()
        tasks: list[asyncio.Task[Any]] = []

        def _loop_run_task(func: Any, *args: Any) -> asyncio.Task[Any]:
            task = loop.create_task(func(*args))
            tasks.append(task)
            return task

        mock_page.run_task = _loop_run_task
        with (
            patch(
                "amplifyp.gui.utils.data_helpers.pick_and_read_file",
                new=AsyncMock(side_effect=mock_pick_and_read_file),
            ),
            patch(
                "amplifyp.gui.views.designer.designer_view_base."
                "threading.Thread",
                side_effect=lambda target, daemon: MagicMock(start=target),
            ),
        ):
            await view._load_designer_2d_click(MagicMock())
            await asyncio.gather(*tasks)

    asyncio.run(_load_and_drain())

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

    # Loading while an analysis is already running must notify, skip file
    # picking and form mutation, and not abort the active analysis.
    active_cancel = threading.Event()
    view._analysis_running = True
    view._cancel_event = active_cancel
    view.form.fwd_dna_input.value = "AAAAAAAAAA"
    with (
        patch(
            "amplifyp.gui.utils.data_helpers.pick_and_read_file",
            new=AsyncMock(side_effect=mock_pick_and_read_file),
        ) as mock_pick,
        patch.object(view, "_show_notification") as mock_notify,
    ):
        asyncio.run(view._load_designer_2d_click(MagicMock()))
        mock_pick.assert_not_awaited()
        mock_notify.assert_called_once_with(
            "Cannot load parameters while an analysis is running."
        )
    assert active_cancel.is_set() is False
    assert view.form.fwd_dna_input.value == "AAAAAAAAAA"

    # Starting an analysis while the file picker is open must also prevent
    # overwriting the form or aborting the newly started run.
    async def _pick_starts_analysis(*args: Any, **kwargs: Any) -> str:
        view._analysis_running = True
        view._cancel_event = active_cancel
        return saved_content

    view._analysis_running = False
    view._cancel_event = None
    with (
        patch(
            "amplifyp.gui.utils.data_helpers.pick_and_read_file",
            new=AsyncMock(side_effect=_pick_starts_analysis),
        ),
        patch.object(view, "_show_notification") as mock_notify,
    ):
        asyncio.run(view._load_designer_2d_click(MagicMock()))
        mock_notify.assert_called_once_with(
            "Cannot load parameters while an analysis is running."
        )
    assert active_cancel.is_set() is False
    assert view.form.fwd_dna_input.value == "AAAAAAAAAA"
    view._analysis_running = False
    view._cancel_event = None


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

    # 8. Card helpers: format_primer_properties range, exception, and copy
    r_deg = settings.calculate_primer_tm_range(Primer("ACGTNNACGT"))
    tm_deg, _ = format_primer_properties(Primer("ACGTNNACGT"), settings)
    assert tm_deg.startswith(f"Tm: {r_deg.max_tm:.1f}°C (")
    assert " to " in tm_deg
    r_approx = settings.calculate_primer_tm_range(Primer("N" * 10))
    tm_approx, _ = format_primer_properties(Primer("N" * 10), settings)
    assert tm_approx.startswith(f"Tm: {r_approx.max_tm:.1f}°C (~")

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
    # flush coroutine on any later show_loading call. Restore _sync_run_task.
    mock_page.run_task = _sync_run_task

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


def test_designer_2d_send_primer_to_input_buttons() -> None:
    """Test send-to-input buttons sit right of copy buttons and update input."""
    from amplifyp.dna import Primer
    from amplifyp.gui.views.designer.designer_card_helpers import (
        build_primer_summary_row,
    )

    mock_page = MagicMock(spec=ft.Page)
    mock_page.run_task = _sync_run_task
    input_data = GUIInput(primers=[{"name": "", "seq": "", "active": False}])
    settings = GUISettings()

    # 1. Helper puts send button right of copy button when callback provided
    sent: list[str] = []
    row_col = build_primer_summary_row(
        "Forward Primer",
        Primer("ATGCGTAC"),
        settings,
        on_send_to_input=lambda: sent.append("clicked"),
    )
    action_row = cast(ft.Row, row_col.controls[1])
    assert len(action_row.controls) == 3
    copy_btn = cast(ft.IconButton, action_row.controls[1])
    send_btn = cast(ft.IconButton, action_row.controls[2])
    assert copy_btn.icon == ft.Icons.COPY
    assert send_btn.icon == ft.Icons.PLAYLIST_ADD
    assert send_btn.on_click is not None
    send_btn.on_click(MagicMock())
    assert sent == ["clicked"]

    # 2. Designer2DView wires both Fwd and Rev send buttons cleanly
    view = Designer2DView(mock_page, input_data, settings)
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

    # Trigger settings change and ensure send buttons remain intact
    card.on_settings_change()

    details_col = cast(ft.Column, card.primer_info_section.content)
    fwd_summary_col = cast(ft.Column, details_col.controls[0])
    rev_summary_col = cast(ft.Column, details_col.controls[1])
    fwd_actions = cast(ft.Row, fwd_summary_col.controls[1])
    rev_actions = cast(ft.Row, rev_summary_col.controls[1])
    assert len(fwd_actions.controls) == 3
    assert len(rev_actions.controls) == 3

    fwd_send_btn = cast(ft.IconButton, fwd_actions.controls[2])
    rev_send_btn = cast(ft.IconButton, rev_actions.controls[2])
    assert fwd_send_btn.icon == ft.Icons.PLAYLIST_ADD
    assert rev_send_btn.icon == ft.Icons.PLAYLIST_ADD
    assert fwd_send_btn.on_click is not None
    assert rev_send_btn.on_click is not None

    with patch.object(view, "_show_notification") as mock_notify:
        fwd_send_btn.on_click(MagicMock())
        rev_send_btn.on_click(MagicMock())
        # Re-clicking existing primer notifies without duplicating
        fwd_send_btn.on_click(MagicMock())
        assert mock_notify.call_count == 3
        assert "Added" in mock_notify.call_args_list[0][0][0]
        assert "already in the Input primer list" in mock_notify.call_args[0][0]

    assert len(input_data.primers) == 2
    assert input_data.primers[0]["seq"] == step.fwd_fwd.primer_1.seq
    assert input_data.primers[1]["seq"] == step.rev_rev.primer_1.seq
    assert all(p["active"] is True for p in input_data.primers)

    # Empty sequence shows error dialog
    with patch(
        "amplifyp.gui.views.designer.designer_view_base.show_error_dialog"
    ) as mock_err:
        view._handle_add_primer("Empty", "   ")
        mock_err.assert_called_once()


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

    # Detached page + no running loop: callback discarded without creating
    # or executing a coroutine on a temporary worker-thread loop.
    called = False
    executed = False

    def _discarded_factory() -> Any:
        nonlocal called
        called = True

        async def _inner() -> None:
            nonlocal executed
            executed = True

        return _inner()

    assert view._schedule_on_event_loop(_discarded_factory) is False
    assert called is False
    assert executed is False
    assert len(view._background_tasks) == 0

    async def _discarded() -> None:
        nonlocal executed
        executed = True

    # Running loop on the current thread: task created on that loop.
    async def _drive() -> None:
        nonlocal executed
        assert view._schedule_on_event_loop(_discarded) is True
        assert len(view._background_tasks) == 1
        await asyncio.sleep(0)
        assert executed is True
        # Allow the done callback scheduled on the next tick to discard it.
        await asyncio.sleep(0)
        assert len(view._background_tasks) == 0

        # Unhandled exception in fallback task is retrieved and logged.
        async def _boom() -> None:
            raise RuntimeError("fallback boom")

        with patch(
            "amplifyp.gui.views.designer.designer_view_base.logger.error"
        ) as mock_log:
            assert view._schedule_on_event_loop(_boom) is True
            assert len(view._background_tasks) == 1
            await asyncio.sleep(0)
            await asyncio.sleep(0)
            assert len(view._background_tasks) == 0
            mock_log.assert_called_once()

        # Cancelled fallback task is discarded cleanly without error logs.
        async def _sleeper() -> None:
            await asyncio.sleep(10)

        with patch(
            "amplifyp.gui.views.designer.designer_view_base.logger.error"
        ) as mock_log:
            assert view._schedule_on_event_loop(_sleeper) is True
            task = next(iter(view._background_tasks))
            task.cancel()
            await asyncio.sleep(0)
            await asyncio.sleep(0)
            assert len(view._background_tasks) == 0
            mock_log.assert_not_called()

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
    settings_space["sequence_separator"] = "Space (' ')"
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
    settings_dash["sequence_separator"] = "Dash ('-')"
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


def test_designer_2d_view_update_ui_refreshes_theme_and_cards() -> None:
    """Test Designer2DView.update_ui refreshes form borders, grid, and cards."""
    from amplifyp.gui.colours import GUIColours
    from amplifyp.gui.views.designer_2d.dismissible_2d_card import (
        Dismissible2DCard,
    )

    mock_page = MagicMock(spec=ft.Page)
    mock_page.run_task = _sync_run_task
    settings = GUISettings()
    view = Designer2DView(mock_page, GUIInput(), settings)

    prev_dark = GUIColours.dark_mode
    try:
        # 1. Empty state placeholder refreshes colour cleanly on theme change
        GUIColours.dark_mode = True
        view.update_ui()
        placeholder = view.results_grid.content_column.controls[1].content
        assert isinstance(placeholder, ft.Text)
        assert placeholder.color == GUIColours.TEXT_ON_SURFACE
        assert view.top_left_container.border == ft.Border.all(
            1, GUIColours.OUTLINE_VARIANT
        )
        assert view.form.fwd_dna_input.border_color == GUIColours.OUTLINE
        assert view.form.filter_dna_checkbox.border == ft.Border.all(
            1, GUIColours.OUTLINE
        )
        assert view.form.error_text.color == GUIColours.ERROR_RED

        # 2. Populate results + open card in Light mode, then switch to Dark,
        # new font_family, Tm method, and separator
        GUIColours.dark_mode = False
        view.form.fwd_dna_input.value = "ATGCGTACGT"
        view.form.fwd_min_len_input.value = "10"
        view.form.rev_dna_input.value = "CGTACGATGC"
        view.form.rev_min_len_input.value = "10"
        with patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.threading.Thread",
            side_effect=lambda target, daemon: MagicMock(start=target),
        ):
            view._run_designer_event()
        assert view._cached_designer is not None
        view._on_grid_step_selected(view._cached_designer.all_steps[0])
        assert len(view._active_cards) == 1
        card = cast(Dismissible2DCard, view._active_cards[0])

        fwd_summary_col = card.primer_info_section.content.controls[0]
        old_tm_text = (
            fwd_summary_col.controls[0].controls[1].controls[0].content.value
        )

        GUIColours.dark_mode = True
        settings["font_family"] = "Consolas"
        settings["tm_method"] = "Lander / Amplify 4"
        settings["sequence_separator"] = "Dash ('-')"
        view.update_ui()

        assert card.font_family == "Consolas"
        assert card.primer_info_section.bgcolor == GUIColours.SURFACE_VARIANT
        updated_fwd_summary = card.primer_info_section.content.controls[0]
        new_tm_text = (
            updated_fwd_summary.controls[0]
            .controls[1]
            .controls[0]
            .content.value
        )
        assert new_tm_text != old_tm_text
        seq_field = updated_fwd_summary.controls[1].controls[0]
        assert isinstance(seq_field, ft.TextField)
        assert seq_field.border_color == GUIColours.OUTLINE
        assert seq_field.text_style is not None
        assert seq_field.text_style.font_family == "Consolas"

        first_diag = (
            card.dimer_subcontainers.controls[0]
            .content.controls[1]
            .content.controls[0]
        )
        spans_text = "".join(span.text for span in first_diag.spans)
        assert "5'-" in spans_text
    finally:
        GUIColours.dark_mode = prev_dark


def test_grid_2d_results_cell_click_o1_selection_toggle() -> None:
    """Verify O(1) cell selection restores prior cell borders accurately."""
    from amplifyp.gui.views.designer_2d.grid_2d_results import Grid2DResultsView

    mock_page = MagicMock(spec=ft.Page)
    settings = GUISettings()
    designer = PrimerDesigner2D(
        DNA("ATGCGTACGT"),
        fwd_min_length=9,
        rev_dna=DNA("CGTACGTACG"),
        rev_min_length=9,
    )
    selected: list[Any] = []
    grid = Grid2DResultsView(
        settings=settings,
        page=mock_page,
        on_select_step_callback=selected.append,
    )
    grid.update_grid(designer)
    keys = list(grid._cell_containers.keys())
    assert len(keys) >= 2
    key_a, key_b = keys[0], keys[1]
    step_a, step_b = designer.all_steps[0], designer.all_steps[1]

    grid._on_cell_click(step_a, key_a)
    assert grid._selected_cell_key == key_a
    assert grid._cell_containers[key_a].border == ft.Border.all(
        3, GUIColours.PRIMARY
    )

    grid._on_cell_click(step_b, key_b)
    assert grid._selected_cell_key == key_b
    assert grid._cell_containers[key_b].border == ft.Border.all(
        3, GUIColours.PRIMARY
    )
    expected_a_border = (
        ft.Border.all(2, GUIColours.SUCCESS_GREEN)
        if key_a in grid._best_cell_keys
        else ft.Border.all(1, GUIColours.OUTLINE_VARIANT)
    )
    assert grid._cell_containers[key_a].border == expected_a_border
    assert len(selected) == 2

    # Re-clicking the already-selected cell keeps its selected border intact
    grid._on_cell_click(step_b, key_b)
    assert grid._selected_cell_key == key_b
    assert grid._cell_containers[key_b].border == ft.Border.all(
        3, GUIColours.PRIMARY
    )
    assert len(selected) == 3

    # Clearing resets tracked selection key cleanly
    grid.clear_grid()
    assert grid._selected_cell_key is None


def test_designer_2d_runs_cooperatively_when_threads_unavailable() -> None:
    """Test 2D designer runs on the event loop in Pyodide / threadless mode."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.run_task = _sync_run_task
    view = Designer2DView(mock_page, GUIInput(), GUISettings())
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGATGC"
    view.form.rev_min_len_input.value = "8"

    # 1. Pyodide/Emscripten runtime skips threading.Thread completely.
    with (
        patch(
            "amplifyp.gui.utils.system.can_use_background_threads",
            return_value=False,
        ),
        patch(
            "amplifyp.gui.views.designer.designer_view_base.threading.Thread",
            side_effect=AssertionError("Thread should not be created"),
        ),
    ):
        view._run_designer_event()

    assert view._cached_designer is not None
    assert len(view._cached_designer.all_steps) == 9
    assert view._analysis_running is False
    assert view._cancel_event is None
    assert view.form.analyse_button.content == "Analyse"

    # 2. Fallback when Thread.start() raises "can't start new thread"
    view._cached_designer = None

    def _failing_thread(*_args: Any, **_kwargs: Any) -> Any:
        def _raise() -> None:
            raise RuntimeError("can't start new thread")

        return type("FailingThread", (), {"start": staticmethod(_raise)})()

    with patch(
        "amplifyp.gui.views.designer.designer_view_base.threading.Thread",
        side_effect=_failing_thread,
    ):
        view._run_designer_event()

    assert view._cached_designer is not None
    assert len(view._cached_designer.all_steps) == 9
    assert view._analysis_running is False

    # 3. Error handling path in cooperative async worker
    with (
        patch(
            "amplifyp.gui.utils.system.can_use_background_threads",
            return_value=False,
        ),
        patch.object(
            PrimerDesigner2D,
            "create_async",
            new=AsyncMock(side_effect=RuntimeError("Async 2D failure")),
        ),
        patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.show_error_dialog"
        ) as mock_err,
    ):
        view._run_designer_event()
        mock_err.assert_called_once()
        assert view._analysis_running is False
        assert view._cancel_event is None

    # 4. Mid-run abort via Abort button click during cooperative async run
    orig_create_async = PrimerDesigner2D.create_async

    async def _aborting_create_async(*args: Any, **kwargs: Any) -> Any:
        orig_progress = kwargs.get("on_progress")

        def _progress(done: int, total: int) -> None:
            if orig_progress is not None:
                orig_progress(done, total)
            if done == 1:
                view._run_designer_event()

        kwargs["on_progress"] = _progress
        kwargs["yield_interval"] = 0.0
        return await orig_create_async(*args, **kwargs)

    with (
        patch(
            "amplifyp.gui.utils.system.can_use_background_threads",
            return_value=False,
        ),
        patch.object(
            PrimerDesigner2D,
            "create_async",
            side_effect=_aborting_create_async,
        ),
    ):
        view._run_designer_event()

    assert view._cached_designer is not None
    assert view._cached_designer.aborted is True
    assert len(view._cached_designer.all_steps) == 1
    assert view.results_grid.progress_tracker.is_active is False
    assert view.results_grid.progress_tracker.loading_body is None
    assert view._analysis_running is False
    assert view._cancel_event is None

    # 5. asyncio.CancelledError cleans up ProgressTracker even after 100% tick
    async def _cancel_after_final_tick(*_args: Any, **kwargs: Any) -> Any:
        on_progress = kwargs.get("on_progress")
        if on_progress is not None:
            on_progress(9, 9)
        raise asyncio.CancelledError()

    with (
        patch(
            "amplifyp.gui.utils.system.can_use_background_threads",
            return_value=False,
        ),
        patch.object(
            PrimerDesigner2D,
            "create_async",
            side_effect=_cancel_after_final_tick,
        ),
    ):
        with pytest.raises(asyncio.CancelledError):
            view._run_designer_event()

    assert view.results_grid.progress_tracker.is_active is False
    assert view.results_grid.progress_tracker.loading_body is None
    assert view._analysis_running is False
    assert view._cancel_event is None
    assert view.form.analyse_button.content == "Analyse"

    # 6. Re-entry guard when running without a cancel event
    view._analysis_running = True
    view._cancel_event = None
    with patch.object(view.form, "validate_and_get_params") as mock_validate:
        view._run_designer_event()
        mock_validate.assert_not_called()
    view._analysis_running = False

    # 7. Detached page + no running loop discards async task and resets view
    detached_page = MagicMock(spec=ft.Page)
    detached_page.run_task.side_effect = RuntimeError("detached")
    detached_view = Designer2DView(detached_page, GUIInput(), GUISettings())
    detached_view.form.fwd_dna_input.value = "ATGCGTACGT"
    detached_view.form.fwd_min_len_input.value = "8"
    detached_view.form.rev_dna_input.value = "CGTACGATGC"
    detached_view.form.rev_min_len_input.value = "8"
    with patch(
        "amplifyp.gui.utils.system.can_use_background_threads",
        return_value=False,
    ):
        detached_view._run_designer_event()

    assert detached_view._analysis_running is False
    assert detached_view._cancel_event is None
    assert detached_view.results_grid.progress_tracker.is_active is False
    assert detached_view.results_grid.progress_tracker.loading_body is None
    assert detached_view.form.analyse_button.content == "Analyse"

    # 8. Error in _on_analysis_success does not trigger _on_analysis_error
    with (
        patch(
            "amplifyp.gui.utils.system.can_use_background_threads",
            return_value=False,
        ),
        patch.object(
            view,
            "_on_analysis_success",
            new=AsyncMock(side_effect=ValueError("UI render failure")),
        ),
        patch(
            "amplifyp.gui.views.designer_2d.designer_2d_view.show_error_dialog"
        ) as mock_err,
    ):
        with pytest.raises(ValueError, match="UI render failure"):
            view._run_designer_event()
        mock_err.assert_not_called()
        assert view._analysis_running is False
        assert view._cancel_event is None
        assert view.results_grid.progress_tracker.loading_body is None
        assert view.form.analyse_button.content == "Analyse"

    # 8b. Mid-render error inside Grid2DResultsView.update_grid still cleans
    # up the loading bar in _finish_analysis_ui.
    with (
        patch(
            "amplifyp.gui.utils.system.can_use_background_threads",
            return_value=False,
        ),
        patch(
            "amplifyp.gui.views.designer_2d.grid_2d_results.designer_2d_colour",
            side_effect=ValueError("colour map boom"),
        ),
        patch.object(
            view.results_grid,
            "clear_grid",
            wraps=view.results_grid.clear_grid,
        ) as mock_clear,
    ):
        with pytest.raises(ValueError, match="colour map boom"):
            view._run_designer_event()
        mock_clear.assert_called_once()
        assert view.results_grid.progress_tracker.loading_body is None
        assert view._analysis_running is False
        assert view._cancel_event is None

    # 9. Shared event loop interleaves _flush_task and concurrent Abort click
    async def _run_shared_loop() -> None:
        loop = asyncio.get_running_loop()
        scheduled_tasks: list[asyncio.Task[Any]] = []
        flush_updates = 0
        started = False
        abort_ready = asyncio.Event()

        loop_page = MagicMock(spec=ft.Page)

        def _loop_run_task(func: Any, *args: Any) -> asyncio.Task[Any]:
            task = loop.create_task(func(*args))
            scheduled_tasks.append(task)
            return task

        loop_page.run_task = _loop_run_task

        loop_view = Designer2DView(loop_page, GUIInput(), GUISettings())
        loop_view.form.fwd_dna_input.value = "ATGCGTACGT"
        loop_view.form.fwd_min_len_input.value = "8"
        loop_view.form.rev_dna_input.value = "CGTACGATGC"
        loop_view.form.rev_min_len_input.value = "8"

        def _count_flush_update() -> None:
            nonlocal flush_updates
            if (
                started
                and loop_view._analysis_running
                and loop_view.results_grid.progress_tracker.loading_body
                is not None
            ):
                flush_updates += 1

        loop_page.update.side_effect = _count_flush_update

        async def _yielding_create_async(*args: Any, **kwargs: Any) -> Any:
            orig_progress = kwargs.get("on_progress")

            def _progress(done: int, total: int) -> None:
                if orig_progress is not None:
                    orig_progress(done, total)
                if done == 2:
                    abort_ready.set()

            kwargs["on_progress"] = _progress
            kwargs["yield_interval"] = 0.0
            return await orig_create_async(*args, **kwargs)

        async def _concurrent_abort() -> None:
            try:
                await asyncio.wait_for(abort_ready.wait(), timeout=5)
            except TimeoutError as exc:
                raise AssertionError(
                    "progress tick done == 2 was not reported"
                ) from exc
            loop_view._run_designer_event()

        with (
            patch(
                "amplifyp.gui.utils.system.can_use_background_threads",
                return_value=False,
            ),
            patch(
                "amplifyp.gui.views.designer.progress_tracker."
                "_FLUSH_INTERVAL_SECONDS",
                0.001,
            ),
            patch.object(
                PrimerDesigner2D,
                "create_async",
                side_effect=_yielding_create_async,
            ),
        ):
            abort_task = loop.create_task(_concurrent_abort())
            loop_view._run_designer_event()
            started = True
            try:
                await asyncio.wait_for(
                    asyncio.gather(abort_task, *scheduled_tasks),
                    timeout=30,
                )
            except BaseException:
                abort_task.cancel()
                for task in scheduled_tasks:
                    task.cancel()
                await asyncio.gather(
                    abort_task, *scheduled_tasks, return_exceptions=True
                )
                raise

        assert loop_view._cached_designer is not None
        assert loop_view._cached_designer.aborted is True
        assert 2 <= len(loop_view._cached_designer.all_steps) < 9
        assert flush_updates >= 1
        assert loop_view._analysis_running is False
        assert loop_view._cancel_event is None
        assert loop_view.results_grid.progress_tracker.loading_body is None

    asyncio.run(_run_shared_loop())


def test_designer_2d_clear_all_guard_and_stale_callbacks() -> None:
    """Test Clear All guard, stale callback ignoring, and threaded teardown."""
    from amplifyp.gui.views.designer.designer_view_base import BaseDesignerView

    mock_page = MagicMock(spec=ft.Page)
    mock_page.run_task = _sync_run_task
    view = Designer2DView(mock_page, GUIInput(), GUISettings())
    view.form.fwd_dna_input.value = "ATGCGTACGT"
    view.form.fwd_min_len_input.value = "8"
    view.form.rev_dna_input.value = "CGTACGATGC"
    view.form.rev_min_len_input.value = "8"

    # 1. _clear_all while analysis is running notifies and does not clear form.
    view._analysis_running = True
    with patch.object(view, "_show_notification") as mock_notify:
        view._clear_all()
        mock_notify.assert_called_once_with(
            "Cannot clear parameters while an analysis is running."
        )
    assert view.form.fwd_dna_input.value == "ATGCGTACGT"
    view._analysis_running = False

    # 2. Stale cancel_event callbacks in _on_analysis_success / _error are
    # ignored without overwriting _cached_designer or showing error dialogues.
    active_event = threading.Event()
    stale_event = threading.Event()
    view._cancel_event = active_event
    view._cached_designer = None
    with patch(
        "amplifyp.gui.views.designer_2d.designer_2d_view.show_error_dialog"
    ) as mock_err:
        asyncio.run(view._on_analysis_success(MagicMock(), stale_event))
        asyncio.run(view._on_analysis_error(RuntimeError("stale"), stale_event))
        mock_err.assert_not_called()
    assert view._cached_designer is None
    view._cancel_event = None

    # 3. Threaded worker keeps _analysis_running True until
    # _on_analysis_finished executes on the event loop, and clears it if
    # scheduling fails (detached).
    queued: list[tuple[Any, tuple[Any, ...]]] = []

    def _queue_task(func: Any, *args: Any) -> None:
        if getattr(func, "__name__", "") == "_flush_task":
            return
        queued.append((func, args))

    mock_page.run_task = _queue_task
    with patch(
        "amplifyp.gui.views.designer.designer_view_base.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        view._run_designer_event()
    # Worker finished, but queued event-loop callbacks have not run yet.
    assert view._analysis_running is True
    assert view._cached_designer is None
    for fn, args in queued:
        asyncio.run(fn(*args))
    assert view._analysis_running is False
    assert view._cached_designer is not None

    # Detached page in threaded worker clears _analysis_running via fallback.
    detached_page = MagicMock(spec=ft.Page)
    detached_page.run_task.side_effect = RuntimeError("detached")
    detached_view = Designer2DView(detached_page, GUIInput(), GUISettings())
    detached_view.form.fwd_dna_input.value = "ATGCGTACGT"
    detached_view.form.fwd_min_len_input.value = "8"
    detached_view.form.rev_dna_input.value = "CGTACGATGC"
    detached_view.form.rev_min_len_input.value = "8"
    with patch(
        "amplifyp.gui.views.designer.designer_view_base.threading.Thread",
        side_effect=lambda target, daemon: MagicMock(start=target),
    ):
        detached_view._run_designer_event()
    assert detached_view._analysis_running is False

    # 4. BaseDesignerView without a form handles abort/finish helpers cleanly.
    bare = BaseDesignerView(mock_page, GUIInput(), GUISettings())
    bare._analysis_running = True
    bare._set_button_abort_mode(True)
    bare._finish_analysis_ui()
    assert bare._analysis_running is False


def test_designer_2d_scroll_and_error_clearing() -> None:
    """Test form scrollability and stale error clearing."""
    mock_page = MagicMock(spec=ft.Page)
    view = Designer2DView(mock_page, GUIInput(), GUISettings())
    form = view.form

    assert view.top_left_container.height == 425
    assert form.scroll == ft.ScrollMode.AUTO

    # 1. Lengthening forward/reverse DNA clears stale min-length overflow errors
    form.fwd_dna_input.value = "ATGC"
    form.fwd_min_len_input.value = "10"
    form.rev_dna_input.value = "CGTA"
    form.rev_min_len_input.value = "10"
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        with pytest.raises(ValueError):
            form.validate_and_get_params()
        assert "Exceeds sequence length" in (form.fwd_min_len_input.error or "")
        assert "Exceeds sequence length" in (form.rev_min_len_input.error or "")

        form.fwd_dna_input.value = "ATGCGTACGTAT"
        form._on_fwd_dna_change(
            MagicMock(spec=ft.ControlEvent, control=form.fwd_dna_input)
        )
        assert form.fwd_min_len_input.error is None

        # Rev Comp also clears rev_min_len_input error if updated seq covers it
        form.rev_dna_input.value = "ATGCGTACGTAT"
        form._on_reverse_complement_click(None)
        assert form.rev_min_len_input.error is None

        # 2. Unrelated min-length errors (empty / non-positive / invalid) remain
        form.fwd_min_len_input.value = "0"
        form.fwd_min_len_input.error = "Must be > 0"
        form.fwd_dna_input.value = "ATGCGTACGTATGC"
        form._on_fwd_dna_change(
            MagicMock(spec=ft.ControlEvent, control=form.fwd_dna_input)
        )
        assert form.fwd_min_len_input.error == "Must be > 0"

        form.rev_min_len_input.value = "abc"
        form.rev_min_len_input.error = "Must be > 0"
        form.rev_dna_input.value = "ATGCGTACGTATGC"
        form._on_rev_dna_change(
            MagicMock(spec=ft.ControlEvent, control=form.rev_dna_input)
        )
        assert form.rev_min_len_input.error == "Must be > 0"

        # 3. Unchecking template filter clears error banner and max field error
        form.filter_dna_checkbox.value = True
        form.max_amplicons_input.error = "Must be > 0"
        form.show_error("No template sequence available for amplicon check.")
        form.filter_dna_checkbox.value = False
        form._on_filter_dna_change(
            MagicMock(spec=ft.ControlEvent, control=form.filter_dna_checkbox)
        )
        assert form.max_amplicons_input.error is None
        assert form.error_text.visible is False
        assert form.error_text.value == ""

    # 4. Detached page RuntimeError in _update_page is ignored cleanly
    mock_page.update.side_effect = RuntimeError("detached")
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        form.fwd_min_len_input.value = "8"
        form.fwd_min_len_input.error = "Exceeds sequence length (4)"
        form.fwd_dna_input.value = "ATGCGTACGT"
        form._on_fwd_dna_change(
            MagicMock(spec=ft.ControlEvent, control=form.fwd_dna_input)
        )
        assert form.fwd_min_len_input.error is None
