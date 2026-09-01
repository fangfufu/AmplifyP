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

"""Comprehensive unit tests for InputView and subcomponents."""

import asyncio
from typing import Any, cast
from unittest.mock import AsyncMock, MagicMock, patch

import flet as ft
import pytest

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.views.input.events import (
    auto_add_empty_row_if_needed,
    handle_field_blur,
    handle_field_focus,
    handle_field_submit,
)
from amplifyp.gui.views.input.input_view import InputView
from amplifyp.gui.views.input.layout import (
    adjust_template_wrap,
    handle_resize,
    on_pan_end,
    on_pan_update,
)
from amplifyp.gui.views.input.primer.clipboard import (
    copy_primers_click,
    parse_primer_clipboard_text,
    paste_primers_click,
)
from amplifyp.gui.views.input.primer.row import PrimerRow
from amplifyp.gui.views.input.template.casing import (
    change_selection_case,
    lower_case_click,
    upper_case_click,
)
from amplifyp.gui.views.input.template.formatter import (
    adjust_wrap_length,
    update_line_numbers,
    validate_bases_per_line,
)


def _setup_test_view() -> tuple[InputView, MagicMock, GUIInput, GUISettings]:
    """Helper to initialise a mock-enabled InputView."""
    page = MagicMock(spec=ft.Page)
    page.width = 1000.0
    page.height = 600.0
    page.overlay = []

    input_data = GUIInput()
    input_data.template = "ATGCGATCGATC"
    input_data.primers = [
        {"name": "P1", "seq": "ATGC", "active": True},
        {"name": "P2", "seq": "GGCC", "active": True},
    ]

    settings = GUISettings()
    settings["show_primer_temperature"] = True
    settings["tm_colour_scheme"] = "SantaLucia"

    view = InputView(page, input_data, settings=settings)
    view.update = MagicMock()
    view.template_input.update = MagicMock()
    view.template_input.template_sequence.update = MagicMock()
    view.template_input.template_sequence_wrapper.update = MagicMock()
    view.template_input.template_sequence_container.update = MagicMock()
    view.template_input.line_numbers_text.update = MagicMock()
    view.template_input.line_numbers_container.update = MagicMock()

    view.primer_input.update = MagicMock()
    view.primer_input.primer_header.update = MagicMock()
    view.primer_input.primers_list.update = MagicMock()

    return view, page, input_data, settings


def test_input_layout_and_resizing() -> None:
    """Test layout pan and resize handlers."""
    view, _page, _, _ = _setup_test_view()

    # 1. on_pan_update with local delta
    ev_drag = MagicMock()
    ev_drag.local_delta = MagicMock(x=50.0)
    on_pan_update(view, ev_drag)
    assert view.right_fraction > 0.0

    # 2. on_pan_end
    on_pan_end(view, MagicMock())
    assert view.template_input.width is None

    # 3. handle_resize
    handle_resize(view, MagicMock())

    # 4. adjust_template_wrap
    adjust_template_wrap(view, update_first=True)
    adjust_template_wrap(view, update_first=False)


def test_primer_layout_manager() -> None:
    """Test PrimerLayoutManager column sizing and caching."""
    view, _, _, _ = _setup_test_view()
    mgr = view.primer_input.layout_manager

    # 1. get_panel_width
    width = mgr.get_panel_width()
    assert width > 0

    # 2. on_primer_divider_pan
    ev_drag = MagicMock()
    ev_drag.local_delta = MagicMock(x=20.0)
    mgr.on_primer_divider_pan(ev_drag)
    assert mgr.owner.name_column_width >= 80.0

    # 3. on_primer_divider_pan_end
    mgr.on_primer_divider_pan_end(MagicMock())
    assert mgr.owner._visible_rows_cache is None

    # 4. adjust_name_column_width
    mgr.adjust_name_column_width(600.0, during_drag=True)
    mgr.adjust_name_column_width(600.0, during_drag=False)


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_primer_clipboard_operations() -> None:
    """Test primer clipboard parsing, copying, and pasting."""
    view, _, _, _ = _setup_test_view()
    primer_input = view.primer_input

    # 1. parse_primer_clipboard_text
    parsed = parse_primer_clipboard_text(
        "P1\tATGC\nP2,GGCC\nTTTT\nPrimerOnly\n\n"
    )
    assert len(parsed) == 4
    assert parsed[0]["name"] == "P1"
    assert parsed[0]["seq"] == "ATGC"
    assert parsed[2]["seq"] == "TTTT"
    assert parsed[3]["name"] == "PrimerOnly"

    # 2. copy_primers_click with selection
    primer_input.selected_indices = {0}
    with patch("flet.Clipboard.set", new_callable=AsyncMock) as mock_set:
        await copy_primers_click(primer_input, None)
        mock_set.assert_called_once()

    # Fallback to focused index
    primer_input.selected_indices = set()
    primer_input.focused_primer_index = 0
    with patch("flet.Clipboard.set", new_callable=AsyncMock) as mock_set:
        await copy_primers_click(primer_input, None)
        mock_set.assert_called_once()

    # No selection / focus
    primer_input.focused_primer_index = None
    await copy_primers_click(primer_input, None)

    # 3. paste_primers_click
    with patch(
        "flet.Clipboard.get",
        new_callable=AsyncMock,
        side_effect=Exception("Clip fail"),
    ):
        await paste_primers_click(primer_input, None)

    with patch("flet.Clipboard.get", new_callable=AsyncMock, return_value=""):
        await paste_primers_click(primer_input, None)

    with patch(
        "flet.Clipboard.get", new_callable=AsyncMock, return_value="   \n\n"
    ):
        await paste_primers_click(primer_input, None)

    with patch(
        "flet.Clipboard.get", new_callable=AsyncMock, return_value="P3\tGGCC"
    ):
        primer_input.input_data.primers = [
            {"name": "", "seq": "", "active": False}
        ]
        await paste_primers_click(primer_input, None)
        assert len(primer_input.input_data.primers) == 1
        assert primer_input.input_data.primers[0]["name"] == "P3"


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_primer_file_manager_operations() -> None:
    """Test loading and saving primers via PrimerFileManager."""
    view, _, input_data, _ = _setup_test_view()
    file_mgr = view.primer_input.file_manager

    # 1. _parse_primers_from_text
    tsv_content = (
        "# Comment\nATGC,P1\nGGCC\tP2\tThirdColumnDescription\n\nInvalid\n"
    )
    parsed = file_mgr._parse_primers_from_text(tsv_content)
    assert len(parsed) == 2
    assert parsed[1]["name"] == "P2 - ThirdColumnDescription"

    # 2. _serialise_primers_to_tsv
    serialised = file_mgr._serialise_primers_to_tsv(input_data.primers)
    assert "ATGC\tP1" in serialised

    # 3. load_primers_click
    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new_callable=AsyncMock,
        return_value=None,
    ):
        await file_mgr.load_primers_click(None)

    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new_callable=AsyncMock,
        return_value="ATGC,P1",
    ):
        input_data.primers = [{"name": "", "seq": "", "active": False}]
        await file_mgr.load_primers_click(None)
        assert input_data.primers[0]["name"] == "P1"

    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new_callable=AsyncMock,
        return_value="# no primers",
    ):
        await file_mgr.load_primers_click(None)

    # 4. save_primers_click
    input_data.primers = []
    await file_mgr.save_primers_click(None)

    input_data.primers = [{"name": "P1", "seq": "ATGC", "active": True}]
    with patch(
        "amplifyp.gui.utils.data_helpers.save_and_write_file",
        new_callable=AsyncMock,
    ):
        await file_mgr.save_primers_click(None)


def test_primer_action_controller_mutations() -> None:
    """Test selections, movements, reordering, and deletions in controller."""
    view, _, _, _ = _setup_test_view()
    act = view.primer_input.action_controller
    primer_input = view.primer_input

    # 1. Single click select / deselect
    mock_field = MagicMock(spec=ft.TextField)
    mock_field.focus = MagicMock(return_value=None)
    act.handle_row_click(0, mock_field)
    assert 0 in primer_input.selected_indices

    act.handle_row_click(0, mock_field)
    assert 0 not in primer_input.selected_indices

    # 2. Double click range selection
    act.handle_row_double_click(0, mock_field)
    assert act._click_a == 0
    act.handle_row_double_click(1, mock_field)
    assert act._click_a is None
    assert {0, 1}.issubset(primer_input.selected_indices)

    # 3. on_add_primer_row
    act.on_add_primer_row(0)
    assert len(view.input_data.primers) == 3

    # 4. move_primers
    primer_input.selected_indices = {1}
    act.move_primers({1}, -1)
    act.move_primers({0}, -1)  # top boundary
    act.move_primers({0}, 1)
    act.move_primers({2}, 1)  # bottom boundary

    # 5. reverse_complement_primers
    act.reverse_complement_primers({0})

    # 6. delete_primers
    act.delete_primers({0})

    # 7. Drag events
    act.handle_drag_start(0, MagicMock())

    ev_drag_up = MagicMock()
    ev_drag_up.local_delta = MagicMock(y=-40.0)
    act.handle_drag_update(0, ev_drag_up)

    ev_drag_down = MagicMock()
    ev_drag_down.local_delta = MagicMock(y=40.0)
    act.handle_drag_update(0, ev_drag_down)

    act.handle_drag_end(0, MagicMock())

    # 8. Header action clicks
    act.header_add_click(None)
    act.header_delete_click(None)
    act.header_up_click(None)
    act.header_down_click(None)


def test_primer_row_components_and_events() -> None:
    """Test PrimerRow selection changes, blur, errors, and Tm updates."""
    view, page, _, settings = _setup_test_view()
    row = view.primer_input.primers_list.controls[0]
    assert isinstance(row, PrimerRow)

    # 1. _on_selection_change
    ev_sel = MagicMock()
    ev_sel.control = row.name_field
    ev_sel.control.selection = MagicMock(base_offset=3)
    row._on_selection_change(ev_sel)
    assert row.name_field.data["cursor_pos"] == 3

    # 2. _on_blur
    ev_blur = MagicMock()
    ev_blur.control = row.name_field
    ev_blur.page = page
    row._on_blur(ev_blur)

    # 3. update_highlight_and_reorder
    row.update_highlight_and_reorder(is_focused=True, is_dup=True)
    row.update_highlight_and_reorder(is_focused=False, is_dup=True)
    row.update_highlight_and_reorder(is_focused=True, is_dup=False)
    row.update_highlight_and_reorder(is_focused=False, is_dup=False)

    # 4. set_error
    row.set_error({"name": "Name error", "seq": "Seq error"})
    row.set_error("Duplicate primer name")
    row.set_error("Other error")
    row.set_error(None)

    # 5. update_tm
    row.update_tm(settings)
    row.seq_field.value = "invalid sequence %%%"
    row.update_tm(settings)
    settings["show_primer_temperature"] = False
    row.update_tm(settings)


def test_template_casing_and_formatter() -> None:
    """Test TemplateInput casing controller, wrap calculation, and lines."""
    view, _, _, _ = _setup_test_view()
    template_input = view.template_input

    # 1. Casing clicks
    template_input.template_sequence.value = "atgcgatc"
    template_input.template_sequence.selection = ft.TextSelection(
        base_offset=0, extent_offset=4
    )
    upper_case_click(template_input, MagicMock())
    assert template_input.template_sequence.value.startswith("ATGC")

    lower_case_click(template_input, MagicMock())
    assert template_input.template_sequence.value.startswith("atgc")

    # Selection invalid / empty
    template_input.template_sequence.selection = None
    change_selection_case(template_input, to_upper=True)

    # 2. adjust_wrap_length & validate_bases_per_line
    assert validate_bases_per_line(template_input, "Auto") == "Auto"
    assert validate_bases_per_line(template_input, "60") == 60
    assert validate_bases_per_line(template_input, "invalid") is None

    template_input.bases_per_line_value_text.value = "Auto"
    wrap_len = adjust_wrap_length(template_input, 600.0, update=False)
    assert wrap_len >= 10

    template_input.bases_per_line_value_text.value = "40"
    wrap_len = adjust_wrap_length(template_input, 600.0, update=False)
    assert wrap_len == 40

    # 3. update_line_numbers
    update_line_numbers(template_input, update=False, gutter_only=True)


def test_input_events_coordination() -> None:
    """Test handle_field_focus, blur, submit, and empty row append."""
    view, _, _, _ = _setup_test_view()
    row = view.primer_input.primers_list.controls[0]
    assert isinstance(row, PrimerRow)

    # 1. handle_field_focus
    ev_focus_name = MagicMock()
    ev_focus_name.control = row.name_field
    handle_field_focus(view, ev_focus_name)

    ev_focus_seq = MagicMock()
    ev_focus_seq.control = row.seq_field
    handle_field_focus(view, ev_focus_seq)

    # 2. handle_field_blur
    ev_blur = MagicMock()
    ev_blur.control = row.seq_field
    handle_field_blur(view, ev_blur)

    # Blur on template sequence
    ev_blur_tmpl = MagicMock()
    ev_blur_tmpl.control = view.template_sequence
    handle_field_blur(view, ev_blur_tmpl)

    # 3. handle_field_submit
    ev_submit = MagicMock()
    ev_submit.control = row.name_field
    handle_field_submit(view, ev_submit)

    # 4. auto_add_empty_row_if_needed
    auto_add_empty_row_if_needed(view, row.seq_field)


def test_primer_info_panel_details() -> None:
    """Test PrimerInfoPanel for zero, single, and multiple selected primers."""
    view, _, _input_data, _ = _setup_test_view()
    primer_input = view.primer_input

    # 1. Zero selection / focus
    primer_input.selected_indices = set()
    primer_input.focused_primer_index = None
    primer_input._update_primer_info_panel()

    # 2. Single selection
    primer_input.selected_indices = {0}
    primer_input.focused_primer_index = 0
    primer_input._update_primer_info_panel()
    assert primer_input.info_header.visible is True

    # 3. Multiple selection (triggers dimer stats calculation)
    primer_input.selected_indices = {0, 1}
    primer_input._update_primer_info_panel()


def test_input_view_state_and_lifecycle() -> None:
    """Test InputView state queries, serialisation, and panel repositioning."""
    view, _, input_data, _ = _setup_test_view()

    # 1. get_primers & get_all_primers_state
    primers = view.get_primers()
    assert len(primers) == 2

    all_state = view.get_all_primers_state()
    assert len(all_state) == 2

    # 2. get_state & set_state
    state = view.get_state()
    assert "primers" in state
    view.set_state(state)

    # 3. reposition_primer_info_panel
    view.reposition_primer_info_panel()
    view._update_row_highlights()
    view._update_primer_info_panel()

    # 4. Clear and delete callbacks
    view._clear_template(None)
    assert input_data.template == ""

    view._clear_primers(None)
    assert len(input_data.primers) == 1

    view.primer_input.selected_indices = {0}
    view._delete_selected_primers(None)


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_template_input_file_operations() -> None:
    """Test TemplateInput file loading, saving, and formatting."""
    from amplifyp.gui.views.input.template.file_manager import (
        load_template_click,
        save_template_click,
    )

    view, _, input_data, _ = _setup_test_view()
    tmpl = view.template_input

    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new_callable=AsyncMock,
        return_value=">seq1\nATGCGATCGATC",
    ):
        await load_template_click(tmpl, MagicMock())
        assert "ATGCGATCGATC" in input_data.template

    # Empty cancel
    with patch(
        "amplifyp.gui.utils.data_helpers.pick_and_read_file",
        new_callable=AsyncMock,
        return_value=None,
    ):
        await load_template_click(tmpl, MagicMock())

    # Empty template save
    input_data.template = ""
    await save_template_click(tmpl, MagicMock())

    # Valid template save
    input_data.template = "ATGCGATCGATC"
    with patch(
        "amplifyp.gui.utils.data_helpers.save_and_write_file",
        new_callable=AsyncMock,
    ):
        await save_template_click(tmpl, MagicMock())

    tmpl._on_copy_click(MagicMock())
    tmpl._handle_menu_select("Auto")
    tmpl._handle_menu_select(60)
    tmpl._upper_case_click(MagicMock())
    tmpl._lower_case_click(MagicMock())
    tmpl._show_notification("Message")


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_all_subcomponent_edge_cases() -> None:
    """Test all remaining branch edge cases across input view submodules."""
    from amplifyp.gui.views.input.primer.validation import validate_primers

    view, page, input_data, settings = _setup_test_view()

    tasks: list[asyncio.Task[Any]] = []

    def mock_run_task(fn: Any, *args: Any) -> None:
        tasks.append(asyncio.create_task(fn(*args)))

    page.run_task = mock_run_task

    # 1. coordinator.py line 60 (non-PrimerRow in primers_list)
    view.primer_input.primers_list.controls.append(ft.Container())
    extracted = view.primer_input.coordinator.extract_primer_data_from_ui()
    assert len(extracted) >= 1
    view.primer_input.primers_list.controls.pop()

    # 2. validation.py line 153 (invalid primer auto-activation check)
    from amplifyp.gui.views.input.primer.validation import (
        reconcile_primer_states,
    )

    ui_pr = [{"name": "P1", "seq": "123INVALID", "active": True}]
    prev_pr = [{"name": "P1", "seq": "123INVALID", "active": False}]
    reconciled = reconcile_primer_states(ui_pr, prev_pr)
    validated = validate_primers(
        reconciled,
        ignore_inactive_name_dup=True,
        ignore_inactive_seq_dup=True,
    )
    assert len(validated) == 1

    # 3. primer/list.py (on_scroll, focused index clamp, non-row controls)
    list_ctrl = view.primer_input.primers_list
    list_ctrl.on_scroll(MagicMock(pixels=150.0))
    assert list_ctrl.scroll_pixels == 150.0

    view.primer_input.focused_primer_index = 999
    list_ctrl.update_list_ui()
    assert view.primer_input.focused_primer_index is None

    # 4. primer/row.py (Tm calculation in __init__ with invalid seq)
    row_bad_tm = PrimerRow(
        idx=0,
        name="Bad",
        seq="NNNNNNNN",
        is_active=True,
        is_dup=False,
        name_error=None,
        seq_error=None,
        font_family="Roboto Mono",
        name_column_width=150.0,
        settings=settings,
        on_change_handler=MagicMock(),
        handle_field_focus=MagicMock(),
        handle_field_blur=MagicMock(),
        handle_field_submit=MagicMock(),
        on_row_click=MagicMock(),
        on_row_double_click=MagicMock(),
        on_divider_pan=MagicMock(),
        on_divider_pan_end=MagicMock(),
        is_focused=True,
        is_last_row=True,
        on_drag_start=MagicMock(),
        on_drag_update=MagicMock(),
        on_drag_end=MagicMock(),
    )
    assert row_bad_tm is not None

    # Drag divider handlers in PrimerRow
    if row_bad_tm.divider.on_pan_update:
        row_bad_tm.divider.on_pan_update(MagicMock())
    if row_bad_tm.divider.on_pan_end:
        row_bad_tm.divider.on_pan_end(MagicMock())

    # 5. primer/action_controller.py
    # (delayed delete, double click out-of-bounds)
    act = view.primer_input.action_controller
    act._click_a = 999  # out of bounds anchor

    act.handle_row_double_click(0, MagicMock(spec=ft.TextField))
    assert act._click_a == 0
    act.handle_row_double_click(0, MagicMock(spec=ft.TextField))  # unselect
    assert act._click_a is None

    # Tm error in PrimerRow __init__
    with patch.object(
        settings, "calculate_primer_tm", side_effect=ValueError("Tm err")
    ):
        row_tm_err = PrimerRow(
            idx=0,
            name="TmErr",
            seq="ATGC",
            is_active=True,
            is_dup=False,
            name_error=None,
            seq_error=None,
            font_family="Roboto Mono",
            name_column_width=150.0,
            settings=settings,
            on_change_handler=MagicMock(),
            handle_field_focus=MagicMock(),
            handle_field_blur=MagicMock(),
            handle_field_submit=MagicMock(),
            on_row_click=MagicMock(),
            on_row_double_click=MagicMock(),
            on_divider_pan=MagicMock(),
            on_divider_pan_end=MagicMock(),
            is_focused=True,
            is_last_row=True,
            on_drag_start=MagicMock(),
            on_drag_update=MagicMock(),
            on_drag_end=MagicMock(),
        )
        assert row_tm_err.tm_text.value == "-"

    # Drag swapping
    input_data.primers = [
        {"name": "P1", "seq": "ATGC", "active": True},
        {"name": "P2", "seq": "GGCC", "active": True},
        {"name": "P3", "seq": "TTTT", "active": True},
    ]
    view.update_ui()
    view.primer_input.selected_indices = {0}
    act.handle_drag_start(0, MagicMock())
    ev_down = MagicMock(local_delta=MagicMock(y=100.0), delta_y=100.0)
    act.handle_drag_update(0, ev_down)
    ev_up = MagicMock(local_delta=MagicMock(y=-100.0), delta_y=-100.0)
    act.handle_drag_update(0, ev_up)
    act.handle_drag_end(0, MagicMock())

    # Delayed delete with mock page run_task
    act.delete_primers({0, 1})
    if tasks:
        await asyncio.gather(*tasks)

    # 6. input_view.py properties & delegates
    assert view.delete_selected_button is not None
    assert view.reverse_complement_button is not None
    assert view.primer_info_panel is not None
    assert view.info_header is not None
    assert view.info_seq_text is not None
    assert view.info_tm_text is not None
    assert view.info_pairs_text is not None
    assert view.info_redundancy_text is not None
    assert view.info_dimer_text is not None

    view.validation_errors = [{"name": "Err", "seq": None}]
    assert len(view.validation_errors) == 1

    view.focused_primer_index = 0
    assert view.focused_primer_index == 0

    # 7. events.py (pasted text handler, submit on template, focus/blur tasks)
    txt_fld = ft.TextField(value="ATGC")
    txt_fld.data = {"idx": 0, "field": "seq"}
    txt_fld.focus = AsyncMock()
    view._handle_pasted_text("P9\tATGCATGC\n", 0, "name", txt_fld)

    # Focus on seq with cursor reset task
    ev_focus_seq = MagicMock(control=txt_fld, page=page)
    handle_field_focus(view, ev_focus_seq)
    if tasks:
        await asyncio.gather(*tasks)

    # Reset currently focused control before blur test
    view._currently_focused_control = view.template_sequence

    # Blur debouncer timer trigger
    on_stop_called = False

    def on_stop(e: Any) -> None:
        nonlocal on_stop_called
        on_stop_called = True

    view.on_stop_editing_callback = on_stop
    view.template_sequence.data = {"idx": 0, "field": "name"}
    handle_field_blur(
        view, MagicMock(control=view.template_sequence, page=page)
    )
    await asyncio.sleep(0.2)

    assert on_stop_called or view._focus_debouncer is not None

    # Submit on template sequence
    tmpl_ev = MagicMock(control=view.template_sequence)
    view._handle_field_submit(tmpl_ev)

    # Pan and resize on input_view
    view.primer_input.primer_header.update = MagicMock()
    view.primer_input.primers_list.update = MagicMock()
    view._on_pan_update(MagicMock(local_delta=MagicMock(x=10.0)))
    view._on_pan_end(MagicMock())
    view._handle_resize(MagicMock())


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_input_additional_branches() -> None:
    """Test additional missing branches in primer input and events."""
    view, page, input_data, settings = _setup_test_view()
    primer_input = view.primer_input

    # 1. primer_input._reposition_info_panel top vs bottom
    settings["primer_info_panel_position"] = "top"
    primer_input._reposition_info_panel()
    settings["primer_info_panel_position"] = "bottom"
    primer_input._reposition_info_panel()

    # 2. _reverse_complement_selected_click
    primer_input.selected_indices = {0}
    primer_input._reverse_complement_selected_click(None)

    # 3. _on_toggle_all_primers tri-state
    primer_input.all_primers_checkbox.value = None
    primer_input._prev_header_checkbox_value = True
    primer_input._on_toggle_all_primers(MagicMock())

    primer_input.all_primers_checkbox.value = None
    primer_input._prev_header_checkbox_value = False
    primer_input._on_toggle_all_primers(MagicMock())

    # Empty primers toggle
    primer_input.input_data.primers = []
    primer_input._on_toggle_all_primers(MagicMock())

    # 4. _copy_primers_click and _paste_primers_click
    with patch(
        "amplifyp.gui.views.input.primer.clipboard.copy_primers_click",
        new_callable=AsyncMock,
    ):
        await primer_input._copy_primers_click(None)

    with patch(
        "amplifyp.gui.views.input.primer.clipboard.paste_primers_click",
        new_callable=AsyncMock,
    ):
        await primer_input._paste_primers_click(None)

    # 5. _load_primers_click and _save_primers_click
    with patch(
        "amplifyp.gui.views.input.primer.file_manager.PrimerFileManager.load_primers_click",
        new_callable=AsyncMock,
    ):
        await primer_input._load_primers_click(None)

    with patch(
        "amplifyp.gui.views.input.primer.file_manager.PrimerFileManager.save_primers_click",
        new_callable=AsyncMock,
    ):
        await primer_input._save_primers_click(None)

    # 6. Primer info panel fixed height when focused_idx is None
    info_panel = primer_input.primer_info_panel
    info_panel.update = MagicMock()
    with patch(
        "amplifyp.gui.views.input.primer.info_panel.PrimerInfoPanel.page",
        new=property(lambda self: page),
    ):
        settings["primer_info_panel_fixed_height"] = True
        info_panel.update_panel(
            focused_idx=None,
            input_data=input_data,
            app_page=page,
            on_update_highlights=MagicMock(),
        )
        assert info_panel.height == 260

        # 7. Primer info panel with empty seq or invalid index
        info_panel.update_panel(
            focused_idx=-1,
            input_data=input_data,
            app_page=page,
            on_update_highlights=MagicMock(),
        )
        input_data.primers = [{"name": "P1", "seq": "", "active": True}]
        info_panel.update_panel(
            focused_idx=0,
            input_data=input_data,
            app_page=page,
            on_update_highlights=MagicMock(),
        )

        # Primer info panel with valid primer having dimers
        input_data.primers = [
            {"name": "D1", "seq": "GCATGCATGCATGCATGC", "active": True},
            {"name": "D2", "seq": "GCATGCATGCATGCATGC", "active": True},
        ]
        info_panel.update_panel(
            focused_idx=0,
            input_data=input_data,
            app_page=page,
            on_update_highlights=MagicMock(),
        )

        # Primer info panel calculation exception handling
        # (outer try/except block)
        with patch(
            "amplifyp.dimer.PrimerDimerGenerator.generate_primer_dimer",
            side_effect=ValueError("Dimer fail"),
        ):
            info_panel.update_panel(
                focused_idx=0,
                input_data=input_data,
                app_page=page,
                on_update_highlights=MagicMock(),
            )

    # 8. Template input text selection clamping in _handle_change
    tmpl = view.template_input
    tmpl.template_sequence.selection = ft.TextSelection(
        base_offset=0, extent_offset=0
    )
    tmpl.template_sequence.value = "ATGC"
    tmpl._handle_change(MagicMock())

    tmpl.template_sequence.selection = ft.TextSelection(
        base_offset=100, extent_offset=100
    )
    tmpl._handle_change(MagicMock())


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_all_remaining_input_branches_to_100_percent() -> None:
    """Test all remaining branches across input views for 100% coverage."""
    view, page, input_data, settings = _setup_test_view()
    tasks: list[asyncio.Task[Any]] = []

    def mock_run_task(fn: Any, *args: Any) -> None:
        tasks.append(asyncio.create_task(fn(*args)))

    page.run_task = mock_run_task

    # 1. input_view.py properties
    assert view.template_circular is not None
    assert view.upper_case_button is not None
    assert view.lower_case_button is not None

    # Callbacks attached in _clear_primers, _delete_selected_primers
    change_called = False

    stop_called = False

    def on_chg(e: Any) -> None:
        nonlocal change_called
        change_called = True

    def on_stp(e: Any) -> None:
        nonlocal stop_called
        stop_called = True

    view.on_change = on_chg
    view.on_stop_editing_callback = on_stp

    view._clear_primers(None)
    assert change_called and stop_called

    change_called = False
    stop_called = False
    view.primer_input.selected_indices = {0}
    view._delete_selected_primers(None)
    assert change_called and stop_called

    change_called = False
    stop_called = False
    view._clear_template(None)
    assert change_called and stop_called

    # get_all_primers_state with empty primer
    input_data.primers = [{"name": "", "seq": "", "active": False}]
    assert len(view.get_all_primers_state()) == 0

    # 2. primer/clipboard.py
    # empty parsed
    from amplifyp.gui.views.input.primer.clipboard import (
        parse_primer_clipboard_text,
    )

    p_empty = parse_primer_clipboard_text("\n   \n")
    assert len(p_empty) == 0

    # valid_selected insertion with selected indices
    view.primer_input.selected_indices = {0}
    with patch(
        "flet.Clipboard.get",
        new_callable=AsyncMock,
        return_value="PNew\tATGCAAAA",
    ):
        await view.primer_input._paste_primers_click(None)

    # 3. primer/file_manager.py
    file_mgr = view.primer_input.file_manager
    # csv.reader error
    parsed_lines = file_mgr._parse_primers_from_text('""\n,\n')
    assert len(parsed_lines) == 0

    with (
        patch(
            "amplifyp.gui.utils.data_helpers.pick_and_read_file",
            new_callable=AsyncMock,
            return_value="some,data",
        ),
        patch.object(
            file_mgr,
            "_parse_primers_from_text",
            side_effect=ValueError("Parse error"),
        ),
    ):
        await file_mgr.load_primers_click(None)

    # 4. primer/layout_manager.py
    lm = view.primer_input.layout_manager
    view.primer_input._visible_rows_cache = []
    lm._cache_visible_rows_if_needed()  # hits early return

    # 5. primer/list.py
    # update with page present
    with patch(
        "amplifyp.gui.views.input.primer.list.PrimerList.page",
        new=property(lambda self: page),
    ):
        view.primer_input.primers_list.update_row_highlights()

    # 6. primer/row.py
    row = view.primer_input.primers_list.controls[0]
    if isinstance(row, PrimerRow):
        # selection without base_offset
        ev_no_base = MagicMock(control=MagicMock(selection=MagicMock(spec=[])))
        row._on_selection_change(ev_no_base)

        # _on_blur with coroutine and error handling
        async def mock_scroll(offset: int = 0) -> None:
            pass

        row.name_scroll.scroll_to = MagicMock(
            side_effect=lambda offset=0: mock_scroll(offset)
        )

        with (
            patch.object(PrimerRow, "page", new=property(lambda self: page)),
            patch.object(ft.ListView, "page", new=property(lambda self: page)),
        ):
            row._on_blur(MagicMock(control=row.name_field, page=page))
            if tasks:
                await asyncio.gather(*tasks)

    # 7. primer/validation.py
    from amplifyp.gui.views.input.primer.validation import (
        reconcile_primer_states,
    )

    val_res = reconcile_primer_states(
        [{"name": "P1", "seq": "ZZZ", "active": True}],
        prev_primers=[{"name": "P1", "seq": "ZZZ", "active": False}],
        auto_activate_new=True,
    )
    assert len(val_res) == 1

    # 8. template/formatter.py & template/input.py
    tmpl = view.template_input
    with (
        patch.object(ft.Container, "page", new=property(lambda self: page)),
        patch.object(ft.Row, "page", new=property(lambda self: page)),
        patch.object(ft.TextField, "page", new=property(lambda self: page)),
    ):
        tmpl.adjust_wrap_length(500.0, update=True)
        tmpl._update_line_numbers(update=True, gutter_only=True)
        tmpl._update_line_numbers(update=True, gutter_only=False)
        tmpl._update_status_bar(update=True)

    # template current_left_width
    tmpl.width = 400.0
    assert tmpl.current_left_width == 400.0
    tmpl.width = None
    assert tmpl.current_left_width > 0

    # template copy empty
    input_data.template = ""
    tmpl.template_sequence.value = ""
    tmpl.template_sequence.selection = None
    tmpl._on_copy_click(MagicMock())

    # template _validate_bases_per_line
    assert tmpl._validate_bases_per_line("50") == 50

    # 9. events.py & action_controller.py edge cases
    # auto_add_empty_row_if_needed
    input_data.primers = [
        {"name": "P1", "seq": "ATGC", "active": True},
    ]
    view.primer_input.validation_errors = [{"name": None, "seq": None}]
    view.update_ui()
    last_row = view.primer_input.primers_list.controls[-1]
    if isinstance(last_row, PrimerRow):
        last_row.seq_field.data = {"idx": 0, "field": "seq"}
        from amplifyp.gui.views.input.events import auto_add_empty_row_if_needed

        auto_add_empty_row_if_needed(view, last_row.seq_field)

    # on_change_handler with newline text submission
    from amplifyp.gui.views.input.events import on_change_handler

    mock_field_nl = ft.TextField(value="ATGC\n")
    mock_field_nl.data = {"idx": 0, "field": "seq"}
    on_change_handler(view, MagicMock(control=mock_field_nl))

    # on_change_handler with tab pasted text
    mock_field_tab = ft.TextField(value="P1\tATGC\nP2\tGGCC\n")
    mock_field_tab.data = {"idx": 0, "field": "seq"}
    on_change_handler(view, MagicMock(control=mock_field_tab))

    # handle_pasted_text empty
    view._handle_pasted_text("", 0, "seq", mock_field_tab)

    # handle_row_click with coroutine focus
    act = view.primer_input.action_controller

    async def mock_async_focus() -> None:
        pass

    mock_fld = MagicMock(spec=ft.TextField)
    mock_fld.focus = MagicMock(side_effect=mock_async_focus)
    act.handle_row_click(0, mock_fld)
    if tasks:
        await asyncio.gather(*tasks)

    # move_primers and reverse_complement_primers empty
    act.move_primers(set(), -1)
    act.reverse_complement_primers(set())
    act._delete_primers_impl(set())
    act._delete_primers_impl({9999})

    # drag start non-contiguous block scan
    view.primer_input.selected_indices = {0, 2}
    act.handle_drag_start(0, MagicMock())
    assert act.drag_block == [0]

    # header_add_click when focused_primer_index is set vs neither
    view.primer_input.selected_indices = set()
    view.primer_input.focused_primer_index = 0
    act.header_add_click(None)

    view.primer_input.focused_primer_index = None
    act.header_add_click(None)

    # primer_input content reset
    view.primer_input.content = None
    view.primer_input._reposition_info_panel()

    # primer_list_container update with page
    with patch.object(ft.Container, "page", new=property(lambda self: page)):
        view.primer_input.update_ui()
        view.primer_input._update_header_buttons_state()

    # 10. action_controller.py drag & selection block scanning
    input_data.primers = [
        {"name": "P0", "seq": "ATGC", "active": True},
        {"name": "P1", "seq": "GGCC", "active": True},
        {"name": "P2", "seq": "TTTT", "active": True},
        {"name": "P3", "seq": "CCCC", "active": True},
    ]
    view.update_ui()
    # Scanning contiguous block surrounding row 2: {1, 2, 3}
    view.primer_input.selected_indices = {1, 2, 3}
    act.handle_drag_start(2, MagicMock())
    assert act.drag_block == [1, 2, 3]

    # Drag update with no drag_block
    delattr(act, "drag_block")
    act.handle_drag_update(0, MagicMock())

    # Drag end when on_change_handler is None
    cast(Any, view.primer_input).on_change_handler = None
    act.handle_drag_end(0, MagicMock())

    # _delete_primers_impl when on_change_handler is None
    act._delete_primers_impl({id(input_data.primers[0])})

    # Restore on_change_handler
    cast(Any, view.primer_input).on_change_handler = MagicMock()

    # Double click on selected and focused row
    view.primer_input.selected_indices = {1}
    view.primer_input.focused_primer_index = 1
    act._click_a = None
    act.handle_row_double_click(1, MagicMock(spec=ft.TextField))
    assert 1 not in view.primer_input.selected_indices

    # 11. primer_input._on_toggle_all_primers
    # when cb_value is None and prev_value is None
    input_data.primers = [{"name": "P1", "seq": "ATGC", "active": True}]
    view.primer_input.all_primers_checkbox.value = None
    view.primer_input._prev_header_checkbox_value = None
    view.primer_input._on_toggle_all_primers(MagicMock())

    # _update_header_buttons_state when primer_header missing
    delattr(view.primer_input, "primer_header")
    view.primer_input._update_header_buttons_state()

    # 12. primer/list.py when non-PrimerRow is encountered during update_list_ui
    input_data.primers = [{"name": "P1", "seq": "ATGC", "active": True}]
    view.primer_input.primers_list.controls = [ft.Container()]
    view.primer_input.primers_list.update_list_ui()

    # 13. primer/row.py divider handlers, errors in __init__, and _on_blur
    row_err_init = PrimerRow(
        idx=0,
        name="ErrP",
        seq="ATGC",
        is_active=True,
        is_dup=False,
        name_error="Name err init",
        seq_error="Seq err init",
        font_family="Roboto Mono",
        name_column_width=150.0,
        settings=settings,
        on_change_handler=MagicMock(),
        handle_field_focus=MagicMock(),
        handle_field_blur=MagicMock(),
        handle_field_submit=MagicMock(),
        on_row_click=MagicMock(),
        on_row_double_click=MagicMock(),
        on_divider_pan=MagicMock(),
        on_divider_pan_end=MagicMock(),
        is_focused=True,
        is_last_row=True,
        on_drag_start=MagicMock(),
        on_drag_update=MagicMock(),
        on_drag_end=MagicMock(),
    )
    assert row_err_init.name_field.error == "Name err init"
    assert row_err_init.seq_field.error == "Seq err init"

    row = view.primer_input.primers_list.controls[0]
    if isinstance(row, PrimerRow):
        # Coroutine _on_blur with async exception inside _do_scroll
        async def failing_scroll(offset: int = 0) -> None:
            raise RuntimeError("Scroll failure")

        row.name_scroll.scroll_to = MagicMock(
            side_effect=lambda offset=0: failing_scroll(offset)
        )

        with (
            patch.object(PrimerRow, "page", new=property(lambda self: page)),
            patch.object(ft.ListView, "page", new=property(lambda self: page)),
            patch.object(
                ft.ListView,
                "update",
                side_effect=RuntimeError("Scroll update err"),
            ),
        ):
            row._on_blur(MagicMock(control=row.name_field, page=page))
            if tasks:
                await asyncio.gather(*tasks)

        # Synchronous _on_blur with update raising RuntimeError
        row.seq_scroll.scroll_to = MagicMock(return_value=None)
        with (
            patch.object(PrimerRow, "page", new=property(lambda self: page)),
            patch.object(ft.ListView, "page", new=property(lambda self: page)),
            patch.object(
                ft.ListView,
                "update",
                side_effect=RuntimeError("Scroll update err"),
            ),
        ):
            row._on_blur(MagicMock(control=row.seq_field, page=None))

    # 14. template/input.py and formatter.py exceptions during update
    tmpl = view.template_input
    parent_mock = MagicMock(right_fraction=0.4)
    with patch.object(
        type(tmpl), "parent", new=property(lambda self: parent_mock)
    ):
        assert tmpl.current_left_width > 0

    with patch.object(ft.Container, "page", new=property(lambda self: page)):
        tmpl.status_bar.update = MagicMock(
            side_effect=RuntimeError("Status bar err")
        )
        tmpl._update_status_bar(update=True)

        tmpl.update = MagicMock(side_effect=RuntimeError("Tmpl update err"))
        from amplifyp.gui.views.input.template.formatter import (
            update_line_numbers,
        )

        update_line_numbers(tmpl, update=True, gutter_only=False)

    # 15. events.py focus and blur timer callback direct execution
    def throw_no_page(self: Any) -> None:
        raise RuntimeError("No page")

    with patch.object(
        InputView,
        "page",
        new=property(throw_no_page),
    ):
        fld_seq = ft.TextField()
        fld_seq.data = {"idx": 0, "field": "seq"}
        ev_focus_err = MagicMock(control=fld_seq, page=None)
        handle_field_focus(view, ev_focus_err)

    # Blur timer_callback execution with page None and page present
    view._currently_focused_control = None
    with patch.object(
        view._focus_debouncer, "trigger", side_effect=lambda cb: cb()
    ):
        with patch.object(InputView, "page", new=property(lambda self: None)):
            handle_field_blur(
                view, MagicMock(control=view.template_sequence, page=None)
            )

        with patch.object(InputView, "page", new=property(lambda self: page)):
            view.on_stop_editing_callback = MagicMock()
            handle_field_blur(
                view, MagicMock(control=view.template_sequence, page=page)
            )
            view.on_stop_editing_callback.assert_called_once()

    # 16. primer/clipboard.py line 68
    from amplifyp.gui.views.input.primer.clipboard import (
        parse_primer_clipboard_text,
    )

    parsed_extra = parse_primer_clipboard_text("part1,part2,part3\n")
    assert len(parsed_extra) == 1

    class CustomLine(str):
        def split(self, *args: Any, **kwargs: Any) -> list[str]:
            return []

        def strip(self, *args: Any, **kwargs: Any) -> "CustomLine":
            return self

    mock_text = MagicMock()
    mock_text.splitlines.return_value = [CustomLine("foo")]
    parsed_empty_split = parse_primer_clipboard_text(mock_text)
    assert len(parsed_empty_split) == 0

    # 17. primer/file_manager.py line 77-78 (StopIteration / csv.Error)
    with patch(
        "amplifyp.gui.views.input.primer.file_manager.csv.reader",
        return_value=iter([]),
    ):
        parsed_err = file_mgr._parse_primers_from_text("P1\tATGC\n")
        assert len(parsed_err) == 0

    # 18. Delayed delete async execution
    with patch.object(
        type(view.primer_input), "page", new=property(lambda self: page)
    ):
        view.primer_input.action_controller.delete_primers({0})
        await asyncio.sleep(0.1)
        if tasks:
            await asyncio.gather(*tasks)

    # 19. primer/action_controller.py all remaining branches
    act = view.primer_input.action_controller
    view.primer_input.app_page = page

    # Coroutine focus task in handle_row_click (when not yet selected)
    view.primer_input.selected_indices = set()
    view.primer_input.focused_primer_index = None

    async def mock_coro_focus() -> None:
        pass

    mock_field_coro = MagicMock(spec=ft.TextField)
    mock_field_coro.focus = MagicMock(side_effect=mock_coro_focus)
    act.handle_row_click(0, mock_field_coro)
    if tasks:
        await asyncio.gather(*tasks)

    # Contiguous block boundary breaking on non-contiguous selected indices
    input_data.primers = [
        {"name": "P0", "seq": "ATGC", "active": True},
        {"name": "P1", "seq": "GGCC", "active": True},
        {"name": "P2", "seq": "TTTT", "active": True},
        {"name": "P3", "seq": "CCCC", "active": True},
        {"name": "P4", "seq": "AAAA", "active": True},
    ]
    view.update_ui()
    view.primer_input.selected_indices = {0, 2, 4}
    act.handle_drag_start(2, MagicMock())
    assert act.drag_block == [2]

    # Drag start when start_idx is not in selected_indices
    view.primer_input.selected_indices = {0}
    act.handle_drag_start(3, MagicMock())
    assert view.primer_input.selected_indices == {3}

    # Drag update boundary checks (at bottom and at top)
    view.primer_input.primers_list.controls = [
        ft.Container() for _ in range(10)
    ]
    act.drag_block = [len(input_data.primers) + 5, 0]
    act.current_drag_y = 100.0
    act.handle_drag_update(0, MagicMock(local_delta=MagicMock(y=10.0)))

    act.drag_block = [1, len(input_data.primers) + 5]
    act.current_drag_y = -100.0
    act.handle_drag_update(0, MagicMock(local_delta=MagicMock(y=-10.0)))

    # Header add click with no selection / focus
    view.primer_input.selected_indices = set()
    view.primer_input.focused_primer_index = None
    act.header_add_click(None)

    # 20. primer/row.py _on_blur exception paths
    row_blur_test = PrimerRow(
        idx=0,
        name="BlurTest",
        seq="ATGC",
        is_active=True,
        is_dup=False,
        name_error=None,
        seq_error=None,
        font_family="Roboto Mono",
        name_column_width=150.0,
        settings=settings,
        on_change_handler=MagicMock(),
        handle_field_focus=MagicMock(),
        handle_field_blur=MagicMock(),
        handle_field_submit=MagicMock(),
        on_row_click=MagicMock(),
        on_row_double_click=MagicMock(),
        on_divider_pan=MagicMock(),
        on_divider_pan_end=MagicMock(),
        is_focused=True,
        is_last_row=True,
        on_drag_start=MagicMock(),
        on_drag_update=MagicMock(),
        on_drag_end=MagicMock(),
    )

    async def failing_coro_scroll(offset: int = 0) -> None:
        raise RuntimeError("Scroll failure")

    captured_coro = None

    def record_task(fn: Any, *args: Any) -> None:
        nonlocal captured_coro
        captured_coro = fn(*args)

    mock_page = MagicMock()
    mock_page.run_task = record_task

    row_blur_test.name_scroll.scroll_to = MagicMock(
        side_effect=lambda offset=0: failing_coro_scroll(offset)
    )

    row_blur_test.name_scroll._Control__page = mock_page
    row_blur_test.name_scroll.update = MagicMock(
        side_effect=RuntimeError("Scroll update err")
    )

    row_blur_test._on_blur(
        MagicMock(control=row_blur_test.name_field, page=mock_page)
    )
    assert captured_coro is not None
    await captured_coro

    # Synchronous scroll where scroll_target.update raises RuntimeError
    row_blur_test.seq_scroll.scroll_to = MagicMock(return_value=None)
    row_blur_test.seq_scroll.update = MagicMock(
        side_effect=RuntimeError("Sync scroll update err")
    )
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        row_blur_test._on_blur(
            MagicMock(control=row_blur_test.seq_field, page=mock_page)
        )
