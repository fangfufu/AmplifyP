# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""Tests for InputView primer editing in-place."""

from unittest.mock import MagicMock

import flet as ft

from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.views.input import InputView


def test_input_view_row_boxes_editing() -> None:
    """Test that primers in InputView render in Row boxes and sync."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [{"name": "P1", "seq": "GCATGCATGC", "active": True}]

    stop_editing_called = False

    def on_stop_editing() -> None:
        nonlocal stop_editing_called
        stop_editing_called = True

    view = InputView(mock_page, input_data, on_stop_editing=on_stop_editing)

    # 1. Verify primer list has correct Row controls (1 row, no trailing row)
    assert len(view.primers_list.controls) == 1
    container = view.primers_list.controls[0]
    assert isinstance(container, ft.Container)

    row = container.content
    assert isinstance(row, ft.Row)

    checkbox_container = row.controls[0]
    assert isinstance(checkbox_container, ft.Container)
    checkbox = checkbox_container.content
    active_divider = row.controls[1]
    name_field = row.controls[2]
    divider = row.controls[3]
    seq_field = row.controls[4]

    assert isinstance(checkbox, ft.Checkbox)
    assert isinstance(active_divider, ft.Container)
    assert isinstance(name_field, ft.TextField)
    assert isinstance(divider, ft.GestureDetector)
    assert isinstance(seq_field, ft.TextField)

    assert checkbox.value is True
    assert name_field.value == "P1"
    assert seq_field.value == "GCATGCATGC"

    # 2. Modify values in the TextFields and verify sync/trigger
    name_field.value = "P1-Modified"
    seq_field.value = "AAAAAAAAAA"

    assert not stop_editing_called
    view._handle_field_submit(MagicMock())

    assert stop_editing_called
    assert input_data.primers[0]["name"] == "P1-Modified"
    assert input_data.primers[0]["seq"] == "AAAAAAAAAA"


def test_input_view_duplicate_warning() -> None:
    """Test that rows with duplicate names or sequences are colored red."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "AAAAAAAAAA", "active": True},
    ]

    view = InputView(mock_page, input_data)

    # Verify no duplicate highlights initially
    assert view.primers_list.controls[0].bgcolor is None
    assert view.primers_list.controls[1].bgcolor is None

    # Change the second primer's name to P1 (duplicate name)
    second_row = view.primers_list.controls[1].content
    second_row.controls[2].value = "P1"
    view.sync_to_state()

    # Both rows should have color warning set to RED_100
    assert view.primers_list.controls[0].bgcolor == ft.Colors.RED_100
    assert view.primers_list.controls[1].bgcolor == ft.Colors.RED_100

    # Resolve duplicate name, introduce duplicate sequence (case-insensitive)
    second_row.controls[2].value = "P2"
    second_row.controls[4].value = "gcatgcatgc"
    view.sync_to_state()

    assert view.primers_list.controls[0].bgcolor == ft.Colors.RED_100
    assert view.primers_list.controls[1].bgcolor == ft.Colors.RED_100


def test_input_view_focus_validation() -> None:
    """Test empty errors are shown only after both fields are focused."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    view = InputView(mock_page, input_data)
    view.update_ui()

    assert len(view.primers_list.controls) == 1
    row = view.primers_list.controls[0]
    checkbox = row.checkbox
    name_field = row.name_field
    seq_field = row.seq_field

    # Initially, no errors are shown even though they are empty
    assert name_field.error is None
    assert seq_field.error is None
    assert checkbox.disabled is True

    # Focus Name field only
    mock_event1 = MagicMock(spec=ft.ControlEvent)
    mock_event1.control = name_field
    view._handle_field_focus(mock_event1)
    view.sync_to_state(rebuild_if_needed=False)
    assert name_field.error is None
    assert seq_field.error is None

    # Focus Sequence field as well (now both touched)
    mock_event2 = MagicMock(spec=ft.ControlEvent)
    mock_event2.control = seq_field
    view._handle_field_focus(mock_event2)
    view.sync_to_state(rebuild_if_needed=False)

    # Now empty validation errors should be populated
    assert name_field.error == "Name cannot be empty"
    assert seq_field.error == "Sequence cannot be empty"
    assert checkbox.disabled is True

    # Fill name and sequence
    name_field.value = "NewPrimer"
    seq_field.value = "ATGATGATG"
    view.sync_to_state(rebuild_if_needed=False)

    # Errors should be cleared and checkbox enabled
    assert name_field.error is None
    assert seq_field.error is None
    assert checkbox.disabled is False


def test_input_view_clear_buttons() -> None:
    """Test that clear primers and template buttons work correctly."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.template = "ATGATGC"
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "AAAAAAAAAA", "active": True},
    ]

    view = InputView(mock_page, input_data)
    view.update_ui()

    assert view.template_sequence.value == "ATGATGC"
    assert len(input_data.primers) == 2  # 2 loaded
    assert view.clear_primers_button.content == "Clear"
    assert view.clear_template_button.content == "Clear"

    # 1. Trigger clear template button
    view._clear_template(MagicMock(spec=ft.ControlEvent))
    assert view.template_sequence.value == ""
    assert input_data.template == ""

    # 2. Trigger clear primers button
    view._clear_primers(MagicMock(spec=ft.ControlEvent))
    assert len(input_data.primers) == 1  # only default empty remaining
    assert input_data.primers[0]["name"] == ""
    assert input_data.primers[0]["seq"] == ""


def test_input_view_primer_info_panel() -> None:
    """Test that the primer info panel displays correct information.

    Verifies behavior when a primer is focused.
    """
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "ATCGGTACTTGTGACGCTAC", "active": True},
        {"name": "P2", "seq": "RACGGTACGTACGTACGTY", "active": True},
    ]

    view = InputView(mock_page, input_data)
    view.update_ui()

    # Initially not visible because focused_primer_index is None
    assert view.primer_info_panel.visible is False

    # Simulate focusing on first primer (index 0)
    mock_event = MagicMock(spec=ft.ControlEvent)
    mock_control = MagicMock()
    mock_control.data = 0
    mock_event.control = mock_control

    view._handle_field_focus(mock_event)

    assert view.focused_primer_index == 0
    assert view.primer_info_panel.visible is True
    assert view.info_header.content.value == "Primer: P1"
    assert view.info_seq_text.value == "20 bp:   ATCGGTACTTGTGACGCTAC"
    assert "Tm =" in view.info_tm_text.value
    assert "10 AT Pairs, 10 GC Pairs, 50.0% AT" in view.info_pairs_text.value
    assert view.info_redundancy_text.value == "No redundant bases."

    # Simulate focusing on second primer (index 1) which has redundant bases
    mock_control.data = 1
    view._handle_field_focus(mock_event)

    assert view.focused_primer_index == 1
    assert view.primer_info_panel.visible is True
    assert view.info_header.content.value == "Primer: P2"
    assert "2 redundant bases" in view.info_redundancy_text.value
    assert "redundancy fold = 4" in view.info_redundancy_text.value

    # Simulate clicking the close button on the info panel
    view.primer_info_panel.close_button.on_click(MagicMock())
    assert view.focused_primer_index is None
    assert view.primer_info_panel.visible is False


def test_input_view_sequence_validation() -> None:
    """Test that invalid sequence characters trigger validation errors."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [{"name": "P1", "seq": "GCATGCATGC", "active": True}]

    view = InputView(mock_page, input_data)
    view.update_ui()

    container = view.primers_list.controls[0]
    row = container.content
    seq_field = row.controls[4]

    # Verify no error initially
    assert seq_field.error is None
    assert container.height == 30

    # Put an invalid character (X) in sequence
    seq_field.value = "GCATGCATGX"
    view.sync_to_state()

    # Re-fetch elements since view rebuilt
    container = view.primers_list.controls[0]
    row = container.content
    seq_field = row.controls[4]

    assert seq_field.error is not None
    assert "contains invalid characters" in seq_field.error
    assert container.height is None  # Row container expanded/autosized

    # Fix the sequence
    seq_field.value = "GCATGCATGC"
    view.sync_to_state()

    container = view.primers_list.controls[0]
    row = container.content
    seq_field = row.controls[4]

    assert seq_field.error is None
    assert container.height == 30


def test_input_view_row_highlighting() -> None:
    """Test focusing a primer row updates the background color."""
    from amplifyp.gui.settings import GUIColors

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "AAAAAAAAAA", "active": True},
    ]

    view = InputView(mock_page, input_data)
    view.update_ui()

    # Initially, no highlights
    assert view.primers_list.controls[0].bgcolor is None
    assert view.primers_list.controls[1].bgcolor is None

    # Focus on the first primer row (index 0)
    mock_event = MagicMock()
    mock_event.control = MagicMock()
    mock_event.control.data = 0
    view._handle_field_focus(mock_event)

    assert view.primers_list.controls[0].bgcolor == GUIColors.SELECTED_ROW_BG
    assert view.primers_list.controls[1].bgcolor is None

    # Focus on the second primer row (index 1)
    mock_event.control.data = 1
    view._handle_field_focus(mock_event)

    assert view.primers_list.controls[0].bgcolor is None
    assert view.primers_list.controls[1].bgcolor == GUIColors.SELECTED_ROW_BG


def test_input_view_row_single_click() -> None:
    """Test that a single click on the row container focuses it."""
    from amplifyp.gui.settings import GUIColors

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [{"name": "P1", "seq": "GCATGCATGC", "active": True}]

    view = InputView(mock_page, input_data)
    view.update_ui()

    container = view.primers_list.controls[0]
    name_field = container.content.controls[2]

    # Mock focus method on TextField
    name_field.focus = MagicMock()

    # Trigger row container click
    container.on_click(MagicMock())

    assert view.focused_primer_index == 0
    assert container.bgcolor == GUIColors.SELECTED_ROW_BG
    name_field.focus.assert_called_once()


def test_input_view_primer_reordering() -> None:
    """Test primer reordering (up/down buttons and list swapping)."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "AAAAA", "active": True},
        {"name": "P2", "seq": "TTTTT", "active": True},
        {"name": "P3", "seq": "GGGGG", "active": True},
    ]

    view = InputView(mock_page, input_data)
    view.update_ui()

    # Verify initial controls length: 3 rows (no trailing empty row)
    assert len(view.primers_list.controls) == 3

    # Initially focused_primer_index is None. Reorder controls
    # should be invisible.
    for i in range(3):
        row = view.primers_list.controls[i].content
        reorder_controls = row.controls[5].content
        assert reorder_controls.visible is False

    # Select/focus Row 1 (P2)
    view.focused_primer_index = 1
    view._update_row_highlights()

    # Now Row 1 reorder controls should be visible.
    row1 = view.primers_list.controls[1].content
    row1_btns = row1.controls[5].content
    assert row1_btns.visible is True
    assert isinstance(row1_btns, ft.Row)
    assert row1_btns.controls[2].disabled is False  # Up (index 2)
    assert row1_btns.controls[3].disabled is False  # Down (index 3)

    # Row 0 and Row 2 reorder controls should be invisible
    assert (
        view.primers_list.controls[0].content.controls[5].content.visible
        is False
    )
    assert (
        view.primers_list.controls[2].content.controls[5].content.visible
        is False
    )

    # Click Up on Row 1 (P2)
    row1_btns.controls[2].on_click(MagicMock())

    # After moving P2 up, the list order should be P2, P1, P3, and focused
    # index should follow P2 (index 0)
    assert input_data.primers[0]["name"] == "P2"
    assert input_data.primers[1]["name"] == "P1"
    assert input_data.primers[2]["name"] == "P3"
    assert view.focused_primer_index == 0

    # In the updated UI, Row 0 (P2) should now have the buttons since it
    # has focus
    row0_btns = view.primers_list.controls[0].content.controls[5].content
    assert row0_btns.visible is True
    assert row0_btns.controls[2].disabled is True  # Up (first row)
    assert row0_btns.controls[3].disabled is False  # Down

    # Click Down on Row 0 (P2)
    row0_btns.controls[3].on_click(MagicMock())

    # P2 should be back at index 1, P1 at index 0, and focused index at 1
    assert input_data.primers[0]["name"] == "P1"
    assert input_data.primers[1]["name"] == "P2"
    assert input_data.primers[2]["name"] == "P3"
    assert view.focused_primer_index == 1

    # Now click Down on P2 (index 1) to swap with P3 (index 2)
    row1_btns = view.primers_list.controls[1].content.controls[5].content
    row1_btns.controls[3].on_click(MagicMock())

    # P2 should now be at index 2, P3 at index 1, and focused index at 2
    assert input_data.primers[0]["name"] == "P1"
    assert input_data.primers[1]["name"] == "P3"
    assert input_data.primers[2]["name"] == "P2"
    assert view.focused_primer_index == 2

    # Verify that the down button is disabled on the last row (index 2)
    row2_btns = view.primers_list.controls[2].content.controls[5].content
    assert (
        row2_btns.controls[3].disabled is True
    )  # Down is disabled on the last row


def test_input_view_block_invalid_primer() -> None:
    """Test invalid primer doesn't block row editing but disables checkbox."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    view = InputView(mock_page, input_data)
    view.update_ui()

    assert len(view.primers_list.controls) == 1
    row = view.primers_list.controls[0]
    checkbox = row.checkbox
    name_field = row.name_field
    seq_field = row.seq_field

    # Fill name and sequence with an invalid character (X)
    name_field.value = "P1"
    seq_field.value = "GCATGCATGX"

    view.sync_to_state(rebuild_if_needed=False)

    # Verify that the checkbox is disabled
    assert checkbox.disabled is True

    # Correct the sequence to be valid (all valid bases)
    seq_field.value = "GCATGCATGC"
    view.sync_to_state(rebuild_if_needed=False)

    # Now verify that the checkbox is enabled but NOT auto-activated (unchecked)
    assert checkbox.disabled is False
    assert checkbox.value is False
    assert len(view.primers_list.controls) == 1


def test_input_view_duplicate_validation_and_enabling() -> None:
    """Test duplicate name/sequence is invalid & cannot be enabled."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "AAAAAAAAAA", "active": True},
    ]

    view = InputView(mock_page, input_data)
    view.update_ui()

    # Verify initially valid and active
    assert view.primer_input.validation_errors[0] == {"name": None, "seq": None}
    assert view.primer_input.validation_errors[1] == {"name": None, "seq": None}
    assert view.primers_list.controls[0].checkbox.disabled is False
    assert view.primers_list.controls[1].checkbox.disabled is False

    # 1. Test Duplicate Name: Make second primer's name "P1"
    view.primers_list.controls[1].name_field.value = "P1"
    view.sync_to_state()

    # Check both are marked as invalid with "Duplicate primer name"
    assert view.primer_input.validation_errors[0] == {
        "name": "Duplicate primer name",
        "seq": None,
    }
    assert view.primer_input.validation_errors[1] == {
        "name": "Duplicate primer name",
        "seq": None,
    }

    # Both must be inactive (active = False) and checkboxes disabled
    assert input_data.primers[0]["active"] is False
    assert input_data.primers[1]["active"] is False
    assert view.primers_list.controls[0].checkbox.disabled is True
    assert view.primers_list.controls[1].checkbox.disabled is True

    # 2. Resolve duplicate name but introduce duplicate sequence
    view.primers_list.controls[1].name_field.value = "P2"
    view.primers_list.controls[1].seq_field.value = "gcatgcatgc"
    view.sync_to_state()

    # Check both are marked invalid with "Duplicate primer sequence"
    assert view.primer_input.validation_errors[0] == {
        "name": None,
        "seq": "Duplicate primer sequence",
    }
    assert view.primer_input.validation_errors[1] == {
        "name": None,
        "seq": "Duplicate primer sequence",
    }
    assert input_data.primers[0]["active"] is False
    assert input_data.primers[1]["active"] is False
    assert view.primers_list.controls[0].checkbox.disabled is True
    assert view.primers_list.controls[1].checkbox.disabled is True


def test_header_checkbox_indeterminate_empty() -> None:
    """Test header checkbox is indeterminate when primer list is empty."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    view = InputView(mock_page, input_data)

    assert view.primer_input.all_primers_checkbox.value is None


def test_header_checkbox_all_active() -> None:
    """Test header checkbox is checked when all primers are active."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "AAAAA", "active": True},
        {"name": "P2", "seq": "TTTTT", "active": True},
    ]
    view = InputView(mock_page, input_data)
    view.sync_to_state()

    assert view.primer_input.all_primers_checkbox.value is True


def test_header_checkbox_all_inactive() -> None:
    """Test header checkbox is unchecked when all primers are inactive."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "AAAAA", "active": False},
        {"name": "P2", "seq": "TTTTT", "active": False},
    ]
    view = InputView(mock_page, input_data)
    view.sync_to_state()

    assert view.primer_input.all_primers_checkbox.value is False


def test_header_checkbox_mixed() -> None:
    """Test header checkbox is indeterminate when some primers are active."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "AAAAA", "active": True},
        {"name": "P2", "seq": "TTTTT", "active": False},
    ]
    view = InputView(mock_page, input_data)
    view.sync_to_state()

    assert view.primer_input.all_primers_checkbox.value is None


def test_header_checkbox_click_all_active() -> None:
    """Test checkbox value is authoritative: True→activate, False→deactivate."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "AAAAA", "active": True},
        {"name": "P2", "seq": "TTTTT", "active": True},
    ]
    view = InputView(mock_page, input_data)
    view.sync_to_state()

    assert view.primer_input.all_primers_checkbox.value is True

    # Simulate Flet cycle: True → None (first click from all-active state)
    view.primer_input.all_primers_checkbox.value = None
    view.primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))

    non_empty = [
        p
        for p in input_data.primers
        if str(p.get("name", "")).strip() or p.get("seq", "").strip()
    ]
    assert all(p["active"] for p in non_empty)

    # Simulate Flet cycle: None → False (second click)
    view.primer_input.all_primers_checkbox.value = False
    view.primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))

    non_empty = [
        p
        for p in input_data.primers
        if str(p.get("name", "")).strip() or p.get("seq", "").strip()
    ]
    assert all(not p["active"] for p in non_empty)


def test_header_checkbox_click_partial() -> None:
    """Test checkbox value is authoritative: True→activate, False→deactivate."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "AAAAA", "active": True},
        {"name": "P2", "seq": "TTTTT", "active": False},
    ]
    view = InputView(mock_page, input_data)
    view.sync_to_state()

    assert view.primer_input.all_primers_checkbox.value is None

    # Simulate Flet cycle: None → True (first click from mixed state)
    view.primer_input.all_primers_checkbox.value = True
    view.primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))

    non_empty = [
        p
        for p in input_data.primers
        if str(p.get("name", "")).strip() or p.get("seq", "").strip()
    ]
    assert all(p["active"] for p in non_empty)

    # Simulate Flet cycle: True → None (second click)
    view.primer_input.all_primers_checkbox.value = None
    view.primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))

    non_empty = [
        p
        for p in input_data.primers
        if str(p.get("name", "")).strip() or p.get("seq", "").strip()
    ]
    assert all(not p["active"] for p in non_empty)


def test_header_checkbox_click_all_inactive() -> None:
    """Test checkbox value is authoritative: True→activate, False→deactivate."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "AAAAA", "active": False},
        {"name": "P2", "seq": "TTTTT", "active": False},
    ]
    view = InputView(mock_page, input_data)
    view.sync_to_state()

    assert view.primer_input.all_primers_checkbox.value is False

    # Simulate Flet cycle: False → None (first click from all-inactive state)
    view.primer_input.all_primers_checkbox.value = None
    view.primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))

    non_empty = [
        p
        for p in input_data.primers
        if str(p.get("name", "")).strip() or p.get("seq", "").strip()
    ]
    assert all(p["active"] for p in non_empty)

    # Simulate Flet cycle: None → True (second click)
    view.primer_input.all_primers_checkbox.value = True
    view.primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))

    non_empty = [
        p
        for p in input_data.primers
        if str(p.get("name", "")).strip() or p.get("seq", "").strip()
    ]
    assert all(p["active"] for p in non_empty)


def test_header_checkbox_tristate_cycle() -> None:
    """Test the Flet tri-state cycle: None→True→None→False→None."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "AAAAA", "active": True},
        {"name": "P2", "seq": "TTTTT", "active": False},
    ]
    view = InputView(mock_page, input_data)
    view.sync_to_state()

    primer_input = view.primer_input
    cb = primer_input.all_primers_checkbox

    # Simulating Flet's tri-state cycle on click:
    # None→True→None→False→None→True...
    # The on_change handler fires AFTER Flet updates the checkbox value.

    # Click 1: Flet cycles None → True, handler sees cb_value=True
    cb.value = True
    primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))
    assert all(
        p["active"] for p in input_data.primers if p.get("seq", "").strip()
    )
    assert primer_input._prev_header_checkbox_value is True

    # Click 2: Flet cycles True → None, handler sees None, prev=True
    # → deactivate
    cb.value = None
    primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))
    assert all(
        not p["active"] for p in input_data.primers if p.get("seq", "").strip()
    )
    assert primer_input._prev_header_checkbox_value is None

    # Click 3: Flet cycles None → False, handler sees False → deactivate
    cb.value = False
    primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))
    assert all(
        not p["active"] for p in input_data.primers if p.get("seq", "").strip()
    )
    assert primer_input._prev_header_checkbox_value is False

    # Click 4: Flet cycles False → None, handler sees None, prev=False
    # → activate
    cb.value = None
    primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))
    assert all(
        p["active"] for p in input_data.primers if p.get("seq", "").strip()
    )
    assert primer_input._prev_header_checkbox_value is None

    # Click 5: Flet cycles None → True, handler sees cb_value=True → activate
    cb.value = True
    primer_input._on_toggle_all_primers(MagicMock(spec=ft.ControlEvent))
    assert all(
        p["active"] for p in input_data.primers if p.get("seq", "").strip()
    )
    assert primer_input._prev_header_checkbox_value is True


def test_template_load_save() -> None:
    """Test loading and saving template using FilePicker."""
    from unittest.mock import AsyncMock, MagicMock, patch

    from amplifyp.gui.user_data import GUIInput

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = True
    input_data = GUIInput()

    # Mock FilePicker class
    mock_file_picker_instance = MagicMock(spec=ft.FilePicker)
    mock_file_picker_instance.save_file = AsyncMock(return_value="template.txt")
    mock_page.file_picker = mock_file_picker_instance

    with patch(
        "amplifyp.gui.views.input.template_file_manager.ft.FilePicker",
        return_value=mock_file_picker_instance,
    ):
        view = InputView(mock_page, input_data)
        view.update_ui()

        # Check buttons exist
        assert hasattr(view.template_input, "save_template_button")
        assert hasattr(view.template_input, "load_template_button")

        # --- TEST SAVE ---
        # 1. Save empty template
        import asyncio

        asyncio.run(
            view.template_input.save_template_button.on_click(
                MagicMock(spec=ft.ControlEvent)
            )
        )
        # Verify snackbar was shown (and save_file was not called because
        # template is empty)
        mock_file_picker_instance.save_file.assert_not_called()

        # 2. Save non-empty template
        input_data.template = "ATGCATGC"
        view.update_ui()
        asyncio.run(
            view.template_input.save_template_button.on_click(
                MagicMock(spec=ft.ControlEvent)
            )
        )
        mock_file_picker_instance.save_file.assert_called_once()
        kwargs = mock_file_picker_instance.save_file.call_args[1]
        assert kwargs["file_name"] == "template.txt"
        assert kwargs["allowed_extensions"] == ["txt"]
        assert kwargs["src_bytes"] == b"ATGCATGC"

        # --- TEST LOAD ---
        mock_file = MagicMock(spec=ft.FilePickerFile)
        mock_file.name = "template.txt"
        mock_file.bytes = b"GGGCCC"
        mock_file.path = None

        mock_file_picker_instance.pick_files = AsyncMock(
            return_value=[mock_file]
        )

        asyncio.run(
            view.template_input.load_template_button.on_click(
                MagicMock(spec=ft.ControlEvent)
            )
        )

        mock_file_picker_instance.pick_files.assert_called_once_with(
            dialog_title="Load",
            allowed_extensions=["txt"],
            file_type=ft.FilePickerFileType.CUSTOM,
            with_data=True,
        )

        assert input_data.template == "GGGCCC"
        assert view.template_sequence.value == "GGGCCC"


def test_input_view_focus_preservation_transition() -> None:
    """Test that focus transitions cancel timer and preserve state."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "AAAAAAAAAA", "active": True},
    ]

    stop_editing_called = False

    def on_stop_editing() -> None:
        nonlocal stop_editing_called
        stop_editing_called = True

    view = InputView(mock_page, input_data, on_stop_editing=on_stop_editing)
    view.update_ui()

    p1_name_field = view.primers_list.controls[0].content.controls[2]
    p2_name_field = view.primers_list.controls[1].content.controls[2]

    # Focus P1
    mock_focus_event_p1 = MagicMock(spec=ft.ControlEvent)
    mock_focus_event_p1.control = p1_name_field
    p1_name_field.data = 0
    view._handle_field_focus(mock_focus_event_p1)
    assert view._currently_focused_control == p1_name_field

    # Focus transitions from P1 to P2
    # Scenario 1: Focus event on P2 occurs, then Blur event on P1 occurs
    mock_focus_event_p2 = MagicMock(spec=ft.ControlEvent)
    mock_focus_event_p2.control = p2_name_field
    p2_name_field.data = 1
    view._handle_field_focus(mock_focus_event_p2)
    assert view._currently_focused_control == p2_name_field

    # Blur event on P1
    mock_blur_event_p1 = MagicMock(spec=ft.ControlEvent)
    mock_blur_event_p1.control = p1_name_field
    p1_name_field.data = 0
    view._handle_field_blur(mock_blur_event_p1)

    # Verify that the debouncer was NOT scheduled/triggered (timer is None)
    assert not stop_editing_called
    assert view._focus_debouncer._timer is None
    assert view._currently_focused_control == p2_name_field

    # Scenario 2: Blur event on P1 occurs, then Focus event on P2 occurs
    view._currently_focused_control = p1_name_field
    # Blur event on P1
    view._handle_field_blur(mock_blur_event_p1)
    # The debouncer should be scheduled (timer is not None)
    assert view._currently_focused_control is None
    assert view._focus_debouncer._timer is not None

    # Focus event on P2 immediately cancels the debouncer (timer becomes None)
    view._handle_field_focus(mock_focus_event_p2)
    assert view._currently_focused_control == p2_name_field
    assert view._focus_debouncer._timer is None
