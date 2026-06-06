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
from amplifyp.gui.views.input_view import InputView


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

    # 1. Verify primer list has correct Row controls (1 + 1 empty)
    assert len(view.primers_list.controls) == 2
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
    view.handle_field_submit(MagicMock())

    assert stop_editing_called
    assert input_data.primers[0]["name"] == "P1-Modified"
    assert input_data.primers[0]["seq"] == "AAAAAAAAAA"


def test_input_view_auto_inactivation_rules() -> None:
    """Test primer becomes inactive if name or sequence is cleared."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [{"name": "P1", "seq": "GCATGCATGC", "active": True}]

    view = InputView(mock_page, input_data)

    # Clear only the name
    container = view.primers_list.controls[0]
    row = container.content
    checkbox_container = row.controls[0]
    checkbox = (
        checkbox_container.content
        if isinstance(checkbox_container, ft.Container)
        else checkbox_container
    )
    name_field = row.controls[2]

    name_field.value = ""
    view.sync_to_state()

    # Verify it is now inactive in the state and checkbox updated
    assert input_data.primers[0]["active"] is False
    assert checkbox.value is False

    # Restore name, clear sequence
    input_data.primers = [{"name": "P1", "seq": "GCATGCATGC", "active": True}]
    view.update_ui()

    container = view.primers_list.controls[0]
    row = container.content
    checkbox_container2 = row.controls[0]
    checkbox2 = (
        checkbox_container2.content
        if isinstance(checkbox_container2, ft.Container)
        else checkbox_container2
    )
    seq_field = row.controls[4]

    seq_field.value = ""
    view.sync_to_state()

    assert input_data.primers[0]["active"] is False
    assert checkbox2.value is False


def test_input_view_deletion_rules() -> None:
    """Test that a primer is deleted if both name and sequence are cleared."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [{"name": "P1", "seq": "GCATGCATGC", "active": True}]

    view = InputView(mock_page, input_data)

    container = view.primers_list.controls[0]
    row = container.content
    name_field = row.controls[2]
    seq_field = row.controls[4]

    name_field.value = ""
    seq_field.value = ""

    view.sync_to_state()

    # Verify that the primer is deleted (excluding the trailing empty row)
    assert len(input_data.primers) == 1  # 1 trailing empty row
    assert input_data.primers[0]["name"] == ""
    assert input_data.primers[0]["seq"] == ""
    assert len(view.primers_list.controls) == 1


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


def test_input_view_trailing_row_activation() -> None:
    """Test that the trailing empty row is inactive by default.

    It should become active once both fields are filled.
    """
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    # Initially has no primers, so we should get exactly 1 trailing empty row
    view = InputView(mock_page, input_data)

    assert len(view.primers_list.controls) == 1
    container = view.primers_list.controls[0]
    row = container.content
    checkbox_container = row.controls[0]
    checkbox = checkbox_container.content
    name_field = row.controls[2]
    seq_field = row.controls[4]

    # 1. By default, trailing empty row is inactive (unchecked)
    assert checkbox.value is False

    # 2. Fill only Name, should remain inactive (unchecked)
    name_field.value = "NewPrimer"
    view.sync_to_state()

    # Get updated controls since UI rebuilt
    assert len(view.primers_list.controls) == 2
    row = view.primers_list.controls[0].content
    checkbox = row.controls[0].content
    name_field = row.controls[2]
    seq_field = row.controls[4]
    assert checkbox.value is False

    # 3. Fill Sequence too, should become active (checked)
    seq_field.value = "ATGATGATG"
    view.sync_to_state()

    # The active status transitions from False to True, and a new
    # trailing empty row will be added
    assert (
        len(input_data.primers) == 2
    )  # The filled one + a new trailing empty one
    assert input_data.primers[0]["name"] == "NewPrimer"
    assert input_data.primers[0]["seq"] == "ATGATGATG"
    assert input_data.primers[0]["active"] is True

    # The new trailing row should be inactive
    assert input_data.primers[1]["name"] == ""
    assert input_data.primers[1]["seq"] == ""
    assert input_data.primers[1]["active"] is False


def test_input_view_clear_buttons() -> None:
    """Test that clear primers button works correctly."""
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
    assert len(input_data.primers) == 3  # 2 loaded + 1 trailing empty

    # 1. Trigger clear primers button
    view.clear_primers(MagicMock(spec=ft.ControlEvent))
    assert len(input_data.primers) == 1  # only trailing empty remaining
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

    view.handle_field_focus(mock_event)

    assert view.focused_primer_index == 0
    assert view.primer_info_panel.visible is True
    assert view.info_header.content.value == "Primer: P1"
    assert view.info_seq_text.value == "20 bp:   ATCGGTACTTGTGACGCTAC"
    assert "Tm =" in view.info_tm_text.value
    assert "10 AT Pairs, 10 GC Pairs, 50.0% AT" in view.info_pairs_text.value
    assert view.info_redundancy_text.value == "No redundant bases."

    # Simulate focusing on second primer (index 1) which has redundant bases
    mock_control.data = 1
    view.handle_field_focus(mock_event)

    assert view.focused_primer_index == 1
    assert view.primer_info_panel.visible is True
    assert view.info_header.content.value == "Primer: P2"
    assert "2 redundant bases" in view.info_redundancy_text.value
    assert "redundancy fold = 4" in view.info_redundancy_text.value
