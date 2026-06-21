# Copyright (C) 2026 Fufu Fang
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

    def on_stop_editing(e: ft.Event | None) -> None:
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
    name_scroll = row.controls[2]
    divider = row.controls[3]
    seq_scroll = row.controls[4]

    assert isinstance(checkbox, ft.Checkbox)
    assert isinstance(active_divider, ft.Container)
    assert isinstance(name_scroll, ft.ListView)
    assert isinstance(divider, ft.GestureDetector)
    assert isinstance(seq_scroll, ft.ListView)

    name_field = container.name_field
    seq_field = container.seq_field
    assert isinstance(name_field, ft.TextField)
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
    """Test that rows with duplicate names or sequences are coloured red."""
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
    second_row = view.primers_list.controls[1]
    second_row.name_field.value = "P1"
    view.sync_to_state()

    # Both rows should have colour warning set to RED_100
    assert view.primers_list.controls[0].bgcolor == ft.Colors.RED_100
    assert view.primers_list.controls[1].bgcolor == ft.Colors.RED_100

    # Resolve duplicate name, introduce duplicate sequence (case-insensitive)
    second_row.name_field.value = "P2"
    second_row.seq_field.value = "gcatgcatgc"
    view.sync_to_state()

    assert view.primers_list.controls[0].bgcolor == ft.Colors.RED_100
    assert view.primers_list.controls[1].bgcolor == ft.Colors.RED_100


def test_input_view_activation_validation() -> None:
    """Test empty errors are shown only when the row is active/selected."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    view = InputView(mock_page, input_data)
    view.update_ui()

    assert len(view.primers_list.controls) == 1
    row = view.primers_list.controls[0]
    checkbox = row.checkbox
    name_field = row.name_field
    seq_field = row.seq_field

    # Initially, no errors are shown and checkbox is enabled
    assert name_field.error is None
    assert seq_field.error is None
    assert checkbox.disabled is False
    assert checkbox.value is False

    # Try to activate the primer by setting checkbox to checked
    checkbox.value = True
    view.sync_to_state(rebuild_if_needed=False)

    # Validation errors populated and checkbox value reverted to False
    assert name_field.error == "Name cannot be empty"
    assert seq_field.error == "Sequence cannot be empty"
    assert checkbox.disabled is False
    assert checkbox.value is False

    # Simulate typing or changing focus to deactivate the error state
    view.sync_to_state(rebuild_if_needed=False)

    # Empty validation errors should be cleared when syncing inactive state
    assert name_field.error is None
    assert seq_field.error is None

    # Fill name and sequence and check the box
    checkbox.value = True
    name_field.value = "NewPrimer"
    seq_field.value = "ATGATGATG"
    view.sync_to_state(rebuild_if_needed=False)

    # Errors should be cleared and checkbox active
    assert name_field.error is None
    assert seq_field.error is None
    assert checkbox.disabled is False
    assert checkbox.value is True


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

    Verifies behaviour when a primer is focused.
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
    seq_field = container.seq_field

    # Verify no error initially
    assert seq_field.error is None
    assert container.height == 30

    # Put an invalid character (X) in sequence
    seq_field.value = "GCATGCATGX"
    view.sync_to_state()

    # Re-fetch elements since view rebuilt
    container = view.primers_list.controls[0]
    seq_field = container.seq_field

    assert seq_field.error is not None
    assert "contains invalid characters" in seq_field.error
    assert container.height is None  # Row container expanded/autosized

    # Fix the sequence
    seq_field.value = "GCATGCATGC"
    view.sync_to_state()

    container = view.primers_list.controls[0]
    seq_field = container.seq_field

    assert seq_field.error is None
    assert container.height == 30


def test_input_view_row_highlighting() -> None:
    """Test focusing a primer row updates the background colour."""
    from amplifyp.gui.colours import GUIColours

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

    assert view.primers_list.controls[0].bgcolor == GUIColours.SELECTED_ROW_BG
    assert view.primers_list.controls[1].bgcolor is None

    # Focus on the second primer row (index 1)
    mock_event.control.data = 1
    view._handle_field_focus(mock_event)

    assert view.primers_list.controls[0].bgcolor is None
    assert view.primers_list.controls[1].bgcolor == GUIColours.SELECTED_ROW_BG


def test_input_view_row_single_click() -> None:
    """Test that a single click on the row container focuses it."""
    from amplifyp.gui.colours import GUIColours

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [{"name": "P1", "seq": "GCATGCATGC", "active": True}]

    view = InputView(mock_page, input_data)
    view.update_ui()

    container = view.primers_list.controls[0]
    name_field = container.name_field

    # Mock focus method on TextField
    name_field.focus = MagicMock()

    # Trigger row container click
    container.on_click(MagicMock())

    assert view.focused_primer_index == 0
    assert container.bgcolor == GUIColours.SELECTED_ROW_BG
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

    # Initially focused_primer_index is None. Reorder controls on header.
    header = view.primer_input.primer_header
    assert header.add_button.disabled is False
    assert header.delete_button.disabled is True
    assert header.up_button.disabled is True
    assert header.down_button.disabled is True

    # Select/focus Row 1 (P2)
    view.focused_primer_index = 1
    view._update_row_highlights()

    # Now header reorder controls should be enabled for Row 1
    assert header.add_button.disabled is False
    assert header.delete_button.disabled is False
    assert header.up_button.disabled is False  # Up enabled
    assert header.down_button.disabled is False  # Down enabled

    # Click Up on Row 1 (P2) via header button
    header.up_button.on_click(MagicMock())
    header = view.primer_input.primer_header

    # After moving P2 up, the list order should be P2, P1, P3, and focused
    # index should follow P2 (index 0)
    assert input_data.primers[0]["name"] == "P2"
    assert input_data.primers[1]["name"] == "P1"
    assert input_data.primers[2]["name"] == "P3"
    assert view.focused_primer_index == 0

    # In the updated UI, header buttons should update for index 0 (first row)
    assert header.add_button.disabled is False
    assert header.delete_button.disabled is False
    assert header.up_button.disabled is True  # Up is disabled on the first row
    assert header.down_button.disabled is False  # Down is enabled

    # Click Down on P2 (index 0) via header button
    header.down_button.on_click(MagicMock())
    header = view.primer_input.primer_header

    # P2 should be back at index 1, P1 at index 0, and focused index at 1
    assert input_data.primers[0]["name"] == "P1"
    assert input_data.primers[1]["name"] == "P2"
    assert input_data.primers[2]["name"] == "P3"
    assert view.focused_primer_index == 1

    # Now click Down on P2 (index 1) to swap with P3 (index 2)
    header.down_button.on_click(MagicMock())
    header = view.primer_input.primer_header

    # P2 should now be at index 2, P3 at index 1, and focused index at 2
    assert input_data.primers[0]["name"] == "P1"
    assert input_data.primers[1]["name"] == "P3"
    assert input_data.primers[2]["name"] == "P2"
    assert view.focused_primer_index == 2

    # Verify that the down button is disabled on the last row (index 2)
    assert (
        header.down_button.disabled is True
    )  # Down is disabled on the last row


def test_input_view_block_invalid_primer() -> None:
    """Test invalid primer does not block row editing.

    It should still allow checking the checkbox.
    """
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

    # Verify that the checkbox is NOT disabled
    # (since it is invalid but not empty)
    assert checkbox.disabled is False


def test_input_view_duplicate_validation_and_enabling() -> None:
    """Test duplicate name/sequence is invalid but checkbox remains enabled."""
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

    # Checkboxes must NOT be disabled, and active status should remain True
    assert input_data.primers[0]["active"] is True
    assert input_data.primers[1]["active"] is True
    assert view.primers_list.controls[0].checkbox.disabled is False
    assert view.primers_list.controls[1].checkbox.disabled is False


def test_app_views_disabled_on_invalid_selected() -> None:
    """Test PCR/Dimer buttons disabled.

    Also tests error banner shows when active primer is invalid.
    """
    import amplifyp.gui.app as gui_app

    mock_page = MagicMock(spec=ft.Page)
    mock_page.window = MagicMock()
    mock_page.controls = []
    mock_page.overlay = []
    mock_page.platform_brightness = "light"
    mock_page.web = False

    # We will invoke gui_app.main to test full integration of buttons
    # and error banner
    gui_app.main(mock_page)

    # Find components from main page controls
    # The view container is controls[1] or similar. Let's find InputView.
    input_view = None
    for control in mock_page.controls:
        if (
            isinstance(control, ft.Container)
            and hasattr(control, "content")
            and isinstance(control.content, InputView)
        ):
            input_view = control.content
            break
        # Let's also check nested content
        if hasattr(control, "controls"):
            for sub in control.controls:
                if isinstance(sub, ft.Container) and isinstance(
                    sub.content, InputView
                ):
                    input_view = sub.content
                    break

    # If not found directly, search overlay or search child controls
    if not input_view:
        for c in mock_page.controls:
            if hasattr(c, "controls"):
                for sub in c.controls:
                    # check if it is InputView or contains it
                    if isinstance(sub, ft.Container) and isinstance(
                        sub.content, ft.ResponsiveRow
                    ):
                        # custom header
                        pass
                    if (
                        isinstance(sub, ft.Container)
                        and hasattr(sub, "content")
                        and isinstance(sub.content, ft.Container)
                        and isinstance(sub.content.content, InputView)
                    ):
                        sub.content.content  # noqa: B018
                    elif isinstance(sub, ft.Container) and isinstance(
                        sub.content, InputView
                    ):
                        sub.content  # noqa: B018

    # Alternatively, let's test it by mocking or directly calling
    # the update_pcr_button_state logic.
    # In app.py, let's look at the closure variables or how it's tested.
    # Actually, we can test it directly on InputView and app logic.
    # Let's check update_pcr_button_state directly or mock the page and buttons.
    # Let's write a targeted test:
    input_data = GUIInput()
    input_data.template = "ATGATGC"
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "AAAAAAAAAA", "active": True},
    ]

    # Create refs for buttons
    pcr_btn = ft.FilledButton("PCR", disabled=False)
    visible_pcr_btn = ft.FilledButton("PCR", disabled=False)
    dimers_btn = ft.FilledButton("Primer Dimers", disabled=False)
    visible_dimers_btn = ft.FilledButton("Primer Dimers", disabled=False)

    # Define update_pcr_button_state with mocked references
    pcr_button_ref = MagicMock()
    pcr_button_ref.current = pcr_btn
    visible_pcr_button_ref = MagicMock()
    visible_pcr_button_ref.current = visible_pcr_btn
    dimers_button_ref = MagicMock()
    dimers_button_ref.current = dimers_btn
    visible_dimers_button_ref = MagicMock()
    visible_dimers_button_ref.current = visible_dimers_btn

    view = InputView(mock_page, input_data)
    view.update_ui()

    # Define the update function like in app.py
    def update_pcr_button_state() -> None:
        has_template = bool(input_data.template.strip())
        active_primers = input_data.get_active_primers()
        has_enough_primers = len(active_primers) >= 2

        has_invalid_selected = False
        for idx, p in enumerate(input_data.primers):
            if p.get("active", False):
                if idx < len(view.primer_input.validation_errors):
                    err = view.primer_input.validation_errors[idx]
                    if err.get("name") or err.get("seq"):
                        has_invalid_selected = True
                        break

        if hasattr(view.primer_input, "error_banner"):
            view.primer_input.error_banner.visible = has_invalid_selected

        pcr_is_enabled = (
            has_template and has_enough_primers and not has_invalid_selected
        )

        if pcr_button_ref.current:
            pcr_button_ref.current.disabled = not pcr_is_enabled
        if visible_pcr_button_ref.current:
            visible_pcr_button_ref.current.disabled = not pcr_is_enabled
        if dimers_button_ref.current:
            dimers_button_ref.current.disabled = (
                len(active_primers) < 1
            ) or has_invalid_selected
        if visible_dimers_button_ref.current:
            visible_dimers_button_ref.current.disabled = (
                len(active_primers) < 1
            ) or has_invalid_selected

    # Verify initially enabled (valid primers)
    update_pcr_button_state()
    assert pcr_btn.disabled is False
    assert dimers_btn.disabled is False
    assert view.primer_input.error_banner.visible is False

    # Make P1 invalid but keep it active
    view.primer_input.primers_list.controls[0].seq_field.value = "GCATGCATGX"
    view.sync_to_state()

    # Run status update
    update_pcr_button_state()

    # Buttons should now be disabled and error banner should be visible
    assert pcr_btn.disabled is True
    assert dimers_btn.disabled is True
    assert view.primer_input.error_banner.visible is True


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


def _run_async(coro: object) -> None:
    """Run an async coroutine in a fresh event loop in a separate thread."""
    import asyncio
    import concurrent.futures

    def _run() -> None:
        asyncio.run(coro)  # type: ignore[arg-type]

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        executor.submit(_run).result()


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
        "amplifyp.gui.views.input.template_input.ft.FilePicker",
        return_value=mock_file_picker_instance,
    ):
        view = InputView(mock_page, input_data)
        view.update_ui()

        # Check buttons exist
        assert hasattr(view.template_input, "save_template_button")
        assert hasattr(view.template_input, "load_template_button")

        # --- TEST SAVE ---
        # 1. Save empty template
        _run_async(
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
        _run_async(
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

        _run_async(
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

    def on_stop_editing(e: ft.Event | None) -> None:
        nonlocal stop_editing_called
        stop_editing_called = True

    view = InputView(mock_page, input_data, on_stop_editing=on_stop_editing)
    view.update_ui()

    p1_name_field = view.primers_list.controls[0].name_field
    p2_name_field = view.primers_list.controls[1].name_field

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


def test_delete_button_disabled_state() -> None:
    """Test delete button disabled state based on active primers."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    # 1. No active primers initially (empty fields lead to inactive)
    input_data.primers = [{"name": "", "seq": "", "active": False}]
    view = InputView(mock_page, input_data)
    view.update_ui()
    assert view.delete_selected_button.content == "Delete"
    assert view.delete_selected_button.visible is not False  # Always visible
    assert view.delete_selected_button.disabled is True

    # 2. Add an active primer, should become enabled (disabled = False)
    input_data.primers = [{"name": "P1", "seq": "GCATGCATGC", "active": True}]
    view.update_ui()
    assert view.delete_selected_button.disabled is False

    # 3. Setting active to False should disable it again
    input_data.primers[0]["active"] = False
    view.update_ui()
    assert view.delete_selected_button.disabled is True


def test_delete_selected_primers() -> None:
    """Test that deleting removes only active primers."""
    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    input_data.primers = [
        {"name": "P1", "seq": "GCATGCATGC", "active": True},
        {"name": "P2", "seq": "AAAAAAAAAA", "active": False},
        {"name": "P3", "seq": "GGGGGGGGGG", "active": True},
    ]
    view = InputView(mock_page, input_data)
    view.update_ui()

    assert len(input_data.primers) == 3
    assert view.delete_selected_button.disabled is False

    # Click delete
    view._delete_selected_primers(MagicMock(spec=ft.ControlEvent))

    # P1 and P3 (active) should be removed, P2 (inactive) should remain
    assert len(input_data.primers) == 1
    assert input_data.primers[0]["name"] == "P2"
    assert input_data.primers[0]["seq"] == "AAAAAAAAAA"
    assert view.delete_selected_button.disabled is True

    # If all remaining are deleted, fallback to a single empty inactive row
    input_data.primers[0]["active"] = True
    view.update_ui()
    assert view.delete_selected_button.disabled is False

    view._delete_selected_primers(MagicMock(spec=ft.ControlEvent))
    assert len(input_data.primers) == 1
    assert input_data.primers[0]["name"] == ""
    assert input_data.primers[0]["seq"] == ""
    assert input_data.primers[0]["active"] is False
    assert view.delete_selected_button.disabled is True


def test_tsv_paste_handling() -> None:
    """Test pasting TSV data and parsing it correctly."""
    from unittest.mock import MagicMock

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()
    view = InputView(mock_page, input_data)
    view.update_ui()

    # Get the text field (e.g. name field of first row)
    name_field = view.primers_list.controls[0].name_field

    # Trigger change handler with a multiline/tab value (paste)
    mock_event = MagicMock(spec=ft.ControlEvent)
    mock_event.control = name_field
    name_field.value = "PrimerA\tATGCTAG\nPrimerB\tCGATCGAT"

    view._on_change_handler(mock_event)

    # Check that both primers were added and not squashed
    assert len(input_data.primers) == 2
    assert input_data.primers[0]["name"] == "PrimerA"
    assert input_data.primers[0]["seq"] == "ATGCTAG"
    assert input_data.primers[1]["name"] == "PrimerB"
    assert input_data.primers[1]["seq"] == "CGATCGAT"


def test_enter_press_handling() -> None:
    """Test that pressing Enter triggers submit and strips the newline."""
    from unittest.mock import MagicMock

    mock_page = MagicMock(spec=ft.Page)
    input_data = GUIInput()

    submit_called = False

    def on_stop_editing(e: ft.Event | None) -> None:
        nonlocal submit_called
        submit_called = True

    view = InputView(mock_page, input_data, on_stop_editing=on_stop_editing)
    view.update_ui()

    # Get the name field
    name_field = view.primers_list.controls[0].name_field

    # Trigger change handler with a value ending with a newline (Enter press)
    mock_event = MagicMock(spec=ft.ControlEvent)
    mock_event.control = name_field
    name_field.value = "PrimerA\n"

    view._on_change_handler(mock_event)

    # Check that the newline was stripped from the field
    assert name_field.value == "PrimerA"
    # Check that the submission callback was invoked
    assert submit_called
    # Check that no new row was parsed/added (length remains 1)
    assert len(input_data.primers) == 1
