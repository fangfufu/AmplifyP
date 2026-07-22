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

"""Event coordination and handling for DNA template and primers input views."""

from typing import Any, cast

import flet as ft  # type: ignore[import-not-found, unused-ignore]

from amplifyp.gui.views.input.primer.row import PrimerRow


def handle_field_focus(input_view: Any, e: ft.Event[ft.TextField]) -> None:
    """Handle focus on input fields to cancel auto-trigger timer."""
    input_view._focus_debouncer.cancel()

    if e.control.data is not None:
        idx = (
            e.control.data["idx"]
            if isinstance(e.control.data, dict)
            else e.control.data
        )

        field = (
            e.control.data["field"]
            if isinstance(e.control.data, dict)
            else None
        )

        if idx not in input_view.primer_input.selected_indices:
            input_view.primer_input.selected_indices = {idx}
            input_view.primer_input._update_delete_button_disabled_state()
        input_view.primer_input.focused_primer_index = idx

        # Set touched status in state
        if 0 <= idx < len(input_view.input_data.primers):
            if field == "name":
                input_view.input_data.primers[idx]["name_touched"] = True
            elif field == "seq":
                input_view.input_data.primers[idx]["seq_touched"] = True
                try:
                    page = e.page or input_view.page
                except RuntimeError:
                    page = None
                if page and not getattr(
                    input_view, "_skip_seq_focus_reset", False
                ):

                    async def set_seq_cursor() -> None:
                        """Reset cursor to start of the text field.

                        This prevents losing focus when re-focusing.
                        """
                        await e.control.focus()
                        e.control.selection = ft.TextSelection(
                            base_offset=0, extent_offset=0
                        )
                        e.control.update()

                    page.run_task(set_seq_cursor)
                input_view._skip_seq_focus_reset = False

        input_view.primer_input._update_row_highlights()
        input_view.primer_input._update_primer_info_panel()
        if input_view.app_page:
            input_view.app_page.update()
    input_view._currently_focused_control = cast(ft.Control, e.control)


def handle_field_blur(input_view: Any, e: ft.Event[ft.TextField]) -> None:
    """Handle blur on input fields to trigger results page after a delay."""
    if (
        input_view._currently_focused_control is not None
        and input_view._currently_focused_control != e.control
    ):
        input_view.sync_to_state(rebuild_if_needed=False)
        return

    input_view._currently_focused_control = None

    input_view.sync_to_state(rebuild_if_needed=False)
    if e.control == input_view.template_sequence:
        input_view._adjust_template_wrap(update_first=True)

    if e.control.data is not None:
        idx = (
            e.control.data["idx"]
            if isinstance(e.control.data, dict)
            else e.control.data
        )
        if idx < len(input_view.primer_input.validation_errors):
            err = input_view.primer_input.validation_errors[idx]
            for row in input_view.primer_input.primers_list.controls:
                if isinstance(row, PrimerRow) and row.idx == idx:
                    row.set_error(err)
                    break
            input_view._auto_add_empty_row_if_needed(
                cast(ft.Control, e.control)
            )
            input_view.app_page.update()

    def timer_callback() -> None:
        """Execute the on-stop-editing callback after a short delay.

        This is triggered by a debouncer to handle cases where the user stops
        interacting with an input field.
        """
        try:
            page = input_view.page
        except RuntimeError:
            return
        if not page:
            return
        if input_view.on_stop_editing_callback:
            input_view.on_stop_editing_callback(None)

    input_view._focus_debouncer.trigger(timer_callback)


def handle_field_submit(input_view: Any, e: ft.Event[ft.TextField]) -> None:
    """Handle submission (Enter key) to immediately trigger results."""
    input_view._focus_debouncer.cancel()
    input_view.sync_to_state()
    if e.control == input_view.template_sequence:
        input_view._adjust_template_wrap(update_first=True)
    input_view._auto_add_empty_row_if_needed(cast(ft.Control, e.control))
    if input_view.app_page:
        input_view.app_page.update()
    if input_view.on_stop_editing_callback:
        input_view.on_stop_editing_callback(None)


def auto_add_empty_row_if_needed(input_view: Any, control: ft.Control) -> None:
    """Append a new empty row if the last row is filled and valid."""
    if (
        control.data is not None
        and isinstance(control.data, dict)
        and control.data.get("field") == "seq"
    ):
        idx = control.data["idx"]
        num_primers = len(input_view.input_data.primers)
        if idx == num_primers - 1:
            p = input_view.input_data.primers[idx]
            if p.get("name", "").strip() and p.get("seq", "").strip():
                if idx < len(input_view.primer_input.validation_errors):
                    err = input_view.primer_input.validation_errors[idx]
                    if not err.get("name") and not err.get("seq"):
                        input_view.input_data.primers.append(
                            {"name": "", "seq": "", "active": False}
                        )
                        input_view.primer_input.update_ui()


def on_change_handler(input_view: Any, e: ft.Event | None) -> None:
    """Handle change in input fields."""
    if (
        e
        and hasattr(e, "control")
        and isinstance(e.control, ft.TextField)
        and e.control.data is not None
    ):
        data = e.control.data

        if isinstance(data, dict) and "idx" in data and "field" in data:
            val = str(e.control.value or "")
            if "\t" in val or "\n" in val:
                non_empty_lines = [
                    line for line in val.splitlines() if line.strip()
                ]
                if "\t" not in val and len(non_empty_lines) <= 1:
                    e.control.value = val.replace("\n", "")
                    input_view._handle_field_submit(e)
                    return
                input_view._handle_pasted_text(
                    val, data["idx"], data["field"], e.control
                )
                return

    input_view.sync_to_state()
    input_view.primer_input._update_primer_info_panel()
    if input_view.on_change:
        input_view.on_change(e)


def handle_pasted_text(
    input_view: Any, text: str, idx: int, field: str, control: ft.TextField
) -> None:
    """Parse pasted text and insert into the primer list."""
    from amplifyp.gui.views.input.primer.clipboard import (
        parse_primer_clipboard_text,
    )

    parsed = parse_primer_clipboard_text(text)
    if not parsed:
        return

    primers = input_view.input_data.primers

    if (
        len(primers) == 1
        and not primers[0].get("name")
        and not primers[0].get("seq")
    ):
        primers.clear()
        idx = 0

    for i, new_p in enumerate(parsed):
        target_idx = idx + i
        if target_idx < len(primers):
            primers[target_idx]["name"] = new_p["name"]
            primers[target_idx]["seq"] = new_p["seq"]
            primers[target_idx]["active"] = False
        else:
            primers.append(new_p)

    input_view.update_ui()
    input_view.sync_to_state(rebuild_if_needed=True, skip_extract=True)
    if input_view.on_change:
        input_view.on_change(None)
