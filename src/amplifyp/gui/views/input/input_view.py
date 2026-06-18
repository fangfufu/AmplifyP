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

"""Input view composing DNA Template input and Primer inputs."""

from typing import Any, cast

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import Debouncer, clean_sequence, format_sequence

from .primer_input import PrimerInput
from .primer_row import PrimerRow
from .template_input import TemplateInput


class InputView(ft.Row):  # type: ignore[misc]
    """Input view composing DNA Template input and Primer inputs."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
        on_change: Any | None = None,
        on_stop_editing: Any | None = None,
    ) -> None:
        """Initialize the InputView.

        Args:
            page: The Flet page instance.
            input_data: The central input data store.
            settings: The GUI settings instance.
            on_change: Callback for input field changes.
            on_stop_editing: Callback for when editing stops.
        """
        super().__init__(
            expand=True, vertical_alignment=ft.CrossAxisAlignment.STRETCH
        )
        self.app_page = page
        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()
        self.on_change = on_change
        self.on_stop_editing_callback = on_stop_editing
        self._focus_debouncer = Debouncer(delay_seconds=0.15)
        self._currently_focused_control: ft.Control | None = None

        self.template_input = TemplateInput(
            page=self.app_page,
            settings=self.settings,
            input_data=self.input_data,
            on_change_handler=self._on_change_handler,
            handle_field_focus=self._handle_field_focus,
            handle_field_blur=self._handle_field_blur,
            handle_field_submit=self._handle_field_submit,
            clear_template_callback=self._clear_template,
        )

        self.primer_input = PrimerInput(
            page=self.app_page,
            settings=self.settings,
            input_data=self.input_data,
            on_change_handler=self._on_change_handler,
            handle_field_focus=self._handle_field_focus,
            handle_field_blur=self._handle_field_blur,
            handle_field_submit=self._handle_field_submit,
            clear_primers_callback=self._clear_primers,
            delete_selected_callback=self._delete_selected_primers,
        )

        self.right_fraction = 0.5

        self.divider = ft.GestureDetector(
            on_pan_update=self._on_pan_update,
            on_pan_end=self._on_pan_end,
            content=ft.Container(
                width=5,
                bgcolor=GUIColours.DIVIDER_GREY,
                border_radius=5,
                margin=ft.Margin.symmetric(horizontal=5),
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
        )

        self.controls = [
            self.template_input,
            self.divider,
            self.primer_input,
        ]

        # Sync initial UI state
        self.update_ui()
        self.app_page.on_resize = self._handle_resize

    # Properties/Delegates for backward compatibility and test access
    @property
    def template_sequence(self) -> ft.TextField:
        """Get the template sequence text field."""
        return self.template_input.template_sequence

    @property
    def template_circular(self) -> ft.Checkbox:
        """Get the circular template checkbox."""
        return self.template_input.template_circular

    @property
    def clear_template_button(self) -> ft.OutlinedButton:
        """Get the clear template button."""
        return self.template_input.clear_template_button

    @property
    def primers_list(self) -> ft.ListView:
        """Get the list of primers view."""
        return self.primer_input.primers_list

    @property
    def clear_primers_button(self) -> ft.OutlinedButton:
        """Get the clear primers button."""
        return self.primer_input.clear_primers_button

    @property
    def delete_selected_button(self) -> ft.OutlinedButton:
        """Get the delete selected primers button."""
        return self.primer_input.delete_selected_button

    @property
    def primer_info_panel(self) -> ft.Card:
        """Get the primer info panel."""
        return self.primer_input.primer_info_panel

    @property
    def info_header(self) -> ft.Container:
        """Get the primer info header container."""
        return self.primer_input.info_header

    @property
    def info_seq_text(self) -> ft.Text:
        """Get the primer info sequence text control."""
        return self.primer_input.info_seq_text

    @property
    def info_tm_text(self) -> ft.Text:
        """Get the primer info Tm text control."""
        return self.primer_input.info_tm_text

    @property
    def info_pairs_text(self) -> ft.Text:
        """Get the primer info AT/GC pairs text control."""
        return self.primer_input.info_pairs_text

    @property
    def info_redundancy_text(self) -> ft.Text:
        """Get the primer info redundancy text control."""
        return self.primer_input.info_redundancy_text

    @property
    def info_dimer_text(self) -> ft.Text:
        """Get the primer info dimer text control."""
        return self.primer_input.info_dimer_text

    @property
    def validation_errors(self) -> list[dict[str, str | None]]:
        """Get the list of validation errors."""
        return self.primer_input.validation_errors

    @validation_errors.setter
    def validation_errors(self, val: list[dict[str, str | None]]) -> None:
        """Set the list of validation errors."""
        self.primer_input.validation_errors = val

    @property
    def focused_primer_index(self) -> int | None:
        """Get the currently focused primer index."""
        return self.primer_input.focused_primer_index

    @focused_primer_index.setter
    def focused_primer_index(self, val: int | None) -> None:
        """Set the currently focused primer index."""
        self.primer_input.focused_primer_index = val

    def _handle_field_focus(self, e: ft.ControlEvent) -> None:
        """Handle focus on input fields to cancel auto-trigger timer.

        Updates the focused primer index, marks field touch status,
        highlights rows, and updates the primer info panel.

        Args:
            e: The Flet control event containing focus information.
        """
        self._focus_debouncer.cancel()

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

            self.primer_input.focused_primer_index = idx

            # Set touched status in state
            if 0 <= idx < len(self.input_data.primers):
                if field == "name":
                    self.input_data.primers[idx]["name_touched"] = True
                elif field == "seq":
                    self.input_data.primers[idx]["seq_touched"] = True

            self.primer_input._update_row_highlights()
            self.primer_input._update_primer_info_panel()
            if self.app_page:
                self.app_page.update()
        self._currently_focused_control = cast(ft.Control, e.control)

    def _handle_field_blur(self, e: ft.ControlEvent) -> None:
        """Handle blur on input fields to trigger results page after a delay.

        Syncs state, applies validation errors to rows, auto-adds empty
        rows if needed, and triggers a debounced callback for results.

        Args:
            e: The Flet control event containing blur information.
        """
        # If focus has moved to another input control, just sync state.
        if (
            self._currently_focused_control is not None
            and self._currently_focused_control != e.control
        ):
            self.sync_to_state(rebuild_if_needed=False)
            return

        self._currently_focused_control = None

        self.sync_to_state(rebuild_if_needed=False)

        if isinstance(e.control, ft.TextField) and e.control.data is not None:
            idx = (
                e.control.data["idx"]
                if isinstance(e.control.data, dict)
                else e.control.data
            )
            if idx < len(self.primer_input.validation_errors):
                err = self.primer_input.validation_errors[idx]
                for row in self.primer_input.primers_list.controls:
                    if isinstance(row, PrimerRow) and row.idx == idx:
                        row.set_error(err)
                        break
                # Auto-add new empty row if valid
                self._auto_add_empty_row_if_needed(cast(ft.Control, e.control))
                self.app_page.update()

        def timer_callback() -> None:
            """Callback for debounced field blur."""
            if not self.page:
                return
            if self.on_stop_editing_callback:
                self.on_stop_editing_callback()

        self._focus_debouncer.trigger(timer_callback)

    def _handle_field_submit(self, e: ft.ControlEvent) -> None:
        """Handle submission (Enter key) to immediately trigger results.

        Cancels any pending debounced callbacks, syncs state, auto-adds
        empty rows if needed, and triggers the stop-editing callback.

        Args:
            e: The Flet control event containing submission information.
        """
        self._focus_debouncer.cancel()
        self.sync_to_state()
        self._auto_add_empty_row_if_needed(cast(ft.Control, e.control))
        if self.app_page:
            self.app_page.update()
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback()

    def _auto_add_empty_row_if_needed(self, control: ft.Control) -> None:
        """Append a new empty row if the last row is filled and valid.

        Checks if the control that triggered this method is a sequence
        field on the last primer row, and if so, verifies that both
        name and sequence are filled and valid. If so, appends a new
        empty primer row.

        Args:
            control: The Flet control that triggered the check.
        """
        if (
            control.data is not None
            and isinstance(control.data, dict)
            and control.data.get("field") == "seq"
        ):
            idx = control.data["idx"]
            num_primers = len(self.input_data.primers)
            if idx == num_primers - 1:
                p = self.input_data.primers[idx]
                if p.get("name", "").strip() and p.get("seq", "").strip():
                    if idx < len(self.primer_input.validation_errors):
                        err = self.primer_input.validation_errors[idx]
                        if not err.get("name") and not err.get("seq"):
                            self.input_data.primers.append(
                                {"name": "", "seq": "", "active": False}
                            )
                            self.primer_input.update_ui()

    def sync_to_state(self, rebuild_if_needed: bool = False) -> None:
        """Sync current UI controls back to the central state.

        Args:
            rebuild_if_needed: Whether to rebuild UI if changes detected.
        """
        self.template_input.sync_to_state()
        self.primer_input.sync_to_state(rebuild_if_needed=rebuild_if_needed)

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state.

        Refreshes both the template input and primer input views to
        reflect the current state.
        """
        self.template_input.update_ui()
        self.primer_input.update_ui()

    def _on_change_handler(self, e: ft.ControlEvent) -> None:
        """Handle change in input fields.

        Syncs state to the central store, updates the primer info panel,
        and triggers the on_change callback.

        Args:
            e: The Flet control event containing change information.
        """
        self.sync_to_state()
        self.primer_input._update_primer_info_panel()
        if self.on_change:
            self.on_change(e)

    def _clear_primers(self, e: ft.ControlEvent) -> None:
        """Clear all primers.

        Resets the primer list to a single empty primer row, clears
        the focused primer index, updates the UI, and triggers callbacks.

        Args:
            e: The Flet control event containing click information.
        """
        self.input_data.primers = [{"name": "", "seq": "", "active": False}]
        self.primer_input.focused_primer_index = None
        self.update_ui()
        if self.on_change:
            self.on_change(e)
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback()

    def _delete_selected_primers(self, e: ft.ControlEvent) -> None:
        """Delete all selected primers.

        Identifies all active (selected) primers and removes them using
        the action controller, then triggers callbacks.

        Args:
            e: The Flet control event containing click information.
        """
        active_indices = {
            i for i, p in enumerate(self.input_data.primers) if p.get("active")
        }
        self.primer_input.action_controller.delete_primers(active_indices)
        if self.on_change:
            self.on_change(e)
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback()

    def _clear_template(self, e: ft.ControlEvent) -> None:
        """Clear the DNA template.

        Clears the template sequence field, syncs state, and triggers
        callbacks.

        Args:
            e: The Flet control event containing click information.
        """
        self.template_input.template_sequence.value = ""
        self.sync_to_state()
        if self.on_change:
            self.on_change(e)
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback()

    def _on_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle resizing the bottom (right) container via the divider.

        Adjusts the width of the template input panel and the name
        column width of primer rows based on drag delta.

        Args:
            e: The Flet drag update event containing delta information.
        """
        page_width = self.app_page.width
        if isinstance(page_width, (int, float)) and page_width > 0:
            delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
            # Calculate current pixel width of the right container
            current_width = page_width * self.right_fraction
            new_width = max(200.0, current_width - delta_x)
            # Ensure the left container stays at least 200px wide
            new_width = min(new_width, page_width - 200.0)

            # Recalculate fractions
            self.right_fraction = new_width / page_width

            # Resize template_input via fixed width to avoid updating
            # the massive primer_input ListView controls tree.
            self.template_input.expand = None
            self.template_input.width = page_width - new_width - 5.0
            self.template_input.update()

            # Adjust the name column width of only the visible rows and header
            self.primer_input.layout_manager.adjust_name_column_width(
                new_width, during_drag=True
            )

    def _on_pan_end(self, e: ft.DragEndEvent) -> None:
        """Handle finishing the drag of the main layout divider.

        Restores responsive expand weights for both panels, performs a
        full rebuild to sync everything at the final width, and updates
        the name column width.

        Args:
            e: The Flet drag end event.
        """
        page_width = self.app_page.width
        if isinstance(page_width, (int, float)) and page_width > 0:
            panel_width = page_width * self.right_fraction

            # Restore responsive expand weights for both panels
            self.template_input.width = None
            self.template_input.expand = int((1.0 - self.right_fraction) * 1000)
            self.primer_input.expand = int(self.right_fraction * 1000)

            # Full rebuild to sync everything at the final width
            self.primer_input.layout_manager.adjust_name_column_width(
                panel_width, during_drag=False
            )
            self.update()

    def _handle_resize(self, e: Any) -> None:
        """Handle page resize to proportionally scale name column.

        Args:
            e: The Flet resize event.
        """
        page_width = self.app_page.width
        if isinstance(page_width, (int, float)) and page_width > 0:
            panel_width = page_width * self.right_fraction
            self.primer_input.layout_manager.adjust_name_column_width(
                panel_width
            )
            self.update()

    def get_primers(self) -> list[dict[str, Any]]:
        """Get the list of active primers.

        Returns:
            A list of dictionaries containing active primer data.
        """
        return self.input_data.get_active_primers()

    def get_all_primers_state(self) -> list[dict[str, Any]]:
        """Get all primers (active and inactive) for serialisation.

        Returns:
            A list of dictionaries containing all primer data with name,
            formatted sequence, and active status.
        """
        primers: list[dict[str, Any]] = []
        for p in self.input_data.primers:
            if (
                not str(p.get("name", "")).strip()
                and not clean_sequence(str(p.get("seq", ""))).strip()
            ):
                continue
            primers.append(
                {
                    "name": p["name"],
                    "seq": format_sequence(p["seq"]),
                    "active": p.get("active", True),
                }
            )
        return primers

    def get_state(self) -> dict[str, Any]:
        """Get the current input data state for serialisation.

        Syncs UI controls to state and returns the central input data
        as a dictionary.

        Returns:
            A dictionary containing the current input data state.
        """
        self.sync_to_state()
        return self.input_data.to_dict()

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current input data from deserialized data.

        Args:
            state: A dictionary containing the input data to restore.
        """
        self.input_data.from_dict(state)
        self.update_ui()
        self.app_page.update()

    def _update_row_highlights(self) -> None:
        """Update background colours of all row containers.

        Delegates to the primer input's row highlight update method.
        """
        self.primer_input._update_row_highlights()

    def _update_primer_info_panel(self) -> None:
        """Update the primer information panel based on the focused primer.

        Delegates to the primer input's info panel update method.
        """
        self.primer_input._update_primer_info_panel()
