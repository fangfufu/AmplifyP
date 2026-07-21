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

"""Input view composing DNA Template input and Primer inputs."""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.data_helpers import clean_sequence, format_sequence
from amplifyp.gui.utils.gui_helpers import Debouncer

from .primer.input import PrimerInput
from .template.input import TemplateInput

logger = logging.getLogger(__name__)


class InputView(ft.Row):  # type: ignore[misc]
    """Input view composing DNA Template input and Primer inputs."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
        on_change: Callable[[ft.Event | None], None] | None = None,
        on_stop_editing: Callable[[ft.Event | None], None] | None = None,
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
    def upper_case_button(self) -> ft.OutlinedButton:
        """Get the upper case button."""
        return self.template_input.upper_case_button

    @property
    def lower_case_button(self) -> ft.OutlinedButton:
        """Get the lower case button."""
        return self.template_input.lower_case_button

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

    def _handle_field_focus(self, e: ft.Event[ft.TextField]) -> None:
        """Handle focus on input fields to cancel auto-trigger timer."""
        from .events import handle_field_focus

        handle_field_focus(self, e)

    def _handle_field_blur(self, e: ft.Event[ft.TextField]) -> None:
        """Handle blur on input fields to trigger results page after a delay."""
        from .events import handle_field_blur

        handle_field_blur(self, e)

    def _handle_field_submit(self, e: ft.Event[ft.TextField]) -> None:
        """Handle submission (Enter key) to immediately trigger results."""
        from .events import handle_field_submit

        handle_field_submit(self, e)

    def _auto_add_empty_row_if_needed(self, control: ft.Control) -> None:
        """Append a new empty row if the last row is filled and valid."""
        from .events import auto_add_empty_row_if_needed

        auto_add_empty_row_if_needed(self, control)

    def sync_to_state(
        self, rebuild_if_needed: bool = False, skip_extract: bool = False
    ) -> None:
        """Sync current UI controls back to the central state.

        Args:
            rebuild_if_needed: Whether to rebuild UI if changes detected.
            skip_extract: Whether to skip UI extraction for primer input.
        """
        self.template_input.sync_to_state()
        self.primer_input.sync_to_state(
            rebuild_if_needed=rebuild_if_needed, skip_extract=skip_extract
        )

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state.

        Refreshes both the template input and primer input views to
        reflect the current state.
        """
        self.template_input.update_ui()
        self.primer_input.update_ui()

        self._adjust_template_wrap(update_first=False)

    def _on_change_handler(self, e: ft.Event | None) -> None:
        """Handle change in input fields."""
        from .events import on_change_handler

        on_change_handler(self, e)

    def _handle_pasted_text(
        self, text: str, idx: int, field: str, control: ft.TextField
    ) -> None:
        """Parse pasted text and insert into the primer list."""
        from .events import handle_pasted_text

        handle_pasted_text(self, text, idx, field, control)

    def _clear_primers(self, e: ft.Event | None) -> None:
        """Clear all primers."""
        self.input_data.primers = [{"name": "", "seq": "", "active": False}]
        self.primer_input.focused_primer_index = None
        self.update_ui()
        if self.on_change:
            self.on_change(e)
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback(None)

    def _delete_selected_primers(self, e: ft.Event | None) -> None:
        """Delete all selected primers."""
        active_indices = self.primer_input.selected_indices
        if active_indices:
            self.primer_input.action_controller.delete_primers(
                active_indices.copy()
            )
        if self.on_change:
            self.on_change(e)
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback(None)

    def _clear_template(self, e: ft.Event | None) -> None:
        """Clear the DNA template."""
        self.template_input.template_sequence.value = ""
        self.sync_to_state()
        if self.on_change:
            self.on_change(e)
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback(None)

    def _on_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle resizing the bottom (right) container via the divider."""
        from .layout import on_pan_update

        on_pan_update(self, e)

    def _on_pan_end(self, e: ft.DragEndEvent) -> None:
        """Handle finishing the drag of the main layout divider."""
        from .layout import on_pan_end

        on_pan_end(self, e)

    def _handle_resize(self, e: ft.PageResizeEvent) -> None:
        """Handle page resize to proportionally scale name column."""
        from .layout import handle_resize

        handle_resize(self, e)

    def _adjust_template_wrap(self, update_first: bool = True) -> None:
        """Adjust the template wrap length based on the available width."""
        from .layout import adjust_template_wrap

        adjust_template_wrap(self, update_first)

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

    def reposition_primer_info_panel(self) -> None:
        """Reposition the primer info panel based on the current setting.

        Delegates to the primer input's reposition method and triggers
        a page update.
        """
        self.primer_input.reposition_info_panel()
