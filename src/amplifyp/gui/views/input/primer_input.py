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

"""Input component for DNA primers list and details."""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import Any, cast

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.sequence import clean_sequence
from amplifyp.gui.utils.ui import NotificationHelper

from .primer_action_controller import PrimerActionController
from .primer_file_manager import PrimerFileManager
from .primer_header import PrimerHeader
from .primer_info_panel import PrimerInfoPanel
from .primer_layout_manager import PrimerLayoutManager
from .primer_list import PrimerList
from .primer_row import PrimerRow
from .primer_toolbar import PrimerToolbar
from .primer_validation import (
    get_duplicate_primer_indices,
    reconcile_primer_states,
    validate_primers,
)

logger = logging.getLogger(__name__)


class PrimerInput(ft.Container):  # type: ignore[misc]
    """Input component for DNA primers list and details."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        input_data: GUIInput,
        on_change_handler: Callable[[ft.Event | None], None],
        handle_field_focus: Callable[[ft.Event[ft.TextField]], None],
        handle_field_blur: Callable[[ft.Event[ft.TextField]], None],
        handle_field_submit: Callable[[ft.Event[ft.TextField]], None],
        clear_primers_callback: Callable[[ft.Event | None], None],
        delete_selected_callback: Callable[[ft.Event | None], None],
    ) -> None:
        """Initialise the PrimerInput component.

        Args:
            page: The Flet page instance for file picker and notifications.
            settings: Application GUI settings instance.
            input_data: Central GUI input state object.
            on_change_handler: Callback invoked when primer content changes.
            handle_field_focus: Callback for field focus events.
            handle_field_blur: Callback for field blur events.
            handle_field_submit: Callback for field submit events.
            clear_primers_callback: Callback to clear all primers.
            delete_selected_callback: Callback to delete selected primers.
        """
        super().__init__(expand=5)
        self.app_page = page
        self.settings = settings
        self.input_data = input_data
        self.on_change_handler = on_change_handler
        self.handle_field_focus = handle_field_focus
        self.handle_field_blur = handle_field_blur
        self.handle_field_submit = handle_field_submit
        self.clear_primers_callback = clear_primers_callback
        self.delete_selected_callback = delete_selected_callback

        self.layout_manager = PrimerLayoutManager(self)
        self.action_controller = PrimerActionController(self)

        self.focused_primer_index: int | None = None
        self.selected_indices: set[int] = set()
        self.validation_errors: list[dict[str, str | None]] = []
        self._prev_header_checkbox_value: bool | None = None
        self._visible_rows_cache: list[PrimerRow] | None = None

        font_family = self.settings.get("font_family", "Roboto Mono")
        self.name_column_width = 150.0

        # Calculate initial name column ratio
        self.name_column_ratio = (
            self.name_column_width / self.layout_manager.get_panel_width()
        )

        # Primer List Component
        self.primers_list = PrimerList(self)

        # Primer Header Component
        self.primer_header = PrimerHeader(
            settings=self.settings,
            on_toggle_all=self._on_toggle_all_primers,
            on_divider_pan=self.layout_manager.on_primer_divider_pan,
            on_divider_pan_end=self.layout_manager.on_primer_divider_pan_end,
            name_column_width=self.name_column_width,
            on_add_primer=self._header_add_click,
            on_delete_primer=self._header_delete_click,
            on_move_primer_up=self._header_up_click,
            on_move_primer_down=self._header_down_click,
        )
        # Compatibility links
        self.all_primers_checkbox = self.primer_header.all_primers_checkbox
        self.primers_header = self.primer_header.header_row
        self.primers_header_container = self.primer_header

        # File Manager Component
        self.file_manager = PrimerFileManager(
            page=self.app_page,
            input_data=self.input_data,
            on_update_ui=self.update_ui,
            on_change_handler=self.on_change_handler,
            show_notification=self._show_notification,
        )

        # Primer Toolbar Component
        self.primer_toolbar = PrimerToolbar(
            on_save=self._save_primers_click,
            on_load=self._load_primers_click,
            on_clear=self.clear_primers_callback,
            on_delete_selected=self.delete_selected_callback,
            on_copy=self._copy_primers_click,
            on_paste=self._paste_primers_click,
        )
        # Compatibility links
        self.clear_primers_button = self.primer_toolbar.clear_button
        self.delete_selected_button = self.primer_toolbar.delete_selected_button

        self.error_message_text = ft.Text(
            value=(
                "PCR and Primer Dimer views are disabled because "
                "one or more selected primers are invalid, or "
                "have duplicated names/sequences."
            ),
            color=ft.Colors.ON_ERROR_CONTAINER,
            weight=ft.FontWeight.BOLD,
            size=13,
        )
        self.error_banner = ft.Container(
            content=self.error_message_text,
            padding=ft.Padding(10, 5, 10, 5),
            bgcolor=ft.Colors.ERROR_CONTAINER,
            border_radius=5,
            visible=False,
        )

        # Primer Info Panel
        self.primer_info_panel = PrimerInfoPanel(
            settings=self.settings, font_family=font_family
        )

        # Primer List Container (stored for dynamic repositioning)
        self.primer_list_container = ft.Container(
            content=ft.Column(
                [
                    self.primers_header_container,
                    self.primers_list,
                ],
                expand=True,
                spacing=0,
            ),
            expand=True,
            border=ft.Border.all(1, GUIColours.OUTLINE),
            border_radius=5,
            padding=0,
        )

        # Title row (stored for dynamic repositioning)
        self.primer_title_row = ft.ResponsiveRow(
            [
                ft.Row(
                    [
                        ft.Text(
                            "Primers",
                            weight=ft.FontWeight.BOLD,
                            no_wrap=True,
                        ),
                        self.primer_toolbar,
                    ],
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    vertical_alignment=ft.CrossAxisAlignment.CENTER,
                    wrap=True,
                ),
            ],
            run_spacing=0,
        )

        self.content = ft.Column(
            [
                self.primer_title_row,
                self.primer_info_panel,
                self.primer_list_container,
                self.error_banner,
            ],
            expand=True,
            spacing=5,
        )

        # Position info panel based on setting
        self._reposition_info_panel()

    def _reposition_info_panel(self) -> None:
        """Reposition the primer info panel based on the current setting.

        Moves the info panel to either the top (after the toolbar row)
        or the bottom (before the error banner) of the primer list area.

        Rebuilds the entire Column so Flet detects the structural change
        and re-renders the new order.
        """
        position = str(
            self.settings.get("primer_info_panel_position", "bottom")
        ).lower()

        if position == "top":
            new_controls: list[ft.Control] = [
                self.primer_title_row,
                self.primer_info_panel,
                self.primer_list_container,
                self.error_banner,
            ]
        else:
            new_controls = [
                self.primer_title_row,
                self.primer_list_container,
                self.primer_info_panel,
                self.error_banner,
            ]

        if isinstance(self.content, ft.Column):
            self.content.controls = new_controls
        else:
            self.content = ft.Column(
                new_controls,
                expand=True,
                spacing=5,
            )
        # Trigger Flet to re-render the changed content
        try:
            self.update()
        except RuntimeError:
            pass

    def reposition_info_panel(self) -> None:
        """Public method to reposition the info panel.

        Call this when the setting changes to ensure the panel is
        immediately repositioned. The internal update is handled by
        _reposition_info_panel.
        """
        self._reposition_info_panel()

    # Delegate properties of the info panel for test and backward compatibility
    @property
    def info_header(self) -> ft.Container:
        """Get the primer info header container."""
        return self.primer_info_panel.info_header

    @property
    def info_seq_text(self) -> ft.Text:
        """Get the primer info sequence text control."""
        return self.primer_info_panel.info_seq_text

    @property
    def info_tm_text(self) -> ft.Text:
        """Get the primer info Tm text control."""
        return self.primer_info_panel.info_tm_text

    @property
    def info_pairs_text(self) -> ft.Text:
        """Get the primer info AT/GC pairs text control."""
        return self.primer_info_panel.info_pairs_text

    @property
    def info_redundancy_text(self) -> ft.Text:
        """Get the primer info redundancy text control."""
        return self.primer_info_panel.info_redundancy_text

    @property
    def info_dimer_text(self) -> ft.Text:
        """Get the primer info dimer text control."""
        return self.primer_info_panel.info_dimer_text

    def _extract_primer_data_from_ui(self) -> list[dict[str, Any]]:
        """Extract primer data from UI controls.

        Returns:
            A list of dicts with keys 'name', 'seq', 'active', 'container',
            and 'checkbox' for each primer row in the UI.
        """
        ui_primers = []
        for row in self.primers_list.controls:
            if not isinstance(row, PrimerRow):
                continue

            name_val = str(row.name_field.value or "").strip()
            seq_val = clean_sequence(str(row.seq_field.value or ""))
            is_active = bool(row.checkbox.value)

            ui_primers.append(
                {
                    "name": name_val,
                    "seq": seq_val,
                    "active": is_active,
                    "container": row,
                    "checkbox": row.checkbox,
                }
            )
        return ui_primers

    def sync_to_state(
        self, rebuild_if_needed: bool = False, skip_extract: bool = False
    ) -> bool:
        """Sync current UI controls back to the central state.

        Args:
            rebuild_if_needed: If True, triggers a full UI rebuild after
                syncing. Otherwise, updates error states in-place.
            skip_extract: If True, skips UI extraction and uses existing
                state directly. Used after paste handling where state is
                already updated.

        Returns:
            True if a UI rebuild was needed and performed.
        """
        if skip_extract:
            ui_primers = self.input_data.primers
        else:
            ui_primers = self._extract_primer_data_from_ui()

        auto_activate_new = self.settings.get(
            "auto_activate_new_valid_primer", False
        )
        primers = reconcile_primer_states(
            ui_primers,
            self.input_data.primers,
            auto_activate_new=auto_activate_new,
        )

        # Update checkbox values in-place on UI controls if they
        # were updated during reconciliation.
        for reconciled_p, ui_p in zip(primers, ui_primers, strict=True):
            checkbox = ui_p.get("checkbox")
            if checkbox:
                checkbox.value = reconciled_p["active"]

        ignore_inactive_name_dup = self.settings.get(
            "ignore_inactive_name_dup_warn", True
        )
        ignore_inactive_seq_dup = self.settings.get(
            "ignore_inactive_seq_dup_warn", True
        )

        dup_indices = get_duplicate_primer_indices(
            ui_primers, ignore_inactive_name_dup, ignore_inactive_seq_dup
        )
        for p in ui_primers:
            container = p.get("container")
            if container is None:
                continue
            c_idx = container.data
            is_dup = c_idx in dup_indices
            new_color = GUIColours.DUPLICATE_BG if is_dup else None
            if container.bgcolor != new_color:
                container.bgcolor = new_color

        # Run background primer construction/validation
        new_validation_errors = validate_primers(
            primers, ignore_inactive_name_dup, ignore_inactive_seq_dup
        )

        self.input_data.primers = primers
        self.validation_errors = new_validation_errors
        if rebuild_if_needed:
            self.update_ui()
        else:
            # Update error status and duplicate highlights in-place on
            # existing controls
            for idx, row in enumerate(self.primers_list.controls):
                if isinstance(row, PrimerRow) and idx < len(
                    new_validation_errors
                ):
                    row.set_error(new_validation_errors[idx])
                    row.update_tm(self.settings)
            self._update_row_highlights()
            self._update_header_checkbox_state()
            self._update_delete_button_disabled_state()

        return rebuild_if_needed

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state.

        Rebuilds the primer list and updates the delete button disabled
        state based on current selections.
        """
        # Recreate the header to make sure controls match settings and indices
        self.primer_header = PrimerHeader(
            settings=self.settings,
            on_toggle_all=self._on_toggle_all_primers,
            on_divider_pan=self.layout_manager.on_primer_divider_pan,
            on_divider_pan_end=self.layout_manager.on_primer_divider_pan_end,
            name_column_width=self.name_column_width,
            on_add_primer=self._header_add_click,
            on_delete_primer=self._header_delete_click,
            on_move_primer_up=self._header_up_click,
            on_move_primer_down=self._header_down_click,
        )
        self.all_primers_checkbox = self.primer_header.all_primers_checkbox
        self.primers_header = self.primer_header.header_row
        self.primers_header_container = self.primer_header

        # Reposition info panel based on current setting
        self._reposition_info_panel()

        # Replace the header in the UI container controls
        inner_column = cast(ft.Column, self.primer_list_container.content)
        inner_column.controls[0] = self.primer_header
        try:
            if self.primer_list_container.page:
                self.primer_list_container.update()
        except RuntimeError:
            logger.debug("Container page detached, skipping update")

        self.primers_list.update_list_ui()
        self._update_delete_button_disabled_state()
        self._update_header_buttons_state()
        try:
            if self.page:
                self.update()
        except RuntimeError:
            pass

    def _update_delete_button_disabled_state(self) -> None:
        """Update disabled state of the delete button based on selection.

        The delete button is enabled only when at least one primer is
        active (selected).
        """
        self.delete_selected_button.disabled = not self.selected_indices

    def _on_toggle_all_primers(self, e: ft.Event[ft.Checkbox]) -> None:
        """Toggle all primers active/inactive based on tri-state checkbox.

        Handles the three checkbox states (True, False, None/tristate)
        to determine the target active state for all primers.

        Args:
            e: The Flet control event triggered by the checkbox change.
        """
        primers = self.input_data.primers
        if not primers:
            return

        cb_value = self.all_primers_checkbox.value
        prev_value = self._prev_header_checkbox_value

        if cb_value is True:
            target_active = True
        elif cb_value is False:
            target_active = False
        else:
            if prev_value is True:
                target_active = False
            elif prev_value is False:
                target_active = True
            else:
                target_active = True

        # Update checkbox values in-place on existing UI controls
        for row in self.primers_list.controls:
            if isinstance(row, PrimerRow) and row.data is not None:
                idx = row.data
                if idx < len(primers):
                    if target_active:
                        if not row.checkbox.disabled:
                            row.checkbox.value = True
                    else:
                        row.checkbox.value = False

        self._prev_header_checkbox_value = cb_value
        self.on_change_handler(e)

    def _update_header_checkbox_state(self) -> None:
        """Update the header checkbox to reflect the current primer states.

        Sets the tri-state checkbox to True (all active), False (none
        active), or None (mixed state) based on the current primer list.
        """
        primers = self.input_data.primers
        non_empty = [
            p
            for p in primers
            if str(p.get("name", "")).strip()
            or clean_sequence(str(p.get("seq", ""))).strip()
        ]
        if not non_empty:
            self.all_primers_checkbox.value = None
        elif all(p.get("active", True) for p in non_empty):
            self.all_primers_checkbox.value = True
        elif all(not p.get("active", True) for p in non_empty):
            self.all_primers_checkbox.value = False
        else:
            self.all_primers_checkbox.value = None
        self._prev_header_checkbox_value = self.all_primers_checkbox.value
        if self.app_page:
            self.app_page.update()

    def _get_duplicate_indices(self) -> set[int]:
        """Find indices of primers with duplicate names or sequences.

        Returns:
            Set of indices corresponding to primers with duplicate names
            or sequences in the central state.
        """
        ignore_inactive_name_dup = self.settings.get(
            "ignore_inactive_name_dup_warn", True
        )
        ignore_inactive_seq_dup = self.settings.get(
            "ignore_inactive_seq_dup_warn", True
        )
        return get_duplicate_primer_indices(
            self.input_data.primers,
            ignore_inactive_name_dup,
            ignore_inactive_seq_dup,
        )

    def _update_row_highlights(self) -> None:
        """Update background colours of all row containers.

        Highlights rows based on selection (focused primer) and
        duplicates (by name or sequence).
        """
        self.primers_list.update_row_highlights()
        self._update_header_buttons_state()

    def _update_primer_info_panel(self) -> None:
        """Update the primer information panel based on the focused primer.

        Displays Tm, sequence details, base counts, redundancy, and
        dimer potential for the currently focused primer row.
        """

        def on_dismiss() -> None:
            """Clear focus when the primer info panel is dismissed."""
            self.focused_primer_index = None
            self._update_row_highlights()
            if self.app_page:
                self.app_page.update()

        self.primer_info_panel.update_panel(
            focused_idx=self.focused_primer_index,
            input_data=self.input_data,
            app_page=self.app_page,
            on_update_highlights=self._update_row_highlights,
            on_dismiss=on_dismiss,
        )

    async def _load_primers_click(self, e: ft.Event | None) -> None:
        """Open file picker to load primers from CSV/TSV file.

        Delegates to the PrimerFileManager component.

        Args:
            e: The Flet control event triggered by the load button click.
        """
        await self.file_manager.load_primers_click(e)

    async def _save_primers_click(self, e: ft.Event | None) -> None:
        """Save active primers to a CSV file.

        Delegates to the PrimerFileManager component.

        Args:
            e: The Flet control event triggered by the save button click.
        """
        await self.file_manager.save_primers_click(e)

    def _show_notification(self, message: str) -> None:
        """Show a notification message.

        Args:
            message: The notification message text to display.
        """
        if not hasattr(self, "_notification_helper"):
            self._notification_helper = NotificationHelper(self.app_page)
        self._notification_helper.show_message(message)

    async def _copy_primers_click(self, e: ft.Event | None) -> None:
        """Copy selected or focused primers to clipboard in TSV format."""
        primers = self.input_data.primers
        selected_primers = (
            [
                primers[i]
                for i in sorted(self.selected_indices)
                if 0 <= i < len(primers)
            ]
            if self.selected_indices
            else []
        )

        # Fallback to focused primer if no selected rows
        if not selected_primers and self.focused_primer_index is not None:
            if 0 <= self.focused_primer_index < len(primers):
                selected_primers = [primers[self.focused_primer_index]]

        if not selected_primers:
            self._show_notification(
                "No primers highlighted or focused to copy."
            )
            return

        lines = []
        for p in selected_primers:
            name = str(p.get("name", "")).strip()
            seq = str(p.get("seq", "")).strip()
            lines.append(f"{name}\t{seq}")

        tsv_text = "\n".join(lines)
        await ft.Clipboard().set(tsv_text)
        self._show_notification(
            f"Copied {len(selected_primers)} primer(s) to clipboard."
        )

    async def _paste_primers_click(self, e: ft.Event | None) -> None:
        """Paste primers from clipboard starting at focused index or end."""
        try:
            clipboard_text = await ft.Clipboard().get()
        except Exception as ex:
            logger.warning("Failed to access clipboard: %s", ex)
            self._show_notification(
                "Unable to read clipboard. Try pasting directly into a "
                "primer field using Ctrl+V."
            )
            return
        if not clipboard_text:
            self._show_notification("Clipboard is empty.")
            return

        from .primer_clipboard import parse_primer_clipboard_text

        parsed = parse_primer_clipboard_text(clipboard_text)
        if not parsed:
            self._show_notification("No valid primers found in clipboard.")
            return

        primers = self.input_data.primers
        valid_selected = {
            i for i in self.selected_indices if 0 <= i < len(primers)
        }
        if valid_selected:
            insert_idx = max(valid_selected) + 1
        else:
            insert_idx = len(primers)

        # Replace single empty row
        if (
            len(primers) == 1
            and not primers[0].get("name")
            and not primers[0].get("seq")
        ):
            primers.clear()
            insert_idx = 0

        for i, new_p in enumerate(parsed):
            new_p["active"] = False
            primers.insert(insert_idx + i, new_p)

        self.selected_indices = set(range(insert_idx, insert_idx + len(parsed)))
        self.focused_primer_index = insert_idx

        self.sync_to_state(rebuild_if_needed=True, skip_extract=True)
        if self.on_change_handler is not None:
            self.on_change_handler(None)

        self._show_notification(f"Pasted {len(parsed)} primer(s).")

    def _update_header_buttons_state(self) -> None:
        """Update enabled/disabled state of header action buttons."""
        if not hasattr(self, "primer_header"):
            return

        num_primers = len(self.input_data.primers)
        has_sel = bool(self.selected_indices)

        self.primer_header.add_button.disabled = False
        self.primer_header.delete_button.disabled = not has_sel

        if has_sel:
            min_idx = min(self.selected_indices)
            max_idx = max(self.selected_indices)
            self.primer_header.up_button.disabled = min_idx == 0
            self.primer_header.down_button.disabled = max_idx == num_primers - 1
        else:
            self.primer_header.up_button.disabled = True
            self.primer_header.down_button.disabled = True

        try:
            if self.primer_header.page:
                self.primer_header.update()
        except RuntimeError:
            logger.debug("Header page detached, skipping update")

    def _header_add_click(self, e: ft.Event | None) -> None:
        """Handle header Add button click."""
        num_primers = len(self.input_data.primers)
        if self.selected_indices:
            idx = min(max(self.selected_indices), num_primers - 1)
        elif self.focused_primer_index is not None:
            idx = min(self.focused_primer_index, num_primers - 1)
        else:
            idx = num_primers - 1

        self.selected_indices = {idx + 1}
        self.focused_primer_index = idx + 1
        self.action_controller.on_add_primer_row(idx)

    def _header_delete_click(self, e: ft.Event | None) -> None:
        """Handle header Delete button click."""
        if self.selected_indices:
            self.action_controller.delete_primers(self.selected_indices.copy())
            self._update_header_buttons_state()

    def _header_up_click(self, e: ft.Event | None) -> None:
        """Handle header Move Up button click."""
        if self.selected_indices and min(self.selected_indices) > 0:
            self.action_controller.move_primers(self.selected_indices, -1)
            self._update_header_buttons_state()

    def _header_down_click(self, e: ft.Event | None) -> None:
        """Handle header Move Down button click."""
        if (
            self.selected_indices
            and max(self.selected_indices) < len(self.input_data.primers) - 1
        ):
            self.action_controller.move_primers(self.selected_indices, 1)
            self._update_header_buttons_state()
