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
from typing import cast

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.data_helpers import clean_sequence
from amplifyp.gui.utils.gui_helpers import NotificationHelper

from .action_controller import PrimerActionController
from .coordinator import PrimerCoordinator
from .file_manager import PrimerFileManager
from .header import PrimerHeader
from .info_panel import PrimerInfoPanel
from .layout_manager import PrimerLayoutManager
from .list import PrimerList
from .row import PrimerRow
from .toolbar import PrimerToolbar

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
        """Initialise the PrimerInput component."""
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
        self.coordinator = PrimerCoordinator(self)

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
            on_add_primer=self.action_controller.header_add_click,
            on_delete_primer=self.action_controller.header_delete_click,
            on_move_primer_up=self.action_controller.header_up_click,
            on_move_primer_down=self.action_controller.header_down_click,
        )
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
            on_reverse_complement=self._reverse_complement_selected_click,
        )
        self.clear_primers_button = self.primer_toolbar.clear_button
        self.delete_selected_button = self.primer_toolbar.delete_selected_button
        self.reverse_complement_button = (
            self.primer_toolbar.reverse_complement_button
        )

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

        # Primer List Container
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

        # Title row
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

        self._reposition_info_panel()

    def _reposition_info_panel(self) -> None:
        """Reposition the primer info panel based on the current setting."""
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
        try:
            self.update()
        except RuntimeError:
            pass

    def reposition_info_panel(self) -> None:
        """Public method to reposition the info panel."""
        self._reposition_info_panel()

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

    def sync_to_state(
        self, rebuild_if_needed: bool = False, skip_extract: bool = False
    ) -> bool:
        """Sync current UI controls back to the central state."""
        return self.coordinator.sync_to_state(rebuild_if_needed, skip_extract)

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        self.primer_header = PrimerHeader(
            settings=self.settings,
            on_toggle_all=self._on_toggle_all_primers,
            on_divider_pan=self.layout_manager.on_primer_divider_pan,
            on_divider_pan_end=self.layout_manager.on_primer_divider_pan_end,
            name_column_width=self.name_column_width,
            on_add_primer=self.action_controller.header_add_click,
            on_delete_primer=self.action_controller.header_delete_click,
            on_move_primer_up=self.action_controller.header_up_click,
            on_move_primer_down=self.action_controller.header_down_click,
        )
        self.all_primers_checkbox = self.primer_header.all_primers_checkbox
        self.primers_header = self.primer_header.header_row
        self.primers_header_container = self.primer_header

        self._reposition_info_panel()
        self._update_primer_info_panel()

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
        """Update disabled state of delete and rev comp buttons."""
        is_disabled = not self.selected_indices
        self.delete_selected_button.disabled = is_disabled
        self.reverse_complement_button.disabled = is_disabled

    def _reverse_complement_selected_click(self, _e: ft.Event | None) -> None:
        """Reverse complement highlighted primers."""
        if self.selected_indices:
            self.action_controller.reverse_complement_primers(
                self.selected_indices.copy()
            )

    def _on_toggle_all_primers(self, e: ft.Event[ft.Checkbox]) -> None:
        """Toggle all primers active/inactive based on tri-state checkbox."""
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
        """Update the header checkbox to reflect the current primer states."""
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
        """Find indices of primers with duplicate names or sequences."""
        return self.coordinator.get_duplicate_indices()

    def _update_row_highlights(self) -> None:
        """Update background colours of all row containers."""
        self.primers_list.update_row_highlights()
        self._update_header_buttons_state()

    def _update_primer_info_panel(self) -> None:
        """Update the primer information panel based on the focused primer."""

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
        """Open file picker to load primers from CSV/TSV file."""
        await self.file_manager.load_primers_click(e)

    async def _save_primers_click(self, e: ft.Event | None) -> None:
        """Save active primers to a CSV file."""
        await self.file_manager.save_primers_click(e)

    def _show_notification(self, message: str) -> None:
        """Show a notification message."""
        if not hasattr(self, "_notification_helper"):
            self._notification_helper = NotificationHelper(self.app_page)
        self._notification_helper.show_message(message)

    async def _copy_primers_click(self, e: ft.Event | None) -> None:
        """Copy selected or focused primers to clipboard in TSV format."""
        from .clipboard import copy_primers_click

        await copy_primers_click(self, e)

    async def _paste_primers_click(self, e: ft.Event | None) -> None:
        """Paste primers from clipboard starting at focused index or end."""
        from .clipboard import paste_primers_click

        await paste_primers_click(self, e)

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
