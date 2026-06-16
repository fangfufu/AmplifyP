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

"""Input component for DNA primers list and details."""

from typing import Any

import flet as ft

from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import NotificationHelper, clean_sequence

from .primer_action_controller import PrimerActionController
from .primer_file_manager import PrimerFileManager
from .primer_header import PrimerHeader
from .primer_info_panel import PrimerInfoPanel
from .primer_layout_manager import PrimerLayoutManager
from .primer_list import PrimerList
from .primer_row import PrimerRow
from .primer_toolbar import PrimerToolbar
from .primer_validation import validate_primers


class PrimerInput(ft.Container):  # type: ignore[misc]
    """Input component for DNA primers list and details."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        input_data: GUIInput,
        on_change_handler: Any,
        handle_field_focus: Any,
        handle_field_blur: Any,
        handle_field_submit: Any,
        clear_primers_callback: Any,
        delete_selected_callback: Any,
    ) -> None:
        """Initialize the PrimerInput component."""
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
        )
        # Compatibility links
        self.clear_primers_button = self.primer_toolbar.clear_button
        self.delete_selected_button = self.primer_toolbar.delete_selected_button

        # Error Banner for selected invalid primers
        self.error_message_text = ft.Text(
            value=(
                "PCR and Primer Dimer views are disabled because "
                "one or more selected (active) primers are invalid."
            ),
            color=GUIColors.ERROR_RED,
            weight=ft.FontWeight.BOLD,
            size=13,
        )
        self.error_banner = ft.Container(
            content=self.error_message_text,
            padding=ft.Padding(10, 5, 10, 5),
            bgcolor=GUIColors.DUPLICATE_BG,
            border_radius=5,
            visible=False,
        )

        # Primer Info Panel
        self.primer_info_panel = PrimerInfoPanel(
            settings=self.settings, font_family=font_family
        )

        self.content = ft.Column(
            [
                ft.ResponsiveRow(
                    [
                        ft.Container(
                            content=ft.Text(
                                "Primers",
                                weight=ft.FontWeight.BOLD,
                                no_wrap=True,
                            ),
                            col={"lg": 3, "md": 3, "sm": 12, "xs": 12},
                            height=32,
                            alignment=ft.Alignment(-1, 0),
                        ),
                        ft.Container(
                            content=self.primer_toolbar,
                            col={"lg": 9, "md": 9, "sm": 12, "xs": 12},
                            alignment=ft.Alignment(1, 0),
                        ),
                    ],
                    run_spacing=0,
                ),
                ft.Container(
                    content=ft.Column(
                        [
                            self.primers_header_container,
                            self.primers_list,
                        ],
                        expand=True,
                        spacing=0,
                    ),
                    expand=True,
                    border=ft.Border.all(1, GUIColors.OUTLINE),
                    border_radius=5,
                    padding=0,
                ),
                self.error_banner,
                self.primer_info_panel,
            ],
            expand=True,
            spacing=5,
        )

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
        """Extract primer data from UI controls."""
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

    def sync_to_state(self, rebuild_if_needed: bool = False) -> bool:
        """Sync current UI controls back to the central state.

        Returns True if a UI rebuild is needed.
        """
        ui_primers = self._extract_primer_data_from_ui()
        primers = []
        for i, p in enumerate(ui_primers):
            prev_p = (
                self.input_data.primers[i]
                if i < len(self.input_data.primers)
                else {}
            )
            # Auto-activate when transitioning from empty to filled
            is_filled = bool(p["name"].strip() and p["seq"].strip())
            is_active = p["active"]
            checkbox = p.get("checkbox")
            if checkbox and checkbox.disabled and is_filled:
                is_active = True
                checkbox.value = True
                checkbox.disabled = False

            primers.append(
                {
                    "name": p["name"],
                    "seq": p["seq"],
                    "active": is_active,
                    "name_touched": prev_p.get("name_touched", False),
                    "seq_touched": prev_p.get("seq_touched", False),
                }
            )

        dup_indices = self._get_duplicate_indices_for_list(ui_primers)
        for p in ui_primers:
            container = p["container"]
            c_idx = container.data
            is_dup = c_idx in dup_indices
            new_color = GUIColors.DUPLICATE_BG if is_dup else None
            if container.bgcolor != new_color:
                container.bgcolor = new_color

        # Run background primer construction/validation
        new_validation_errors = validate_primers(primers)

        # Force active = False in state for empty fields
        for p in primers:
            is_empty = not p["name"].strip() or not p["seq"].strip()
            if is_empty:
                p["active"] = False

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
            self._update_row_highlights()
            self._update_header_checkbox_state()
            self._update_delete_button_disabled_state()

        return rebuild_if_needed

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        self.primers_list.update_list_ui()
        self._update_delete_button_disabled_state()

    def _update_delete_button_disabled_state(self) -> None:
        """Update disabled state of the delete button based on selection.

        Selection means that the primer is active.
        """
        has_selected = any(
            p.get("active", False) for p in self.input_data.primers
        )
        self.delete_selected_button.disabled = not has_selected
        if self.delete_selected_button.parent:
            self.delete_selected_button.update()

    def _on_toggle_all_primers(self, e: Any) -> None:
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
        if self.on_change_handler:
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
        if self.app_page:
            self.app_page.update()

    def _get_duplicate_indices_for_list(
        self, primers: list[dict[str, Any]]
    ) -> set[int]:
        """Find and return indices of duplicate primers by name/sequence."""
        names_count: dict[str, int] = {}
        seqs_count: dict[str, int] = {}
        for p in primers:
            n_lower = str(p.get("name", "")).strip().lower()
            s_lower = clean_sequence(str(p.get("seq", ""))).lower()
            if n_lower:
                names_count[n_lower] = names_count.get(n_lower, 0) + 1
            if s_lower:
                seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1

        dup_indices = set()
        for idx, p in enumerate(primers):
            n_lower = str(p.get("name", "")).strip().lower()
            s_lower = clean_sequence(str(p.get("seq", ""))).lower()
            if (n_lower and names_count.get(n_lower, 0) > 1) or (
                s_lower and seqs_count.get(s_lower, 0) > 1
            ):
                # If the item has a container reference, use container.data,
                # otherwise use idx
                c = p.get("container")
                c_idx = (
                    c.data if (c is not None and hasattr(c, "data")) else idx
                )
                dup_indices.add(c_idx)
        return dup_indices

    def _get_duplicate_indices(self) -> set[int]:
        """Find indices of primers with duplicate names or sequences."""
        return self._get_duplicate_indices_for_list(self.input_data.primers)

    def _update_row_highlights(self) -> None:
        """Update background colors of all row containers.

        Highlights rows based on selection and duplicates.
        """
        self.primers_list.update_row_highlights()

    def _update_primer_info_panel(self) -> None:
        """Update the primer information panel based on the focused primer."""

        def on_dismiss() -> None:
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

    async def _load_primers_click(self, e: ft.ControlEvent) -> None:
        """Open file picker to load primers from CSV/TSV file."""
        await self.file_manager.load_primers_click(e)

    async def _save_primers_click(self, e: ft.ControlEvent) -> None:
        """Save active primers to a CSV file."""
        await self.file_manager.save_primers_click(e)

    def _show_notification(self, message: str) -> None:
        """Show a notification message."""
        if not hasattr(self, "_notification_helper"):
            self._notification_helper = NotificationHelper(self.app_page)
        self._notification_helper.show_message(message)
