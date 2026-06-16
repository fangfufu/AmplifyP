# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""Input component for DNA primers list and details."""

from typing import Any

import flet as ft

from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import NotificationHelper, clean_sequence

from .primer_file_manager import PrimerFileManager
from .primer_header import PrimerHeader
from .primer_info_panel import PrimerInfoPanel
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

        self.focused_primer_index: int | None = None
        self.validation_errors: list[dict[str, str | None]] = []
        self._prev_header_checkbox_value: bool | None = None
        self._visible_rows_cache: list[PrimerRow] | None = None

        font_family = self.settings.get("font_family", "Roboto Mono")
        self.name_column_width = 150.0

        # Calculate initial name column ratio
        self.name_column_ratio = self.name_column_width / self.get_panel_width()

        # Primer List Component
        self.primers_list = PrimerList(self)

        # Primer Header Component
        self.primer_header = PrimerHeader(
            settings=self.settings,
            on_toggle_all=self._on_toggle_all_primers,
            on_divider_pan=self._on_primer_divider_pan,
            on_divider_pan_end=self._on_primer_divider_pan_end,
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
        )
        # Compatibility links
        self.clear_primers_button = self.primer_toolbar.clear_button
        self.save_primers_button = self.primer_toolbar.save_button
        self.load_primers_button = self.primer_toolbar.load_button

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

    def _handle_row_click(self, idx: int, name_edit: ft.TextField) -> None:
        """Handle clicking on the row container.

        Selects the row and focuses the name field.
        """
        import inspect

        self.focused_primer_index = idx
        self._update_row_highlights()
        self._update_primer_info_panel()
        res = name_edit.focus()
        if inspect.iscoroutine(res):
            if self.app_page:

                async def do_focus() -> None:
                    await res

                self.app_page.run_task(do_focus)

    def _on_add_primer_row(self, idx: int) -> None:
        """Add a new empty primer row immediately below the row at idx."""
        self.sync_to_state(rebuild_if_needed=False)
        self.input_data.primers.insert(
            idx + 1, {"name": "", "seq": "", "active": False}
        )
        self.update_ui()
        if self.on_change_handler:
            self.on_change_handler(None)

    def _move_primer(self, idx: int, direction: int) -> None:
        """Move primer at idx up (direction=-1) or down (direction=1)."""
        self.sync_to_state(rebuild_if_needed=False)
        primers = self.input_data.primers
        target_idx = idx + direction

        # Swap only if both indices are valid filled primers.
        if 0 <= idx < len(primers) and 0 <= target_idx < len(primers):
            primers[idx], primers[target_idx] = (
                primers[target_idx],
                primers[idx],
            )
            if self.focused_primer_index == idx:
                self.focused_primer_index = target_idx
            elif self.focused_primer_index == target_idx:
                self.focused_primer_index = idx

            self.update_ui()
            if self.on_change_handler:
                self.on_change_handler(None)

    def _delete_primer(self, idx: int) -> None:
        """Delete primer at idx."""
        self.sync_to_state(rebuild_if_needed=False)
        primers = self.input_data.primers
        if 0 <= idx < len(primers):
            primers.pop(idx)
            if self.focused_primer_index == idx:
                self.focused_primer_index = None
            elif (
                self.focused_primer_index is not None
                and self.focused_primer_index > idx
            ):
                self.focused_primer_index -= 1

            self.update_ui()
            if self.on_change_handler:
                self.on_change_handler(None)

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
            primers.append(
                {
                    "name": p["name"],
                    "seq": p["seq"],
                    "active": p["active"],
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

        # Force active = False in state for any primer with errors or
        # empty fields
        for idx, p in enumerate(primers):
            err = new_validation_errors[idx]
            is_empty = not p["name"].strip() or not p["seq"].strip()
            if err.get("name") or err.get("seq") or is_empty:
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

        return rebuild_if_needed

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        self.primers_list.update_list_ui()

    def get_panel_width(self) -> float:
        """Get the current width of the primer input panel."""
        page_width = (
            self.app_page.width
            if (self.app_page and self.app_page.width)
            else 800.0
        )
        parent_view = getattr(self.on_change_handler, "__self__", None)
        right_fraction = (
            getattr(parent_view, "right_fraction", 0.5) if parent_view else 0.5
        )
        return page_width * right_fraction

    def _cache_visible_rows_if_needed(self) -> None:
        """Cache the list of visible rows based on viewport scroll state."""
        if self._visible_rows_cache is not None:
            return

        self._visible_rows_cache = []
        scroll_y = self.primers_list.scroll_pixels
        viewport_h = self.primers_list.viewport_dimension
        if self.app_page and self.app_page.height:
            viewport_h = max(viewport_h, float(self.app_page.height))
        current_y = 0.0
        for row in self.primers_list.controls:
            if isinstance(row, PrimerRow):
                row_h = (
                    30.0
                    if not (row.name_field.error or row.seq_field.error)
                    else 50.0
                )
                row_top = current_y
                row_bottom = current_y + row_h

                # Check if row is visible (plus 60px buffer above/below)
                if (row_bottom >= scroll_y - 60.0) and (
                    row_top <= scroll_y + viewport_h + 60.0
                ):
                    self._visible_rows_cache.append(row)

                current_y += row_h

    def _on_primer_divider_pan(self, e: ft.DragUpdateEvent) -> None:
        """Handle dragging the vertical divider between Name and Sequence."""
        delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
        panel_width = self.get_panel_width()

        # Subtract space for other UI components in the row:
        # - Checkbox container (55)
        # - Two dividers (4 + 4 = 8)
        # - Control container when focused (108)
        # - Minimum space to display "Sequence" text and edit cursor (100)
        max_name_width = max(80.0, panel_width - 271.0)

        target_width = max(
            80.0, min(max_name_width, self.name_column_width + delta_x)
        )

        self.name_column_width = target_width
        self.name_column_ratio = self.name_column_width / panel_width

        # Update the width of the Name header control
        self.primers_header.controls[2].width = self.name_column_width
        self.primer_header.update()

        self._cache_visible_rows_if_needed()
        visible_rows = self._visible_rows_cache or []

        # Update and render only the name fields of the visible rows directly
        for row in visible_rows:
            row.name_field.width = self.name_column_width
            row.name_field.update()

    def _on_primer_divider_pan_end(self, e: ft.DragEndEvent) -> None:
        """Handle finishing the drag of the vertical divider."""
        # Clear the visible rows cache
        self._visible_rows_cache = None

        # Ensure the final exact width is applied to header and all rows in sync
        self.primers_header.controls[2].width = self.name_column_width
        self.primer_header.update()

        for row in self.primers_list.controls:
            if isinstance(row, PrimerRow):
                row.name_field.width = self.name_column_width
        self.primers_list.update()

    def adjust_name_column_width(
        self, new_panel_width: float, during_drag: bool = False
    ) -> None:
        """Scale name column width proportionally based on panel width."""
        max_name_width = max(80.0, new_panel_width - 271.0)
        target_width = max(
            80.0,
            min(max_name_width, new_panel_width * self.name_column_ratio),
        )

        self.name_column_width = target_width

        # Apply the new width to header
        if (
            hasattr(self, "primers_header")
            and self.primers_header
            and len(self.primers_header.controls) > 2
        ):
            self.primers_header.controls[2].width = self.name_column_width
            self.primer_header.update()

        if during_drag:
            self._cache_visible_rows_if_needed()
            visible_rows = self._visible_rows_cache or []

            # Update/render name fields of visible rows directly
            for row in visible_rows:
                row.name_field.width = self.name_column_width
                row.name_field.update()
        else:
            # Clear cache and update all rows
            self._visible_rows_cache = None
            for row in self.primers_list.controls:
                if isinstance(row, PrimerRow):
                    row.name_field.width = self.name_column_width
            self.primers_list.update()

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
