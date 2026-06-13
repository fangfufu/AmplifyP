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

        font_family = self.settings.get("font_family", "Roboto Mono")
        self.name_column_width = 150.0

        # Primer List Component
        self.primers_list = PrimerList(self)

        # Primer Header Component
        self.primer_header = PrimerHeader(
            settings=self.settings,
            on_toggle_all=self._on_toggle_all_primers,
            on_divider_pan=self._on_primer_divider_pan,
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

    def _move_primer(self, idx: int, direction: int) -> None:
        """Move primer at idx up (direction=-1) or down (direction=1)."""
        parent_view = getattr(self.on_change_handler, "__self__", None)
        if parent_view is not None:
            if getattr(parent_view, "_focus_debouncer", None) is not None:
                parent_view._focus_debouncer.cancel()
            parent_view._skip_blur_timer = True

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
        parent_view = getattr(self.on_change_handler, "__self__", None)
        if parent_view is not None:
            if getattr(parent_view, "_focus_debouncer", None) is not None:
                parent_view._focus_debouncer.cancel()
            parent_view._skip_blur_timer = True

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

    def sync_to_state(self, rebuild_if_needed: bool = True) -> bool:
        """Sync current UI controls back to the central state.

        Returns True if a UI rebuild is needed.
        """
        ui_primers = self._extract_primer_data_from_ui()
        filtered_ui_primers, should_rebuild = self._apply_activation_rules(
            ui_primers
        )

        dup_indices = self._get_duplicate_indices_for_list(filtered_ui_primers)
        primers = []
        for p in filtered_ui_primers:
            container = p["container"]
            c_idx = container.data
            is_dup = c_idx in dup_indices
            new_color = GUIColors.DUPLICATE_BG if is_dup else None
            if container.bgcolor != new_color:
                container.bgcolor = new_color

            primers.append(
                {
                    "name": p["name"],
                    "seq": p["seq"],
                    "active": p["active"],
                }
            )

        # Run background primer construction/validation
        new_validation_errors = self.validate_primers(primers)

        # Rebuild if the number of controls in the UI doesn't match the
        # expected rows (filtered + 1 trailing empty row)
        expected_rows = len(filtered_ui_primers) + 1
        if len(self.primers_list.controls) != expected_rows:
            should_rebuild = True

        self.input_data.primers = primers
        self.validation_errors = new_validation_errors
        if should_rebuild and rebuild_if_needed:
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

        return should_rebuild

    def validate_primers(
        self, primers: list[dict[str, Any]]
    ) -> list[dict[str, str | None]]:
        """Validate a list of primers, detecting format and duplicate errors."""
        names_count: dict[str, int] = {}
        seqs_count: dict[str, int] = {}
        for p in primers:
            n_lower = str(p.get("name", "")).strip().lower()
            s_lower = clean_sequence(str(p.get("seq", ""))).lower()
            if n_lower:
                names_count[n_lower] = names_count.get(n_lower, 0) + 1
            if s_lower:
                seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1

        errors = []
        for p in primers:
            name_val = p.get("name", "")
            seq_val = p.get("seq", "")
            seq_err = PrimerRow.validate(name_val, seq_val)
            name_err = None

            n_lower = str(name_val).strip().lower()
            s_lower = clean_sequence(str(seq_val)).lower()

            if not seq_err and s_lower and seqs_count.get(s_lower, 0) > 1:
                seq_err = "Duplicate primer sequence"
            if n_lower and names_count.get(n_lower, 0) > 1:
                name_err = "Duplicate primer name"

            errors.append({"name": name_err, "seq": seq_err})
        return errors

    def _apply_activation_rules(
        self,
        ui_primers: list[dict[str, Any]],
    ) -> tuple[list[dict[str, Any]], bool]:
        """Apply auto-activation, auto-inactivation and deletion rules."""
        from amplifyp.dna import Primer

        prev_primers = self.input_data.primers
        filtered_primers = []
        should_rebuild = False

        # Precompute counts for duplicate checks in activation logic
        names_count: dict[str, int] = {}
        seqs_count: dict[str, int] = {}
        for p in ui_primers:
            n_lower = str(p.get("name", "")).strip().lower()
            s_lower = clean_sequence(str(p.get("seq", ""))).lower()
            if n_lower:
                names_count[n_lower] = names_count.get(n_lower, 0) + 1
            if s_lower:
                seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1

        for i, p in enumerate(ui_primers):
            name_val = p["name"]
            seq_val = p["seq"]
            is_active = p["active"]
            checkbox = p.get("checkbox")

            # Deletion rule: if both are empty, delete it
            if not name_val and not seq_val:
                continue

            prev_primer = prev_primers[i] if i < len(prev_primers) else None
            prev_name = (
                str(prev_primer.get("name", "")).strip() if prev_primer else ""
            )
            prev_seq = (
                clean_sequence(str(prev_primer.get("seq", "")))
                if prev_primer
                else ""
            )

            prev_seq_invalid = False
            if prev_seq:
                try:
                    Primer(sequence=prev_seq, name=prev_name)
                except ValueError:
                    prev_seq_invalid = True

            # Check if current primer would be duplicate/invalid
            n_lower = name_val.strip().lower()
            s_lower = clean_sequence(seq_val).lower()
            is_dup = (n_lower and names_count.get(n_lower, 0) > 1) or (
                s_lower and seqs_count.get(s_lower, 0) > 1
            )

            # Auto-activation rule
            if (not prev_name or not prev_seq or prev_seq_invalid) and (
                name_val and seq_val
            ):
                is_valid = True
                try:
                    Primer(sequence=seq_val, name=name_val)
                except ValueError:
                    is_valid = False
                if is_valid and not is_dup:
                    is_active = True
                    if checkbox is not None:
                        checkbox.value = True
                    should_rebuild = True

            # Auto-inactivation rule if fields are cleared or invalid
            is_valid = True
            if name_val or seq_val:
                try:
                    Primer(sequence=seq_val, name=name_val)
                except ValueError:
                    is_valid = False
                if is_dup:
                    is_valid = False

            if not name_val or not seq_val or not is_valid:
                if is_active:
                    is_active = False
                    if checkbox is not None:
                        checkbox.value = False
                    should_rebuild = True

            p["active"] = is_active
            filtered_primers.append(p)

        return filtered_primers, should_rebuild

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        self.primers_list.update_list_ui()

    def _on_primer_divider_pan(self, e: ft.DragUpdateEvent) -> None:
        """Handle dragging the vertical divider between Name and Sequence."""
        delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
        self.name_column_width = max(
            80.0, min(300.0, self.name_column_width + delta_x)
        )
        # Update the width of the Name header control
        self.primers_header.controls[2].width = self.name_column_width
        # Update the width of all Name TextFields in the list controls
        for row in self.primers_list.controls:
            if isinstance(row, PrimerRow):
                row.name_field.width = self.name_column_width
        self.app_page.update()

    def _on_toggle_all_primers(self, e: Any) -> None:
        """Toggle all primers active/inactive based on tri-state checkbox."""
        primers = self.input_data.primers
        if not primers:
            return
        non_empty = [
            p
            for p in primers
            if str(p.get("name", "")).strip()
            or clean_sequence(str(p.get("seq", ""))).strip()
        ]
        if not non_empty:
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

        for p in non_empty:
            p["active"] = target_active

        # Update checkbox values in-place on existing controls to avoid
        # rebuilding the list
        for row in self.primers_list.controls:
            if isinstance(row, PrimerRow) and row.data is not None:
                idx = row.data
                if idx < len(primers):
                    p = primers[idx]
                    is_non_empty = (
                        str(p.get("name", "")).strip()
                        or clean_sequence(str(p.get("seq", ""))).strip()
                    )
                    if is_non_empty and not row.checkbox.disabled:
                        row.checkbox.value = target_active

        self._prev_header_checkbox_value = cb_value
        self._update_header_checkbox_state()
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
