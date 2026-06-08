# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Input view for DNA template and primers."""

import csv
import io
import threading
from typing import Any

import flet as ft

from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import clean_sequence, format_sequence


class TemplateInput(ft.Container):  # type: ignore[misc]
    """Input component for DNA template sequence."""

    def __init__(
        self,
        settings: GUISettings,
        input_data: GUIInput,
        on_change_handler: Any,
        handle_field_focus: Any,
        handle_field_blur: Any,
        handle_field_submit: Any,
        clear_template_callback: Any,
    ) -> None:
        """Initialize the TemplateInput component."""
        super().__init__(expand=5)
        self.settings = settings
        self.input_data = input_data
        font_family = self.settings.get("font_family", "Roboto Mono")

        self.template_sequence = ft.TextField(
            multiline=True,
            max_lines=None,
            border=ft.InputBorder.NONE,
            text_align=ft.TextAlign.LEFT,
            hint_text="Enter DNA sequence here...",
            content_padding=0,
            on_change=on_change_handler,
            on_focus=handle_field_focus,
            on_blur=handle_field_blur,
            on_submit=handle_field_submit,
            text_style=ft.TextStyle(font_family=font_family),
        )
        self.template_circular = ft.Checkbox(
            label="Circular Template",
            value=False,
            on_change=on_change_handler,
        )
        self.circular_container = ft.Container(
            content=self.template_circular,
            border=ft.Border.all(1, GUIColors.OUTLINE),
            border_radius=5,
            padding=ft.Padding(10, 0, 10, 0),
            height=32,
            alignment=ft.Alignment(0, 0),
        )

        self.clear_template_button = ft.OutlinedButton(
            "Clear",
            icon=ft.Icons.DELETE_OUTLINE,
            tooltip="Clear Template",
            on_click=clear_template_callback,
            height=32,
        )

        self.content = ft.Column(
            [
                ft.Row(
                    [
                        ft.Text("Template Sequence", weight=ft.FontWeight.BOLD),
                        ft.Row(
                            [
                                self.circular_container,
                                self.clear_template_button,
                            ],
                            spacing=10,
                        ),
                    ],
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    height=40,
                ),
                ft.Container(
                    content=ft.ListView(
                        [self.template_sequence],
                        expand=True,
                        scroll=ft.ScrollMode.ALWAYS,
                    ),
                    expand=True,
                    border=ft.Border.all(1, GUIColors.OUTLINE),
                    border_radius=5,
                    padding=0,
                ),
            ],
            expand=True,
            spacing=5,
        )

    def sync_to_state(self) -> None:
        """Sync template text field to the central state."""
        self.template_sequence.value = self.template_sequence.value or ""
        self.input_data.template = clean_sequence(
            str(self.template_sequence.value)
        )
        self.input_data.template_circular = bool(self.template_circular.value)

    def update_ui(self) -> None:
        """Update template UI elements to match central state."""
        font_family = self.settings.get("font_family", "Roboto Mono")
        self.template_sequence.text_style = ft.TextStyle(
            font_family=font_family
        )
        self.template_sequence.value = self.input_data.template
        self.template_circular.value = self.input_data.template_circular


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
        self.validation_errors: list[str | None] = []

        font_family = self.settings.get("font_family", "Roboto Mono")
        self.name_column_width = 150.0
        self.primers_list = ft.ListView(
            expand=True, spacing=0, padding=0, scroll=ft.ScrollMode.ALWAYS
        )
        self.all_primers_checkbox = ft.Checkbox(
            value=None,
            tristate=True,
            on_change=self._on_toggle_all_primers,
        )
        self.primers_header = ft.Row(
            [
                ft.Container(
                    content=self.all_primers_checkbox,
                    width=55,
                    alignment=ft.Alignment(0, 0),
                ),
                ft.Container(
                    width=4,
                    bgcolor=GUIColors.DIVIDER_GREY,
                    margin=0,
                    height=36,
                ),
                ft.Container(
                    content=ft.Text(
                        "Name",
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_small", 12),
                    ),
                    width=self.name_column_width,
                    padding=ft.Padding(5, 0, 0, 0),
                    alignment=ft.Alignment(-1, 0),
                ),
                ft.GestureDetector(
                    on_pan_update=self._on_primer_divider_pan,
                    content=ft.Container(
                        width=4,
                        bgcolor=GUIColors.DIVIDER_GREY,
                        margin=0,
                        height=36,
                    ),
                    mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
                ),
                ft.Container(
                    content=ft.Text(
                        "Sequence",
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_small", 12),
                    ),
                    expand=True,
                    padding=ft.Padding(5, 0, 0, 0),
                    alignment=ft.Alignment(-1, 0),
                ),
            ],
            alignment=ft.MainAxisAlignment.START,
            height=36,
            spacing=0,
        )
        self.primers_header_container = ft.Container(
            content=self.primers_header,
            padding=0,
            border=ft.Border(bottom=ft.BorderSide(4, GUIColors.DIVIDER_GREY)),
            height=40,
        )

        self.clear_primers_button = ft.OutlinedButton(
            "Clear",
            icon=ft.Icons.DELETE_OUTLINE,
            tooltip="Clear All Primers",
            on_click=self.clear_primers_callback,
            height=32,
        )

        self.save_primers_button = ft.FilledTonalButton(
            "Save Primers",
            icon=ft.Icons.FILE_DOWNLOAD,
            on_click=self._save_primers_click,
            tooltip="Save active primers to CSV",
            height=32,
        )
        self.load_primers_button = ft.FilledTonalButton(
            "Load Primers",
            icon=ft.Icons.FILE_OPEN,
            on_click=self._load_primers_click,
            tooltip="Load primers from CSV/TSV",
            height=32,
        )

        # Primer Info Panel UI Components
        self.info_header = ft.Container(
            content=ft.Text(
                "Primer: -",
                weight=ft.FontWeight.BOLD,
                size=self.settings.get("font_size_default", 14),
                color=GUIColors.DIAGRAM_BLACK,
                selectable=True,
            ),
            bgcolor=GUIColors.INFO_HEADER_BG,
            padding=ft.Padding(10, 5, 10, 5),
            alignment=ft.Alignment(-1, 0),
        )

        self.info_seq_text = ft.Text(
            "",
            font_family=font_family,
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_tm_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_pairs_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_redundancy_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_dimer_text = ft.Text(
            "",
            color=GUIColors.ERROR_RED,
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )

        self.primer_info_panel = ft.Container(
            content=ft.Column(
                [
                    self.info_header,
                    ft.Container(
                        content=ft.Column(
                            [
                                self.info_seq_text,
                                self.info_tm_text,
                                self.info_pairs_text,
                                self.info_redundancy_text,
                                self.info_dimer_text,
                            ],
                            spacing=3,
                        ),
                        padding=ft.Padding(10, 5, 10, 10),
                    ),
                ],
                spacing=0,
            ),
            border=ft.Border.all(1, GUIColors.OUTLINE),
            border_radius=5,
            visible=False,
        )

        self.content = ft.Column(
            [
                ft.Row(
                    [
                        ft.Text("Primers", weight=ft.FontWeight.BOLD),
                        ft.Row(
                            [
                                self.save_primers_button,
                                self.load_primers_button,
                                self.clear_primers_button,
                            ],
                            spacing=8,
                        ),
                    ],
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    height=40,
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

    def _handle_row_click(self, idx: int, name_edit: ft.TextField) -> None:
        """Handle clicking on the row container.

        Selects the row and focuses the name field.
        """
        self.focused_primer_index = idx
        self._update_row_highlights()
        self._update_primer_info_panel()
        name_edit.focus()

    def _move_primer(self, idx: int, direction: int) -> None:
        """Move primer at idx up (direction=-1) or down (direction=1)."""
        parent_view = getattr(self.on_change_handler, "__self__", None)
        if parent_view is not None:
            if getattr(parent_view, "_focus_timer", None) is not None:
                parent_view._focus_timer.cancel()
                parent_view._focus_timer = None
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

    def _extract_primer_data_from_ui(self) -> list[dict[str, Any]]:
        """Extract primer data from UI controls."""
        ui_primers = []
        for container in self.primers_list.controls:
            if not isinstance(container, ft.Container):
                continue
            row = container.content
            if not isinstance(row, ft.Row) or len(row.controls) < 5:
                continue

            checkbox_control = row.controls[0]
            checkbox = (
                checkbox_control.content
                if isinstance(checkbox_control, ft.Container)
                else checkbox_control
            )
            name_tf = row.controls[2]
            seq_ctrl = row.controls[4]
            seq_tf = (
                seq_ctrl.controls[0]
                if isinstance(seq_ctrl, ft.Stack) and seq_ctrl.controls
                else seq_ctrl
            )

            name_val = str(name_tf.value or "").strip()
            seq_val = clean_sequence(str(seq_tf.value or ""))
            is_active = bool(checkbox.value)

            ui_primers.append(
                {
                    "name": name_val,
                    "seq": seq_val,
                    "active": is_active,
                    "container": container,
                    "checkbox": checkbox,
                }
            )
        return ui_primers

    def _apply_activation_rules(
        self, ui_primers: list[dict[str, Any]]
    ) -> tuple[list[dict[str, Any]], bool]:
        """Apply auto-activation, auto-inactivation and deletion rules.

        Returns (filtered_ui_primers, rebuild_needed).
        """
        filtered_primers = []
        should_rebuild = False

        for i, p in enumerate(ui_primers):
            name_val = p["name"]
            seq_val = p["seq"]
            is_active = p["active"]
            checkbox = p["checkbox"]

            # Deletion rule: if both are empty, delete it
            if not name_val and not seq_val:
                should_rebuild = True
                continue

            prev_primer = (
                self.input_data.primers[i]
                if i < len(self.input_data.primers)
                else None
            )
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
                    from amplifyp.dna import Primer

                    Primer(sequence=prev_seq, name=prev_name)
                except ValueError:
                    prev_seq_invalid = True

            # Auto-activation rule: if it was previously incomplete or invalid
            # and both fields are now filled
            if (not prev_name or not prev_seq or prev_seq_invalid) and (
                name_val and seq_val
            ):
                # Verify that the sequence is valid before auto-activating
                is_valid = True
                try:
                    from amplifyp.dna import Primer

                    Primer(sequence=seq_val, name=name_val)
                except ValueError:
                    is_valid = False
                if is_valid:
                    is_active = True
                    checkbox.value = True
                    should_rebuild = True

            # Auto-inactivation rule: if either is empty, set active to False
            if not name_val or not seq_val:
                if is_active:
                    is_active = False
                    checkbox.value = False
                    should_rebuild = True

            p["active"] = is_active
            filtered_primers.append(p)

        return filtered_primers, should_rebuild

    def _detect_and_mark_duplicates(
        self, ui_primers: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        """Detect duplicates and mark container colors.

        Returns primers_list for state.
        """
        names_count: dict[str, int] = {}
        seqs_count: dict[str, int] = {}
        for p in ui_primers:
            n_lower = p["name"].lower()
            s_lower = p["seq"].lower()
            if n_lower:
                names_count[n_lower] = names_count.get(n_lower, 0) + 1
            if s_lower:
                seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1

        primers = []
        for p in ui_primers:
            container = p["container"]
            n_lower = p["name"].lower()
            s_lower = p["seq"].lower()

            is_dup = (n_lower and names_count.get(n_lower, 0) > 1) or (
                s_lower and seqs_count.get(s_lower, 0) > 1
            )

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
        return primers

    def sync_to_state(self, rebuild_if_needed: bool = True) -> bool:
        """Sync current UI controls back to the central state.

        Returns True if a UI rebuild is needed.
        """
        ui_primers = self._extract_primer_data_from_ui()
        filtered_ui_primers, should_rebuild = self._apply_activation_rules(
            ui_primers
        )
        primers = self._detect_and_mark_duplicates(filtered_ui_primers)

        # Run background primer construction/validation
        new_validation_errors = []
        for p in primers:
            name_val = p["name"]
            seq_val = p["seq"]
            error_message = None
            if seq_val:
                try:
                    from amplifyp.dna import Primer

                    Primer(sequence=seq_val, name=name_val)
                except ValueError as ex:
                    error_message = str(ex)
            new_validation_errors.append(error_message)

        if getattr(self, "validation_errors", []) != new_validation_errors:
            should_rebuild = True

        # If the number of non-empty primers changed, we rebuild
        if len(primers) != len(
            [
                p
                for p in self.input_data.primers
                if p.get("name", "").strip() or p.get("seq", "").strip()
            ]
        ):
            should_rebuild = True

        self.input_data.primers = primers
        self.validation_errors = new_validation_errors
        self._update_header_checkbox_state()

        if should_rebuild and rebuild_if_needed:
            self.update_ui()

        return should_rebuild

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        font_family = self.settings.get("font_family", "Roboto Mono")
        self.primers_list.controls.clear()

        # Filter out all empty primers from state first
        clean_primers = [
            p
            for p in self.input_data.primers
            if str(p.get("name", "")).strip()
            or clean_sequence(str(p.get("seq", ""))).strip()
        ]

        # Only append a trailing empty row if there are no invalid primers
        any_invalid = False
        for p in clean_primers:
            name_val = p.get("name", "")
            seq_val = p.get("seq", "")
            if seq_val:
                try:
                    from amplifyp.dna import Primer

                    Primer(sequence=seq_val, name=name_val)
                except ValueError:
                    any_invalid = True
                    break

        if not any_invalid:
            clean_primers.append({"name": "", "seq": "", "active": False})
        self.input_data.primers = clean_primers

        dup_indices = self._get_duplicate_indices()

        self.validation_errors = []
        for idx, p in enumerate(self.input_data.primers):
            name_val = p["name"]
            seq_val = p["seq"]
            is_active = p.get("active", True)

            error_message = None
            if seq_val:
                try:
                    from amplifyp.dna import Primer

                    Primer(sequence=seq_val, name=name_val)
                except ValueError as ex:
                    error_message = str(ex)
            self.validation_errors.append(error_message)

            is_dup = idx in dup_indices

            checkbox = ft.Checkbox(
                value=is_active,
                on_change=self.on_change_handler,
            )
            checkbox_container = ft.Container(
                content=checkbox,
                width=55,
                height=30,
                alignment=ft.Alignment(0, 0),
            )
            name_edit = ft.TextField(
                value=name_val,
                hint_text="New Primer Name",
                dense=True,
                content_padding=ft.Padding(5, 0, 0, 0),
                height=30,
                width=self.name_column_width,
                border=ft.InputBorder.NONE,
                data=idx,
                on_focus=self.handle_field_focus,
                on_blur=self.handle_field_blur,
                on_submit=self.handle_field_submit,
            )
            seq_edit = ft.TextField(
                value=seq_val,
                hint_text="New Primer Sequence",
                dense=True,
                content_padding=ft.Padding(
                    5,
                    0,
                    60
                    if idx == self.focused_primer_index
                    and idx < len(self.input_data.primers) - 1
                    else 0,
                    0,
                ),
                height=30 if not error_message else None,
                border=ft.InputBorder.NONE,
                text_style=ft.TextStyle(font_family=font_family),
                data=idx,
                on_focus=self.handle_field_focus,
                on_blur=self.handle_field_blur,
                on_submit=self.handle_field_submit,
            )
            if error_message:
                seq_edit.error = error_message
            divider = ft.GestureDetector(
                on_pan_update=self._on_primer_divider_pan,
                content=ft.Container(
                    width=4,
                    bgcolor=GUIColors.DIVIDER_GREY,
                    margin=0,
                    height=30,
                ),
                mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
            )

            active_divider = ft.Container(
                width=4,
                bgcolor=GUIColors.DIVIDER_GREY,
                margin=0,
                height=30,
            )

            seq_stack = ft.Stack(
                [
                    seq_edit,
                ],
                expand=True,
            )

            if idx < len(self.input_data.primers) - 1:
                up_button = ft.IconButton(
                    icon=ft.Icons.ARROW_UPWARD,
                    icon_size=16,
                    width=24,
                    height=24,
                    padding=0,
                    tooltip="Move Up",
                    disabled=(idx == 0),
                    on_click=lambda e, idx=idx: self._move_primer(idx, -1),
                )
                down_button = ft.IconButton(
                    icon=ft.Icons.ARROW_DOWNWARD,
                    icon_size=16,
                    width=24,
                    height=24,
                    padding=0,
                    tooltip="Move Down",
                    disabled=(idx == len(self.input_data.primers) - 2),
                    on_click=lambda e, idx=idx: self._move_primer(idx, 1),
                )
                reorder_controls = ft.Row(
                    [up_button, down_button],
                    spacing=2,
                    alignment=ft.MainAxisAlignment.CENTER,
                )
                reorder_container = ft.Container(
                    content=reorder_controls,
                    right=10,
                    top=3,
                    visible=(idx == self.focused_primer_index),
                )
                seq_stack.controls.append(reorder_container)

            row = ft.Row(
                [
                    checkbox_container,
                    active_divider,
                    name_edit,
                    divider,
                    seq_stack,
                ],
                spacing=0,
                vertical_alignment=ft.CrossAxisAlignment.START,
            )

            row_container = ft.Container(
                content=row,
                bgcolor=GUIColors.DUPLICATE_BG if is_dup else None,
                padding=0,
                height=30 if not error_message else None,
                data=idx,
                on_click=lambda e, idx=idx, name_edit=name_edit: (
                    self._handle_row_click(idx, name_edit)
                ),
            )
            self.primers_list.controls.append(row_container)

        self._update_row_highlights()
        self._update_primer_info_panel()
        self._update_header_checkbox_state()

    def _on_primer_divider_pan(self, e: ft.DragUpdateEvent) -> None:
        """Handle dragging the vertical divider between Name and Sequence."""
        delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
        self.name_column_width = max(
            80.0, min(300.0, self.name_column_width + delta_x)
        )
        # Update the width of the Name header control
        self.primers_header.controls[2].width = self.name_column_width
        # Update the width of all Name TextFields in the list controls
        for container in self.primers_list.controls:
            if isinstance(container, ft.Container):
                row = container.content
                if isinstance(row, ft.Row) and len(row.controls) >= 3:
                    name_tf = row.controls[2]
                    if isinstance(name_tf, ft.TextField):
                        name_tf.width = self.name_column_width
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
        all_active = all(p.get("active", True) for p in non_empty)
        target_active = not all_active
        for p in non_empty:
            p["active"] = target_active
        self.update_ui()

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

    def _get_duplicate_indices(self) -> set[int]:
        """Find indices of primers with duplicate names or sequences."""
        names_count: dict[str, int] = {}
        seqs_count: dict[str, int] = {}
        for p in self.input_data.primers:
            n_lower = str(p.get("name", "")).strip().lower()
            s_lower = clean_sequence(str(p.get("seq", ""))).lower()
            if n_lower:
                names_count[n_lower] = names_count.get(n_lower, 0) + 1
            if s_lower:
                seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1

        dup_indices = set()
        for idx, p in enumerate(self.input_data.primers):
            n_lower = str(p.get("name", "")).strip().lower()
            s_lower = clean_sequence(str(p.get("seq", ""))).lower()
            if (n_lower and names_count.get(n_lower, 0) > 1) or (
                s_lower and seqs_count.get(s_lower, 0) > 1
            ):
                dup_indices.add(idx)
        return dup_indices

    def _update_row_highlights(self) -> None:
        """Update background colors of all row containers.

        Highlights rows based on selection and duplicates.
        """
        dup_indices = self._get_duplicate_indices()

        for container in self.primers_list.controls:
            if (
                isinstance(container, ft.Container)
                and container.data is not None
            ):
                c_idx = container.data
                is_dup = c_idx in dup_indices

                if c_idx == self.focused_primer_index:
                    container.bgcolor = GUIColors.SELECTED_ROW_BG
                elif is_dup:
                    container.bgcolor = GUIColors.DUPLICATE_BG
                else:
                    container.bgcolor = None

                # Update reorder container visibility and text padding
                # dynamically.
                row = container.content
                if isinstance(row, ft.Row) and len(row.controls) >= 5:
                    seq_stack = row.controls[4]
                    if (
                        isinstance(seq_stack, ft.Stack)
                        and len(seq_stack.controls) > 1
                    ):
                        reorder_container = seq_stack.controls[1]
                        seq_edit = seq_stack.controls[0]
                        is_focused = c_idx == self.focused_primer_index
                        reorder_container.visible = is_focused
                        if isinstance(seq_edit, ft.TextField):
                            seq_edit.content_padding = ft.Padding(
                                5, 0, 60 if is_focused else 0, 0
                            )

    def _update_primer_info_panel(self) -> None:
        """Update the primer information panel based on the focused primer."""
        if self.focused_primer_index is None:
            self.primer_info_panel.visible = False
            self._update_row_highlights()
            self.app_page.update()
            return

        try:
            primers = self.input_data.primers
            if (
                self.focused_primer_index < 0
                or self.focused_primer_index >= len(primers)
            ):
                self.primer_info_panel.visible = False
                self.app_page.update()
                return

            p_data = primers[self.focused_primer_index]
            name_val = p_data.get("name", "").strip()
            seq_val = clean_sequence(p_data.get("seq", ""))

            if not seq_val:
                self.primer_info_panel.visible = False
                self.app_page.update()
                return

            from amplifyp.dimer import PrimerDimerGenerator
            from amplifyp.dna import Primer

            primer_obj = Primer(sequence=seq_val, name=name_val)

            # Header
            header_text = (
                f"Primer: {name_val}" if name_val else f"Primer: {seq_val}"
            )
            self.info_header.content.value = header_text

            # Sequence details
            self.info_seq_text.value = (
                f"{len(primer_obj)} bp:   {primer_obj.seq}"
            )

            # Melting temperature
            tm = self.settings.calculate_primer_tm(primer_obj)
            self.info_tm_text.value = f"Tm = {tm:.2f}°C"

            # Base counts / percentage
            self.info_pairs_text.value = (
                f"{primer_obj.count_at()} AT Pairs, "
                f"{primer_obj.count_cg()} GC Pairs, "
                f"{primer_obj.ratio_at() * 100:.1f}% AT"
            )

            # Redundancy
            if primer_obj.redundant_base_count == 0:
                self.info_redundancy_text.value = "No redundant bases."
            else:
                self.info_redundancy_text.value = (
                    f"{primer_obj.redundant_base_count} redundant bases "
                    f"(redundancy fold = {primer_obj.redundancy_fold})."
                )

            # Dimer potential check
            pd_settings = self.settings.get_primer_dimer_settings()
            generator = PrimerDimerGenerator(settings=pd_settings)

            # Check self-dimer and cross-dimer potential
            # against all active primers
            active_primers = self.input_data.get_active_primers()
            max_dimer = None

            # Always check self-dimer
            res_self = generator.generate_primer_dimer(primer_obj, primer_obj)
            if (
                res_self.overlap > pd_settings.min_overlap
                and res_self.quality > pd_settings.threshold
            ):
                max_dimer = res_self

            for ap in active_primers:
                ap_name = ap.get("name", "").strip()
                ap_seq = clean_sequence(ap.get("seq", ""))
                if not ap_seq:
                    continue
                try:
                    ap_obj = Primer(sequence=ap_seq, name=ap_name)
                    res = generator.generate_primer_dimer(primer_obj, ap_obj)
                    if (
                        res.overlap > pd_settings.min_overlap
                        and res.quality > pd_settings.threshold
                    ):
                        if max_dimer is None or res.quality > max_dimer.quality:
                            max_dimer = res
                except ValueError:
                    continue

            if max_dimer is not None:
                self.info_dimer_text.value = (
                    "Potential Primer Dimer with quality = "
                    f"{max_dimer.quality:.1f} "
                    f"and overlap = {max_dimer.overlap}"
                )
                self.info_dimer_text.visible = True
            else:
                self.info_dimer_text.value = ""
                self.info_dimer_text.visible = False

            self.primer_info_panel.visible = True

        except Exception:
            self.primer_info_panel.visible = False

        self.app_page.update()

    async def _load_primers_click(self, e: ft.ControlEvent) -> None:
        """Open file picker to load primers from CSV/TSV file."""
        try:
            files = await ft.FilePicker().pick_files(
                dialog_title="Load Primers",
                allowed_extensions=["csv", "tsv", "txt"],
                file_type=ft.FilePickerFileType.CUSTOM,
                with_data=True,
            )
            if not files:
                return

            file = files[0]
            if file.bytes is not None:
                content = file.bytes.decode("utf-8")
            else:
                if not file.path:
                    self._show_snackbar("Error: Could not read file content.")
                    return
                with open(file.path, encoding="utf-8") as f:
                    content = f.read()

            loaded_count = 0
            for line in content.strip().splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                if "\t" in line:
                    delimiter = "\t"
                else:
                    delimiter = ","

                parts = line.split(delimiter, 1)
                if len(parts) != 2:
                    continue

                name = parts[0].strip()
                seq = parts[1].strip()

                if not name or not seq:
                    continue

                try:
                    from amplifyp.dna import Primer

                    Primer(sequence=seq, name=name)
                except ValueError:
                    continue

                self.input_data.primers.append(
                    {
                        "name": name,
                        "seq": seq,
                        "active": True,
                    }
                )
                loaded_count += 1

            if loaded_count > 0:
                self.update_ui()
                if self.on_change_handler:
                    self.on_change_handler(None)
                self._show_snackbar(f"Loaded {loaded_count} primer(s).")
            else:
                self._show_snackbar("No valid primers found in file.")

        except Exception as ex:
            self._show_snackbar(f"Error loading file: {ex}")

    async def _save_primers_click(self, e: ft.ControlEvent) -> None:
        """Save active primers to a CSV file."""
        active_primers = [
            p for p in self.input_data.primers if p.get("active", True)
        ]
        if not active_primers:
            self._show_snackbar("No active primers to save.")
            return

        output = io.StringIO()
        writer = csv.writer(output)
        for p in active_primers:
            writer.writerow([p.get("name", ""), p.get("seq", "")])

        csv_content = output.getvalue()
        output.close()

        try:
            file_path = await ft.FilePicker().save_file(
                dialog_title="Save Primers",
                file_name="primers.csv",
                allowed_extensions=["csv"],
                file_type=ft.FilePickerFileType.CUSTOM,
                src_bytes=csv_content.encode("utf-8"),
            )
            if self.app_page.web:
                self._show_snackbar("Primers ready for download!")
            else:
                if file_path is None:
                    return
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(csv_content)
                self._show_snackbar(f"Saved {len(active_primers)} primer(s).")
        except Exception as ex:
            self._show_snackbar(f"Error saving file: {ex}")

    def _show_snackbar(self, message: str) -> None:
        """Show a snackbar message."""
        self.app_page.overlay.append(ft.SnackBar(ft.Text(message), open=True))
        self.app_page.update()


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
        """Initialize the InputView."""
        super().__init__(
            expand=True, vertical_alignment=ft.CrossAxisAlignment.STRETCH
        )
        self.app_page = page
        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()
        self.on_change = on_change
        self.on_stop_editing_callback = on_stop_editing
        self._focus_timer: threading.Timer | None = None

        self.template_input = TemplateInput(
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
        )

        self.right_fraction = 0.5

        self.divider = ft.GestureDetector(
            on_pan_update=self._on_pan_update,
            content=ft.Container(
                width=5,
                bgcolor=GUIColors.DIVIDER_GREY,
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
    def primer_info_panel(self) -> ft.Container:
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
    def validation_errors(self) -> list[str | None]:
        """Get the list of validation errors."""
        return self.primer_input.validation_errors

    @validation_errors.setter
    def validation_errors(self, val: list[str | None]) -> None:
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

    # Control layout compatibility properties
    @property
    def top_container(self) -> ft.Container:
        """Get the template input container (layout backward compatibility)."""
        return self.template_input

    @property
    def bottom_container(self) -> ft.Container:
        """Get the primer input container (layout backward compatibility)."""
        return self.primer_input

    def _handle_row_click(self, idx: int, name_edit: ft.TextField) -> None:
        """Delegate row click handling to PrimerInput."""
        self.primer_input._handle_row_click(idx, name_edit)

    def _handle_field_focus(self, e: ft.ControlEvent) -> None:
        """Handle focus on input fields to cancel auto-trigger timer."""
        if self._focus_timer is not None:
            self._focus_timer.cancel()
            self._focus_timer = None

        if e.control.data is not None:
            self.primer_input.focused_primer_index = e.control.data
            self.primer_input._update_row_highlights()
            self.primer_input._update_primer_info_panel()

    def _handle_field_blur(self, e: ft.ControlEvent) -> None:
        """Handle blur on input fields to trigger results page after a delay."""
        if getattr(self, "_skip_blur_timer", False):
            self._skip_blur_timer = False
            return

        self.sync_to_state(rebuild_if_needed=False)

        # If the control that blurred has a validation error,
        # update its display immediately.
        if isinstance(e.control, ft.TextField) and e.control.data is not None:
            idx = e.control.data
            if idx < len(self.primer_input.validation_errors):
                err = self.primer_input.validation_errors[idx]
                e.control.error = err
                e.control.height = 30 if not err else None
                for container in self.primer_input.primers_list.controls:
                    if (
                        isinstance(container, ft.Container)
                        and container.data == idx
                    ):
                        container.height = 30 if not err else None
                        break
                self.app_page.update()

        if self._focus_timer is not None:
            self._focus_timer.cancel()
            self._focus_timer = None

        def timer_callback() -> None:
            if not self.page:
                return
            self.sync_to_state(rebuild_if_needed=True)
            if self.on_stop_editing_callback:
                self.on_stop_editing_callback()

        self._focus_timer = threading.Timer(0.15, timer_callback)
        self._focus_timer.daemon = True
        try:
            self._focus_timer.start()
        except RuntimeError:
            # Pyodide (WebAssembly) may not support starting OS threads when
            # SharedArrayBuffer / cross-origin isolation is unavailable.
            # Fall back to a direct synchronous invocation of the callback.
            self._focus_timer = None
            timer_callback()

    def _handle_field_submit(self, e: ft.ControlEvent) -> None:
        """Handle submission (Enter key) to immediately trigger results."""
        if self._focus_timer is not None:
            self._focus_timer.cancel()
            self._focus_timer = None
        self.sync_to_state()
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback()

    def will_unmount(self) -> None:
        """Clean up when the view is unmounted."""
        if self._focus_timer is not None:
            self._focus_timer.cancel()
            self._focus_timer = None

    def sync_to_state(self, rebuild_if_needed: bool = True) -> None:
        """Sync current UI controls back to the central state."""
        self.template_input.sync_to_state()
        self.primer_input.sync_to_state(rebuild_if_needed=rebuild_if_needed)

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        self.template_input.update_ui()
        self.primer_input.update_ui()

    def _on_change_handler(self, e: ft.ControlEvent) -> None:
        """Handle change in input fields."""
        self.sync_to_state()
        self.primer_input._update_primer_info_panel()
        if self.on_change:
            self.on_change(e)

    def _clear_primers(self, e: ft.ControlEvent) -> None:
        """Clear all primers."""
        self.input_data.primers = []
        self.primer_input.focused_primer_index = None
        self.update_ui()
        if self.on_change:
            self.on_change(e)
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback()

    def _clear_template(self, e: ft.ControlEvent) -> None:
        """Clear the DNA template."""
        self.template_input.template_sequence.value = ""
        self.sync_to_state()
        if self.on_change:
            self.on_change(e)
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback()

    def _on_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle resizing the bottom (right) container via the divider."""
        page_width = self.app_page.width
        if isinstance(page_width, (int, float)) and page_width > 0:
            delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
            # Calculate current pixel width of the right container
            current_width = page_width * self.right_fraction
            new_width = max(200.0, current_width - delta_x)
            # Ensure the left container stays at least 200px wide
            new_width = min(new_width, page_width - 200.0)

            # Recalculate fractions and set relative expand weights
            self.right_fraction = new_width / page_width
            self.primer_input.expand = int(self.right_fraction * 1000)
            self.template_input.expand = int((1.0 - self.right_fraction) * 1000)
            self.app_page.update()

    def get_template(self) -> str:
        """Get the current template sequence."""
        return self.input_data.template

    def is_circular(self) -> bool:
        """Check if the template is circular."""
        return self.input_data.template_circular

    def get_primers(self) -> list[dict[str, Any]]:
        """Get the list of active primers."""
        return self.input_data.get_active_primers()

    def get_all_primers_state(self) -> list[dict[str, Any]]:
        """Get all primers (active and inactive) for serialization."""
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
        """Get the current input data state for serialization."""
        self.sync_to_state()
        return self.input_data.to_dict()

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current input data from deserialized data."""
        self.input_data.from_dict(state)
        self.update_ui()
        self.app_page.update()

    def _update_row_highlights(self) -> None:
        """Update background colors of all row containers."""
        self.primer_input._update_row_highlights()

    def _update_primer_info_panel(self) -> None:
        """Update the primer information panel based on the focused primer."""
        self.primer_input._update_primer_info_panel()
