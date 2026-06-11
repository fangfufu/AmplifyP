# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""ListView for rendering, styling, and resizing primers."""

from typing import Any

import flet as ft

from amplifyp.gui.settings import GUIColors
from amplifyp.gui.util import clean_sequence

from .primer_row import PrimerRow


class PrimerList(ft.ListView):  # type: ignore[misc]
    """ListView for rendering, styling, and resizing primers."""

    def __init__(
        self,
        primer_input: Any,
    ) -> None:
        """Initialize the PrimerList."""
        super().__init__(
            expand=True, spacing=0, padding=0, scroll=ft.ScrollMode.ALWAYS
        )
        self.primer_input = primer_input

    def update_list_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        font_family = self.primer_input.settings.get(
            "font_family", "Roboto Mono"
        )
        self.controls.clear()

        # Filter out all empty primers from state first
        clean_primers = [
            p
            for p in self.primer_input.input_data.primers
            if str(p.get("name", "")).strip()
            or clean_sequence(str(p.get("seq", ""))).strip()
        ]

        # Always append a trailing empty row
        clean_primers.append({"name": "", "seq": "", "active": False})
        self.primer_input.input_data.primers = clean_primers

        # Precompute duplicate indices & name/seq counts for error reporting
        names_count: dict[str, int] = {}
        seqs_count: dict[str, int] = {}
        for p in self.primer_input.input_data.primers:
            n_lower = str(p.get("name", "")).strip().lower()
            s_lower = clean_sequence(str(p.get("seq", ""))).lower()
            if n_lower:
                names_count[n_lower] = names_count.get(n_lower, 0) + 1
            if s_lower:
                seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1

        dup_indices = self.primer_input._get_duplicate_indices()

        self.primer_input.validation_errors = []
        num_primers = len(self.primer_input.input_data.primers)
        for idx, p in enumerate(self.primer_input.input_data.primers):
            name_val = p["name"]
            seq_val = p["seq"]
            is_active = p.get("active", True)

            error_message = PrimerRow.validate(name_val, seq_val)
            if not error_message:
                n_lower = name_val.strip().lower()
                s_lower = clean_sequence(seq_val).lower()
                if n_lower and names_count.get(n_lower, 0) > 1:
                    error_message = "Duplicate primer name"
                elif s_lower and seqs_count.get(s_lower, 0) > 1:
                    error_message = "Duplicate primer sequence"

            self.primer_input.validation_errors.append(error_message)

            is_dup = idx in dup_indices
            is_focused = idx == self.primer_input.focused_primer_index
            is_last_row = idx == num_primers - 1
            is_penultimate_row = idx == num_primers - 2

            row = PrimerRow(
                idx=idx,
                name=name_val,
                seq=seq_val,
                is_active=is_active,
                is_dup=is_dup,
                error_message=error_message,
                font_family=font_family,
                name_column_width=self.primer_input.name_column_width,
                on_change_handler=self.primer_input.on_change_handler,
                handle_field_focus=self.primer_input.handle_field_focus,
                handle_field_blur=self.primer_input.handle_field_blur,
                handle_field_submit=self.primer_input.handle_field_submit,
                on_row_click=self.primer_input._handle_row_click,
                on_move_primer=self.primer_input._move_primer,
                on_divider_pan=self.primer_input._on_primer_divider_pan,
                is_focused=is_focused,
                is_last_row=is_last_row,
                is_penultimate_row=is_penultimate_row,
            )
            self.controls.append(row)

        self.update_row_highlights()
        self.primer_input._update_primer_info_panel()
        self.primer_input._update_header_checkbox_state()

    def update_row_highlights(self) -> None:
        """Update background colors of all row containers.

        Highlights rows based on selection and duplicates.
        """
        dup_indices = self.primer_input._get_duplicate_indices()

        for container in self.controls:
            if (
                isinstance(container, ft.Container)
                and container.data is not None
            ):
                c_idx = container.data
                is_dup = c_idx in dup_indices

                if c_idx == self.primer_input.focused_primer_index:
                    container.bgcolor = GUIColors.SELECTED_ROW_BG
                elif is_dup:
                    container.bgcolor = GUIColors.DUPLICATE_BG
                else:
                    container.bgcolor = None

                # Update reorder container visibility & text padding.
                row = container.content
                if isinstance(row, ft.Row) and len(row.controls) >= 5:
                    seq_stack = row.controls[4]
                    if (
                        isinstance(seq_stack, ft.Stack)
                        and len(seq_stack.controls) > 1
                    ):
                        reorder_container = seq_stack.controls[1]
                        seq_edit = seq_stack.controls[0]
                        is_focused = (
                            c_idx == self.primer_input.focused_primer_index
                        )
                        reorder_container.visible = is_focused
                        if isinstance(seq_edit, ft.TextField):
                            seq_edit.content_padding = ft.Padding(
                                5, 0, 60 if is_focused else 0, 0
                            )
