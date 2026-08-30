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

"""Form component and input validation for 1D Primer Designer View."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.dna import DNADirection
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.data_helpers import clean_sequence
from amplifyp.gui.views.designer import BaseDesignerForm, create_field_container


class Designer1DForm(BaseDesignerForm):
    """Form controls and validation logic for 1D Primer Designer."""

    def __init__(
        self,
        settings: GUISettings,
        on_submit_callback: Callable[[], Any],
        on_clear_error_callback: Callable[[Any], None] | None = None,
        on_save_callback: Callable[[ft.ControlEvent], Any] | None = None,
        on_load_callback: Callable[[ft.ControlEvent], Any] | None = None,
        on_clear_all_callback: Callable[[ft.ControlEvent], Any] | None = None,
    ) -> None:
        """Initialise the Designer1DForm."""
        super().__init__(
            settings=settings,
            on_submit_callback=on_submit_callback,
            on_clear_error_callback=on_clear_error_callback,
            on_save_callback=on_save_callback,
            on_load_callback=on_load_callback,
            on_clear_all_callback=on_clear_all_callback,
            spacing=8,
        )

        # 1D-specific input fields
        self.dna_input = ft.TextField(
            hint_text="e.g. ATGCGTACGT...",
            expand=True,
            multiline=False,
            autofocus=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._on_dna_change,
        )
        self.length_display = ft.TextField(
            value="0",
            read_only=True,
            width=90,
            text_align=ft.TextAlign.CENTER,
            border_color=GUIColours.OUTLINE,
            content_padding=ft.Padding(8, 4, 8, 4),
        )
        self.min_len_input = ft.TextField(
            hint_text="e.g. 18",
            value="",
            expand=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )

        self.controls = [
            self._build_header_container("1D Truncation Parameters"),
            ft.Row(
                [
                    create_field_container(
                        "Candidate Primer Sequence", self.dna_input, expand=True
                    ),
                    create_field_container(
                        "Length (nt)",
                        self.length_display,
                        expand=False,
                        width=90,
                    ),
                ],
                spacing=8,
            ),
            self._build_filter_row(
                extra_controls=[
                    create_field_container(
                        "Min Length (nt)", self.min_len_input, expand=True
                    ),
                ]
            ),
            self.error_text,
        ]

    def _on_dna_change(self, e: ft.ControlEvent) -> None:
        """Update length counter and clear error when DNA input changes."""
        dna_raw = self.dna_input.value or ""
        cleaned = clean_sequence(dna_raw)
        self.length_display.value = str(len(cleaned))
        self._clear_field_error(e)

    def clear_errors(self) -> None:
        """Clear all field error indicators and general error message."""
        super().clear_errors()
        self.dna_input.error = None
        self.min_len_input.error = None

    def validate_and_get_params(
        self,
    ) -> tuple[str, int, DNADirection, float | None, int | None] | None:
        """Validate input fields and return parsed parameters tuple if valid."""
        self.clear_errors()

        dna_raw = (self.dna_input.value or "").strip()
        if not dna_raw:
            self.show_field_error(
                self.dna_input, "Please enter a valid DNA sequence."
            )
            return None

        clean_seq = clean_sequence(dna_raw)
        if not clean_seq:
            self.show_field_error(
                self.dna_input, "DNA sequence contains no valid nucleotides."
            )
            return None

        min_len_raw = (self.min_len_input.value or "").strip()
        if not min_len_raw:
            self.show_field_error(
                self.min_len_input, "Minimum length is required."
            )
            return None

        if not min_len_raw.isdigit():
            self.show_field_error(
                self.min_len_input, "Minimum length must be a positive integer."
            )
            return None

        min_length = int(min_len_raw)
        if min_length <= 0:
            self.show_field_error(
                self.min_len_input, "Minimum length must be greater than 0."
            )
            return None

        if min_length > len(clean_seq):
            self.show_field_error(
                self.min_len_input,
                f"Minimum length ({min_length}) cannot exceed sequence length "
                f"({len(clean_seq)}).",
            )
            return None

        mode = DNADirection.FWD
        threshold, q_valid = self.validate_max_quality(int_only=True)
        if not q_valid:
            return None

        max_overlap, o_valid = self.validate_max_overlap()
        if not o_valid:
            return None

        return clean_seq, min_length, mode, threshold, max_overlap
