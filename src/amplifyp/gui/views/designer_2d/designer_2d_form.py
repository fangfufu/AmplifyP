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

"""Form component and input validation for 2D Primer Designer View."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.dna import DNA, DNADirection
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.data_helpers import clean_sequence
from amplifyp.gui.views.designer import BaseDesignerForm, create_field_container
from amplifyp.primer_designer_2d import FilterMetric


class Designer2DForm(BaseDesignerForm):
    """Form controls and validation logic for 2D Primer Designer."""

    def __init__(
        self,
        settings: GUISettings,
        on_submit_callback: Callable[[], Any],
        on_clear_error_callback: Callable[[Any], None] | None = None,
        on_save_callback: Callable[[ft.ControlEvent], Any] | None = None,
        on_load_callback: Callable[[ft.ControlEvent], Any] | None = None,
        on_clear_all_callback: Callable[[ft.ControlEvent], Any] | None = None,
    ) -> None:
        """Initialise the Designer2DForm."""
        super().__init__(
            settings=settings,
            on_submit_callback=on_submit_callback,
            on_clear_error_callback=on_clear_error_callback,
            on_save_callback=on_save_callback,
            on_load_callback=on_load_callback,
            on_clear_all_callback=on_clear_all_callback,
            spacing=8,
        )

        # 2D-specific input fields
        self.fwd_dna_input = ft.TextField(
            hint_text="e.g. ATGCGTACGT...",
            expand=True,
            multiline=False,
            autofocus=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.fwd_min_len_input = ft.TextField(
            hint_text="e.g. 18",
            value="",
            width=160,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.rev_dna_input = ft.TextField(
            hint_text="e.g. CGTACGATGC...",
            expand=True,
            multiline=False,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.rev_min_len_input = ft.TextField(
            hint_text="e.g. 18",
            value="",
            width=160,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )

        self.controls = [
            self._build_header_container("2D Truncation Parameters"),
            ft.Row(
                [
                    create_field_container(
                        "Forward Candidate Primer Sequence",
                        self.fwd_dna_input,
                        expand=True,
                    ),
                    create_field_container(
                        "Fwd Min Length (nt)",
                        self.fwd_min_len_input,
                        expand=False,
                        width=160,
                    ),
                ],
                spacing=8,
            ),
            ft.Row(
                [
                    create_field_container(
                        "Reverse Candidate Primer Sequence",
                        self.rev_dna_input,
                        expand=True,
                    ),
                    create_field_container(
                        "Rev Min Length (nt)",
                        self.rev_min_len_input,
                        expand=False,
                        width=160,
                    ),
                ],
                spacing=8,
            ),
            self._build_filter_row(),
            self.error_text,
        ]

    # --- Property aliases for backwards compatibility ---
    @property
    def quality_filter_input(self) -> ft.TextField:
        """Alias for max_quality_input."""
        return self.max_quality_input

    @property
    def overlap_filter_input(self) -> ft.TextField:
        """Alias for max_overlap_input."""
        return self.max_overlap_input

    def clear_errors(self) -> None:
        """Clear all field error indicators and general error message."""
        super().clear_errors()
        self.fwd_dna_input.error = None
        self.fwd_min_len_input.error = None
        self.rev_dna_input.error = None
        self.rev_min_len_input.error = None

    def _validate_primer_input(
        self,
        dna_input: ft.TextField,
        min_len_input: ft.TextField,
        name: str,
    ) -> tuple[str, int, bool]:
        """Validate sequence and min length for a candidate primer."""
        has_error = False
        raw_seq = dna_input.value or ""
        clean_seq = clean_sequence(raw_seq)
        if not clean_seq:
            dna_input.error = (
                f"{name} candidate primer sequence cannot be empty"
            )
            has_error = True

        min_len_raw = (min_len_input.value or "").strip()
        min_len = 0
        min_len_valid = False
        if not min_len_raw:
            min_len_input.error = "Must be > 0"
            has_error = True
        else:
            try:
                min_len = int(min_len_raw)
                if min_len <= 0:
                    raise ValueError
                min_len_valid = True
            except ValueError:
                min_len_input.error = "Must be > 0"
                has_error = True

        if min_len_valid and clean_seq and len(clean_seq) < min_len:
            min_len_input.error = f"Exceeds sequence length ({len(clean_seq)})"
            has_error = True

        return clean_seq, min_len, not has_error

    def validate_and_get_params(
        self,
    ) -> tuple[DNA, int, DNA, int, float | None, int | None, FilterMetric]:
        """Validate input fields and return 2D primer designer parameters.

        Returns:
            Tuple of (fwd_dna, fwd_min_length, rev_dna, rev_min_length,
                threshold, max_overlap, filter_metric).

        Raises:
            ValueError: If any input field contains invalid data.
        """
        self.clear_errors()

        fwd_seq, fwd_min_len, fwd_valid = self._validate_primer_input(
            self.fwd_dna_input, self.fwd_min_len_input, "Forward"
        )
        rev_seq, rev_min_len, rev_valid = self._validate_primer_input(
            self.rev_dna_input, self.rev_min_len_input, "Reverse"
        )

        threshold, q_valid = self.validate_max_quality()
        if not q_valid:
            self.max_quality_input.error = "Must be >= 0"

        max_overlap, o_valid = self.validate_max_overlap()
        if not o_valid:
            self.max_overlap_input.error = "Must be >= 0"

        if not (fwd_valid and rev_valid and q_valid and o_valid):
            try:
                if self.page:
                    self.page.update()
            except RuntimeError:
                pass
            raise ValueError("Input validation failed")

        filter_metric = FilterMetric.MAX
        fwd_dna = DNA(fwd_seq, direction=DNADirection.FWD)
        rev_dna = DNA(rev_seq, direction=DNADirection.REV)

        return (
            fwd_dna,
            fwd_min_len,
            rev_dna,
            rev_min_len,
            threshold,
            max_overlap,
            filter_metric,
        )
