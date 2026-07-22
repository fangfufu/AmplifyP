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

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.dna import DNA, DNADirection
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.data_helpers import clean_sequence
from amplifyp.primer_designer_2d import FilterMetric


class Designer2DForm(ft.Column):  # type: ignore[misc]
    """Form controls and validation logic for 2D Primer Designer."""

    def __init__(
        self,
        settings: GUISettings,
        on_submit_callback: Callable[[], Any],
        on_clear_error_callback: Callable[[ft.ControlEvent], None]
        | None = None,
    ) -> None:
        """Initialise the Designer2DForm."""
        super().__init__(spacing=8)
        self.settings = settings
        self.on_submit_callback = on_submit_callback
        self.on_clear_error_callback = on_clear_error_callback

        # Input fields
        self.fwd_dna_input = ft.TextField(
            label="Forward DNA Sequence",
            hint_text="e.g. ATGCGTACGT...",
            expand=True,
            multiline=False,
            autofocus=True,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.fwd_min_len_input = ft.TextField(
            label="Fwd Min Len (bp)",
            value="18",
            width=120,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.rev_dna_input = ft.TextField(
            label="Reverse DNA Sequence",
            hint_text="e.g. CGTACGATGC...",
            expand=True,
            multiline=False,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.rev_min_len_input = ft.TextField(
            label="Rev Min Len (bp)",
            value="18",
            width=120,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )

        pd_settings = self.settings.get_primer_dimer_settings()

        self.quality_filter_input = ft.TextField(
            label="Quality Filter",
            hint_text="Unconstrained if empty",
            value=f"{pd_settings.threshold:.1f}",
            expand=True,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.overlap_filter_input = ft.TextField(
            label="Overlap Filter (bp)",
            hint_text="Unconstrained if empty",
            value=str(pd_settings.min_overlap),
            expand=True,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.filter_metric_dropdown = ft.Dropdown(
            label="Metric",
            options=[
                ft.dropdown.Option("MAX", "Max"),
                ft.dropdown.Option("MEAN", "Mean"),
            ],
            value="MAX",
            width=110,
        )
        self.analyse_button = ft.FilledButton(
            "Analyse",
            icon=ft.Icons.PLAY_ARROW,
            tooltip="Run 2D Primer Truncation Analysis",
            on_click=self._on_submit_event,
        )
        self.error_text = ft.Text(
            "", color=GUIColours.ERROR_RED, visible=False, size=12
        )

        self.controls = [
            ft.Text(
                "2D Truncation Parameters",
                weight=ft.FontWeight.BOLD,
                size=self.settings.get("font_size_subheader", 16),
            ),
            ft.Row([self.fwd_dna_input, self.fwd_min_len_input], spacing=8),
            ft.Row([self.rev_dna_input, self.rev_min_len_input], spacing=8),
            ft.Row(
                [
                    self.quality_filter_input,
                    self.overlap_filter_input,
                    self.filter_metric_dropdown,
                    ft.Container(
                        content=self.analyse_button,
                        margin=ft.Margin.only(left=8),
                    ),
                ],
                alignment=ft.MainAxisAlignment.START,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=8,
            ),
            self.error_text,
        ]

    def _on_submit_event(self, e: ft.ControlEvent) -> None:
        """Handle submit/click events from form controls."""
        self.on_submit_callback()

    def _clear_field_error(self, e: ft.ControlEvent) -> None:
        """Clear error text when user edits an input field."""
        ctrl = getattr(e, "control", None)
        if isinstance(ctrl, ft.TextField) and ctrl.error:
            ctrl.error = None
            try:
                if self.page:
                    self.page.update()
            except RuntimeError:
                pass
        if self.on_clear_error_callback:
            self.on_clear_error_callback(e)

    def clear_errors(self) -> None:
        """Clear all field error indicators and general error message."""
        self.fwd_dna_input.error = None
        self.fwd_min_len_input.error = None
        self.rev_dna_input.error = None
        self.rev_min_len_input.error = None
        self.quality_filter_input.error = None
        self.overlap_filter_input.error = None
        self.error_text.visible = False
        self.error_text.value = ""

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

        # Clean and validate Forward DNA sequence
        raw_fwd_seq = self.fwd_dna_input.value or ""
        cleaned_fwd_seq = clean_sequence(raw_fwd_seq)
        if not cleaned_fwd_seq:
            self.fwd_dna_input.error = "Forward DNA sequence cannot be empty"
            self.error_text.value = "Forward DNA sequence cannot be empty"
            self.error_text.visible = True
            raise ValueError("Forward DNA sequence cannot be empty")

        # Validate Forward Min Length
        try:
            fwd_min_len = int(self.fwd_min_len_input.value or "18")
            if fwd_min_len <= 0:
                raise ValueError
        except ValueError as err:
            self.fwd_min_len_input.error = "Must be > 0"
            self.error_text.value = (
                "Forward minimum length must be a positive integer"
            )
            self.error_text.visible = True
            raise ValueError(
                "Forward min length must be a positive integer"
            ) from err

        if len(cleaned_fwd_seq) < fwd_min_len:
            self.fwd_min_len_input.error = (
                f"Exceeds sequence length ({len(cleaned_fwd_seq)})"
            )
            self.error_text.value = (
                f"Forward minimum length ({fwd_min_len}) cannot exceed "
                f"sequence length ({len(cleaned_fwd_seq)})"
            )
            self.error_text.visible = True
            raise ValueError("Forward min length exceeds sequence length")

        # Clean and validate Reverse DNA sequence
        raw_rev_seq = self.rev_dna_input.value or ""
        cleaned_rev_seq = clean_sequence(raw_rev_seq)
        if not cleaned_rev_seq:
            self.rev_dna_input.error = "Reverse DNA sequence cannot be empty"
            self.error_text.value = "Reverse DNA sequence cannot be empty"
            self.error_text.visible = True
            raise ValueError("Reverse DNA sequence cannot be empty")

        # Validate Reverse Min Length
        try:
            rev_min_len = int(self.rev_min_len_input.value or "18")
            if rev_min_len <= 0:
                raise ValueError
        except ValueError as err:
            self.rev_min_len_input.error = "Must be > 0"
            self.error_text.value = (
                "Reverse minimum length must be a positive integer"
            )
            self.error_text.visible = True
            raise ValueError(
                "Reverse min length must be a positive integer"
            ) from err

        if len(cleaned_rev_seq) < rev_min_len:
            self.rev_min_len_input.error = (
                f"Exceeds sequence length ({len(cleaned_rev_seq)})"
            )
            self.error_text.value = (
                f"Reverse minimum length ({rev_min_len}) cannot exceed "
                f"sequence length ({len(cleaned_rev_seq)})"
            )
            self.error_text.visible = True
            raise ValueError("Reverse min length exceeds sequence length")

        # Quality filter
        threshold: float | None = None
        q_str = (self.quality_filter_input.value or "").strip()
        if q_str:
            try:
                threshold = float(q_str)
                if threshold < 0:
                    raise ValueError
            except ValueError as err:
                self.quality_filter_input.error = "Must be >= 0"
                self.error_text.value = "Quality filter must be non-negative"
                self.error_text.visible = True
                raise ValueError("Quality filter must be non-negative") from err

        # Overlap filter
        max_overlap: int | None = None
        o_str = (self.overlap_filter_input.value or "").strip()
        if o_str:
            try:
                max_overlap = int(o_str)
                if max_overlap < 0:
                    raise ValueError
            except ValueError as err:
                self.overlap_filter_input.error = "Must be >= 0"
                self.error_text.value = "Overlap filter must be non-negative"
                self.error_text.visible = True
                raise ValueError("Overlap filter must be non-negative") from err

        # Filter metric
        metric_str = (self.filter_metric_dropdown.value or "MAX").lower()
        filter_metric = (
            FilterMetric.MEAN if metric_str == "mean" else FilterMetric.MAX
        )

        fwd_dna = DNA(cleaned_fwd_seq, direction=DNADirection.FWD)
        rev_dna = DNA(cleaned_rev_seq, direction=DNADirection.REV)

        return (
            fwd_dna,
            fwd_min_len,
            rev_dna,
            rev_min_len,
            threshold,
            max_overlap,
            filter_metric,
        )
