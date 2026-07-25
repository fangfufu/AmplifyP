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
        on_clear_error_callback: Callable[[Any], None] | None = None,
        on_save_callback: Callable[[ft.ControlEvent], Any] | None = None,
        on_load_callback: Callable[[ft.ControlEvent], Any] | None = None,
    ) -> None:
        """Initialise the Designer2DForm."""
        super().__init__(spacing=8)
        self.settings = settings
        self.on_submit_callback = on_submit_callback
        self.on_clear_error_callback = on_clear_error_callback

        # Save and Load buttons
        self.save_button = ft.FilledTonalButton(
            "Save",
            icon=ft.Icons.SAVE,
            tooltip="Save parameters to YAML",
            on_click=on_save_callback,
            height=32,
        )
        self.load_button = ft.FilledTonalButton(
            "Load",
            icon=ft.Icons.UPLOAD_FILE,
            tooltip="Load parameters from YAML",
            on_click=on_load_callback,
            height=32,
        )

        # Input fields
        self.fwd_dna_input = ft.TextField(
            label="Forward Candidate Primer Sequence",
            hint_text="e.g. ATGCGTACGT...",
            expand=True,
            multiline=False,
            autofocus=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.fwd_min_len_input = ft.TextField(
            label="Fwd Min Len (bp)",
            value="18",
            width=160,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.rev_dna_input = ft.TextField(
            label="Reverse Candidate Primer Sequence",
            hint_text="e.g. CGTACGATGC...",
            expand=True,
            multiline=False,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.rev_min_len_input = ft.TextField(
            label="Rev Min Len (bp)",
            value="18",
            width=160,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )

        pd_settings = self.settings.get_primer_dimer_settings()

        self.quality_filter_input = ft.TextField(
            label="Quality Filter",
            hint_text="Unconstrained if empty",
            value=f"{pd_settings.threshold:.1f}",
            expand=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.overlap_filter_input = ft.TextField(
            label="Overlap Filter (bp)",
            hint_text="Unconstrained if empty",
            value="",
            expand=True,
            border_color=GUIColours.OUTLINE,
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
            width=150,
            border_color=GUIColours.OUTLINE,
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
            ft.Container(
                content=ft.Row(
                    [
                        ft.Text(
                            "2D Truncation Parameters",
                            weight=ft.FontWeight.BOLD,
                            size=self.settings.get("font_size_subheader", 16),
                        ),
                        ft.Row(
                            [self.load_button, self.save_button],
                            spacing=4,
                            tight=True,
                        ),
                    ],
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    vertical_alignment=ft.CrossAxisAlignment.CENTER,
                ),
                margin=ft.Margin.only(bottom=6),
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

    def _on_submit_event(self, e: Any) -> None:
        """Handle submit/click events from form controls."""
        self.on_submit_callback()

    def _clear_field_error(self, e: Any) -> None:
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
        has_error = False

        # Clean and validate Forward DNA sequence
        raw_fwd_seq = self.fwd_dna_input.value or ""
        cleaned_fwd_seq = clean_sequence(raw_fwd_seq)
        if not cleaned_fwd_seq:
            self.fwd_dna_input.error = (
                "Forward candidate primer sequence cannot be empty"
            )
            has_error = True

        # Validate Forward Min Length
        fwd_min_len = 18
        fwd_min_len_valid = False
        try:
            fwd_min_len = int(self.fwd_min_len_input.value or "18")
            if fwd_min_len <= 0:
                raise ValueError
            fwd_min_len_valid = True
        except ValueError:
            self.fwd_min_len_input.error = "Must be > 0"
            has_error = True

        if (
            fwd_min_len_valid
            and cleaned_fwd_seq
            and len(cleaned_fwd_seq) < fwd_min_len
        ):
            self.fwd_min_len_input.error = (
                f"Exceeds sequence length ({len(cleaned_fwd_seq)})"
            )
            has_error = True

        # Clean and validate Reverse DNA sequence
        raw_rev_seq = self.rev_dna_input.value or ""
        cleaned_rev_seq = clean_sequence(raw_rev_seq)
        if not cleaned_rev_seq:
            self.rev_dna_input.error = (
                "Reverse candidate primer sequence cannot be empty"
            )
            has_error = True

        # Validate Reverse Min Length
        rev_min_len = 18
        rev_min_len_valid = False
        try:
            rev_min_len = int(self.rev_min_len_input.value or "18")
            if rev_min_len <= 0:
                raise ValueError
            rev_min_len_valid = True
        except ValueError:
            self.rev_min_len_input.error = "Must be > 0"
            has_error = True

        if (
            rev_min_len_valid
            and cleaned_rev_seq
            and len(cleaned_rev_seq) < rev_min_len
        ):
            self.rev_min_len_input.error = (
                f"Exceeds sequence length ({len(cleaned_rev_seq)})"
            )
            has_error = True

        # Quality filter
        threshold: float | None = None
        q_str = (self.quality_filter_input.value or "").strip()
        if q_str:
            try:
                threshold = float(q_str)
                if threshold < 0:
                    raise ValueError
            except ValueError:
                self.quality_filter_input.error = "Must be >= 0"
                has_error = True

        # Overlap filter
        max_overlap: int | None = None
        o_str = (self.overlap_filter_input.value or "").strip()
        if o_str:
            try:
                max_overlap = int(o_str)
                if max_overlap < 0:
                    raise ValueError
            except ValueError:
                self.overlap_filter_input.error = "Must be >= 0"
                has_error = True

        if has_error:
            try:
                if self.page:
                    self.page.update()
            except RuntimeError:
                pass
            raise ValueError("Input validation failed")

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
