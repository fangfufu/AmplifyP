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

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.dna import DNADirection
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.data_helpers import clean_sequence


class Designer1DForm(ft.Column):  # type: ignore[misc]
    """Form controls and validation logic for 1D Primer Designer."""

    def __init__(
        self,
        settings: GUISettings,
        on_submit_callback: Callable[[], Any],
        on_clear_error_callback: Callable[[ft.ControlEvent], None]
        | None = None,
        on_save_callback: Callable[[ft.ControlEvent], Any] | None = None,
        on_load_callback: Callable[[ft.ControlEvent], Any] | None = None,
    ) -> None:
        """Initialise the Designer1DForm."""
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
        self.dna_input = ft.TextField(
            label="DNA Sequence",
            hint_text="e.g. ATGCGTACGT...",
            expand=True,
            multiline=False,
            autofocus=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.min_len_input = ft.TextField(
            label="Min Length (bp)",
            value="18",
            expand=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.mode_dropdown = ft.Dropdown(
            label="Direction",
            options=[
                ft.dropdown.Option("FWD", "Forward"),
                ft.dropdown.Option("REV", "Reverse"),
            ],
            value="FWD",
            expand=True,
            border_color=GUIColours.OUTLINE,
        )
        pd_settings = self.settings.get_primer_dimer_settings()

        self.max_quality_input = ft.TextField(
            label="Max Quality",
            hint_text="Unconstrained if empty",
            value=f"{pd_settings.threshold:.1f}",
            expand=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.max_overlap_input = ft.TextField(
            label="Max Overlap (bp)",
            hint_text="Unconstrained if empty",
            value=str(pd_settings.min_overlap),
            expand=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.analyse_button = ft.FilledButton(
            "Analyse",
            icon=ft.Icons.PLAY_ARROW,
            tooltip="Run 1D Primer Truncation Analysis",
            on_click=self._on_submit_event,
        )
        self.error_text = ft.Text(
            "", color=GUIColours.ERROR_RED, visible=False, size=12
        )

        self.controls = [
            ft.Row(
                [
                    ft.Text(
                        "1D Truncation Parameters",
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
            ft.Row([self.dna_input]),
            ft.Row(
                [self.min_len_input, self.mode_dropdown],
                alignment=ft.MainAxisAlignment.START,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=14,
            ),
            ft.Row(
                [
                    self.max_quality_input,
                    self.max_overlap_input,
                    ft.Container(
                        content=self.analyse_button,
                        margin=ft.Margin.only(left=24),
                    ),
                ],
                alignment=ft.MainAxisAlignment.START,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=14,
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
        self.dna_input.error = None
        self.min_len_input.error = None
        self.max_quality_input.error = None
        self.max_overlap_input.error = None
        self.error_text.visible = False
        self.error_text.value = ""

    def show_field_error(self, field: ft.TextField, message: str) -> None:
        """Set error text on a specific field and display general error."""
        field.error = message
        self.show_error(message)
        try:
            if self.page:
                self.page.update()
        except RuntimeError:
            pass

    def show_error(self, message: str) -> None:
        """Display validation error message."""
        self.error_text.value = message
        self.error_text.visible = True

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

        mode_val = self.mode_dropdown.value or "FWD"
        mode = DNADirection.FWD if mode_val == "FWD" else DNADirection.REV

        max_q_raw = (self.max_quality_input.value or "").strip()
        threshold: float | None = None
        if max_q_raw:
            try:
                threshold = float(max_q_raw)
                if threshold < 0:
                    self.show_field_error(
                        self.max_quality_input,
                        "Max Quality must be a non-negative number.",
                    )
                    return None
            except ValueError:
                self.show_field_error(
                    self.max_quality_input,
                    "Max Quality must be a valid number.",
                )
                return None

        max_overlap_raw = (self.max_overlap_input.value or "").strip()
        max_overlap: int | None = None
        if max_overlap_raw:
            if not max_overlap_raw.isdigit():
                self.show_field_error(
                    self.max_overlap_input,
                    "Max Overlap must be a non-negative integer.",
                )
                return None
            max_overlap = int(max_overlap_raw)
            if max_overlap < 0:
                self.show_field_error(
                    self.max_overlap_input,
                    "Max Overlap must be a non-negative integer.",
                )
                return None

        return clean_seq, min_length, mode, threshold, max_overlap
