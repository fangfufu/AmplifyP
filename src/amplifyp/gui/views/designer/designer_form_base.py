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

"""Base form component and input validation for primer designer views."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings


def create_field_container(
    label_text: str,
    field: ft.Control,
    expand: bool | int | None = True,
    width: int | None = None,
) -> ft.Column:
    """Create a column with a fixed header label above an input control.

    Args:
        label_text: Text displayed above the control.
        field: Input control (e.g. TextField).
        expand: Expand property for the container. Defaults to True.
        width: Optional fixed width in pixels.

    Returns:
        Flet Column containing label and input field.
    """
    return ft.Column(
        [
            ft.Text(
                label_text,
                size=12,
                weight=ft.FontWeight.W_500,
                color=GUIColours.TEXT_ON_SURFACE,
            ),
            field,
        ],
        spacing=2,
        expand=expand,
        width=width,
    )


class BaseDesignerForm(ft.Column):  # type: ignore[misc]
    """Base form providing standard designer controls and validation."""

    def __init__(
        self,
        settings: GUISettings,
        on_submit_callback: Callable[[], Any],
        on_clear_error_callback: Callable[[Any], None] | None = None,
        on_save_callback: Callable[[ft.ControlEvent], Any] | None = None,
        on_load_callback: Callable[[ft.ControlEvent], Any] | None = None,
        on_clear_all_callback: Callable[[ft.ControlEvent], Any] | None = None,
        spacing: int = 8,
    ) -> None:
        """Initialise BaseDesignerForm.

        Args:
            settings: GUI settings object.
            on_submit_callback: Callback when user submits or analyses.
            on_clear_error_callback: Callback when input changes.
            on_save_callback: Callback for save button.
            on_load_callback: Callback for load button.
            on_clear_all_callback: Callback for clear all button.
            spacing: Vertical spacing between form rows.
        """
        super().__init__(spacing=spacing)
        self.settings = settings
        self.on_submit_callback = on_submit_callback
        self.on_clear_error_callback = on_clear_error_callback

        # Action buttons
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
        self.clear_all_button = ft.FilledTonalButton(
            "Clear All",
            icon=ft.Icons.CLEAR_ALL,
            tooltip="Clear all parameters and results",
            on_click=on_clear_all_callback,
            height=32,
        )

        # Common filter inputs
        self.max_quality_input = ft.TextField(
            hint_text="Unconstrained if empty",
            value="",
            expand=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )
        self.max_overlap_input = ft.TextField(
            hint_text="Unconstrained if empty",
            value="",
            expand=True,
            border_color=GUIColours.OUTLINE,
            on_submit=self._on_submit_event,
            on_change=self._clear_field_error,
        )

        # Analyse button & error text
        self.analyse_button = ft.FilledButton(
            "Analyse",
            icon=ft.Icons.PLAY_ARROW,
            tooltip="Run Primer Truncation Analysis",
            on_click=self._on_submit_event,
        )
        self.error_text = ft.Text(
            "", color=GUIColours.ERROR_RED, visible=False, size=12
        )

    def _build_header_container(self, title: str) -> ft.Container:
        """Build standard form header row with title and action buttons.

        Args:
            title: Header title text.

        Returns:
            Flet Container with title and action buttons.
        """
        return ft.Container(
            content=ft.Row(
                [
                    ft.Text(
                        title,
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_subheader", 16),
                    ),
                    ft.Row(
                        [
                            self.load_button,
                            self.save_button,
                            self.clear_all_button,
                        ],
                        spacing=4,
                        tight=True,
                    ),
                ],
                alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
            ),
            margin=ft.Margin.only(bottom=6),
        )

    def _build_filter_row(
        self, extra_controls: list[ft.Control] | None = None
    ) -> ft.Row:
        """Build the standard bottom filter row with filter controls.

        Args:
            extra_controls: Optional additional controls to prepend or append.

        Returns:
            Flet Row containing filter controls.
        """
        row_controls: list[ft.Control] = []
        if extra_controls:
            row_controls.extend(extra_controls)

        row_controls.extend(
            [
                create_field_container(
                    "Max Quality", self.max_quality_input, expand=True
                ),
                create_field_container(
                    "Max Overlap (bp)", self.max_overlap_input, expand=True
                ),
                ft.Container(
                    content=self.analyse_button,
                    margin=ft.Margin.only(top=18, left=8),
                ),
            ]
        )

        return ft.Row(
            row_controls,
            alignment=ft.MainAxisAlignment.START,
            vertical_alignment=ft.CrossAxisAlignment.START,
            spacing=8,
        )

    def _on_submit_event(self, e: Any) -> None:
        """Handle submit/click events from form controls."""
        self.on_submit_callback()

    def _clear_field_error(self, e: Any) -> None:
        """Clear error highlight when user edits an input field."""
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
        self.max_quality_input.error = None
        self.max_overlap_input.error = None
        self.error_text.visible = False
        self.error_text.value = ""

    def show_field_error(self, field: ft.TextField, message: str) -> None:
        """Set error text on a specific input field.

        Args:
            field: The target TextField control.
            message: The error message to display.
        """
        field.error = message
        try:
            if self.page:
                self.page.update()
        except RuntimeError:
            pass

    def show_error(self, message: str) -> None:
        """Display general validation or execution error message.

        Args:
            message: The error message to display.
        """
        self.error_text.value = message
        self.error_text.visible = True

    def validate_max_quality(
        self, int_only: bool = False
    ) -> tuple[float | None, bool]:
        """Validate max quality input field.

        Args:
            int_only: Whether quality value must be an integer. Defaults to
                False.

        Returns:
            Tuple of (quality_value_or_none, is_valid).
        """
        raw = (self.max_quality_input.value or "").strip()
        if not raw:
            return None, True
        if int_only:
            try:
                int_val = int(raw)
                if int_val < 0:
                    self.show_field_error(
                        self.max_quality_input,
                        "Max Quality must be a non-negative integer.",
                    )
                    return None, False
                return float(int_val), True
            except ValueError:
                self.show_field_error(
                    self.max_quality_input,
                    "Max Quality must be an integer.",
                )
                return None, False
        try:
            val = float(raw)
            if val < 0:
                self.show_field_error(
                    self.max_quality_input,
                    "Max Quality must be non-negative.",
                )
                return None, False
            return val, True
        except ValueError:
            self.show_field_error(
                self.max_quality_input,
                "Max Quality must be an integer.",
            )
            return None, False

    def validate_max_overlap(self) -> tuple[int | None, bool]:
        """Validate max overlap input field.

        Returns:
            Tuple of (overlap_value_or_none, is_valid).
        """
        raw = (self.max_overlap_input.value or "").strip()
        if not raw:
            return None, True
        try:
            val = int(raw)
            if val < 0:
                self.show_field_error(
                    self.max_overlap_input,
                    "Max Overlap must be non-negative.",
                )
                return None, False
            return val, True
        except ValueError:
            self.show_field_error(
                self.max_overlap_input,
                "Max Overlap must be a non-negative integer.",
            )
            return None, False
