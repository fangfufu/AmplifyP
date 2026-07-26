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

"""Input component for DNA template sequence."""

from __future__ import annotations

import logging
from collections.abc import Callable

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.data_helpers import clean_sequence
from amplifyp.gui.utils.gui_helpers import NotificationHelper

from .casing import change_selection_case
from .file_manager import load_template_click, save_template_click
from .formatter import adjust_wrap_length, update_line_numbers

logger = logging.getLogger(__name__)


class TemplateInput(ft.Container):  # type: ignore[misc]
    """Input component for DNA template sequence."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        input_data: GUIInput,
        on_change_handler: Callable[[ft.Event | None], None],
        handle_field_focus: Callable[[ft.Event[ft.TextField]], None],
        handle_field_blur: Callable[[ft.Event[ft.TextField]], None],
        handle_field_submit: Callable[[ft.Event[ft.TextField]], None],
        clear_template_callback: Callable[[ft.Event | None], None],
    ) -> None:
        """Initialise the TemplateInput component."""
        super().__init__(expand=5)
        self.app_page = page
        self.settings = settings
        self.input_data = input_data
        self.on_change_handler = on_change_handler

        # Line numbers / character count gutter
        self.line_numbers_text = ft.TextField(
            value="0",
            dense=True,
            multiline=True,
            read_only=True,
            text_align=ft.TextAlign.RIGHT,
            border=ft.InputBorder.NONE,
            bgcolor=ft.Colors.TRANSPARENT,
            hover_color=ft.Colors.TRANSPARENT,
            focused_bgcolor=ft.Colors.TRANSPARENT,
            content_padding=ft.Padding(0, 10, 0, 10),
            can_request_focus=False,
        )
        self.line_numbers_container = ft.Container(
            content=self.line_numbers_text,
            bgcolor=GUIColours.GUTTER_BG,
            padding=0,
            alignment=ft.Alignment(1, -1),
        )

        self.on_change_handler = on_change_handler
        self.handle_field_focus = handle_field_focus
        self.handle_field_blur = handle_field_blur
        self._is_focused = False
        self._cleaned_len = 0
        self._last_left_width: float | None = None
        self._last_active_selection: ft.TextSelection | None = None
        self._current_selected_cleaned_text: str = ""

        self.template_sequence = ft.TextField(
            dense=True,
            multiline=True,
            expand=False,
            hint_text="Enter DNA sequence here...",
            border=ft.InputBorder.NONE,
            bgcolor=ft.Colors.TRANSPARENT,
            hover_color=ft.Colors.TRANSPARENT,
            focused_bgcolor=ft.Colors.TRANSPARENT,
            content_padding=ft.Padding(0, 10, 10, 10),
            on_change=self._handle_change,
            on_focus=self._handle_focus,
            on_blur=self._handle_blur,
            on_submit=handle_field_submit,
            on_selection_change=self._handle_selection_change,
        )

        self.status_text = ft.Text(
            value="Insertion Point After Base: 0",
            size=12,
        )

        self.copy_selection_floating_button = ft.GestureDetector(
            content=ft.Container(
                content=ft.Row(
                    [
                        ft.Icon(
                            ft.Icons.COPY,
                            size=14,
                            color=GUIColours.WHITE,
                        ),
                        ft.Text(
                            "Copy Selection",
                            size=12,
                            weight=ft.FontWeight.BOLD,
                            color=GUIColours.WHITE,
                        ),
                    ],
                    spacing=5,
                    alignment=ft.MainAxisAlignment.CENTER,
                    tight=True,
                ),
                padding=ft.Padding(8, 4, 10, 4),
                bgcolor=GUIColours.PRIMARY,
                border_radius=5,
                shadow=ft.BoxShadow(
                    spread_radius=1, blur_radius=4, color="#60000000"
                ),
            ),
            on_tap=self._on_copy_selection_click,
            mouse_cursor=ft.MouseCursor.CLICK,
            top=10,
            right=20,
            visible=False,
        )
        self.copy_selection_button = self.copy_selection_floating_button

        self.bases_per_line_label = ft.Text(
            value="Bases per line:",
            size=12,
        )
        self.bases_per_line_value_text = ft.Text(
            value=str(self.settings.get("template_bases_per_line", "Auto")),
            size=12,
        )
        self.bases_per_line_value_container = ft.Container(
            content=self.bases_per_line_value_text,
            padding=ft.Padding(5, 2, 5, 2),
            border=ft.Border.all(1, GUIColours.OUTLINE),
            border_radius=3,
            alignment=ft.Alignment(0, 0),
        )
        self.bases_per_line_menu = ft.PopupMenuButton(
            content=self.bases_per_line_value_container,
            items=[
                ft.PopupMenuItem(
                    content="Auto",
                    on_click=lambda e: self._handle_menu_select("Auto"),
                )
            ]
            + [
                ft.PopupMenuItem(
                    content=str(val),
                    on_click=lambda e, v=val: self._handle_menu_select(v),
                )
                for val in range(10, 110, 10)
            ],
            padding=0,
        )
        self.bases_per_line_container = ft.Container(
            content=ft.Row(
                [
                    self.bases_per_line_label,
                    self.bases_per_line_menu,
                ],
                spacing=5,
                alignment=ft.MainAxisAlignment.START,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
            ),
            height=24,
            alignment=ft.Alignment(0, 0),
        )

        self.status_bar = ft.Container(
            content=ft.Row(
                [
                    self.status_text,
                    ft.Container(expand=True),
                    self.bases_per_line_container,
                ],
                alignment=ft.MainAxisAlignment.START,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
                height=28,
                spacing=5,
            ),
            bgcolor=GUIColours.INFO_HEADER_BG,
            padding=ft.Padding(10, 0, 10, 0),
            height=28,
        )

        self.template_sequence_wrapper = self.template_sequence

        self.template_sequence_container = ft.Row(
            [self.template_sequence_wrapper],
            spacing=0,
            vertical_alignment=ft.CrossAxisAlignment.START,
            expand=True,
            scroll=ft.Scrollbar(
                orientation=ft.ScrollbarOrientation.TOP,
                thumb_visibility=True,
            ),
        )

        self.sequence_layout = ft.Row(
            [
                self.line_numbers_container,
                self.template_sequence_container,
            ],
            spacing=0,
            vertical_alignment=ft.CrossAxisAlignment.START,
        )

        self.sequence_layout_stack = ft.Stack(
            [
                self.sequence_layout,
                self.copy_selection_floating_button,
            ],
            expand=True,
        )

        self.template_circular = ft.Checkbox(
            label="Circular",
            value=False,
            on_change=on_change_handler,
        )
        self.circular_container = ft.Container(
            content=self.template_circular,
            border=ft.Border.all(1, GUIColours.OUTLINE),
            border_radius=5,
            padding=ft.Padding(0, 0, 0, 0),
            height=32,
            width=115,
        )

        self.save_template_button = ft.FilledTonalButton(
            "Save",
            icon=ft.Icons.FILE_DOWNLOAD,
            on_click=self._save_template_click,
            tooltip="Save template to TXT",
            height=32,
        )

        self.load_template_button = ft.FilledTonalButton(
            "Load",
            icon=ft.Icons.FILE_OPEN,
            on_click=self._load_template_click,
            tooltip="Load template from TXT",
            height=32,
        )

        self.upper_case_button = ft.OutlinedButton(
            "Upper",
            icon=ft.Icons.TEXT_FIELDS,
            tooltip="Convert selected bases to upper case",
            on_click=self._upper_case_click,
            height=32,
        )

        self.lower_case_button = ft.OutlinedButton(
            "Lower",
            icon=ft.Icons.TEXT_FIELDS,
            tooltip="Convert selected bases to lower case",
            on_click=self._lower_case_click,
            height=32,
        )

        self.casing_group = ft.Row(
            [self.upper_case_button, self.lower_case_button],
            spacing=10,
            tight=True,
        )

        self.copy_template_button = ft.OutlinedButton(
            "Copy",
            icon=ft.Icons.COPY,
            tooltip="Copy template sequence without linebreaks",
            on_click=self._on_copy_click,
            height=32,
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
                ft.ResponsiveRow(
                    [
                        ft.Container(
                            content=ft.Row(
                                [
                                    self.circular_container,
                                    ft.Text(
                                        "Template Sequence",
                                        weight=ft.FontWeight.BOLD,
                                        no_wrap=True,
                                    ),
                                ],
                                spacing=10,
                                tight=True,
                                wrap=True,
                            ),
                            col={
                                "xl": 4,
                                "lg": 4,
                                "md": 12,
                                "sm": 12,
                                "xs": 12,
                            },
                        ),
                        ft.Container(
                            content=ft.Row(
                                [
                                    self.load_template_button,
                                    self.save_template_button,
                                    self.clear_template_button,
                                    self.copy_template_button,
                                    self.casing_group,
                                ],
                                spacing=10,
                                tight=True,
                                wrap=True,
                            ),
                            col={
                                "xl": 8,
                                "lg": 8,
                                "md": 12,
                                "sm": 12,
                                "xs": 12,
                            },
                            alignment=ft.Alignment(1, 0),
                        ),
                    ],
                ),
                ft.Stack(
                    [
                        ft.Container(
                            content=ft.Column(
                                [
                                    ft.ListView(
                                        [self.sequence_layout],
                                        expand=True,
                                        scroll=ft.ScrollMode.ALWAYS,
                                    ),
                                    self.status_bar,
                                ],
                                spacing=0,
                                horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                            ),
                            expand=True,
                            border=ft.Border.all(1, GUIColours.OUTLINE),
                            border_radius=5,
                            padding=0,
                        ),
                        self.copy_selection_floating_button,
                    ],
                    expand=True,
                ),
            ],
            expand=True,
            spacing=5,
        )
        self.update_ui()

    # ------------------------------------------------------------------ #
    # File management                                                      #
    # ------------------------------------------------------------------ #

    async def _load_template_click(self, e: ft.Event) -> None:
        """Open file picker to load template sequence from a TXT file."""
        await load_template_click(self, e)

    async def _save_template_click(self, e: ft.Event) -> None:
        """Save template sequence to a TXT file."""
        await save_template_click(self, e)

    # ------------------------------------------------------------------ #
    # Case conversion & Copying                                            #
    # ------------------------------------------------------------------ #

    def _on_copy_click(self, e: ft.Event) -> None:
        """Copy selected template text (or full sequence) without linebreaks."""
        sel = self.template_sequence.selection
        raw_text = self.template_sequence.value or ""
        if sel is not None and sel.is_valid and sel.start != sel.end:
            selected_text = raw_text[sel.start : sel.end]
        else:
            selected_text = raw_text

        cleaned_text = clean_sequence(selected_text)
        if not cleaned_text:
            self._show_notification("Nothing to copy!")
            return

        if self.app_page:
            if getattr(self.app_page, "web", False) and hasattr(
                self.app_page, "run_javascript"
            ):
                import json

                escaped_text = json.dumps(cleaned_text)
                self.app_page.run_javascript(  # pyright: ignore[reportAttributeAccessIssue]
                    f"navigator.clipboard.writeText({escaped_text});"
                )
            else:
                import pyperclip

                pyperclip.copy(cleaned_text)

        self._show_notification("Copied to clipboard!")

    def _on_copy_selection_click(self, e: ft.Event) -> None:
        """Copy selected template sequence without linebreaks."""
        text_to_copy = self._current_selected_cleaned_text

        if not text_to_copy:
            sel = self.template_sequence.selection
            if (
                sel is None
                or not getattr(sel, "is_valid", False)
                or getattr(sel, "base_offset", -1) < 0
                or getattr(sel, "base_offset", None)
                == getattr(sel, "extent_offset", None)
            ):
                sel = self._last_active_selection

            raw_text = self.template_sequence.value or ""
            if sel is not None and getattr(sel, "is_valid", False):
                base = getattr(sel, "base_offset", None)
                extent = getattr(sel, "extent_offset", None)
                if (
                    base is not None
                    and extent is not None
                    and base >= 0
                    and extent >= 0
                    and base != extent
                ):
                    start_idx = min(base, extent)
                    end_idx = max(base, extent)
                    text_to_copy = clean_sequence(raw_text[start_idx:end_idx])
                elif (
                    hasattr(sel, "start")
                    and hasattr(sel, "end")
                    and sel.start is not None
                    and sel.end is not None
                    and sel.start != sel.end
                ):
                    start_idx = min(sel.start, sel.end)
                    end_idx = max(sel.start, sel.end)
                    text_to_copy = clean_sequence(raw_text[start_idx:end_idx])
                else:
                    text_to_copy = clean_sequence(raw_text)
            else:
                text_to_copy = clean_sequence(raw_text)

        if not text_to_copy:
            self._show_notification("Nothing to copy!")
            return

        if self.app_page:
            try:
                self.app_page.set_clipboard(text_to_copy)
            except Exception as ex:
                logger.debug("Flet set_clipboard fallback error: %s", ex)
            if getattr(self.app_page, "web", False) and hasattr(
                self.app_page, "run_javascript"
            ):
                import json

                escaped_text = json.dumps(text_to_copy)
                self.app_page.run_javascript(
                    f"navigator.clipboard.writeText({escaped_text});"
                )
            else:
                try:
                    import pyperclip

                    pyperclip.copy(text_to_copy)
                except Exception as ex:
                    logger.debug("pyperclip copy fallback error: %s", ex)

        self._show_notification("Copied selection to clipboard!")

    def _upper_case_click(self, e: ft.Event) -> None:
        """Handle upper case button click."""
        change_selection_case(self, to_upper=True)

    def _lower_case_click(self, e: ft.Event) -> None:
        """Handle lower case button click."""
        change_selection_case(self, to_upper=False)

    # ------------------------------------------------------------------ #
    # Formatting / gutter                                                  #
    # ------------------------------------------------------------------ #

    @property
    def current_left_width(self) -> float:
        """Get the current available width of the left panel."""
        if isinstance(self.width, int | float) and self.width > 0:
            return float(self.width)
        right_fraction = 0.5
        if self.parent:
            right_fraction = getattr(self.parent, "right_fraction", 0.5)
        page_width = (
            self.app_page.width
            if (self.app_page and self.app_page.width)
            else 800
        )
        return float(page_width * (1.0 - right_fraction))

    def adjust_wrap_length(self, left_width: float, update: bool = True) -> int:
        """Adjust wrap length based on selected width or Auto."""
        return adjust_wrap_length(self, left_width, update)

    def _update_line_numbers(
        self, update: bool = True, gutter_only: bool = False
    ) -> None:
        """Update the line numbers gutter based on current template sequence."""
        update_line_numbers(self, update=update, gutter_only=gutter_only)

    def _validate_bases_per_line(
        self, val_str: str | None = None
    ) -> int | str | None:
        """Validate bases per line, enforcing 10..100 or Auto."""
        from .formatter import validate_bases_per_line

        return validate_bases_per_line(self, val_str)

    def _handle_menu_select(self, val: int | str) -> None:
        """Handle selection of bases per line from popup menu."""
        self.settings["template_bases_per_line"] = val
        self.settings.save_to_local(self.app_page)
        self.bases_per_line_value_text.value = str(val)
        self.adjust_wrap_length(self.current_left_width)

    # ------------------------------------------------------------------ #
    # Notification                                                         #
    # ------------------------------------------------------------------ #

    def _show_notification(self, message: str) -> None:
        """Show a notification message."""
        if not hasattr(self, "_notification_helper"):
            self._notification_helper = NotificationHelper(self.app_page)
        self._notification_helper.show_message(message)

    # ------------------------------------------------------------------ #
    # State synchronisation                                                #
    # ------------------------------------------------------------------ #

    def sync_to_state(self) -> None:
        """Sync template text field to the central state."""
        self.template_sequence.value = self.template_sequence.value or ""
        self.input_data.template = clean_sequence(
            str(self.template_sequence.value)
        )
        self._cleaned_len = len(self.input_data.template)
        self.input_data.template_circular = bool(self.template_circular.value)
        self._update_line_numbers(update=False)
        self._update_status_bar(update=False)

    def _sync_dark_mode(self) -> None:
        """Synchronise GUIColours dark_mode state with page settings."""
        dark_setting = str(self.settings.get("dark_mode", False)).lower()
        if dark_setting == "system":
            brightness = str(
                getattr(self.app_page, "platform_brightness", "")
            ).lower()
            GUIColours.dark_mode = brightness == "dark"
        else:
            GUIColours.dark_mode = dark_setting not in (
                "false",
                "0",
                "no",
                "none",
                "",
            )

    def update_ui(self) -> None:
        """Update template controls, line numbers gutter, and status bar."""
        self._sync_dark_mode()
        font_size = self.settings.get("font_size_default", 14)
        font_family = self.settings.get("font_family", "Roboto Mono")

        self.template_sequence.text_style = ft.TextStyle(
            font_family=font_family,
            size=font_size,
            height=1.5,
        )
        self.line_numbers_text.text_style = ft.TextStyle(
            font_family=font_family,
            color=GUIColours.TEXT_ON_SURFACE,
            size=font_size,
            height=1.5,
        )
        self.line_numbers_container.bgcolor = GUIColours.GUTTER_BG

        self.template_sequence.value = self.input_data.template
        self._cleaned_len = len(self.input_data.template)
        self.template_circular.value = self.input_data.template_circular

        # Update status bar style
        self.status_text.style = ft.TextStyle(
            font_family=font_family,
            color=GUIColours.TEXT_ON_SURFACE,
            size=12,
        )
        self.bases_per_line_label.style = ft.TextStyle(
            font_family=font_family,
            color=GUIColours.TEXT_ON_SURFACE,
            size=12,
        )
        self.bases_per_line_value_text.value = str(
            self.settings.get("template_bases_per_line", "Auto")
        )
        self.bases_per_line_value_text.style = ft.TextStyle(
            font_family=font_family,
            color=GUIColours.TEXT_ON_SURFACE,
            size=12,
        )

        self.status_bar.bgcolor = GUIColours.INFO_HEADER_BG
        self.status_bar.border = ft.Border(
            top=ft.BorderSide(1, GUIColours.OUTLINE)
        )

        self._update_line_numbers(update=False)
        self._update_status_bar(None, update=False)

    # ------------------------------------------------------------------ #
    # Event handlers                                                       #
    # ------------------------------------------------------------------ #

    def _handle_change(self, e: ft.ControlEvent) -> None:
        """Handle template text changes, updating gutter line numbers."""
        sel = self.template_sequence.selection
        clean_idx = 0
        has_selection = sel is not None and sel.is_valid
        if sel is not None and sel.is_valid:
            raw_text = self.template_sequence.value or ""
            prefix = raw_text[: sel.base_offset]
            clean_idx = len(clean_sequence(prefix))

        raw_val = self.template_sequence.value or ""
        cleaned = clean_sequence(raw_val)
        self.input_data.template = cleaned
        self._cleaned_len = len(cleaned)

        left_width = self.current_left_width
        wrap_length = self.adjust_wrap_length(left_width, update=False)

        if has_selection:
            if clean_idx <= 0:
                new_pos = 0
            elif clean_idx >= len(cleaned):
                num_newlines = max(0, (len(cleaned) - 1) // wrap_length)
                new_pos = clean_idx + num_newlines
            else:
                num_newlines = clean_idx // wrap_length
                new_pos = clean_idx + num_newlines

            self.template_sequence.selection = ft.TextSelection(
                base_offset=new_pos, extent_offset=new_pos
            )

        try:
            self.update()
        except (RuntimeError, AssertionError):
            pass

        self.on_change_handler(e)

    def _handle_focus(self, e: ft.Event[ft.TextField]) -> None:
        """Handle template focus, updating status bar text selection info."""
        self._is_focused = True
        self._update_status_bar(self.template_sequence.selection)
        self.handle_field_focus(e)

    def _handle_blur(self, e: ft.Event[ft.TextField]) -> None:
        """Handle template focus loss, resetting status bar selection info."""
        self._is_focused = False
        self._update_status_bar(None)
        self.handle_field_blur(e)

    def _handle_selection_change(self, e: ft.ControlEvent) -> None:
        """Handle selection change event on the template text field."""
        sel = getattr(e, "selection", None) or self.template_sequence.selection
        if sel is not None and getattr(sel, "is_valid", False):
            base = getattr(sel, "base_offset", None)
            extent = getattr(sel, "extent_offset", None)
            if base is None or extent is None:
                base = getattr(sel, "start", None)
                extent = getattr(sel, "end", None)

            if (
                base is not None
                and extent is not None
                and base >= 0
                and extent >= 0
                and base != extent
            ):
                self._last_active_selection = sel
                start_idx = min(base, extent)
                end_idx = max(base, extent)
                raw_text = self.template_sequence.value or ""
                cleaned = clean_sequence(raw_text[start_idx:end_idx])
                if cleaned:
                    self._current_selected_cleaned_text = cleaned
        self._update_status_bar(sel)

    def _count_bases(self, prefix: str) -> int:
        """Count the number of biological bases in a prefix string."""
        return len(clean_sequence(prefix))

    def _update_status_bar(
        self,
        selection: ft.TextSelection | None = None,
        update: bool = True,
    ) -> None:
        """Update the status bar text based on text selection and focus."""
        sel = (
            selection
            if (selection is not None and selection.is_valid)
            else self.template_sequence.selection
        )
        if (sel is None or not sel.is_valid or sel.start == sel.end) and (
            self._last_active_selection is not None
            and self._last_active_selection.is_valid
            and self._last_active_selection.start
            != self._last_active_selection.end
        ):
            sel = self._last_active_selection

        has_selection = (
            sel is not None and sel.is_valid and sel.start != sel.end
        ) or bool(self._current_selected_cleaned_text)

        if not self._is_focused:
            self.status_text.value = f"Total Bases: {self._cleaned_len}"
        else:
            if (
                has_selection
                and sel is not None
                and sel.is_valid
                and sel.start != sel.end
            ):
                raw_text = self.template_sequence.value or ""
                bases_before = self._count_bases(raw_text[: sel.start])
                bases_total = self._count_bases(raw_text[: sel.end])

                if bases_total > bases_before:
                    self.status_text.value = (
                        f"Selected Bases: {bases_before + 1} - {bases_total}"
                    )
                else:
                    self.status_text.value = (
                        f"Insertion Point After Base: {bases_before}"
                    )
            elif sel is not None and sel.is_valid:
                raw_text = self.template_sequence.value or ""
                bases_before = self._count_bases(raw_text[: sel.start])
                self.status_text.value = (
                    f"Insertion Point After Base: {bases_before}"
                )
            else:
                self.status_text.value = (
                    f"Insertion Point After Base: {self._cleaned_len}"
                )

        if (
            has_selection
            and sel is not None
            and sel.is_valid
            and sel.start != sel.end
        ):
            raw_text = self.template_sequence.value or ""
            text_up_to_end = raw_text[: sel.end]
            lines = text_up_to_end.split("\n")
            line_idx = max(0, len(lines) - 1)
            col_idx = len(lines[-1]) if lines else 0

            font_size = max(1, self.settings.get("font_size_default", 14))
            line_height = font_size * 1.5
            char_width = font_size * 0.70

            max_digits = len(str(max(1, len(self.input_data.template))))
            gutter_width = 20 + max_digits * char_width

            top_pos = 5.0 + (line_idx * line_height)
            left_pos = gutter_width + (col_idx * char_width) + 10.0

            max_left = max(gutter_width + 50.0, self.current_left_width - 150.0)
            left_pos = min(left_pos, max_left)

            self.copy_selection_floating_button.top = top_pos
            self.copy_selection_floating_button.left = left_pos
            self.copy_selection_floating_button.right = None
            self.copy_selection_button.visible = True
        elif has_selection:
            self.copy_selection_button.visible = True
        else:
            self.copy_selection_button.visible = False

        if update:
            try:
                if self.copy_selection_floating_button.page:
                    self.copy_selection_floating_button.update()
            except (RuntimeError, AssertionError):
                pass
            try:
                if self.sequence_layout_stack.page:
                    self.sequence_layout_stack.update()
            except (RuntimeError, AssertionError):
                pass
            try:
                page = self.status_bar.page
            except RuntimeError:
                page = None
            if page:
                try:
                    self.status_bar.update()
                except (RuntimeError, AssertionError):
                    pass

        if update:
            try:
                if self.copy_selection_floating_button.page:
                    self.copy_selection_floating_button.update()
            except (RuntimeError, AssertionError):
                pass
            try:
                if self.sequence_layout_stack.page:
                    self.sequence_layout_stack.update()
            except (RuntimeError, AssertionError):
                pass
            try:
                page = self.status_bar.page
            except RuntimeError:
                page = None
            if page:
                try:
                    self.status_bar.update()
                except (RuntimeError, AssertionError):
                    pass
