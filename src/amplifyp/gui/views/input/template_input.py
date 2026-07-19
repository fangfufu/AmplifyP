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

from collections.abc import Callable

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.sequence import clean_sequence, format_sequence
from amplifyp.gui.utils.ui import NotificationHelper


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

        self.bases_per_line_label = ft.Text(
            value="Bases per line:",
            size=12,
        )
        self.bases_per_line_value_text = ft.Text(
            value=str(self.settings.get("template_bases_per_line", 50)),
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
            padding=ft.Padding(10, 0, 10, 0),
            height=28,
        )

        self.template_sequence_wrapper = ft.Container(
            content=self.template_sequence,
            padding=0,
            margin=0,
        )

        self.template_sequence_container = ft.Row(
            [self.template_sequence_wrapper],
            spacing=0,
            vertical_alignment=ft.CrossAxisAlignment.START,
            expand=True,
            scroll=ft.Scrollbar(
                orientation=ft.ScrollbarOrientation.TOP,
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
            ],
            expand=True,
            spacing=5,
        )
        self.update_ui()

    async def _load_template_click(self, e: ft.Event) -> None:
        """Open file picker to load template sequence from a TXT file.

        Args:
            e: The Flet control event triggered by the load button click.
        """
        from amplifyp.gui.utils.io import pick_and_read_file

        content = await pick_and_read_file(
            page=self.app_page,
            dialog_title="Load",
            allowed_extensions=["txt"],
            show_notification=self._show_notification,
        )
        if content is None:
            return

        self.input_data.template = content
        self.update_ui()
        self.on_change_handler(None)
        self._show_notification("Template loaded successfully.")

    async def _save_template_click(self, e: ft.Event) -> None:
        """Save template sequence to a TXT file.

        Args:
            e: The Flet control event triggered by the save button click.
        """
        template_content = self.input_data.template
        if not template_content.strip():
            self._show_notification("No template to save.")
            return

        from amplifyp.gui.utils.io import save_and_write_file

        await save_and_write_file(
            page=self.app_page,
            dialog_title="Save",
            file_name="template.txt",
            allowed_extensions=["txt"],
            content=template_content,
            show_notification=self._show_notification,
            success_message_desktop="Template saved successfully.",
            success_message_web="Template ready for download!",
        )

    def _upper_case_click(self, e: ft.Event) -> None:
        """Handle upper case button click."""
        self._change_selection_case(to_upper=True)

    def _lower_case_click(self, e: ft.Event) -> None:
        """Handle lower case button click."""
        self._change_selection_case(to_upper=False)

    def _change_selection_case(self, to_upper: bool) -> None:
        """Convert selected bases in template to upper/lower case."""
        sel = self.template_sequence.selection
        if not sel or not sel.is_valid or sel.start == sel.end:
            self._show_notification("Please select sequence text first.")
            return

        raw_val = self.template_sequence.value or ""
        start, end = sel.start, sel.end

        selected_text = raw_val[start:end]
        modified_text = (
            selected_text.upper() if to_upper else selected_text.lower()
        )
        new_val = raw_val[:start] + modified_text + raw_val[end:]

        self.template_sequence.value = new_val
        self.input_data.template = clean_sequence(new_val)
        self._cleaned_len = len(self.input_data.template)

        self._update_line_numbers(update=False)
        self.template_sequence.selection = ft.TextSelection(
            base_offset=sel.base_offset, extent_offset=sel.extent_offset
        )
        try:
            self.template_sequence.update()
            self.update()
        except (RuntimeError, AssertionError):
            pass
        self.on_change_handler(None)

    def _show_notification(self, message: str) -> None:
        """Show a notification message.

        Args:
            message: The message to display in the notification.
        """
        if not hasattr(self, "_notification_helper"):
            self._notification_helper = NotificationHelper(self.app_page)
        self._notification_helper.show_message(message)

    def sync_to_state(self) -> None:
        """Sync template text field to the central state.

        Reads the current UI values and writes them into the central
        ``GUIInput`` state object.
        """
        self.template_sequence.value = self.template_sequence.value or ""
        self.input_data.template = clean_sequence(
            str(self.template_sequence.value)
        )
        self._cleaned_len = len(self.input_data.template)
        self.input_data.template_circular = bool(self.template_circular.value)
        self._update_line_numbers(update=False)
        self._update_status_bar(update=False)

    def update_ui(self) -> None:
        """Update template UI elements to match central state.

        Applies values from the central ``GUIInput`` state object to the
        template text field and circular checkbox controls.
        """
        font_family = self.settings.get("font_family", "Roboto Mono")
        font_size = self.settings.get("font_size_default", 14)

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
            self.settings.get("template_bases_per_line", 50)
        )
        self.bases_per_line_value_text.style = ft.TextStyle(
            font_family=font_family,
            color=GUIColours.TEXT_ON_SURFACE,
            size=12,
        )

        self.status_bar.bgcolor = GUIColours.GUTTER_BG
        self.status_bar.border = ft.Border(
            top=ft.BorderSide(1, GUIColours.OUTLINE)
        )

        self._update_line_numbers(update=False)
        self._update_status_bar(None, update=False)

    @property
    def current_left_width(self) -> float:
        """Get the current available width of the left panel."""
        if isinstance(self.width, (int, float)) and self.width > 0:
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
        self._last_left_width = left_width
        self.sequence_layout.width = max(100.0, left_width - 15.0)

        font_size = max(1, self.settings.get("font_size_default", 14))
        char_width = font_size * 0.70

        # Calculate dynamic gutter width based on template digits
        template_len = len(self.input_data.template)
        max_digits = len(str(max(1, template_len)))
        gutter_width = 20 + max_digits * char_width

        wrap_setting = self._validate_bases_per_line()
        if wrap_setting == "Auto":
            available_width = left_width - gutter_width - 100
            max_fit = int(available_width / char_width)
            wrap_length = (max_fit // 10) * 10
            wrap_length = max(10, min(100, wrap_length))
        else:
            wrap_length = wrap_setting if isinstance(wrap_setting, int) else 50

        field_available_width = max(100.0, left_width - gutter_width - 35.0)

        # Disable autowrapping at window edge, enable horizontal scroll
        target_width = max(
            field_available_width + 100.0, wrap_length * char_width + 100.0
        )
        self.template_sequence_wrapper.width = target_width
        self.template_sequence.width = target_width
        self.template_sequence.expand = False

        # Update TextField content with new wrapping
        self.template_sequence.value = format_sequence(
            self.input_data.template, wrap_length
        )
        self._update_line_numbers(update=update)

        try:
            if self.page:
                self.template_sequence_wrapper.update()
                self.template_sequence.update()
                self.template_sequence_container.update()
        except (AssertionError, RuntimeError):
            pass

        return wrap_length

    def _update_line_numbers(
        self, update: bool = True, gutter_only: bool = False
    ) -> None:
        """Update the line numbers gutter based on current template sequence."""
        text = self.template_sequence.value or ""
        lines = text.split("\n")
        line_indices = []
        current_idx = 0
        for line in lines:
            line_indices.append(str(current_idx))
            current_idx += len(clean_sequence(line))

        self.line_numbers_text.value = "\n".join(line_indices)

        # Set dynamic gutter width
        font_size = max(1, self.settings.get("font_size_default", 14))
        char_width = font_size * 0.75
        # Calculate from live text excluding whitespaces instead of stored state
        template_len = current_idx
        max_digits = len(str(max(1, template_len)))
        gutter_width = 20 + max_digits * char_width
        self.line_numbers_container.width = gutter_width

        if update:
            try:
                page = self.page
            except RuntimeError:
                page = None
            if page:
                try:
                    if gutter_only:
                        self.line_numbers_text.update()
                        self.line_numbers_container.update()
                    else:
                        self.update()
                except (RuntimeError, AssertionError):
                    pass

    def _validate_bases_per_line(
        self, val_str: str | None = None
    ) -> int | str | None:
        """Validate bases per line, enforcing 10..100 or Auto.

        Args:
            val_str: Optional string value to validate. If None, reads from
                the current UI selection.

        Returns:
            The validated integer value, 'Auto', or None if invalid.
        """
        if val_str is None:
            val_str = (self.bases_per_line_value_text.value or "").strip()
        if val_str.lower() == "auto":
            return "Auto"
        try:
            val_int = int(val_str.strip())
            if 10 <= val_int <= 100 and val_int % 10 == 0:
                return val_int
        except ValueError:
            pass
        return None

    def _handle_menu_select(self, val: int | str) -> None:
        """Handle selection of bases per line from popup menu."""
        self.settings["template_bases_per_line"] = val
        self.settings.save_to_local(self.app_page)
        self.bases_per_line_value_text.value = str(val)
        self.adjust_wrap_length(self.current_left_width)

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
        self._update_status_bar(getattr(e, "selection", None))

    def _count_bases(self, prefix: str) -> int:
        """Count the number of biological bases in a prefix string."""
        return len(clean_sequence(prefix))

    def _update_status_bar(
        self, selection: ft.TextSelection | None = None, update: bool = True
    ) -> None:
        """Update the status bar text based on text selection and focus."""
        if not self._is_focused:
            self.status_text.value = f"Total Bases: {self._cleaned_len}"
        else:
            sel = (
                selection
                if (selection is not None and selection.is_valid)
                else self.template_sequence.selection
            )
            if sel is not None and sel.is_valid:
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
            else:
                self.status_text.value = (
                    f"Insertion Point After Base: {self._cleaned_len}"
                )

        if update:
            try:
                page = self.status_bar.page
            except RuntimeError:
                page = None
            if page:
                try:
                    self.status_bar.update()
                except (RuntimeError, AssertionError):
                    pass
