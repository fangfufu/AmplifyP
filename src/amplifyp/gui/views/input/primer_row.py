# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""A single row representing a primer in the list."""

from typing import Any

import flet as ft

from amplifyp.gui.settings import GUIColors


class PrimerRow(ft.Container):  # type: ignore[misc]
    """A single row representing a primer in the list."""

    def __init__(
        self,
        idx: int,
        name: str,
        seq: str,
        is_active: bool,
        is_dup: bool,
        name_error: str | None,
        seq_error: str | None,
        font_family: str,
        name_column_width: float,
        on_change_handler: Any,
        handle_field_focus: Any,
        handle_field_blur: Any,
        handle_field_submit: Any,
        on_row_click: Any,
        on_move_primer: Any,
        on_delete_primer: Any,
        on_divider_pan: Any,
        on_divider_pan_end: Any,
        is_focused: bool,
        is_last_row: bool,
        is_penultimate_row: bool,
    ) -> None:
        """Initialize the PrimerRow."""
        has_err = bool(name_error or seq_error)
        super().__init__(
            data=idx,
            bgcolor=GUIColors.DUPLICATE_BG if is_dup else None,
            padding=0,
            height=30 if not has_err else None,
        )
        self.idx = idx
        self.checkbox = ft.Checkbox(
            value=is_active if not has_err else False,
            on_change=on_change_handler,
            disabled=has_err,
            visible=not is_last_row,
        )
        self.checkbox_container = ft.Container(
            content=self.checkbox,
            width=55,
            height=30,
            alignment=ft.Alignment(0, 0),
        )
        self.name_field = ft.TextField(
            value=name,
            hint_text="New Primer Name",
            dense=True,
            content_padding=ft.Padding(5, 0, 0, 0),
            height=30 if not name_error else None,
            width=name_column_width,
            border=ft.InputBorder.NONE,
            data=idx,
            on_focus=handle_field_focus,
            on_blur=handle_field_blur,
            on_submit=handle_field_submit,
        )
        self.seq_field = ft.TextField(
            value=seq,
            hint_text="New Primer Sequence",
            dense=True,
            content_padding=ft.Padding(5, 0, 5, 0),
            height=30 if not seq_error else None,
            border=ft.InputBorder.NONE,
            text_style=ft.TextStyle(font_family=font_family),
            data=idx,
            expand=True,
            on_focus=handle_field_focus,
            on_blur=handle_field_blur,
            on_submit=handle_field_submit,
        )
        if name_error:
            self.name_field.error = name_error
        if seq_error:
            self.seq_field.error = seq_error

        self.divider = ft.GestureDetector(
            on_pan_update=on_divider_pan,
            on_pan_end=on_divider_pan_end,
            content=ft.Container(
                width=4,
                bgcolor=GUIColors.DIVIDER_GREY,
                margin=0,
                height=30,
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
        )

        self.active_divider = ft.Container(
            width=4,
            bgcolor=GUIColors.DIVIDER_GREY,
            margin=0,
            height=30,
        )

        self.reorder_controls = None
        self.control_container = None
        if not is_last_row:
            up_button = ft.IconButton(
                icon=ft.Icons.ARROW_UPWARD,
                icon_size=16,
                width=24,
                height=24,
                padding=0,
                tooltip="Move Up",
                disabled=(idx == 0),
                on_click=lambda e: on_move_primer(idx, -1),
            )
            down_button = ft.IconButton(
                icon=ft.Icons.ARROW_DOWNWARD,
                icon_size=16,
                width=24,
                height=24,
                padding=0,
                tooltip="Move Down",
                disabled=is_penultimate_row,
                on_click=lambda e: on_move_primer(idx, 1),
            )
            delete_button = ft.IconButton(
                icon=ft.Icons.DELETE_OUTLINE,
                icon_size=16,
                width=24,
                height=24,
                padding=0,
                tooltip="Delete Primer",
                on_click=lambda e: on_delete_primer(idx),
            )
            self.reorder_controls = ft.Row(
                [delete_button, up_button, down_button],
                spacing=2,
                alignment=ft.MainAxisAlignment.CENTER,
            )
            self.control_container = ft.Container(
                content=self.reorder_controls,
                width=82 if is_focused else 0,
                height=30,
                alignment=ft.Alignment(0, 0),
            )
            self.reorder_controls.visible = is_focused
        else:
            self.control_container = ft.Container(width=82 if is_focused else 0)

        self.content = ft.Row(
            [
                self.checkbox_container,
                self.active_divider,
                self.name_field,
                self.divider,
                self.seq_field,
                self.control_container,
            ],
            spacing=0,
            vertical_alignment=ft.CrossAxisAlignment.START,
        )
        self.on_click = lambda e: on_row_click(idx, self.name_field)

    def update_highlight_and_reorder(
        self, is_focused: bool, is_dup: bool
    ) -> None:
        """Update the background color and reorder buttons layout."""
        if is_focused:
            self.bgcolor = GUIColors.SELECTED_ROW_BG
        elif is_dup:
            self.bgcolor = GUIColors.DUPLICATE_BG
        else:
            self.bgcolor = None  # type: ignore[assignment]

        if self.control_container is not None:
            self.control_container.width = 82 if is_focused else 0
        if self.reorder_controls is not None:
            self.reorder_controls.visible = is_focused

    def set_error(self, err: dict[str, str | None] | str | None) -> None:
        """Set or clear the error message.

        Also adjusts the sequence field, name field and container height.
        """
        if isinstance(err, dict):
            name_error = err.get("name")
            seq_error = err.get("seq")
        else:
            # Compatibility: if string or None is passed
            if err == "Duplicate primer name":
                name_error = err
                seq_error = None
            else:
                name_error = None
                seq_error = err

        self.name_field.error = name_error
        self.name_field.height = 30 if not name_error else None
        self.seq_field.error = seq_error
        self.seq_field.height = 30 if not seq_error else None

        has_err = bool(name_error or seq_error)
        self.height = 30 if not has_err else None
        self.checkbox.disabled = has_err
        if has_err:
            self.checkbox.value = False

    @staticmethod
    def validate(name: str, seq: str) -> str | None:
        """Validate a primer sequence and name.

        Returns an error message if the primer is invalid.
        """
        if not seq:
            return None
        try:
            from amplifyp.dna import Primer

            Primer(sequence=seq, name=name)
            return None
        except ValueError as ex:
            return str(ex)
