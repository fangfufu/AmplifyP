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

"""Toolbar component for PrimerInput containing managing actions."""

from __future__ import annotations

from collections.abc import Awaitable, Callable

import flet as ft


class PrimerToolbar(ft.Row):  # type: ignore[misc]
    """Toolbar containing actions for managing primers (Save, Load, Clear)."""

    def __init__(
        self,
        on_save: Callable[[ft.Event | None], None | Awaitable[None]],
        on_load: Callable[[ft.Event | None], None | Awaitable[None]],
        on_clear: Callable[[ft.Event | None], None | Awaitable[None]],
        on_delete_selected: Callable[[ft.Event | None], None | Awaitable[None]],
        on_copy: Callable[[ft.Event | None], None | Awaitable[None]],
        on_paste: Callable[[ft.Event | None], None | Awaitable[None]],
    ) -> None:
        """Initialise the PrimerToolbar."""
        self.save_button = ft.FilledTonalButton(
            "Save",
            icon=ft.Icons.FILE_DOWNLOAD,
            on_click=on_save,
            tooltip="Save primers to CSV",
            height=32,
        )
        self.load_button = ft.FilledTonalButton(
            "Load",
            icon=ft.Icons.FILE_OPEN,
            on_click=on_load,
            tooltip="Load primers from CSV/TSV",
            height=32,
        )
        self.clear_button = ft.OutlinedButton(
            "Clear",
            icon=ft.Icons.DELETE_OUTLINE,
            tooltip="Clear All Primers",
            on_click=on_clear,
            height=32,
        )
        self.delete_selected_button = ft.OutlinedButton(
            "Delete",
            icon=ft.Icons.DELETE_SWEEP,
            tooltip="Delete Highlighted Primers",
            on_click=on_delete_selected,
            height=32,
            disabled=True,
        )
        self.copy_button = ft.OutlinedButton(
            "Copy",
            icon=ft.Icons.CONTENT_COPY,
            tooltip="Copy Highlighted Primers to Clipboard",
            on_click=on_copy,
            height=32,
        )
        self.paste_button = ft.OutlinedButton(
            "Paste",
            icon=ft.Icons.CONTENT_PASTE,
            tooltip="Paste Primers from Clipboard",
            on_click=on_paste,
            height=32,
        )
        self.file_group = ft.Row(
            [self.load_button, self.save_button],
            spacing=10,
            tight=True,
        )
        self.clipboard_group = ft.Row(
            [self.copy_button, self.paste_button],
            spacing=10,
            tight=True,
        )
        self.edit_group = ft.Row(
            [self.delete_selected_button, self.clear_button],
            spacing=10,
            tight=True,
        )
        super().__init__(
            [
                self.file_group,
                self.clipboard_group,
                self.edit_group,
            ],
            spacing=10,
            alignment=ft.MainAxisAlignment.END,
            tight=True,
            wrap=True,
        )
