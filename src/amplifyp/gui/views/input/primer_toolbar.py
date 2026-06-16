# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""Toolbar containing actions for managing primers (Save, Load, Clear)."""

from typing import Any

import flet as ft


class PrimerToolbar(ft.Row):  # type: ignore[misc]
    """Toolbar containing actions for managing primers (Save, Load, Clear)."""

    def __init__(
        self,
        on_save: Any,
        on_load: Any,
        on_clear: Any,
        on_delete_selected: Any,
    ) -> None:
        """Initialize the PrimerToolbar."""
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
            tooltip="Delete Selected Primers",
            on_click=on_delete_selected,
            height=32,
            disabled=True,
        )
        super().__init__(
            [
                self.load_button,
                self.save_button,
                self.delete_selected_button,
                self.clear_button,
            ],
            spacing=10,
            tight=True,
            wrap=True,
        )
