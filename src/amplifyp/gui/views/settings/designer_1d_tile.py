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

"""Designer1DTile expansion tile component for settings view."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.gui_helpers import BorderedCheckbox


class Designer1DTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Designer 1D settings."""

    def __init__(
        self,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.ControlEvent], None],
        header_size: int,
    ) -> None:
        """Initialise the Designer1DTile."""
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        self.show_tm_checkbox = BorderedCheckbox(
            label="Show melting temperature (Tm) on cards",
            value=self.settings.get("designer_1d_show_tm", True),
            on_change=self.on_change_handler,
        )

        self.show_pct_at_checkbox = BorderedCheckbox(
            label="Show % AT on cards",
            value=self.settings.get("designer_1d_show_pct_at", False),
            on_change=self.on_change_handler,
        )

        self.settings_map["designer_1d_show_tm"] = self.show_tm_checkbox
        self.settings_map["designer_1d_show_pct_at"] = self.show_pct_at_checkbox

        super().__init__(
            title=ft.Text(
                "Designer 1D",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Column(
                        [
                            ft.Row(
                                [
                                    ft.Container(
                                        content=ft.Column(
                                            [
                                                self.show_tm_checkbox,
                                                self.show_pct_at_checkbox,
                                            ],
                                            spacing=15,
                                            horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                        ),
                                        width=700,
                                    ),
                                ],
                                alignment=ft.MainAxisAlignment.CENTER,
                            ),
                        ],
                        spacing=15,
                        horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )

    @property
    def set_designer_1d_show_tm(self) -> BorderedCheckbox:
        """Get the show melting temperature checkbox."""
        return self.show_tm_checkbox

    @property
    def set_designer_1d_show_pct_at(self) -> BorderedCheckbox:
        """Get the show % AT checkbox."""
        return self.show_pct_at_checkbox
