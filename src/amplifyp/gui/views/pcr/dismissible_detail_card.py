# Copyright (C) 2026 Fufu Fang
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

"""DismissibleDetailCard base component for PCR view detail cards."""

from collections.abc import Callable
from typing import Any

import flet as ft


class DismissibleDetailCard(ft.Card):  # type: ignore[misc]
    """Base card containing title, close button, and dismiss logic."""

    def __init__(
        self,
        card_id: str,
        title: str,
        settings: Any,
        dismiss_callback: Callable[[ft.Card], None],
        body_controls: list[ft.Control],
    ) -> None:
        """Initialise the DismissibleDetailCard."""
        super().__init__()
        self._card_id = card_id

        def remove_card(e: Any) -> None:
            dismiss_callback(self)

        self.content = ft.Container(
            padding=10,
            content=ft.Column(
                [
                    ft.Row(
                        [
                            ft.Text(
                                title,
                                weight=ft.FontWeight.BOLD,
                                size=settings.get("font_size_subheader", 16),
                                selectable=True,
                            ),
                            ft.IconButton(
                                icon=ft.Icons.CLOSE,
                                icon_size=18,
                                tooltip="Dismiss",
                                on_click=remove_card,
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    ),
                    *body_controls,
                ]
            ),
        )
