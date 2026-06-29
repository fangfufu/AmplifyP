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

import flet as ft

from amplifyp.gui.settings import GUISettings


class DismissibleDetailCard(ft.Card):  # type: ignore[misc]
    """Base card containing title, close button, and dismiss logic."""

    def __init__(
        self,
        card_id: str,
        title: str,
        settings: GUISettings,
        dismiss_callback: Callable[[ft.Card], None],
        body_controls: list[ft.Control],
    ) -> None:
        """Initialise the DismissibleDetailCard.

        Args:
            card_id: Unique identifier string for this card.
            title: The card title text displayed in the header.
            settings: Application GUI settings instance for font sizes.
            dismiss_callback: Callback invoked when the close button is
                clicked. Receives the card instance as argument.
            body_controls: List of Flet control objects to display in
                the card body below the header.
        """
        super().__init__()
        self._card_id = card_id

        def remove_card(e: ft.Event) -> None:
            """Callback to trigger the dismiss callback."""
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
