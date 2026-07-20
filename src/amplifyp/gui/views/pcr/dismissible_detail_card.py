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

"""DismissibleDetailCard base component for PCR view detail cards."""

from collections.abc import Callable, Sequence

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
        body_controls: Sequence[ft.Control],
        title_controls: Sequence[ft.Control] | None = None,
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
            title_controls: Optional list of Flet controls to add to the
                title row, placed between the title and close button.
        """
        super().__init__()
        self._card_id = card_id

        def remove_card(e: ft.Event) -> None:
            """Callback to trigger the dismiss callback."""
            dismiss_callback(self)

        title_row_controls: list[ft.Control] = [
            ft.Container(
                content=ft.Text(
                    title,
                    weight=ft.FontWeight.BOLD,
                    size=settings.get("font_size_subheader", 16),
                    selectable=True,
                ),
                expand=True,
            ),
        ]
        if title_controls:
            title_row_controls.extend(title_controls)
        title_row_controls.append(
            ft.IconButton(
                icon=ft.Icons.CLOSE,
                icon_size=18,
                tooltip="Dismiss",
                on_click=remove_card,
            )
        )

        self.content = ft.Container(
            padding=10,
            content=ft.Column(
                [
                    ft.Row(
                        title_row_controls,
                        alignment=ft.MainAxisAlignment.END,
                    ),
                    *body_controls,
                ]
            ),
        )
