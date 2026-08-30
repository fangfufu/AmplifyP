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

"""PrimerItemCard component for generated primers list."""

from collections.abc import Callable

import flet as ft

from amplifyp.dimer import PrimerDimer
from amplifyp.dna import DNADirection
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings


class PrimerItemCard(ft.Card):  # type: ignore[misc]
    """Card item representing a generated primer in the list view."""

    def __init__(
        self,
        dimer: PrimerDimer,
        step_index: int,
        mode: DNADirection,
        settings: GUISettings,
        on_select_callback: Callable[[PrimerDimer, int], None],
    ) -> None:
        """Initialise PrimerItemCard."""
        font_size_default = settings.get("font_size_default", 14)
        font_size_header = settings.get("font_size_header", 18)

        is_reverse = mode == DNADirection.REV
        seq_align = ft.TextAlign.RIGHT if is_reverse else ft.TextAlign.LEFT
        seq = dimer.primer_1.seq
        length = len(seq)
        overlap_str = f"Overlap: {dimer.overlap} bp"

        card_container = ft.Container(
            content=ft.Row(
                [
                    ft.Column(
                        [
                            ft.Text(
                                f"{length} nt",
                                weight=ft.FontWeight.BOLD,
                                size=font_size_header,
                                text_align=ft.TextAlign.LEFT,
                            ),
                            ft.Container(
                                content=ft.TextField(
                                    value=seq,
                                    read_only=True,
                                    dense=True,
                                    expand=True,
                                    text_size=font_size_default,
                                    text_style=ft.TextStyle(
                                        font_family="Roboto Mono"
                                    ),
                                    text_align=seq_align,
                                    content_padding=ft.Padding(6, 4, 6, 4),
                                ),
                                expand=True,
                                alignment=(
                                    ft.Alignment(1, 0)
                                    if is_reverse
                                    else ft.Alignment(-1, 0)
                                ),
                            ),
                        ],
                        spacing=2,
                        expand=True,
                        horizontal_alignment=ft.CrossAxisAlignment.START,
                    ),
                    ft.Column(
                        [
                            ft.Container(
                                content=ft.Text(
                                    f"Quality: {round(dimer.quality)}",
                                    weight=ft.FontWeight.BOLD,
                                    size=font_size_default,
                                    color=GUIColours.DIAGRAM_BLACK,
                                ),
                                bgcolor=GUIColours.SELECTED_ROW_BG,
                                padding=ft.Padding(6, 3, 6, 3),
                                border_radius=4,
                            ),
                            ft.Container(
                                content=ft.Text(
                                    overlap_str,
                                    weight=ft.FontWeight.BOLD,
                                    size=font_size_default,
                                    color=GUIColours.DIAGRAM_BLACK,
                                ),
                                bgcolor=GUIColours.SELECTED_ROW_BG,
                                padding=ft.Padding(6, 3, 6, 3),
                                border_radius=4,
                            ),
                        ],
                        spacing=4,
                        horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                ],
                alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
            ),
            padding=10,
            ink=True,
            on_click=lambda _ev: on_select_callback(dimer, step_index),
        )

        super().__init__(content=card_container)
