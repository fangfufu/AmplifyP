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

"""Dimer Card component."""

import flet as ft

from amplifyp.dimer import PrimerDimer
from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.util import (
    create_overlapped_sequence_view,
)


class DimerCard(ft.Card):  # type: ignore[misc]
    """Card representing a single dimer with sequence alignment."""

    def __init__(
        self,
        d: PrimerDimer,
        settings: GUISettings,
        font_family: str = "Roboto Mono",
    ) -> None:
        """Initialize the DimerCard."""
        super().__init__()
        self.d = d
        self.settings = settings
        self.font_family = font_family

        # Fetch font sizes from settings, fallback to defaults
        font_size_subheader = self.settings.get("font_size_subheader", 16)
        font_size_small = self.settings.get("font_size_small", 12)
        font_size_default = self.settings.get("font_size_default", 14)

        header = self._build_card_header(font_size_subheader, font_size_small)
        diagram = self._build_alignment_diagram(font_size_default)

        self.content = ft.Container(
            padding=15,
            content=ft.Column(
                [
                    header,
                    ft.Container(height=8),
                    diagram,
                ]
            ),
        )

    def _build_alignment_diagram(self, font_size_default: int) -> ft.Container:
        """Build the visual alignment container control for a dimer."""
        seq1 = self.d.primer_1.seq
        seq2 = self.d.primer_2.seq

        middle_str = self.d.binding_strength_str

        # Build visually aligned lines.
        p2_line = f"5'-{seq2}-3'"
        mid_line = " " * (3 + self.d.p1_pos) + middle_str
        p1_line = f"{' ' * self.d.p1_pos}3'-{seq1[::-1]}-5'"

        # Create visual alignment stack using generic helper
        diagram_stack = create_overlapped_sequence_view(
            p2_line,
            mid_line,
            p1_line,
            font_family=self.font_family,
            font_size=font_size_default,
        )
        return ft.Container(
            content=ft.Row(
                [diagram_stack],
                scroll=ft.ScrollMode.ALWAYS,
            ),
            padding=12,
            border_radius=6,
            border=ft.Border.all(1, GUIColors.OUTLINE_VARIANT),
            height=110,
        )

    def _build_card_header(
        self, font_size_subheader: int, font_size_small: int
    ) -> ft.Row:
        """Build the header row of the dimer card."""
        p1_name = self.d.primer_1.name
        p2_name = self.d.primer_2.name
        seq1 = self.d.primer_1.seq
        seq2 = self.d.primer_2.seq

        is_self = p1_name == p2_name and seq1 == seq2
        dimer_title = (
            f"{p1_name} self-dimer" if is_self else f"{p1_name} vs {p2_name}"
        )
        quality_text = f"Quality: {self.d.quality:.1f}"

        return ft.Row(
            [
                ft.Text(
                    dimer_title,
                    weight=ft.FontWeight.BOLD,
                    size=font_size_subheader,
                    selectable=True,
                ),
                ft.Row(
                    [
                        ft.Container(
                            content=ft.Text(
                                f"Overlap: {self.d.overlap} bp",
                                weight=ft.FontWeight.BOLD,
                                color=GUIColors.TEXT_ON_SURFACE,
                                size=font_size_small,
                            ),
                            bgcolor=GUIColors.CONTAINER_HIGHEST,
                            padding=ft.Padding(8, 4, 8, 4),
                            border_radius=4,
                        ),
                        ft.Container(
                            content=ft.Text(
                                quality_text,
                                weight=ft.FontWeight.BOLD,
                                color=GUIColors.TEXT_ON_SURFACE,
                                size=font_size_small,
                            ),
                            bgcolor=GUIColors.CONTAINER_HIGHEST,
                            padding=ft.Padding(8, 4, 8, 4),
                            border_radius=4,
                        ),
                    ],
                    spacing=8,
                    vertical_alignment=ft.CrossAxisAlignment.CENTER,
                ),
            ],
            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            vertical_alignment=ft.CrossAxisAlignment.CENTER,
        )
