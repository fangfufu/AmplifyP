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
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
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
        """Initialise the DimerCard.

        Args:
            d: The PrimerDimer object containing sequence alignment data.
            settings: Application GUI settings instance.
            font_family: Font family for sequence display.
        """
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
                    diagram,
                ],
                spacing=8,
            ),
        )

    def _build_alignment_diagram(self, font_size_default: int) -> ft.Container:
        """Build the visual alignment container control for a dimer.

        Creates a three-line text visualisation showing the two primer
        sequences aligned at their binding interface with strength
        indicators in the middle line.

        Args:
            font_size_default: Font size for the sequence display.

        Returns:
            A Container with the alignment diagram and border styling.
        """
        seq1 = self.d.primer_1.seq
        seq2 = self.d.primer_2.seq

        middle_str = self.d.binding_strength_str

        # Build visually aligned lines.
        p2_label = f"{self.d.primer_2.name}: "
        p1_label = f"{self.d.primer_1.name}: "
        label_width = max(len(p2_label), len(p1_label))

        p2_prefix = p2_label.ljust(label_width)
        p1_prefix = p1_label.ljust(label_width)

        p2_line = f"{p2_prefix}5'-{seq2}-3'"
        mid_line = f"{' ' * label_width}{' ' * (3 + self.d.p1_pos)}{middle_str}"
        p1_line = f"{p1_prefix}{' ' * self.d.p1_pos}3'-{seq1[::-1]}-5'"

        # Create visual alignment stack using generic helper
        diagram_stack = create_overlapped_sequence_view(
            p2_line,
            mid_line,
            p1_line,
            font_family=self.font_family,
            font_size=font_size_default,
            is_dimer=True,
        )
        return ft.Container(
            content=ft.Row(
                [diagram_stack],
                scroll=ft.ScrollMode.ALWAYS,
            ),
            padding=12,
            border_radius=6,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            height=110,
        )

    def _build_card_header(
        self, font_size_subheader: int, font_size_small: int
    ) -> ft.Column:
        """Build the header row of the dimer card.

        Displays the dimer title (primer names or self-dimer label),
        overlap length, and quality score in styled containers.

        Args:
            font_size_subheader: Font size for the dimer title.
            font_size_small: Font size for the overlap and quality labels.

        Returns:
            A Column containing the title and metric containers.
        """
        p1_name = self.d.primer_1.name
        p2_name = self.d.primer_2.name
        seq1 = self.d.primer_1.seq
        seq2 = self.d.primer_2.seq

        is_self = p1_name == p2_name and seq1 == seq2
        dimer_title = (
            f"{p1_name} self-dimer" if is_self else f"{p1_name} vs {p2_name}"
        )
        quality_text = f"Quality: {self.d.quality:.1f}"

        return ft.Column(
            [
                ft.Row(
                    [
                        ft.Text(
                            dimer_title,
                            weight=ft.FontWeight.BOLD,
                            size=font_size_subheader,
                            selectable=True,
                        ),
                    ],
                    alignment=ft.MainAxisAlignment.START,
                ),
                ft.Row(
                    [
                        ft.Container(
                            content=ft.Text(
                                f"Overlap: {self.d.overlap} bp",
                                weight=ft.FontWeight.BOLD,
                                color=GUIColours.DIAGRAM_BLACK,
                                size=font_size_small,
                            ),
                            bgcolor=GUIColours.SELECTED_ROW_BG,
                            padding=ft.Padding(8, 4, 8, 4),
                            border_radius=4,
                        ),
                        ft.Container(
                            content=ft.Text(
                                quality_text,
                                weight=ft.FontWeight.BOLD,
                                color=GUIColours.DIAGRAM_BLACK,
                                size=font_size_small,
                            ),
                            bgcolor=GUIColours.SELECTED_ROW_BG,
                            padding=ft.Padding(8, 4, 8, 4),
                            border_radius=4,
                        ),
                    ],
                    spacing=8,
                    vertical_alignment=ft.CrossAxisAlignment.CENTER,
                    alignment=ft.MainAxisAlignment.START,
                ),
            ],
            spacing=4,
        )
