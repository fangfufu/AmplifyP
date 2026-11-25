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

"""Dimer Card component."""

import flet as ft

from amplifyp.dimer import PrimerDimer
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.data_helpers import create_overlapped_sequence_view


class DimerCard(ft.Card):  # type: ignore[misc]
    """Card representing a single dimer with sequence alignment."""

    def __init__(
        self,
        d: PrimerDimer,
        settings: GUISettings,
        font_family: str = "Roboto Mono",
        show_names: bool = True,
    ) -> None:
        """Initialise the DimerCard.

        Args:
            d: The PrimerDimer object containing sequence alignment data.
            settings: Application GUI settings instance.
            font_family: Font family for sequence display.
            show_names: Whether to show primer names.
        """
        super().__init__()
        self.d = d
        self.settings = settings
        self.font_family = font_family
        self.show_names = show_names

        # Fetch font sizes from settings, fallback to defaults
        font_size_subheader = self.settings.get("font_size_subheader", 16)
        font_size_small = self.settings.get("font_size_small", 12)
        font_size_default = self.settings.get("font_size_default", 14)

        header = self._build_card_header(font_size_subheader, font_size_small)
        diagram = self._build_alignment_diagram(font_size_default)

        self.content = ft.Container(
            padding=ft.Padding(16, 16, 16, 16),
            content=ft.Column(
                [
                    header,
                    diagram,
                ],
                spacing=12,
            ),
        )

    def _build_alignment_diagram(self, font_size_default: int) -> ft.Container:
        """Build the visual alignment container control for a dimer.

        Creates a visual alignment showing the two primer sequences and
        their binding interface, optionally including their names.

        Args:
            font_size_default: Font size for the sequence display.

        Returns:
            A Container with the alignment diagram and border styling.
        """
        sep_setting = (
            self.settings.get("sequence_separator", "Dash ('-')")
            if hasattr(self.settings, "get")
            else " "
        )
        sep = "-" if str(sep_setting) in ("Dash ('-')", "-") else " "
        alignment = PrimerDimer.format_alignment(
            self.d,
            separator=sep,
            include_names=self.show_names,
        )

        # Create visual alignment stack using generic helper
        diagram_stack = create_overlapped_sequence_view(
            alignment.top_seq_line,
            alignment.bonds_line,
            alignment.bottom_seq_line,
            font_family=self.font_family,
            font_size=font_size_default,
            is_dimer=True,
            top_name_line=alignment.top_name_line,
            bottom_name_line=alignment.bottom_name_line,
        )
        return ft.Container(
            content=ft.Row(
                [diagram_stack],
                scroll=ft.ScrollMode.ALWAYS,
            ),
            padding=ft.Padding(14, 10, 14, 10),
            border_radius=6,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            height=150 if self.show_names else 120,
        )

    def _build_card_header(
        self, font_size_subheader: int, font_size_small: int
    ) -> ft.Row:
        """Build the header row of the dimer card.

        Displays the dimer title (primer names or self-dimer label),
        overlap length, and quality score in styled containers.

        Args:
            font_size_subheader: Font size for the dimer title.
            font_size_small: Font size for the overlap and quality labels.

        Returns:
            A Row containing the title and metric containers.
        """
        p1_name = self.d.primer_1.name
        p2_name = self.d.primer_2.name

        is_self = self.d.is_self_dimer
        if not self.show_names and is_self:
            dimer_title = "Self-dimer"
        else:
            dimer_title = (
                f"{p1_name} self-dimer"
                if is_self
                else f"{p1_name} vs {p2_name}"
            )
        quality_text = f"Quality: {round(self.d.quality)}"

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
                                quality_text,
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
                                f"Overlap: {self.d.overlap} bp",
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
                ),
            ],
            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            vertical_alignment=ft.CrossAxisAlignment.CENTER,
        )
