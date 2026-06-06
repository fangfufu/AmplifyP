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

"""Primer Dimer View for the Flet application."""

import traceback

import flet as ft

from amplifyp.dimer import PrimerDimer, PrimerDimerGenerator
from amplifyp.dna import Primer
from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import create_overlapped_sequence_view, show_error_dialog


class PrimerDimerView(ft.Column):  # type: ignore[misc]
    """View to show calculated primer dimer results."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
    ) -> None:
        """Initialize the PrimerDimerView."""
        super().__init__(expand=True)
        self.app_page = page
        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()

        self.result_list = ft.ListView(
            expand=True, spacing=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.controls = [
            ft.Container(content=self.result_list, expand=True),
        ]

    def _build_alignment_diagram(
        self, d: PrimerDimer, font_family: str, font_size_default: int
    ) -> ft.Container:
        """Build the visual alignment container control for a primer dimer."""
        seq1 = d.primer_1.seq
        seq2 = d.primer_2.seq

        middle_str = d.binding_strength_str

        # Build visually aligned lines.
        p2_line = f"5'-{seq2}-3'"
        mid_line = " " * (3 + d.p1_pos) + middle_str
        p1_line = f"{' ' * d.p1_pos}3'-{seq1[::-1]}-5'"

        # Create visual alignment stack using generic helper
        diagram_stack = create_overlapped_sequence_view(
            p2_line,
            mid_line,
            p1_line,
            font_family=font_family,
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
        self, d: PrimerDimer, font_size_subheader: int, font_size_small: int
    ) -> ft.Row:
        """Build the header row of the primer dimer card."""
        p1_name = d.primer_1.name
        p2_name = d.primer_2.name
        seq1 = d.primer_1.seq
        seq2 = d.primer_2.seq

        is_self = p1_name == p2_name and seq1 == seq2
        dimer_title = (
            f"{p1_name} self-dimer" if is_self else f"{p1_name} vs {p2_name}"
        )
        quality_text = f"Quality: {d.quality:.1f}"

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
                                f"Overlap: {d.overlap} bp",
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

    def _create_primer_dimer_card(
        self, d: PrimerDimer, font_family: str = "Roboto Mono"
    ) -> ft.Card:
        """Create a Flet Card showing visually aligned primer dimer results."""
        # Fetch font sizes from settings, fallback to defaults
        font_size_subheader = self.settings.get("font_size_subheader", 16)
        font_size_small = self.settings.get("font_size_small", 12)
        font_size_default = self.settings.get("font_size_default", 14)

        header = self._build_card_header(
            d, font_size_subheader, font_size_small
        )
        diagram = self._build_alignment_diagram(
            d, font_family, font_size_default
        )

        return ft.Card(
            content=ft.Container(
                padding=15,
                content=ft.Column(
                    [
                        header,
                        ft.Container(height=8),
                        diagram,
                    ]
                ),
            )
        )

    def run_analysis(self) -> bool:
        """Run primer dimer analysis and populate the UI."""
        self.result_list.controls.clear()
        success = True
        try:
            pd_settings = self.settings.get_primer_dimer_settings()
            generator = PrimerDimerGenerator(settings=pd_settings)

            primers = self.input_data.get_active_primers()
            for p in primers:
                name = p["name"]
                seq = p["seq"]
                generator.add_primer(Primer(sequence=seq, name=name))

            generator.analyse_primers()

            dimers = generator.primer_dimers

            if not dimers:
                self.result_list.controls.append(
                    ft.Container(
                        content=ft.Text(
                            "No primer dimers detected above threshold.",
                            size=self.settings.get("font_size_subheader", 16),
                            italic=True,
                            text_align=ft.TextAlign.CENTER,
                        ),
                        alignment=ft.Alignment(0, 0),
                        padding=20,
                    )
                )
            else:
                font_family = self.settings.get("font_family", "Roboto Mono")
                for d in dimers:
                    card = self._create_primer_dimer_card(
                        d,
                        font_family=font_family,
                    )
                    self.result_list.controls.append(card)
        except Exception as ex:
            self.result_list.controls.append(
                ft.Text(
                    f"Error running analysis: {ex}\n{traceback.format_exc()}",
                    color=GUIColors.ERROR_RED,
                )
            )
            show_error_dialog(self.app_page, "Error running analysis", str(ex))
            success = False
        self.app_page.update()
        return success
