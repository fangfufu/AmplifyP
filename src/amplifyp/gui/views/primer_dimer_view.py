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
import flet.canvas as cv

from amplifyp.dimer import PrimerDimerGenerator
from amplifyp.dna import Primer
from amplifyp.gui.state import GUIState


class PrimerDimerView(ft.Column):  # type: ignore[misc]
    """View to show calculated primer dimer results."""

    def __init__(self, page: ft.Page, state: GUIState | None = None) -> None:
        """Initialize the PrimerDimerView."""
        super().__init__(expand=True)
        self.app_page = page
        self.state = state if state is not None else GUIState()

        self.result_list = ft.ListView(expand=True, spacing=10)
        self.controls = [
            ft.Container(content=self.result_list, expand=True),
        ]

    def run_analysis(self) -> None:
        """Run primer dimer analysis and populate the UI."""
        self.result_list.controls.clear()
        try:
            pd_settings = self.state.get_primer_dimer_settings()
            generator = PrimerDimerGenerator(settings=pd_settings)

            primers = self.state.get_active_primers()
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
                            size=16,
                            italic=True,
                            text_align=ft.TextAlign.CENTER,
                        ),
                        alignment=ft.Alignment(0, 0),
                        padding=20,
                    )
                )
            else:
                for d in dimers:
                    p1_name = d.primer_1.name
                    p2_name = d.primer_2.name
                    seq1 = d.primer_1.seq.upper()
                    seq2 = d.primer_2.seq.upper()

                    middle_str = d.binding_strength_str

                    # Build visually aligned lines.
                    p2_line = f"5'-{seq2}-3'"
                    mid_line = " " * (3 + d.p1_pos) + middle_str
                    p1_line = f"{' ' * d.p1_pos}3'-{seq1[::-1]}-5'"

                    # Ensure all three lines have the same total width so that
                    # the overlap region aligns vertically across all lines.
                    target_width = max(
                        len(p2_line), len(mid_line), len(p1_line)
                    )
                    p2_line = p2_line.ljust(target_width)
                    mid_line = mid_line.ljust(target_width)
                    p1_line = p1_line.ljust(target_width)

                    is_self = p1_name == p2_name and seq1 == seq2
                    dimer_title = (
                        f"{p1_name} self-dimer"
                        if is_self
                        else f"{p1_name} vs {p2_name}"
                    )

                    bg_color = (
                        ft.Colors.BLUE_ACCENT
                        if d.quality >= 100
                        else ft.Colors.ORANGE_ACCENT
                    )
                    quality_text = f"Quality: {d.quality:.1f}"

                    # Build Flet Canvas shapes for robust character alignment
                    shapes = []
                    char_width = 11.5  # Width of each character cell in pixels
                    x_offset = 12  # Left margin padding
                    y_p2 = 12  # Y coord for top primer line
                    y_mid = 38  # Y coord for middle bonds line
                    y_p1 = 64  # Y coord for bottom primer line

                    # Draw each line character-by-character
                    for idx, char in enumerate(p2_line):
                        if char != " ":
                            shapes.append(
                                cv.Text(
                                    x_offset + idx * char_width,
                                    y_p2,
                                    char,
                                    style=ft.TextStyle(
                                        font_family="Roboto Mono",
                                        size=14,
                                        weight=ft.FontWeight.BOLD,
                                        color=ft.Colors.ON_SURFACE,
                                    ),
                                )
                            )

                    for idx, char in enumerate(mid_line):
                        if char != " ":
                            shapes.append(
                                cv.Text(
                                    x_offset + idx * char_width,
                                    y_mid,
                                    char,
                                    style=ft.TextStyle(
                                        font_family="Roboto Mono",
                                        size=14,
                                        weight=ft.FontWeight.BOLD,
                                        color=ft.Colors.GREEN_400,
                                    ),
                                )
                            )

                    for idx, char in enumerate(p1_line):
                        if char != " ":
                            shapes.append(
                                cv.Text(
                                    x_offset + idx * char_width,
                                    y_p1,
                                    char,
                                    style=ft.TextStyle(
                                        font_family="Roboto Mono",
                                        size=14,
                                        weight=ft.FontWeight.BOLD,
                                        color=ft.Colors.ON_SURFACE,
                                    ),
                                )
                            )

                    # Dynamic width based on the number of characters
                    canvas_width = target_width * char_width + x_offset * 2
                    diagram_canvas = cv.Canvas(
                        shapes=shapes,
                        height=90,
                        width=canvas_width,
                    )
                    diagram_canvas.pointer_events = "none"

                    # Create a selectable text layer behind the canvas
                    # that perfectly matches the character spacing
                    selectable_text = ft.Text(
                        f"{p2_line}\n{mid_line}\n{p1_line}",
                        font_family="Roboto Mono",
                        size=14,
                        selectable=True,
                        color=ft.Colors.TRANSPARENT,
                    )
                    selectable_text.line_height = 1.85
                    selectable_text.left = x_offset - 2

                    diagram_stack = ft.Stack(
                        [
                            diagram_canvas,
                            selectable_text,
                        ],
                        width=canvas_width,
                        height=90,
                    )

                    self.result_list.controls.append(
                        ft.Card(
                            content=ft.Container(
                                padding=15,
                                content=ft.Column(
                                    [
                                        ft.Row(
                                            [
                                                ft.Text(
                                                    dimer_title,
                                                    weight=ft.FontWeight.BOLD,
                                                    size=16,
                                                    selectable=True,
                                                ),
                                                ft.Container(
                                                    content=ft.Text(
                                                        quality_text,
                                                        weight=ft.FontWeight.BOLD,
                                                        color=ft.Colors.WHITE,
                                                        size=12,
                                                    ),
                                                    bgcolor=bg_color,
                                                    padding=ft.Padding(
                                                        8, 4, 8, 4
                                                    ),
                                                    border_radius=4,
                                                ),
                                            ],
                                            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                                        ),
                                        ft.Text(
                                            f"Overlap length: {d.overlap} bp",
                                            size=13,
                                            color=ft.Colors.ON_SURFACE_VARIANT,
                                            selectable=True,
                                        ),
                                        ft.Container(height=8),
                                        ft.Container(
                                            content=ft.Row(
                                                [diagram_stack],
                                                scroll=ft.ScrollMode.AUTO,
                                            ),
                                            bgcolor=ft.Colors.SURFACE_CONTAINER_HIGHEST,
                                            padding=12,
                                            border_radius=6,
                                            border=ft.Border.all(
                                                1, ft.Colors.OUTLINE_VARIANT
                                            ),
                                            height=110,
                                        ),
                                    ]
                                ),
                            )
                        )
                    )
        except Exception as ex:
            self.result_list.controls.append(
                ft.Text(
                    f"Error running analysis: {ex}\n{traceback.format_exc()}",
                    color=ft.Colors.RED,
                )
            )
        self.app_page.update()
