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

"""Result View for the Flet application."""

import traceback
from typing import TYPE_CHECKING, Any

import flet as ft
import flet.canvas as cv

from amplifyp.dna import DNA, DNAType, Primer
from amplifyp.pcr import PCR

if TYPE_CHECKING:
    pass


class ResultView(ft.Column):  # type: ignore[misc]
    """Result view for rendering PCR custom execution targets."""

    def __init__(
        self,
        page: ft.Page,
        state_or_input: Any = None,
        settings_view: Any = None,
    ) -> None:
        """Initialize the ResultView."""
        super().__init__(expand=True)
        self.app_page = page

        # Flexible initialization for decoupling and compatibility
        from amplifyp.gui.state import GUIState

        from .input_view import InputView

        if isinstance(state_or_input, GUIState):
            self.state = state_or_input
            self.input_view = None
            self.settings_view = settings_view
        elif isinstance(state_or_input, InputView):
            self.input_view = state_or_input
            self.settings_view = settings_view
            self.state = state_or_input.state
        else:
            self.state = GUIState()
            self.input_view = None
            self.settings_view = None

        self.result_list = ft.ListView(expand=True, spacing=10)
        self.diagram_canvas = cv.Canvas(
            expand=True,
            content=ft.Container(height=300),
            shapes=[],
        )
        self.diagram_container = ft.Container(
            content=self.diagram_canvas,
            visible=False,
            border=ft.Border.all(1, ft.Colors.OUTLINE),
            border_radius=5,
            padding=10,
            height=300,
        )

        self.controls = [
            self.diagram_container,
            ft.Container(content=self.result_list, expand=True),
        ]

    def run_pcr(self) -> None:
        """Execute the PCR simulation and update the UI."""
        self.result_list.controls.clear()
        self.diagram_canvas.shapes.clear()
        self.diagram_container.visible = False
        try:
            from amplifyp.gui.util import clean_sequence

            # Read parameters directly from the central state
            clean_template = clean_sequence(self.state.template)
            t_type = (
                DNAType.CIRCULAR
                if self.state.template_circular
                else DNAType.LINEAR
            )
            template_dna = DNA(clean_template, dna_type=t_type)

            rep_settings = self.state.get_replication_settings()
            pcr = PCR(template_dna, settings=rep_settings)

            primers = self.state.get_active_primers()
            for p in primers:
                name = p["name"]
                seq = clean_sequence(p["seq"])
                pcr.add_primer(Primer(sequence=seq, name=name))

            num_amplicons = pcr.predict_amplicons()

            if num_amplicons == 0:
                self.result_list.controls.append(ft.Text("No amplicons found."))
            else:
                self.diagram_container.visible = True

                # Canvas drawing constants based on Amplify4's `makePlot`
                # Drawing configuration
                v_target = 50.0  # Y position of target baseline
                h_margin = 20.0  # X padding
                c_width = 800.0  # Canvas effective drawing width
                t_width = c_width - (2.0 * h_margin)
                target_length = len(template_dna)

                if target_length > 0:
                    self.diagram_canvas.shapes.append(
                        cv.Text(h_margin, 10, "1", style=ft.TextStyle(size=12))
                    )
                    self.diagram_canvas.shapes.append(
                        cv.Text(
                            c_width - h_margin,
                            10,
                            str(target_length),
                            style=ft.TextStyle(size=12),
                            text_align=ft.TextAlign.RIGHT,
                        )
                    )

                    # Target Baseline
                    self.diagram_canvas.shapes.append(
                        cv.Path(
                            [
                                cv.Path.MoveTo(h_margin, v_target),
                                cv.Path.LineTo(c_width - h_margin, v_target),
                            ],
                            paint=ft.Paint(
                                style=ft.PaintingStyle.STROKE,
                                stroke_width=2.5,
                            ),
                        )
                    )

                    # Ticks every 100bp
                    tick_paint = ft.Paint(
                        style=ft.PaintingStyle.STROKE, stroke_width=1.0
                    )
                    for tick_pos in [0.0, 0.5, 1.0]:
                        x_pos = h_margin + (tick_pos * t_width)
                        self.diagram_canvas.shapes.append(
                            cv.Path(
                                [
                                    cv.Path.MoveTo(x_pos, v_target - 5),
                                    cv.Path.LineTo(x_pos, v_target + 5),
                                ],
                                paint=tick_paint,
                            )
                        )

                # Draw amplicons
                v_frag_start = v_target + 40
                v_frag_step = 40

                for idx, amp in enumerate(pcr.amplicons):
                    start_idx = amp.start.index
                    end_idx = amp.end.index

                    # X coords
                    x_start = (
                        h_margin + (start_idx / target_length * t_width)
                        if target_length
                        else h_margin
                    )
                    x_end = (
                        h_margin + (end_idx / target_length * t_width)
                        if target_length
                        else (c_width - h_margin)
                    )

                    y_pos = v_frag_start + (idx * v_frag_step)

                    # Amplicon Bar
                    amp_paint = ft.Paint(
                        color=ft.Colors.BLUE, style=ft.PaintingStyle.FILL
                    )
                    bar_height = 8.0
                    if amp.circular:
                        self.diagram_canvas.shapes.append(
                            cv.Path(
                                [
                                    cv.Path.MoveTo(h_margin, y_pos),
                                    cv.Path.LineTo(c_width - h_margin, y_pos),
                                    cv.Path.LineTo(
                                        c_width - h_margin,
                                        y_pos + bar_height,
                                    ),
                                    cv.Path.LineTo(
                                        h_margin, y_pos + bar_height
                                    ),
                                    cv.Path.Close(),
                                ],
                                paint=amp_paint,
                            )
                        )
                        label_x = h_margin + (t_width / 2.0)
                    else:
                        self.diagram_canvas.shapes.append(
                            cv.Path(
                                [
                                    cv.Path.MoveTo(x_start, y_pos),
                                    cv.Path.LineTo(x_end, y_pos),
                                    cv.Path.LineTo(x_end, y_pos + bar_height),
                                    cv.Path.LineTo(x_start, y_pos + bar_height),
                                    cv.Path.Close(),
                                ],
                                paint=amp_paint,
                            )
                        )
                        label_x = x_start + ((x_end - x_start) / 2.0)

                    # Amplicon Length Text
                    self.diagram_canvas.shapes.append(
                        cv.Text(
                            label_x,
                            y_pos + bar_height + 5,
                            f"{len(amp.product)}bp",
                            style=ft.TextStyle(size=11, color=ft.Colors.BLUE),
                            text_align=ft.TextAlign.CENTER,
                        )
                    )

                    # Card Details
                    self.result_list.controls.append(
                        ft.Card(
                            content=ft.Container(
                                padding=10,
                                content=ft.Column(
                                    [
                                        ft.Text(
                                            f"Length: {len(amp.product)} bp",
                                            weight=ft.FontWeight.BOLD,
                                            size=16,
                                        ),
                                        ft.Text(
                                            "Forward Primer: "
                                            f"{amp.fwd_origin.name} "
                                            f"(Start: {amp.start.index})"
                                        ),
                                        ft.Text(
                                            "Reverse Primer: "
                                            f"{amp.rev_origin.name} "
                                            f"(End: {amp.end.index})"
                                        ),
                                        ft.Text(
                                            "Quality Score: "
                                            f"{amp.q_score:.2f} - "
                                            f"{amp.q_score_report_str()}"
                                        ),
                                        ft.Text(
                                            "Amplified Sequence:",
                                            weight=ft.FontWeight.BOLD,
                                        ),
                                        ft.TextField(
                                            value=str(amp.product.seq),
                                            read_only=True,
                                            multiline=True,
                                            min_lines=3,
                                            max_lines=8,
                                            expand=True,
                                            text_style=ft.TextStyle(
                                                font_family="monospace"
                                            ),
                                        ),
                                    ]
                                ),
                            )
                        )
                    )
        except Exception as ex:
            self.result_list.controls.append(
                ft.Text(
                    f"Error: {ex}\n{traceback.format_exc()}",
                    color=ft.Colors.RED,
                )
            )

        self.app_page.update()
