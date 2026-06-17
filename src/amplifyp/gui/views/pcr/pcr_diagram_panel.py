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

"""Diagram panel widget for rendering PCR execution targets."""

from collections.abc import Callable
from typing import Any

import flet as ft
import flet.canvas as cv

from amplifyp.gui.settings import MAX_AMPLICONS_RENDER, GUIColors, GUISettings
from amplifyp.pcr import PCR

from .amplicon_drawing import DrawnAmplicon
from .pcr_layout import PCRLayoutSolver
from .primer_drawing import DrawnPrimer

__all__ = ["PCRDrawingPanel"]


class PCRDrawingPanel(ft.Column):  # type: ignore[misc]
    """Custom control encapsulating the PCR diagram canvas, stack, and resizer.

    This wraps the Canvas and GestureDetector elements.
    """

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        on_primer_click: Callable[[str, int, Any, Any], None],
        on_amplicon_click: Callable[[Any], None],
    ) -> None:
        """Initialize the PCRDrawingPanel."""
        super().__init__(spacing=0, tight=True)
        self.app_page = page
        self.settings = settings
        self.on_primer_click = on_primer_click
        self.on_amplicon_click = on_amplicon_click

        self.diagram_canvas = cv.Canvas(
            expand=True,
            shapes=[],
        )
        self.diagram_stack = ft.Stack(
            controls=[self.diagram_canvas],
        )
        self.diagram_scrollable = ft.Column(
            controls=[self.diagram_stack],
            scroll=ft.ScrollMode.ALWAYS,
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
        )
        self.diagram_container = ft.Container(
            content=self.diagram_scrollable,
            visible=False,
            border=ft.Border.all(1, GUIColors.OUTLINE),
            border_radius=5,
            padding=10,
            height=300,
        )
        self.divider = ft.GestureDetector(
            on_pan_update=self._on_pan_update,
            content=ft.Container(
                height=5,
                bgcolor=GUIColors.DIVIDER_GREY,
                border_radius=5,
                margin=ft.Margin.symmetric(vertical=5),
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_UP_DOWN,
            visible=False,
        )

        self.controls = [
            self.diagram_container,
            self.divider,
        ]

    def _on_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle vertical resizing of the diagram container."""
        delta_y = getattr(e.local_delta, "y", 0.0) if e.local_delta else 0.0
        self.diagram_container.height = max(
            150.0, float(self.diagram_container.height or 300.0) + delta_y
        )
        self.app_page.update()

    def reset_ui(self) -> None:
        """Reset the PCR view diagram canvas shapes and controls."""
        self.diagram_canvas.shapes.clear()
        self.diagram_stack.controls.clear()
        self.diagram_stack.controls.append(self.diagram_canvas)
        self.diagram_container.visible = False
        self.divider.visible = False

    def render_diagram(self, pcr: PCR) -> None:
        """Perform coordinates calculations and render diagram elements.

        This draws baseline, primers, and amplicons.
        """
        amplicons = pcr.amplicons
        num_amplicons = len(amplicons)
        if num_amplicons == 0:
            self.reset_ui()
            return

        self.diagram_container.visible = True
        self.divider.visible = True

        # Sort amplicons by ascending q_score (lower is better quality)
        amplicons = sorted(amplicons, key=lambda a: a.q_score)
        num_amplicons = len(amplicons)

        # Limit to top amplicons if there are too many
        if num_amplicons > MAX_AMPLICONS_RENDER:
            amplicons = amplicons[:MAX_AMPLICONS_RENDER]
            num_amplicons = MAX_AMPLICONS_RENDER

        target_length = len(pcr.template)

        fwd_bindings, rev_bindings = PCRLayoutSolver.collect_primer_bindings(
            pcr, amplicons
        )

        v_target, h_margin, c_width, t_width, v_frag_start = (
            PCRLayoutSolver.calculate_canvas_dimensions(
                target_length,
                num_amplicons,
                fwd_bindings,
                rev_bindings,
                self.app_page.width,
            )
        )

        # Apply dimensions directly to UI controls
        self.diagram_canvas.width = c_width
        self.diagram_canvas.height = v_frag_start + num_amplicons * 35 + 30.0
        self.diagram_stack.height = self.diagram_canvas.height

        if target_length > 0:
            self._draw_template_baseline(
                v_target, h_margin, c_width, t_width, target_length
            )

        self._draw_primers(
            fwd_bindings,
            rev_bindings,
            target_length,
            t_width,
            h_margin,
            v_target,
        )

        self._draw_amplicons(
            pcr,
            target_length,
            t_width,
            h_margin,
            v_target,
            c_width,
            amplicons=amplicons,
            v_frag_start=v_frag_start,
        )

    def _draw_template_baseline(
        self,
        v_target: float,
        h_margin: float,
        c_width: float,
        t_width: float,
        target_length: int,
    ) -> None:
        """Draw boundary lines, baseline, ticks, and texts."""
        # Draw vertical boundary lines at start and end of template
        boundary_paint = ft.Paint(
            color=GUIColors.DIAGRAM_BLACK,
            style=ft.PaintingStyle.STROKE,
            stroke_width=1.0,
        )
        self.diagram_canvas.shapes.append(
            cv.Path(
                [
                    cv.Path.MoveTo(h_margin, v_target),
                    cv.Path.LineTo(h_margin, v_target - 65),
                ],
                paint=boundary_paint,
            )
        )
        self.diagram_canvas.shapes.append(
            cv.Path(
                [
                    cv.Path.MoveTo(c_width - h_margin, v_target),
                    cv.Path.LineTo(c_width - h_margin, v_target - 65),
                ],
                paint=boundary_paint,
            )
        )

        # Draw 1 and target_length text
        self.diagram_canvas.shapes.append(
            cv.Text(
                h_margin,
                v_target - 85,
                "1",
                style=ft.TextStyle(
                    size=self.settings.get("font_size_map_baseline", 16),
                    weight=ft.FontWeight.BOLD,
                    color=GUIColors.DIAGRAM_BLACK,
                ),
            )
        )
        self.diagram_canvas.shapes.append(
            cv.Text(
                c_width - h_margin,
                v_target - 85,
                str(target_length),
                style=ft.TextStyle(
                    size=self.settings.get("font_size_map_baseline", 16),
                    weight=ft.FontWeight.BOLD,
                    color=GUIColors.DIAGRAM_BLACK,
                ),
                alignment=ft.Alignment(1.0, -1.0),
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
                    color=GUIColors.DIAGRAM_BLACK,
                    style=ft.PaintingStyle.STROKE,
                    stroke_width=2.5,
                ),
            )
        )

        # Ticks dynamically scaled
        tick_interval = 100
        if target_length > 10000:
            tick_interval = 1000
        elif target_length > 5000:
            tick_interval = 500

        tick_paint = ft.Paint(
            color=GUIColors.DIAGRAM_BLACK,
            style=ft.PaintingStyle.STROKE,
            stroke_width=1.0,
        )
        for bp in range(0, target_length + 1, tick_interval):
            x_pos = h_margin + (bp / target_length * t_width)
            self.diagram_canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_pos, v_target - 3),
                        cv.Path.LineTo(x_pos, v_target + 3),
                    ],
                    paint=tick_paint,
                )
            )

    def _draw_primers(
        self,
        fwd_bindings: dict[tuple[int, str], tuple[str, float, Any, Any]],
        rev_bindings: dict[tuple[int, str], tuple[str, float, Any, Any]],
        target_length: int,
        t_width: float,
        h_margin: float,
        v_target: float,
    ) -> None:
        """Draw forward and reverse primers using DrawnPrimer instances."""
        fwd_shifted = PCRLayoutSolver.calculate_shifted_x(
            fwd_bindings, target_length, t_width, h_margin
        )
        rev_shifted = PCRLayoutSolver.calculate_shifted_x(
            rev_bindings, target_length, t_width, h_margin
        )

        def make_click(
            n_val: str, idx_val: int, c_val: Any, v_val: Any
        ) -> Callable[[], None]:
            return lambda: self.on_primer_click(n_val, idx_val, c_val, v_val)

        for (start_idx, name), (
            name_val,
            S,
            fwd_conf,
            fwd_var,
        ) in fwd_bindings.items():
            drawn = DrawnPrimer(
                name=name_val,
                index=start_idx,
                conf=fwd_conf,
                var=fwd_var,
                S=S,
                target_length=target_length,
                t_width=t_width,
                h_margin=h_margin,
                v_target=v_target,
                settings=self.settings,
                on_click=make_click(name_val, start_idx, fwd_conf, fwd_var),
                x_shifted=fwd_shifted.get((start_idx, name)),
            )
            drawn.draw(self.diagram_canvas, self.diagram_stack)

        for (end_idx, name), (
            name_val,
            S,
            rev_conf,
            rev_var,
        ) in rev_bindings.items():
            drawn = DrawnPrimer(
                name=name_val,
                index=end_idx,
                conf=rev_conf,
                var=rev_var,
                S=S,
                target_length=target_length,
                t_width=t_width,
                h_margin=h_margin,
                v_target=v_target,
                settings=self.settings,
                on_click=make_click(name_val, end_idx, rev_conf, rev_var),
                x_shifted=rev_shifted.get((end_idx, name)),
            )
            drawn.draw(self.diagram_canvas, self.diagram_stack)

    def _draw_amplicons(
        self,
        pcr: PCR,
        target_length: int,
        t_width: float,
        h_margin: float,
        v_target: float,
        c_width: float,
        amplicons: list[Any] | None = None,
        v_frag_start: float | None = None,
    ) -> None:
        """Draw amplicons using DrawnAmplicon instances."""
        if amplicons is None:
            amplicons = pcr.amplicons
        for idx, amp in enumerate(amplicons):
            drawn = DrawnAmplicon(
                amp=amp,
                idx=idx,
                target_length=target_length,
                t_width=t_width,
                h_margin=h_margin,
                v_target=v_target,
                c_width=c_width,
                settings=self.settings,
                on_click=lambda a=amp: self.on_amplicon_click(a),  # type: ignore[misc]
                v_frag_start=v_frag_start,
            )
            drawn.draw(self.diagram_canvas, self.diagram_stack)
