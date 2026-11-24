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

"""Diagram panel widget for rendering PCR execution targets."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

import flet as ft
import flet.canvas as cv

from amplifyp.amplicon import AmpliconRanking, sort_amplicons
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import MAX_AMPLICONS_RENDER, GUISettings
from amplifyp.pcr import PCR

if TYPE_CHECKING:
    from amplifyp.amplicon import Amplicon
    from amplifyp.repliconf import DirIdx, Repliconf

from .amplicon_drawing import DrawnAmplicon
from .pcr_layout import PCRLayoutSolver
from .primer_drawing import DrawnPrimer

__all__ = ["PCRDrawingPanel", "sort_amplicons"]

_AUTORESIZE_VIEWPORT_FRACTIONS: dict[str, float] = {
    "autoresize_1_4": 1.0 / 4.0,
    "autoresize": 1.0 / 3.0,
    "autoresize_1_2": 1.0 / 2.0,
    "autoresize_2_3": 2.0 / 3.0,
}
_NAV_BAR_HEIGHT = 100.0


class PCRDrawingPanel(ft.Column):  # type: ignore[misc]
    """Custom control encapsulating the PCR diagram canvas, stack, and resizer.

    This wraps the Canvas and GestureDetector elements.
    """

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        on_primer_click: Callable[[str, int, Repliconf, DirIdx], None],
        on_amplicon_click: Callable[[Amplicon], None],
    ) -> None:
        """Initialise the PCRDrawingPanel.

        Args:
            page: The Flet page instance for UI updates.
            settings: Application GUI settings instance.
            on_primer_click: Callback invoked when a primer is clicked.
                Receives (name, index, conf, var) arguments.
            on_amplicon_click: Callback invoked when an amplicon is clicked.
                Receives the amplicon object as argument.
        """
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
            border=ft.Border.all(1, GUIColours.OUTLINE),
            border_radius=5,
            padding=10,
            height=300,
        )
        self.divider = ft.GestureDetector(
            on_pan_update=self._on_pan_update,
            on_pan_end=self._on_pan_end,
            content=ft.Container(
                height=5,
                bgcolor=GUIColours.DIVIDER_GREY,
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

    def _get_viewport_height(self) -> float:
        """Return available viewport height excluding the top navigation bar."""
        page_h = getattr(self.app_page, "height", None)
        if isinstance(page_h, (int, float)) and page_h > 0:
            return max(0.0, float(page_h) - _NAV_BAR_HEIGHT)
        return 800.0 - _NAV_BAR_HEIGHT

    def _apply_container_height(self) -> None:
        """Size the diagram container according to the divider behaviour."""
        viewport_height = self._get_viewport_height()
        behaviour = str(
            self.settings.get("pcr_divider_behaviour", "autoresize")
        )
        if behaviour == "save_proportionality":
            ratio = self.settings._safe_float("pcr_map_height_ratio", 0.375)
            max_height = max(150.0, viewport_height - 150.0)
            self.diagram_container.height = max(
                150.0, min(max_height, viewport_height * ratio)
            )
            return

        # Default: autoresize to fit content, capped at viewport fraction.
        fraction = _AUTORESIZE_VIEWPORT_FRACTIONS.get(behaviour, 1.0 / 3.0)
        content_height = float(self.diagram_canvas.height or 300.0) + 22.0
        max_auto_height = max(150.0, viewport_height * fraction)
        self.diagram_container.height = max(
            150.0, min(content_height, max_auto_height)
        )

    def _on_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle vertical resizing of the diagram container.

        Args:
            e: The Flet drag update event containing delta y.
        """
        delta_y = getattr(e.local_delta, "y", 0.0) if e.local_delta else 0.0
        viewport_height = self._get_viewport_height()
        max_height = max(150.0, viewport_height - 150.0)
        new_height = max(
            150.0,
            min(
                max_height,
                float(self.diagram_container.height or 300.0) + delta_y,
            ),
        )
        self.diagram_container.height = new_height
        if (
            self.settings.get("pcr_divider_behaviour", "autoresize")
            == "save_proportionality"
            and viewport_height > 0
        ):
            self.settings["pcr_map_height_ratio"] = new_height / viewport_height
        try:
            if self.page:
                self.update()
        except RuntimeError:
            pass

    def _on_pan_end(self, _e: ft.DragEndEvent) -> None:
        """Persist map height proportionality when drag finishes."""
        if (
            self.app_page
            and self.settings.get("pcr_divider_behaviour", "autoresize")
            == "save_proportionality"
        ):
            self.settings.save_to_local(
                self.app_page, keys=("pcr_map_height_ratio",)
            )

    def reset_ui(self) -> None:
        """Reset the PCR view diagram canvas shapes and controls.

        Clears all shapes, hides the diagram container and divider,
        and resets the stack to contain only the canvas.
        """
        self.diagram_canvas.shapes = []
        self.diagram_stack.controls.clear()
        self.diagram_stack.controls.append(self.diagram_canvas)
        self.diagram_container.visible = False
        self.divider.visible = False
        self._sorted_amplicons: list[Amplicon] = []

    def render_diagram(self, pcr: PCR) -> None:
        """Perform coordinates calculations and render diagram elements.

        This draws baseline, primers, and amplicons.
        Sorts amplicons according to the configured ranking mode and limits
        rendering to MAX_AMPLICONS_RENDER if there are too many.

        Args:
            pcr: The PCR simulation instance containing amplicons and
                template data.
        """
        amplicons = pcr.amplicons
        num_amplicons = len(amplicons)

        if amplicons:
            ranking_mode = str(
                self.settings.get(
                    "pcr_amplicon_ranking",
                    AmpliconRanking.POSITION_THEN_LENGTH,
                )
            )
            try:
                amplicons = sort_amplicons(amplicons, ranking_mode)
            except ValueError:
                amplicons = sort_amplicons(
                    amplicons, AmpliconRanking.QUALITY_SCORE
                )
            num_amplicons = len(amplicons)

            # Limit to top amplicons if there are too many
            if num_amplicons > MAX_AMPLICONS_RENDER:
                amplicons = amplicons[:MAX_AMPLICONS_RENDER]
                num_amplicons = MAX_AMPLICONS_RENDER

        self._sorted_amplicons = amplicons
        target_length = len(pcr.template)

        fwd_bindings, rev_bindings = PCRLayoutSolver.collect_primer_bindings(
            pcr, amplicons
        )

        if not fwd_bindings and not rev_bindings:
            self.reset_ui()
            return

        self.diagram_container.border = ft.Border.all(1, GUIColours.OUTLINE)
        self.diagram_container.visible = True
        self.divider.visible = True

        label_orientation = str(
            self.settings.get("pcr_primer_label_orientation", "Horizontal")
        )
        font_size_primer = self.settings._safe_int("font_size_map_primer", 13)

        v_target, h_margin, c_width, t_width, v_frag_start = (
            PCRLayoutSolver.calculate_canvas_dimensions(
                target_length,
                num_amplicons,
                fwd_bindings,
                rev_bindings,
                self.app_page.width,
                label_orientation=label_orientation,
                font_size_primer=font_size_primer,
            )
        )

        # Apply dimensions directly to UI controls
        self.diagram_canvas.width = c_width
        self.diagram_canvas.height = v_frag_start + num_amplicons * 35 + 30.0
        self.diagram_stack.height = self.diagram_canvas.height
        self._apply_container_height()

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
            c_width=c_width,
            label_orientation=label_orientation,
            font_size_primer=font_size_primer,
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

        # Force a new list reference for canvas shapes to trigger
        # Flet change detection
        self.diagram_canvas.shapes = list(self.diagram_canvas.shapes)

    def _draw_template_baseline(
        self,
        v_target: float,
        h_margin: float,
        c_width: float,
        t_width: float,
        target_length: int,
    ) -> None:
        """Draw boundary lines, baseline, ticks, and texts.

        Draws vertical boundary lines at start/end, position labels (1
        and target_length), the main baseline, and dynamic tick marks
        scaled based on template length.

        Args:
            v_target: Vertical position of the baseline.
            h_margin: Horizontal margin in pixels.
            c_width: Total canvas width in pixels.
            t_width: Template drawing width in pixels.
            target_length: Total length of the template in base pairs.
        """
        # Draw vertical boundary lines at start and end of template
        boundary_paint = ft.Paint(
            color=GUIColours.DIAGRAM_BLACK,
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
                    color=GUIColours.DIAGRAM_BLACK,
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
                    color=GUIColours.DIAGRAM_BLACK,
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
                    color=GUIColours.DIAGRAM_BLACK,
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
            color=GUIColours.DIAGRAM_BLACK,
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
        fwd_bindings: dict[
            tuple[int, str], tuple[str, float, Repliconf, DirIdx]
        ],
        rev_bindings: dict[
            tuple[int, str], tuple[str, float, Repliconf, DirIdx]
        ],
        target_length: int,
        t_width: float,
        h_margin: float,
        v_target: float,
        c_width: float = 800.0,
        label_orientation: str = "Horizontal",
        font_size_primer: int = 13,
    ) -> None:
        """Draw forward and reverse primers using DrawnPrimer instances.

        Calculates shifted x-positions to prevent overlap, then creates
        and draws DrawnPrimer objects for each binding site.

        Args:
            fwd_bindings: Forward primer binding data dict.
            rev_bindings: Reverse primer binding data dict.
            target_length: Total length of the template in base pairs.
            t_width: Template drawing width in pixels.
            h_margin: Horizontal margin in pixels.
            v_target: Vertical position of the baseline.
            c_width: Total canvas width in pixels.
            label_orientation: Primer label orientation ('Horizontal', '45°',
                or 'Vertical').
            font_size_primer: Font size in pixels for primer labels.
        """
        fwd_shifted = PCRLayoutSolver.calculate_shifted_x(
            fwd_bindings, target_length, t_width, h_margin
        )
        rev_shifted = PCRLayoutSolver.calculate_shifted_x(
            rev_bindings, target_length, t_width, h_margin
        )

        fwd_label_layout: dict[tuple[int, str], tuple[float, int, float]] = {}
        rev_label_layout: dict[tuple[int, str], tuple[float, int, float]] = {}
        if label_orientation == "Horizontal":
            fwd_label_layout, _ = (
                PCRLayoutSolver.calculate_horizontal_label_layout(
                    fwd_bindings,
                    fwd_shifted,
                    c_width,
                    font_size_primer,
                    is_fwd=True,
                )
            )
            rev_label_layout, _ = (
                PCRLayoutSolver.calculate_horizontal_label_layout(
                    rev_bindings,
                    rev_shifted,
                    c_width,
                    font_size_primer,
                    is_fwd=False,
                )
            )

        def make_click(
            n_val: str, idx_val: int, c_val: Repliconf, v_val: DirIdx
        ) -> Callable[[], None]:
            """Return a lambda that calls the on_primer_click callback."""
            return lambda: self.on_primer_click(n_val, idx_val, c_val, v_val)

        for (start_idx, name), (
            name_val,
            S,
            fwd_conf,
            fwd_var,
        ) in fwd_bindings.items():
            label_left, label_row, label_width = fwd_label_layout.get(
                (start_idx, name), (None, 0, None)
            )
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
                label_left=label_left,
                label_row=label_row,
                label_width=label_width,
            )
            drawn.draw(self.diagram_canvas, self.diagram_stack)

        for (end_idx, name), (
            name_val,
            s,
            rev_conf,
            rev_var,
        ) in rev_bindings.items():
            label_left, label_row, label_width = rev_label_layout.get(
                (end_idx, name), (None, 0, None)
            )
            drawn = DrawnPrimer(
                name=name_val,
                index=end_idx,
                conf=rev_conf,
                var=rev_var,
                S=s,
                target_length=target_length,
                t_width=t_width,
                h_margin=h_margin,
                v_target=v_target,
                settings=self.settings,
                on_click=make_click(name_val, end_idx, rev_conf, rev_var),
                x_shifted=rev_shifted.get((end_idx, name)),
                label_left=label_left,
                label_row=label_row,
                label_width=label_width,
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
        amplicons: list[Amplicon] | None = None,
        v_frag_start: float | None = None,
    ) -> None:
        """Draw amplicons using DrawnAmplicon instances.

        Args:
            pcr: The PCR simulation instance.
            target_length: Total length of the template in base pairs.
            t_width: Template drawing width in pixels.
            h_margin: Horizontal margin in pixels.
            v_target: Vertical position of the baseline.
            c_width: Total canvas width in pixels.
            amplicons: List of amplicon objects to draw, or None to use
                pcr.amplicons.
            v_frag_start: Vertical start position for fragment rows, or
                None to use a default.
        """
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
