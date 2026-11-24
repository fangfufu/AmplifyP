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

"""Primer drawing class and context card helpers for the PCRView."""

import math
from collections.abc import Callable
from typing import TYPE_CHECKING

import flet as ft
import flet.canvas as cv

from amplifyp.dna import DNA, DNADirection
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.repliconf import DirIdx, Repliconf

from .dismissible_detail_card import DismissibleDetailCard

if TYPE_CHECKING:
    from amplifyp.origin import ReplicationOrigin


class DrawnPrimer:
    """Class representing a primer drawn on the PCR diagram canvas."""

    def __init__(
        self,
        name: str,
        index: int,
        conf: Repliconf,
        var: DirIdx,
        S: float,
        target_length: int,
        t_width: float,
        h_margin: float,
        v_target: float,
        settings: GUISettings,
        on_click: Callable[[], None],
        x_shifted: float | None = None,
        label_left: float | None = None,
        label_row: int = 0,
        label_width: float | None = None,
    ) -> None:
        """Initialise the DrawnPrimer.

        Args:
            name: The primer name for display.
            index: The binding position index on the template.
            conf: The Repliconf configuration object for this primer.
            var: The variable containing direction information.
            S: Triangle size factor based on quality score.
            target_length: Total length of the template in base pairs.
            t_width: Template drawing width in pixels.
            h_margin: Horizontal margin in pixels.
            v_target: Vertical position of the baseline.
            settings: Application GUI settings instance.
            on_click: Callback invoked when the primer is clicked.
            x_shifted: Optional overridden x-coordinate for shifted
                positioning to prevent overlap.
            label_left: Optional overridden x-coordinate for horizontal
                label placement.
            label_row: Row index for stacking horizontal labels without
                overlap.
            label_width: Optional precomputed horizontal label width.
        """
        self.name = name
        self.index = index
        self.conf = conf
        self.var = var
        self.S = S
        self.target_length = target_length
        self.t_width = t_width
        self.h_margin = h_margin
        self.v_target = v_target
        self.settings = settings
        self.on_click = on_click
        self.x_shifted = x_shifted
        self.label_left = label_left
        self.label_row = label_row
        self.label_width = label_width

        # Calculate coordinates
        self.x_pos = (
            h_margin + (index / target_length * t_width)
            if target_length
            else h_margin
        )
        self.direction = var.direction

    def _add_label(self, stack: ft.Stack, x_render: float) -> None:
        """Add a horizontal, 45-degree, or vertical primer label to stack."""
        orientation = str(
            self.settings.get("pcr_primer_label_orientation", "Horizontal")
        )
        font_size = self.settings._safe_int("font_size_map_primer", 13)
        is_fwd = self.direction == DNADirection.FWD
        colour = GUIColours.FWD_PRIMER if is_fwd else GUIColours.REV_LABEL

        if orientation in ("Vertical", "45°"):
            top = (
                self.v_target - 25 - self.S - font_size
                if is_fwd
                else self.v_target + 25 + self.S - 3
            )
            c_width = self.t_width + 2.0 * self.h_margin
            width = (
                self.label_width
                if self.label_width is not None
                else max(20.0, len(self.name) * (font_size * 0.72) + 4.0)
            )
            if orientation == "45°" and not is_fwd:
                left = x_render - width
                angle = -math.pi / 4
                alignment = ft.Alignment(1, 0)
                text_control = ft.Text(
                    self.name,
                    color=colour,
                    size=font_size,
                    weight=ft.FontWeight.BOLD,
                    width=width,
                    text_align=ft.TextAlign.RIGHT,
                )
            elif orientation == "45°":
                left = x_render
                angle = -math.pi / 4
                alignment = ft.Alignment(-1, 0)
                text_control = ft.Text(
                    self.name,
                    color=colour,
                    size=font_size,
                    weight=ft.FontWeight.BOLD,
                )
            else:
                half_h = font_size / 2.0
                left = max(
                    4.0 + half_h,
                    min(c_width - 4.0 - half_h, x_render),
                )
                angle = -math.pi / 2 if is_fwd else math.pi / 2
                alignment = ft.Alignment(-1, 0)
                text_control = ft.Text(
                    self.name,
                    color=colour,
                    size=font_size,
                    weight=ft.FontWeight.BOLD,
                )
            stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=lambda _: self.on_click(),
                    content=text_control,
                    left=left,
                    top=top,
                    rotate=ft.Rotate(angle, alignment=alignment),
                )
            )
            return

        row_height = font_size + 5.0
        half_text_h = (font_size + 3.0) / 2.0
        width = (
            self.label_width
            if self.label_width is not None
            else max(20.0, len(self.name) * (font_size * 0.72) + 4.0)
        )
        if is_fwd:
            left = (
                self.label_left
                if self.label_left is not None
                else x_render + self.S / 2.0 + 4.0
            )
            top = (
                self.v_target
                - 25
                - self.S / 2.0
                - half_text_h
                - self.label_row * row_height
            )
        else:
            left = (
                self.label_left
                if self.label_left is not None
                else x_render - self.S / 2.0 - 4.0 - width
            )
            top = (
                self.v_target
                + 25
                + self.S / 2.0
                - half_text_h
                + self.label_row * row_height
            )
        if left < x_render:
            text_control = ft.Text(
                self.name,
                color=colour,
                size=font_size,
                weight=ft.FontWeight.BOLD,
                width=width,
                text_align=ft.TextAlign.RIGHT,
            )
        else:
            text_control = ft.Text(
                self.name,
                color=colour,
                size=font_size,
                weight=ft.FontWeight.BOLD,
            )
        stack.controls.append(
            ft.GestureDetector(
                mouse_cursor=ft.MouseCursor.CLICK,
                on_tap=lambda _: self.on_click(),
                content=text_control,
                left=left,
                top=top,
            )
        )

    def draw(self, canvas: cv.Canvas, stack: ft.Stack) -> None:
        """Draw the primer indicator elements onto the canvas and stack.

        Draws connector lines (straight or bent if shifted), direction-
        specific triangles (right-pointing blue for forward, left-pointing
        red for reverse), rotated labels, and click overlays.

        Args:
            canvas: The Flet canvas to draw shapes on.
            stack: The Flet stack to add labels and gesture detectors to.
        """
        x_render = self.x_shifted if self.x_shifted is not None else self.x_pos

        orientation = str(
            self.settings.get("pcr_primer_label_orientation", "Horizontal")
        )
        font_size = self.settings._safe_int("font_size_map_primer", 13)
        row_height = font_size + 5.0
        row_offset = (
            self.label_row * row_height if orientation == "Horizontal" else 0.0
        )

        if self.direction == DNADirection.FWD:
            tri_bottom = self.v_target - 25 - row_offset
            tri_mid = tri_bottom - self.S / 2
            tri_top = tri_bottom - self.S
            # Draw FWD primers (blue, float above baseline, pointing right)
            if abs(x_render - self.x_pos) > 0.1:
                # Bent leader line:
                # 1. vertical up
                # 2. diagonal to x_render
                # 3. vertical up to triangle
                canvas.shapes.append(
                    cv.Path(
                        [
                            cv.Path.MoveTo(self.x_pos, self.v_target),
                            cv.Path.LineTo(self.x_pos, self.v_target - 10),
                            cv.Path.LineTo(x_render, self.v_target - 20),
                            cv.Path.LineTo(x_render, tri_mid),
                        ],
                        paint=ft.Paint(
                            color=GUIColours.FWD_PRIMER,
                            style=ft.PaintingStyle.STROKE,
                            stroke_width=1.0,
                        ),
                    )
                )
            else:
                # Straight connector line
                canvas.shapes.append(
                    cv.Path(
                        [
                            cv.Path.MoveTo(self.x_pos, self.v_target),
                            cv.Path.LineTo(self.x_pos, tri_mid),
                        ],
                        paint=ft.Paint(
                            color=GUIColours.FWD_PRIMER,
                            style=ft.PaintingStyle.STROKE,
                            stroke_width=1.0,
                        ),
                    )
                )
            # Right-pointing triangle of size S:
            canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_render + self.S / 2, tri_mid),
                        cv.Path.LineTo(x_render - self.S / 2, tri_top),
                        cv.Path.LineTo(x_render - self.S / 2, tri_bottom),
                        cv.Path.Close(),
                    ],
                    paint=ft.Paint(
                        color=GUIColours.FWD_PRIMER,
                        style=ft.PaintingStyle.FILL,
                    ),
                )
            )
            # Add label in the stack
            self._add_label(stack, x_render)
            # Click overlay
            stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=lambda _: self.on_click(),
                    content=ft.Container(
                        bgcolor=GUIColours.TRANSPARENT,
                        width=20,
                        height=25 + self.S + row_offset,
                    ),
                    left=x_render - 10,
                    top=tri_top,
                )
            )
        else:
            tri_top = self.v_target + 25 + row_offset
            tri_mid = tri_top + self.S / 2
            tri_bottom = tri_top + self.S
            # Draw REV primers (red, float below baseline, pointing left)
            if abs(x_render - self.x_pos) > 0.1:
                # Bent leader line:
                # 1. vertical down
                # 2. diagonal to x_render
                # 3. vertical down to triangle
                canvas.shapes.append(
                    cv.Path(
                        [
                            cv.Path.MoveTo(self.x_pos, self.v_target),
                            cv.Path.LineTo(self.x_pos, self.v_target + 10),
                            cv.Path.LineTo(x_render, self.v_target + 20),
                            cv.Path.LineTo(x_render, tri_mid),
                        ],
                        paint=ft.Paint(
                            color=GUIColours.REV_PRIMER,
                            style=ft.PaintingStyle.STROKE,
                            stroke_width=1.0,
                        ),
                    )
                )
            else:
                # Straight connector line
                canvas.shapes.append(
                    cv.Path(
                        [
                            cv.Path.MoveTo(self.x_pos, self.v_target),
                            cv.Path.LineTo(self.x_pos, tri_mid),
                        ],
                        paint=ft.Paint(
                            color=GUIColours.REV_PRIMER,
                            style=ft.PaintingStyle.STROKE,
                            stroke_width=1.0,
                        ),
                    )
                )
            # Left-pointing triangle of size S:
            canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_render - self.S / 2, tri_mid),
                        cv.Path.LineTo(x_render + self.S / 2, tri_top),
                        cv.Path.LineTo(x_render + self.S / 2, tri_bottom),
                        cv.Path.Close(),
                    ],
                    paint=ft.Paint(
                        color=GUIColours.REV_PRIMER,
                        style=ft.PaintingStyle.FILL,
                    ),
                )
            )
            # Add label in the stack
            self._add_label(stack, x_render)
            # Click overlay
            stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=lambda _: self.on_click(),
                    content=ft.Container(
                        bgcolor=GUIColours.TRANSPARENT,
                        width=20,
                        height=25 + self.S + row_offset,
                    ),
                    left=x_render - 10,
                    top=self.v_target,
                )
            )


def get_template_substring(template: DNA, start: int, length: int) -> str:
    """Get substring of template DNA, supporting circular wrapping.

    Args:
        template: The DNA template object.
        start: The starting index for the substring.
        length: The number of bases to extract.

    Returns:
        A string of the requested DNA substring, with '-' for out-of-
        bounds positions on linear templates.
    """
    return template.padded_slice(start, length)


def format_context_lines(
    primer_name: str,
    padded_idx: int,
    conf: Repliconf,
    origin: "ReplicationOrigin",
    L: int,
    N: int,
    direction: DNADirection,
    improved_visualisation: bool = False,
    settings: GUISettings | None = None,
    separator: str | None = None,
) -> tuple[str, str, str]:
    """Format context alignment lines for binding site context map.

    Constructs multi-line text visualisation showing the primer binding
    site with upstream/downstream context, pipe markers, arrows, and
    optional complementary strand display.

    Args:
        primer_name: The name of the primer.
        padded_idx: The padded binding index on the template.
        conf: The Repliconf configuration object.
        origin: The replication origin object.
        L: Length of the primer.
        N: Length of the template.
        direction: The DNA direction (FWD or REV).
        improved_visualisation: If True, includes complementary strand
            display for forward primers.
        settings: Optional GUI settings to retrieve separator configuration.
        separator: Optional explicit separator override (' ' or '-').

    Returns:
        A tuple of (top_line, primer_line, bottom_line) strings for
        the context map visualisation.
    """
    _ = (L, N)
    if separator is not None:
        sep = "-" if (separator == "Dash ('-')" or separator == "-") else " "
    elif settings is not None and hasattr(settings, "get"):
        sep_setting = settings.get(
            "sequence_separator",
            settings.get("dimer_sequence_separator", "Space (' ')"),
        )
        sep = (
            "-" if (sep_setting == "Dash ('-')" or sep_setting == "-") else " "
        )
    else:
        sep = " "

    alignment = Repliconf.format_context_alignment(
        conf,
        DirIdx(direction, padded_idx),
        primer_name=primer_name,
        origin=origin,
        improved_visualisation=improved_visualisation,
        separator=sep,
    )
    return alignment.top_block, alignment.primer_line, alignment.bottom_block


class ReplicationContextCard(DismissibleDetailCard):
    """Card showing visually aligned replication context map."""

    def __init__(
        self,
        primer_name: str,
        padded_idx: int,
        conf: Repliconf,
        var: DirIdx,
        settings: GUISettings,
        dismiss_callback: Callable[[ft.Card], None],
    ) -> None:
        """Initialise the ReplicationContextCard.

        Displays primability, stability, quality scores, and a visually
        aligned replication context map for the specified primer binding
        site.

        Args:
            primer_name: The name of the primer.
            padded_idx: The padded binding index on the template.
            conf: The Repliconf configuration object.
            var: The variable containing direction information.
            settings: Application GUI settings instance.
            dismiss_callback: Callback invoked when the card is dismissed.
        """
        self._primer_name = primer_name
        self._padded_idx = padded_idx
        self._conf = conf
        self._var = var
        self._dismiss_callback = dismiss_callback

        card_id = f"context_{primer_name}_{padded_idx}"

        origin = conf.origin(var)

        primer_type = "" if var.direction == DNADirection.FWD else "Reverse"

        from amplifyp.gui.utils.data_helpers import (
            create_overlapped_sequence_view,
        )

        sep_setting = settings.get(
            "sequence_separator",
            settings.get("dimer_sequence_separator", "Space (' ')"),
        )
        sep = (
            "-" if (sep_setting == "Dash ('-')" or sep_setting == "-") else " "
        )
        alignment = Repliconf.format_context_alignment(
            conf,
            var,
            primer_name=primer_name,
            origin=origin,
            improved_visualisation=bool(
                settings.get("improved_visualisation", False)
            ),
            separator=sep,
        )

        font_family = settings.get("font_family", "Roboto Mono")
        font_size = settings.get("font_size_default", 14)
        diagram_text = create_overlapped_sequence_view(
            top_line=alignment.top_block,
            mid_line=alignment.primer_line,
            bottom_line=alignment.bonds_line,
            font_family=font_family,
            font_size=font_size,
            complement_line=alignment.complement_line,
            template_line=alignment.template_line,
        )

        if origin.settings.amplify4_compatibility_mode:
            prim_str = f"{int(origin.primability * 100)}%"
            stab_str = f"{int(origin.stability * 100)}%"
        else:
            prim_str = f"{origin.primability:.3f}"
            stab_str = f"{origin.stability:.3f}"

        font_size_small = settings.get("font_size_small", 12)
        title_controls = [
            ft.Row(
                [
                    ft.Container(
                        content=ft.Text(
                            f"Primeability: {prim_str}",
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
                            f"Stability: {stab_str}",
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
                            f"Quality: {origin.quality:.3f}",
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
                alignment=ft.MainAxisAlignment.END,
            ),
        ]
        body_controls = [
            ft.Container(
                content=ft.Row(
                    [diagram_text],
                    scroll=ft.ScrollMode.ALWAYS,
                ),
                padding=12,
                border_radius=6,
                border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            ),
        ]

        title = (
            f"Context Map - {primer_name} ({primer_type}) Binding Site"
            if primer_type
            else f"Context Map - {primer_name} Binding Site"
        )
        super().__init__(
            card_id=card_id,
            title=title,
            settings=settings,
            dismiss_callback=dismiss_callback,
            body_controls=body_controls,
            title_controls=title_controls,
        )

    def rebuild(self, settings: GUISettings) -> "ReplicationContextCard":
        """Rebuild this card with current theme colours and settings."""
        return ReplicationContextCard(
            primer_name=self._primer_name,
            padded_idx=self._padded_idx,
            conf=self._conf,
            var=self._var,
            settings=settings,
            dismiss_callback=self._dismiss_callback,
        )
