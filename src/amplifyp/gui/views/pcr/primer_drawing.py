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

"""Primer drawing class and context card helpers for the PCRView."""

from collections.abc import Callable
from typing import Any

import flet as ft
import flet.canvas as cv

from amplifyp.dna import DNA, DNADirection, DNAType
from amplifyp.gui.colors import GUIColors

from .dismissible_detail_card import DismissibleDetailCard


class DrawnPrimer:
    """Class representing a primer drawn on the PCR diagram canvas."""

    def __init__(
        self,
        name: str,
        index: int,
        conf: Any,
        var: Any,
        S: float,
        target_length: int,
        t_width: float,
        h_margin: float,
        v_target: float,
        settings: Any,
        on_click: Callable[[], None],
        x_shifted: float | None = None,
    ) -> None:
        """Initialize the DrawnPrimer."""
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

        # Calculate coordinates
        self.x_pos = (
            h_margin + (index / target_length * t_width)
            if target_length
            else h_margin
        )
        self.direction = var.direction

    def draw(self, canvas: cv.Canvas, stack: ft.Stack) -> None:
        """Draw the primer indicator elements onto the canvas and stack."""
        x_render = self.x_shifted if self.x_shifted is not None else self.x_pos

        if self.direction == DNADirection.FWD:
            # Draw FWD primers (blue, float above baseline, pointing down)
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
                            cv.Path.LineTo(x_render, self.v_target - 25),
                        ],
                        paint=ft.Paint(
                            color=GUIColors.FWD_PRIMER,
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
                            cv.Path.LineTo(self.x_pos, self.v_target - 25),
                        ],
                        paint=ft.Paint(
                            color=GUIColors.FWD_PRIMER,
                            style=ft.PaintingStyle.STROKE,
                            stroke_width=1.0,
                        ),
                    )
                )
            # Down-pointing triangle of size S:
            # Tip: (x_render, self.v_target - 25)
            # Top-Left: (x_render - S/2, self.v_target - 25 - S)
            # Top-Right: (x_render + S/2, self.v_target - 25 - S)
            canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_render, self.v_target - 25),
                        cv.Path.LineTo(
                            x_render - self.S / 2, self.v_target - 25 - self.S
                        ),
                        cv.Path.LineTo(
                            x_render + self.S / 2, self.v_target - 25 - self.S
                        ),
                        cv.Path.Close(),
                    ],
                    paint=ft.Paint(
                        color=GUIColors.FWD_PRIMER,
                        style=ft.PaintingStyle.FILL,
                    ),
                )
            )
            # Add rotated label in the stack
            stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=lambda _: self.on_click(),
                    content=ft.Text(
                        self.name,
                        color=GUIColors.FWD_PRIMER,
                        size=self.settings.get("font_size_map_primer", 13),
                        weight=ft.FontWeight.BOLD,
                    ),
                    left=x_render,
                    top=self.v_target - 25 - self.S - 13,
                    rotate=ft.Rotate(-1.5708, alignment=ft.Alignment(-1, 0)),
                )
            )
            # Click overlay
            stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=lambda _: self.on_click(),
                    content=ft.Container(
                        bgcolor=GUIColors.TRANSPARENT,
                        width=20,
                        height=25 + self.S,
                    ),
                    left=x_render - 10,
                    top=self.v_target - 25 - self.S,
                )
            )
        else:
            # Draw REV primers (red, float below baseline, pointing up)
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
                            cv.Path.LineTo(x_render, self.v_target + 25),
                        ],
                        paint=ft.Paint(
                            color=GUIColors.REV_PRIMER,
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
                            cv.Path.LineTo(self.x_pos, self.v_target + 25),
                        ],
                        paint=ft.Paint(
                            color=GUIColors.REV_PRIMER,
                            style=ft.PaintingStyle.STROKE,
                            stroke_width=1.0,
                        ),
                    )
                )
            # Up-pointing triangle of size S:
            # Tip: (x_render, self.v_target + 25)
            # Bottom-Left: (x_render - S/2, self.v_target + 25 + S)
            # Bottom-Right: (x_render + S/2, self.v_target + 25 + S)
            canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_render, self.v_target + 25),
                        cv.Path.LineTo(
                            x_render - self.S / 2, self.v_target + 25 + self.S
                        ),
                        cv.Path.LineTo(
                            x_render + self.S / 2, self.v_target + 25 + self.S
                        ),
                        cv.Path.Close(),
                    ],
                    paint=ft.Paint(
                        color=GUIColors.REV_PRIMER,
                        style=ft.PaintingStyle.FILL,
                    ),
                )
            )
            # Add rotated label in the stack
            stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=lambda _: self.on_click(),
                    content=ft.Text(
                        self.name,
                        color=GUIColors.REV_LABEL,
                        size=self.settings.get("font_size_map_primer", 13),
                        weight=ft.FontWeight.BOLD,
                    ),
                    left=x_render,
                    top=self.v_target + 25 + self.S - 3,
                    rotate=ft.Rotate(1.5708, alignment=ft.Alignment(-1, 0)),
                )
            )
            # Click overlay
            stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=lambda _: self.on_click(),
                    content=ft.Container(
                        bgcolor=GUIColors.TRANSPARENT,
                        width=20,
                        height=25 + self.S,
                    ),
                    left=x_render - 10,
                    top=self.v_target,
                )
            )


def get_template_substring(template: DNA, start: int, length: int) -> str:
    """Get substring of template DNA, supporting circular wrapping."""
    N_len = len(template)
    if template.type == DNAType.CIRCULAR:
        return "".join(
            template.seq[i % N_len] for i in range(start, start + length)
        )
    else:
        return "".join(
            template.seq[i] if 0 <= i < N_len else "-"
            for i in range(start, start + length)
        )


def format_context_lines(
    primer_name: str,
    padded_idx: int,
    conf: Any,
    origin: Any,
    L: int,
    N: int,
    direction: DNADirection,
    improved_visualisation: bool = False,
) -> tuple[str, str, str]:
    """Format context alignment lines for binding site context map."""
    if direction == DNADirection.FWD:
        start_genomic = (padded_idx - L) % N
        primer_display_seq = conf.primer.seq
        primer_label = f"{primer_name} (Forward)"
        primer_ends = ("5'", "3'")
        strength_display = origin.binding_strength_str[::-1]
    else:
        start_genomic = padded_idx % N
        primer_display_seq = conf.primer.seq[::-1]
        primer_label = f"{primer_name} (Reverse)"
        primer_ends = ("3'", "5'")
        strength_display = origin.binding_strength_str

    end_genomic = start_genomic + L
    start_num_str = str(((start_genomic) % N) + 1)
    end_num_str = str(((end_genomic - 1) % N) + 1)

    # Construct primer line:
    primer_line = (
        f"{primer_label:<29}"
        f"{primer_ends[0]}-{primer_display_seq}-{primer_ends[1]}"
    )

    # Construct strength line:
    bonds_line = f"{' ' * 12}{' ' * 20}{strength_display}"

    # Construct pipes and arrows:
    pipe_line = f"{' ' * 12}{' ' * 20}|{' ' * (L - 2)}|"
    arrow_line = f"{' ' * 12}{' ' * 20}V{' ' * (L - 2)}V"

    # Construct numbers:
    num_line_chars = [" "] * (12 + 20 + L + 20)
    col1 = 12 + 20 - len(start_num_str) // 2
    for idx, char in enumerate(start_num_str):
        num_line_chars[col1 + idx] = char
    col2 = (12 + 20 + L - 1) - len(end_num_str) // 2
    for idx, char in enumerate(end_num_str):
        num_line_chars[col2 + idx] = char
    top_line = "".join(num_line_chars).rstrip()

    # Combined top line:
    top_line = f"{top_line}\n{pipe_line}\n{arrow_line}"

    upstream_seq = get_template_substring(conf.template, start_genomic - 20, 20)
    binding_seq = get_template_substring(conf.template, start_genomic, L)
    downstream_seq = get_template_substring(
        conf.template, start_genomic + L, 20
    )

    context_line_prefix = "Context  "
    if improved_visualisation and direction == DNADirection.FWD:
        from amplifyp.dna import GLOBAL_COMPLEMENT_TABLE

        comp_upstream = upstream_seq.translate(GLOBAL_COMPLEMENT_TABLE)
        comp_binding = binding_seq.translate(GLOBAL_COMPLEMENT_TABLE)
        comp_downstream = downstream_seq.translate(GLOBAL_COMPLEMENT_TABLE)
        comp_line = (
            f"{' ' * 9}3'-{comp_upstream}{comp_binding}{comp_downstream}-5'"
        )
        bottom_line = (
            f"{bonds_line}\n"
            f"{comp_line}\n"
            f"{context_line_prefix}5'-{upstream_seq}{binding_seq}{downstream_seq}-3'"
        )
    else:
        bottom_line = (
            f"{bonds_line}\n"
            f"{context_line_prefix}5'-{upstream_seq}{binding_seq}{downstream_seq}-3'"
        )

    return top_line, primer_line, bottom_line


class ReplicationContextCard(DismissibleDetailCard):
    """Card showing visually aligned replication context map."""

    def __init__(
        self,
        primer_name: str,
        padded_idx: int,
        conf: Any,
        var: Any,
        settings: Any,
        dismiss_callback: Callable[[ft.Card], None],
    ) -> None:
        """Initialize the ReplicationContextCard."""
        card_id = f"context_{primer_name}_{padded_idx}"

        origin = conf.origin(var)
        if origin is None:
            body_controls = [ft.Text("Error: Replication origin not found")]
            super().__init__(
                card_id=card_id,
                title="Error",
                settings=settings,
                dismiss_callback=dismiss_callback,
                body_controls=body_controls,
            )
            return

        primer_type = (
            "Forward Primer"
            if var.direction == DNADirection.FWD
            else "Reverse Primer"
        )
        L = len(conf.primer)
        N = len(conf.template)

        from amplifyp.gui.util import create_overlapped_sequence_view

        top_line, primer_line, bottom_line = format_context_lines(
            primer_name=primer_name,
            padded_idx=padded_idx,
            conf=conf,
            origin=origin,
            L=L,
            N=N,
            direction=var.direction,
            improved_visualisation=bool(
                settings.get("improved_visualisation", False)
            ),
        )

        font_family = settings.get("font_family", "Roboto Mono")
        font_size = settings.get("font_size_default", 14)
        diagram_text = create_overlapped_sequence_view(
            top_line=top_line,
            mid_line=primer_line,
            bottom_line=bottom_line,
            font_family=font_family,
            font_size=font_size,
        )

        if origin.settings.amplify4_compatibility_mode:
            prim_str = f"{int(origin.primability * 100)}%"
            stab_str = f"{int(origin.stability * 100)}%"
        else:
            prim_str = f"{origin.primability:.3f}"
            stab_str = f"{origin.stability:.3f}"

        font_size_small = settings.get("font_size_small", 12)
        body_controls = [
            ft.Row(
                [
                    ft.Container(
                        content=ft.Text(
                            f"Primeability: {prim_str}",
                            weight=ft.FontWeight.BOLD,
                            color=GUIColors.DIAGRAM_BLACK,
                            size=font_size_small,
                        ),
                        bgcolor=GUIColors.SELECTED_ROW_BG,
                        padding=ft.Padding(8, 4, 8, 4),
                        border_radius=4,
                    ),
                    ft.Container(
                        content=ft.Text(
                            f"Stability: {stab_str}",
                            weight=ft.FontWeight.BOLD,
                            color=GUIColors.DIAGRAM_BLACK,
                            size=font_size_small,
                        ),
                        bgcolor=GUIColors.SELECTED_ROW_BG,
                        padding=ft.Padding(8, 4, 8, 4),
                        border_radius=4,
                    ),
                    ft.Container(
                        content=ft.Text(
                            f"Quality: {origin.quality:.4f}",
                            weight=ft.FontWeight.BOLD,
                            color=GUIColors.DIAGRAM_BLACK,
                            size=font_size_small,
                        ),
                        bgcolor=GUIColors.SELECTED_ROW_BG,
                        padding=ft.Padding(8, 4, 8, 4),
                        border_radius=4,
                    ),
                ],
                spacing=8,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
                alignment=ft.MainAxisAlignment.START,
            ),
            ft.Container(
                content=ft.Row(
                    [diagram_text],
                    scroll=ft.ScrollMode.ALWAYS,
                ),
                padding=12,
                border_radius=6,
                border=ft.Border.all(1, GUIColors.OUTLINE_VARIANT),
            ),
        ]

        super().__init__(
            card_id=card_id,
            title=(f"Context Map - {primer_name} ({primer_type}) Binding Site"),
            settings=settings,
            dismiss_callback=dismiss_callback,
            body_controls=body_controls,
        )
