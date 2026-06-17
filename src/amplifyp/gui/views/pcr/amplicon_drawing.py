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

"""Amplicon drawing class and detail card helpers for the PCRView."""

from collections.abc import Callable
from typing import Any

import flet as ft
import flet.canvas as cv

from amplifyp.gui.colors import GUIColors

from .dismissible_detail_card import DismissibleDetailCard


class DrawnAmplicon:
    """Amplicon bar drawn on the PCR diagram canvas."""

    def __init__(
        self,
        amp: Any,
        idx: int,
        target_length: int,
        t_width: float,
        h_margin: float,
        v_target: float,
        c_width: float,
        settings: Any,
        on_click: Callable[[], None],
        v_frag_start: float | None = None,
    ) -> None:
        """Initialize the DrawnAmplicon."""
        self.amp = amp
        self.idx = idx
        self.target_length = target_length
        self.t_width = t_width
        self.h_margin = h_margin
        self.v_target = v_target
        self.c_width = c_width
        self.settings = settings
        self.on_click = on_click

        self.start_idx = amp.start.index
        self.end_idx = amp.end.index

        # Calculate coordinates
        self.x_start = (
            h_margin + (self.start_idx / target_length * t_width)
            if target_length
            else h_margin
        )
        self.x_end = (
            h_margin + (self.end_idx / target_length * t_width)
            if target_length
            else (c_width - h_margin)
        )

        if v_frag_start is None:
            v_frag_start = v_target + 40.0
        v_frag_step = 35
        self.y_pos = v_frag_start + (idx * v_frag_step)

        # Dynamic thickness based on the quality score
        if amp.q_score < 300:
            self.bar_height = 8.0
        elif amp.q_score < 700:
            self.bar_height = 5.5
        elif amp.q_score < 1500:
            self.bar_height = 3.5
        elif amp.q_score < 4000:
            self.bar_height = 2.0
        else:
            self.bar_height = 1.0

    def draw(self, canvas: cv.Canvas, stack: ft.Stack) -> None:
        """Draw the amplicon bar, text label, and gesture detector overlay."""
        # Amplicon Bar (Black, filling path)
        amp_paint = ft.Paint(
            color=GUIColors.DIAGRAM_BLACK,
            style=ft.PaintingStyle.FILL,
        )
        if self.amp.circular:
            canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(self.h_margin, self.y_pos),
                        cv.Path.LineTo(
                            self.c_width - self.h_margin, self.y_pos
                        ),
                        cv.Path.LineTo(
                            self.c_width - self.h_margin,
                            self.y_pos + self.bar_height,
                        ),
                        cv.Path.LineTo(
                            self.h_margin, self.y_pos + self.bar_height
                        ),
                        cv.Path.Close(),
                    ],
                    paint=amp_paint,
                )
            )
            label_x = self.h_margin + (self.t_width / 2.0)
        else:
            canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(self.x_start, self.y_pos),
                        cv.Path.LineTo(self.x_end, self.y_pos),
                        cv.Path.LineTo(
                            self.x_end, self.y_pos + self.bar_height
                        ),
                        cv.Path.LineTo(
                            self.x_start, self.y_pos + self.bar_height
                        ),
                        cv.Path.Close(),
                    ],
                    paint=amp_paint,
                )
            )
            label_x = self.x_start + ((self.x_end - self.x_start) / 2.0)

        # Amplicon Length Text (just number, black, centered)
        canvas.shapes.append(
            cv.Text(
                label_x,
                self.y_pos + self.bar_height + 5,
                str(len(self.amp.product)),
                style=ft.TextStyle(
                    size=self.settings.get("font_size_map_amplicon", 13),
                    color=GUIColors.DIAGRAM_BLACK,
                ),
                alignment=ft.Alignment(0.0, -1.0),
            )
        )

        # Click overlay for the amplicon
        amp_width = max(
            10.0,
            (self.x_end - self.x_start)
            if not self.amp.circular
            else (self.c_width - 2 * self.h_margin),
        )
        amp_left = self.x_start if not self.amp.circular else self.h_margin
        stack.controls.append(
            ft.GestureDetector(
                mouse_cursor=ft.MouseCursor.CLICK,
                on_tap=lambda _: self.on_click(),
                content=ft.Container(
                    bgcolor=GUIColors.TRANSPARENT,
                    width=amp_width,
                    height=20 + self.bar_height,
                ),
                left=amp_left,
                top=self.y_pos - 3,
            )
        )


class AmpliconDetailCard(DismissibleDetailCard):
    """Card displaying details of a single amplicon."""

    def __init__(
        self,
        amp: Any,
        settings: Any,
        dismiss_callback: Callable[[ft.Card], None],
    ) -> None:
        """Initialize the AmpliconDetailCard."""
        card_id = (
            f"amplicon_{amp.fwd_origin.name}_{amp.rev_origin.name}_"
            f"{amp.start.index}_{amp.end.index}"
        )

        full_seq = str(amp.product.seq)
        fwd_len = len(amp.fwd_origin.seq)
        rev_len = len(amp.rev_origin.seq)
        mid_len = len(amp.product) - (fwd_len + rev_len)

        if len(full_seq) >= fwd_len + rev_len:
            fwd_part = full_seq[:fwd_len]
            mid_part = full_seq[fwd_len:-rev_len]
            rev_part = full_seq[-rev_len:]
        else:
            fwd_part = full_seq
            mid_part = ""
            rev_part = ""

        from amplifyp.gui.util import _resolve_font_family

        font_family = settings.get("font_family", "Roboto Mono")
        resolved = _resolve_font_family(font_family)
        sequence_text = ft.Text(
            spans=[
                ft.TextSpan(
                    fwd_part,
                    style=ft.TextStyle(
                        color=GUIColors.FWD_PRIMER,
                        weight=ft.FontWeight.BOLD,
                    ),
                ),
                ft.TextSpan(
                    mid_part,
                    style=ft.TextStyle(
                        color=GUIColors.TEXT_ON_SURFACE,
                    ),
                ),
                ft.TextSpan(
                    rev_part,
                    style=ft.TextStyle(
                        color=GUIColors.REV_LABEL,
                        weight=ft.FontWeight.BOLD,
                    ),
                ),
            ],
            font_family=resolved,
            size=settings.get("font_size_body", 13),
            selectable=True,
        )

        body_controls = [
            ft.Text(
                spans=[
                    ft.TextSpan("▶(primer: "),
                    ft.TextSpan(
                        amp.fwd_origin.name,
                        style=ft.TextStyle(
                            color=GUIColors.FWD_PRIMER,
                            weight=ft.FontWeight.BOLD,
                        ),
                    ),
                    ft.TextSpan(f") — {mid_len} bp — (primer: "),
                    ft.TextSpan(
                        amp.rev_origin.name,
                        style=ft.TextStyle(
                            color=GUIColors.REV_LABEL,
                            weight=ft.FontWeight.BOLD,
                        ),
                    ),
                    ft.TextSpan(
                        f")◀      Q = {amp.q_score:.1f} "
                        f"({amp.q_score_report_str(verbose=True)})"
                    ),
                ],
                selectable=True,
                size=settings.get("font_size_body", 13),
            ),
            ft.Text(
                "Amplified Sequence:",
                weight=ft.FontWeight.BOLD,
                selectable=True,
            ),
            ft.Container(
                content=sequence_text,
                padding=12,
                border_radius=6,
                border=ft.Border.all(1, GUIColors.OUTLINE_VARIANT),
            ),
        ]

        super().__init__(
            card_id=card_id,
            title=f"Amplicon: {len(amp.product)} bp",
            settings=settings,
            dismiss_callback=dismiss_callback,
            body_controls=body_controls,
        )
