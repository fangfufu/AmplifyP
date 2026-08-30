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

"""QualityBarChart component for 1D primer designer view."""

from collections.abc import Callable, Sequence

import flet as ft

from amplifyp.dimer import PrimerDimer
from amplifyp.gui.colours import GUIColours


class QualityBarChart:
    """Vertical quality bar chart rendering for 1D self-dimers."""

    @staticmethod
    def build_chart(
        dimers: Sequence[PrimerDimer],
        container_height: float,
        on_primer_selected: Callable[[PrimerDimer, int], None],
    ) -> ft.Control:
        """Build vertical quality bar chart for generated self-dimers."""
        if not dimers:
            return ft.Container(
                content=ft.Text(
                    "No analysis results yet. Enter sequence and "
                    "click Analyse.",
                    size=12,
                    italic=True,
                    color=GUIColours.TEXT_ON_SURFACE,
                ),
                alignment=ft.Alignment(0, 0),
                expand=True,
            )

        max_quality = max((d.quality for d in dimers), default=1.0)
        max_quality = max(max_quality, 1.0)

        max_bar_h = max(10.0, container_height - 98.0)
        bar_stack_height = max_bar_h + 20.0

        bar_controls: list[ft.Control] = []

        for step_idx, dimer in enumerate(dimers):
            primer_seq = dimer.primer_1.seq
            primer_len = len(primer_seq)
            q_val = dimer.quality

            bar_h = max(6.0, (q_val / max_quality) * max_bar_h)

            bar_with_score = ft.Column(
                [
                    ft.Text(
                        f"{round(q_val)}",
                        size=11,
                        weight=ft.FontWeight.BOLD,
                    ),
                    ft.Container(
                        width=22,
                        height=bar_h,
                        bgcolor=GUIColours.PRIMARY,
                        border_radius=4,
                        tooltip=(
                            f"Step {step_idx + 1}: {primer_len} nt\n"
                            f"Quality: {round(q_val)}"
                        ),
                    ),
                ],
                height=bar_stack_height,
                alignment=ft.MainAxisAlignment.END,
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=2,
            )

            bar_column = ft.Column(
                [
                    bar_with_score,
                    ft.Text(f"{primer_len} nt", size=11),
                ],
                alignment=ft.MainAxisAlignment.END,
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=2,
            )

            bar_item = ft.Container(
                content=bar_column,
                padding=ft.Padding(6, 4, 6, 26),
                ink=True,
                on_click=lambda _ev, d=dimer, idx=step_idx: on_primer_selected(
                    d, idx
                ),
            )
            bar_controls.append(bar_item)

        return ft.Row(
            bar_controls,
            scroll=ft.ScrollMode.ALWAYS,
            alignment=ft.MainAxisAlignment.START,
            vertical_alignment=ft.CrossAxisAlignment.END,
            spacing=4,
        )
