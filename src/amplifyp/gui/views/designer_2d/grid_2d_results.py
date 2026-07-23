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

"""Grid component for displaying 2D quality and overlap results matrix."""

from collections.abc import Callable

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.primer_designer_2d import (
    FilterMetric,
    PrimerDesigner2D,
    PrimerDimers2D,
)


class Grid2DResultsView(ft.Container):  # type: ignore[misc]
    """Grid matrix displaying 2D quality and overlap results."""

    def __init__(
        self,
        settings: GUISettings,
        on_select_step_callback: Callable[[PrimerDimers2D], None],
    ) -> None:
        """Initialise the Grid2DResultsView."""
        super().__init__(expand=True)
        self.settings = settings
        self.on_select_step_callback = on_select_step_callback

        self._selected_step: PrimerDimers2D | None = None
        self._cell_containers: dict[tuple[int, int], ft.Container] = {}

        self.content_column = ft.Column(
            [
                ft.Text(
                    "2D Truncation Results Grid",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_subheader", 16),
                ),
                ft.Container(
                    content=ft.Text(
                        "Run 2D analysis to view quality and overlap matrix.",
                        italic=True,
                        size=self.settings.get("font_size_small", 12),
                        color=GUIColours.TEXT_ON_SURFACE,
                    ),
                    expand=True,
                    alignment=ft.Alignment(0, 0),
                ),
            ],
            spacing=6,
            expand=True,
            scroll=ft.ScrollMode.ALWAYS,
        )
        self.content = self.content_column

    def update_grid(self, designer: PrimerDesigner2D) -> None:
        """Populate and render the 2D matrix grid.

        Args:
            designer: Computed PrimerDesigner2D analysis instance.
        """
        self._selected_step = None
        self._cell_containers.clear()

        steps = designer.all_steps
        if not steps:
            self.content_column.controls = [
                ft.Text(
                    "2D Truncation Results Grid",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_subheader", 16),
                ),
                ft.Container(
                    content=ft.Text(
                        "No valid 2D truncation combinations matched "
                        "the active filters.",
                        italic=True,
                        size=self.settings.get("font_size_small", 12),
                        color=GUIColours.ERROR_RED,
                    ),
                    expand=True,
                    alignment=ft.Alignment(0, 0),
                ),
            ]
            try:
                if self.page:
                    self.page.update()
            except RuntimeError:
                pass
            return

        # Map steps by (fwd_len, rev_len)
        step_map: dict[tuple[int, int], PrimerDimers2D] = {}
        fwd_lengths_set: set[int] = set()
        rev_lengths_set: set[int] = set()

        for step in steps:
            fwd_len = len(step.fwd_fwd.primer_1.seq)
            rev_len = len(step.rev_rev.primer_1.seq)
            step_map[(fwd_len, rev_len)] = step
            fwd_lengths_set.add(fwd_len)
            rev_lengths_set.add(rev_len)

        fwd_lengths = sorted(fwd_lengths_set, reverse=True)
        rev_lengths = sorted(rev_lengths_set, reverse=True)

        use_max = designer.filter_metric == FilterMetric.MAX
        font_small = self.settings.get("font_size_small", 11)

        # Header Row
        header_cells: list[ft.Control] = [
            ft.Container(
                content=ft.Text(
                    "Rev \\ Fwd",
                    weight=ft.FontWeight.BOLD,
                    size=font_small,
                    color=GUIColours.PRIMARY,
                ),
                width=80,
                height=36,
                alignment=ft.Alignment(0, 0),
                bgcolor=GUIColours.SURFACE_VARIANT,
                border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
                border_radius=4,
            )
        ]

        for f_len in fwd_lengths:
            header_cells.append(
                ft.Container(
                    content=ft.Text(
                        f"{f_len} bp",
                        weight=ft.FontWeight.BOLD,
                        size=font_small,
                    ),
                    width=72,
                    height=36,
                    alignment=ft.Alignment(0, 0),
                    bgcolor=GUIColours.SURFACE_VARIANT,
                    border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
                    border_radius=4,
                )
            )

        grid_rows: list[ft.Control] = [ft.Row(header_cells, spacing=4)]

        # Grid Data Rows
        for r_len in rev_lengths:
            row_cells: list[ft.Control] = [
                ft.Container(
                    content=ft.Text(
                        f"{r_len} bp",
                        weight=ft.FontWeight.BOLD,
                        size=font_small,
                    ),
                    width=80,
                    height=46,
                    alignment=ft.Alignment(0, 0),
                    bgcolor=GUIColours.SURFACE_VARIANT,
                    border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
                    border_radius=4,
                )
            ]

            for f_len in fwd_lengths:
                step = step_map.get((f_len, r_len))
                if step is None:
                    row_cells.append(
                        ft.Container(
                            content=ft.Text(
                                "N/A",
                                size=font_small,
                                color=GUIColours.TEXT_ON_SURFACE,
                                opacity=0.5,
                            ),
                            width=72,
                            height=46,
                            alignment=ft.Alignment(0, 0),
                            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
                            border_radius=4,
                        )
                    )
                    continue

                q_val = step.max_quality if use_max else step.mean_quality
                o_val = step.max_overlap if use_max else step.mean_overlap
                o_str = f"{o_val}" if use_max else f"{o_val:.1f}"

                # Top half: Quality, Bottom half: Overlap
                cell_content = ft.Column(
                    [
                        ft.Container(
                            content=ft.Text(
                                f"Q: {q_val:.1f}",
                                weight=ft.FontWeight.BOLD,
                                size=font_small,
                                color=GUIColours.PRIMARY,
                            ),
                            alignment=ft.Alignment(0, 0),
                            expand=True,
                        ),
                        ft.Divider(height=1, color=GUIColours.OUTLINE_VARIANT),
                        ft.Container(
                            content=ft.Text(
                                f"O: {o_str}",
                                size=font_small - 1,
                                color=GUIColours.DIAGRAM_BLACK,
                            ),
                            alignment=ft.Alignment(0, 0),
                            expand=True,
                        ),
                    ],
                    spacing=0,
                    alignment=ft.MainAxisAlignment.CENTER,
                )

                def _make_click_handler(
                    st: PrimerDimers2D, key: tuple[int, int]
                ) -> Callable[[ft.ControlEvent], None]:
                    return lambda e: self._on_cell_click(st, key)

                cell_container = ft.Container(
                    content=cell_content,
                    width=72,
                    height=46,
                    padding=2,
                    border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
                    border_radius=4,
                    ink=True,
                    on_click=_make_click_handler(step, (f_len, r_len)),
                    tooltip=(
                        f"Forward: {f_len} bp | Reverse: {r_len} bp\n"
                        f"{'Max' if use_max else 'Mean'} Quality: {q_val:.1f}\n"
                        f"{'Max' if use_max else 'Mean'} Overlap: {o_str} bp"
                    ),
                )
                self._cell_containers[(f_len, r_len)] = cell_container
                row_cells.append(cell_container)

            grid_rows.append(ft.Row(row_cells, spacing=4))

        # Metric Legend Header
        metric_label = "Max" if use_max else "Mean"
        self.content_column.controls = [
            ft.Row(
                [
                    ft.Text(
                        "2D Truncation Results Grid",
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_subheader", 16),
                    ),
                    ft.Container(
                        content=ft.Text(
                            f"Metric: {metric_label} Quality / Overlap",
                            size=font_small,
                            weight=ft.FontWeight.BOLD,
                            color=GUIColours.PRIMARY,
                        ),
                        bgcolor=GUIColours.SELECTED_ROW_BG,
                        padding=ft.Padding(6, 2, 6, 2),
                        border_radius=4,
                    ),
                ],
                alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            ),
            ft.Row(
                [
                    ft.Column(
                        grid_rows,
                        spacing=4,
                        scroll=ft.ScrollMode.ALWAYS,
                    )
                ],
                scroll=ft.ScrollMode.ALWAYS,
                expand=True,
            ),
        ]

        try:
            if self.page:
                self.page.update()
        except RuntimeError:
            pass

    def _on_cell_click(
        self, step: PrimerDimers2D, key: tuple[int, int]
    ) -> None:
        """Handle selection of a grid cell.

        Args:
            step: Selected PrimerDimers2D instance.
            key: (fwd_len, rev_len) tuple key for the clicked cell container.
        """
        # Highlight selected cell container
        for cell_key, container in self._cell_containers.items():
            if cell_key == key:
                container.border = ft.Border.all(2, GUIColours.PRIMARY)
                container.bgcolor = GUIColours.SELECTED_ROW_BG
            else:
                container.border = ft.Border.all(1, GUIColours.OUTLINE_VARIANT)
                container.bgcolor = None

        self._selected_step = step
        self.on_select_step_callback(step)

        try:
            if self.page:
                self.page.update()
        except RuntimeError:
            pass
