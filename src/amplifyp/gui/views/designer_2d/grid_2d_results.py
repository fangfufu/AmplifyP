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

from amplifyp.gui.colours import (
    GUIColours,
    designer_2d_colour,
    get_text_contrast_colour,
)
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
        self._cell_bg_colours: dict[tuple[int, int], str | None] = {}
        self._best_cell_keys: set[tuple[int, int]] = set()

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
            spacing=2,
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
        self._cell_bg_colours.clear()
        self._best_cell_keys.clear()

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
        scheme = self.settings.get("designer_2d_colour_scheme", "None")

        qualities = [
            (s.max_quality if use_max else s.mean_quality) for s in steps
        ]
        min_q = min(qualities) if qualities else 0.0
        max_q = max(qualities) if qualities else 1.0

        # Header Row
        header_cells: list[ft.Control] = [
            ft.Container(
                content=ft.Text(
                    "Rev \\ Fwd",
                    weight=ft.FontWeight.BOLD,
                    size=font_small,
                    color=GUIColours.PRIMARY,
                ),
                width=76,
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
                    width=52,
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
                    width=76,
                    height=36,
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
                            width=52,
                            height=36,
                            alignment=ft.Alignment(0, 0),
                            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
                            border_radius=4,
                        )
                    )
                    continue

                q_val = step.max_quality if use_max else step.mean_quality
                o_val = step.max_overlap if use_max else step.mean_overlap
                o_str = f"{o_val}" if use_max else f"{o_val:.1f}"

                is_best = abs(q_val - min_q) < 1e-6
                if is_best:
                    self._best_cell_keys.add((f_len, r_len))

                bg_col = designer_2d_colour(q_val, min_q, max_q, scheme)
                self._cell_bg_colours[(f_len, r_len)] = bg_col
                text_col = (
                    get_text_contrast_colour(bg_col)
                    if bg_col is not None
                    else GUIColours.PRIMARY
                )

                cell_content = ft.Text(
                    f"{q_val:.1f}",
                    weight=ft.FontWeight.BOLD,
                    size=font_small,
                    color=text_col,
                    text_align=ft.TextAlign.CENTER,
                )

                def _make_click_handler(
                    st: PrimerDimers2D, key: tuple[int, int]
                ) -> Callable[[ft.ControlEvent], None]:
                    return lambda e: self._on_cell_click(st, key)

                border = (
                    ft.Border.all(2, GUIColours.SUCCESS_GREEN)
                    if is_best
                    else ft.Border.all(1, GUIColours.OUTLINE_VARIANT)
                )

                tooltip_prefix = (
                    "★ Best Quality (Lowest Score)\n" if is_best else ""
                )

                cell_container = ft.Container(
                    content=cell_content,
                    width=52,
                    height=36,
                    alignment=ft.Alignment(0, 0),
                    padding=2,
                    bgcolor=bg_col,
                    border=border,
                    border_radius=4,
                    ink=True,
                    on_click=_make_click_handler(step, (f_len, r_len)),
                    tooltip=(
                        f"{tooltip_prefix}"
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
        header_badges: list[ft.Control] = [
            ft.Container(
                content=ft.Text(
                    f"Metric: {metric_label} Quality",
                    size=font_small,
                    weight=ft.FontWeight.BOLD,
                    color=GUIColours.PRIMARY,
                ),
                bgcolor=GUIColours.SELECTED_ROW_BG,
                padding=ft.Padding(8, 2, 12, 2),
                border_radius=4,
            ),
            ft.Container(
                content=ft.Text(
                    f"★ Best Quality: {min_q:.1f}",
                    size=font_small,
                    weight=ft.FontWeight.BOLD,
                    color=GUIColours.SUCCESS_GREEN,
                ),
                bgcolor=GUIColours.SELECTED_ROW_BG,
                padding=ft.Padding(8, 2, 12, 2),
                border_radius=4,
            ),
        ]
        if scheme != "None":
            header_badges.append(
                ft.Container(
                    content=ft.Text(
                        f"Colour Map: {scheme} ({max_q:.1f} - {min_q:.1f})",
                        size=font_small,
                        weight=ft.FontWeight.BOLD,
                        color=GUIColours.PRIMARY,
                    ),
                    bgcolor=GUIColours.SELECTED_ROW_BG,
                    padding=ft.Padding(8, 2, 12, 2),
                    border_radius=4,
                )
            )

        self.content_column.controls = [
            ft.Row(
                [
                    ft.Text(
                        "2D Truncation Results Grid",
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_subheader", 16),
                    ),
                    ft.Container(
                        content=ft.Row(header_badges, spacing=6),
                        padding=ft.Padding(0, 0, 12, 0),
                    ),
                ],
                alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            ),
            ft.Container(
                content=ft.Row(
                    [
                        ft.Container(
                            content=ft.Column(
                                grid_rows,
                                spacing=4,
                                scroll=ft.ScrollMode.ALWAYS,
                            ),
                            padding=ft.Padding(0, 12, 0, 0),
                        ),
                    ],
                    scroll=ft.Scrollbar(
                        orientation=ft.ScrollbarOrientation.TOP,
                        thumb_visibility=True,
                    ),
                    expand=True,
                ),
                padding=ft.Padding(0, 0, 0, 0),
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
            bg_col = self._cell_bg_colours.get(cell_key)
            is_best = cell_key in self._best_cell_keys
            if cell_key == key:
                container.border = ft.Border.all(3, GUIColours.PRIMARY)
                container.bgcolor = (
                    bg_col if bg_col is not None else GUIColours.SELECTED_ROW_BG
                )
            else:
                container.border = (
                    ft.Border.all(2, GUIColours.SUCCESS_GREEN)
                    if is_best
                    else ft.Border.all(1, GUIColours.OUTLINE_VARIANT)
                )
                container.bgcolor = bg_col

        self._selected_step = step
        self.on_select_step_callback(step)

        try:
            if self.page:
                self.page.update()
        except RuntimeError:
            pass
