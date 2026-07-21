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

"""Layout manager for PrimerInput component."""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import flet as ft

from .row import PrimerRow

if TYPE_CHECKING:
    from .input import PrimerInput


class PrimerLayoutManager:
    """Handles layout math, column resizing, and row visibility caching."""

    def __init__(self, owner: PrimerInput) -> None:
        """Initialise the PrimerLayoutManager."""
        self.owner = owner

    def get_panel_width(self) -> float:
        """Get the current width of the primer input panel."""
        page_width = (
            self.owner.app_page.width
            if (self.owner.app_page and self.owner.app_page.width)
            else 800.0
        )
        parent_view = getattr(self.owner.on_change_handler, "__self__", None)
        right_fraction = (
            getattr(parent_view, "right_fraction", 0.5) if parent_view else 0.5
        )
        return page_width * right_fraction

    def _cache_visible_rows_if_needed(self) -> None:
        """Cache the list of visible rows based on viewport scroll state."""
        if self.owner._visible_rows_cache is not None:
            return

        self.owner._visible_rows_cache = []
        scroll_y = self.owner.primers_list.scroll_pixels
        viewport_h = self.owner.primers_list.viewport_dimension
        if self.owner.app_page and self.owner.app_page.height:
            viewport_h = max(viewport_h, float(self.owner.app_page.height))
        current_y = 0.0

        for row in self.owner.primers_list.controls:
            if isinstance(row, PrimerRow):
                row_h = (
                    30.0
                    if not (row.name_field.error or row.seq_field.error)
                    else 55.0
                )
                row_top = current_y
                row_bottom = current_y + row_h

                if (row_bottom >= scroll_y - 60.0) and (
                    row_top <= scroll_y + viewport_h + 60.0
                ):
                    self.owner._visible_rows_cache.append(row)

                current_y += row_h

    def on_primer_divider_pan(self, e: ft.DragUpdateEvent) -> None:
        """Handle dragging the vertical divider between Name and Sequence."""
        delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
        panel_width = self.get_panel_width()

        show_temp = self.owner.settings.get("show_primer_temperature", False)
        extra_space = 54.0 if show_temp else 0.0
        max_name_width = max(80.0, panel_width - 271.0 - extra_space)

        target_width = max(
            80.0,
            min(
                max_name_width,
                self.owner.name_column_width + delta_x,
            ),
        )

        self.owner.name_column_width = target_width
        self.owner.name_column_ratio = (
            self.owner.name_column_width / panel_width
        )

        cast(
            ft.Container, self.owner.primers_header.controls[3]
        ).width = self.owner.name_column_width
        self.owner.primer_header.update()

        self._cache_visible_rows_if_needed()
        visible_rows = self.owner._visible_rows_cache or []

        for row in visible_rows:
            row.name_scroll.width = self.owner.name_column_width
        self.owner.primers_list.update()

    def on_primer_divider_pan_end(self, _e: ft.DragEndEvent) -> None:
        """Handle finishing the drag of the vertical divider."""
        self.owner._visible_rows_cache = None

        cast(
            ft.Container, self.owner.primers_header.controls[3]
        ).width = self.owner.name_column_width
        self.owner.primer_header.update()

        for row in self.owner.primers_list.controls:
            if isinstance(row, PrimerRow):
                row.name_scroll.width = self.owner.name_column_width
        self.owner.primers_list.update()

    def adjust_name_column_width(
        self, new_panel_width: float, during_drag: bool = False
    ) -> None:
        """Scale name column width proportionally based on panel width."""
        show_temp = self.owner.settings.get("show_primer_temperature", False)
        extra_space = 54.0 if show_temp else 0.0
        max_name_width = max(80.0, new_panel_width - 271.0 - extra_space)
        target_width = max(
            80.0,
            min(
                max_name_width,
                new_panel_width * self.owner.name_column_ratio,
            ),
        )

        self.owner.name_column_width = target_width

        if (
            hasattr(self.owner, "primers_header")
            and self.owner.primers_header
            and len(self.owner.primers_header.controls) > 3
        ):
            cast(
                ft.Container, self.owner.primers_header.controls[3]
            ).width = self.owner.name_column_width
            self.owner.primer_header.update()

        if during_drag:
            self._cache_visible_rows_if_needed()
            visible_rows = self.owner._visible_rows_cache or []

            for row in visible_rows:
                row.name_scroll.width = self.owner.name_column_width
            self.owner.primers_list.update()
        else:
            self.owner._visible_rows_cache = None
            for row in self.owner.primers_list.controls:
                if isinstance(row, PrimerRow):
                    row.name_scroll.width = self.owner.name_column_width
            self.owner.primers_list.update()
