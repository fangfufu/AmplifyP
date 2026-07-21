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

"""Layout and resizing handlers for DNA template and primers input views."""

from typing import Any

import flet as ft  # type: ignore[import-not-found, unused-ignore]


def on_pan_update(input_view: Any, e: ft.DragUpdateEvent) -> None:
    """Handle resizing the bottom (right) container via the divider."""
    page_width = input_view.app_page.width
    if isinstance(page_width, (int, float)) and page_width > 0:
        delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
        current_width = page_width * input_view.right_fraction
        new_width = max(200.0, current_width - delta_x)
        new_width = min(new_width, page_width - 200.0)

        input_view.right_fraction = new_width / page_width

        input_view.template_input.expand = None
        left_width = page_width - new_width - 5.0
        input_view.template_input.width = left_width
        input_view.template_input.update()

        input_view.primer_input.layout_manager.adjust_name_column_width(
            new_width, during_drag=True
        )


def on_pan_end(input_view: Any, e: ft.DragEndEvent) -> None:
    """Handle finishing the drag of the main layout divider."""
    page_width = input_view.app_page.width
    if isinstance(page_width, (int, float)) and page_width > 0:
        panel_width = page_width * input_view.right_fraction

        input_view.template_input.width = None
        input_view.template_input.expand = int(
            (1.0 - input_view.right_fraction) * 1000
        )
        input_view.primer_input.expand = int(input_view.right_fraction * 1000)

        input_view.primer_input.layout_manager.adjust_name_column_width(
            panel_width, during_drag=False
        )
        input_view._adjust_template_wrap(update_first=False)
        input_view.update()


def handle_resize(input_view: Any, e: ft.PageResizeEvent) -> None:
    """Handle page resize to proportionally scale name column."""
    page_width = input_view.app_page.width
    if isinstance(page_width, (int, float)) and page_width > 0:
        panel_width = page_width * input_view.right_fraction
        input_view.primer_input.layout_manager.adjust_name_column_width(
            panel_width
        )
        input_view._adjust_template_wrap(update_first=False)
        input_view.update()


def adjust_template_wrap(input_view: Any, update_first: bool = True) -> None:
    """Adjust the template wrap length based on the available width."""
    page_width = input_view.app_page.width
    if isinstance(page_width, (int, float)) and page_width > 0:
        left_width = page_width * (1.0 - input_view.right_fraction)
        if update_first:
            input_view.template_input.adjust_wrap_length(
                left_width, update=True
            )
        else:
            input_view.template_input.adjust_wrap_length(
                left_width, update=False
            )
