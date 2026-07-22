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

"""Case conversion helpers for DNA template sequence input."""

from __future__ import annotations

from typing import TYPE_CHECKING

import flet as ft

from amplifyp.gui.utils.data_helpers import clean_sequence

if TYPE_CHECKING:
    from .input import TemplateInput


def upper_case_click(template_input: TemplateInput, _e: ft.Event) -> None:
    """Handle upper case button click."""
    change_selection_case(template_input, to_upper=True)


def lower_case_click(template_input: TemplateInput, _e: ft.Event) -> None:
    """Handle lower case button click."""
    change_selection_case(template_input, to_upper=False)


def change_selection_case(
    template_input: TemplateInput, to_upper: bool
) -> None:
    """Convert selected bases in template to upper/lower case."""
    sel = template_input.template_sequence.selection
    if not sel or not sel.is_valid or sel.start == sel.end:
        template_input._show_notification("Please select sequence text first.")
        return

    raw_val = template_input.template_sequence.value or ""
    start, end = sel.start, sel.end

    selected_text = raw_val[start:end]
    modified_text = selected_text.upper() if to_upper else selected_text.lower()
    new_val = raw_val[:start] + modified_text + raw_val[end:]

    template_input.template_sequence.value = new_val
    template_input.input_data.template = clean_sequence(new_val)
    template_input._cleaned_len = len(template_input.input_data.template)

    template_input._update_line_numbers(update=False)
    template_input.template_sequence.selection = ft.TextSelection(
        base_offset=sel.base_offset, extent_offset=sel.extent_offset
    )
    try:
        template_input.template_sequence.update()
        template_input.update()
    except (RuntimeError, AssertionError):
        pass
    template_input.on_change_handler(None)
