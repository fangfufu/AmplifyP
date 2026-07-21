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

"""Sequence formatting and line number gutter helpers for template input."""

from __future__ import annotations

from typing import TYPE_CHECKING

from amplifyp.gui.utils.data_helpers import clean_sequence, format_sequence

if TYPE_CHECKING:
    from .input import TemplateInput


def adjust_wrap_length(
    template_input: TemplateInput, left_width: float, update: bool = True
) -> int:
    """Adjust wrap length based on selected width or Auto.

    Args:
        template_input: The parent TemplateInput component.
        left_width: The available width for the left panel in pixels.
        update: Whether to trigger a UI update after adjusting.

    Returns:
        The computed wrap length (number of bases per line).
    """
    template_input._last_left_width = left_width
    template_input.sequence_layout.width = max(100.0, left_width - 15.0)

    font_size = max(1, template_input.settings.get("font_size_default", 14))
    char_width = font_size * 0.70

    # Calculate dynamic gutter width based on template digits
    template_len = len(template_input.input_data.template)
    max_digits = len(str(max(1, template_len)))
    gutter_width = 20 + max_digits * char_width

    wrap_setting = validate_bases_per_line(template_input)
    if wrap_setting == "Auto":
        available_width = left_width - gutter_width - 100
        max_fit = int(available_width / char_width)
        wrap_length = (max_fit // 10) * 10
        wrap_length = max(10, min(100, wrap_length))
    else:
        wrap_length = wrap_setting if isinstance(wrap_setting, int) else 50

    field_available_width = max(100.0, left_width - gutter_width - 35.0)

    # Disable autowrapping at window edge, enable horizontal scroll
    target_width = max(
        field_available_width + 100.0, wrap_length * char_width + 100.0
    )
    template_input.template_sequence_wrapper.width = target_width
    template_input.template_sequence.width = target_width
    template_input.template_sequence.expand = False

    # Update TextField content with new wrapping
    template_input.template_sequence.value = format_sequence(
        template_input.input_data.template, wrap_length
    )
    update_line_numbers(template_input, update=update)

    try:
        if template_input.page:
            template_input.template_sequence_wrapper.update()
            template_input.template_sequence.update()
            template_input.template_sequence_container.update()
    except (AssertionError, RuntimeError):
        pass

    return wrap_length


def update_line_numbers(
    template_input: TemplateInput,
    update: bool = True,
    gutter_only: bool = False,
) -> None:
    """Update the line numbers gutter based on current template sequence.

    Args:
        template_input: The parent TemplateInput component.
        update: Whether to trigger a UI update after adjusting.
        gutter_only: If True, only update the gutter text and container.
    """
    text = template_input.template_sequence.value or ""
    lines = text.split("\n")
    line_indices = []
    current_idx = 0
    for line in lines:
        line_indices.append(str(current_idx))
        current_idx += len(clean_sequence(line))

    template_input.line_numbers_text.value = "\n".join(line_indices)

    # Set dynamic gutter width
    font_size = max(1, template_input.settings.get("font_size_default", 14))
    char_width = font_size * 0.75
    template_len = current_idx
    max_digits = len(str(max(1, template_len)))
    gutter_width = 20 + max_digits * char_width
    template_input.line_numbers_container.width = gutter_width

    if update:
        try:
            page = template_input.page
        except RuntimeError:
            page = None
        if page:
            try:
                if gutter_only:
                    template_input.line_numbers_text.update()
                    template_input.line_numbers_container.update()
                else:
                    template_input.update()
            except (RuntimeError, AssertionError):
                pass


def validate_bases_per_line(
    template_input: TemplateInput, val_str: str | None = None
) -> int | str | None:
    """Validate bases per line, enforcing 10..100 or Auto.

    Args:
        template_input: The parent TemplateInput component.
        val_str: Optional string value to validate. If None, reads from UI.

    Returns:
        ``"Auto"``, an integer multiple of 10 between 10 and 100, or None.
    """
    if val_str is None:
        val_str = (template_input.bases_per_line_value_text.value or "").strip()
    if val_str.lower() == "auto":
        return "Auto"
    try:
        val_int = int(val_str.strip())
        if 10 <= val_int <= 100 and val_int % 10 == 0:
            return val_int
    except ValueError:
        pass
    return None
