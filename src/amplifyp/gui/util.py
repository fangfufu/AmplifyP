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

"""Utility functions for sequence handling and state serialization."""

import flet as ft
import flet.canvas as cv
import yaml


def clean_sequence(seq: str) -> str:
    """Clean sequence of escaped and standard whitespaces."""
    if not seq:
        return ""
    clean = str(seq).replace("\\n", "").replace("\\t", "").replace("\\r", "")
    return "".join(clean.split()).upper()


def format_sequence(seq: str, wrap_length: int = 80) -> str:
    """Format sequence into lines of specified length."""
    clean = clean_sequence(seq)
    return "\n".join(
        [clean[i : i + wrap_length] for i in range(0, len(clean), wrap_length)]
    )


def serialize_state(state: dict[str, object]) -> str:
    """Serialize state dict to YAML string, handling multiline strings."""

    def multiline_presenter(dumper: yaml.Dumper, data: str) -> yaml.ScalarNode:
        if "\n" in data:
            return dumper.represent_scalar(
                "tag:yaml.org,2002:str", data, style="|"
            )
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    class _StateDumper(yaml.Dumper):
        pass

    _StateDumper.add_representer(str, multiline_presenter)
    return yaml.dump(state, Dumper=_StateDumper, sort_keys=False)


def create_overlapped_sequence_view(
    top_line: str,
    mid_line: str,
    bottom_line: str,
) -> ft.Stack:
    """Create a Flet Stack showing visually aligned sequences.

    Draws character-by-character on a Canvas to prevent rendering offset
    issues, with a transparent selectable text layer behind it for copy-paste.
    """
    # Ensure all three lines have the same total width
    target_width = max(len(top_line), len(mid_line), len(bottom_line))
    top_line = top_line.ljust(target_width)
    mid_line = mid_line.ljust(target_width)
    bottom_line = bottom_line.ljust(target_width)

    # Build Flet Canvas shapes for robust character alignment
    shapes = []
    char_width = 11.5  # Width of each character cell in pixels
    x_offset = 12  # Left margin padding
    y_top = 12  # Y coord for top line
    y_mid = 38  # Y coord for middle line
    y_bottom = 64  # Y coord for bottom line

    # Draw each line character-by-character
    for idx, char in enumerate(top_line):
        if char != " ":
            shapes.append(
                cv.Text(
                    x_offset + idx * char_width,
                    y_top,
                    char,
                    style=ft.TextStyle(
                        font_family="Roboto Mono",
                        size=14,
                        weight=ft.FontWeight.BOLD,
                        color=ft.Colors.ON_SURFACE,
                    ),
                )
            )

    for idx, char in enumerate(mid_line):
        if char != " ":
            shapes.append(
                cv.Text(
                    x_offset + idx * char_width,
                    y_mid,
                    char,
                    style=ft.TextStyle(
                        font_family="Roboto Mono",
                        size=14,
                        weight=ft.FontWeight.BOLD,
                        color=ft.Colors.GREEN_400,
                    ),
                )
            )

    for idx, char in enumerate(bottom_line):
        if char != " ":
            shapes.append(
                cv.Text(
                    x_offset + idx * char_width,
                    y_bottom,
                    char,
                    style=ft.TextStyle(
                        font_family="Roboto Mono",
                        size=14,
                        weight=ft.FontWeight.BOLD,
                        color=ft.Colors.ON_SURFACE,
                    ),
                )
            )

    # Dynamic width based on the number of characters
    canvas_width = target_width * char_width + x_offset * 2
    diagram_canvas = cv.Canvas(
        shapes=shapes,
        height=90,
        width=canvas_width,
    )
    diagram_canvas.pointer_events = "none"

    # Create a selectable text layer behind the canvas
    # that perfectly matches the character spacing
    selectable_text = ft.Text(
        f"{top_line}\n{mid_line}\n{bottom_line}",
        font_family="Roboto Mono",
        size=14,
        selectable=True,
        color=ft.Colors.TRANSPARENT,
    )
    selectable_text.line_height = 1.85
    selectable_text.left = x_offset - 2

    return ft.Stack(
        [
            diagram_canvas,
            selectable_text,
        ],
        width=canvas_width,
        height=90,
    )
