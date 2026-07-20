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

"""Sequence handling utility functions for GUI."""

import flet as ft


def clean_sequence(seq: str) -> str:
    r"""Clean sequence of escaped and standard whitespaces.

    Removes escaped newlines, tabs, carriage returns, and all standard
    whitespace characters from the sequence.

    Args:
        seq: The sequence string to clean.

    Returns:
        The cleaned sequence with all whitespace removed.
    """
    if not seq:
        return ""
    clean = str(seq).replace("\\n", "").replace("\\t", "").replace("\\r", "")
    return "".join(clean.split())


# Font families that are guaranteed to exist on all platforms.
_SYSTEM_MONOSPACE_FONTS = (
    "monospace",
    "Courier New",
    "Courier",
    "Consolas",
    "Lucida Console",
    "DejaVu Sans Mono",
    "Ubuntu Mono",
    "Liberation Mono",
)


def _resolve_font_family(font_family: str) -> str:
    """Return *font_family* if it is a known system font, else "Roboto Mono".

    In Flet desktop mode the ``page.fonts`` mapping only covers custom
    font-files that the application ships.  When a font name that is *not*
    registered with ``page.fonts`` is used inside an ``ft.Text`` the text
    (and its ``TextSpans``) can fail to render silently.  This helper
    guards against that by falling back to the default ``"Roboto Mono"``
    family which is shipped with this app.
    """
    if not font_family:
        return "Roboto Mono"
    ff = font_family.strip()
    if ff.lower() in {f.lower() for f in _SYSTEM_MONOSPACE_FONTS}:
        return ff
    if ff.lower() == "roboto mono":
        return "Roboto Mono"
    return ff


def format_sequence(seq: str, wrap_length: int = 80) -> str:
    """Format sequence into lines of specified length.

    First cleans the sequence, then wraps it into lines of the given
    length.

    Args:
        seq: The sequence string to format.
        wrap_length: Maximum number of characters per line.

    Returns:
        The formatted sequence with newlines at wrap boundaries.
    """
    clean = clean_sequence(seq)
    return "\n".join(
        [clean[i : i + wrap_length] for i in range(0, len(clean), wrap_length)]
    )


def create_overlapped_sequence_view(
    top_line: str,
    mid_line: str,
    bottom_line: str,
    font_family: str = "Roboto Mono",
    font_size: int = 14,
    is_dimer: bool = False,
    top_name_line: str = "",
    bottom_name_line: str = "",
) -> ft.Text:
    """Create a Flet Text control showing visually aligned sequences.

    Uses TextSpans for the visual representation with colour-coded
    lines. For dimer view, displays five lines (top primer name, top
    sequence, middle bonds, bottom sequence, bottom primer name).
    For context map view, displays five lines (coordinates, primer row,
    bonds, complement, template).

    Args:
        top_line: The top line content (coordinates/arrows for context map,
            or top sequence for dimer view).
        mid_line: The middle line content (primer row for context map,
            or bonds line for dimer view).
        bottom_line: The bottom line content (bonds + template for context
            map, or bottom sequence for dimer view).
        font_family: The font family to use for the text.
        font_size: The font size in pixels.
        is_dimer: Whether this is a dimer view (True) or context map (False).
        top_name_line: The top primer name for dimer view.
        bottom_name_line: The bottom primer name for dimer view.

    Returns:
        A Flet Text control with styled TextSpans.
    """
    from amplifyp.gui.colours import GUIColours

    resolved = _resolve_font_family(font_family)

    if is_dimer:
        spans = []
        if top_name_line:
            spans.append(
                ft.TextSpan(
                    f"{top_name_line}\n",
                    style=ft.TextStyle(
                        color=GUIColours.PURPLE,
                        weight=ft.FontWeight.BOLD,
                    ),
                )
            )
        spans.extend(
            [
                ft.TextSpan(
                    f"{top_line}\n",
                    style=ft.TextStyle(
                        color=GUIColours.TEXT_ON_SURFACE,
                        weight=ft.FontWeight.BOLD,
                    ),
                ),
                ft.TextSpan(
                    f"{mid_line}\n",
                    style=ft.TextStyle(
                        color=GUIColours.FWD_PRIMER,
                        weight=ft.FontWeight.BOLD,
                    ),
                ),
                ft.TextSpan(
                    f"{bottom_line}" + ("\n" if bottom_name_line else ""),
                    style=ft.TextStyle(
                        color=GUIColours.TEXT_ON_SURFACE,
                        weight=ft.FontWeight.BOLD,
                    ),
                ),
            ]
        )
        if bottom_name_line:
            spans.append(
                ft.TextSpan(
                    bottom_name_line,
                    style=ft.TextStyle(
                        color=GUIColours.PURPLE,
                        weight=ft.FontWeight.BOLD,
                    ),
                )
            )
    else:
        # Context Map:
        # top_line = coordinates / arrows (black)
        # mid_line = primer row (black)
        # bottom_line = bonds (blue) + template (black)
        bottom_parts = bottom_line.split("\n")
        if len(bottom_parts) == 3:
            bonds_line = bottom_parts[0]
            comp_line = bottom_parts[1]
            template_line = bottom_parts[2]
        elif len(bottom_parts) >= 2:
            bonds_line = bottom_parts[0]
            comp_line = ""
            template_line = "\n".join(bottom_parts[1:])
        else:
            bonds_line = bottom_line
            comp_line = ""
            template_line = ""

        comp_span_text = f"{comp_line}\n" if comp_line else ""

        spans = [
            ft.TextSpan(
                f"{top_line}\n",
                style=ft.TextStyle(
                    color=GUIColours.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                f"{mid_line}\n",
                style=ft.TextStyle(
                    color=GUIColours.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                f"{bonds_line}\n",
                style=ft.TextStyle(
                    color=GUIColours.FWD_PRIMER,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                comp_span_text,
                style=ft.TextStyle(
                    color=GUIColours.MUTED_GREY,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
            ft.TextSpan(
                template_line,
                style=ft.TextStyle(
                    color=GUIColours.TEXT_ON_SURFACE,
                    weight=ft.FontWeight.BOLD,
                ),
            ),
        ]

    return ft.Text(
        spans=spans,
        font_family=resolved,
        size=font_size,
        selectable=True,
    )
