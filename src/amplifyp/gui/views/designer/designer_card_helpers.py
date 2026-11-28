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

"""Shared card helper components and badge builders for designer views."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.utils.gui_helpers import copy_text_to_clipboard

if TYPE_CHECKING:
    from amplifyp.dna import Primer
    from amplifyp.gui.settings import GUISettings


def create_badge(
    text: str,
    bg_colour: str | None = None,
    text_colour: str = GUIColours.DIAGRAM_BLACK,
    font_size: int = 12,
    padding: ft.Padding | None = None,
) -> ft.Container:
    """Create a styled property badge container.

    Args:
        text: Badge label text.
        bg_colour: Background colour. Defaults to SELECTED_ROW_BG if None.
        text_colour: Badge text colour. Defaults to DIAGRAM_BLACK.
        font_size: Font size in pixels. Defaults to 12.
        padding: Custom padding. Defaults to (8, 4, 8, 4).

    Returns:
        A styled Flet Container with rounded corners.
    """
    return ft.Container(
        content=ft.Text(
            text,
            weight=ft.FontWeight.BOLD,
            color=text_colour,
            size=font_size,
        ),
        bgcolor=bg_colour or GUIColours.SELECTED_ROW_BG,
        padding=padding or ft.Padding(8, 4, 8, 4),
        border_radius=4,
    )


def format_primer_properties(
    primer: Primer, settings: GUISettings
) -> tuple[str, str]:
    """Calculate and format melting temperature and % AT for a primer.

    Args:
        primer: Primer object.
        settings: Application GUI settings for Tm calculation.

    Returns:
        Tuple of (tm_text, pct_at_text).
    """
    pct_at = primer.ratio_at() * 100.0
    try:
        tm_val, tm_range = settings.resolve_primer_tm_display(primer)
        if tm_range is not None:
            prefix = "~" if tm_range.is_approximate else ""
            tm_text = (
                f"Tm: {tm_val:.1f}°C ({prefix}{tm_range.min_tm:.1f} to "
                f"{tm_range.max_tm:.1f}°C)"
            )
        else:
            tm_text = f"Tm: {tm_val:.1f}°C"
    except (KeyError, ValueError, RuntimeError):
        tm_text = "Tm: N/A"

    pct_at_text = f"% AT: {pct_at:.1f}%"
    return tm_text, pct_at_text


def build_sequence_action_row(
    seq_control: ft.Control,
    seq: str,
    label: str,
    on_send_to_input: Callable[[], None] | None = None,
    vertical_alignment: ft.CrossAxisAlignment | None = None,
) -> tuple[ft.Row, ft.IconButton, ft.IconButton | None]:
    """Build a row with a sequence control, copy button, and send button.

    Args:
        seq_control: Control displaying the primer sequence.
        seq: Primer sequence string copied to clipboard when clicked.
        label: Primer label used in button tooltips.
        on_send_to_input: Optional callback to send this primer to the Input
            primer list.
        vertical_alignment: Optional vertical alignment for the row.

    Returns:
        Tuple of ``(row, copy_button, send_button)`` where ``send_button`` is
        None when ``on_send_to_input`` is omitted.
    """

    def _copy_seq(e: ft.ControlEvent | None) -> None:
        page = getattr(e, "page", None) if e is not None else None
        copy_text_to_clipboard(page, seq)

    copy_btn = ft.IconButton(
        icon=ft.Icons.COPY,
        icon_size=16,
        tooltip=f"Copy {label} sequence",
        on_click=_copy_seq,
    )

    action_row_controls: list[ft.Control] = [seq_control, copy_btn]
    send_btn: ft.IconButton | None = None
    if on_send_to_input is not None:

        def _send_to_input(_e: ft.ControlEvent | None) -> None:
            on_send_to_input()

        send_btn = ft.IconButton(
            icon=ft.Icons.PLAYLIST_ADD,
            icon_size=16,
            tooltip=f"Send {label} to Input primer list",
            on_click=_send_to_input,
        )
        action_row_controls.append(send_btn)

    row = (
        ft.Row(
            action_row_controls,
            spacing=4,
            vertical_alignment=vertical_alignment,
        )
        if vertical_alignment is not None
        else ft.Row(action_row_controls, spacing=4)
    )
    return row, copy_btn, send_btn


def build_primer_summary_row(
    label: str,
    primer: Primer,
    settings: GUISettings,
    font_size_small: int = 12,
    on_send_to_input: Callable[[], None] | None = None,
) -> ft.Column:
    """Build a primer sequence display row with copy button and property badges.

    Args:
        label: Primer label description (e.g. 'Forward Primer').
        primer: Primer instance.
        settings: Application GUI settings.
        font_size_small: Font size for small text and badges.
        on_send_to_input: Optional callback to send this primer to the Input
            primer list.

    Returns:
        Flet Column containing property badges, read-only sequence field, and
        copy button (plus optional send-to-input button).
    """
    tm_text, pct_at_text = format_primer_properties(primer, settings)
    font_family = str(settings.get("font_family", "Roboto Mono"))

    seq_field = ft.TextField(
        value=primer.seq,
        read_only=True,
        expand=True,
        dense=True,
        text_style=ft.TextStyle(font_family=font_family),
        border_color=GUIColours.OUTLINE,
        content_padding=ft.Padding(8, 4, 8, 4),
    )

    action_row, _copy_btn, _send_btn = build_sequence_action_row(
        seq_control=seq_field,
        seq=primer.seq,
        label=label,
        on_send_to_input=on_send_to_input,
    )

    return ft.Column(
        [
            ft.Row(
                [
                    ft.Text(
                        f"{label} ({len(primer.seq)} nt):",
                        weight=ft.FontWeight.BOLD,
                        size=font_size_small,
                    ),
                    ft.Row(
                        [
                            create_badge(
                                tm_text,
                                font_size=font_size_small,
                                padding=ft.Padding(6, 3, 6, 3),
                            ),
                            create_badge(
                                pct_at_text,
                                font_size=font_size_small,
                                padding=ft.Padding(6, 3, 6, 3),
                            ),
                        ],
                        spacing=8,
                    ),
                ],
                alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            ),
            action_row,
        ],
        spacing=4,
    )
