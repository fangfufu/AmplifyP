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

from typing import TYPE_CHECKING

import flet as ft

from amplifyp.gui.colours import GUIColours

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
        tm_val = settings.calculate_primer_tm(primer)
        tm_text = f"Tm: {tm_val:.1f}°C"
    except (KeyError, ValueError, RuntimeError):
        tm_text = "Tm: N/A"

    pct_at_text = f"% AT: {pct_at:.1f}%"
    return tm_text, pct_at_text


def build_primer_summary_row(
    label: str,
    primer: Primer,
    settings: GUISettings,
    font_size_small: int = 12,
) -> ft.Column:
    """Build a primer sequence display row with copy button and property badges.

    Args:
        label: Primer label description (e.g. 'Forward Primer').
        primer: Primer instance.
        settings: Application GUI settings.
        font_size_small: Font size for small text and badges.

    Returns:
        Flet Column containing property badges, read-only sequence field, and
        copy button.
    """
    tm_text, pct_at_text = format_primer_properties(primer, settings)

    seq_field = ft.TextField(
        value=primer.seq,
        read_only=True,
        expand=True,
        dense=True,
        border_color=GUIColours.OUTLINE,
        content_padding=ft.Padding(8, 4, 8, 4),
    )

    async def _copy_seq_async() -> None:
        await ft.Clipboard().set(primer.seq)

    def _copy_seq(e: ft.ControlEvent) -> None:
        try:
            if e.page:
                e.page.run_task(_copy_seq_async)
        except RuntimeError:
            pass

    copy_btn = ft.IconButton(
        icon=ft.Icons.COPY,
        icon_size=16,
        tooltip=f"Copy {label} sequence",
        on_click=_copy_seq,
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
            ft.Row(
                [
                    seq_field,
                    copy_btn,
                ],
                spacing=4,
            ),
        ],
        spacing=4,
    )
