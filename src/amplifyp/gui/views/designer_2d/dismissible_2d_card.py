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

"""Dismissible2DCard component for 2D primer designer view."""

from __future__ import annotations

from collections.abc import Callable
from typing import cast

import flet as ft

from amplifyp.dimer import PrimerDimer
from amplifyp.dna import Primer
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.data_helpers import create_overlapped_sequence_view
from amplifyp.gui.views.designer.designer_card_helpers import (
    build_primer_summary_row,
    create_badge,
)
from amplifyp.gui.views.pcr.dismissible_detail_card import DismissibleDetailCard
from amplifyp.primer_designer_2d import PrimerDimers2D


class Dismissible2DCard(DismissibleDetailCard):
    """Dismissible card displaying 2D primer pair details."""

    def __init__(
        self,
        card_id: str,
        step: PrimerDimers2D,
        settings: GUISettings,
        dismiss_callback: Callable[[ft.Card], None],
        font_family: str = "Roboto Mono",
        on_run_pcr_callback: (
            Callable[[str, str, str, str], None] | None
        ) = None,
        on_add_primer_callback: Callable[[str, str], None] | None = None,
        amplicon_count: int | None = None,
    ) -> None:
        """Initialise the Dismissible2DCard.

        Args:
            card_id: Unique identifier for this card.
            step: PrimerDimers2D object containing the 4 dimer alignments.
            settings: GUI settings object.
            dismiss_callback: Callback invoked when dismissed.
            font_family: Sequence alignment font family.
            on_run_pcr_callback: Optional callback to execute PCR with this
                primer pair.
            on_add_primer_callback: Optional callback to send a primer from
                this card to the Input primer list.
            amplicon_count: Optional number of amplicons predicted on template.
        """
        self.step = step
        self.settings = settings
        self.font_family = font_family
        self.on_run_pcr_callback = on_run_pcr_callback
        self.on_add_primer_callback = on_add_primer_callback
        self.amplicon_count = (
            amplicon_count
            if amplicon_count is not None
            else step.amplicon_count
        )

        fwd_p = step.fwd_fwd.primer_1
        rev_p = step.rev_rev.primer_1
        fwd_len = len(fwd_p.seq)
        rev_len = len(rev_p.seq)
        self._fwd_p = fwd_p
        self._rev_p = rev_p

        title = f"2D Primer Pair (Forward: {fwd_len} nt, Reverse: {rev_len} nt)"

        font_size_small = settings.get("font_size_small", 12)
        font_size_default = settings.get("font_size_default", 14)

        def _on_pcr_click(e: ft.ControlEvent) -> None:
            if self.on_run_pcr_callback:
                self.on_run_pcr_callback(
                    fwd_p.seq,
                    f"2D Fwd ({fwd_len} nt)",
                    rev_p.seq,
                    f"2D Rev ({rev_len} nt)",
                )

        self.pcr_button = ft.FilledTonalButton(
            "Run PCR",
            icon=ft.Icons.PLAY_ARROW,
            height=28,
            tooltip=(
                f"Run PCR using template with {fwd_len} nt and "
                f"{rev_len} nt primers"
            ),
            on_click=_on_pcr_click,
        )
        title_controls = self._build_title_controls(font_size_small)

        # Primer details section (Forward & Reverse)
        self.primer_info_section = self._build_primer_details_section(
            fwd_p=fwd_p,
            rev_p=rev_p,
            font_size_small=font_size_small,
        )
        self._primer_divider = ft.Divider(
            height=1, color=GUIColours.OUTLINE_VARIANT
        )

        # Dimer pair subcontainers
        self.dimer_subcontainers = self._build_dimer_subcontainers(
            font_size_default=font_size_default,
            font_size_small=font_size_small,
        )

        card_body = [
            self.primer_info_section,
            self._primer_divider,
            self.dimer_subcontainers,
        ]

        super().__init__(
            card_id=card_id,
            title=title,
            settings=settings,
            dismiss_callback=dismiss_callback,
            body_controls=card_body,
            title_controls=title_controls,
        )

    def _build_title_controls(self, font_size_small: int) -> list[ft.Control]:
        """Build header metric badges and Run PCR button."""
        title_controls: list[ft.Control] = [
            create_badge(
                f"Max Quality: {round(self.step.max_quality)}",
                font_size=font_size_small,
            ),
            create_badge(
                f"Max Overlap: {self.step.max_overlap} bp",
                font_size=font_size_small,
            ),
        ]
        if self.amplicon_count is not None:
            title_controls.append(
                create_badge(
                    f"Amplicons: {self.amplicon_count}",
                    font_size=font_size_small,
                )
            )
        title_controls.append(self.pcr_button)
        return title_controls

    def on_settings_change(self) -> None:
        """Update header badges, primer details, alignment cards, and font."""
        self.font_family = self.settings.get("font_family", self.font_family)
        font_size_small = self.settings.get("font_size_small", 12)
        font_size_default = self.settings.get("font_size_default", 14)

        container = cast(ft.Container, self.content)
        card_col = cast(ft.Column, container.content)
        title_row = cast(ft.Row, card_col.controls[0])
        title_row.controls = [
            title_row.controls[0],
            *self._build_title_controls(font_size_small),
            title_row.controls[-1],
        ]

        new_primer_info = self._build_primer_details_section(
            fwd_p=self._fwd_p,
            rev_p=self._rev_p,
            font_size_small=font_size_small,
        )
        self.primer_info_section.content = new_primer_info.content
        self.primer_info_section.bgcolor = new_primer_info.bgcolor
        self._primer_divider.color = GUIColours.OUTLINE_VARIANT

        new_subcontainers = self._build_dimer_subcontainers(
            font_size_default=font_size_default,
            font_size_small=font_size_small,
        )
        self.dimer_subcontainers.controls = new_subcontainers.controls
        try:
            if self.page:
                self.page.update()
        except RuntimeError:
            pass

    def update_ui(self) -> None:
        """Update card contents when settings change."""
        self.on_settings_change()

    def _build_primer_details_section(
        self, fwd_p: Primer, rev_p: Primer, font_size_small: int
    ) -> ft.Container:
        """Build the top section displaying sequence and properties for Fwd/Rev.

        Args:
            fwd_p: Forward primer object.
            rev_p: Reverse primer object.
            font_size_small: Font size for property badges.

        Returns:
            Container with formatted rows for forward and reverse primers.
        """
        fwd_name = f"2D Fwd ({len(fwd_p.seq)} nt)"
        rev_name = f"2D Rev ({len(rev_p.seq)} nt)"

        def _send_fwd() -> None:
            if self.on_add_primer_callback:
                self.on_add_primer_callback(fwd_name, fwd_p.seq)

        def _send_rev() -> None:
            if self.on_add_primer_callback:
                self.on_add_primer_callback(rev_name, rev_p.seq)

        fwd_col = build_primer_summary_row(
            "Forward Primer",
            fwd_p,
            self.settings,
            font_size_small,
            on_send_to_input=_send_fwd if self.on_add_primer_callback else None,
        )
        rev_col = build_primer_summary_row(
            "Reverse Primer",
            rev_p,
            self.settings,
            font_size_small,
            on_send_to_input=_send_rev if self.on_add_primer_callback else None,
        )

        return ft.Container(
            content=ft.Column(
                [
                    fwd_col,
                    rev_col,
                ],
                spacing=10,
            ),
            padding=8,
            bgcolor=GUIColours.SURFACE_VARIANT,
            border_radius=6,
        )

    def _build_dimer_subcontainers(
        self, font_size_default: int, font_size_small: int
    ) -> ft.Column:
        """Build the subcontainers for each primer dimer pair alignment.

        Args:
            font_size_default: Default sequence font size.
            font_size_small: Small text font size.

        Returns:
            Column containing styled dimer pair alignment subcontainers.
        """
        dimer_pairs: list[tuple[str, PrimerDimer]] = [
            ("Forward Self-Dimer (Fwd-Fwd)", self.step.fwd_fwd),
            ("Reverse Self-Dimer (Rev-Rev)", self.step.rev_rev),
            ("Forward-Reverse Cross-Dimer (Fwd-Rev)", self.step.fwd_rev),
        ]

        show_rev_fwd = self.settings.get("designer_2d_show_rev_fwd", False)
        if isinstance(show_rev_fwd, str):
            show_rev_fwd = show_rev_fwd.lower() in ("true", "1", "yes")
        if show_rev_fwd:
            dimer_pairs.append(
                ("Reverse-Forward Cross-Dimer (Rev-Fwd)", self.step.rev_fwd)
            )

        subcontainers: list[ft.Control] = []

        sep_setting = (
            self.settings.get("sequence_separator", "Dash ('-')")
            if hasattr(self.settings, "get")
            else "-"
        )

        sep = "-" if str(sep_setting) in ("Dash ('-')", "-") else " "

        for label, dimer in dimer_pairs:
            alignment = PrimerDimer.format_alignment(
                dimer, separator=sep, include_names=False
            )

            diagram_stack = create_overlapped_sequence_view(
                alignment.top_seq_line,
                alignment.bonds_line,
                alignment.bottom_seq_line,
                font_family=self.font_family,
                font_size=font_size_default,
                is_dimer=True,
            )

            diagram_container = ft.Container(
                content=ft.Row(
                    [diagram_stack],
                    scroll=ft.ScrollMode.ALWAYS,
                    expand=True,
                ),
                padding=8,
                alignment=ft.Alignment(-1, 0),
                bgcolor=GUIColours.DIAGRAM_BG,
                border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
                border_radius=4,
            )

            header_row = ft.Row(
                [
                    ft.Text(
                        label,
                        weight=ft.FontWeight.BOLD,
                        size=font_size_small + 1,
                        color=GUIColours.PRIMARY,
                    ),
                    ft.Row(
                        [
                            create_badge(
                                f"Quality: {round(dimer.quality)}",
                                font_size=font_size_small,
                                padding=ft.Padding(6, 2, 6, 2),
                            ),
                            create_badge(
                                f"Overlap: {dimer.overlap} bp",
                                font_size=font_size_small,
                                padding=ft.Padding(6, 2, 6, 2),
                            ),
                        ],
                        spacing=6,
                    ),
                ],
                alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            )

            subcontainer = ft.Container(
                content=ft.Column(
                    [
                        header_row,
                        diagram_container,
                    ],
                    spacing=6,
                ),
                padding=8,
                border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
                border_radius=6,
            )
            subcontainers.append(subcontainer)

        return ft.Column(subcontainers, spacing=8)
