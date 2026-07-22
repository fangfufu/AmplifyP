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

from collections.abc import Callable

import flet as ft

from amplifyp.dimer import PrimerDimer
from amplifyp.dna import Primer
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.data_helpers import create_overlapped_sequence_view
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
    ) -> None:
        """Initialise the Dismissible2DCard.

        Args:
            card_id: Unique identifier for this card.
            step: PrimerDimers2D object containing the 4 dimer alignments.
            settings: GUI settings object.
            dismiss_callback: Callback invoked when dismissed.
            font_family: Sequence alignment font family.
        """
        self.step = step
        self.settings = settings
        self.font_family = font_family

        fwd_p = step.fwd_fwd.primer_1
        rev_p = step.rev_rev.primer_1
        fwd_len = len(fwd_p.seq)
        rev_len = len(rev_p.seq)

        title = f"2D Primer Pair (Forward: {fwd_len} bp, Reverse: {rev_len} bp)"

        font_size_small = settings.get("font_size_small", 12)
        font_size_default = settings.get("font_size_default", 14)

        # Mean overlap across all 4 dimers
        mean_overlap = step.mean_overlap

        def _make_badge(
            text: str, bg_colour: str | None = None
        ) -> ft.Container:
            return ft.Container(
                content=ft.Text(
                    text,
                    weight=ft.FontWeight.BOLD,
                    color=GUIColours.DIAGRAM_BLACK,
                    size=font_size_small,
                ),
                bgcolor=bg_colour or GUIColours.SELECTED_ROW_BG,
                padding=ft.Padding(8, 4, 8, 4),
                border_radius=4,
            )

        # Right hand side title metric boxes
        title_controls = [
            _make_badge(f"Max Quality: {step.max_quality:.1f}"),
            _make_badge(f"Mean Quality: {step.mean_quality:.1f}"),
            _make_badge(f"Max Overlap: {step.max_overlap} bp"),
            _make_badge(f"Mean Overlap: {mean_overlap:.1f} bp"),
        ]

        # Primer details section (Forward & Reverse)
        primer_info_section = self._build_primer_details_section(
            fwd_p=fwd_p,
            rev_p=rev_p,
            font_size_small=font_size_small,
        )

        # 4 Dimer pair subcontainers
        subcontainers_section = self._build_4_dimer_subcontainers(
            font_size_default=font_size_default,
            font_size_small=font_size_small,
        )

        body_controls = [
            primer_info_section,
            ft.Divider(height=1, color=GUIColours.OUTLINE_VARIANT),
            subcontainers_section,
        ]

        super().__init__(
            card_id=card_id,
            title=title,
            settings=settings,
            dismiss_callback=dismiss_callback,
            body_controls=body_controls,
            title_controls=title_controls,
        )

    def _build_primer_details_section(
        self, fwd_p: Primer, rev_p: Primer, font_size_small: int
    ) -> ft.Container:
        """Build section displaying Forward and Reverse primer details.

        Args:
            fwd_p: Forward Primer object.
            rev_p: Reverse Primer object.
            font_size_small: Font size for badges.

        Returns:
            Container control with read-only sequence fields and metric badges.
        """

        def _make_badge(text: str) -> ft.Container:
            return ft.Container(
                content=ft.Text(
                    text,
                    weight=ft.FontWeight.BOLD,
                    color=GUIColours.DIAGRAM_BLACK,
                    size=font_size_small,
                ),
                bgcolor=GUIColours.SELECTED_ROW_BG,
                padding=ft.Padding(6, 3, 6, 3),
                border_radius=4,
            )

        def _build_primer_row(label: str, primer: Primer) -> ft.Column:
            at_pairs = primer.count_at()
            gc_pairs = primer.count_cg()
            pct_at = primer.ratio_at() * 100.0

            try:
                tm_val = self.settings.calculate_primer_tm(primer)
                tm_text = f"Tm: {tm_val:.1f}°C"
            except (KeyError, ValueError, RuntimeError):
                tm_text = "Tm: N/A"

            seq_field = ft.TextField(
                value=primer.seq,
                read_only=True,
                expand=True,
                dense=True,
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
                                f"{label} ({len(primer.seq)} bp):",
                                weight=ft.FontWeight.BOLD,
                                size=font_size_small,
                            ),
                            _make_badge(tm_text),
                            _make_badge(f"AT Pairs: {at_pairs}"),
                            _make_badge(f"GC Pairs: {gc_pairs}"),
                            _make_badge(f"% AT: {pct_at:.1f}%"),
                        ],
                        alignment=ft.MainAxisAlignment.START,
                        spacing=8,
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

        fwd_col = _build_primer_row("Forward Primer", fwd_p)
        rev_col = _build_primer_row("Reverse Primer", rev_p)

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

    def _build_4_dimer_subcontainers(
        self, font_size_default: int, font_size_small: int
    ) -> ft.Column:
        """Build the 4 subcontainers for each primer dimer pair alignment.

        Args:
            font_size_default: Default sequence font size.
            font_size_small: Small text font size.

        Returns:
            Column containing 4 styled dimer pair alignment subcontainers.
        """
        dimer_pairs: list[tuple[str, PrimerDimer]] = [
            ("Forward Self-Dimer (Fwd-Fwd)", self.step.fwd_fwd),
            ("Reverse Self-Dimer (Rev-Rev)", self.step.rev_rev),
            ("Forward-Reverse Cross-Dimer (Fwd-Rev)", self.step.fwd_rev),
            ("Reverse-Forward Cross-Dimer (Rev-Fwd)", self.step.rev_fwd),
        ]

        subcontainers: list[ft.Control] = []

        for label, dimer in dimer_pairs:
            seq1 = dimer.primer_1.seq
            seq2 = dimer.primer_2.seq
            middle_str = dimer.binding_strength_str

            p2_line = f"5'-{seq2}-3'"
            mid_line = f"{' ' * (3 + dimer.p1_pos)}{middle_str}"
            p1_line = f"{' ' * dimer.p1_pos}3'-{seq1[::-1]}-5'"

            diagram_stack = create_overlapped_sequence_view(
                p2_line,
                mid_line,
                p1_line,
                font_family=self.font_family,
                font_size=font_size_default,
            )

            diagram_container = ft.Container(
                content=diagram_stack,
                padding=8,
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
                            ft.Container(
                                content=ft.Text(
                                    f"Quality: {dimer.quality:.1f}",
                                    weight=ft.FontWeight.BOLD,
                                    size=font_size_small,
                                ),
                                bgcolor=GUIColours.SELECTED_ROW_BG,
                                padding=ft.Padding(6, 2, 6, 2),
                                border_radius=4,
                            ),
                            ft.Container(
                                content=ft.Text(
                                    f"Overlap: {dimer.overlap} bp",
                                    weight=ft.FontWeight.BOLD,
                                    size=font_size_small,
                                ),
                                bgcolor=GUIColours.SELECTED_ROW_BG,
                                padding=ft.Padding(6, 2, 6, 2),
                                border_radius=4,
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
