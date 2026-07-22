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

"""DismissibleSelfDimerCard component for 1D primer designer view."""

from collections.abc import Callable

import flet as ft

from amplifyp.dimer import PrimerDimer
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.data_helpers import create_overlapped_sequence_view
from amplifyp.gui.views.pcr.dismissible_detail_card import DismissibleDetailCard


class DismissibleSelfDimerCard(DismissibleDetailCard):
    """Dismissible card displaying 1D self-dimer details."""

    def __init__(
        self,
        card_id: str,
        dimer: PrimerDimer,
        settings: GUISettings,
        dismiss_callback: Callable[[ft.Card], None],
        font_family: str = "Roboto Mono",
        step_index: int | None = None,
    ) -> None:
        """Initialise the DismissibleSelfDimerCard.

        Args:
            card_id: Unique identifier for this card.
            dimer: PrimerDimer object containing self-dimer alignment.
            settings: GUI settings object.
            dismiss_callback: Callback invoked when dismissed.
            font_family: Sequence alignment font family.
            step_index: Optional 0-indexed step number in truncation sequence.
        """
        self.dimer = dimer
        self.settings = settings
        self.font_family = font_family
        self.step_index = step_index

        font_size_small = settings.get("font_size_small", 12)
        font_size_default = settings.get("font_size_default", 14)

        primer_seq = dimer.primer_1.seq
        step_label = (
            f"Step {step_index + 1}: " if step_index is not None else ""
        )
        title = f"{step_label}Self-dimer ({len(primer_seq)} bp)"

        title_controls = [
            ft.Container(
                content=ft.Text(
                    f"Overlap: {dimer.overlap} bp",
                    weight=ft.FontWeight.BOLD,
                    color=GUIColours.DIAGRAM_BLACK,
                    size=font_size_small,
                ),
                bgcolor=GUIColours.SELECTED_ROW_BG,
                padding=ft.Padding(8, 4, 8, 4),
                border_radius=4,
            ),
            ft.Container(
                content=ft.Text(
                    f"Quality: {dimer.quality:.1f}",
                    weight=ft.FontWeight.BOLD,
                    color=GUIColours.DIAGRAM_BLACK,
                    size=font_size_small,
                ),
                bgcolor=GUIColours.SELECTED_ROW_BG,
                padding=ft.Padding(8, 4, 8, 4),
                border_radius=4,
            ),
        ]

        diagram = self._build_alignment_diagram(font_size_default)

        super().__init__(
            card_id=card_id,
            title=title,
            settings=settings,
            dismiss_callback=dismiss_callback,
            body_controls=[diagram],
            title_controls=title_controls,
        )

    def _build_alignment_diagram(self, font_size_default: int) -> ft.Container:
        """Build the visual sequence alignment diagram container control."""
        seq1 = self.dimer.primer_1.seq
        seq2 = self.dimer.primer_2.seq

        middle_str = self.dimer.binding_strength_str

        p2_line = f"5'-{seq2}-3'"
        mid_line = f"{' ' * (3 + self.dimer.p1_pos)}{middle_str}"
        p1_line = f"{' ' * self.dimer.p1_pos}3'-{seq1[::-1]}-5'"

        diagram_stack = create_overlapped_sequence_view(
            p2_line,
            mid_line,
            p1_line,
            font_family=self.font_family,
            font_size=font_size_default,
            is_dimer=True,
        )
        return ft.Container(
            content=ft.Row(
                [diagram_stack],
                scroll=ft.ScrollMode.ALWAYS,
            ),
            padding=ft.Padding(12, 4, 12, 4),
            border_radius=6,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            height=85,
        )
