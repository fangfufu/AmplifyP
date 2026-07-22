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
from typing import cast

import flet as ft

from amplifyp.dimer import PrimerDimer
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.views.dimer.dimer_card import DimerCard
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

        dimer_card = DimerCard(
            d=dimer,
            settings=settings,
            font_family=font_family,
            show_names=False,
        )

        primer_seq = dimer.primer_1.seq
        title = f"Self-dimer ({len(primer_seq)} bp)"

        font_size_subheader = settings.get("font_size_subheader", 16)
        font_size_small = settings.get("font_size_small", 12)
        font_size_default = settings.get("font_size_default", 14)

        header_row = dimer_card._build_card_header(
            font_size_subheader, font_size_small
        )
        metric_controls = cast(ft.Row, header_row.controls[1])

        super().__init__(
            card_id=card_id,
            title=title,
            settings=settings,
            dismiss_callback=dismiss_callback,
            body_controls=[
                dimer_card._build_alignment_diagram(font_size_default)
            ],
            title_controls=list(metric_controls.controls),
        )
