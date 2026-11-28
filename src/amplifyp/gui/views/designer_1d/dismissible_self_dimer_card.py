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
from amplifyp.gui.views.designer.designer_card_helpers import (
    create_badge,
    format_primer_properties,
)
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
        on_run_pcr_callback: Callable[[str, str], None] | None = None,
        origin_count: int | None = None,
    ) -> None:
        """Initialise the DismissibleSelfDimerCard.

        Args:
            card_id: Unique identifier for this card.
            dimer: PrimerDimer object containing self-dimer alignment.
            settings: GUI settings object.
            dismiss_callback: Callback invoked when dismissed.
            font_family: Sequence alignment font family.
            step_index: Optional 0-indexed step number in truncation sequence.
            on_run_pcr_callback: Optional callback to execute PCR with this
                primer.
            origin_count: Optional number of template binding sites.
        """
        self.dimer = dimer
        self.settings = settings
        self.font_family = font_family
        self.step_index = step_index
        self.on_run_pcr_callback = on_run_pcr_callback
        self.origin_count = origin_count

        primer_seq = dimer.primer_1.seq
        title = f"Self-dimer ({len(primer_seq)} nt)"

        def _on_pcr_click(e: ft.ControlEvent) -> None:
            if self.on_run_pcr_callback:
                self.on_run_pcr_callback(
                    primer_seq, f"1D Primer ({len(primer_seq)} nt)"
                )

        self.pcr_button = ft.FilledTonalButton(
            "Run PCR",
            icon=ft.Icons.PLAY_ARROW,
            height=28,
            tooltip=f"Run PCR using template with {len(primer_seq)} nt primer",
            on_click=_on_pcr_click,
        )

        self._badges_container = ft.Container(
            content=self._build_badges_row(),
        )
        self._alignment_container = ft.Container(
            content=self._build_alignment_diagram(),
        )

        super().__init__(
            card_id=card_id,
            title=title,
            settings=settings,
            dismiss_callback=dismiss_callback,
            body_controls=[
                self._badges_container,
                self._alignment_container,
            ],
            title_controls=[self.pcr_button],
        )

    @staticmethod
    def _is_enabled(val: object) -> bool:
        """Coerce bool or string setting values to a boolean."""
        if isinstance(val, str):
            return val.lower() not in ("false", "0", "no", "")
        return bool(val)

    def _build_badge_controls(self) -> list[ft.Control]:
        """Build metric badges for the card from current settings."""
        dimer_card = DimerCard(
            d=self.dimer,
            settings=self.settings,
            font_family=self.font_family,
            show_names=False,
        )
        font_size_subheader = self.settings._safe_int("font_size_subheader", 16)
        font_size_small = self.settings._safe_int("font_size_small", 12)
        header_row = dimer_card._build_card_header(
            font_size_subheader, font_size_small
        )
        metric_controls = cast(ft.Row, header_row.controls[1])

        tm_text, pct_at_text = format_primer_properties(
            self.dimer.primer_1, self.settings
        )
        badge_controls: list[ft.Control] = list(metric_controls.controls)
        if self._is_enabled(self.settings.get("designer_1d_show_tm", True)):
            badge_controls.append(
                create_badge(tm_text, font_size=font_size_small)
            )
        if self._is_enabled(
            self.settings.get("designer_1d_show_pct_at", False)
        ):
            badge_controls.append(
                create_badge(pct_at_text, font_size=font_size_small)
            )
        if self.origin_count is not None:
            badge_controls.append(
                create_badge(
                    f"Sites: {self.origin_count}", font_size=font_size_small
                )
            )
        return badge_controls

    def _build_badges_row(self) -> ft.Row:
        """Build a wrapping row of property badges."""
        return ft.Row(
            self._build_badge_controls(),
            wrap=True,
            spacing=6,
            run_spacing=4,
            vertical_alignment=ft.CrossAxisAlignment.CENTER,
        )

    def _build_alignment_diagram(self) -> ft.Control:
        """Build self-dimer sequence alignment diagram from current settings."""
        dimer_card = DimerCard(
            d=self.dimer,
            settings=self.settings,
            font_family=self.font_family,
            show_names=False,
        )
        font_size_default = self.settings._safe_int("font_size_default", 14)
        return dimer_card._build_alignment_diagram(font_size_default)

    def on_settings_change(self) -> None:
        """Refresh header badges, alignment diagram, and font from settings."""
        self.font_family = str(
            self.settings.get("font_family", self.font_family)
        )
        self._badges_container.content = self._build_badges_row()
        self._alignment_container.content = self._build_alignment_diagram()
        try:
            if self.page:
                self.page.update()
        except RuntimeError:
            pass

    def update_ui(self) -> None:
        """Alias for on_settings_change."""
        self.on_settings_change()
