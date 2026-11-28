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

"""Info panel component for DNA primers."""

from __future__ import annotations

import logging
from collections.abc import Callable

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.data_helpers import clean_sequence

logger = logging.getLogger(__name__)


class PrimerInfoPanel(ft.Card):  # type: ignore[misc]
    """Panel to display Tm, redundancy, and dimer info for a primer."""

    def __init__(self, settings: GUISettings, font_family: str) -> None:
        """Initialise the PrimerInfoPanel."""
        super().__init__()
        self.settings = settings
        self.font_family = font_family
        self._on_dismiss_callback: Callable[[], None] | None = None

        self.info_header_text = ft.Text(
            "Primer: -",
            weight=ft.FontWeight.BOLD,
            size=self.settings.get("font_size_subheader", 16),
            color=GUIColours.TEXT_ON_SURFACE,
            selectable=True,
        )
        self.info_header = ft.Container(
            content=self.info_header_text,
            alignment=ft.Alignment(-1, 0),
        )

        self.close_button = ft.IconButton(
            icon=ft.Icons.CLOSE,
            icon_size=16,
            tooltip="Dismiss",
            on_click=self._on_close_click,
            padding=0,
            height=20,
            width=20,
        )

        self.info_seq_text = ft.Text(
            "",
            font_family=font_family,
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_stats_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_dimer_text = ft.Text(
            "",
            color=GUIColours.ERROR_RED,
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_dimer_card_container = ft.Container(visible=False)

        self.info_body_column = ft.Column(
            [
                self.info_seq_text,
                self.info_stats_text,
                self.info_dimer_text,
                self.info_dimer_card_container,
            ],
            spacing=1,
        )

        self.content = ft.Container(
            padding=ft.Padding(10, 2, 10, 2),
            content=ft.Column(
                [
                    ft.Row(
                        [
                            self.info_header,
                            self.close_button,
                        ],
                        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    ),
                    self.info_body_column,
                ],
                spacing=2,
            ),
        )
        self.visible = False

    def _apply_height_mode(self, fixed_height_enabled: bool) -> None:
        """Configure height and body scrolling based on fixed-height mode."""
        if fixed_height_enabled:
            raw_height = (
                self.settings.get("primer_info_panel_height", 260) or 260
            )
            try:
                saved_height = float(raw_height)
            except (OverflowError, TypeError, ValueError):
                saved_height = 260.0
            page_height: float | None = None
            try:
                if self.page and self.page.height:
                    page_height = float(self.page.height)
            except RuntimeError:
                pass
            if page_height is not None:
                max_height = max(120.0, page_height - 200.0)
                self.height = max(80.0, min(max_height, saved_height))
            else:
                self.height = saved_height
            self.info_body_column.scroll = ft.ScrollMode.AUTO
            self.info_body_column.expand = True
        else:
            self.height = None
            self.info_body_column.scroll = None
            self.info_body_column.expand = False

    def _on_close_click(self, e: ft.Event) -> None:
        """Handle close button click."""
        self._show_empty_or_hide()
        self._update_safe()
        if self._on_dismiss_callback:
            self._on_dismiss_callback()

    def _show_empty_or_hide(self) -> None:
        """Show an empty placeholder in fixed-height mode, or hide the panel."""
        if self.settings.get("primer_info_panel_fixed_height", False):
            self.info_header_text.value = "Primer: -"
            self.info_seq_text.value = ""
            self.info_stats_text.value = ""
            self.info_dimer_text.value = ""
            self.info_dimer_text.visible = False
            self.info_dimer_card_container.content = None
            self.info_dimer_card_container.visible = False
            self.close_button.visible = False
            self._apply_height_mode(True)
            self.visible = True
        else:
            self.close_button.visible = True
            self.visible = False
            self._apply_height_mode(False)
            self.info_dimer_card_container.visible = False
            self.info_dimer_card_container.content = None

    def update_panel(
        self,
        focused_idx: int | None,
        input_data: GUIInput,
        app_page: ft.Page,
        on_update_highlights: Callable[[], None],
        on_dismiss: Callable[[], None] | None = None,
    ) -> None:
        """Update the primer information panel based on the focused primer."""
        self._on_dismiss_callback = on_dismiss
        self.font_family = str(
            self.settings.get("font_family", self.font_family)
        )
        self.info_header_text.color = GUIColours.TEXT_ON_SURFACE
        self.info_seq_text.font_family = self.font_family
        self.info_dimer_text.color = GUIColours.ERROR_RED
        if focused_idx is None:
            self._show_empty_or_hide()
            on_update_highlights()
            self._update_safe()
            return

        try:
            primers = input_data.primers
            if focused_idx < 0 or focused_idx >= len(primers):
                self._show_empty_or_hide()
                self._update_safe()
                return

            p_data = primers[focused_idx]
            name_val = p_data.get("name", "").strip()
            seq_val = clean_sequence(p_data.get("seq", ""))

            if not seq_val:
                self._show_empty_or_hide()
                self._update_safe()
                return

            from amplifyp.dimer import PrimerDimerGenerator
            from amplifyp.dna import Primer

            primer_obj = Primer(sequence=seq_val, name=name_val)

            header_text = (
                f"Primer: {name_val}" if name_val else f"Primer: {seq_val}"
            )
            self.info_header_text.value = header_text
            self.info_seq_text.value = (
                f"{len(primer_obj)} nt:   {primer_obj.seq}"
            )

            try:
                tm, tm_range = self.settings.resolve_primer_tm_display(
                    primer_obj
                )
                if tm_range is not None:
                    label = (
                        "approx. range" if tm_range.is_approximate else "range"
                    )
                    tm_str = (
                        f"Tm = {tm:.2f}°C (max; {label} "
                        f"{tm_range.min_tm:.2f} to {tm_range.max_tm:.2f}°C)"
                    )
                else:
                    tm_str = f"Tm = {tm:.2f}°C"
            except (KeyError, ValueError):
                tm_str = "Tm = N/A"

            at_str = f"{primer_obj.ratio_at() * 100:.1f}% AT"

            if primer_obj.redundant_base_count == 0:
                redundancy_str = "No redundant bases."
            else:
                redundancy_str = (
                    f"{primer_obj.redundant_base_count} redundant bases "
                    f"(redundancy fold = {primer_obj.redundancy_fold})."
                )

            self.info_stats_text.value = f"{tm_str}, {at_str}, {redundancy_str}"

            pd_settings = self.settings.get_primer_dimer_settings()
            generator = PrimerDimerGenerator(settings=pd_settings)

            max_dimer = None
            res_self = generator.generate_primer_dimer(primer_obj, primer_obj)
            if (
                res_self.overlap > pd_settings.min_overlap
                and res_self.quality > pd_settings.threshold
            ):
                max_dimer = res_self

            if max_dimer is not None:
                from amplifyp.gui.views.dimer.dimer_card import DimerCard

                self.info_dimer_text.value = ""
                self.info_dimer_text.visible = False
                self.info_dimer_card_container.content = DimerCard(
                    max_dimer, self.settings, self.font_family, show_names=False
                )
                self.info_dimer_card_container.visible = True
            else:
                self.info_dimer_text.value = ""
                self.info_dimer_text.visible = False
                self.info_dimer_card_container.content = None
                self.info_dimer_card_container.visible = False

            self._apply_height_mode(
                bool(self.settings.get("primer_info_panel_fixed_height", False))
            )
            self.close_button.visible = True
            self.visible = True

        except (ValueError, AttributeError, ArithmeticError):
            logger.debug("Failed to calculate primer info, hiding panel")
            self._show_empty_or_hide()

        self._update_safe()

    def _update_safe(self) -> None:
        """Safely update the panel if attached to a page."""
        try:
            if self.page:
                self.update()
        except RuntimeError:
            logger.debug("Info panel detached, skipping update")
