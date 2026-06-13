# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""Flet component displaying details of a selected DNA primer."""

from typing import Any

import flet as ft

from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import clean_sequence


class PrimerInfoPanel(ft.Card):  # type: ignore[misc]
    """Panel to display Tm, redundancy, and dimer info for a primer."""

    def __init__(self, settings: GUISettings, font_family: str) -> None:
        """Initialize the PrimerInfoPanel."""
        super().__init__()
        self.settings = settings
        self._on_dismiss_callback: Any = None

        self.info_header_text = ft.Text(
            "Primer: -",
            weight=ft.FontWeight.BOLD,
            size=self.settings.get("font_size_subheader", 16),
            color=GUIColors.TEXT_ON_SURFACE,
            selectable=True,
        )
        self.info_header = ft.Container(
            content=self.info_header_text,
            alignment=ft.Alignment(-1, 0),
        )

        self.close_button = ft.IconButton(
            icon=ft.Icons.CLOSE,
            icon_size=18,
            tooltip="Dismiss",
            on_click=self._on_close_click,
        )

        self.info_seq_text = ft.Text(
            "",
            font_family=font_family,
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_tm_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_pairs_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_redundancy_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )
        self.info_dimer_text = ft.Text(
            "",
            color=GUIColors.ERROR_RED,
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )

        self.content = ft.Container(
            padding=10,
            content=ft.Column(
                [
                    ft.Row(
                        [
                            self.info_header,
                            self.close_button,
                        ],
                        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                    ),
                    ft.Column(
                        [
                            self.info_seq_text,
                            self.info_tm_text,
                            self.info_pairs_text,
                            self.info_redundancy_text,
                            self.info_dimer_text,
                        ],
                        spacing=3,
                    ),
                ],
                spacing=10,
            ),
        )
        self.visible = False

    def _on_close_click(self, e: Any) -> None:
        """Handle close button click: clear focus and hide panel."""
        self.visible = False
        if self._on_dismiss_callback:
            self._on_dismiss_callback()

    def update_panel(
        self,
        focused_idx: int | None,
        input_data: GUIInput,
        app_page: ft.Page,
        on_update_highlights: Any,
        on_dismiss: Any = None,
    ) -> None:
        """Update the primer information panel based on the focused primer."""
        self._on_dismiss_callback = on_dismiss
        if focused_idx is None:
            self.visible = False
            on_update_highlights()
            app_page.update()
            return

        try:
            primers = input_data.primers
            if focused_idx < 0 or focused_idx >= len(primers):
                self.visible = False
                app_page.update()
                return

            p_data = primers[focused_idx]
            name_val = p_data.get("name", "").strip()
            seq_val = clean_sequence(p_data.get("seq", ""))

            if not seq_val:
                self.visible = False
                app_page.update()
                return

            from amplifyp.dimer import PrimerDimerGenerator
            from amplifyp.dna import Primer

            primer_obj = Primer(sequence=seq_val, name=name_val)

            # Header
            header_text = (
                f"Primer: {name_val}" if name_val else f"Primer: {seq_val}"
            )
            self.info_header_text.value = header_text

            # Sequence details
            self.info_seq_text.value = (
                f"{len(primer_obj)} bp:   {primer_obj.seq}"
            )

            # Melting temperature
            tm = self.settings.calculate_primer_tm(primer_obj)
            self.info_tm_text.value = f"Tm = {tm:.2f}°C"

            # Base counts / percentage
            self.info_pairs_text.value = (
                f"{primer_obj.count_at()} AT Pairs, "
                f"{primer_obj.count_cg()} GC Pairs, "
                f"{primer_obj.ratio_at() * 100:.1f}% AT"
            )

            # Redundancy
            if primer_obj.redundant_base_count == 0:
                self.info_redundancy_text.value = "No redundant bases."
            else:
                self.info_redundancy_text.value = (
                    f"{primer_obj.redundant_base_count} redundant bases "
                    f"(redundancy fold = {primer_obj.redundancy_fold})."
                )

            # Dimer potential check
            pd_settings = self.settings.get_primer_dimer_settings()
            generator = PrimerDimerGenerator(settings=pd_settings)

            # Check self-dimer potential
            max_dimer = None

            # Always check self-dimer
            res_self = generator.generate_primer_dimer(primer_obj, primer_obj)
            if (
                res_self.overlap > pd_settings.min_overlap
                and res_self.quality > pd_settings.threshold
            ):
                max_dimer = res_self

            if max_dimer is not None:
                self.info_dimer_text.value = (
                    "Potential Primer Dimer with quality = "
                    f"{max_dimer.quality:.1f} "
                    f"and overlap = {max_dimer.overlap}"
                )
                self.info_dimer_text.visible = True
            else:
                self.info_dimer_text.value = ""
                self.info_dimer_text.visible = False

            self.visible = True

        except Exception:
            self.visible = False

        app_page.update()
