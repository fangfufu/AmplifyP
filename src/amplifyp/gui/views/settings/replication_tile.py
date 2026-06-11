# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""ReplicationTile component for Flet settings view."""

from typing import Any

import flet as ft

from amplifyp.gui.views.settings.score_table import ScoreTable


class ReplicationTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Origin of Replication settings."""

    def __init__(
        self,
        settings: Any,
        settings_map: dict[str, Any],
        on_change_handler: Any,
        header_size: int,
        font_size_default: int,
        font_size_micro: int,
        font_size_table_header: int,
    ) -> None:
        """Initialize the ReplicationTile."""
        from amplifyp.dna import Nucleotides

        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        self.set_primability_cutoff = ft.TextField(
            label="Primability Cutoff",
            value="0.8",
            on_change=self.on_change_handler,
        )
        self.set_stability_cutoff = ft.TextField(
            label="Stability Cutoff",
            value="0.4",
            on_change=self.on_change_handler,
        )
        self.set_amp4_compat = ft.Checkbox(
            label="Amplify4 Compatibility Mode",
            value=False,
            on_change=self.on_change_handler,
        )

        self.settings_map["primability_cutoff"] = self.set_primability_cutoff
        self.settings_map["stability_cutoff"] = self.set_stability_cutoff
        self.settings_map["amp4_compat"] = self.set_amp4_compat

        col_headers = [c for c in Nucleotides.TEMPLATE if c != Nucleotides.GAP]

        # Map initialization (must happen before building ScoreTable)
        from amplifyp.gui.util import initialize_score_fields

        initialize_score_fields(
            settings_map=self.settings_map,
            prefix="bp_score",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=col_headers,
            on_change_handler=self.on_change_handler,
            font_size=font_size_default,
        )

        self.bp_table = ScoreTable(
            label="Base Pair Weights",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=col_headers,
            row_label="Primer",
            col_label="Template",
            prefix="bp_score",
            settings_map=self.settings_map,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
        )

        super().__init__(
            title=ft.Text(
                "Origin of Replication Settings",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Row(
                        [
                            self.bp_table,
                            ft.VerticalDivider(),
                            ft.Container(
                                content=ft.Column(
                                    [
                                        ft.Text(
                                            "Parameters",
                                            weight=ft.FontWeight.BOLD,
                                            size=self.settings.get(
                                                "font_size_default", 14
                                            ),
                                        ),
                                        self.set_primability_cutoff,
                                        self.set_stability_cutoff,
                                        self.set_amp4_compat,
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=220,
                            ),
                        ],
                        vertical_alignment=ft.CrossAxisAlignment.START,
                        alignment=ft.MainAxisAlignment.CENTER,
                        scroll=ft.ScrollMode.ALWAYS,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )
