# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""DimerTile component for Flet settings view."""

from typing import Any

import flet as ft

from amplifyp.gui.views.settings.score_table import ScoreTable


class DimerTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Dimer settings."""

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
        """Initialize the DimerTile."""
        from amplifyp.dna import Nucleotides

        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        self.set_pd_min_overlap = ft.TextField(
            label="Min Overlap",
            value="3",
            on_change=self.on_change_handler,
        )
        self.set_pd_threshold = ft.TextField(
            label="Threshold",
            value="60.0",
            on_change=self.on_change_handler,
        )

        self.settings_map["pd_min_overlap"] = self.set_pd_min_overlap
        self.settings_map["pd_threshold"] = self.set_pd_threshold

        # Map initialization (must happen before building ScoreTable)
        from amplifyp.gui.util import initialize_score_fields

        initialize_score_fields(
            settings_map=self.settings_map,
            prefix="pd_score",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=list(Nucleotides.PRIMER),
            on_change_handler=self.on_change_handler,
            font_size=font_size_default,
        )

        self.pd_table = ScoreTable(
            label="Primer Dimer Weights",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=list(Nucleotides.PRIMER),
            row_label="Primer",
            col_label="Primer",
            prefix="pd_score",
            settings_map=self.settings_map,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
        )

        super().__init__(
            title=ft.Text(
                "Primer Dimer Settings",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Row(
                        [
                            self.pd_table,
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
                                        self.set_pd_min_overlap,
                                        self.set_pd_threshold,
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
