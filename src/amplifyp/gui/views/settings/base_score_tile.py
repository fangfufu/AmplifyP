# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Base score settings tile for Flet settings view."""

from typing import Any

import flet as ft

from amplifyp.gui.views.settings.score_table import ScoreTable


class BaseScoreTile(ft.ExpansionTile):  # type: ignore[misc]
    """Base settings tile with a ScoreTable and parameter controls."""

    def __init__(
        self,
        settings: Any,
        settings_map: dict[str, Any],
        on_change_handler: Any,
        header_size: int,
        font_size_default: int,
        font_size_micro: int,
        font_size_table_header: int,
        title: str,
        score_table_label: str,
        score_table_prefix: str,
        row_headers: list[str],
        col_headers: list[str],
        row_label: str,
        col_label: str,
        parameter_controls: list[ft.Control],
    ) -> None:
        """Initialize the BaseScoreTile."""
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler

        # Map initialization (must happen before building ScoreTable)
        from amplifyp.gui.util import initialize_score_fields

        initialize_score_fields(
            settings_map=self.settings_map,
            prefix=score_table_prefix,
            row_headers=row_headers,
            col_headers=col_headers,
            on_change_handler=self.on_change_handler,
            font_size=font_size_default,
        )

        self.score_table = ScoreTable(
            label=score_table_label,
            row_headers=row_headers,
            col_headers=col_headers,
            row_label=row_label,
            col_label=col_label,
            prefix=score_table_prefix,
            settings_map=self.settings_map,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
        )

        super().__init__(
            title=ft.Text(
                title,
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Row(
                        [
                            self.score_table,
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
                                        *parameter_controls,
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
