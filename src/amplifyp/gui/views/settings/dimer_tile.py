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

from amplifyp.gui.views.settings.base_score_tile import BaseScoreTile


class DimerTile(BaseScoreTile):
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

        self.set_pd_min_overlap = ft.TextField(
            label="Min Overlap",
            value="3",
            on_change=on_change_handler,
        )
        self.set_pd_threshold = ft.TextField(
            label="Threshold",
            value="60.0",
            on_change=on_change_handler,
        )

        settings_map["pd_min_overlap"] = self.set_pd_min_overlap
        settings_map["pd_threshold"] = self.set_pd_threshold

        super().__init__(
            settings=settings,
            settings_map=settings_map,
            on_change_handler=on_change_handler,
            header_size=header_size,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
            title="Primer Dimer Settings",
            score_table_label="Primer Dimer Weights",
            score_table_prefix="pd_score",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=list(Nucleotides.PRIMER),
            row_label="Primer",
            col_label="Primer",
            parameter_controls=[
                self.set_pd_min_overlap,
                self.set_pd_threshold,
            ],
        )
