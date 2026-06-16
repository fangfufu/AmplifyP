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

from amplifyp.gui.views.settings.base_score_tile import BaseScoreTile


class ReplicationTile(BaseScoreTile):
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

        self.set_primability_cutoff = ft.TextField(
            label="Primability Cutoff",
            value="0.8",
            on_change=on_change_handler,
        )
        self.set_stability_cutoff = ft.TextField(
            label="Stability Cutoff",
            value="0.4",
            on_change=on_change_handler,
        )
        self.set_amp4_compat = ft.Checkbox(
            label="Amplify4 Compatibility Mode",
            value=False,
            on_change=on_change_handler,
        )

        settings_map["primability_cutoff"] = self.set_primability_cutoff
        settings_map["stability_cutoff"] = self.set_stability_cutoff
        settings_map["amp4_compat"] = self.set_amp4_compat

        col_headers = [c for c in Nucleotides.TEMPLATE if c != Nucleotides.GAP]

        super().__init__(
            settings=settings,
            settings_map=settings_map,
            on_change_handler=on_change_handler,
            header_size=header_size,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
            title="Origin of Replication Settings",
            score_table_label="Base Pair Weights",
            score_table_prefix="bp_score",
            row_headers=list(Nucleotides.PRIMER),
            col_headers=col_headers,
            row_label="Primer",
            col_label="Template",
            parameter_controls=[
                self.set_primability_cutoff,
                self.set_stability_cutoff,
                self.set_amp4_compat,
            ],
        )
