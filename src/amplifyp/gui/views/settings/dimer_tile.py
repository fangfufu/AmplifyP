# Copyright (C) 2026 Fufu Fang
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

"""DimerTile component for Flet settings view."""

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.views.settings.base_score_tile import BaseScoreTile


class DimerTile(BaseScoreTile):
    """Expansion tile for Dimer settings."""

    def __init__(
        self,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.Event | None], None],
        header_size: int,
        font_size_default: int,
        font_size_micro: int,
        font_size_table_header: int,
    ) -> None:
        """Initialise the DimerTile.

        Args:
            settings: The settings object.
            settings_map: A dictionary mapping setting keys to UI
                components for population and retrieval.
            on_change_handler: The handler to call when a setting changes.
            header_size: The size of the expansion tile header text.
            font_size_default: Default font size for text elements.
            font_size_micro: Micro font size for small labels.
            font_size_table_header: Font size for table header cells.
        """
        from amplifyp.dna import Nucleotides

        self.set_pd_min_overlap = ft.TextField(
            label="Min Overlap",
            value="3",
            on_change=on_change_handler,
            border_color=GUIColours.OUTLINE,
        )
        self.set_pd_threshold = ft.TextField(
            label="Threshold",
            value="60.0",
            on_change=on_change_handler,
            border_color=GUIColours.OUTLINE,
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
