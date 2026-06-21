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

"""Dimer View for the Flet application."""

import logging
import traceback
from typing import TYPE_CHECKING, Any

import flet as ft

from amplifyp.dimer import PrimerDimerGenerator
from amplifyp.dna import Primer
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import MAX_DIMERS_RENDER, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import clean_sequence, show_error_dialog
from amplifyp.gui.views.dimer.dimer_card import DimerCard

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from amplifyp.dimer import PrimerDimer


class DimerView(ft.Column):  # type: ignore[misc]
    """View to show calculated dimer results."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
    ) -> None:
        """Initialise the DimerView."""
        super().__init__(expand=True)
        self.app_page = page
        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()
        self._cached_dimers: list[PrimerDimer] | None = None
        self._cached_state_key: tuple[dict[str, Any], dict[str, Any]] | None = (
            None
        )

        self.result_list = ft.ListView(
            expand=True, spacing=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.controls = [
            ft.Container(content=self.result_list, expand=True),
        ]

    def run_analysis(self) -> bool:
        """Run dimer analysis and populate the UI."""
        self.result_list.controls.clear()
        success = True
        try:
            current_state_key = (
                self.input_data.to_dict(),
                self.settings.to_dict(),
            )
            if (
                self._cached_state_key == current_state_key
                and self._cached_dimers is not None
            ):
                dimers = self._cached_dimers
            else:
                pd_settings = self.settings.get_primer_dimer_settings()
                generator = PrimerDimerGenerator(settings=pd_settings)
                primers = self.input_data.get_active_primers()
                for p in primers:
                    name = p["name"]
                    seq = clean_sequence(p["seq"])
                    generator.add_primer(Primer(sequence=seq, name=name))

                generator.analyse_primers()
                dimers = generator.primer_dimers
                self._cached_dimers = dimers
                self._cached_state_key = current_state_key

            if not dimers:
                self.result_list.controls.append(
                    ft.Container(
                        content=ft.Text(
                            "No primer dimers detected above threshold.",
                            size=self.settings.get("font_size_subheader", 16),
                            italic=True,
                            text_align=ft.TextAlign.CENTER,
                        ),
                        alignment=ft.Alignment(0, 0),
                        padding=20,
                    )
                )
            else:
                num_dimers = len(dimers)
                display_dimers = dimers
                if num_dimers > MAX_DIMERS_RENDER:
                    display_dimers = dimers[:MAX_DIMERS_RENDER]
                    self.result_list.controls.append(
                        ft.Container(
                            content=ft.Text(
                                f"Warning: {num_dimers} primer dimers "
                                "detected. Only the top "
                                f"{MAX_DIMERS_RENDER} strongest binding "
                                "dimers are displayed to prevent "
                                "UI freeze.",
                                color=GUIColours.ERROR_RED,
                                weight=ft.FontWeight.BOLD,
                            ),
                            padding=10,
                        )
                    )
                font_family = self.settings.get("font_family", "Roboto Mono")
                for d in display_dimers:
                    card = DimerCard(
                        d,
                        self.settings,
                        font_family=font_family,
                    )
                    self.result_list.controls.append(card)
        except (OSError, ValueError, RuntimeError) as ex:
            logger.exception("Dimer analysis failed: %s", ex)
            self.result_list.controls.append(
                ft.Text(
                    f"Error running analysis: {ex}\n{traceback.format_exc()}",
                    color=GUIColours.ERROR_RED,
                )
            )
            show_error_dialog(self.app_page, "Error running analysis", str(ex))
            success = False
        self.app_page.update()
        return success
