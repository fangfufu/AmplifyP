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

"""Settings View for the Flet application."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any

import flet as ft

from amplifyp.gui.settings import GUISettings

if TYPE_CHECKING:
    from amplifyp.gui.util import BorderedCheckbox
from amplifyp.gui.views.settings.appearance_tile import AppearanceTile
from amplifyp.gui.views.settings.dimer_tile import DimerTile
from amplifyp.gui.views.settings.replication_tile import ReplicationTile
from amplifyp.gui.views.settings.tm_tile import TmTile
from amplifyp.settings import ReplicationSettings


class SettingsView(ft.ListView):  # type: ignore[misc]
    """Settings view for configuring PCR replication settings."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings | None = None,
        on_change: Callable[[ft.ControlEvent], None] | None = None,
        on_apply: Callable[[ft.ControlEvent | None], None] | None = None,
        on_reset: Callable[[ft.ControlEvent | None], None] | None = None,
    ) -> None:
        """Initialise the SettingsView."""
        super().__init__(
            expand=True, spacing=20, padding=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.app_page = page
        self.settings = settings if settings is not None else GUISettings()
        self.on_change = on_change
        self.on_apply = on_apply
        self.on_reset = on_reset

        self.settings_map: dict[str, Any] = {}

        header_size = self.settings.get("font_size_header", 18)
        font_size_default = self.settings.get("font_size_default", 14)
        font_size_micro = self.settings.get("font_size_micro", 10)
        font_size_table_header = self.settings.get("font_size_table_header", 15)

        # Initialise sub-control tiles
        self.replication_tile = ReplicationTile(
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,
            header_size=header_size,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
        )

        self.tm_tile = TmTile(
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,
            header_size=header_size,
        )

        self.dimer_tile = DimerTile(
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,
            header_size=header_size,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
        )

        self.appearance_tile = AppearanceTile(
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,
            header_size=header_size,
        )

        self.controls = [
            self.replication_tile,
            self.tm_tile,
            self.dimer_tile,
            self.appearance_tile,
            ft.Divider(),
            self._build_action_buttons(),
        ]

        # Sync initial UI state
        self.update_ui()

    # Backwards compatibility properties for unit tests
    @property
    def set_primability_cutoff(self) -> ft.TextField:
        """Get the primability cutoff text field."""
        return self.replication_tile.set_primability_cutoff

    @property
    def set_stability_cutoff(self) -> ft.TextField:
        """Get the stability cutoff text field."""
        return self.replication_tile.set_stability_cutoff

    @property
    def set_amp4_compat(self) -> ft.Checkbox | BorderedCheckbox:
        """Get the amplify4 compatibility mode checkbox."""
        return self.replication_tile.set_amp4_compat

    @property
    def set_tm_dna_conc(self) -> ft.TextField:
        """Get the DNA conc text field."""
        return self.tm_tile.set_tm_dna_conc

    @property
    def set_tm_dnap_conc(self) -> ft.TextField:
        """Get the DNA Pol conc text field."""
        return self.tm_tile.set_tm_dnap_conc

    @property
    def set_tm_mono_salt(self) -> ft.TextField:
        """Get the monovalent salt conc text field."""
        return self.tm_tile.set_tm_mono_salt

    @property
    def set_tm_div_salt(self) -> ft.TextField:
        """Get the divalent salt conc text field."""
        return self.tm_tile.set_tm_div_salt

    @property
    def set_tm_dNTP_conc(self) -> ft.TextField:
        """Get the dNTP conc text field."""
        return self.tm_tile.set_tm_dNTP_conc

    @property
    def set_tm_method(self) -> ft.Dropdown:
        """Get the Tm method calculation dropdown."""
        return self.tm_tile.set_tm_method

    @property
    def set_pd_min_overlap(self) -> ft.TextField:
        """Get the dimer min overlap text field."""
        return self.dimer_tile.set_pd_min_overlap

    @property
    def set_pd_threshold(self) -> ft.TextField:
        """Get the dimer threshold text field."""
        return self.dimer_tile.set_pd_threshold

    @property
    def set_font_family(self) -> ft.Dropdown:
        """Get the font family dropdown."""
        return self.appearance_tile.set_font_family

    @property
    def set_colour_deficient(self) -> ft.Checkbox:
        """Get the colour deficient mode checkbox."""
        return self.appearance_tile.set_colour_deficient

    @property
    def set_improved_visualisation(self) -> ft.Checkbox | BorderedCheckbox:
        """Get the improved visualisation mode checkbox."""
        return self.appearance_tile.set_improved_visualisation

    @property
    def set_show_primer_temperature(self) -> ft.Checkbox | BorderedCheckbox:
        """Get the show primer temperature checkbox."""
        return self.tm_tile.set_show_primer_temperature

    @property
    def set_tm_colour_scheme(self) -> ft.Dropdown:
        """Get the Tm colour scheme dropdown."""
        return self.tm_tile.set_tm_colour_scheme

    def _build_action_buttons(self) -> ft.Row:
        """Build the Action buttons Row (Apply & Reset)."""
        return ft.Row(
            list[ft.Control](
                [
                    ft.FilledButton(
                        "Apply",
                        icon=ft.Icons.DONE,
                        on_click=self._on_apply_handler,
                    ),
                    ft.OutlinedButton(
                        "Reset to Default",
                        icon=ft.Icons.RESTORE,
                        on_click=self._on_reset_handler,
                    ),
                ]
            ),
            alignment=ft.MainAxisAlignment.END,
            spacing=10,
        )

    def sync_to_state(self) -> None:
        """Sync current settings UI controls back to the central state."""
        for k, v in self.settings_map.items():
            self.settings[k] = v.value
        self.appearance_tile.sync_colour_scheme_to_settings()

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central settings."""
        from amplifyp.gui.util import BorderedCheckbox

        for k, v in self.settings.items():
            if k in self.settings_map:
                if isinstance(
                    self.settings_map[k], (ft.Checkbox, BorderedCheckbox)
                ):
                    self.settings_map[k].value = bool(v)
                else:
                    self.settings_map[k].value = str(v)
        self.appearance_tile.update_colour_scheme_dropdown()

    def _on_change_handler(self, e: ft.ControlEvent) -> None:
        """Handle change in settings fields."""
        # Sync the triggering control's value first (Dropdown .value may not
        # be updated yet when on_select fires in Flet).
        ctrl = getattr(e, "control", None)
        if ctrl is not None:
            for k, v in self.settings_map.items():
                if v is ctrl:
                    self.settings[k] = ctrl.value
                    break
        self.sync_to_state()
        if self.on_change:
            self.on_change(e)

    def _on_apply_handler(self, *args: Any) -> None:
        """Handle apply button click."""
        self.sync_to_state()
        if self.on_apply:
            self.on_apply(args[0] if args else None)

    def _on_reset_handler(self, *args: Any) -> None:
        """Handle reset to default button click."""
        from amplifyp.dna import Nucleotides
        from amplifyp.settings import (
            DEFAULT_BASE_PAIR_WEIGHTS,
            DEFAULT_PRIMABILITY_CUTOFF,
            DEFAULT_PRIMER_DIMER_OVERLAP,
            DEFAULT_PRIMER_DIMER_THRESHOLD,
            DEFAULT_PRIMER_DIMER_WEIGHTS,
            DEFAULT_STABILITY_CUTOFF,
            GLOBAL_TM_SETTINGS,
        )

        reset_dict = {
            "primability_cutoff": str(DEFAULT_PRIMABILITY_CUTOFF),
            "stability_cutoff": str(DEFAULT_STABILITY_CUTOFF),
            "amp4_compat": False,
            "tm_method": "SantaLucia 1998 / Owczarzy 2008 (Default)",
            "tm_dna_conc": str(GLOBAL_TM_SETTINGS.dna_conc),
            "tm_dnap_conc": str(GLOBAL_TM_SETTINGS.dnap_conc),
            "tm_mono_salt": str(GLOBAL_TM_SETTINGS.monovalent_salt_conc),
            "tm_div_salt": str(GLOBAL_TM_SETTINGS.divalent_salt_conc),
            "tm_dNTP_conc": str(GLOBAL_TM_SETTINGS.dnTP_conc),
            "pd_min_overlap": str(DEFAULT_PRIMER_DIMER_OVERLAP),
            "pd_threshold": str(DEFAULT_PRIMER_DIMER_THRESHOLD),
            "font_family": "Roboto Mono",
            "colour_deficient": False,
            "dark_mode": "system",
            "improved_visualisation": True,
            "show_primer_temperature": False,
            "tm_colour_scheme": "None",
        }

        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.TEMPLATE:
                if c_char == Nucleotides.GAP:
                    continue
                key = f"bp_score_{r_char}_{c_char}"
                reset_dict[key] = str(DEFAULT_BASE_PAIR_WEIGHTS[r_char, c_char])

        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.PRIMER:
                key = f"pd_score_{r_char}_{c_char}"
                reset_dict[key] = str(
                    DEFAULT_PRIMER_DIMER_WEIGHTS[r_char, c_char]
                )

        self.settings.from_dict(reset_dict)
        self.update_ui()
        self.app_page.update()
        if self.on_reset:
            self.on_reset(args[0] if args else None)

    def get_replication_settings(self) -> ReplicationSettings:
        """Get the current settings as a ReplicationSettings object."""
        self.sync_to_state()
        return self.settings.get_replication_settings()

    def get_state(self) -> dict[str, Any]:
        """Get the current settings for serialisation."""
        self.sync_to_state()
        return self.settings.to_dict()

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current settings from deserialized data."""
        self.settings.from_dict(state)
        self.update_ui()
        self.app_page.update()
