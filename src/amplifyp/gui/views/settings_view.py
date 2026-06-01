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

from typing import Any

import flet as ft

from amplifyp.gui.state import GUIState
from amplifyp.settings import ReplicationSettings


class SettingsView(ft.ListView):  # type: ignore[misc]
    """Settings view for configuring PCR replication settings."""

    def __init__(self, page: ft.Page, state: GUIState | None = None) -> None:
        """Initialize the SettingsView."""
        super().__init__(expand=True, spacing=20, padding=10)
        self.app_page = page
        self.state = state if state is not None else GUIState()

        # Settings State
        # Replication Settings
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

        # TM Settings
        self.set_tm_dna_conc = ft.TextField(
            label="DNA Conc (nM)",
            value="50.0",
            on_change=self.on_change_handler,
        )
        self.set_tm_dnap_conc = ft.TextField(
            label="DNA Pol Conc", value="0.0", on_change=self.on_change_handler
        )
        self.set_tm_mono_salt = ft.TextField(
            label="Monovalent Salt Conc (mM)",
            value="50.0",
            on_change=self.on_change_handler,
        )
        self.set_tm_div_salt = ft.TextField(
            label="Divalent Salt Conc (mM)",
            value="1.5",
            on_change=self.on_change_handler,
        )
        self.set_tm_dNTP_conc = ft.TextField(
            label="dNTP Conc (mM)",
            value="0.0",
            on_change=self.on_change_handler,
        )

        # Amplify4 TM Settings
        self.set_amp4tm_dna_conc = ft.TextField(
            label="Amplify4 DNA Conc (nM)",
            value="50.0",
            on_change=self.on_change_handler,
        )
        self.set_amp4tm_mono_salt = ft.TextField(
            label="Amplify4 Monovalent Salt Conc (mM)",
            value="50.0",
            on_change=self.on_change_handler,
        )

        # Primer Dimer Settings
        self.set_pd_min_overlap = ft.TextField(
            label="Min Overlap", value="3", on_change=self.on_change_handler
        )
        self.set_pd_threshold = ft.TextField(
            label="Threshold", value="60.0", on_change=self.on_change_handler
        )

        self.settings_map = {
            "primability_cutoff": self.set_primability_cutoff,
            "stability_cutoff": self.set_stability_cutoff,
            "amp4_compat": self.set_amp4_compat,
            "tm_dna_conc": self.set_tm_dna_conc,
            "tm_dnap_conc": self.set_tm_dnap_conc,
            "tm_mono_salt": self.set_tm_mono_salt,
            "tm_div_salt": self.set_tm_div_salt,
            "tm_dNTP_conc": self.set_tm_dNTP_conc,
            "amp4tm_dna_conc": self.set_amp4tm_dna_conc,
            "amp4tm_mono_salt": self.set_amp4tm_mono_salt,
            "pd_min_overlap": self.set_pd_min_overlap,
            "pd_threshold": self.set_pd_threshold,
        }

        self.controls = [
            ft.Text("Replication Settings", weight=ft.FontWeight.BOLD, size=18),
            self.set_primability_cutoff,
            self.set_stability_cutoff,
            self.set_amp4_compat,
            ft.Divider(),
            ft.Text("TM Settings", weight=ft.FontWeight.BOLD, size=18),
            self.set_tm_dna_conc,
            self.set_tm_dnap_conc,
            self.set_tm_mono_salt,
            self.set_tm_div_salt,
            self.set_tm_dNTP_conc,
            ft.Divider(),
            ft.Text("Amplify4 TM Settings", weight=ft.FontWeight.BOLD, size=18),
            self.set_amp4tm_dna_conc,
            self.set_amp4tm_mono_salt,
            ft.Divider(),
            ft.Text(
                "Primer Dimer Settings", weight=ft.FontWeight.BOLD, size=18
            ),
            self.set_pd_min_overlap,
            self.set_pd_threshold,
        ]

        # Sync initial UI state
        self.update_ui()

    def sync_to_state(self) -> None:
        """Sync current settings UI controls back to the central state."""
        for k, v in self.settings_map.items():
            self.state.settings[k] = v.value

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        for k, v in self.state.settings.items():
            if k in self.settings_map:
                if isinstance(self.settings_map[k], ft.Checkbox):
                    self.settings_map[k].value = bool(v)
                else:
                    self.settings_map[k].value = str(v)

    def on_change_handler(self, e: ft.ControlEvent) -> None:
        """Handle change in settings fields."""
        self.sync_to_state()
        if self.state.has_run_pcr:
            self.state.results_outdated = True

    def get_replication_settings(self) -> ReplicationSettings:
        """Get the current settings as a ReplicationSettings object."""
        self.sync_to_state()
        return self.state.get_replication_settings()

    def get_state(self) -> dict[str, Any]:
        """Get the current state for serialization."""
        self.sync_to_state()
        return dict(self.state.settings)

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current state from deserialized data."""
        self.state.from_dict(state)
        self.update_ui()
        self.app_page.update()
