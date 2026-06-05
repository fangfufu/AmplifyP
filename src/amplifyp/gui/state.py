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

"""Centralized GUI state store class."""

from typing import TYPE_CHECKING, Any

from amplifyp.gui.util import clean_sequence, format_sequence

if TYPE_CHECKING:
    from amplifyp.settings import PrimerDimerSettings, ReplicationSettings


class GUIState:
    """Central store for all states of the GUI."""

    def __init__(self) -> None:
        """Initialize default GUI state."""
        from amplifyp.settings import (
            DEFAULT_PRIMABILITY_CUTOFF,
            DEFAULT_PRIMER_DIMER_OVERLAP,
            DEFAULT_PRIMER_DIMER_THRESHOLD,
            DEFAULT_STABILITY_CUTOFF,
            GLOBAL_AMPLIFY4_TM_SETTINGS,
            GLOBAL_TM_SETTINGS,
        )

        self.template: str = ""
        self.template_circular: bool = False
        self.primers: list[dict[str, Any]] = []
        self.settings: dict[str, Any] = {
            "primability_cutoff": str(DEFAULT_PRIMABILITY_CUTOFF),
            "stability_cutoff": str(DEFAULT_STABILITY_CUTOFF),
            "amp4_compat": False,
            "tm_dna_conc": str(GLOBAL_TM_SETTINGS.dna_conc),
            "tm_dnap_conc": str(GLOBAL_TM_SETTINGS.dnap_conc),
            "tm_mono_salt": str(GLOBAL_TM_SETTINGS.monovalent_salt_conc),
            "tm_div_salt": str(GLOBAL_TM_SETTINGS.divalent_salt_conc),
            "tm_dNTP_conc": str(GLOBAL_TM_SETTINGS.dnTP_conc),
            "amp4tm_dna_conc": str(GLOBAL_AMPLIFY4_TM_SETTINGS.dna_conc),
            "amp4tm_mono_salt": str(
                GLOBAL_AMPLIFY4_TM_SETTINGS.monovalent_salt_conc
            ),
            "pd_min_overlap": str(DEFAULT_PRIMER_DIMER_OVERLAP),
            "pd_threshold": str(DEFAULT_PRIMER_DIMER_THRESHOLD),
        }
        self.results_outdated: bool = False
        self.has_run_pcr: bool = False

    def _safe_float(self, key: str, default: float) -> float:
        val = self.settings.get(key)
        if val is None or val == "":
            return default
        try:
            return float(val)
        except (ValueError, TypeError):
            return default

    def _safe_int(self, key: str, default: int) -> int:
        val = self.settings.get(key)
        if val is None or val == "":
            return default
        try:
            return int(val)
        except (ValueError, TypeError):
            return default

    def get_replication_settings(self) -> "ReplicationSettings":
        """Get ReplicationSettings from the central state settings."""
        from amplifyp.settings import ReplicationSettings

        return ReplicationSettings(
            primability_cutoff=self._safe_float("primability_cutoff", 0.8),
            stability_cutoff=self._safe_float("stability_cutoff", 0.4),
            amplify4_compatibility_mode=self._safe_int("amp4_compat", 0) != 0,
        )

    def get_primer_dimer_settings(self) -> "PrimerDimerSettings":
        """Get PrimerDimerSettings from the central state settings."""
        from amplifyp.settings import PrimerDimerSettings

        return PrimerDimerSettings(
            min_overlap=self._safe_int("pd_min_overlap", 3),
            threshold=self._safe_float("pd_threshold", 60.0),
        )

    def get_active_primers(self) -> list[dict[str, Any]]:
        """Get the active primers."""
        return [
            p
            for p in self.primers
            if p.get("active", True)
            and str(p.get("name") or "").strip()
            and clean_sequence(str(p.get("seq") or "")).strip()
        ]

    def to_dict(self) -> dict[str, Any]:
        """Convert state to a dictionary for serialization."""
        formatted_primers = []
        for p in self.primers:
            formatted_primers.append(
                {
                    "name": str(p.get("name", "")),
                    "seq": format_sequence(str(p.get("seq", ""))),
                    "active": bool(p.get("active", True)),
                }
            )

        return {
            "template": format_sequence(self.template),
            "template_circular": self.template_circular,
            "primers": formatted_primers,
            "settings": dict(self.settings),
        }

    def from_dict(self, state_dict: dict[str, Any]) -> None:
        """Update state from a dictionary."""
        if "template" in state_dict:
            self.template = clean_sequence(state_dict["template"])
        if "template_circular" in state_dict:
            self.template_circular = bool(state_dict["template_circular"])
        if "primers" in state_dict:
            self._load_primers(state_dict["primers"])
        if "settings" in state_dict:
            self._load_settings(state_dict["settings"])

    def _load_primers(self, primers: list[dict[str, Any]]) -> None:
        """Load primer data from a list of primer dicts."""
        self.primers = []
        for p in primers:
            p_seq = p.get("seq", "")
            seq_str = "".join(p_seq) if isinstance(p_seq, list) else p_seq
            self.primers.append(
                {
                    "name": str(p.get("name", "")),
                    "seq": clean_sequence(seq_str),
                    "active": bool(p.get("active", True)),
                }
            )

    def _load_settings(self, settings: dict[str, Any]) -> None:
        """Load settings from a dict, coercing types to match self.settings."""
        for k in self.settings:
            if k in settings:
                v = settings[k]
                self.settings[k] = (
                    bool(v) if isinstance(self.settings[k], bool) else str(v)
                )
