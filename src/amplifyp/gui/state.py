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

import flet as ft

from amplifyp.gui.util import clean_sequence, format_sequence

if TYPE_CHECKING:
    from amplifyp.settings import PrimerDimerSettings, ReplicationSettings


class GUIColors:
    """Centralized semantic color constants for the GUI."""

    TEXT_ON_SURFACE = ft.Colors.ON_SURFACE
    SUCCESS_GREEN = ft.Colors.GREEN_400
    CONTAINER_HIGHEST = ft.Colors.SURFACE_CONTAINER_HIGHEST
    OUTLINE_VARIANT = ft.Colors.OUTLINE_VARIANT
    OUTLINE = ft.Colors.OUTLINE
    ERROR_RED = ft.Colors.RED
    DIVIDER_GREY = ft.Colors.GREY_400
    DUPLICATE_BG = ft.Colors.RED_100
    DIAGRAM_BLACK = ft.Colors.BLACK
    FWD_PRIMER = ft.Colors.BLUE_800
    REV_PRIMER = ft.Colors.RED_ACCENT_700
    REV_LABEL = ft.Colors.RED_800
    TRANSPARENT = ft.Colors.TRANSPARENT


class GUISettings:
    """Encapsulates configuration settings for the GUI."""

    def __init__(self, settings_dict: dict[str, Any] | None = None) -> None:
        """Initialize GUISettings with optional dictionary or defaults."""
        if settings_dict is not None:
            self._settings = dict(settings_dict)
        else:
            from amplifyp.settings import (
                DEFAULT_PRIMABILITY_CUTOFF,
                DEFAULT_PRIMER_DIMER_OVERLAP,
                DEFAULT_PRIMER_DIMER_THRESHOLD,
                DEFAULT_STABILITY_CUTOFF,
                GLOBAL_AMPLIFY4_TM_SETTINGS,
                GLOBAL_TM_SETTINGS,
            )

            self._settings = {
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
                "font_family": "Roboto Mono",
            }

    def __getitem__(self, key: str) -> Any:
        """Get a setting value by key."""
        return self._settings[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """Set a setting value by key."""
        self._settings[key] = value

    def get(self, key: str, default: Any = None) -> Any:
        """Get a setting value by key with an optional default."""
        return self._settings.get(key, default)

    def items(self) -> Any:
        """Get the settings key-value items."""
        return self._settings.items()

    def keys(self) -> Any:
        """Get the settings keys."""
        return self._settings.keys()

    def __iter__(self) -> Any:
        """Iterate over settings keys."""
        return iter(self._settings)

    def __contains__(self, key: str) -> bool:
        """Check if a key exists in settings."""
        return key in self._settings

    def _safe_float(self, key: str, default: float) -> float:
        val = self._settings.get(key)
        if val is None or val == "":
            return default
        try:
            return float(val)
        except (ValueError, TypeError):
            return default

    def _safe_int(self, key: str, default: int) -> int:
        val = self._settings.get(key)
        if val is None or val == "":
            return default
        try:
            return int(val)
        except (ValueError, TypeError):
            return default

    def get_replication_settings(self) -> "ReplicationSettings":
        """Get ReplicationSettings from the central settings."""
        from amplifyp.settings import ReplicationSettings

        return ReplicationSettings(
            primability_cutoff=self._safe_float("primability_cutoff", 0.8),
            stability_cutoff=self._safe_float("stability_cutoff", 0.4),
            amplify4_compatibility_mode=self._safe_int("amp4_compat", 0) != 0,
        )

    def get_primer_dimer_settings(self) -> "PrimerDimerSettings":
        """Get PrimerDimerSettings from the central settings."""
        from amplifyp.settings import PrimerDimerSettings

        return PrimerDimerSettings(
            min_overlap=self._safe_int("pd_min_overlap", 3),
            threshold=self._safe_float("pd_threshold", 60.0),
        )


class GUIInput:
    """Encapsulates user input including DNA template and primers list."""

    def __init__(
        self,
        template: str = "",
        template_circular: bool = False,
        primers: list[dict[str, Any]] | None = None,
    ) -> None:
        """Initialize GUIInput."""
        self.template: str = template
        self.template_circular: bool = template_circular
        self.primers: list[dict[str, Any]] = (
            primers if primers is not None else []
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


class GUIState:
    """Central store for all states of the GUI."""

    def __init__(self) -> None:
        """Initialize default GUI state."""
        self.input = GUIInput()
        self._settings = GUISettings()
        self.results_outdated: bool = False
        self.has_run_pcr: bool = False

    @property
    def template(self) -> str:
        """Get template DNA sequence."""
        return self.input.template

    @template.setter
    def template(self, value: str) -> None:
        """Set template DNA sequence."""
        self.input.template = value

    @property
    def template_circular(self) -> bool:
        """Get whether template is circular."""
        return self.input.template_circular

    @template_circular.setter
    def template_circular(self, value: bool) -> None:
        """Set whether template is circular."""
        self.input.template_circular = value

    @property
    def primers(self) -> list[dict[str, Any]]:
        """Get primers list."""
        return self.input.primers

    @primers.setter
    def primers(self, value: list[dict[str, Any]]) -> None:
        """Set primers list."""
        self.input.primers = value

    @property
    def settings(self) -> GUISettings:
        """Get GUISettings."""
        return self._settings

    @settings.setter
    def settings(self, value: Any) -> None:
        """Set GUISettings, coercing dictionaries to GUISettings."""
        if isinstance(value, GUISettings):
            self._settings = value
        elif isinstance(value, dict):
            self._settings = GUISettings(value)
        else:
            raise TypeError("settings must be a GUISettings or dict")

    def get_replication_settings(self) -> "ReplicationSettings":
        """Get ReplicationSettings from the central settings."""
        return self.settings.get_replication_settings()

    def get_primer_dimer_settings(self) -> "PrimerDimerSettings":
        """Get PrimerDimerSettings from the central settings."""
        return self.settings.get_primer_dimer_settings()

    def get_active_primers(self) -> list[dict[str, Any]]:
        """Get the active primers."""
        return self.input.get_active_primers()

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
            "settings": dict(self.settings._settings),
        }

    def from_dict(self, state_dict: dict[str, Any]) -> None:
        """Update state from a dictionary."""
        if "template" in state_dict:
            self.template = clean_sequence(state_dict["template"])
        if "template_circular" in state_dict:
            self.template_circular = bool(state_dict["template_circular"])
        if "primers" in state_dict:
            self.input._load_primers(state_dict["primers"])
        if "settings" in state_dict:
            self._load_settings(state_dict["settings"])

    def _load_settings(self, settings: dict[str, Any]) -> None:
        """Load settings from a dict, coercing types to match self.settings."""
        for k in self.settings._settings:
            if k in settings:
                v = settings[k]
                self.settings._settings[k] = (
                    bool(v)
                    if isinstance(self.settings._settings[k], bool)
                    else str(v)
                )
