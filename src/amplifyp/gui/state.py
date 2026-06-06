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


class _GUIColorsMeta(type):
    """Metaclass for GUIColors to support dynamic color shifting."""

    @property
    def color_deficient_mode(cls) -> bool:
        """Get color deficient mode status."""
        return cls._color_deficient_mode

    @color_deficient_mode.setter
    def color_deficient_mode(cls, value: bool) -> None:
        """Set color deficient mode status."""
        cls._color_deficient_mode = value

    @property
    def TEXT_ON_SURFACE(cls) -> str:
        """Get standard text color on surface."""
        return str(ft.Colors.ON_SURFACE.value)

    @property
    def SUCCESS_GREEN(cls) -> str:
        """Get color-blind friendly success color."""
        # Blue is highly distinguishable for most common
        # red-green color blindness.
        return (
            str(ft.Colors.BLUE_400.value)
            if cls._color_deficient_mode
            else str(ft.Colors.GREEN_400.value)
        )

    @property
    def CONTAINER_HIGHEST(cls) -> str:
        """Get container highest color."""
        return str(ft.Colors.SURFACE_CONTAINER_HIGHEST.value)

    @property
    def OUTLINE_VARIANT(cls) -> str:
        """Get outline variant color."""
        return str(ft.Colors.OUTLINE_VARIANT.value)

    @property
    def OUTLINE(cls) -> str:
        """Get outline color."""
        return str(ft.Colors.OUTLINE.value)

    @property
    def ERROR_RED(cls) -> str:
        """Get color-blind friendly error color."""
        # Use Orange/Vermilion instead of Red in color deficient mode.
        return (
            str(ft.Colors.ORANGE_800.value)
            if cls._color_deficient_mode
            else str(ft.Colors.RED.value)
        )

    @property
    def DIVIDER_GREY(cls) -> str:
        """Get divider grey color."""
        return str(ft.Colors.GREY_400.value)

    @property
    def DUPLICATE_BG(cls) -> str:
        """Get duplicate warning background color."""
        return (
            str(ft.Colors.ORANGE_100.value)
            if cls._color_deficient_mode
            else str(ft.Colors.RED_100.value)
        )

    @property
    def DIAGRAM_BLACK(cls) -> str:
        """Get diagram black color."""
        return str(ft.Colors.BLACK.value)

    @property
    def FWD_PRIMER(cls) -> str:
        """Get forward primer color."""
        # Sky blue / clear blue for forward primer in color deficient mode
        return (
            str(ft.Colors.BLUE_600.value)
            if cls._color_deficient_mode
            else str(ft.Colors.BLUE_800.value)
        )

    @property
    def REV_PRIMER(cls) -> str:
        """Get reverse primer color."""
        # Vermilion / orange-red for reverse primer in color deficient mode
        # (instead of red-accent)
        return (
            str(ft.Colors.ORANGE_700.value)
            if cls._color_deficient_mode
            else str(ft.Colors.RED_ACCENT_700.value)
        )

    @property
    def REV_LABEL(cls) -> str:
        """Get reverse primer label color."""
        return (
            str(ft.Colors.ORANGE_900.value)
            if cls._color_deficient_mode
            else str(ft.Colors.RED_800.value)
        )

    @property
    def TRANSPARENT(cls) -> str:
        """Get transparent color."""
        return str(ft.Colors.TRANSPARENT.value)


class GUIColors(metaclass=_GUIColorsMeta):
    """Centralized semantic color constants for the GUI."""

    _color_deficient_mode = False


class GUISettings:
    """Encapsulates configuration settings for the GUI."""

    def __init__(self, settings_dict: dict[str, Any] | None = None) -> None:
        """Initialize GUISettings with optional dictionary or defaults."""
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

        self._settings: dict[str, Any] = {
            "primability_cutoff": str(DEFAULT_PRIMABILITY_CUTOFF),
            "stability_cutoff": str(DEFAULT_STABILITY_CUTOFF),
            "amp4_compat": False,
            "tm_method": "Amplify P (Default)",
            "tm_dna_conc": str(GLOBAL_TM_SETTINGS.dna_conc),
            "tm_dnap_conc": str(GLOBAL_TM_SETTINGS.dnap_conc),
            "tm_mono_salt": str(GLOBAL_TM_SETTINGS.monovalent_salt_conc),
            "tm_div_salt": str(GLOBAL_TM_SETTINGS.divalent_salt_conc),
            "tm_dNTP_conc": str(GLOBAL_TM_SETTINGS.dnTP_conc),
            "pd_min_overlap": str(DEFAULT_PRIMER_DIMER_OVERLAP),
            "pd_threshold": str(DEFAULT_PRIMER_DIMER_THRESHOLD),
            "font_family": "Roboto Mono",
            "color_deficient": False,
            "font_size_map_baseline": 16,
            "font_size_map_primer": 13,
            "font_size_map_amplicon": 13,
            "font_size_header": 18,
            "font_size_subheader": 16,
            "font_size_body": 13,
            "font_size_small": 12,
            "font_size_default": 14,
        }

        # Initialize base-pair weights
        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.TEMPLATE:
                if c_char == Nucleotides.GAP:
                    continue
                key = f"bp_score_{r_char}_{c_char}"
                self._settings[key] = str(
                    DEFAULT_BASE_PAIR_WEIGHTS[r_char, c_char]
                )

        # Initialize primer-dimer weights
        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.PRIMER:
                key = f"pd_score_{r_char}_{c_char}"
                self._settings[key] = str(
                    DEFAULT_PRIMER_DIMER_WEIGHTS[r_char, c_char]
                )

        if settings_dict is not None:
            self._settings.update(settings_dict)
            if "color_deficient" in settings_dict:
                val = settings_dict["color_deficient"]
                if isinstance(val, str):
                    val = val.lower() in ("true", "1", "yes")
                GUIColors.color_deficient_mode = bool(val)

    def __getitem__(self, key: str) -> Any:
        """Get a setting value by key."""
        return self._settings[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """Set a setting value by key."""
        self._settings[key] = value
        if key == "color_deficient":
            val = value
            if isinstance(val, str):
                val = val.lower() in ("true", "1", "yes")
            GUIColors.color_deficient_mode = bool(val)

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
        from amplifyp.dna import Nucleotides
        from amplifyp.settings import BasePairWeightsTbl, ReplicationSettings

        bp_weights = []
        for r_char in Nucleotides.PRIMER:
            row_vals = []
            for c_char in Nucleotides.TEMPLATE:
                if c_char == Nucleotides.GAP:
                    continue
                key = f"bp_score_{r_char}_{c_char}"
                val = self._safe_float(key, 0.0)
                row_vals.append(val)
            bp_weights.append(row_vals)

        base_pair_scores = BasePairWeightsTbl(
            row=Nucleotides.PRIMER,
            col=Nucleotides.TEMPLATE,
            weight=bp_weights,
        )

        return ReplicationSettings(
            primability_cutoff=self._safe_float("primability_cutoff", 0.8),
            stability_cutoff=self._safe_float("stability_cutoff", 0.4),
            amplify4_compatibility_mode=self._safe_int("amp4_compat", 0) != 0,
            base_pair_scores=base_pair_scores,
        )

    def get_primer_dimer_settings(self) -> "PrimerDimerSettings":
        """Get PrimerDimerSettings from the central settings."""
        from amplifyp.dna import Nucleotides
        from amplifyp.settings import BasePairWeightsTbl, PrimerDimerSettings

        pd_weights = []
        for r_char in Nucleotides.PRIMER:
            row_vals = []
            for c_char in Nucleotides.PRIMER:
                key = f"pd_score_{r_char}_{c_char}"
                val = self._safe_float(key, 0.0)
                row_vals.append(val)
            pd_weights.append(row_vals)

        weights = BasePairWeightsTbl(
            row=Nucleotides.PRIMER,
            col=Nucleotides.PRIMER,
            weight=pd_weights,
        )

        return PrimerDimerSettings(
            min_overlap=self._safe_int("pd_min_overlap", 3),
            threshold=self._safe_float("pd_threshold", 60.0),
            weights=weights,
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
        GUIColors.color_deficient_mode = bool(
            self.settings._settings.get("color_deficient", False)
        )
