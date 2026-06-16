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

"""Centralized GUI settings and color definitions."""

import os
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import flet as ft
import yaml

from amplifyp.settings import (
    DEFAULT_BASE_PAIR_WEIGHTS,
    DEFAULT_PRIMABILITY_CUTOFF,
    DEFAULT_PRIMER_DIMER_OVERLAP,
    DEFAULT_PRIMER_DIMER_THRESHOLD,
    DEFAULT_PRIMER_DIMER_WEIGHTS,
    DEFAULT_STABILITY_CUTOFF,
    GLOBAL_TM_SETTINGS,
)

# Maximum number of amplicons to draw on the PCR result diagram to
# prevent UI freeze.
MAX_AMPLICONS_RENDER = 100

# Maximum number of primer dimer cards to display in the UI to
# prevent UI freeze.
MAX_DIMERS_RENDER = 100

if TYPE_CHECKING:
    from amplifyp.settings import (
        PrimerDimerSettings,
        ReplicationSettings,
        TMSettings,
    )


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
    def dark_mode(cls) -> bool:
        """Get dark mode status."""
        return cls._dark_mode

    @dark_mode.setter
    def dark_mode(cls, value: bool) -> None:
        """Set dark mode status."""
        cls._dark_mode = value

    @property
    def TEXT_ON_SURFACE(cls) -> str:
        """Get standard text color on surface."""
        return cast(str, ft.Colors.ON_SURFACE)

    @property
    def SUCCESS_GREEN(cls) -> str:
        """Get color-blind friendly success color."""
        # Blue is highly distinguishable for most common
        # red-green color blindness.
        return cast(
            str,
            ft.Colors.BLUE_400
            if cls._color_deficient_mode
            else ft.Colors.GREEN_400,
        )

    @property
    def OUTLINE_VARIANT(cls) -> str:
        """Get outline variant color."""
        return cast(str, ft.Colors.OUTLINE_VARIANT)

    @property
    def OUTLINE(cls) -> str:
        """Get outline color."""
        return cast(str, ft.Colors.OUTLINE)

    @property
    def ERROR_RED(cls) -> str:
        """Get color-blind friendly error color."""
        # Use Orange/Vermilion instead of Red in color deficient mode.
        return cast(
            str,
            ft.Colors.ORANGE_800
            if cls._color_deficient_mode
            else ft.Colors.RED,
        )

    @property
    def DIVIDER_GREY(cls) -> str:
        """Get divider grey color."""
        return cast(str, ft.Colors.GREY_400)

    @property
    def DUPLICATE_BG(cls) -> str:
        """Get duplicate warning background color."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.ORANGE_900
                if cls._color_deficient_mode
                else ft.Colors.RED_900,
            )
        return cast(
            str,
            ft.Colors.ORANGE_100
            if cls._color_deficient_mode
            else ft.Colors.RED_100,
        )

    @property
    def SELECTED_ROW_BG(cls) -> str:
        """Get selected/focused row background color."""
        return cast(
            str,
            ft.Colors.BLUE_900 if cls._dark_mode else ft.Colors.BLUE_50,
        )

    @property
    def DIAGRAM_BLACK(cls) -> str:
        """Get diagram black/white color depending on dark mode."""
        return cast(
            str,
            ft.Colors.WHITE if cls._dark_mode else ft.Colors.BLACK,
        )

    @property
    def INFO_HEADER_BG(cls) -> str:
        """Get background color for info header."""
        return cast(
            str,
            ft.Colors.GREY_800 if cls._dark_mode else ft.Colors.GREY_200,
        )

    @property
    def FWD_PRIMER(cls) -> str:
        """Get forward primer color."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.BLUE_300
                if cls._color_deficient_mode
                else ft.Colors.BLUE_400,
            )
        # Sky blue / clear blue for forward primer in color deficient mode
        return cast(
            str,
            ft.Colors.BLUE_600
            if cls._color_deficient_mode
            else ft.Colors.BLUE_800,
        )

    @property
    def REV_PRIMER(cls) -> str:
        """Get reverse primer color."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.ORANGE_300
                if cls._color_deficient_mode
                else ft.Colors.RED_ACCENT_200,
            )
        # Vermilion / orange-red for reverse primer in color deficient mode
        # (instead of red-accent)
        return cast(
            str,
            ft.Colors.ORANGE_700
            if cls._color_deficient_mode
            else ft.Colors.RED_ACCENT_700,
        )

    @property
    def REV_LABEL(cls) -> str:
        """Get reverse primer label color."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.ORANGE_200
                if cls._color_deficient_mode
                else ft.Colors.RED_300,
            )
        return cast(
            str,
            ft.Colors.ORANGE_900
            if cls._color_deficient_mode
            else ft.Colors.RED_800,
        )

    @property
    def TRANSPARENT(cls) -> str:
        """Get transparent color."""
        return cast(str, ft.Colors.TRANSPARENT)


class GUIColors(metaclass=_GUIColorsMeta):
    """Centralized semantic color constants for the GUI."""

    _color_deficient_mode = False
    _dark_mode = False


class GUISettings:
    """Encapsulates configuration settings for the GUI."""

    def __init__(self, settings_dict: dict[str, Any] | None = None) -> None:
        """Initialize GUISettings with optional dictionary or defaults."""
        from amplifyp.dna import Nucleotides

        self._settings: dict[str, Any] = {
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
            "color_deficient": False,
            "dark_mode": "system",
            "font_size_map_baseline": 16,
            "font_size_map_primer": 13,
            "font_size_map_amplicon": 13,
            "font_size_header": 18,
            "font_size_subheader": 16,
            "font_size_body": 13,
            "font_size_small": 12,
            "font_size_default": 14,
            "font_size_micro": 10,
            "font_size_table_header": 15,
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
            if "dark_mode" in settings_dict:
                val = settings_dict["dark_mode"]
                if isinstance(val, str):
                    val = val.lower() in ("true", "1", "yes")
                GUIColors.dark_mode = bool(val)

    def __getitem__(self, key: str) -> Any:
        """Get a setting value by key."""
        return self._settings[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """Set a setting value by key."""
        if key in self._settings:
            orig_val = self._settings[key]
            if isinstance(orig_val, bool):
                if isinstance(value, str):
                    value = value.lower() in ("true", "1", "yes")
                else:
                    value = bool(value)
            elif isinstance(orig_val, int):
                try:
                    value = int(value)
                except (ValueError, TypeError):
                    pass
            elif isinstance(orig_val, float):
                try:
                    value = float(value)
                except (ValueError, TypeError):
                    pass
        self._settings[key] = value
        if key == "color_deficient":
            val = value
            if isinstance(val, str):
                val = val.lower() in ("true", "1", "yes")
            GUIColors.color_deficient_mode = bool(val)
        elif key == "dark_mode":
            val = value
            if isinstance(val, str) and val.lower() == "system":
                # Handled by apply_theme based on brightness
                pass
            else:
                if isinstance(val, str):
                    val = val.lower() in ("true", "1", "yes")
                GUIColors.dark_mode = bool(val)

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
            primability_cutoff=self._safe_float(
                "primability_cutoff", DEFAULT_PRIMABILITY_CUTOFF
            ),
            stability_cutoff=self._safe_float(
                "stability_cutoff", DEFAULT_STABILITY_CUTOFF
            ),
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
            min_overlap=self._safe_int(
                "pd_min_overlap", DEFAULT_PRIMER_DIMER_OVERLAP
            ),
            threshold=self._safe_float(
                "pd_threshold", DEFAULT_PRIMER_DIMER_THRESHOLD
            ),
            weights=weights,
        )

    def get_tm_settings(self) -> "TMSettings":
        """Get TMSettings from the central settings."""
        from amplifyp.settings import GLOBAL_TM_SETTINGS, TMSettings

        return TMSettings(
            dna_conc=self._safe_float(
                "tm_dna_conc", GLOBAL_TM_SETTINGS.dna_conc
            ),
            dnap_conc=self._safe_float(
                "tm_dnap_conc", GLOBAL_TM_SETTINGS.dnap_conc
            ),
            monovalent_salt_conc=self._safe_float(
                "tm_mono_salt", GLOBAL_TM_SETTINGS.monovalent_salt_conc
            ),
            divalent_salt_conc=self._safe_float(
                "tm_div_salt", GLOBAL_TM_SETTINGS.divalent_salt_conc
            ),
            dnTP_conc=self._safe_float(
                "tm_dNTP_conc", GLOBAL_TM_SETTINGS.dnTP_conc
            ),
        )

    def calculate_primer_tm(self, primer: Any) -> float:
        """Calculate the melting temperature of a primer based on settings."""
        from amplifyp.melting import (
            calculate_tm_lander_amplify4,
            calculate_tm_santalucia_1998_owczarzy_2008,
        )

        tm_method = self.get(
            "tm_method", "SantaLucia 1998 / Owczarzy 2008 (Default)"
        )
        tm_settings = self.get_tm_settings()
        if tm_method == "Lander / Amplify 4":
            return calculate_tm_lander_amplify4(primer, tm_settings)
        return calculate_tm_santalucia_1998_owczarzy_2008(primer, tm_settings)

    def to_dict(self) -> dict[str, Any]:
        """Convert settings to a dictionary for serialization."""
        return dict(self._settings)

    def from_dict(self, settings_dict: dict[str, Any]) -> None:
        """Update settings from a dictionary."""
        for k in self._settings:
            if k in settings_dict:
                v = settings_dict[k]
                if k == "dark_mode":
                    if str(v).lower() == "system":
                        self._settings[k] = "system"
                    else:
                        if isinstance(v, str):
                            self._settings[k] = v.lower() in (
                                "true",
                                "1",
                                "yes",
                            )
                        else:
                            self._settings[k] = bool(v)
                elif isinstance(self._settings[k], bool):
                    self._settings[k] = bool(v)
                elif isinstance(self._settings[k], int):
                    try:
                        self._settings[k] = int(v)
                    except (ValueError, TypeError):
                        pass
                elif isinstance(self._settings[k], float):
                    try:
                        self._settings[k] = float(v)
                    except (ValueError, TypeError):
                        pass
                else:
                    self._settings[k] = str(v)
        GUIColors.color_deficient_mode = bool(
            self._settings.get("color_deficient", False)
        )
        dark_mode_val = self._settings.get("dark_mode", False)
        if str(dark_mode_val).lower() != "system":
            GUIColors.dark_mode = bool(dark_mode_val)

    def _get_config_path(self) -> Path:
        """Get the OS-specific path for user settings configuration."""
        if sys.platform.startswith("win"):
            appdata = os.environ.get("APPDATA")
            if appdata:
                return Path(appdata) / "AmplifyP" / "settings.yaml"
            return (
                Path(os.path.expanduser("~"))
                / "AppData"
                / "Roaming"
                / "AmplifyP"
                / "settings.yaml"
            )
        elif sys.platform.startswith("darwin"):
            home = os.environ.get("HOME") or os.path.expanduser("~")
            return (
                Path(home)
                / "Library"
                / "Application Support"
                / "AmplifyP"
                / "settings.yaml"
            )
        else:
            xdg_config = os.environ.get("XDG_CONFIG_HOME")
            if xdg_config:
                return Path(xdg_config) / "amplifyp" / "settings.yaml"
            home = os.environ.get("HOME") or os.path.expanduser("~")
            return Path(home) / ".config" / "amplifyp" / "settings.yaml"

    def load_from_local(self, page: ft.Page) -> None:
        """Load settings from local storage.

        On web platforms, uses page.client_storage.
        On desktop platforms, reads from OS-specific YAML file.
        """
        if getattr(page, "web", False):
            if (
                not hasattr(page, "client_storage")
                or page.client_storage is None
            ):
                return
            # Web browser path
            settings_dict = {}
            for k in self._settings.keys():
                storage_key = f"amplifyp.settings.{k}"
                if page.client_storage.contains_key(storage_key):
                    settings_dict[k] = page.client_storage.get(storage_key)
            if settings_dict:
                self.from_dict(settings_dict)
        else:
            # Desktop path
            config_path = self._get_config_path()
            if config_path.exists():
                try:
                    with open(config_path, encoding="utf-8") as f:
                        data = yaml.safe_load(f)
                    if isinstance(data, dict):
                        self.from_dict(data)
                except Exception as e:
                    print(f"Error loading settings from {config_path}: {e}")

    def save_to_local(self, page: ft.Page) -> None:
        """Save settings to local storage.

        On web platforms, uses page.client_storage.
        On desktop platforms, writes to OS-specific YAML file.
        """
        if getattr(page, "web", False):
            if (
                not hasattr(page, "client_storage")
                or page.client_storage is None
            ):
                return
            # Web browser path
            for k, v in self._settings.items():
                storage_key = f"amplifyp.settings.{k}"
                page.client_storage.set(storage_key, v)
        else:
            # Desktop path
            config_path = self._get_config_path()
            try:
                config_path.parent.mkdir(parents=True, exist_ok=True)
                with open(config_path, "w", encoding="utf-8") as f:
                    yaml.safe_dump(self.to_dict(), f, default_flow_style=False)
            except Exception as e:
                print(f"Error saving settings to {config_path}: {e}")
