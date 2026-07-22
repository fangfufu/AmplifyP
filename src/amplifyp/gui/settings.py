# Copyright (C) 2026 AmplifyP Contributors
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

"""Centralised GUI settings and configuration."""

import logging
import os
import sys
from collections.abc import ItemsView, Iterator, KeysView
from pathlib import Path
from typing import TYPE_CHECKING, Any

import flet as ft
import yaml

from amplifyp.gui.colours import GUIColours
from amplifyp.settings import (
    DEFAULT_BASE_PAIR_WEIGHTS,
    DEFAULT_PRIMABILITY_CUTOFF,
    DEFAULT_PRIMER_DIMER_OVERLAP,
    DEFAULT_PRIMER_DIMER_THRESHOLD,
    DEFAULT_PRIMER_DIMER_WEIGHTS,
    DEFAULT_STABILITY_CUTOFF,
    GLOBAL_TM_SETTINGS,
)

if TYPE_CHECKING:
    from amplifyp.dna import Primer
    from amplifyp.settings import (
        PrimerDimerSettings,
        ReplicationSettings,
        TMSettings,
    )

logger = logging.getLogger(__name__)

# Maximum number of amplicons to draw on the PCR result diagram to
# prevent UI freeze.
MAX_AMPLICONS_RENDER = 100

# Maximum number of primer dimer cards to display in the UI to
# prevent UI freeze.
MAX_DIMERS_RENDER = 100


class GUISettings:
    """Encapsulates configuration settings for the GUI."""

    def __init__(self, settings_dict: dict[str, Any] | None = None) -> None:
        """Initialize GUISettings with optional dictionary or defaults.

        Args:
            settings_dict: Optional dictionary of initial setting values.
        """
        from amplifyp.dna import Nucleotides

        self._cached_tm_settings: TMSettings | None = None
        self._settings: dict[str, Any] = {
            "primability_cutoff": str(DEFAULT_PRIMABILITY_CUTOFF),
            "stability_cutoff": str(DEFAULT_STABILITY_CUTOFF),
            "amp4_compat": False,
            "tm_method": "SantaLucia 1998 / Owczarzy 2008 (Default)",
            "tm_dna_conc": str(GLOBAL_TM_SETTINGS.dna_conc),
            "tm_dnap_conc": str(GLOBAL_TM_SETTINGS.dnap_conc),
            "tm_mono_salt": str(GLOBAL_TM_SETTINGS.monovalent_salt_conc),
            "tm_div_salt": str(GLOBAL_TM_SETTINGS.divalent_salt_conc),
            "tm_dNTP_conc": str(GLOBAL_TM_SETTINGS.dntp_conc),
            "pd_min_overlap": str(DEFAULT_PRIMER_DIMER_OVERLAP),
            "pd_threshold": str(DEFAULT_PRIMER_DIMER_THRESHOLD),
            "font_family": "Roboto Mono",
            "colour_deficient": False,
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
            "improved_visualisation": True,
            "show_primer_temperature": False,
            "ignore_inactive_name_dup_warn": True,
            "ignore_inactive_seq_dup_warn": True,
            "tm_colour_scheme": "None",
            "log_level_amplifyp": "INFO",
            "log_level_flet": "INFO",
            "log_console_enabled": True,
            "log_file_enabled": True,
            "log_file_path": "(Default)",
            "log_rotation_enabled": True,
            "log_rotation_max_bytes": 5242880,
            "version_checking_frequency": "Once per Month",
            "last_version_check_timestamp": 0.0,
            "auto_reload_on_startup": True,
            "template_fixed_width": False,
            "template_bases_per_line": "Auto",
            "primer_info_panel_position": "bottom",
            "primer_info_panel_fixed_height": False,
            "auto_activate_new_valid_primer": False,
        }

        # Initialise base-pair weights
        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.TEMPLATE:
                if c_char == Nucleotides.GAP:
                    continue
                key = f"bp_score_{r_char}_{c_char}"
                self._settings[key] = str(
                    DEFAULT_BASE_PAIR_WEIGHTS[r_char, c_char]
                )

        # Initialise primer-dimer weights
        for r_char in Nucleotides.PRIMER:
            for c_char in Nucleotides.PRIMER:
                key = f"pd_score_{r_char}_{c_char}"
                self._settings[key] = str(
                    DEFAULT_PRIMER_DIMER_WEIGHTS[r_char, c_char]
                )

        if settings_dict is not None:
            self._settings.update(settings_dict)
            if "colour_deficient" in settings_dict:
                val = settings_dict["colour_deficient"]
                if isinstance(val, str):
                    val = val.lower() in ("true", "1", "yes")
                GUIColours.colour_deficient_mode = bool(val)
            if "dark_mode" in settings_dict:
                val = settings_dict["dark_mode"]
                if isinstance(val, str):
                    val = val.lower() in ("true", "1", "yes")
                GUIColours.dark_mode = bool(val)

    def __getitem__(self, key: str) -> Any:
        """Get a setting value by key.

        Args:
            key: The setting key to retrieve.

        Returns:
            The setting value.

        Raises:
            KeyError: If the key does not exist.
        """
        return self._settings[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """Set a setting value by key with type coercion.

        Automatically converts string values to match the expected type
        of the setting (bool, int, float, or str).

        Args:
            key: The setting key to update.
            value: The new value to set.
        """
        if key in self._settings:
            orig_val = self._settings[key]
            if key == "dark_mode":
                if isinstance(value, str):
                    if value.lower() == "system":
                        value = "system"
                    else:
                        value = value.lower() in ("true", "1", "yes")
                else:
                    value = bool(value)
            elif isinstance(orig_val, bool):
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
        self._cached_tm_settings = None
        if key == "colour_deficient":
            val = value
            if isinstance(val, str):
                val = val.lower() in ("true", "1", "yes")
            GUIColours.colour_deficient_mode = bool(val)
        elif key == "dark_mode":
            val = value
            if isinstance(val, str) and val.lower() == "system":
                # Handled by apply_theme based on brightness
                pass
            else:
                if isinstance(val, str):
                    val = val.lower() in ("true", "1", "yes")
                GUIColours.dark_mode = bool(val)

    def get(self, key: str, default: Any = None) -> Any:
        """Get a setting value by key with an optional default.

        Args:
            key: The setting key to retrieve.
            default: The default value if the key does not exist.

        Returns:
            The setting value, or the default if not found.
        """
        return self._settings.get(key, default)

    def items(self) -> ItemsView[str, Any]:
        """Get the settings key-value items.

        Returns:
            A view of the setting key-value pairs.
        """
        return self._settings.items()

    def keys(self) -> KeysView[str]:
        """Get the settings keys.

        Returns:
            A view of the setting keys.
        """
        return self._settings.keys()

    def __iter__(self) -> Iterator[str]:
        """Iterate over settings keys.

        Returns:
            An iterator over the setting keys.
        """
        return iter(self._settings)

    def __contains__(self, key: str) -> bool:
        """Check if a key exists in settings.

        Args:
            key: The setting key to check.

        Returns:
            True if the key exists, False otherwise.
        """
        return key in self._settings

    def _safe_float(self, key: str, default: float) -> float:
        """Retrieve a setting as float, returning a default if conversion fails.

        Args:
            key (str): The setting key.
            default (float): The default value to return.

        Returns:
            float: The setting value as a float, or the default value.
        """
        val = self._settings.get(key)
        if val is None or val == "":
            return default
        try:
            return float(val)
        except (ValueError, TypeError):
            return default

    def _safe_int(self, key: str, default: int) -> int:
        """Retrieve a setting as int, returning a default if conversion fails.

        Args:
            key (str): The setting key.
            default (int): The default value to return.

        Returns:
            int: The setting value as an integer, or the default value.
        """
        val = self._settings.get(key)
        if val is None or val == "":
            return default
        try:
            return int(val)
        except (ValueError, TypeError):
            return default

    def get_replication_settings(self) -> "ReplicationSettings":
        """Get ReplicationSettings from the central settings.

        Constructs a ReplicationSettings object using the current GUI
        settings values for primability/stability cutoffs, base-pair
        weights, and Amplify4 compatibility mode.

        Returns:
            A ReplicationSettings instance configured from GUI settings.
        """
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
        """Get PrimerDimerSettings from the central settings.

        Constructs a PrimerDimerSettings object using the current GUI
        settings values for minimum overlap, threshold, and base-pair
        weights.

        Returns:
            A PrimerDimerSettings instance configured from GUI settings.
        """
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
        """Get TMSettings from the central settings.

        Constructs a TMSettings object using the current GUI settings
        values for DNA concentration, salt concentrations, and dNTP
        concentration used in melting temperature calculations.

        Returns:
            A TMSettings instance configured from GUI settings.
        """
        from amplifyp.settings import TMSettings

        if self._cached_tm_settings is not None:
            return self._cached_tm_settings

        self._cached_tm_settings = TMSettings(
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
            dntp_conc=self._safe_float(
                "tm_dNTP_conc", GLOBAL_TM_SETTINGS.dntp_conc
            ),
        )
        return self._cached_tm_settings

    def calculate_primer_tm(self, primer: "Primer") -> float:
        """Calculate the melting temperature of a primer based on settings.

        Uses the configured TM method (SantaLucia 1998 / Owczarzy 2008
        or Lander / Amplify 4) and current thermodynamic settings to
        compute the melting temperature.

        Args:
            primer: The primer object to calculate Tm for.

        Returns:
            The melting temperature as a float.
        """
        from amplifyp.errors import InsufficientThermodynamicDataError
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
        try:
            return calculate_tm_santalucia_1998_owczarzy_2008(
                primer, tm_settings
            )
        except InsufficientThermodynamicDataError:
            return calculate_tm_lander_amplify4(primer, tm_settings)

    def to_dict(self) -> dict[str, Any]:
        """Convert settings to a dictionary for serialisation.

        Returns:
            A copy of the internal settings dictionary.
        """
        return dict(self._settings)

    def from_dict(self, settings_dict: dict[str, Any]) -> None:
        """Update settings from a dictionary.

        Applies values from the dictionary to matching keys, performing
        type coercion as needed. Also updates GUIColours state for
        colour_deficient and dark_mode settings.

        Args:
            settings_dict: Dictionary of setting key-value pairs to apply.
        """
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
        GUIColours.colour_deficient_mode = bool(
            self._settings.get("colour_deficient", False)
        )
        dark_mode_val = self._settings.get("dark_mode", False)
        if str(dark_mode_val).lower() != "system":
            GUIColours.dark_mode = bool(dark_mode_val)
        self._cached_tm_settings = None

    def _get_config_path(self) -> Path:
        """Get the OS-specific path for user settings configuration.

        Returns:
            Path object pointing to the settings.yaml file location:
            - Windows: %APPDATA%/AmplifyP/settings.yaml
            - macOS: ~/Library/Application Support/AmplifyP/settings.yaml
            - Linux: $XDG_CONFIG_HOME/amplifyp/settings.yaml or
                ~/.config/amplifyp/settings.yaml
        """
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
            storage = getattr(page, "client_storage", None)
            if storage is not None:
                # Web browser path
                settings_dict = {}
                for k in self._settings.keys():
                    storage_key = f"amplifyp.settings.{k}"
                    if storage.contains_key(storage_key):
                        settings_dict[k] = storage.get(storage_key)
                if settings_dict:
                    self.from_dict(settings_dict)
            return

        # Desktop path
        config_path = self._get_config_path()
        if config_path.exists():
            try:
                with open(config_path, encoding="utf-8") as f:
                    data = yaml.safe_load(f)
                if isinstance(data, dict):
                    self.from_dict(data)
            except (OSError, yaml.YAMLError, ValueError) as e:
                logger.exception(
                    "Error loading settings from %s: %s", config_path, e
                )

    def save_to_local(self, page: ft.Page) -> None:
        """Save settings to local storage.

        On web platforms, uses page.client_storage.
        On desktop platforms, writes to OS-specific YAML file.
        """
        if getattr(page, "web", False):
            storage = getattr(page, "client_storage", None)
            if storage is not None:
                # Web browser path
                for k, v in self._settings.items():
                    storage_key = f"amplifyp.settings.{k}"
                    storage.set(storage_key, v)
            return

        # Desktop path
        config_path = self._get_config_path()
        try:
            config_path.parent.mkdir(parents=True, exist_ok=True)
            with open(config_path, "w", encoding="utf-8") as f:
                yaml.safe_dump(self.to_dict(), f, default_flow_style=False)
        except (OSError, yaml.YAMLError, ValueError) as e:
            logger.exception("Error saving settings to %s: %s", config_path, e)
