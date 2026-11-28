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

from __future__ import annotations

import inspect
import logging
from collections.abc import ItemsView, Iterator, KeysView
from pathlib import Path
from typing import TYPE_CHECKING, Any

import flet as ft
import yaml

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.os_paths import get_settings_yaml_path
from amplifyp.io import safe_dump_yaml, safe_load_yaml
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
    from amplifyp.melting import TmMethod, TmRange
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
        self._cached_replication_settings: ReplicationSettings | None = None
        self._cached_primer_dimer_settings: PrimerDimerSettings | None = None
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
            "sequence_separator": "Dash ('-')",
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
            "tm_colour_scheme": "None",
            "designer_2d_colour_scheme": "Blue-Orange",
            "designer_2d_show_rev_fwd": False,
            "log_level_amplifyp": "INFO",
            "log_level_flet": "INFO",
            "log_console_enabled": True,
            "log_file_enabled": True,
            "log_file_path": "(Default)",
            "log_rotation_enabled": True,
            "log_rotation_max_bytes": 5242880,
            "version_checking_frequency": "Once per Day",
            "last_version_check_timestamp": 0.0,
            "auto_reload_on_startup": True,
            "template_fixed_width": False,
            "template_bases_per_line": "Auto",
            "primer_info_panel_position": "bottom",
            "primer_info_panel_fixed_height": False,
            "primer_info_panel_height": 260,
            "auto_activate_new_valid_primer": False,
            "designer_1d_show_tm": True,
            "designer_1d_show_pct_at": False,
            "pcr_amplicon_ranking": "Position, then length",
            "pcr_relative_q_score": True,
            "pcr_primer_label_orientation": "Horizontal",
            "pcr_divider_behaviour": "autoresize_3_4",
            "pcr_map_height_ratio": 0.375,
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

        # Initialise lengthwise match and run weights
        from amplifyp.settings import DEFAULT_MATCH_WEIGHTS, DEFAULT_RUN_WEIGHTS

        self._settings["match_weight_default"] = str(
            DEFAULT_MATCH_WEIGHTS.default_weight
        )
        for idx, val in DEFAULT_MATCH_WEIGHTS.overrides.items():
            self._settings[f"match_weight_{idx}"] = str(val)

        self._settings["run_weight_default"] = str(
            DEFAULT_RUN_WEIGHTS.default_weight
        )
        for idx, val in DEFAULT_RUN_WEIGHTS.overrides.items():
            self._settings[f"run_weight_{idx}"] = str(val)

        if settings_dict is not None:
            self.from_dict(settings_dict)
            for k, v in settings_dict.items():
                if (
                    k not in self._settings
                    and k
                    not in (
                        "amplify4_compatibility_mode",
                        "dimer_sequence_separator",
                    )
                    and not k.startswith(("match_weight_", "run_weight_"))
                ):
                    self._settings[k] = v

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

    def _invalidate_caches(self) -> None:
        """Invalidate cached core settings objects after mutations."""
        self._cached_tm_settings = None
        self._cached_replication_settings = None
        self._cached_primer_dimer_settings = None

    def __delitem__(self, key: str) -> None:
        """Delete a setting key."""
        del self._settings[key]
        self._invalidate_caches()

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
        self._invalidate_caches()
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

    def get_replication_settings(self) -> ReplicationSettings:
        """Get ReplicationSettings from the central settings.

        Constructs a ReplicationSettings object using the current GUI
        settings values for primability/stability cutoffs, base-pair
        weights, and Amplify4 compatibility mode.

        Note:
            The returned instance is cached and shared until this
            ``GUISettings`` object is modified. Callers must treat it as
            read-only.

        Returns:
            A ReplicationSettings instance configured from GUI settings.
        """
        if self._cached_replication_settings is not None:
            return self._cached_replication_settings

        from amplifyp.settings import ReplicationSettings

        normalised = dict(self._settings)
        normalised["amplify4_compatibility_mode"] = normalised.get(
            "amp4_compat", False
        )
        self._cached_replication_settings = ReplicationSettings.from_dict(
            normalised
        )
        return self._cached_replication_settings

    def get_primer_dimer_settings(self) -> PrimerDimerSettings:
        """Get PrimerDimerSettings from the central settings.

        Constructs a PrimerDimerSettings object using the current GUI
        settings values for minimum overlap, threshold, and base-pair
        weights.

        Note:
            The returned instance is cached and shared until this
            ``GUISettings`` object is modified. Callers must treat it as
            read-only.

        Returns:
            A PrimerDimerSettings instance configured from GUI settings.
        """
        if self._cached_primer_dimer_settings is not None:
            return self._cached_primer_dimer_settings

        from amplifyp.settings import PrimerDimerSettings

        self._cached_primer_dimer_settings = PrimerDimerSettings.from_dict(
            self._settings
        )
        return self._cached_primer_dimer_settings

    def get_tm_settings(self) -> TMSettings:
        """Get TMSettings from the central settings.

        Constructs a TMSettings object using the current GUI settings
        values for DNA concentration, salt concentrations, and dNTP
        concentration used in melting temperature calculations.

        Note:
            The returned instance is cached and shared until this
            ``GUISettings`` object is modified. Callers must treat it as
            read-only.

        Returns:
            A TMSettings instance configured from GUI settings.
        """
        from amplifyp.settings import TMSettings

        if self._cached_tm_settings is not None:
            return self._cached_tm_settings

        self._cached_tm_settings = TMSettings.from_dict(self._settings)
        return self._cached_tm_settings

    def _get_tm_method(self) -> TmMethod:
        """Resolve the configured Tm calculation method."""
        from amplifyp.melting import TmMethod

        tm_method_raw = str(
            self.get("tm_method", TmMethod.SANTALUCIA_1998_OWCZARZY_2008)
        )
        try:
            return TmMethod(tm_method_raw)
        except ValueError:
            return TmMethod.SANTALUCIA_1998_OWCZARZY_2008

    def calculate_primer_tm(self, primer: Primer) -> float:
        """Calculate the melting temperature of a primer based on settings.

        Uses the configured TM method (SantaLucia 1998 / Owczarzy 2008
        or Lander / Amplify 4) and current thermodynamic settings to
        compute the melting temperature.

        Args:
            primer: The primer object to calculate Tm for.

        Returns:
            The melting temperature as a float.
        """
        from amplifyp.melting import calculate_tm

        return calculate_tm(
            primer,
            settings=self.get_tm_settings(),
            method=self._get_tm_method(),
            fallback_to_amplify4=True,
        )

    def calculate_primer_tm_range(self, primer: Primer) -> TmRange:
        """Calculate min/mean/max melting temperatures for a primer.

        Args:
            primer: The primer object to calculate Tm range for.

        Returns:
            A ``TmRange`` instance across IUPAC variants.
        """
        from amplifyp.melting import calculate_tm_range

        return calculate_tm_range(
            primer,
            settings=self.get_tm_settings(),
            method=self._get_tm_method(),
            fallback_to_amplify4=True,
        )

    def resolve_primer_tm_display(
        self, primer: Primer
    ) -> tuple[float, TmRange | None]:
        """Return headline display Tm (max for degenerate primers) and range.

        Args:
            primer: The primer object to evaluate.

        Returns:
            Tuple of ``(display_tm, tm_range)``, where ``tm_range`` is
            populated only when multiple IUPAC variants were evaluated.
        """
        tm = self.calculate_primer_tm(primer)
        if primer.redundant_base_count > 0 and len(primer) >= 2:
            tm_range = self.calculate_primer_tm_range(primer)
            if tm_range.variant_count > 1:
                return tm_range.max_tm, tm_range
        return tm, None

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
        try:
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
                        if isinstance(v, str):
                            self._settings[k] = v.lower() in (
                                "true",
                                "1",
                                "yes",
                            )
                        else:
                            self._settings[k] = bool(v)
                    elif isinstance(self._settings[k], int):
                        try:
                            self._settings[k] = int(v)
                        except (ValueError, TypeError, OverflowError):
                            pass
                    elif isinstance(self._settings[k], float):
                        try:
                            self._settings[k] = float(v)
                        except (ValueError, TypeError, OverflowError):
                            pass
                    else:
                        self._settings[k] = str(v)
            for prefix in ("match_weight_", "run_weight_"):
                valid_entries: dict[str, str] = {}
                for incoming_k, incoming_v in settings_dict.items():
                    if not incoming_k.startswith(prefix):
                        continue
                    suffix = incoming_k[len(prefix) :]
                    if suffix == "default":
                        valid_entries[incoming_k] = str(incoming_v)
                    elif suffix.isascii() and suffix.isdigit():
                        canonical_k = f"{prefix}{int(suffix)}"
                        if (
                            canonical_k not in valid_entries
                            or incoming_k == canonical_k
                        ):
                            valid_entries[canonical_k] = str(incoming_v)
                if valid_entries:
                    for existing_k in [
                        k for k in self._settings if k.startswith(prefix)
                    ]:
                        del self._settings[existing_k]
                    self._settings.update(valid_entries)
            if (
                "amplify4_compatibility_mode" in settings_dict
                and "amp4_compat" not in settings_dict
            ):
                from amplifyp.settings import coerce_bool

                self._settings["amp4_compat"] = coerce_bool(
                    settings_dict["amplify4_compatibility_mode"],
                    bool(self._settings.get("amp4_compat", False)),
                )
            if (
                "dimer_sequence_separator" in settings_dict
                and "sequence_separator" not in settings_dict
            ):
                legacy_sep = settings_dict["dimer_sequence_separator"]
                if legacy_sep is not None:
                    self._settings["sequence_separator"] = str(legacy_sep)
            GUIColours.colour_deficient_mode = bool(
                self._settings.get("colour_deficient", False)
            )
            dark_mode_val = self._settings.get("dark_mode", False)
            if str(dark_mode_val).lower() != "system":
                GUIColours.dark_mode = bool(dark_mode_val)
        finally:
            self._invalidate_caches()

    def reset_to_defaults(self) -> None:
        """Reset all settings to their factory default values."""
        defaults = GUISettings()
        self._settings = dict(defaults._settings)
        self._invalidate_caches()
        GUIColours.colour_deficient_mode = bool(
            self._settings.get("colour_deficient", False)
        )
        dark_mode_val = self._settings.get("dark_mode", False)
        if str(dark_mode_val).lower() == "system":
            GUIColours.dark_mode = False
        else:
            GUIColours.dark_mode = bool(dark_mode_val)

    def _get_config_path(self) -> Path:
        """Get the OS-specific path for user settings configuration.

        Returns:
            Path object pointing to the settings.yaml file location:
            - Windows: %APPDATA%/AmplifyP/settings.yaml
            - macOS: ~/Library/Application Support/AmplifyP/settings.yaml
            - Linux: $XDG_CONFIG_HOME/amplifyp/settings.yaml or
                ~/.config/amplifyp/settings.yaml
        """
        return get_settings_yaml_path()

    async def load_from_local_async(self, page: ft.Page) -> None:
        """Asynchronously load settings from web SharedPreferences."""
        prefs = ft.SharedPreferences()
        settings_dict: dict[str, Any] = {}
        storage_prefix = "amplifyp.settings."
        if inspect.iscoroutinefunction(getattr(prefs, "get_keys", None)):
            all_keys = await prefs.get_keys(storage_prefix)
            for full_key in all_keys:
                if full_key.startswith(storage_prefix):
                    k = full_key[len(storage_prefix) :]
                    settings_dict[k] = await prefs.get(full_key)
        else:
            for k in self._settings:
                storage_key = f"{storage_prefix}{k}"
                if await prefs.contains_key(storage_key):
                    settings_dict[k] = await prefs.get(storage_key)
            for prefix in ("match_weight_", "run_weight_"):
                idx = 0
                while await prefs.contains_key(
                    f"{storage_prefix}{prefix}{idx}"
                ):
                    k = f"{prefix}{idx}"
                    if k not in settings_dict:
                        settings_dict[k] = await prefs.get(
                            f"{storage_prefix}{k}"
                        )
                    idx += 1
        if settings_dict:
            self.from_dict(settings_dict)

    def load_from_local(self, page: ft.Page) -> None:
        """Load settings from local storage.

        On web platforms, uses ft.SharedPreferences via page.run_task.
        On desktop platforms, reads from OS-specific YAML file.
        """
        if getattr(page, "web", False):
            if hasattr(page, "run_task"):
                page.run_task(self.load_from_local_async, page)
            return

        # Desktop path
        config_path = self._get_config_path()
        if config_path.exists():
            try:
                with open(config_path, encoding="utf-8") as f:
                    data = safe_load_yaml(f)
                if isinstance(data, dict):
                    self.from_dict(data)
            except (OSError, yaml.YAMLError, ValueError) as e:
                logger.exception(
                    "Error loading settings from %s: %s", config_path, e
                )

    def save_to_local(
        self, page: ft.Page, keys: tuple[str, ...] | list[str] | None = None
    ) -> Any:
        """Save settings to local storage.

        On web platforms, uses ft.SharedPreferences via page.run_task.
        On desktop platforms, writes to OS-specific YAML file.
        """
        if getattr(page, "web", False):
            if hasattr(page, "run_task"):
                target_keys = (
                    keys if keys is not None else tuple(self._settings.keys())
                )
                snapshot = {
                    k: self._settings[k]
                    for k in target_keys
                    if k in self._settings
                }
                prev_task = getattr(self, "_web_save_task", None)

                async def _save_web() -> None:
                    if prev_task is not None and (
                        not hasattr(prev_task, "done") or not prev_task.done()
                    ):
                        try:
                            await prev_task
                        except Exception as e:
                            logger.debug("Previous web save failed: %s", e)
                    prefs = ft.SharedPreferences()
                    storage_prefix = "amplifyp.settings."
                    for k, v in snapshot.items():
                        storage_key = f"{storage_prefix}{k}"
                        await prefs.set(storage_key, v)
                    if (
                        keys is None
                        and inspect.iscoroutinefunction(
                            getattr(prefs, "get_keys", None)
                        )
                        and inspect.iscoroutinefunction(
                            getattr(prefs, "remove", None)
                        )
                    ):
                        existing_keys = await prefs.get_keys(storage_prefix)
                        for full_key in existing_keys:
                            k = full_key[len(storage_prefix) :]
                            if k not in snapshot and (
                                k.startswith(("match_weight_", "run_weight_"))
                                or k == "dimer_sequence_separator"
                            ):
                                await prefs.remove(full_key)

                task = page.run_task(_save_web)
                self._web_save_task = task
                return task
            return None

        # Desktop path
        config_path = self._get_config_path()
        try:
            config_path.parent.mkdir(parents=True, exist_ok=True)
            with open(config_path, "w", encoding="utf-8") as f:
                safe_dump_yaml(self.to_dict(), f, default_flow_style=False)
        except (OSError, yaml.YAMLError, ValueError) as e:
            logger.exception("Error saving settings to %s: %s", config_path, e)
