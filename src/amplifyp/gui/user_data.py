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

"""Centralized GUI state and user data orchestrator."""

from typing import TYPE_CHECKING, Any

from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.util import clean_sequence, format_sequence

if TYPE_CHECKING:
    from amplifyp.settings import PrimerDimerSettings, ReplicationSettings


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
