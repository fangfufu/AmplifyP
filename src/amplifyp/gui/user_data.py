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

"""Centralised GUI state and user data orchestrator."""

from typing import Any

from amplifyp.gui.utils.data_helpers import clean_sequence, format_sequence


class GUIInput:
    """Encapsulates user input including DNA template and primers list."""

    def __init__(
        self,
        template: str = "",
        template_circular: bool = False,
        primers: list[dict[str, Any]] | None = None,
    ) -> None:
        """Initialise GUIInput."""
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

    def to_dict(self) -> dict[str, Any]:
        """Convert input data to a dictionary for serialisation."""
        return {
            "template": format_sequence(self.template),
            "template_circular": self.template_circular,
            "primers": [
                {
                    "name": str(p.get("name", "")),
                    "seq": format_sequence(str(p.get("seq", ""))),
                    "active": bool(p.get("active", True)),
                }
                for p in self.primers
            ],
        }

    def from_dict(self, state_dict: dict[str, Any]) -> None:
        """Update input data from a dictionary."""
        if "template" in state_dict:
            self.template = clean_sequence(state_dict["template"])
        if "template_circular" in state_dict:
            self.template_circular = bool(state_dict["template_circular"])
        if "primers" in state_dict:
            self._load_primers(state_dict["primers"])
