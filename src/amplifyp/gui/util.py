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

"""Utility functions for sequence handling and state serialization."""

import yaml


def clean_sequence(seq: str) -> str:
    """Clean sequence of escaped and standard whitespaces."""
    if not seq:
        return ""
    clean = str(seq).replace("\\n", "").replace("\\t", "").replace("\\r", "")
    return "".join(clean.split()).upper()


def format_sequence(seq: str, wrap_length: int = 80) -> str:
    """Format sequence into lines of specified length."""
    clean = clean_sequence(seq)
    return "\n".join(
        [clean[i : i + wrap_length] for i in range(0, len(clean), wrap_length)]
    )


def serialize_state(state: dict[str, object]) -> str:
    """Serialize state dict to YAML string, handling multiline strings."""

    def multiline_presenter(dumper: yaml.Dumper, data: str) -> yaml.ScalarNode:
        if "\n" in data:
            return dumper.represent_scalar(
                "tag:yaml.org,2002:str", data, style="|"
            )
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    yaml.add_representer(str, multiline_presenter)
    return yaml.dump(state, sort_keys=False)
