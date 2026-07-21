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

"""GUI utility function re-exports."""

import yaml

from amplifyp.gui.utils.data_helpers import (
    _read_file,
    _resolve_font_family,
    _write_file,
    clean_sequence,
    create_overlapped_sequence_view,
    format_sequence,
    pick_and_read_file,
    save_and_write_file,
)
from amplifyp.gui.utils.gui_helpers import (
    BorderedCheckbox,
    Debouncer,
    NotificationHelper,
    initialise_score_fields,
    show_error_dialog,
)
from amplifyp.gui.utils.system import get_full_sha, get_git_sha, get_version


def serialise_state(state: dict[str, object]) -> str:
    """Serialise state dict to YAML string, handling multiline strings.

    Uses a custom YAML dumper that represents multiline strings using
    the '|' block style for better readability.

    Args:
        state: The state dictionary to serialise.

    Returns:
        A YAML string representation of the state.
    """

    def multiline_presenter(dumper: yaml.Dumper, data: str) -> yaml.ScalarNode:
        """Represent multiline strings using the '|' block style in YAML."""
        if "\n" in data:
            return dumper.represent_scalar(
                "tag:yaml.org,2002:str", data, style="|"
            )
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    class _StateDumper(yaml.Dumper):
        """Custom YAML dumper using multiline representer for strings."""

    _StateDumper.add_representer(str, multiline_presenter)
    return yaml.dump(state, Dumper=_StateDumper, sort_keys=False)


__all__ = [
    "BorderedCheckbox",
    "Debouncer",
    "NotificationHelper",
    "_read_file",
    "_resolve_font_family",
    "_write_file",
    "clean_sequence",
    "create_overlapped_sequence_view",
    "format_sequence",
    "get_full_sha",
    "get_git_sha",
    "get_version",
    "initialise_score_fields",
    "pick_and_read_file",
    "save_and_write_file",
    "serialise_state",
    "show_error_dialog",
]
