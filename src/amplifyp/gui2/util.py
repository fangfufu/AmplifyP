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

from amplifyp.gui2.utils.data_helpers import (
    _read_file,
    _resolve_font_family,
    _write_file,
    clean_sequence,
    format_sequence,
    pick_and_read_file,
    save_and_write_file,
    serialise_state,
)
from amplifyp.gui2.utils.gui_helpers import (
    BorderedCheckbox,
    Debouncer,
    NotificationHelper,
    initialise_score_fields,
    show_error_dialog,
)
from amplifyp.gui2.utils.system import get_full_sha, get_git_sha, get_version


def serialise_state_wrapper(state: dict[str, object]) -> str:
    """Serialise state dict to YAML string, handling multiline strings."""
    return serialise_state(state)


__all__ = [
    "BorderedCheckbox",
    "Debouncer",
    "NotificationHelper",
    "_read_file",
    "_resolve_font_family",
    "_write_file",
    "clean_sequence",
    "format_sequence",
    "get_full_sha",
    "get_git_sha",
    "get_version",
    "initialise_score_fields",
    "pick_and_read_file",
    "save_and_write_file",
    "serialise_state",
    "serialise_state_wrapper",
    "show_error_dialog",
]
