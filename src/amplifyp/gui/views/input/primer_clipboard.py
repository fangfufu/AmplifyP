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

"""Helper function for parsing primer clipboard text."""

from typing import Any

from amplifyp.gui.utils.sequence import clean_sequence


def parse_primer_clipboard_text(text: str) -> list[dict[str, Any]]:
    """Parse pasted text containing primers.

    Args:
        text: The raw clipboard or pasted text.

    Returns:
        A list of dicts with keys 'name', 'seq', and 'active'.
    """
    parsed = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            name = parts[0].strip()
            seq = clean_sequence(parts[1])
        elif len(parts) == 1:
            subparts = line.split(",")
            if len(subparts) >= 2:
                name = subparts[0].strip()
                seq = clean_sequence(subparts[1])
            else:
                val = line.strip()
                cleaned = clean_sequence(val)
                is_seq = False
                if cleaned:
                    is_seq = all(
                        c in "ACGTRYSWKMBDHVNacgtryswkmbdhvn" for c in cleaned
                    )
                if is_seq:
                    name = ""
                    seq = cleaned
                else:
                    name = val
                    seq = ""
        else:
            continue
        parsed.append({"name": name, "seq": seq, "active": False})
    return parsed
