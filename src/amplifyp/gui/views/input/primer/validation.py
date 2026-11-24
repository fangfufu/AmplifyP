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

"""Validation helpers for DNA primers."""

from __future__ import annotations

from typing import Any

from amplifyp.dna import clean_sequence
from amplifyp.io import (
    count_active_primer_names_and_sequences as _count_names_and_sequences,
)
from amplifyp.io import (
    get_duplicate_primer_indices as _core_get_duplicate_primer_indices,
)
from amplifyp.io import validate_primer, validate_primers

__all__ = [
    "_count_names_and_sequences",
    "get_duplicate_primer_indices",
    "reconcile_primer_states",
    "validate_primer",
    "validate_primers",
]


def get_duplicate_primer_indices(
    primers: list[dict[str, Any]],
) -> set[int]:
    """Find and return indices of duplicate active primers by name/sequence."""
    raw_dup_indices = _core_get_duplicate_primer_indices(primers)
    dup_indices: set[int] = set()
    for idx in raw_dup_indices:
        c = primers[idx].get("container")
        c_idx = (
            c.data
            if (c is not None and getattr(c, "data", None) is not None)
            else idx
        )
        dup_indices.add(c_idx)
    return dup_indices


def reconcile_primer_states(
    ui_primers: list[dict[str, Any]],
    prev_primers: list[dict[str, Any]],
    auto_activate_new: bool = False,
) -> list[dict[str, Any]]:
    """Reconcile UI primer state with the previous central state."""
    primers = []
    for i, p in enumerate(ui_primers):
        prev_p = prev_primers[i] if i < len(prev_primers) else {}
        p_name = str(p.get("name") or "")
        p_seq = str(p.get("seq") or "")
        is_filled = bool(p_name.strip() and clean_sequence(p_seq))
        was_active = prev_p.get("active", False)
        is_active = p.get("active", True)

        show_empty_errors = bool(p.get("show_empty_errors", False))

        if auto_activate_new and is_filled and not was_active:
            name_err, seq_err = validate_primer(
                p_name, p_seq, show_empty_errors=False
            )
            if name_err is None and seq_err is None:
                is_active = True
                show_empty_errors = False

        primers.append(
            {
                "name": p_name,
                "seq": p_seq,
                "active": is_active,
                "show_empty_errors": show_empty_errors,
                "name_touched": prev_p.get("name_touched", False),
                "seq_touched": prev_p.get("seq_touched", False),
            }
        )
    return primers
