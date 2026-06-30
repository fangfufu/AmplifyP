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

"""Validation logic for DNA primers."""

from typing import Any

from amplifyp.dna import Primer
from amplifyp.gui.utils.sequence import clean_sequence


def validate_primer(
    name: str, seq: str, show_empty_errors: bool = False
) -> tuple[str | None, str | None]:
    """Validate a single primer sequence and name.

    Checks for empty fields and validates the DNA sequence using the
    Primer class. Returns error messages for invalid fields.

    Args:
        name: The primer name to validate.
        seq: The primer DNA sequence to validate.
        show_empty_errors: If True, reports empty name/sequence as errors.

    Returns:
        A tuple of (name_error, seq_error) where each element is an
        error message string or None if valid.
    """
    name_err = None
    seq_err = None

    if not name.strip() and show_empty_errors:
        name_err = "Name cannot be empty"

    if not seq.strip() and show_empty_errors:
        seq_err = "Sequence cannot be empty"
    else:
        try:
            Primer(sequence=seq, name=name)
        except ValueError as ex:
            seq_err = str(ex)

    return name_err, seq_err


def _count_names_and_sequences(
    primers: list[dict[str, Any]],
) -> tuple[dict[str, int], dict[str, int]]:
    """Count occurrences of primer names and sequences."""
    names_count: dict[str, int] = {}
    seqs_count: dict[str, int] = {}
    for p in primers:
        n_lower = str(p.get("name", "")).strip().lower()
        s_lower = clean_sequence(str(p.get("seq", ""))).lower()
        if n_lower:
            names_count[n_lower] = names_count.get(n_lower, 0) + 1
        if s_lower:
            seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1
    return names_count, seqs_count


def get_duplicate_primer_indices(primers: list[dict[str, Any]]) -> set[int]:
    """Find and return indices of duplicate primers by name/sequence.

    Args:
        primers: List of primer dicts to check for duplicates.

    Returns:
        Set of indices corresponding to primers with duplicate names
        or sequences.
    """
    names_count, seqs_count = _count_names_and_sequences(primers)

    dup_indices = set()
    for idx, p in enumerate(primers):
        n_lower = str(p.get("name", "")).strip().lower()
        s_lower = clean_sequence(str(p.get("seq", ""))).lower()
        if (n_lower and names_count.get(n_lower, 0) > 1) or (
            s_lower and seqs_count.get(s_lower, 0) > 1
        ):
            c = p.get("container")
            c_idx = c.data if (c is not None and hasattr(c, "data")) else idx
            dup_indices.add(c_idx)
    return dup_indices


def reconcile_primer_states(
    ui_primers: list[dict[str, Any]], prev_primers: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Reconcile UI primer state with the previous central state.

    Handles auto-activation transitions, touched flags, and empty
    error visibility.

    Args:
        ui_primers: List of primer dicts extracted from UI.
        prev_primers: List of primer dicts representing the previous
            central state.

    Returns:
        A list of reconciled primer dicts.
    """
    primers = []
    for i, p in enumerate(ui_primers):
        prev_p = prev_primers[i] if i < len(prev_primers) else {}
        is_filled = bool(p["name"].strip() and p["seq"].strip())
        was_empty = (
            not prev_p.get("name", "").strip()
            or not prev_p.get("seq", "").strip()
        )
        was_active = prev_p.get("active", False)
        is_active = p["active"]

        show_empty_errors = prev_p.get("show_empty_errors", False)
        if is_active and not is_filled:
            is_active = False
            show_empty_errors = True
        elif not is_active:
            show_empty_errors = False

        if is_filled and was_empty and not was_active:
            is_active = True
            show_empty_errors = False

        primers.append(
            {
                "name": p["name"],
                "seq": p["seq"],
                "active": is_active,
                "show_empty_errors": show_empty_errors,
                "name_touched": prev_p.get("name_touched", False),
                "seq_touched": prev_p.get("seq_touched", False),
            }
        )
    return primers


def validate_primers(
    primers: list[dict[str, Any]],
) -> list[dict[str, str | None]]:
    """Validate a list of primers, detecting format and duplicate errors.

    Validates each primer's name and sequence, then checks for
    duplicate names and sequences across the entire list.

    Args:
        primers: List of primer dicts with 'name', 'seq', 'name_touched',
            and 'seq_touched' keys.

    Returns:
        List of error dicts, one per primer, with 'name' and 'seq' keys
        containing error message strings or None.
    """
    names_count, seqs_count = _count_names_and_sequences(primers)

    errors = []
    for p in primers:
        name_val = str(p.get("name", "")).strip()
        seq_val = clean_sequence(str(p.get("seq", "")))

        show_empty_errors = bool(p.get("show_empty_errors", False))
        name_err, seq_err = validate_primer(
            name_val, seq_val, show_empty_errors=show_empty_errors
        )

        n_lower = name_val.lower()
        s_lower = seq_val.lower()

        if not seq_err and s_lower and seqs_count.get(s_lower, 0) > 1:
            seq_err = "Duplicate primer sequence"
        if not name_err and n_lower and names_count.get(n_lower, 0) > 1:
            name_err = "Duplicate primer name"

        errors.append({"name": name_err, "seq": seq_err})
    return errors
