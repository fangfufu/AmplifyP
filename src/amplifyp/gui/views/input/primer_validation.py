# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Validation logic for DNA primers."""

from typing import Any

from amplifyp.gui.util import clean_sequence


def validate_primer(
    name: str, seq: str, show_empty_errors: bool = False
) -> tuple[str | None, str | None]:
    """Validate a single primer sequence and name.

    Returns a tuple of (name_error, seq_error).
    """
    name_err = None
    seq_err = None

    if not name.strip():
        if show_empty_errors:
            name_err = "Name cannot be empty"

    if not seq.strip():
        if show_empty_errors:
            seq_err = "Sequence cannot be empty"
    else:
        try:
            from amplifyp.dna import Primer

            Primer(sequence=seq, name=name)
        except ValueError as ex:
            seq_err = str(ex)

    return name_err, seq_err


def validate_primers(
    primers: list[dict[str, Any]],
) -> list[dict[str, str | None]]:
    """Validate a list of primers, detecting format and duplicate errors."""
    names_count: dict[str, int] = {}
    seqs_count: dict[str, int] = {}
    for p in primers:
        n_lower = str(p.get("name", "")).strip().lower()
        s_lower = clean_sequence(str(p.get("seq", ""))).lower()
        if n_lower:
            names_count[n_lower] = names_count.get(n_lower, 0) + 1
        if s_lower:
            seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1

    errors = []
    for p in primers:
        name_val = str(p.get("name", "")).strip()
        seq_val = clean_sequence(str(p.get("seq", "")))

        show_empty_errors = bool(
            p.get("name_touched", False) and p.get("seq_touched", False)
        )
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
