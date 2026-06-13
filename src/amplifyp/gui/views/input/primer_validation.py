# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Validation logic for DNA primers."""

from typing import Any

from amplifyp.gui.util import clean_sequence


def validate_primer(name: str, seq: str) -> str | None:
    """Validate a single primer sequence and name.

    Returns an error message if the primer is invalid, otherwise None.
    """
    if not seq:
        return None
    try:
        from amplifyp.dna import Primer

        Primer(sequence=seq, name=name)
        return None
    except ValueError as ex:
        return str(ex)


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
        name_val = p.get("name", "")
        seq_val = p.get("seq", "")
        seq_err = validate_primer(name_val, seq_val)
        name_err = None

        n_lower = str(name_val).strip().lower()
        s_lower = clean_sequence(str(seq_val)).lower()

        if not seq_err and s_lower and seqs_count.get(s_lower, 0) > 1:
            seq_err = "Duplicate primer sequence"
        if n_lower and names_count.get(n_lower, 0) > 1:
            name_err = "Duplicate primer name"

        errors.append({"name": name_err, "seq": seq_err})
    return errors
