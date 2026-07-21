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

"""Clipboard interaction helpers for DNA primers input."""

from __future__ import annotations

import logging
from typing import Any

import flet as ft

from amplifyp.gui.utils.data_helpers import clean_sequence

logger = logging.getLogger(__name__)


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


async def copy_primers_click(primer_input: Any, e: ft.Event | None) -> None:
    """Copy selected or focused primers to clipboard in TSV format."""
    primers = primer_input.input_data.primers
    selected_primers = (
        [
            primers[i]
            for i in sorted(primer_input.selected_indices)
            if 0 <= i < len(primers)
        ]
        if primer_input.selected_indices
        else []
    )

    # Fallback to focused primer if no selected rows
    if not selected_primers and primer_input.focused_primer_index is not None:
        if 0 <= primer_input.focused_primer_index < len(primers):
            selected_primers = [primers[primer_input.focused_primer_index]]

    if not selected_primers:
        primer_input._show_notification(
            "No primers highlighted or focused to copy."
        )
        return

    lines = []
    for p in selected_primers:
        name = str(p.get("name", "")).strip()
        seq = str(p.get("seq", "")).strip()
        lines.append(f"{name}\t{seq}")

    tsv_text = "\n".join(lines)
    await ft.Clipboard().set(tsv_text)
    primer_input._show_notification(
        f"Copied {len(selected_primers)} primer(s) to clipboard."
    )


async def paste_primers_click(primer_input: Any, e: ft.Event | None) -> None:
    """Paste primers from clipboard starting at focused index or end."""
    try:
        clipboard_text = await ft.Clipboard().get()
    except Exception as ex:
        logger.warning("Failed to access clipboard: %s", ex)
        primer_input._show_notification(
            "Unable to read clipboard. Try pasting directly into a "
            "primer field using Ctrl+V."
        )
        return
    if not clipboard_text:
        primer_input._show_notification("Clipboard is empty.")
        return

    parsed = parse_primer_clipboard_text(clipboard_text)
    if not parsed:
        primer_input._show_notification("No valid primers found in clipboard.")
        return

    primers = primer_input.input_data.primers
    valid_selected = {
        i for i in primer_input.selected_indices if 0 <= i < len(primers)
    }
    if valid_selected:
        insert_idx = max(valid_selected) + 1
    else:
        insert_idx = len(primers)

    # Replace single empty row
    if (
        len(primers) == 1
        and not primers[0].get("name")
        and not primers[0].get("seq")
    ):
        primers.clear()
        insert_idx = 0

    for i, new_p in enumerate(parsed):
        new_p["active"] = False
        primers.insert(insert_idx + i, new_p)

    primer_input.selected_indices = set(
        range(insert_idx, insert_idx + len(parsed))
    )
    primer_input.focused_primer_index = insert_idx

    primer_input.sync_to_state(rebuild_if_needed=True, skip_extract=True)
    if primer_input.on_change_handler is not None:
        primer_input.on_change_handler(None)

    primer_input._show_notification(f"Pasted {len(parsed)} primer(s).")
