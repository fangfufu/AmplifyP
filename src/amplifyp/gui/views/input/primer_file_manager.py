# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""Handles loading and saving of primers using file pickers."""

import csv
import io
from typing import Any

import flet as ft

from amplifyp.gui.user_data import GUIInput


class PrimerFileManager:
    """Handles loading and saving of primers using file pickers."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput,
        on_update_ui: Any,
        on_change_handler: Any,
        show_snackbar: Any,
    ) -> None:
        """Initialize the PrimerFileManager."""
        self.app_page = page
        self.input_data = input_data
        self.on_update_ui = on_update_ui
        self.on_change_handler = on_change_handler
        self.show_snackbar = show_snackbar

    def _parse_primers_from_text(self, content: str) -> list[dict[str, Any]]:
        """Parse primers from CSV/TSV content."""
        parsed_primers = []
        for line in content.strip().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            delimiter = "\t" if "\t" in line else ","
            reader = csv.reader([line], delimiter=delimiter)
            try:
                row = next(reader)
            except (StopIteration, csv.Error):
                continue

            if not row or len(row) < 2:
                continue

            seq = row[0].strip()
            name = row[1].strip()

            if len(row) >= 3:
                extra = row[2].strip()
                if extra:
                    name = f"{name} - {extra}"

            if not name and not seq:
                continue

            parsed_primers.append(
                {
                    "name": name,
                    "seq": seq,
                    "active": False,
                }
            )
        return parsed_primers

    def _serialize_primers_to_tsv(self, primers: list[dict[str, Any]]) -> str:
        """Serialize primers list to a TSV string."""
        output = io.StringIO()
        writer = csv.writer(output, delimiter="\t")
        for p in primers:
            writer.writerow([p.get("seq", ""), p.get("name", "")])

        tsv_content = output.getvalue()
        output.close()
        return tsv_content

    async def load_primers_click(self, e: ft.ControlEvent) -> None:
        """Open file picker to load primers from CSV/TSV file."""
        try:
            files = await ft.FilePicker().pick_files(
                dialog_title="Load",
                allowed_extensions=["csv", "tsv", "txt"],
                file_type=ft.FilePickerFileType.CUSTOM,
                with_data=True,
            )
            if not files:
                return

            file = files[0]
            if file.bytes is not None:
                content = file.bytes.decode("utf-8")
            else:
                if not file.path:
                    self.show_snackbar("Error: Could not read file content.")
                    return
                with open(file.path, encoding="utf-8") as f:
                    content = f.read()

            parsed = self._parse_primers_from_text(content)
            for p in parsed:
                self.input_data.primers.append(p)

            if len(parsed) > 0:
                self.on_update_ui()
                if self.on_change_handler:
                    self.on_change_handler(None)
                self.show_snackbar(f"Loaded {len(parsed)} primer(s).")
            else:
                self.show_snackbar("No valid primers found in file.")

        except Exception as ex:
            self.show_snackbar(f"Error loading file: {ex}")

    async def save_primers_click(self, e: ft.ControlEvent) -> None:
        """Save primers to a TSV file."""
        primers_to_save = [
            p
            for p in self.input_data.primers
            if str(p.get("name", "")).strip() or str(p.get("seq", "")).strip()
        ]
        if not primers_to_save:
            self.show_snackbar("No primers to save.")
            return

        tsv_content = self._serialize_primers_to_tsv(primers_to_save)

        try:
            file_path = await ft.FilePicker().save_file(
                dialog_title="Save",
                file_name="primers.tsv",
                allowed_extensions=["tsv"],
                file_type=ft.FilePickerFileType.CUSTOM,
                src_bytes=tsv_content.encode("utf-8"),
            )
            if self.app_page.web:
                self.show_snackbar("Primers ready for download!")
            else:
                if file_path is None:
                    return
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(tsv_content)
                self.show_snackbar(f"Saved {len(primers_to_save)} primer(s).")
        except Exception as ex:
            self.show_snackbar(f"Error saving file: {ex}")
