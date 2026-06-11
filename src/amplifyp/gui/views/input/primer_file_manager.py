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

            if "\t" in line:
                delimiter = "\t"
            else:
                delimiter = ","

            parts = line.split(delimiter, 1)
            if len(parts) != 2:
                continue

            name = parts[0].strip()
            seq = parts[1].strip()

            if not name and not seq:
                continue

            parsed_primers.append(
                {
                    "name": name,
                    "seq": seq,
                    "active": True,
                }
            )
        return parsed_primers

    def _serialize_primers_to_csv(self, primers: list[dict[str, Any]]) -> str:
        """Serialize primers list to a CSV string."""
        output = io.StringIO()
        writer = csv.writer(output)
        for p in primers:
            writer.writerow([p.get("name", ""), p.get("seq", "")])

        csv_content = output.getvalue()
        output.close()
        return csv_content

    async def load_primers_click(self, e: ft.ControlEvent) -> None:
        """Open file picker to load primers from CSV/TSV file."""
        try:
            files = await ft.FilePicker().pick_files(
                dialog_title="Load Primers",
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
        """Save primers to a CSV file."""
        primers_to_save = [
            p
            for p in self.input_data.primers
            if str(p.get("name", "")).strip() or str(p.get("seq", "")).strip()
        ]
        if not primers_to_save:
            self.show_snackbar("No primers to save.")
            return

        csv_content = self._serialize_primers_to_csv(primers_to_save)

        try:
            file_path = await ft.FilePicker().save_file(
                dialog_title="Save Primers",
                file_name="primers.csv",
                allowed_extensions=["csv"],
                file_type=ft.FilePickerFileType.CUSTOM,
                src_bytes=csv_content.encode("utf-8"),
            )
            if self.app_page.web:
                self.show_snackbar("Primers ready for download!")
            else:
                if file_path is None:
                    return
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(csv_content)
                self.show_snackbar(f"Saved {len(primers_to_save)} primer(s).")
        except Exception as ex:
            self.show_snackbar(f"Error saving file: {ex}")
