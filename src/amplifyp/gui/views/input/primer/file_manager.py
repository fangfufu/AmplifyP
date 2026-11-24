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

"""Handles loading and saving of primers using file pickers."""

import csv
from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.user_data import GUIInput
from amplifyp.io import parse_primers_from_text, serialise_primers_to_tsv


class PrimerFileManager:
    """Handles loading and saving of primers using file pickers."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput,
        on_update_ui: Callable[[], None],
        on_change_handler: Callable[[ft.Event[ft.Control] | None], None],
        show_notification: Callable[[str], None],
    ) -> None:
        """Initialise the PrimerFileManager.

        Args:
            page: The Flet page instance for file picker operations.
            input_data: Central GUI input state containing primer list.
            on_update_ui: Callback to refresh the primer UI after loading.
            on_change_handler: Callback to notify of primer state changes.
            show_notification: Callback to display notification messages.
        """
        self.app_page = page
        self.input_data = input_data
        self.on_update_ui = on_update_ui
        self.on_change_handler = on_change_handler
        self.show_notification = show_notification

    def _parse_primers_from_text(self, content: str) -> list[dict[str, Any]]:
        """Parse primers from CSV/TSV content."""
        return parse_primers_from_text(content, csv_reader=csv.reader)

    def _serialise_primers_to_tsv(self, primers: list[dict[str, Any]]) -> str:
        """Serialise primers list to a TSV string."""
        return serialise_primers_to_tsv(primers)

    async def load_primers_click(self, _e: ft.Event | None) -> None:
        """Open file picker to load primers from CSV/TSV file.

        Reads the selected file, parses primers, appends them to the
        current primer list, and updates the UI.

        Args:
            _e: The Flet control event (unused).
        """
        from amplifyp.gui.utils.data_helpers import pick_and_read_file

        content = await pick_and_read_file(
            page=self.app_page,
            dialog_title="Load",
            allowed_extensions=["csv", "tsv", "txt"],
            show_notification=self.show_notification,
        )
        if content is None:
            return

        try:
            parsed = self._parse_primers_from_text(content)
            if len(parsed) > 0:
                # Overwrite if list currently contains only a single empty
                # primer
                if (
                    len(self.input_data.primers) == 1
                    and not str(
                        self.input_data.primers[0].get("name", "")
                    ).strip()
                    and not str(
                        self.input_data.primers[0].get("seq", "")
                    ).strip()
                ):
                    self.input_data.primers.clear()

                for p in parsed:
                    self.input_data.primers.append(p)

                self.on_update_ui()
                self.on_change_handler(None)
                self.show_notification(f"Loaded {len(parsed)} primer(s).")
            else:
                self.show_notification("No valid primers found in file.")
        except (OSError, ValueError, csv.Error) as ex:
            self.show_notification(f"Error parsing primers: {ex}")

    async def save_primers_click(self, _e: ft.Event | None) -> None:
        """Save primers to a TSV file.

        Opens a file picker for the user to choose a save location,
        then writes the primer list as tab-separated values.

        Args:
            _e: The Flet control event (unused).
        """
        primers_to_save = [
            p
            for p in self.input_data.primers
            if str(p.get("name", "")).strip() or str(p.get("seq", "")).strip()
        ]
        if not primers_to_save:
            self.show_notification("No primers to save.")
            return

        tsv_content = self._serialise_primers_to_tsv(primers_to_save)

        from amplifyp.gui.utils.data_helpers import save_and_write_file

        await save_and_write_file(
            page=self.app_page,
            dialog_title="Save",
            file_name="primers.tsv",
            allowed_extensions=["tsv"],
            content=tsv_content,
            show_notification=self.show_notification,
            success_message_desktop=f"Saved {len(primers_to_save)} primer(s).",
            success_message_web="Primers ready for download!",
        )
