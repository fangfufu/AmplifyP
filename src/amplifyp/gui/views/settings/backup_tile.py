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

"""BackupTile component for Flet settings view."""

import logging
from collections.abc import Callable

import flet as ft
import yaml

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.util import (
    NotificationHelper,
    pick_and_read_file,
    save_and_write_file,
    serialise_state,
)

logger = logging.getLogger(__name__)


class BackupTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for Backup settings (saving and loading configuration)."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        sync_to_state_callback: Callable[[], None],
        update_ui_callback: Callable[[], None],
        on_change_callback: Callable[[ft.ControlEvent | None], None] | None,
        header_size: int,
        font_size_default: int,
    ) -> None:
        """Initialise the BackupTile."""
        self.app_page = page
        self.settings = settings
        self.sync_to_state_callback = sync_to_state_callback
        self.update_ui_callback = update_ui_callback
        self.on_change_callback = on_change_callback
        self.notification_helper = NotificationHelper(page)
        self.filepicker_open = False

        self.save_button = ft.FilledButton(
            "Save Settings",
            icon=ft.Icons.SAVE,
            on_click=self._save_settings_async,
            tooltip="Save settings to a YAML file",
        )

        self.load_button = ft.OutlinedButton(
            "Load Settings",
            icon=ft.Icons.UPLOAD_FILE,
            on_click=self._load_settings_async,
            tooltip="Load settings from a YAML file",
        )

        super().__init__(
            title=ft.Text(
                "Backup",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Column(
                        [
                            ft.Text(
                                "Export or import all application settings "
                                "exposed in this view.",
                                size=font_size_default,
                            ),
                            ft.Container(height=5),
                            ft.Row(
                                [
                                    self.save_button,
                                    self.load_button,
                                ],
                                spacing=15,
                                alignment=ft.MainAxisAlignment.START,
                            ),
                        ],
                        spacing=10,
                        horizontal_alignment=ft.CrossAxisAlignment.START,
                    ),
                    padding=ft.Padding(20, 20, 20, 20),
                )
            ],
        )

    async def _save_settings_async(self, e: ft.ControlEvent) -> None:
        """Save settings to YAML configuration file."""
        if self.filepicker_open:
            return
        self.filepicker_open = True
        try:
            self.sync_to_state_callback()
            settings_dict: dict[str, object] = {
                "settings": self.settings.to_dict(),
            }
            yaml_str = serialise_state(settings_dict)

            await save_and_write_file(
                page=self.app_page,
                dialog_title="Save Settings",
                file_name="amplify_settings.yaml",
                allowed_extensions=["yaml", "yml"],
                content=yaml_str,
                show_notification=self.notification_helper.show_message,
                success_message_desktop="Settings saved successfully!",
                success_message_web="Settings ready for download!",
            )
        except (OSError, ValueError) as ex:
            self.notification_helper.show_message(
                f"Error saving settings: {ex}"
            )
        finally:
            self.filepicker_open = False

    async def _load_settings_async(self, e: ft.ControlEvent) -> None:
        """Load settings from YAML configuration file."""
        if self.filepicker_open:
            return
        self.filepicker_open = True
        try:
            content = await pick_and_read_file(
                page=self.app_page,
                dialog_title="Load Settings",
                allowed_extensions=["yaml", "yml"],
                show_notification=self.notification_helper.show_message,
            )
            if content is None:
                return

            parsed_state = yaml.safe_load(content)
            if not isinstance(parsed_state, dict):
                self.notification_helper.show_message(
                    "Error: Invalid settings file format."
                )
                return

            settings_dict = (
                parsed_state.get("settings")
                if "settings" in parsed_state
                else parsed_state
            )

            if not isinstance(settings_dict, dict):
                self.notification_helper.show_message(
                    "Error: Invalid settings data format."
                )
                return

            self.settings.from_dict(settings_dict)
            self.settings.save_to_local(self.app_page)

            # Sync values back to the UI controls
            self.update_ui_callback()

            # Apply loaded changes (theme, notify views, update PCR
            # button, etc.)
            if self.on_change_callback:
                self.on_change_callback(None)

            self.notification_helper.show_message(
                "Settings loaded successfully!"
            )
        except (OSError, ValueError, yaml.YAMLError) as ex:
            logger.exception("Error loading settings:")
            self.notification_helper.show_message(
                f"Error loading settings: {ex}"
            )
        finally:
            self.filepicker_open = False
