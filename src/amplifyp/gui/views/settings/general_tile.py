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

"""GeneralTile component for Flet settings view."""

from __future__ import annotations

import asyncio
import logging
import time
from collections.abc import Callable
from typing import Any

import flet as ft
import yaml

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.util import serialise_state
from amplifyp.gui.utils.io import pick_and_read_file, save_and_write_file
from amplifyp.gui.utils.ui import BorderedCheckbox, NotificationHelper

logger = logging.getLogger(__name__)


class GeneralTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for General settings."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.ControlEvent | None], None],
        header_size: int,
        font_size_default: int,
        sync_to_state_callback: Callable[[], None],
        update_ui_callback: Callable[[], None],
        on_update_found: Callable[[str], None] | None = None,
    ) -> None:
        """Initialise the GeneralTile.

        Args:
            page: The Flet page object.
            settings: The settings object.
            settings_map: A dictionary mapping setting keys to UI
                components for population and retrieval.
            on_change_handler: The handler to call when a setting changes.
            header_size: The size of the header text.
            font_size_default: The default font size for body text.
            sync_to_state_callback: Callback to sync UI to settings state.
            update_ui_callback: Callback to update UI controls.
            on_update_found: Callback when an update is found.
        """
        self.app_page = page
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler
        self.sync_to_state_callback = sync_to_state_callback
        self.update_ui_callback = update_ui_callback
        self.on_update_found = on_update_found
        self.notification_helper = NotificationHelper(page)
        self.filepicker_open = False

        # --- Appearance settings ---
        self.set_font_family = ft.Dropdown(
            label="Font Family",
            options=[
                ft.dropdown.Option("Roboto Mono"),
                ft.dropdown.Option("Courier New"),
                ft.dropdown.Option("Consolas"),
                ft.dropdown.Option("monospace"),
            ],
            width=500,
            on_select=self.on_change_handler,
            border_color=GUIColours.OUTLINE,
        )

        self.set_colour_scheme = ft.Dropdown(
            label="Colour Scheme",
            options=[
                ft.dropdown.Option("Light"),
                ft.dropdown.Option("Dark"),
                ft.dropdown.Option("Light (Colour Deficient Friendly)"),
                ft.dropdown.Option("Dark (Colour Deficient Friendly)"),
                ft.dropdown.Option("System"),
                ft.dropdown.Option("System (Colour Deficient Friendly)"),
            ],
            width=500,
            on_select=self._on_colour_scheme_change,
            border_color=GUIColours.OUTLINE,
        )

        self._dummy_colour_deficient = ft.Checkbox(visible=False)

        # --- General settings ---
        self.auto_reload_checkbox = BorderedCheckbox(
            label="Automatically reload last template and primers on startup",
            on_change=self.on_change_handler,
        )
        self.settings_map["auto_reload_on_startup"] = self.auto_reload_checkbox

        # --- Updates settings ---
        self.set_version_checking_frequency = ft.Dropdown(
            label="Version Checking Frequency",
            options=[
                ft.dropdown.Option("At Startup"),
                ft.dropdown.Option("Once per Day"),
                ft.dropdown.Option("Once per Week"),
                ft.dropdown.Option("Once per Month"),
                ft.dropdown.Option("Disabled"),
            ],
            width=500,
            on_select=self.on_change_handler,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            border_color=GUIColours.OUTLINE,
        )
        self.check_button = ft.OutlinedButton(
            "Check for Updates",
            icon=ft.Icons.REFRESH,
            on_click=self._on_manual_check_click,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )
        self.status_text = ft.Text(
            "",
            size=13,
            italic=True,
            color=GUIColours.MUTED_GREY,
        )
        self.settings_map["version_checking_frequency"] = (
            self.set_version_checking_frequency
        )

        # --- Backup settings ---
        self.save_button = ft.FilledButton(
            "Save Settings",
            icon=ft.Icons.SAVE,
            on_click=self._save_settings_async,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="Save settings to a YAML file",
        )
        self.load_button = ft.OutlinedButton(
            "Load Settings",
            icon=ft.Icons.UPLOAD_FILE,
            on_click=self._load_settings_async,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="Load settings from a YAML file",
        )

        self.settings_map["font_family"] = self.set_font_family
        self.settings_map["colour_deficient"] = self._dummy_colour_deficient

        super().__init__(
            title=ft.Text(
                "General Settings",
                weight=ft.FontWeight.BOLD,
                size=header_size,
            ),
            expanded_cross_axis_alignment=ft.CrossAxisAlignment.STRETCH,
            controls=[
                ft.Container(
                    content=ft.Row(
                        [
                            ft.Container(
                                content=ft.Column(
                                    [
                                        ft.Text(
                                            "Appearance",
                                            weight=ft.FontWeight.BOLD,
                                            size=font_size_default + 2,
                                        ),
                                        self.set_font_family,
                                        self.set_colour_scheme,
                                        self._dummy_colour_deficient,
                                        ft.Divider(),
                                        ft.Text(
                                            "Autosave",
                                            weight=ft.FontWeight.BOLD,
                                            size=font_size_default + 2,
                                        ),
                                        self.auto_reload_checkbox,
                                        ft.Divider(),
                                        ft.Text(
                                            "Updates",
                                            weight=ft.FontWeight.BOLD,
                                            size=font_size_default + 2,
                                        ),
                                        ft.Text(
                                            "Configure how AmplifyP checks "
                                            "for new releases on GitHub.",
                                            size=13,
                                            color=GUIColours.TEXT_ON_SURFACE,
                                        ),
                                        self.set_version_checking_frequency,
                                        ft.Row(
                                            [
                                                self.check_button,
                                                self.status_text,
                                            ],
                                            spacing=15,
                                            vertical_alignment=ft.CrossAxisAlignment.CENTER,
                                        ),
                                        ft.Divider(),
                                        ft.Text(
                                            "Backup Settings",
                                            weight=ft.FontWeight.BOLD,
                                            size=font_size_default + 2,
                                        ),
                                        ft.Text(
                                            "Export or import all application "
                                            "settings exposed in this view.",
                                            size=font_size_default,
                                        ),
                                        ft.Row(
                                            [
                                                self.save_button,
                                                self.load_button,
                                            ],
                                            spacing=15,
                                            alignment=ft.MainAxisAlignment.CENTER,
                                        ),
                                    ],
                                    spacing=15,
                                    horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
                                ),
                                width=500,
                            ),
                        ],
                        alignment=ft.MainAxisAlignment.CENTER,
                    ),
                    padding=ft.Padding(0, 20, 0, 10),
                )
            ],
        )

    @property
    def set_auto_reload_on_startup(self) -> BorderedCheckbox:
        """Get the auto reload checkbox."""
        return self.auto_reload_checkbox

    def _on_manual_check_click(self, e: ft.ControlEvent) -> None:
        """Handle manual update check trigger."""
        self.app_page.run_task(self.perform_manual_check)

    async def perform_manual_check(self) -> None:
        """Asynchronously run update check and update UI."""
        from amplifyp import __version__ as current_version
        from amplifyp.gui.utils.version_check import (
            fetch_latest_release_version,
            is_newer_version,
        )

        self.check_button.disabled = True
        self.status_text.value = "Checking for updates..."
        self.status_text.color = GUIColours.MUTED_GREY
        self.app_page.update()

        loop = asyncio.get_running_loop()
        latest_tag = await loop.run_in_executor(
            None, fetch_latest_release_version
        )

        self.check_button.disabled = False

        if latest_tag is None:
            self.status_text.value = (
                "Could not check for updates. "
                "Please check your network connection."
            )
            self.status_text.color = GUIColours.ERROR_RED
        else:
            self.settings["last_version_check_timestamp"] = float(time.time())
            self.settings.save_to_local(self.app_page)

            if is_newer_version(latest_tag, current_version):
                self.status_text.value = (
                    f"New version {latest_tag} is available!"
                )
                self.status_text.color = GUIColours.SUCCESS_GREEN
                if self.on_update_found:
                    self.on_update_found(latest_tag)
            else:
                self.status_text.value = "AmplifyP is up to date."
                self.status_text.color = GUIColours.SUCCESS_GREEN

        self.app_page.update()

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
            if self.on_change_handler is not None:
                self.on_change_handler(None)

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

    def update_ui(self) -> None:
        """Sync component with settings state."""
        self.set_version_checking_frequency.value = self.settings.get(
            "version_checking_frequency", "Once per Month"
        )

    @property
    def set_colour_deficient(self) -> ft.Checkbox:
        """Get colour deficient checkbox for backwards compatibility."""
        return self._dummy_colour_deficient

    def _on_colour_scheme_change(self, e: ft.ControlEvent) -> None:
        """Handle colour scheme dropdown change.

        Syncs the selected colour scheme to settings and triggers the
        change handler.

        Args:
            e: The Flet control event triggered by the dropdown change.
        """
        self.sync_colour_scheme_to_settings()
        self.on_change_handler(e)

    def update_colour_scheme_dropdown(self) -> None:
        """Update the colour scheme dropdown value based on settings.

        Reads the current dark_mode and colour_deficient settings and
        sets the dropdown to the matching display option.
        """
        dark = self.settings.get("dark_mode", False)
        deficient = bool(self.settings.get("colour_deficient", False))
        if str(dark).lower() == "system":
            if deficient:
                self.set_colour_scheme.value = (
                    "System (Colour Deficient Friendly)"
                )
            else:
                self.set_colour_scheme.value = "System"
        elif bool(dark):
            if deficient:
                self.set_colour_scheme.value = (
                    "Dark (Colour Deficient Friendly)"
                )
            else:
                self.set_colour_scheme.value = "Dark"
        else:
            if deficient:
                self.set_colour_scheme.value = (
                    "Light (Colour Deficient Friendly)"
                )
            else:
                self.set_colour_scheme.value = "Light"

    def sync_colour_scheme_to_settings(self) -> None:
        """Sync the colour scheme dropdown selection back to settings.

        Parses the dropdown value and updates dark_mode and
        colour_deficient settings accordingly, including the hidden
        checkbox state.
        """
        val = self.set_colour_scheme.value
        if val == "Dark":
            self.settings["dark_mode"] = True
            self.settings["colour_deficient"] = False
            self._dummy_colour_deficient.value = False
        elif val == "Dark (Colour Deficient Friendly)":
            self.settings["dark_mode"] = True
            self.settings["colour_deficient"] = True
            self._dummy_colour_deficient.value = True
        elif val == "Light (Colour Deficient Friendly)":
            self.settings["dark_mode"] = False
            self.settings["colour_deficient"] = True
            self._dummy_colour_deficient.value = True
        elif val == "System":
            self.settings["dark_mode"] = "system"
            self.settings["colour_deficient"] = False
            self._dummy_colour_deficient.value = False
        elif val == "System (Colour Deficient Friendly)":
            self.settings["dark_mode"] = "system"
            self.settings["colour_deficient"] = True
            self._dummy_colour_deficient.value = True
        else:  # Light
            self.settings["dark_mode"] = False
            self.settings["colour_deficient"] = False
            self._dummy_colour_deficient.value = False
