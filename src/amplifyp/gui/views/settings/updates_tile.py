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

"""UpdatesTile component for Flet settings view."""

import asyncio
import time
from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings


class UpdatesTile(ft.ExpansionTile):  # type: ignore[misc]
    """Expansion tile for updates and version checking settings."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        settings_map: dict[str, Any],
        on_change_handler: Callable[[ft.ControlEvent | None], None],
        header_size: int,
        on_update_found: Callable[[str], None] | None = None,
    ) -> None:
        """Initialise the UpdatesTile."""
        self._page = page
        self.settings = settings
        self.settings_map = settings_map
        self.on_change_handler = on_change_handler
        self.on_update_found = on_update_found

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

        super().__init__(
            title=ft.Text(
                "Updates Settings",
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

    def _on_manual_check_click(self, e: ft.ControlEvent) -> None:
        """Handle manual update check trigger."""
        self._page.run_task(self.perform_manual_check)

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
        self._page.update()

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
            self.settings.save_to_local(self._page)

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

        self._page.update()

    def update_ui(self) -> None:
        """Sync component with settings state."""
        self.set_version_checking_frequency.value = self.settings.get(
            "version_checking_frequency", "Once per Month"
        )
