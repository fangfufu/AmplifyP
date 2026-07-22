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

"""Application header navigation control."""

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.utils.system import get_version


class AppHeader(ft.Column):  # type: ignore[misc]
    """Application header component with navigation and state buttons."""

    def __init__(
        self,
        settings: GUISettings,
        on_switch_input: Callable[[ft.ControlEvent], None],
        on_switch_settings: Callable[[ft.ControlEvent], None],
        on_switch_about: Callable[[ft.ControlEvent], None],
        on_pcr_click: Callable[[ft.ControlEvent], None],
        on_dimers_click: Callable[[ft.ControlEvent], None],
        on_save: Callable[[ft.ControlEvent], Any],
        on_load: Callable[[ft.ControlEvent], Any],
        pcr_button_ref: ft.Ref[ft.FilledButton],
        dimers_button_ref: ft.Ref[ft.FilledButton],
        visible_pcr_button_ref: ft.Ref[ft.FilledButton],
        visible_dimers_button_ref: ft.Ref[ft.FilledButton],
        on_clear_all: Callable[[ft.ControlEvent], Any] | None = None,
        on_switch_designer: Callable[[ft.ControlEvent], None] | None = None,
    ) -> None:
        """Initialise the AppHeader navigation component."""
        super().__init__(
            spacing=8, horizontal_alignment=ft.CrossAxisAlignment.START
        )
        self.settings = settings

        # AppBar buttons (for test compatibility)
        input_button = ft.FilledButton(
            "Input",
            icon=ft.Icons.INPUT,
            on_click=on_switch_input,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="Input",
        )
        input_button.tooltip = "Input"

        pcr_button = ft.FilledButton(
            "PCR",
            ref=pcr_button_ref,
            on_click=on_pcr_click,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            disabled=True,
            icon=ft.Icons.ANALYTICS,
            tooltip="PCR",
        )
        pcr_button.tooltip = "PCR"

        dimers_button = ft.FilledButton(
            "Primer Dimers",
            ref=dimers_button_ref,
            on_click=on_dimers_click,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            disabled=True,
            icon=ft.Icons.COMPARE_ARROWS,
            tooltip="Primer Dimers",
        )
        dimers_button.tooltip = "Primer Dimers"

        designer_button = ft.FilledButton(
            "1D Designer",
            icon=ft.Icons.TUNE,
            on_click=on_switch_designer,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="1D Primer Designer",
        )
        designer_button.tooltip = "1D Primer Designer"

        settings_button = ft.FilledButton(
            "Settings",
            icon=ft.Icons.SETTINGS,
            on_click=on_switch_settings,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="Settings",
        )
        settings_button.tooltip = "Settings"

        about_button = ft.FilledButton(
            "About",
            icon=ft.Icons.INFO,
            on_click=on_switch_about,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="About",
        )
        about_button.tooltip = "About"

        save_btn_control = ft.FilledButton(
            "Save all",
            icon=ft.Icons.SAVE,
            tooltip="Save all",
            on_click=on_save,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        clear_btn_control = ft.FilledButton(
            "Clear all",
            icon=ft.Icons.DELETE,
            tooltip="Clear all",
            on_click=on_clear_all,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        load_btn_control = ft.FilledButton(
            "Load all",
            icon=ft.Icons.UPLOAD_FILE,
            tooltip="Load all",
            on_click=on_load,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        self.appbar_actions = [
            input_button,
            pcr_button,
            dimers_button,
            designer_button,
            settings_button,
            about_button,
            clear_btn_control,
            save_btn_control,
            load_btn_control,
        ]

        # Visible navigation controls (Responsive header)
        visible_input_button = ft.FilledButton(
            "Input",
            icon=ft.Icons.INPUT,
            on_click=on_switch_input,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="Input",
        )
        visible_input_button.tooltip = "Input"

        visible_pcr_button = ft.FilledButton(
            "PCR",
            ref=visible_pcr_button_ref,
            on_click=on_pcr_click,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            disabled=True,
            icon=ft.Icons.ANALYTICS,
            tooltip="PCR",
        )
        visible_pcr_button.tooltip = "PCR"

        visible_dimers_button = ft.FilledButton(
            "Primer Dimers",
            ref=visible_dimers_button_ref,
            on_click=on_dimers_click,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            disabled=True,
            icon=ft.Icons.COMPARE_ARROWS,
            tooltip="Primer Dimers",
        )
        visible_dimers_button.tooltip = "Primer Dimers"

        visible_designer_button = ft.FilledButton(
            "1D Designer",
            icon=ft.Icons.TUNE,
            on_click=on_switch_designer,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="1D Primer Designer",
        )
        visible_designer_button.tooltip = "1D Primer Designer"

        visible_settings_button = ft.FilledButton(
            "Settings",
            icon=ft.Icons.SETTINGS,
            on_click=on_switch_settings,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="Settings",
        )
        visible_settings_button.tooltip = "Settings"

        visible_about_button = ft.FilledButton(
            "About",
            icon=ft.Icons.INFO,
            on_click=on_switch_about,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="About",
        )
        visible_about_button.tooltip = "About"

        self.visible_save_btn_control = ft.FilledButton(
            "Save all",
            icon=ft.Icons.SAVE,
            tooltip="Save all",
            on_click=on_save,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        self.visible_clear_btn_control = ft.FilledButton(
            "Clear all",
            icon=ft.Icons.DELETE,
            tooltip="Clear all",
            on_click=on_clear_all,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        self.visible_load_btn_control = ft.FilledButton(
            "Load all",
            icon=ft.Icons.UPLOAD_FILE,
            tooltip="Load all",
            on_click=on_load,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        self.visible_header_divider = ft.Container(
            width=1,
            height=20,
            bgcolor=GUIColours.OUTLINE,
        )

        app_version = get_version()
        self.version_text = ft.Text(
            app_version,
            size=14,
            color=GUIColours.TEXT_ON_SURFACE,
            opacity=0.5,
            weight=ft.FontWeight.W_400,
            selectable=True,
        )

        self.controls = [
            ft.Row(
                [
                    ft.Image(
                        src="/images/favicon.png",
                        height=32,
                        fit=ft.BoxFit.CONTAIN,
                    ),
                    ft.Text(
                        "AmplifyP",
                        size=20,
                        weight=ft.FontWeight.BOLD,
                    ),
                    ft.Container(width=12),
                    self.version_text,
                ],
                spacing=8,
                tight=True,
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
            ),
            ft.Container(
                content=ft.Row(
                    [
                        ft.Row(
                            [
                                visible_input_button,
                                visible_pcr_button,
                                visible_dimers_button,
                                visible_designer_button,
                                visible_settings_button,
                                visible_about_button,
                            ],
                            spacing=10,
                            tight=True,
                        ),
                        self.visible_header_divider,
                        ft.Row(
                            [
                                self.visible_clear_btn_control,
                                self.visible_save_btn_control,
                                self.visible_load_btn_control,
                            ],
                            spacing=10,
                            tight=True,
                        ),
                    ],
                    spacing=10,
                    tight=True,
                    wrap=True,
                    vertical_alignment=ft.CrossAxisAlignment.CENTER,
                ),
            ),
        ]

    def set_update_available(self, new_version: str) -> None:
        """Update the version text to show that a new version is available."""
        current_version = get_version()
        self.version_text.value = (
            f"{current_version} (Update {new_version} available!)"
        )
        self.version_text.color = GUIColours.UPDATE_AVAILABLE_COLOUR
        self.version_text.opacity = 1.0
        self.version_text.tooltip = "Click to open download page"
        self.version_text.cursor = ft.MouseCursor.CLICK  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        self.version_text.on_click = lambda e: self.page.launch_url(  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            "https://github.com/fangfufu/AmplifyP/releases"
        )
        self.version_text.update()
