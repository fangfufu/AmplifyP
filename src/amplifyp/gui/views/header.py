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
        on_clear_all: Callable[[ft.ControlEvent], Any] | None = None,
        on_switch_designer: Callable[[ft.ControlEvent], None] | None = None,
        on_switch_designer_2d: Callable[[ft.ControlEvent], None] | None = None,
    ) -> None:
        """Initialise the AppHeader navigation component."""
        super().__init__(
            spacing=8, horizontal_alignment=ft.CrossAxisAlignment.START
        )
        self.settings = settings

        self.input_button = ft.FilledButton(
            "Input",
            icon=ft.Icons.INPUT,
            on_click=on_switch_input,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="Input",
        )
        self.input_button.tooltip = "Input"

        self.pcr_button = ft.FilledButton(
            "PCR",
            ref=pcr_button_ref,
            on_click=on_pcr_click,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            disabled=True,
            icon=ft.Icons.ANALYTICS,
            tooltip="PCR",
        )
        self.pcr_button.tooltip = "PCR"

        self.dimers_button = ft.FilledButton(
            "Primer Dimers",
            ref=dimers_button_ref,
            on_click=on_dimers_click,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            disabled=True,
            icon=ft.Icons.COMPARE_ARROWS,
            tooltip="Primer Dimers",
        )
        self.dimers_button.tooltip = "Primer Dimers"

        self.designer_button = ft.FilledButton(
            "Designer 1D",
            icon=ft.Icons.TUNE,
            on_click=on_switch_designer,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="1D Primer Designer",
        )
        self.designer_button.tooltip = "1D Primer Designer"

        self.designer_2d_button = ft.FilledButton(
            "Designer 2D",
            icon=ft.Icons.GRID_ON,
            on_click=on_switch_designer_2d,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="2D Primer Designer",
        )
        self.designer_2d_button.tooltip = "2D Primer Designer"

        self.settings_button = ft.FilledButton(
            "Settings",
            icon=ft.Icons.SETTINGS,
            on_click=on_switch_settings,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="Settings",
        )
        self.settings_button.tooltip = "Settings"

        self.about_button = ft.FilledButton(
            "About",
            icon=ft.Icons.INFO,
            on_click=on_switch_about,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            tooltip="About",
        )
        self.about_button.tooltip = "About"

        self.nav_buttons: list[ft.FilledButton] = [
            self.input_button,
            self.pcr_button,
            self.dimers_button,
            self.designer_button,
            self.designer_2d_button,
            self.settings_button,
            self.about_button,
        ]
        self._active_button: ft.FilledButton = self.input_button
        self._view_to_button: dict[Any, ft.FilledButton] = {}
        self._name_to_button: dict[str, ft.FilledButton] = {
            "input": self.input_button,
            "pcr": self.pcr_button,
            "dimer": self.dimers_button,
            "dimers": self.dimers_button,
            "primer_dimers": self.dimers_button,
            "primer dimers": self.dimers_button,
            "designer": self.designer_button,
            "designer_1d": self.designer_button,
            "designer 1d": self.designer_button,
            "designer_2d": self.designer_2d_button,
            "designer 2d": self.designer_2d_button,
            "settings": self.settings_button,
            "about": self.about_button,
        }
        self.set_active_button(self.input_button)

        self.save_btn_control = ft.OutlinedButton(
            "Save all",
            icon=ft.Icons.SAVE,
            tooltip="Save all",
            on_click=on_save,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        self.clear_btn_control = ft.OutlinedButton(
            "Clear all",
            icon=ft.Icons.DELETE,
            tooltip="Clear all",
            on_click=on_clear_all,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        self.load_btn_control = ft.OutlinedButton(
            "Load all",
            icon=ft.Icons.UPLOAD_FILE,
            tooltip="Load all",
            on_click=on_load,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        self.header_divider = ft.Container(
            width=1,
            height=20,
            bgcolor=GUIColours.OUTLINE,
        )

        self.appbar_actions = [
            self.input_button,
            self.pcr_button,
            self.dimers_button,
            self.designer_button,
            self.designer_2d_button,
            self.settings_button,
            self.about_button,
            self.clear_btn_control,
            self.save_btn_control,
            self.load_btn_control,
        ]

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
                                self.input_button,
                                self.pcr_button,
                                self.dimers_button,
                                self.designer_button,
                                self.designer_2d_button,
                                self.settings_button,
                                self.about_button,
                            ],
                            spacing=10,
                            tight=True,
                        ),
                        self.header_divider,
                        ft.Row(
                            [
                                self.clear_btn_control,
                                self.save_btn_control,
                                self.load_btn_control,
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

    @property
    def active_style(self) -> ft.ButtonStyle:
        """Create button style for the active view navigation button."""
        return ft.ButtonStyle(
            bgcolor={
                ft.ControlState.DEFAULT: GUIColours.NAV_ACTIVE_BG,
                ft.ControlState.DISABLED: GUIColours.DISABLED_BG,
            },
            color={
                ft.ControlState.DEFAULT: GUIColours.NAV_ACTIVE_FG,
                ft.ControlState.DISABLED: GUIColours.DISABLED_FG,
            },
        )

    @property
    def inactive_style(self) -> ft.ButtonStyle:
        """Create button style for inactive view navigation buttons."""
        return ft.ButtonStyle(
            bgcolor={
                ft.ControlState.DEFAULT: GUIColours.NAV_INACTIVE_BG,
                ft.ControlState.DISABLED: GUIColours.DISABLED_BG,
            },
            color={
                ft.ControlState.DEFAULT: GUIColours.NAV_INACTIVE_FG,
                ft.ControlState.DISABLED: GUIColours.DISABLED_FG,
            },
        )

    @property
    def active_button(self) -> ft.FilledButton:
        """Get the currently highlighted navigation button."""
        return self._active_button

    def set_active_button(self, target_button: ft.FilledButton) -> None:
        """Set the active navigation button and apply highlighting."""
        self._active_button = target_button
        active_style = self.active_style
        inactive_style = self.inactive_style
        for btn in self.nav_buttons:
            btn.style = active_style if btn is target_button else inactive_style

    def register_view_buttons(
        self, view_map: dict[Any, ft.FilledButton]
    ) -> None:
        """Register mapping from view controls to navigation buttons."""
        self._view_to_button.update(view_map)

    def set_active_view(self, view: Any) -> None:
        """Update active view highlight based on view control or identifier."""
        target_btn: ft.FilledButton | None = None
        if view in self.nav_buttons:
            target_btn = view
        elif view in self._view_to_button:
            target_btn = self._view_to_button[view]
        elif isinstance(view, str):
            target_btn = self._name_to_button.get(view.lower().strip())

        if target_btn is not None:
            self.set_active_button(target_btn)

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
