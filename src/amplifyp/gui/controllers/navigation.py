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

"""Navigation controller for view routing and header orchestration."""

from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours


class NavigationManager:
    """Manages view switching, header setup, and resize event routing."""

    def __init__(self, controller: Any) -> None:
        """Initialize NavigationManager with a reference to the controller."""
        self.controller = controller

    def setup_navigation_controls(self) -> None:
        """Configure navigation controls for the main application window.

        Creates and sets up the AppBar buttons and the visible top header
        buttons (Input, PCR, Primer Dimers, Settings, Save, Load).
        """
        from amplifyp.gui.views.header import AppHeader

        ctrl = self.controller
        ctrl.header = AppHeader(
            settings=ctrl.settings,
            on_switch_input=lambda e: self.switch_view(e, ctrl.input_view),
            on_switch_settings=lambda e: self.switch_view(
                e, ctrl.settings_view
            ),
            on_switch_about=lambda e: self.switch_view(e, ctrl.about_view),
            on_pcr_click=self.on_pcr_click,
            on_dimers_click=self.on_dimers_click,
            on_switch_designer=lambda e: self.switch_view(
                e, ctrl.designer_view
            ),
            on_save=ctrl.save_state,
            on_load=ctrl.load_state,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            on_clear_all=ctrl.clear_all,
            pcr_button_ref=ctrl.pcr_button_ref,
            dimers_button_ref=ctrl.dimers_button_ref,
        )

        # Store aliases for direct accesses
        ctrl.save_btn_control = ctrl.header.save_btn_control
        ctrl.clear_btn_control = ctrl.header.clear_btn_control
        ctrl.load_btn_control = ctrl.header.load_btn_control
        ctrl.header_divider = ctrl.header.header_divider

        # Configure page appbar
        ctrl.page.appbar = ft.AppBar(
            visible=False,
            actions=ctrl.header.appbar_actions,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        )

        ctrl.header_container = ft.Container(
            content=ctrl.header,
            padding=ft.Padding(16, 8, 16, 8),
            bgcolor=GUIColours.SURFACE,
        )

        ctrl.page.add(
            ft.Divider(height=1, thickness=1),
            ctrl.view_container,
        )
        ctrl.page.controls.insert(0, ctrl.header_container)
        ctrl.page.on_resize = ctrl.input_view._handle_resize
        ctrl.page.update()
        # After the first update, platform_brightness is populated.
        # Re-apply theme and refresh views to resolve dynamic colours correctly.
        ctrl.apply_theme()
        ctrl.input_view.update_ui()
        ctrl.page.update()

    def switch_view(self, _e: ft.Event[ft.Control], view: ft.Control) -> None:
        """Switch the main view container to display a different view.

        Updates the container content and configures resize handlers
        appropriate for the target view.

        Args:
            _e: The event that triggered the view switch (unused).
            view: The Flet control to display as the new view.
        """
        ctrl = self.controller
        if view == ctrl.input_view and ctrl.input_view_dirty:
            ctrl.input_view.update_ui()
            ctrl.input_view_dirty = False

        ctrl.view_container.content = view
        is_input = view == ctrl.input_view
        ctrl.save_btn_control.visible = is_input
        ctrl.clear_btn_control.visible = is_input
        ctrl.load_btn_control.visible = is_input
        ctrl.header_divider.visible = is_input

        if view == ctrl.input_view:
            ctrl.page.on_resize = ctrl.input_view._handle_resize
        elif view == ctrl.pcr_view:
            ctrl.page.on_resize = ctrl.pcr_view._handle_resize  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        else:
            ctrl.page.on_resize = None

        ctrl.page.update()

    def on_pcr_click(self, e: ft.ControlEvent) -> None:
        """Handle PCR click: switch view then run PCR.

        The view is switched first so the diagram canvas renders
        while the PCR view is the active content.  This avoids
        Flet's diff algorithm marking canvas shapes as 'already
        sent' before the view becomes visible.
        """
        ctrl = self.controller
        ctrl.update_pcr_button_state(sync=True)
        self.switch_view(e, ctrl.pcr_view)
        if not ctrl.pcr_view.run_pcr():
            self.switch_view(e, ctrl.input_view)

    def on_dimers_click(self, e: ft.ControlEvent) -> None:
        """Handle dimers click: switch view then run analysis."""
        ctrl = self.controller
        ctrl.update_pcr_button_state(sync=True)
        self.switch_view(e, ctrl.dimers_view)
        if not ctrl.dimers_view.run_analysis():
            self.switch_view(e, ctrl.input_view)
