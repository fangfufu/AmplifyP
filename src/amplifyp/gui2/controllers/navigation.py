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


class NavigationManager:
    """Manages view switching, header setup, and resize event routing."""

    def __init__(self, controller: Any) -> None:
        """Initialize NavigationManager with a reference to the controller."""
        self.controller = controller

    def setup_navigation_controls(self) -> None:
        """Configure navigation controls for the main application window."""
        from amplifyp.gui2.views.header import AppHeader

        ctrl = self.controller
        ctrl.header = AppHeader(
            settings=ctrl.settings,
            on_switch_input=lambda: self.switch_view(ctrl.input_view),
            on_switch_settings=lambda: self.switch_view(ctrl.settings_view),
            on_switch_about=lambda: self.switch_view(ctrl.about_view),
            on_pcr_click=self.on_pcr_click,
            on_dimers_click=self.on_dimers_click,
            on_switch_designer=lambda: self.switch_view(ctrl.designer_view),
            on_switch_designer_2d=lambda: self.switch_view(
                ctrl.designer_2d_view
            ),
            on_save=lambda: self._save_state(),
            on_load=lambda: self._load_state(),
            on_clear_all=lambda: self._clear_all(),
            pcr_button_ref=ctrl.pcr_button_ref,
            dimers_button_ref=ctrl.dimers_button_ref,
        )

        ctrl.save_btn_control = ctrl.header.save_btn_control
        ctrl.clear_btn_control = ctrl.header.clear_btn_control
        ctrl.load_btn_control = ctrl.header.load_btn_control
        ctrl.header_divider = ctrl.header.header_divider

        ctrl.header_container = ctrl.header

    def _save_state(self) -> None:
        """Save state via controller."""
        self.controller.save_state()

    def _load_state(self) -> None:
        """Load state via controller."""
        self.controller.load_state()

    def _clear_all(self) -> None:
        """Clear all via controller."""
        self.controller.clear_all()

    def switch_view(self, view: Any) -> None:
        """Switch the main view container to display a different view.

        Args:
            view: The widget to display as the new view.
        """
        ctrl = self.controller
        if view == ctrl.input_view and ctrl.input_view_dirty:
            ctrl.input_view.update_ui()
            ctrl.input_view_dirty = False

        ctrl.view_container.setCurrentWidget(view)
        is_input = view == ctrl.input_view
        ctrl.save_btn_control.setVisible(is_input)
        ctrl.clear_btn_control.setVisible(is_input)
        ctrl.load_btn_control.setVisible(is_input)
        ctrl.header_divider.setVisible(is_input)

    def on_pcr_click(self) -> None:
        """Handle PCR click: switch view then run PCR."""
        ctrl = self.controller
        ctrl.update_pcr_button_state(sync=True)
        self.switch_view(ctrl.pcr_view)
        if not ctrl.pcr_view.run_pcr():
            self.switch_view(ctrl.input_view)

    def on_dimers_click(self) -> None:
        """Handle dimers click: switch view then run analysis."""
        ctrl = self.controller
        ctrl.update_pcr_button_state(sync=True)
        self.switch_view(ctrl.dimers_view)
        if not ctrl.dimers_view.run_analysis():
            self.switch_view(ctrl.input_view)
