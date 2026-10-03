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

from collections.abc import Callable
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.utils.gui_helpers import show_error_dialog


class NavigationManager:
    """Manages view switching, header setup, and resize event routing."""

    def __init__(self, controller: Any) -> None:
        """Initialise NavigationManager with a reference to the controller."""
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
            on_switch_designer_2d=lambda e: self.switch_view(
                e, ctrl.designer_2d_view
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

        # Register views for active button highlighting
        ctrl.header.register_view_buttons(
            {
                ctrl.input_view: ctrl.header.input_button,
                ctrl.pcr_view: ctrl.header.pcr_button,
                ctrl.dimers_view: ctrl.header.dimers_button,
                ctrl.designer_view: ctrl.header.designer_button,
                ctrl.designer_2d_view: ctrl.header.designer_2d_button,
                ctrl.settings_view: ctrl.header.settings_button,
                ctrl.about_view: ctrl.header.about_button,
            }
        )
        ctrl.header.set_active_view(ctrl.input_view)

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

        if hasattr(ctrl, "header") and ctrl.header is not None:
            ctrl.header.set_active_view(view)

        if view == ctrl.input_view:
            ctrl.page.on_resize = ctrl.input_view._handle_resize
        elif view == ctrl.pcr_view:
            ctrl.page.on_resize = ctrl.pcr_view._handle_resize  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
        else:
            ctrl.page.on_resize = None

        ctrl.page.update()

    def _validate_primers_for_action(
        self, action_name: str
    ) -> tuple[str, str] | None:
        """Validate active primers before navigating to an analysis view.

        Args:
            action_name: The action description for the error message
                (e.g., 'PCR' or 'primer dimer analysis').

        Returns:
            A tuple of (title, message) if validation fails, or None if valid.
        """
        ctrl = self.controller
        active_primers = [
            p for p in ctrl.input_data.primers if p.get("active", False)
        ]
        if not active_primers:
            return (
                "Primers Required",
                f"Please select at least one primer before running "
                f"{action_name}.",
            )

        primers_valid = True
        if ctrl.input_view is not None and hasattr(
            ctrl.input_view, "primer_input"
        ):
            res = ctrl.input_view.primer_input.validate_for_run()
            if isinstance(res, bool):
                primers_valid = res
            elif hasattr(ctrl.input_view.primer_input, "validation_errors"):
                val_errs = ctrl.input_view.primer_input.validation_errors
                for idx, p in enumerate(ctrl.input_data.primers):
                    if p.get("active", False) and idx < len(val_errs):
                        err = val_errs[idx]
                        if err and (err.get("name") or err.get("seq")):
                            primers_valid = False
                            break

        if not primers_valid:
            return (
                "Invalid Primers",
                "One or more selected primers are invalid, have "
                "empty names/sequences, or have duplicate "
                "names/sequences.",
            )
        return None

    def _validate_for_pcr(self) -> tuple[str, str] | None:
        """Validate state before navigating to the PCR view.

        Returns:
            A tuple of (title, message) if validation fails, or None if valid.
        """
        ctrl = self.controller
        has_template = bool(ctrl.input_data.template.strip())
        if not has_template:
            return (
                "Template Required",
                "Please enter a DNA template in the Input view before "
                "running PCR.",
            )
        return self._validate_primers_for_action("PCR")

    def _validate_for_dimers(self) -> tuple[str, str] | None:
        """Validate state before navigating to the Primer Dimers view.

        Returns:
            A tuple of (title, message) if validation fails, or None if valid.
        """
        return self._validate_primers_for_action("primer dimer analysis")

    def _handle_analysis_click(
        self,
        e: ft.ControlEvent,
        validator: Callable[[], tuple[str, str] | None],
        target_view: ft.Control,
        run_action: Callable[[], bool],
    ) -> None:
        """Validate input, switch view, and run analysis.

        Args:
            e: The event that triggered the click.
            validator: Validation callback returning (title, message) or None.
            target_view: The view to transition to.
            run_action: The analysis callback to execute.
        """
        ctrl = self.controller
        ctrl.update_pcr_button_state(sync=True)
        validation_error = validator()
        if validation_error:
            if (
                ctrl.input_view is not None
                and ctrl.view_container.content != ctrl.input_view
            ):
                self.switch_view(e, ctrl.input_view)
            title, message = validation_error
            show_error_dialog(ctrl.page, title, message)
            return

        if ctrl.input_view is not None and hasattr(
            ctrl.input_view, "primer_input"
        ):
            ctrl.input_view.primer_input.reset_validation_mode()

        self.switch_view(e, target_view)
        if not run_action():
            self.switch_view(e, ctrl.input_view)

    def on_pcr_click(self, e: ft.ControlEvent) -> None:
        """Handle PCR click: validate input, switch view, then run PCR.

        The view is switched only if validation passes, ensuring canvas
        shapes render while the PCR view is active.
        """
        self._handle_analysis_click(
            e=e,
            validator=self._validate_for_pcr,
            target_view=self.controller.pcr_view,
            run_action=self.controller.pcr_view.run_pcr,
        )

    def on_dimers_click(self, e: ft.ControlEvent) -> None:
        """Handle dimers click: validate, switch view, and run analysis."""
        self._handle_analysis_click(
            e=e,
            validator=self._validate_for_dimers,
            target_view=self.controller.dimers_view,
            run_action=self.controller.dimers_view.run_analysis,
        )
