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

"""GUI controller for orchestrating views, state, and page events."""

import asyncio
import logging
import time
from typing import Any, cast

import flet as ft  # type: ignore[import-not-found, unused-ignore]
import yaml

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import serialise_state
from amplifyp.gui.utils.ui import NotificationHelper
from amplifyp.gui.views import (
    AboutView,
    DimerView,
    InputView,
    PCRView,
    SettingsView,
)

logger = logging.getLogger(__name__)


class GUIController:
    """Manages GUI state, event handlers, views and main orchestration."""

    def __init__(
        self,
        page: ft.Page,
        state_file: str | None = None,
        auto_close: bool = False,
    ) -> None:
        """Initialize the GUIController.

        Args:
            page: The Flet page instance for the application.
            state_file: Optional path to a YAML state file to load on startup.
            auto_close: If True, quit automatically after rendering is complete.
        """
        self.page = page
        self.state_file = state_file
        self.auto_close = auto_close
        self.input_data = GUIInput()
        self.settings = GUISettings()
        self.filepicker_open = False
        self._confirm_dialog = None

        # Refs for buttons (backward compatibility / lookup)
        self.pcr_button_ref = ft.Ref[ft.FilledButton]()
        self.dimers_button_ref = ft.Ref[ft.FilledButton]()
        self.visible_pcr_button_ref = ft.Ref[ft.FilledButton]()
        self.visible_dimers_button_ref = ft.Ref[ft.FilledButton]()

        # Views placeholders
        self.input_view: InputView = cast(InputView, None)
        self.settings_view: SettingsView = cast(SettingsView, None)
        self.about_view: AboutView = cast(AboutView, None)
        self.pcr_view: PCRView = cast(PCRView, None)
        self.dimers_view: DimerView = cast(DimerView, None)
        self.view_container: ft.Container = cast(ft.Container, None)
        self.header_container: ft.Container = cast(ft.Container, None)

        # UI Control placeholders
        self.visible_save_btn_control: ft.FilledButton = cast(
            ft.FilledButton, None
        )
        self.visible_load_btn_control: ft.FilledButton = cast(
            ft.FilledButton, None
        )
        self.visible_header_divider: ft.Container = cast(ft.Container, None)
        self.notification_helper: NotificationHelper = cast(
            NotificationHelper, None
        )
        self.input_view_dirty = False

    def initialise(self) -> None:
        """Configure page setup, window events, views, and custom layout."""
        self._configure_page_and_window()

        self.settings.load_from_local(self.page)
        self.page.on_platform_brightness_change = (
            self.on_platform_brightness_change
        )
        self.apply_theme()

        # Instantiate views
        self.input_view = InputView(
            self.page,
            self.input_data,
            self.settings,
            on_change=lambda e: self.update_pcr_button_state(sync=False),
            on_stop_editing=lambda e: self.update_pcr_button_state(sync=False),
        )
        self.settings_view = SettingsView(
            self.page,
            self.settings,
            on_change=self.on_settings_change,
            on_reset=self.on_settings_change,
            on_update_found=self.on_update_found,
        )
        self.pcr_view = PCRView(self.page, self.input_data, self.settings)
        self.dimers_view = DimerView(self.page, self.input_data, self.settings)
        self.about_view = AboutView(self.page, self.settings)

        self.notification_helper = NotificationHelper(self.page)

        # Load state and handle auto-close asynchronously once page is ready
        if self.state_file:
            self.page.run_task(self._restore_state_and_auto_close_async)

        # Main view container
        self.view_container = ft.Container(
            content=self.input_view, expand=True, padding=10
        )

        # Header controls & routing buttons setup
        self._setup_navigation_controls()

        # Check for updates in the background on startup
        self.page.run_task(self.check_updates_async)

    def _configure_page_and_window(self) -> None:
        """Set up page properties, window styling, and event handlers."""
        self.page.overlay.clear()
        self.page.title = "AmplifyP"
        self.page.vertical_alignment = ft.MainAxisAlignment.START
        self.page.fonts = {"Roboto Mono": "fonts/RobotoMono-Regular.ttf"}
        self.page.padding = 0
        self.page.spacing = 0
        self.page.window.icon = "/images/icon.png"

        # Handle close / reload warnings
        if self.page.web:
            if hasattr(self.page, "run_javascript"):
                self.page.run_javascript(
                    """
                    window.addEventListener('beforeunload', (event) => {
                        event.preventDefault();
                        event.returnValue = '';
                    });
                    """
                )
        else:
            self.page.window.prevent_close = False
            if not self.auto_close:
                self.page.window.prevent_close = True
                self._confirm_dialog = None
                self.page.window.on_event = self.on_window_event

    def _setup_navigation_controls(self) -> None:
        """Configure navigation controls for the main application window.

        Creates and sets up the AppBar buttons and the visible top header
        buttons (Input, PCR, Primer Dimers, Settings, Save, Load).
        """
        from amplifyp.gui.views.header import AppHeader

        self.header = AppHeader(
            settings=self.settings,
            on_switch_input=lambda e: self.switch_view(e, self.input_view),
            on_switch_settings=lambda e: self.switch_view(
                e, self.settings_view
            ),
            on_switch_about=lambda e: self.switch_view(e, self.about_view),
            on_pcr_click=self.on_pcr_click,
            on_dimers_click=self.on_dimers_click,
            on_save=self.save_state,
            on_load=self.load_state,
            pcr_button_ref=self.pcr_button_ref,
            dimers_button_ref=self.dimers_button_ref,
            visible_pcr_button_ref=self.visible_pcr_button_ref,
            visible_dimers_button_ref=self.visible_dimers_button_ref,
        )

        # Store aliases for backward compatibility or direct accesses
        self.visible_save_btn_control = self.header.visible_save_btn_control
        self.visible_load_btn_control = self.header.visible_load_btn_control
        self.visible_header_divider = self.header.visible_header_divider

        # Configure page appbar
        self.page.appbar = ft.AppBar(
            visible=False,
            actions=self.header.appbar_actions,
        )

        self.header_container = ft.Container(
            content=self.header,
            padding=ft.Padding(16, 8, 16, 8),
            bgcolor=GUIColours.SURFACE,
        )

        self.page.add(
            ft.Divider(height=1, thickness=1),
            self.view_container,
        )
        self.page.controls.insert(0, self.header_container)
        self.page.on_resize = self.input_view._handle_resize
        self.page.update()

    def apply_theme(self) -> None:
        """Apply theme settings (light/dark/system mode) to the page."""
        dark_mode_setting = self.settings.get("dark_mode", False)
        is_dark = False
        if str(dark_mode_setting).lower() == "system":
            self.page.theme_mode = ft.ThemeMode.SYSTEM
            self.page.bg_color = None
            is_dark = str(self.page.platform_brightness).lower() == "dark"
        elif bool(dark_mode_setting) and str(dark_mode_setting).lower() not in (
            "false",
            "0",
            "no",
        ):
            self.page.theme_mode = ft.ThemeMode.DARK
            self.page.bg_color = None
            is_dark = True
        else:
            self.page.theme_mode = ft.ThemeMode.LIGHT
            self.page.bg_color = GUIColours.WHITE
            is_dark = False
        GUIColours.dark_mode = is_dark
        if hasattr(self, "header_container") and self.header_container:
            self.header_container.bgcolor = GUIColours.SURFACE

    def on_platform_brightness_change(
        self, e: ft.ControlEvent | None = None
    ) -> None:
        """Handle system brightness shifts."""
        self.apply_theme()
        active_view = self.view_container.content
        if active_view == self.input_view:
            self.input_view.update_ui()
        else:
            self.input_view_dirty = True

        if active_view == self.settings_view:
            self.settings_view.update_ui()

        if active_view == self.pcr_view:
            self.pcr_view.run_pcr(keep_cards=True)
        elif active_view == self.dimers_view:
            self.dimers_view.run_analysis()
        self.page.update()

    def on_pcr_click(self, e: ft.ControlEvent) -> None:
        """Handle PCR click: switch view then run PCR.

        The view is switched first so the diagram canvas renders
        while the PCR view is the active content.  This avoids
        Flet's diff algorithm marking canvas shapes as 'already
        sent' before the view becomes visible.
        """
        self.update_pcr_button_state(sync=True)
        self.switch_view(e, self.pcr_view)
        if not self.pcr_view.run_pcr():
            self.switch_view(e, self.input_view)

    def on_dimers_click(self, e: ft.ControlEvent) -> None:
        """Handle dimers click: switch view then run analysis."""
        self.update_pcr_button_state(sync=True)
        self.switch_view(e, self.dimers_view)
        if not self.dimers_view.run_analysis():
            self.switch_view(e, self.input_view)

    def update_pcr_button_state(self, sync: bool = True) -> None:
        """Enable PCR and dimers buttons only if input is valid."""
        if sync:
            self.input_view.sync_to_state()

        has_template = bool(self.input_data.template.strip())
        active_primers = self.input_data.get_active_primers()
        has_enough_primers = len(active_primers) >= 1

        # Check if any selected (active) primer has validation errors
        # or duplicates
        has_invalid_selected = False
        for idx, p in enumerate(self.input_data.primers):
            if p.get("active", False) and idx < len(
                self.input_view.primer_input.validation_errors
            ):
                err = self.input_view.primer_input.validation_errors[idx]
                if err.get("name") or err.get("seq"):
                    has_invalid_selected = True
                    break

        if hasattr(self.input_view.primer_input, "error_banner"):
            self.input_view.primer_input.error_banner.visible = (
                has_invalid_selected
            )

        pcr_is_enabled = (
            has_template and has_enough_primers and not has_invalid_selected
        )

        btn = self.pcr_button_ref.current
        if btn:
            btn.disabled = not pcr_is_enabled
            btn.text = "PCR"

        visible_btn = self.visible_pcr_button_ref.current
        if visible_btn:
            visible_btn.disabled = not pcr_is_enabled
            visible_btn.text = "PCR"

        dimers_btn = self.dimers_button_ref.current
        if dimers_btn:
            dimers_btn.disabled = (
                len(active_primers) < 1
            ) or has_invalid_selected

        visible_dimers_btn = self.visible_dimers_button_ref.current
        if visible_dimers_btn:
            visible_dimers_btn.disabled = (
                len(active_primers) < 1
            ) or has_invalid_selected

        self.page.update()

    def on_settings_change(self, e: ft.ControlEvent | None = None) -> None:
        """Handle settings changes from the settings view.

        Applies theme, updates PCR button state, and persists settings.

        Args:
            e: The Flet control event triggering the change.
        """
        self.apply_theme()

        # Only update the active view immediately to prevent lag!
        active_view = self.view_container.content
        if active_view == self.input_view:
            self.input_view.update_ui()
        else:
            self.input_view_dirty = True

        if active_view == self.settings_view:
            self.settings_view.update_ui()

        self.update_pcr_button_state()
        self.settings.save_to_local(self.page)

        # Only redraw/re-simulate active views
        if active_view == self.pcr_view:
            self.pcr_view.run_pcr(keep_cards=True)
        elif active_view == self.dimers_view:
            self.dimers_view.run_analysis()

        self.page.update()

    def _restore_state_from_file(self, path: str) -> None:
        """Restore app state from a YAML file on startup.

        Args:
            path: Path to the YAML state file.
        """
        try:
            with open(path, encoding="utf-8") as f:
                content = f.read()

            parsed_state = yaml.safe_load(content)

            if not isinstance(parsed_state, dict):
                logger.warning("Invalid state file format, ignoring.")
                return

            self._apply_parsed_state(parsed_state)
            logger.info("State loaded successfully from %s", path)
        except (OSError, ValueError, yaml.YAMLError):
            logger.exception("Error loading state file '%s'", path)

    def _apply_parsed_state(
        self, parsed_state: dict[str, Any], ignore_settings: bool = False
    ) -> None:
        """Apply parsed YAML state to the application.

        Args:
            parsed_state: Parsed YAML dict containing input and settings.
            ignore_settings: If True, settings are not applied.
        """
        if "input" in parsed_state:
            self.input_data.from_dict(parsed_state["input"])
        else:
            self.input_data.from_dict(parsed_state)
        if not ignore_settings and "settings" in parsed_state:
            self.settings.from_dict(parsed_state["settings"])
            self.settings.save_to_local(self.page)
        self.apply_theme()
        self.input_view.update_ui()
        self.settings_view.update_ui()
        self.update_pcr_button_state()
        self.page.update()

    async def save_state(self, e: ft.Event[ft.Control]) -> None:
        """Save app state to YAML configuration file."""
        if self.filepicker_open:
            return
        self.filepicker_open = True
        try:
            self.input_view.sync_to_state()
            combined: dict[str, object] = {
                "input": self.input_data.to_dict(),
            }
            yaml_str = serialise_state(combined)

            from amplifyp.gui.utils.io import save_and_write_file

            await save_and_write_file(
                page=self.page,
                dialog_title="Save all",
                file_name="amplify_gui_state.yaml",
                allowed_extensions=["yaml", "yml"],
                content=yaml_str,
                show_notification=self.notification_helper.show_message,
                success_message_desktop="State saved successfully!",
                success_message_web="State ready for download!",
            )
        except (OSError, ValueError) as ex:
            self.notification_helper.show_message(f"Error saving state: {ex}")
        finally:
            self.filepicker_open = False

    async def load_state(self, e: ft.Event[ft.Control]) -> None:
        """Load app state from YAML configuration file."""
        if self.filepicker_open:
            return
        self.filepicker_open = True
        try:
            from amplifyp.gui.utils.io import pick_and_read_file

            content = await pick_and_read_file(
                page=self.page,
                dialog_title="Load all",
                allowed_extensions=["yaml", "yml"],
                show_notification=self.notification_helper.show_message,
            )
            if content is None:
                return

            parsed_state = yaml.safe_load(content)

            if not isinstance(parsed_state, dict):
                self.notification_helper.show_message(
                    "Error: Invalid state file format."
                )
                return

            self._apply_parsed_state(parsed_state, ignore_settings=True)
            self.notification_helper.show_message("State loaded successfully!")
        except (OSError, ValueError, yaml.YAMLError) as ex:
            logger.exception("Error loading state:")
            self.notification_helper.show_message(f"Error loading state: {ex}")
        finally:
            self.filepicker_open = False

    def switch_view(self, _e: ft.Event[ft.Control], view: ft.Control) -> None:
        """Switch the main view container to display a different view.

        Updates the container content and configures resize handlers
        appropriate for the target view.

        Args:
            _e: The event that triggered the view switch (unused).
            view: The Flet control to display as the new view.
        """
        if view == self.input_view and self.input_view_dirty:
            self.input_view.update_ui()
            self.input_view_dirty = False

        self.view_container.content = view
        is_input = view == self.input_view
        self.visible_save_btn_control.visible = is_input
        self.visible_load_btn_control.visible = is_input
        self.visible_header_divider.visible = is_input

        if view == self.input_view:
            self.page.on_resize = self.input_view._handle_resize
        elif view == self.pcr_view:
            self.page.on_resize = self.pcr_view._handle_resize
        else:
            self.page.on_resize = None

        self.page.update()

    def confirm_dismiss(self, e: ft.ControlEvent) -> None:
        """Close confirmation dialog."""
        dialog = self._confirm_dialog
        if dialog:
            dialog.open = False
            self.page.update()

    async def confirm_exit_async(self) -> None:
        """Asynchronously destroy the application window.

        Used as a task handler for safe window closure on desktop.
        """
        try:
            await self.page.window.destroy()
        except RuntimeError:
            logger.debug("Window already closed, skipping destroy")

    def confirm_exit(self, e: ft.ControlEvent) -> None:
        """Launch the async window destruction task.

        Args:
            e: The Flet control event triggering the exit.
        """
        self.page.run_task(self.confirm_exit_async)

    async def _restore_state_and_auto_close_async(self) -> None:
        """Restore state from file and run auto-close sequence if requested."""
        # Yield to let the page finish initial rendering and attach controls
        await asyncio.sleep(0)
        if self.state_file:
            self._restore_state_from_file(self.state_file)
        if self.auto_close and self.state_file:
            await self._auto_close_and_quit_delayed()

    async def _auto_close_and_quit_delayed(
        self, _event: ft.ControlEvent | None = None
    ) -> None:
        """Run PCR/dimer analysis then quit for performance regression testing.

        Yields to the event loop to let initial render complete, runs analysis,
        then destroys the window automatically.

        Args:
            _event: Unused event parameter for task compatibility.
        """
        try:
            # Yield to event loop to let the initial page render complete
            await asyncio.sleep(0)

            self.update_pcr_button_state(sync=False)

            pcr_btn = self.pcr_button_ref.current
            if pcr_btn and not pcr_btn.disabled:
                self.pcr_view.run_pcr()

            dimers_btn = self.dimers_button_ref.current
            if dimers_btn and not dimers_btn.disabled:
                self.dimers_view.run_analysis()

            self.page.update()
            # Give Flet/Flutter rendering engine time to finish the pass
            await asyncio.sleep(1)

            await self.confirm_exit_async()
        except Exception:
            logger.exception("Error during auto-close sequence")
            try:
                await self.confirm_exit_async()
            except RuntimeError:
                pass

    def on_window_event(self, e: ft.WindowEvent) -> None:
        """Handle desktop window events, showing close confirmation dialog.

        When the user attempts to close the application window, this method
        displays a confirmation dialog to prevent accidental data loss.

        Args:
            e: The window event containing close information.
        """
        if (
            e.data == "close"
            or getattr(e, "type", None) == ft.WindowEventType.CLOSE
        ):
            dialog = self._confirm_dialog
            if not dialog:
                dialog = ft.AlertDialog(
                    modal=True,
                    title=ft.Text("Confirm Exit"),
                    content=ft.Text(
                        "Are you sure you want to close AmplifyP? "
                        "Unsaved changes will be lost."
                    ),
                    actions=[
                        ft.TextButton("Yes", on_click=self.confirm_exit),
                        ft.TextButton("No", on_click=self.confirm_dismiss),
                    ],
                    actions_alignment=ft.MainAxisAlignment.END,
                )
                self._confirm_dialog = dialog

            if dialog not in self.page.overlay:
                self.page.overlay.append(dialog)
            dialog.open = True
            self.page.update()

    def on_update_found(self, latest_version: str) -> None:
        """Update header version text when a new version is found."""
        if hasattr(self, "header") and self.header:
            self.header.set_update_available(latest_version)

    async def check_updates_async(self) -> None:
        """Run update checking asynchronously without blocking main thread."""
        from amplifyp import __version__ as current_version
        from amplifyp.gui.utils.version_check import (
            fetch_latest_release_version,
            is_newer_version,
            should_check_for_updates,
        )

        frequency = self.settings.get(
            "version_checking_frequency", "Once per Month"
        )
        try:
            last_check = float(
                self.settings.get("last_version_check_timestamp", 0.0)
            )
        except (TypeError, ValueError):
            last_check = 0.0
        current_time = float(time.time())

        if not should_check_for_updates(frequency, last_check, current_time):
            return

        loop = asyncio.get_running_loop()
        latest_tag = await loop.run_in_executor(
            None, fetch_latest_release_version
        )

        if latest_tag is not None:
            # Update last check timestamp
            self.settings["last_version_check_timestamp"] = current_time
            self.settings.save_to_local(self.page)

            if is_newer_version(latest_tag, current_version):
                self.on_update_found(latest_tag)
