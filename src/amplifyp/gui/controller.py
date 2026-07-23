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
from pathlib import Path
from typing import Any, cast

import flet as ft  # type: ignore[import-not-found, unused-ignore]

from amplifyp.gui.controllers import (
    NavigationManager,
    ThemeManager,
    UpdateManager,
)
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.gui_helpers import NotificationHelper
from amplifyp.gui.views import (
    AboutView,
    Designer2DView,
    DimerView,
    InputView,
    PCRView,
    PrimerDesignerView,
    SettingsView,
)
from amplifyp.gui.views.settings.primer_list_tile import PrimerListTile

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
        self._clear_dialog = None

        # Refs for buttons
        self.pcr_button_ref = ft.Ref[ft.FilledButton]()
        self.dimers_button_ref = ft.Ref[ft.FilledButton]()

        # Views placeholders
        self.input_view: InputView = cast(InputView, None)
        self.settings_view: SettingsView = cast(SettingsView, None)
        self.about_view: AboutView = cast(AboutView, None)
        self.pcr_view: PCRView = cast(PCRView, None)
        self.dimers_view: DimerView = cast(DimerView, None)
        self.designer_view: PrimerDesignerView = cast(PrimerDesignerView, None)
        self.designer_2d_view: Designer2DView = cast(Designer2DView, None)
        self.view_container: ft.Container = cast(ft.Container, None)
        self.header_container: ft.Container = cast(ft.Container, None)

        # UI Control placeholders
        self.save_btn_control: ft.FilledButton = cast(ft.FilledButton, None)
        self.clear_btn_control: ft.FilledButton = cast(ft.FilledButton, None)
        self.load_btn_control: ft.FilledButton = cast(ft.FilledButton, None)
        self.header_divider: ft.Container = cast(ft.Container, None)
        self.notification_helper: NotificationHelper = cast(
            NotificationHelper, None
        )
        self.input_view_dirty = False

        # Sub-controllers
        self._theme_manager = ThemeManager(self)
        self._nav_manager = NavigationManager(self)
        self._update_manager = UpdateManager(self)

    def initialise(self) -> None:
        """Configure page setup, window events, views, and custom layout."""
        self._configure_page_and_window()

        self.settings.load_from_local(self.page)
        self.page.on_platform_brightness_change = (
            self._theme_manager.on_platform_brightness_change
        )
        self.apply_theme()
        self.page.on_keyboard_event = self._on_keyboard_event

        def handle_input_change(e: ft.ControlEvent | None) -> None:
            """Handle input changes to update button states.

            Args:
                e: The Flet control event.
            """
            self.update_pcr_button_state(sync=False)

        def handle_stop_editing(e: ft.ControlEvent | None) -> None:
            """Save state if on web when editing is stopped.

            Args:
                e: The Flet control event.
            """
            if self.page.web:
                self.save_last_state()

        # Instantiate views
        self.input_view = InputView(
            self.page,
            self.input_data,
            self.settings,
            on_change=handle_input_change,
            on_stop_editing=handle_stop_editing,
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
        self.designer_view = PrimerDesignerView(
            self.page, self.input_data, self.settings
        )
        self.designer_2d_view = Designer2DView(
            self.page, self.input_data, self.settings
        )
        self.about_view = AboutView(self.page, self.settings)

        self.notification_helper = NotificationHelper(self.page)

        # Load state and handle auto-close asynchronously once page is ready
        if self.state_file:
            self.page.run_task(self._restore_state_and_auto_close_async)
        else:
            self.page.run_task(self._load_last_state_async)

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
                self.page.run_javascript(  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
                    """
                    window.addEventListener('beforeunload', (event) => {
                        event.preventDefault();
                        event.returnValue = '';
                    });
                    """
                )
        else:
            self.page.window.width = 1280
            self.page.window.height = 720
            self.page.window.prevent_close = False
            if not self.auto_close:
                self.page.window.prevent_close = True
                self._confirm_dialog = None
                self.page.window.on_event = self.on_window_event

    def _setup_navigation_controls(self) -> None:
        """Configure navigation controls for the main application window."""
        self._nav_manager.setup_navigation_controls()

    def apply_theme(self) -> None:
        """Apply theme settings (light/dark/system mode) to the page."""
        self._theme_manager.apply_theme()

    def on_platform_brightness_change(
        self, e: ft.ControlEvent | None = None
    ) -> None:
        """Handle system brightness shifts."""
        self._theme_manager.on_platform_brightness_change(e)

    def switch_view(self, _e: ft.Event[ft.Control], view: ft.Control) -> None:
        """Switch the main view container to display a different view.

        Updates the container content and configures resize handlers
        appropriate for the target view.

        Args:
            _e: The event that triggered the view switch (unused).
            view: The Flet control to display as the new view.
        """
        self._nav_manager.switch_view(_e, view)

    def on_pcr_click(self, e: ft.ControlEvent) -> None:
        """Handle PCR click: switch view then run PCR.

        The view is switched first so the diagram canvas renders
        while the PCR view is the active content.  This avoids
        Flet's diff algorithm marking canvas shapes as 'already
        sent' before the view becomes visible.
        """
        self._nav_manager.on_pcr_click(e)

    def on_dimers_click(self, e: ft.ControlEvent) -> None:
        """Handle dimers click: switch view then run analysis."""
        self._nav_manager.on_dimers_click(e)

    def update_pcr_button_state(
        self, sync: bool = True, update_page: bool = True
    ) -> None:
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
            btn.text = "PCR"  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]

        dimers_btn = self.dimers_button_ref.current
        if dimers_btn:
            dimers_btn.disabled = (
                len(active_primers) < 1
            ) or has_invalid_selected

        if update_page:
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
            # Immediately reposition the info panel if its position changed
            if e is not None:
                ctrl = getattr(e, "control", None)
                if ctrl is not None:
                    try:
                        tile = getattr(
                            self.settings_view, "primer_list_tile", None
                        )
                        if (
                            tile is not None
                            and isinstance(tile, PrimerListTile)
                            and tile.set_primer_info_panel_position is ctrl
                        ):
                            self.input_view.reposition_primer_info_panel()
                            return
                    except Exception as ex:
                        logger.debug("Failed to check primer list tile: %s", ex)
        else:
            self.input_view_dirty = True

        if active_view == self.settings_view:
            self.settings_view.update_ui()

        self.update_pcr_button_state(update_page=False)
        self.settings.save_to_local(self.page)

        # Only redraw/re-simulate active views
        if active_view == self.pcr_view:
            self.pcr_view.run_pcr(keep_cards=True)
        elif active_view == self.dimers_view:
            self.dimers_view.run_analysis()
        else:
            self.page.update()

    def _get_last_state_path(self) -> Path:
        """Get the OS-specific path for the last saved GUI state.

        Returns:
            Path object pointing to the last_state.yaml file location.
        """
        from amplifyp.gui.utils.data_helpers import get_last_state_path

        return get_last_state_path(self)

    def save_last_state(self) -> None:
        """Save the last template and primers to local/platform storage."""
        from amplifyp.gui.utils.data_helpers import save_last_state

        save_last_state(self)

    def load_last_state(self) -> None:
        """Load the last template and primers from local/platform storage."""
        from amplifyp.gui.utils.data_helpers import load_last_state

        load_last_state(self)

    async def _load_last_state_async(self) -> None:
        """Asynchronously load the last template and primers from storage."""
        # Yield to let the page finish initial rendering and attach controls
        await asyncio.sleep(0)
        self.load_last_state()

    def _restore_state_from_file(self, path: str) -> None:
        """Restore app state from a YAML file on startup.

        Args:
            path: Path to the YAML state file.
        """
        from amplifyp.gui.utils.data_helpers import restore_state_from_file

        restore_state_from_file(self, path)

    def _apply_parsed_state(
        self, parsed_state: dict[str, Any], ignore_settings: bool = False
    ) -> None:
        """Apply parsed YAML state to the application.

        Args:
            parsed_state: Parsed YAML dict containing input and settings.
            ignore_settings: If True, settings are not applied.
        """
        from amplifyp.gui.utils.data_helpers import apply_parsed_state

        apply_parsed_state(self, parsed_state, ignore_settings)

    async def save_state(self, e: ft.Event[ft.Control]) -> None:
        """Save app state to YAML configuration file."""
        from amplifyp.gui.utils.data_helpers import save_state

        await save_state(self, e)

    async def load_state(self, e: ft.Event[ft.Control]) -> None:
        """Load app state from YAML configuration file."""
        from amplifyp.gui.utils.data_helpers import load_state

        await load_state(self, e)

    def confirm_dismiss(self, e: ft.ControlEvent) -> None:
        """Close close confirmation dialogue."""
        from amplifyp.gui.utils.system import confirm_dismiss

        confirm_dismiss(self, e)

    async def confirm_exit_async(self) -> None:
        """Asynchronously destroy the application window."""
        from amplifyp.gui.utils.system import confirm_exit_async

        await confirm_exit_async(self)

    def confirm_exit(self, e: ft.ControlEvent) -> None:
        """Launch the async window destruction task."""
        from amplifyp.gui.utils.system import confirm_exit

        confirm_exit(self, e)

    async def _restore_state_and_auto_close_async(self) -> None:
        """Restore state from file and run auto-close sequence if requested."""
        from amplifyp.gui.utils.system import (
            restore_state_and_auto_close_async,
        )

        await restore_state_and_auto_close_async(self)

    async def _auto_close_and_quit_delayed(
        self, _event: ft.ControlEvent | None = None
    ) -> None:
        """Run PCR/dimer analysis then quit for regression testing."""
        from amplifyp.gui.utils.system import auto_close_and_quit_delayed

        await auto_close_and_quit_delayed(self, _event)

    def on_window_event(self, e: ft.WindowEvent) -> None:
        """Handle desktop window events, showing close confirmation dialog."""
        from amplifyp.gui.utils.system import on_window_event

        on_window_event(self, e)

    def on_update_found(self, latest_version: str) -> None:
        """Update header version text when a new version is found."""
        self._update_manager.on_update_found(latest_version)

    async def check_updates_async(self) -> None:
        """Run update checking asynchronously without blocking main thread."""
        await self._update_manager.check_updates_async()

    def _confirm_clear(self, _ev: ft.ControlEvent) -> None:
        """Show the clear confirmation dialogue.

        Args:
            _ev: The Flet control event.
        """
        from amplifyp.gui.utils.data_helpers import confirm_clear

        confirm_clear(self, _ev)

    def _dismiss_clear(self, _ev: ft.ControlEvent) -> None:
        """Dismiss the clear confirmation dialogue.

        Args:
            _ev: The Flet control event.
        """
        from amplifyp.gui.utils.data_helpers import dismiss_clear

        dismiss_clear(self, _ev)

    def clear_all(self, e: ft.ControlEvent) -> None:
        """Show a confirmation dialogue before clearing inputs.

        Clears all template sequences and primers if confirmed.
        """
        from amplifyp.gui.utils.data_helpers import clear_all

        clear_all(self, e)

    def _on_keyboard_event(self, e: ft.KeyboardEvent) -> None:
        """Handle global keyboard events for primer navigation.

        Arrow Up/Down navigates between primer rows.
        Arrow Left/Right navigates between name and sequence fields
        in the same row.
        """
        from amplifyp.gui.utils.gui_helpers import handle_keyboard_event

        handle_keyboard_event(self, e)
