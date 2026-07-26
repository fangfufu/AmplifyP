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

"""GUI controller for orchestrating views, state, and window events."""

from __future__ import annotations

import asyncio
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

from PySide6.QtCore import QTimer
from PySide6.QtGui import QCloseEvent, QFontDatabase
from PySide6.QtWidgets import (
    QMainWindow,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from amplifyp.gui2.controllers import (
    NavigationManager,
    ThemeManager,
    UpdateManager,
)
from amplifyp.gui2.settings import GUISettings
from amplifyp.gui2.user_data import GUIInput

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


class GUIController:
    """Manages GUI state, event handlers, views and main orchestration."""

    def __init__(
        self,
        state_file: str | None = None,
        auto_close: bool = False,
        export_screenshots: bool = False,
        screenshots_dir: str | None = None,
        window_width: int | None = None,
        window_height: int | None = None,
    ) -> None:
        """Initialize the GUIController.

        Args:
            state_file: Optional path to a YAML state file to load on startup.
            auto_close: If True, quit automatically after rendering is complete.
            export_screenshots: Save PNG screenshots of views.
            screenshots_dir: Target directory for saved PNG screenshots.
            window_width: Optional application window width in pixels.
            window_height: Optional application window height in pixels.
        """
        self.state_file = state_file
        self.auto_close = auto_close
        self.export_screenshots = export_screenshots
        self.screenshots_dir = (
            Path(screenshots_dir) if screenshots_dir else None
        )
        self.window_width = window_width
        self.window_height = window_height
        self.input_data = GUIInput()
        self.settings = GUISettings()
        self.filepicker_open = False
        self._confirm_dialog = None
        self._clear_dialog = None
        self.input_view_dirty = False

        # Views (populated in initialise)
        self.input_view: Any = None
        self.settings_view: Any = None
        self.about_view: Any = None
        self.pcr_view: Any = None
        self.dimers_view: Any = None
        self.designer_view: Any = None
        self.designer_2d_view: Any = None
        self.current_view: Any = None

        # UI controls
        self.pcr_button_ref: Any = None
        self.dimers_button_ref: Any = None
        self.header: Any = None
        self.header_container: Any = None
        self.save_btn_control: Any = None
        self.clear_btn_control: Any = None
        self.load_btn_control: Any = None
        self.header_divider: Any = None
        self.view_container: QStackedWidget | None = None
        self.main_window: QMainWindow | None = None

        # Sub-controllers
        self._theme_manager = ThemeManager(self)
        self._nav_manager = NavigationManager(self)
        self._update_manager = UpdateManager(self)

    def initialise(self) -> None:
        """Configure window setup, views, and layout."""
        self.settings.load_from_local()
        # Theme will be applied after header is created

        # Create main window
        self.main_window = QMainWindow()
        self.main_window.setWindowTitle("AmplifyP")
        self.main_window.setMinimumSize(800, 600)
        self.main_window.resize(
            self.window_width or 1280,
            self.window_height or 720,
        )

        # Load Roboto Mono font
        self._load_font()

        # Instantiate views
        self.input_view = self._create_input_view()
        self.settings_view = self._create_settings_view()
        self.pcr_view = self._create_pcr_view()
        self.dimers_view = self._create_dimer_view()
        self.designer_view = self._create_designer_view()
        self.designer_2d_view = self._create_designer_2d_view()
        self.about_view = self._create_about_view()

        # View container (stacked widget)
        self.view_container = QStackedWidget()
        self.view_container.addWidget(self.input_view)
        self.current_view = self.input_view

        # Set up navigation
        self._nav_manager.setup_navigation_controls()

        # Apply theme after header is created
        self.apply_theme()

        # Set initial view (shows/hides save/clear/load buttons)
        self._nav_manager.switch_view(self.input_view)

        # Central widget
        central_widget = QWidget()
        layout = QVBoxLayout(central_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.header_container)
        layout.addWidget(self.view_container)
        self.main_window.setCentralWidget(central_widget)

        # Handle window close
        self.main_window.closeEvent = self._on_close_event

        # Load state
        if self.state_file:
            QTimer.singleShot(0, lambda: self._restore_state_and_auto_close())
        else:
            QTimer.singleShot(0, self._load_last_state)

        # Check for updates in background
        QTimer.singleShot(
            0, lambda: asyncio.create_task(self.check_updates_async())
        )

        self.main_window.show()

    def _load_font(self) -> None:
        """Load the Roboto Mono font from assets."""
        import os

        font_path = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "assets",
            "fonts",
            "RobotoMono-Regular.ttf",
        )
        if os.path.exists(font_path):
            font_id = QFontDatabase.addApplicationFont(font_path)
            if font_id != -1:
                font_family = QFontDatabase.applicationFontFamilies(font_id)[0]
                from PySide6.QtGui import QFont

                font = QFont(font_family)
                strategy = getattr(
                    QFont.StyleStrategy, "PreferNoAntialias", None
                )
                if strategy is not None:
                    font.setStyleStrategy(strategy)  # type: ignore[arg-type]
                self.main_window.setFont(font)  # type: ignore[union-attr]

    def apply_theme(self) -> None:
        """Apply theme settings (light/dark/system mode)."""
        self._theme_manager.apply_theme()

    def switch_view(self, view: Any) -> None:
        """Switch the main view container."""
        self._nav_manager.switch_view(view)
        self.current_view = view

    def update_pcr_button_state(
        self, sync: bool = True, update_page: bool = True
    ) -> None:
        """Enable PCR and dimers buttons only if input is valid."""
        if sync:
            self.input_view.sync_to_state()

        has_template = bool(self.input_data.template.strip())
        active_primers = self.input_data.get_active_primers()
        has_enough_primers = len(active_primers) >= 1

        has_invalid_selected = False
        for idx, p in enumerate(self.input_data.primers):
            if p.get("active", False) and idx < len(
                getattr(self.input_view.primer_input, "validation_errors", [])
            ):
                err = self.input_view.primer_input.validation_errors[idx]
                if err.get("name") or err.get("seq"):
                    has_invalid_selected = True
                    break

        pcr_is_enabled = (
            has_template and has_enough_primers and not has_invalid_selected
        )

        btn = self.pcr_button_ref
        if btn is not None:
            btn.setEnabled(pcr_is_enabled)

        dimers_btn = self.dimers_button_ref
        if dimers_btn is not None:
            dimers_btn.setEnabled(
                len(active_primers) >= 1 and not has_invalid_selected
            )

    def on_settings_change(self) -> None:
        """Handle settings changes from the settings view."""
        self.apply_theme()
        self.update_pcr_button_state(update_page=False)
        self.settings.save_to_local()

        active_view = self.current_view
        if active_view == self.pcr_view:
            self.pcr_view.run_pcr(keep_cards=True)
        elif active_view == self.dimers_view:
            self.dimers_view.run_analysis()

    def save_last_state(self) -> None:
        """Save the last template and primers to local storage."""
        from amplifyp.gui2.utils.data_helpers import save_last_state

        save_last_state(self)

    def load_last_state(self) -> None:
        """Load the last template and primers from local storage."""
        from amplifyp.gui2.utils.data_helpers import load_last_state

        load_last_state(self)

    def _load_last_state(self) -> None:
        """Asynchronously load the last state."""
        self.load_last_state()

    def _restore_state_from_file(self, path: str) -> None:
        """Restore app state from a YAML file on startup."""
        from amplifyp.gui2.utils.data_helpers import restore_state_from_file

        restore_state_from_file(self, path)

    def _restore_state_and_auto_close(self) -> None:
        """Restore state from file and run auto-close sequence."""
        if self.state_file:
            self._restore_state_from_file(self.state_file)
        if getattr(self, "export_screenshots", False):
            self._export_screenshots()
        elif getattr(self, "auto_close", False) and self.state_file:
            self._auto_close_and_quit()

    def _export_screenshots(self) -> None:
        """Export screenshots of all views."""
        from amplifyp.gui2.utils.system import capture_all_views_async

        asyncio.get_event_loop().run_until_complete(
            capture_all_views_async(self)
        )

    def _auto_close_and_quit(self) -> None:
        """Run PCR/dimer analysis then quit."""
        asyncio.get_event_loop().run_until_complete(
            self._auto_close_and_quit_async()
        )

    async def _auto_close_and_quit_async(self) -> None:
        """Async version of auto-close and quit."""
        await asyncio.sleep(0)
        self.update_pcr_button_state(sync=False)
        if self.pcr_view:
            self.pcr_view.run_pcr()
        if self.dimers_view:
            self.dimers_view.run_analysis()
        await asyncio.sleep(1)
        await self.confirm_exit_async()

    async def confirm_exit_async(self) -> None:
        """Asynchronously destroy the application window."""
        from amplifyp.gui2.utils.system import confirm_exit_async

        await confirm_exit_async(self)

    def save_state(self) -> bool:
        """Save app state to YAML configuration file."""
        from amplifyp.gui2.utils.data_helpers import save_state

        return save_state(self, self._show_notification)

    def load_state(self) -> bool:
        """Load app state from YAML configuration file."""
        from amplifyp.gui2.utils.data_helpers import load_state

        return load_state(self, self._show_notification)

    def clear_all(self) -> None:
        """Show a confirmation dialogue before clearing inputs."""
        from amplifyp.gui2.utils.data_helpers import clear_all

        clear_all(self)

    def _on_close_event(self, event: QCloseEvent) -> None:
        """Handle window close event."""
        if not self.auto_close:
            from amplifyp.gui2.utils.system import on_window_event

            accepted = on_window_event(self)
            if not accepted:
                event.ignore()
                return
        event.accept()

    def on_update_found(self, latest_version: str) -> None:
        """Update header version text when a new version is found."""
        self._update_manager.on_update_found(latest_version)

    async def check_updates_async(self) -> None:
        """Run update checking asynchronously."""
        await self._update_manager.check_updates_async()

    def _show_notification(self, message: str) -> None:
        """Show a notification message."""
        from PySide6.QtWidgets import QMessageBox

        QMessageBox.information(
            self.main_window,
            "AmplifyP",
            message,
            QMessageBox.StandardButton.Ok,
        )

    # ------------------------------------------------------------------
    # View factories - stub implementations, replaced by Phase 2+ views
    # ------------------------------------------------------------------

    def _create_input_view(self) -> QWidget:
        """Create the input view."""
        from amplifyp.gui2.views.input_view import InputView

        view = InputView(
            self.input_data,
            self.settings,
            on_change=self._on_input_change,
        )
        return view

    def _create_settings_view(self) -> QWidget:
        """Create the settings view."""
        from amplifyp.gui2.views.settings_view import SettingsView

        view = SettingsView(
            self.settings,
            on_change=self.on_settings_change,
            on_reset=self.on_settings_change,
            on_update_found=self.on_update_found,
        )
        return view

    def _create_pcr_view(self) -> QWidget:
        """Create the PCR view."""
        from amplifyp.gui2.views.pcr_view import PCRView

        view = PCRView(self.input_data, self.settings)
        return view

    def _create_dimer_view(self) -> QWidget:
        """Create the dimer view."""
        from amplifyp.gui2.views.dimer_view import DimerView

        view = DimerView(self.input_data, self.settings)
        return view

    def _create_designer_view(self) -> QWidget:
        """Create the 1D designer view."""
        from amplifyp.gui2.views.designer_1d_view import Designer1DView

        view = Designer1DView(self.input_data, self.settings)
        return view

    def _create_designer_2d_view(self) -> QWidget:
        """Create the 2D designer view."""
        from amplifyp.gui2.views.designer_2d_view import Designer2DView

        view = Designer2DView(self.input_data, self.settings)
        return view

    def _create_about_view(self) -> QWidget:
        """Create the about view."""
        from amplifyp.gui2.views.about_view import AboutView

        view = AboutView(self.settings)
        return view

    def _on_input_change(self) -> None:
        """Handle input changes."""
        self.update_pcr_button_state(sync=False)
