# Copyright (C) 2026 Fufu Fang
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

import traceback
from typing import cast

import flet as ft
import yaml

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import NotificationHelper, get_version, serialize_state
from amplifyp.gui.views import (
    DimerView,
    InputView,
    PCRView,
    SettingsView,
)


class GUIController:
    """Manages GUI state, event handlers, views and main orchestration."""

    def __init__(self, page: ft.Page) -> None:
        """Initialize the GUIController."""
        self.page = page
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
        self.pcr_view: PCRView = cast(PCRView, None)
        self.dimers_view: DimerView = cast(DimerView, None)
        self.view_container: ft.Container = cast(ft.Container, None)

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

    def initialize(self) -> None:
        """Configure page setup, window events, views, and custom layout."""
        self.page.overlay.clear()
        self.page.title = "AmplifyP"
        self.page.vertical_alignment = ft.MainAxisAlignment.START
        self.page.fonts = {"Roboto Mono": "fonts/RobotoMono-Regular.ttf"}
        self.page.padding = 0
        self.page.spacing = 0
        self.page.window.icon = "images/icon.png"

        # Handle close / reload warnings
        if self.page.web:
            if hasattr(self.page, "run_javascript"):
                self.page.run_javascript("""
                    window.addEventListener('beforeunload', (event) => {
                        event.preventDefault();
                        event.returnValue = '';
                    });
                """)
        else:
            self.page.window.prevent_close = False
            self.page.window.prevent_close = True
            self.page._confirm_dialog = None
            self.page.window.on_event = self.on_window_event

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
            on_stop_editing=lambda: self.update_pcr_button_state(sync=False),
        )
        self.settings_view = SettingsView(
            self.page,
            self.settings,
            on_change=self.on_settings_change,
            on_apply=self.run_apply_settings,
            on_reset=self.run_apply_settings,
        )
        self.pcr_view = PCRView(self.page, self.input_data, self.settings)
        self.dimers_view = DimerView(self.page, self.input_data, self.settings)

        self.notification_helper = NotificationHelper(self.page)

        # Main view container
        self.view_container = ft.Container(
            content=self.input_view, expand=True, padding=10
        )

        # Header controls & routing buttons setup
        self._setup_navigation_controls()

    def _setup_navigation_controls(self) -> None:
        # AppBar buttons (for test compatibility)
        input_button = ft.FilledButton(
            "Input",
            icon=ft.Icons.INPUT,
            on_click=lambda e: self.switch_view(e, self.input_view),
            tooltip="Input",
        )
        input_button.content_description = "Input"

        pcr_button = ft.FilledButton(
            "PCR",
            ref=self.pcr_button_ref,
            on_click=self.on_pcr_click,
            disabled=True,
            icon=ft.Icons.ANALYTICS,
            tooltip="PCR",
        )
        pcr_button.content_description = "PCR"

        dimers_button = ft.FilledButton(
            "Primer Dimers",
            ref=self.dimers_button_ref,
            on_click=self.on_dimers_click,
            disabled=True,
            icon=ft.Icons.COMPARE_ARROWS,
            tooltip="Primer Dimers",
        )
        dimers_button.content_description = "Primer Dimers"

        settings_button = ft.FilledButton(
            "Settings",
            icon=ft.Icons.SETTINGS,
            on_click=lambda e: self.switch_view(e, self.settings_view),
            tooltip="Settings",
        )
        settings_button.content_description = "Settings"

        save_btn_control = ft.FilledButton(
            "Save all",
            icon=ft.Icons.SAVE,
            tooltip="Save all",
            on_click=self.save_state,
        )
        save_btn_control.content_description = "Save all"

        load_btn_control = ft.FilledButton(
            "Load all",
            icon=ft.Icons.UPLOAD_FILE,
            tooltip="Load all",
            on_click=self.load_state,
        )
        load_btn_control.content_description = "Load all"

        self.page.appbar = ft.AppBar(
            visible=False,
            actions=[
                input_button,
                pcr_button,
                dimers_button,
                settings_button,
                save_btn_control,
                load_btn_control,
            ],
        )

        # Visible navigation controls (Responsive header)
        visible_input_button = ft.FilledButton(
            "Input",
            icon=ft.Icons.INPUT,
            on_click=lambda e: self.switch_view(e, self.input_view),
            tooltip="Input",
        )
        visible_input_button.content_description = "Input"

        visible_pcr_button = ft.FilledButton(
            "PCR",
            ref=self.visible_pcr_button_ref,
            on_click=self.on_pcr_click,
            disabled=True,
            icon=ft.Icons.ANALYTICS,
            tooltip="PCR",
        )
        visible_pcr_button.content_description = "PCR"

        visible_dimers_button = ft.FilledButton(
            "Primer Dimers",
            ref=self.visible_dimers_button_ref,
            on_click=self.on_dimers_click,
            disabled=True,
            icon=ft.Icons.COMPARE_ARROWS,
            tooltip="Primer Dimers",
        )
        visible_dimers_button.content_description = "Primer Dimers"

        visible_settings_button = ft.FilledButton(
            "Settings",
            icon=ft.Icons.SETTINGS,
            on_click=lambda e: self.switch_view(e, self.settings_view),
            tooltip="Settings",
        )
        visible_settings_button.content_description = "Settings"

        self.visible_save_btn_control = ft.FilledButton(
            "Save all",
            icon=ft.Icons.SAVE,
            tooltip="Save all",
            on_click=self.save_state,
        )
        self.visible_save_btn_control.content_description = "Save all"

        self.visible_load_btn_control = ft.FilledButton(
            "Load all",
            icon=ft.Icons.UPLOAD_FILE,
            tooltip="Load all",
            on_click=self.load_state,
        )
        self.visible_load_btn_control.content_description = "Load all"

        self.visible_header_divider = ft.Container(
            width=1,
            height=20,
            bgcolor=GUIColours.OUTLINE,
        )

        app_version = get_version()
        version_text = ft.Text(
            app_version,
            size=14,
            color=GUIColours.TEXT_ON_SURFACE,
            opacity=0.5,
            weight=ft.FontWeight.W_400,
            selectable=True,
        )

        header_container = ft.Container(
            content=ft.ResponsiveRow(
                [
                    ft.Container(
                        content=ft.Row(
                            [
                                ft.Image(
                                    src="images/favicon.png",
                                    height=32,
                                    fit="contain",
                                ),
                                ft.Text(
                                    "AmplifyP",
                                    size=20,
                                    weight=ft.FontWeight.BOLD,
                                ),
                                ft.Container(width=12),
                                version_text,
                            ],
                            spacing=8,
                            tight=True,
                            vertical_alignment=ft.CrossAxisAlignment.CENTER,
                        ),
                        col={"lg": 4, "md": 12, "sm": 12, "xs": 12},
                        alignment=ft.Alignment(-1, 0),
                    ),
                    ft.Container(
                        content=ft.Row(
                            [
                                visible_input_button,
                                visible_pcr_button,
                                visible_dimers_button,
                                visible_settings_button,
                                self.visible_header_divider,
                                self.visible_save_btn_control,
                                self.visible_load_btn_control,
                            ],
                            spacing=10,
                            tight=True,
                            wrap=True,
                            vertical_alignment=ft.CrossAxisAlignment.CENTER,
                        ),
                        col={"lg": 8, "md": 12, "sm": 12, "xs": 12},
                        alignment=ft.Alignment(1, 0),
                    ),
                ],
                vertical_alignment=ft.CrossAxisAlignment.CENTER,
            ),
            padding=ft.Padding(16, 8, 16, 8),
            bgcolor=GUIColours.SURFACE,
        )

        self.page.add(
            ft.Divider(height=1, thickness=1),
            self.view_container,
        )
        self.page.controls.insert(0, header_container)
        self.page.on_resize = self.input_view._handle_resize

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

    def on_platform_brightness_change(self, e: ft.ControlEvent) -> None:
        """Handle system brightness shifts."""
        self.apply_theme()
        self.page.update()

    def on_pcr_click(self, e: ft.ControlEvent) -> None:
        """Handle PCR click: run PCR and switch view if successful."""
        if self.pcr_view.run_pcr():
            self.switch_view(e, self.pcr_view)
            if self.pcr_button_ref.current:
                self.pcr_button_ref.current.text = "PCR"
            if self.visible_pcr_button_ref.current:
                self.visible_pcr_button_ref.current.text = "PCR"
            self.page.update()

    def on_dimers_click(self, e: ft.ControlEvent) -> None:
        """Handle dimers click: run analysis and switch view if successful."""
        if self.dimers_view.run_analysis():
            self.switch_view(e, self.dimers_view)
            self.page.update()

    def update_pcr_button_state(self, sync: bool = True) -> None:
        """Enable PCR and dimers buttons only if input is valid."""
        if sync:
            self.input_view.sync_to_state()
        has_template = bool(self.input_data.template.strip())
        active_primers = self.input_data.get_active_primers()
        has_enough_primers = len(active_primers) >= 2

        # Check if any selected (active) primer has validation errors
        has_invalid_selected = False
        for idx, p in enumerate(self.input_data.primers):
            if p.get("active", False):
                if idx < len(self.input_view.primer_input.validation_errors):
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

    def on_settings_change(self, e: ft.ControlEvent) -> None:
        """Handle settings changes."""
        self.apply_theme()
        self.update_pcr_button_state()
        self.settings.save_to_local(self.page)

    def run_apply_settings(self, e: ft.ControlEvent) -> None:
        """Apply the settings updates."""
        self.apply_theme()
        self.update_pcr_button_state()
        self.settings.save_to_local(self.page)

    async def save_state(self, e: ft.ControlEvent) -> None:
        """Save app state to YAML configuration file."""
        if self.filepicker_open:
            return
        self.filepicker_open = True
        try:
            self.input_view.sync_to_state()
            self.settings_view.sync_to_state()
            combined: dict[str, object] = {
                "input": self.input_data.to_dict(),
                "settings": self.settings.to_dict(),
            }
            yaml_str = serialize_state(combined)

            from amplifyp.gui.util import save_and_write_file

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
        except Exception as ex:
            self.notification_helper.show_message(f"Error saving state: {ex}")
        finally:
            self.filepicker_open = False

    async def load_state(self, e: ft.ControlEvent) -> None:
        """Load app state from YAML configuration file."""
        if self.filepicker_open:
            return
        self.filepicker_open = True
        try:
            from amplifyp.gui.util import pick_and_read_file

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

            if "input" in parsed_state:
                self.input_data.from_dict(parsed_state["input"])
            else:
                # Legacy format: input data at top level
                self.input_data.from_dict(parsed_state)
            if "settings" in parsed_state:
                self.settings.from_dict(parsed_state["settings"])
                self.settings.save_to_local(self.page)
            self.apply_theme()
            self.input_view.update_ui()
            self.settings_view.update_ui()
            self.update_pcr_button_state()
            self.notification_helper.show_message("State loaded successfully!")
        except Exception as ex:
            print("LOAD STATE ERROR:", traceback.format_exc())
            self.notification_helper.show_message(f"Error loading state: {ex}")
        finally:
            self.filepicker_open = False

    def switch_view(self, e: ft.ControlEvent, view: ft.Control) -> None:
        """Switch routing container views."""
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
        dialog = getattr(self.page, "_confirm_dialog", None)
        if dialog:
            dialog.open = False
            self.page.update()

    async def confirm_exit_async(self) -> None:
        """Destroy window helper."""
        try:
            await self.page.window.destroy()
        except Exception:  # noqa: S110
            pass

    def confirm_exit(self, e: ft.ControlEvent) -> None:
        """Exit task launcher."""
        self.page.run_task(self.confirm_exit_async)

    def on_window_event(self, e: ft.WindowEvent) -> None:
        """Handle desktop close warning dialogues."""
        if (
            e.data == "close"
            or getattr(e, "type", None) == ft.WindowEventType.CLOSE
        ):
            dialog = getattr(self.page, "_confirm_dialog", None)
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
                self.page._confirm_dialog = dialog

            if dialog not in self.page.overlay:
                self.page.overlay.append(dialog)
            dialog.open = True
            self.page.update()
