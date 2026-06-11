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

"""Main Flet application logic."""

import flet as ft
import yaml

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import serialize_state
from amplifyp.gui.views import (
    InputView,
    PrimerDimerView,
    ResultView,
    SettingsView,
)


def main(page: ft.Page) -> None:
    """Main entry point for the Flet application."""
    page.title = "AmplifyP"
    page.vertical_alignment = ft.MainAxisAlignment.START
    page.fonts = {"Roboto Mono": "fonts/RobotoMono-Regular.ttf"}
    page.padding = 0
    page.spacing = 0
    page.window.icon = "images/icon.png"

    # Handle close / reload warnings
    if page.web:
        if hasattr(page, "run_javascript"):
            page.run_javascript("""
                window.addEventListener('beforeunload', (event) => {
                    event.preventDefault();
                    event.returnValue = '';
                });
            """)
    else:
        page.window.prevent_close = True

        def confirm_dismiss(e: ft.ControlEvent) -> None:
            confirm_dialog.open = False
            page.update()

        def confirm_exit(e: ft.ControlEvent) -> None:
            page.run_task(page.window.destroy)

        confirm_dialog = ft.AlertDialog(
            modal=True,
            title=ft.Text("Confirm Exit"),
            content=ft.Text(
                "Are you sure you want to close AmplifyP? "
                "Unsaved changes will be lost."
            ),
            actions=[
                ft.TextButton("Yes", on_click=confirm_exit),
                ft.TextButton("No", on_click=confirm_dismiss),
            ],
            actions_alignment=ft.MainAxisAlignment.END,
        )
        page.overlay.append(confirm_dialog)

        def on_window_event(e: ft.WindowEvent) -> None:
            if (
                e.data == "close"
                or getattr(e, "type", None) == ft.WindowEventType.CLOSE
            ):
                confirm_dialog.open = True
                page.update()

        page.window.on_event = on_window_event

    # Centralize state storage
    input_data = GUIInput()
    settings = GUISettings()
    has_run_pcr = False
    pcr_outdated = False
    pcr_button_ref = ft.Ref[ft.FilledButton]()
    dimers_button_ref = ft.Ref[ft.FilledButton]()

    def apply_theme() -> None:
        """Apply theme settings (light/dark mode) to the page."""
        dark_mode = settings.get("dark_mode", False)
        if dark_mode:
            page.theme_mode = ft.ThemeMode.DARK
            page.bg_color = None
        else:
            page.theme_mode = ft.ThemeMode.LIGHT
            page.bg_color = ft.Colors.WHITE

    # Apply initial theme
    apply_theme()

    def on_pcr_click(e: ft.ControlEvent) -> None:
        """Handle PCR click: run PCR and switch view if successful."""
        nonlocal has_run_pcr, pcr_outdated
        if result_view.run_pcr():
            switch_view(e, result_view)
            has_run_pcr = True
            pcr_outdated = False
            if pcr_button_ref.current:
                pcr_button_ref.current.text = "PCR"
            page.update()

    def on_dimers_click(e: ft.ControlEvent) -> None:
        """Handle dimers click: run analysis and switch view if successful."""
        if dimers_view.run_analysis():
            switch_view(e, dimers_view)
            page.update()

    def update_pcr_button_state() -> None:
        """Enable PCR and dimers buttons only if input is valid."""
        nonlocal pcr_outdated
        input_view.sync_to_state()
        has_template = bool(input_data.template.strip())
        has_primers = len(input_data.get_active_primers()) > 0
        is_enabled = has_template and has_primers

        btn = pcr_button_ref.current
        if btn:
            btn.disabled = not is_enabled

            # Set outdated if enabled and inputs change AFTER a first run
            if is_enabled and has_run_pcr:
                pcr_outdated = True

            label = "PCR *" if (pcr_outdated and is_enabled) else "PCR"
            btn.text = label

        dimers_btn = dimers_button_ref.current
        if dimers_btn:
            dimers_btn.disabled = not has_primers

        page.update()

    def run_analysis_in_background() -> None:
        nonlocal has_run_pcr, pcr_outdated
        result_view.run_pcr()
        dimers_view.run_analysis()
        has_run_pcr = True
        pcr_outdated = False
        if pcr_button_ref.current:
            pcr_button_ref.current.text = "PCR"
        page.update()

    def on_settings_change(e: ft.ControlEvent) -> None:
        apply_theme()
        update_pcr_button_state()

    def run_apply_settings(e: ft.ControlEvent) -> None:
        apply_theme()
        run_analysis_in_background()

    input_view = InputView(
        page,
        input_data,
        settings,
        on_change=lambda e: update_pcr_button_state(),
        on_stop_editing=update_pcr_button_state,
    )
    settings_view = SettingsView(
        page,
        settings,
        on_change=on_settings_change,
        on_apply=run_apply_settings,
        on_reset=run_apply_settings,
    )
    result_view = ResultView(page, input_data, settings)
    dimers_view = PrimerDimerView(page, input_data, settings)

    # Save and Load State
    snack_bar = ft.SnackBar(ft.Text(""), open=False)
    page.overlay.append(snack_bar)

    def show_snackbar(message: str) -> None:
        snack_bar.content = ft.Text(message)
        snack_bar.open = True
        page.update()

    async def save_state(e: ft.ControlEvent) -> None:
        input_view.sync_to_state()
        settings_view.sync_to_state()
        combined: dict[str, object] = {
            "input": input_data.to_dict(),
            "settings": settings.to_dict(),
        }
        yaml_str = serialize_state(combined)

        try:
            file_path = await ft.FilePicker().save_file(
                dialog_title="Save all",
                file_name="amplify_gui_state.yaml",
                allowed_extensions=["yaml", "yml"],
                file_type=ft.FilePickerFileType.CUSTOM,
                src_bytes=yaml_str.encode("utf-8"),
            )
            if page.web:
                show_snackbar("State ready for download!")
            else:
                if file_path is None:
                    return
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(yaml_str)
                show_snackbar("State saved successfully!")
        except Exception as ex:
            show_snackbar(f"Error saving state: {ex}")

    async def load_state(e: ft.ControlEvent) -> None:
        try:
            files = await ft.FilePicker().pick_files(
                dialog_title="Load all",
                allowed_extensions=["yaml", "yml"],
                file_type=ft.FilePickerFileType.CUSTOM,
                with_data=True,
            )
            if not files:
                return

            file = files[0]
            if file.bytes is not None:
                content = file.bytes.decode("utf-8")
                parsed_state = yaml.safe_load(content)
            else:
                if not file.path:
                    show_snackbar("Error: Could not read file content.")
                    return
                with open(file.path, encoding="utf-8") as f:
                    parsed_state = yaml.safe_load(f)

            if not isinstance(parsed_state, dict):
                show_snackbar("Error: Invalid state file format.")
                return

            if "input" in parsed_state:
                input_data.from_dict(parsed_state["input"])
            else:
                # Legacy format: input data at top level
                input_data.from_dict(parsed_state)
            if "settings" in parsed_state:
                settings.from_dict(parsed_state["settings"])
            apply_theme()
            input_view.update_ui()
            settings_view.update_ui()
            update_pcr_button_state()
            show_snackbar("State loaded successfully!")
        except Exception as ex:
            import traceback

            tb = traceback.format_exc()
            print("LOAD STATE ERROR:", tb)
            show_snackbar(f"Error loading state: {ex}")

    view_container = ft.Container(content=input_view, expand=True, padding=10)

    def switch_view(e: ft.ControlEvent, view: ft.Control) -> None:
        view_container.content = view
        page.update()

    # AppBar buttons
    input_button = ft.FilledButton(
        "Input",
        icon=ft.Icons.INPUT,
        on_click=lambda e: switch_view(e, input_view),
        tooltip="Input",
    )
    input_button.content_description = "Input"

    pcr_button = ft.FilledButton(
        "PCR",
        ref=pcr_button_ref,
        on_click=on_pcr_click,
        disabled=True,
        icon=ft.Icons.ANALYTICS,
        tooltip="PCR",
    )
    pcr_button.content_description = "PCR"

    dimers_button = ft.FilledButton(
        "Primer Dimers",
        ref=dimers_button_ref,
        on_click=on_dimers_click,
        disabled=True,
        icon=ft.Icons.COMPARE_ARROWS,
        tooltip="Primer Dimers",
    )
    dimers_button.content_description = "Primer Dimers"

    settings_button = ft.FilledButton(
        "Settings",
        icon=ft.Icons.SETTINGS,
        on_click=lambda e: switch_view(e, settings_view),
        tooltip="Settings",
    )
    settings_button.content_description = "Settings"

    save_btn_control = ft.FilledButton(
        "Save all",
        icon=ft.Icons.SAVE,
        tooltip="Save all",
        on_click=save_state,
    )
    save_btn_control.content_description = "Save all"

    load_btn_control = ft.FilledButton(
        "Load all",
        icon=ft.Icons.UPLOAD_FILE,
        tooltip="Load all",
        on_click=load_state,
    )
    load_btn_control.content_description = "Load all"

    page.appbar = ft.AppBar(
        title=ft.Row(
            [
                ft.Image(src="images/favicon.png", height=32, fit="contain"),
                ft.Text("AmplifyP", size=20, weight=ft.FontWeight.BOLD),
            ],
            alignment=ft.MainAxisAlignment.START,
            vertical_alignment=ft.CrossAxisAlignment.CENTER,
            spacing=8,
        ),
        elevation_on_scroll=0,
        actions=[
            input_button,
            ft.Container(width=16),
            pcr_button,
            ft.Container(width=16),
            dimers_button,
            ft.Container(width=16),
            settings_button,
            ft.Container(width=16),
            ft.VerticalDivider(),
            ft.Container(width=16),
            save_btn_control,
            ft.Container(width=16),
            load_btn_control,
            ft.Container(width=20),
        ],
    )

    page.add(
        ft.Divider(height=1, thickness=1),
        view_container,
    )
