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

DIMERS_LABEL = "Primer Dimers"

STATE_FILE = "amplify_gui_state.yaml"


def main(page: ft.Page) -> None:
    """Main entry point for the Flet application."""
    page.title = "AmplifyP"
    page.vertical_alignment = ft.MainAxisAlignment.START
    page.fonts = {"Roboto Mono": "fonts/RobotoMono-Regular.ttf"}

    # Centralize state storage
    input_data = GUIInput()
    settings = GUISettings()
    has_run_pcr = False
    results_outdated = False
    results_button_ref = ft.Ref[ft.FilledButton]()
    dimers_button_ref = ft.Ref[ft.FilledButton]()

    def on_results_click(e: ft.ControlEvent) -> None:
        """Handle results button click: switch view and run PCR."""
        nonlocal has_run_pcr, results_outdated
        switch_view(e, result_view)
        result_view.run_pcr()
        has_run_pcr = True
        results_outdated = False
        if results_button_ref.current:
            results_button_ref.current.text = "Results"
        page.update()

    def on_dimers_click(e: ft.ControlEvent) -> None:
        """Handle dimers button click: switch view and run analysis."""
        switch_view(e, dimers_view)
        dimers_view.run_analysis()
        page.update()

    results_button = ft.FilledButton(
        "Results",
        ref=results_button_ref,
        on_click=on_results_click,
        disabled=True,
        icon=ft.Icons.ANALYTICS,
    )

    dimers_button = ft.FilledButton(
        DIMERS_LABEL,
        ref=dimers_button_ref,
        on_click=on_dimers_click,
        disabled=True,
        icon=ft.Icons.COMPARE_ARROWS,
        tooltip=DIMERS_LABEL,
    )
    dimers_button.content_description = DIMERS_LABEL

    def update_results_button_state() -> None:
        """Enable results and dimers buttons only if input is valid."""
        nonlocal has_run_pcr, results_outdated
        input_view.sync_to_state()
        has_template = bool(input_data.template.strip())
        has_primers = len(input_data.get_active_primers()) > 0
        is_enabled = has_template and has_primers

        btn = results_button_ref.current
        if btn:
            btn.disabled = not is_enabled

            # Set outdated if enabled and inputs change AFTER a first run
            if is_enabled and has_run_pcr:
                results_outdated = True

            label = (
                "Results *" if (results_outdated and is_enabled) else "Results"
            )
            btn.text = label

        dimers_btn = dimers_button_ref.current
        if dimers_btn:
            dimers_btn.disabled = not has_primers

        page.update()

    def run_analysis_in_background() -> None:
        nonlocal has_run_pcr, results_outdated
        result_view.run_pcr()
        dimers_view.run_analysis()
        has_run_pcr = True
        results_outdated = False
        if results_button_ref.current:
            results_button_ref.current.text = "Results"
        page.update()

    def run_results_tab(e: ft.ControlEvent) -> None:
        switch_view(e, result_view)
        run_analysis_in_background()

    def run_apply_settings(e: ft.ControlEvent) -> None:
        run_analysis_in_background()

    input_view = InputView(
        page,
        input_data,
        settings,
        on_change=lambda e: update_results_button_state(),
        on_stop_editing=update_results_button_state,
    )
    settings_view = SettingsView(
        page,
        settings,
        on_change=lambda e: update_results_button_state(),
        on_apply=run_apply_settings,
        on_reset=run_apply_settings,
    )
    result_view = ResultView(page, input_data, settings)
    dimers_view = PrimerDimerView(page, input_data, settings)

    # Save and Load State
    def show_snackbar(message: str) -> None:
        page.overlay.append(ft.SnackBar(ft.Text(message), open=True))
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
                dialog_title="Save State",
                file_name=STATE_FILE,
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
                dialog_title="Load State",
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
            input_view.update_ui()
            settings_view.update_ui()
            update_results_button_state()
            show_snackbar("State loaded successfully!")
        except Exception as ex:
            import traceback

            tb = traceback.format_exc()
            print("LOAD STATE ERROR:", tb)
            show_snackbar(f"Error loading state: {ex}")

    view_container = ft.Container(content=input_view, expand=True)

    def switch_view(e: ft.ControlEvent, view: ft.Control) -> None:
        view_container.content = view
        page.update()

    save_btn_control = ft.FilledButton(
        "Save State",
        icon=ft.Icons.SAVE,
        tooltip="Save State",
        on_click=save_state,
    )
    save_btn_control.content_description = "Save State"

    load_btn_control = ft.FilledButton(
        "Load State",
        icon=ft.Icons.UPLOAD_FILE,
        tooltip="Load State",
        on_click=load_state,
    )
    load_btn_control.content_description = "Load State"

    page.appbar = ft.AppBar(
        title=ft.Text("AmplifyP"),
        elevation_on_scroll=0,
        actions=[
            ft.FilledButton(
                "Input",
                icon=ft.Icons.INPUT,
                on_click=lambda e: switch_view(e, input_view),
            ),
            ft.Container(width=16),
            results_button,
            ft.Container(width=16),
            dimers_button,
            ft.Container(width=16),
            ft.FilledButton(
                "Settings",
                icon=ft.Icons.SETTINGS,
                on_click=lambda e: switch_view(e, settings_view),
            ),
            ft.Container(width=16),
            ft.VerticalDivider(),
            ft.Container(width=16),
            save_btn_control,
            ft.Container(width=16),
            load_btn_control,
            ft.Container(width=20),
        ],
    )

    page.add(view_container)
