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

from amplifyp.gui.views import InputView, ResultView, SettingsView

STATE_FILE = "amplify_gui_state.yaml"


def _serialize_state(state: dict[str, object]) -> str:
    """Serialize state dict to YAML string."""

    def multiline_presenter(dumper: yaml.Dumper, data: str) -> yaml.ScalarNode:
        if "\n" in data:
            return dumper.represent_scalar(
                "tag:yaml.org,2002:str", data, style="|"
            )
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    yaml.add_representer(str, multiline_presenter)
    return yaml.dump(state, sort_keys=False)


def main(page: ft.Page) -> None:
    """Main entry point for the Flet application."""
    page.title = "AmplifyP"
    page.vertical_alignment = ft.MainAxisAlignment.START
    results_outdated = [False]
    has_run_pcr = [False]
    results_button_ref = ft.Ref[ft.TextButton]()

    def on_results_click(e: ft.ControlEvent) -> None:
        """Handle results button click: switch view and run PCR."""
        switch_view(e, result_view)
        result_view.run_pcr()
        has_run_pcr[0] = True
        results_outdated[0] = False
        if results_button_ref.current:
            results_button_ref.current.content = ft.Text("Results")
        page.update()

    results_button = ft.TextButton(
        "Results",
        ref=results_button_ref,
        on_click=on_results_click,
        disabled=True,
    )

    def update_results_button_state() -> None:
        """Enable results button only if input is valid."""
        has_template = bool(input_view.get_template().strip())
        has_primers = len(input_view.get_primers()) > 0
        is_enabled = has_template and has_primers

        btn = results_button_ref.current
        if not btn:
            return

        btn.disabled = not is_enabled

        # Set outdated if enabled and inputs change AFTER a first run
        if is_enabled and has_run_pcr[0]:
            results_outdated[0] = True

        label = (
            "Results *" if (results_outdated[0] and is_enabled) else "Results"
        )
        btn.content = ft.Text(label)
        page.update()

    input_view = InputView(
        page, on_change=lambda e: update_results_button_state()
    )
    settings_view = SettingsView(page)
    result_view = ResultView(page, input_view, settings_view)

    # Save and Load State
    def show_snackbar(message: str) -> None:
        page.overlay.append(ft.SnackBar(ft.Text(message), open=True))
        page.update()

    async def save_state(e: ft.ControlEvent) -> None:
        state = input_view.get_state()
        state["settings"] = settings_view.get_state()
        yaml_str = _serialize_state(state)

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

            input_view.set_state(parsed_state)
            settings_view.set_state(parsed_state)
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

    save_btn_control = ft.IconButton(
        ft.Icons.SAVE,
        tooltip="Save State",
        on_click=save_state,
    )
    save_btn_control.content_description = "Save State"

    load_btn_control = ft.IconButton(
        ft.Icons.UPLOAD_FILE,
        tooltip="Load State",
        on_click=load_state,
    )
    load_btn_control.content_description = "Load State"

    page.appbar = ft.AppBar(
        title=ft.Text("AmplifyP"),
        actions=[
            ft.TextButton(
                "Input", on_click=lambda e: switch_view(e, input_view)
            ),
            results_button,
            ft.TextButton(
                "Settings", on_click=lambda e: switch_view(e, settings_view)
            ),
            ft.VerticalDivider(),
            save_btn_control,
            load_btn_control,
        ],
    )

    page.add(view_container)
