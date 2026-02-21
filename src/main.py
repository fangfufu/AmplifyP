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

"""Main Flet application entry point."""

import flet as ft
import yaml

from amplifyp.gui.views import InputView, ResultView, SettingsView

STATE_FILE = "amplify_gui_state.yaml"


def main(page: ft.Page) -> None:
    """Main entry point for the Flet application."""
    page.title = "AmplifyP"
    page.vertical_alignment = ft.MainAxisAlignment.START

    input_view = InputView(page)
    settings_view = SettingsView(page)
    result_view = ResultView(page, input_view, settings_view)

    # Save and Load State
    def show_snackbar(message: str) -> None:
        page.overlay.append(ft.SnackBar(ft.Text(message), open=True))
        page.update()

    async def save_state(e: ft.ControlEvent) -> None:
        file_path = await ft.FilePicker().save_file(
            dialog_title="Save State",
            file_name="amplify_gui_state.yaml",
            allowed_extensions=["yaml", "yml"],
        )
        if file_path is None:
            return

        state = input_view.get_state()
        state["settings"] = settings_view.get_state()

        def multiline_presenter(
            dumper: yaml.Dumper, data: str
        ) -> yaml.ScalarNode:
            if "\n" in data:
                return dumper.represent_scalar(
                    "tag:yaml.org,2002:str", data, style="|"
                )
            return dumper.represent_scalar("tag:yaml.org,2002:str", data)

        yaml.add_representer(str, multiline_presenter)

        try:
            with open(file_path, "w") as f:
                yaml.dump(state, f, sort_keys=False)
            show_snackbar("State saved successfully!")
        except Exception as ex:
            show_snackbar(f"Error saving state: {ex}")

    async def load_state(e: ft.ControlEvent) -> None:
        import asyncio
        import os

        upload_event = asyncio.Event()
        upload_error = [None]

        def on_upload(e: ft.FilePickerUploadEvent) -> None:
            if e.error:
                upload_error[0] = e.error
                upload_event.set()
            elif e.progress == 1.0:
                upload_event.set()

        file_picker = ft.FilePicker(on_upload=on_upload)
        page.overlay.append(file_picker)
        page.update()

        files = await file_picker.pick_files(
            dialog_title="Load State",
            allowed_extensions=["yaml", "yml"],
        )
        if files is None or len(files) == 0:
            page.overlay.remove(file_picker)
            page.update()
            return

        try:
            file = files[0]
            file_path_str = file.path

            if page.web and not file_path_str:
                upload_name = file.name
                upload_url = page.get_upload_url(upload_name, 60)
                await file_picker.upload(
                    [
                        ft.FilePickerUploadFile(
                            upload_name, upload_url=upload_url
                        )
                    ]
                )
                await upload_event.wait()

                if upload_error[0]:
                    show_snackbar(f"Error uploading file: {upload_error[0]}")
                    return

                possible_paths = [
                    upload_name,
                    f"uploads/{upload_name}",
                    f"/tmp/{upload_name}",  # noqa: S108
                    f"assets/uploads/{upload_name}",
                ]
                file_path_str = next(
                    (p for p in possible_paths if os.path.exists(p)), None
                )

            if not file_path_str:
                show_snackbar("Error: Could not locate the selected file.")
                return

            with open(file_path_str) as f:
                state = yaml.safe_load(f)

            input_view.set_state(state)
            settings_view.set_state(state)

            show_snackbar("State loaded successfully!")
        except FileNotFoundError:
            show_snackbar("No saved state found.")
        except Exception as ex:
            show_snackbar(f"Error loading state: {ex}")
        finally:
            page.overlay.remove(file_picker)
            page.update()

    view_container = ft.Container(content=input_view, expand=True)

    def switch_view(e: ft.ControlEvent, view: ft.Control) -> None:
        view_container.content = view
        page.update()

    page.appbar = ft.AppBar(
        title=ft.Text("AmplifyP"),
        actions=[
            ft.TextButton(
                "Input", on_click=lambda e: switch_view(e, input_view)
            ),
            ft.TextButton(
                "Results", on_click=lambda e: switch_view(e, result_view)
            ),
            ft.TextButton(
                "Settings", on_click=lambda e: switch_view(e, settings_view)
            ),
            ft.VerticalDivider(),
            ft.IconButton(
                ft.Icons.SAVE, tooltip="Save State", on_click=save_state
            ),
            ft.IconButton(
                ft.Icons.UPLOAD_FILE, tooltip="Load State", on_click=load_state
            ),
        ],
    )

    page.add(view_container)


if __name__ == "__main__":  # pragma: no cover
    ft.run(main)
