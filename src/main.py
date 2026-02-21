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

from typing import cast

import flet as ft
import yaml

from amplifyp.gui.views import InputView, ResultView, SettingsView

STATE_FILE = "amplify_gui_state.yaml"
UPLOAD_DIR = "uploads"


def _is_pyodide() -> bool:
    """Check if running in a Pyodide environment."""
    try:
        import pyodide  # noqa: F401

        return True
    except ImportError:
        return False


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


async def _read_file_via_js() -> str | None:
    """Read a file using JavaScript FileReader API (for Pyodide web).

    Opens a browser file input dialog and reads the selected file's
    text content using the JavaScript FileReader API. Returns the file
    content as a string, or None if the user cancels.
    """
    import asyncio

    from pyodide.code import run_js

    loop = asyncio.get_event_loop()
    future: asyncio.Future[str | None] = loop.create_future()

    def resolve(content: object) -> None:
        if not future.done():
            future.set_result(str(content) if content else None)

    run_js(
        """
        (resolve_fn) => {
            const input = document.createElement('input');
            input.type = 'file';
            input.accept = '.yaml,.yml';
            input.onchange = (event) => {
                const file = event.target.files[0];
                if (!file) { resolve_fn(null); return; }
                const reader = new FileReader();
                reader.onload = (e) => resolve_fn(e.target.result);
                reader.onerror = () => resolve_fn(null);
                reader.readAsText(file);
            };
            input.click();
        }
    """
    )(resolve)

    return await future


async def _save_file_via_js(content: str, filename: str) -> None:
    """Save a text file via JavaScript Blob API (for Pyodide web).

    Creates a Blob from the content, generates an object URL,
    attaches it to a temporary anchor element with a download
    attribute, and clicks it to trigger a browser download.
    """
    from pyodide.code import run_js

    run_js(
        """
        (content, filename) => {
            const blob = new Blob([content], {type: 'text/yaml'});
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }
    """
    )(content, filename)


async def _load_via_upload(
    page: ft.Page,
    show_snackbar: object,
) -> dict[str, object] | None:
    """Load state via Flet upload (for server-side web).

    Picks a file, uploads it to the server's upload directory,
    then reads and parses the YAML content. Returns the parsed
    state dict, or None if the user cancels or an error occurs.
    """
    import asyncio
    import os
    from collections.abc import Callable

    _show = cast(Callable[[str], None], show_snackbar)

    upload_complete = asyncio.Event()
    upload_error: list[str | None] = [None]

    def on_upload(e: ft.FilePickerUploadEvent) -> None:
        if e.error:
            upload_error[0] = e.error
            upload_complete.set()
        elif e.progress == 1.0:
            upload_complete.set()

    file_picker = ft.FilePicker(on_upload=on_upload)
    files = await file_picker.pick_files(
        dialog_title="Load State",
        allowed_extensions=["yaml", "yml"],
    )
    if files is None or len(files) == 0:
        return None

    file = files[0]
    upload_name = file.name
    upload_url = page.get_upload_url(upload_name, 60)
    await file_picker.upload(
        [
            ft.FilePickerUploadFile(
                name=upload_name,
                upload_url=upload_url,
            )
        ]
    )
    await upload_complete.wait()

    if upload_error[0]:
        _show(f"Upload error: {upload_error[0]}")
        return None

    # Flet resolves upload_dir relative to the script directory,
    # so we must do the same to find the uploaded file.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, UPLOAD_DIR, upload_name)
    with open(file_path) as f:
        return yaml.safe_load(f)  # type: ignore[no-any-return]


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
        state = input_view.get_state()
        state["settings"] = settings_view.get_state()
        yaml_str = _serialize_state(state)

        if page.web and _is_pyodide():
            # Pyodide: trigger download via JS Blob API.
            await _save_file_via_js(yaml_str, STATE_FILE)
            show_snackbar("State saved successfully!")
        elif page.web:
            # Server-side web: write to assets directory
            # (served as static files) for download.
            import os

            script_dir = os.path.dirname(os.path.abspath(__file__))
            assets_dir = os.path.join(script_dir, "assets")
            os.makedirs(assets_dir, exist_ok=True)
            save_path = os.path.join(assets_dir, STATE_FILE)
            with open(save_path, "w") as f:
                f.write(yaml_str)
            await ft.UrlLauncher().launch_url(f"/{STATE_FILE}")
            show_snackbar("State saved successfully!")
        else:
            file_path = await ft.FilePicker().save_file(
                dialog_title="Save State",
                file_name=STATE_FILE,
                allowed_extensions=["yaml", "yml"],
            )
            if file_path is None:
                return
            try:
                with open(file_path, "w") as f:
                    f.write(yaml_str)
                show_snackbar("State saved successfully!")
            except Exception as ex:
                show_snackbar(f"Error saving state: {ex}")

    async def load_state(e: ft.ControlEvent) -> None:
        try:
            if page.web and _is_pyodide():
                # Pyodide static site: use JS FileReader
                # since upload and file.path are unavailable.
                content = await _read_file_via_js()
                if content is None:
                    return
                state = yaml.safe_load(content)
            elif page.web:
                # Server-side web (flet run --web): pick
                # file then upload to server temp storage.
                state = await _load_via_upload(page, show_snackbar)
                if state is None:
                    return
            else:
                files = await ft.FilePicker().pick_files(
                    dialog_title="Load State",
                    allowed_extensions=["yaml", "yml"],
                )
                if files is None or len(files) == 0:
                    return
                file_path_str = files[0].path
                if not file_path_str:
                    show_snackbar("Error: Could not locate the file.")
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
                ft.Icons.UPLOAD_FILE,
                tooltip="Load State",
                on_click=load_state,
            ),
        ],
    )

    page.add(view_container)


if __name__ == "__main__":  # pragma: no cover
    import os
    import secrets

    os.environ.setdefault("FLET_SECRET_KEY", secrets.token_hex(32))
    ft.run(main, upload_dir=UPLOAD_DIR)
