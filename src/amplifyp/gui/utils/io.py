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

"""File input/output and dialogue helpers."""

import asyncio
from collections.abc import Callable

import flet as ft


def _read_file(path: str) -> str:
    """Synchronous file read helper for use with asyncio.to_thread."""
    with open(path, encoding="utf-8") as f:
        return f.read()


def _write_file(path: str, content: str) -> None:
    """Synchronous file write helper for use with asyncio.to_thread."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)


async def pick_and_read_file(
    page: ft.Page,
    dialog_title: str,
    allowed_extensions: list[str],
    show_notification: Callable[[str], None],
) -> str | None:
    """Open a file picker to load a file, and read its text content."""
    file_picker = ft.FilePicker()
    page.services.append(file_picker)
    page.update()
    try:
        files = await file_picker.pick_files(
            dialog_title=dialog_title,
            allowed_extensions=allowed_extensions,
            file_type=ft.FilePickerFileType.CUSTOM,
            with_data=True,
        )
        if not files:
            return None

        file = files[0]
        if file.bytes is not None:
            return file.bytes.decode("utf-8")  # type: ignore[no-any-return]
        else:
            if not file.path:
                show_notification("Error: Could not read file content.")
                return None
            content = await asyncio.to_thread(_read_file, file.path)
            return content
    except (OSError, ValueError) as ex:
        show_notification(f"Error loading file: {ex}")
        return None
    finally:
        if file_picker in page.services:
            page.services.remove(file_picker)
            page.update()


async def save_and_write_file(
    page: ft.Page,
    dialog_title: str,
    file_name: str,
    allowed_extensions: list[str],
    content: str,
    show_notification: Callable[[str], None],
    success_message_desktop: str = "Saved successfully!",
    success_message_web: str = "Ready for download!",
) -> bool:
    """Save content using the file picker, supporting both Web and Desktop."""
    file_picker = ft.FilePicker()
    page.services.append(file_picker)
    page.update()
    try:
        file_path = await file_picker.save_file(
            dialog_title=dialog_title,
            file_name=file_name,
            allowed_extensions=allowed_extensions,
            file_type=ft.FilePickerFileType.CUSTOM,
            src_bytes=content.encode("utf-8"),
        )
        if page.web:
            show_notification(success_message_web)
            return True
        else:
            if file_path is None:
                return False
            await asyncio.to_thread(_write_file, file_path, content)
            show_notification(success_message_desktop)
            return True
    except OSError as ex:
        show_notification(f"Error saving file: {ex}")
        return False
    finally:
        if file_picker in page.services:
            page.services.remove(file_picker)
            page.update()
