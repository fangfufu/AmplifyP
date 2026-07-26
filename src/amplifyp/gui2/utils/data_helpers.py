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

"""Consolidated data, file I/O, state persistence, and sequence utilities."""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

import yaml
from PySide6.QtWidgets import QFileDialog, QMessageBox


def serialise_state(state: dict[str, object]) -> str:
    """Serialise state dict to YAML string, handling multiline strings."""

    def multiline_presenter(dumper: yaml.Dumper, data: str) -> yaml.ScalarNode:
        if "\n" in data:
            return dumper.represent_scalar(
                "tag:yaml.org,2002:str", data, style="|"
            )
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    class _StateDumper(yaml.Dumper):
        pass

    _StateDumper.add_representer(str, multiline_presenter)
    return yaml.dump(state, Dumper=_StateDumper, sort_keys=False)


logger = logging.getLogger(__name__)


# ==============================================================================
# Sequence Handling Utilities
# ==============================================================================


def clean_sequence(seq: str) -> str:
    """Clean sequence of escaped and standard whitespaces."""
    if not seq:
        return ""
    clean = str(seq).replace("\\n", "").replace("\\t", "").replace("\\r", "")
    return "".join(clean.split())


_SYSTEM_MONOSPACE_FONTS = (
    "monospace",
    "Courier New",
    "Courier",
    "Consolas",
    "Lucida Console",
    "DejaVu Sans Mono",
    "Ubuntu Mono",
    "Liberation Mono",
)


def _resolve_font_family(font_family: str) -> str:
    """Return *font_family* if it is a known system font, else "Roboto Mono"."""
    if not font_family:
        return "Roboto Mono"
    ff = font_family.strip()
    if ff.lower() in {f.lower() for f in _SYSTEM_MONOSPACE_FONTS}:
        return ff
    if ff.lower() == "roboto mono":
        return "Roboto Mono"
    return ff


def format_sequence(seq: str, wrap_length: int = 80) -> str:
    """Format sequence into lines of specified length."""
    clean = clean_sequence(seq)
    return "\n".join(
        [clean[i : i + wrap_length] for i in range(0, len(clean), wrap_length)]
    )


# ==============================================================================
# File I/O Helpers
# ==============================================================================


def _read_file(path: str) -> str:
    """Synchronous file read helper."""
    with open(path, encoding="utf-8") as f:
        return f.read()


def _write_file(path: str, content: str) -> None:
    """Synchronous file write helper."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)


def pick_and_read_file(
    parent_widget: Any,
    dialog_title: str,
    allowed_extensions: list[str],
    show_notification: Callable[[str], None],
) -> str | None:
    """Open a file dialog to load a file, and read its text content."""
    if allowed_extensions:
        exts = " ".join(f"*.{ext}" for ext in allowed_extensions)
        filter_str = f"{dialog_title} ({exts})"
    else:
        filter_str = "All Files (*)"

    file_path, _ = QFileDialog.getOpenFileName(
        parent_widget, dialog_title, "", filter_str
    )
    if not file_path:
        return None

    try:
        with open(file_path, encoding="utf-8") as f:
            return f.read()
    except (OSError, ValueError) as ex:
        show_notification(f"Error loading file: {ex}")
        return None


def save_and_write_file(
    parent_widget: Any,
    dialog_title: str,
    file_name: str,
    allowed_extensions: list[str],
    content: str,
    show_notification: Callable[[str], None],
) -> bool:
    """Save content using a file dialog."""
    if allowed_extensions:
        exts = " ".join(f"*.{ext}" for ext in allowed_extensions)
        filter_str = f"{dialog_title} ({exts})"
    else:
        filter_str = f"{dialog_title} (*.*)"

    file_path, _ = QFileDialog.getSaveFileName(
        parent_widget,
        dialog_title,
        file_name,
        filter_str,
    )
    if not file_path:
        return False

    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        show_notification("Saved successfully!")
        return True
    except OSError as ex:
        show_notification(f"Error saving file: {ex}")
        return False


# ==============================================================================
# State Persistence & Restoration
# ==============================================================================


def get_last_state_path(settings: Any) -> Path:
    """Get the OS-specific path for the last saved GUI state."""
    config_path = settings._get_config_path()
    return config_path.parent / "last_state.yaml"


def save_last_state(controller: Any) -> None:
    """Save the last template and primers to local storage."""
    if not controller.settings.get("auto_reload_on_startup", True):
        return

    controller.input_view.sync_to_state()
    state_dict = {
        "input": controller.input_data.to_dict(),
    }

    path = get_last_state_path(controller.settings)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        yaml_str = serialise_state(state_dict)
        with open(path, "w", encoding="utf-8") as f:
            f.write(yaml_str)
    except Exception as e:
        logger.exception("Error saving last state to %s: %s", path, e)


def load_last_state(controller: Any) -> None:
    """Load the last template and primers from local storage."""
    if not controller.settings.get("auto_reload_on_startup", True):
        return

    path = get_last_state_path(controller.settings)
    state_dict = None

    if path.exists():
        try:
            with open(path, encoding="utf-8") as f:
                content = f.read()
            state_dict = yaml.safe_load(content)
        except Exception as e:
            logger.exception("Error loading last state from %s: %s", path, e)

    if state_dict and isinstance(state_dict, dict):
        if "input" not in state_dict:
            state_dict = {"input": state_dict}
        apply_parsed_state(controller, state_dict, ignore_settings=True)


def restore_state_from_file(controller: Any, path: str) -> None:
    """Restore app state from a YAML file on startup."""
    try:
        with open(path, encoding="utf-8") as f:
            content = f.read()

        parsed_state = yaml.safe_load(content)

        if not isinstance(parsed_state, dict):
            logger.warning("Invalid state file format, ignoring.")
            return

        apply_parsed_state(controller, parsed_state)
        logger.info("State loaded successfully from %s", path)
    except (OSError, ValueError, yaml.YAMLError):
        logger.exception("Error loading state file '%s'", path)


def apply_parsed_state(
    controller: Any, parsed_state: dict[str, Any], ignore_settings: bool = False
) -> None:
    """Apply parsed YAML state to the application."""
    if "input" in parsed_state:
        controller.input_data.from_dict(parsed_state["input"])
    else:
        controller.input_data.from_dict(parsed_state)
    if not ignore_settings and "settings" in parsed_state:
        controller.settings.from_dict(parsed_state["settings"])
        controller.settings.save_to_local()
    controller.apply_theme()
    controller.input_view.update_ui()
    controller.settings_view.update_ui()
    controller.update_pcr_button_state(update_page=False)
    if hasattr(controller, "main_window"):
        controller.main_window.update()


def save_state(
    controller: Any, show_notification: Callable[[str], None]
) -> bool:
    """Save app state to YAML configuration file."""
    if getattr(controller, "filepicker_open", False):
        return False
    controller.filepicker_open = True
    try:
        controller.input_view.sync_to_state()
        combined: dict[str, object] = {
            "input": controller.input_data.to_dict(),
        }
        yaml_str = serialise_state(combined)

        return save_and_write_file(
            parent_widget=controller.main_window,
            dialog_title="Save all",
            file_name="amplify_gui_state.yaml",
            allowed_extensions=["yaml", "yml"],
            content=yaml_str,
            show_notification=show_notification,
        )
    except (OSError, ValueError) as ex:
        show_notification(f"Error saving state: {ex}")
        return False
    finally:
        controller.filepicker_open = False


def load_state(
    controller: Any, show_notification: Callable[[str], None]
) -> bool:
    """Load app state from YAML configuration file."""
    if getattr(controller, "filepicker_open", False):
        return False
    controller.filepicker_open = True
    try:
        content = pick_and_read_file(
            parent_widget=controller.main_window,
            dialog_title="Load all",
            allowed_extensions=["yaml", "yml"],
            show_notification=show_notification,
        )
        if content is None:
            return False

        parsed_state = yaml.safe_load(content)

        if not isinstance(parsed_state, dict):
            show_notification("Error: Invalid state file format.")
            return False

        apply_parsed_state(controller, parsed_state, ignore_settings=True)
        show_notification("State loaded successfully!")
        return True
    except (OSError, ValueError, yaml.YAMLError) as ex:
        logger.exception("Error loading state:")
        show_notification(f"Error loading state: {ex}")
        return False
    finally:
        controller.filepicker_open = False


def confirm_clear(controller: Any) -> None:
    """Clear all inputs and reset state."""
    controller.input_data.template = ""
    controller.input_data.template_circular = False
    controller.input_data.primers = [{"name": "", "seq": "", "active": False}]
    if hasattr(controller.input_view, "primer_input"):
        if hasattr(controller.input_view.primer_input, "focused_primer_index"):
            controller.input_view.primer_input.focused_primer_index = None
        if hasattr(controller.input_view.primer_input, "selected_indices"):
            controller.input_view.primer_input.selected_indices.clear()
    controller.input_view.update_ui()
    controller.update_pcr_button_state(sync=False, update_page=False)
    save_last_state(controller)


def clear_all(controller: Any) -> bool:
    """Show a confirmation dialogue before clearing inputs."""
    reply = QMessageBox.question(
        controller.main_window,
        "Confirm Clear",
        "Are you sure you want to clear all template sequences\nand primers?",
        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        QMessageBox.StandardButton.No,
    )
    if reply == QMessageBox.StandardButton.Yes:
        confirm_clear(controller)
        return True
    return False
