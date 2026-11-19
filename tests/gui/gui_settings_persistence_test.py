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

"""Tests for GUI settings persistence on Desktop and Web platforms."""

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, mock_open, patch

import flet as ft
import pytest
import yaml

from amplifyp.gui.settings import GUISettings


def test_get_config_path() -> None:
    """Test that OS-specific config paths are resolved correctly."""
    settings = GUISettings()

    # Test Windows
    with (
        patch("sys.platform", "win32"),
        patch.dict(
            "os.environ",
            {"APPDATA": "C:\\Users\\test\\AppData\\Roaming"},
            clear=True,
        ),
    ):
        path = settings._get_config_path()
        assert "AppData" in str(path)
        assert "AmplifyP" in str(path)
        assert path.name == "settings.yaml"

    # Test macOS
    with (
        patch("sys.platform", "darwin"),
        patch.dict("os.environ", {"HOME": "/Users/test"}, clear=True),
    ):
        path = settings._get_config_path()
        assert (
            path.as_posix()
            == "/Users/test/Library/Application Support/AmplifyP/settings.yaml"
        )

    # Test Linux with XDG_CONFIG_HOME
    with (
        patch("sys.platform", "linux"),
        patch.dict(
            "os.environ", {"XDG_CONFIG_HOME": "/custom/config"}, clear=True
        ),
    ):
        path = settings._get_config_path()
        assert path.as_posix() == "/custom/config/amplifyp/settings.yaml"

    # Test Linux without XDG_CONFIG_HOME
    with (
        patch("sys.platform", "linux"),
        patch.dict("os.environ", {"HOME": "/home/test"}, clear=True),
    ):
        path = settings._get_config_path()
        assert path.as_posix() == "/home/test/.config/amplifyp/settings.yaml"


def test_desktop_save_to_local() -> None:
    """Test saving settings to local YAML file on desktop."""
    settings = GUISettings()
    settings["primability_cutoff"] = "0.75"

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False

    mock_path = MagicMock(spec=Path)
    mock_path.parent = MagicMock()

    m_open = mock_open()
    with (
        patch.object(settings, "_get_config_path", return_value=mock_path),
        patch("builtins.open", m_open),
    ):
        settings.save_to_local(mock_page)

    mock_path.parent.mkdir.assert_called_once_with(parents=True, exist_ok=True)

    # Filter out external background file opens (e.g. xonsh history json)
    config_open_call = None
    for call in m_open.call_args_list:
        if call[0] and call[0][0] == mock_path:
            config_open_call = call
            break
    assert config_open_call is not None, (
        f"Config file path {mock_path} was not opened"
    )
    assert config_open_call[0][1] == "w"
    assert config_open_call[1].get("encoding") == "utf-8"

    # Check if YAML was dumped (contains primability_cutoff)
    written_data = "".join(
        call[0][0] for call in m_open.return_value.write.call_args_list
    )
    parsed = yaml.safe_load(written_data)
    assert parsed["primability_cutoff"] == "0.75"


def test_desktop_load_from_local() -> None:
    """Test loading settings from local YAML file on desktop."""
    settings = GUISettings()
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False

    mock_path = MagicMock(spec=Path)
    mock_path.exists.return_value = True

    yaml_content = "primability_cutoff: '0.65'\ndark_mode: true"
    m_open = mock_open(read_data=yaml_content)

    with (
        patch.object(settings, "_get_config_path", return_value=mock_path),
        patch("builtins.open", m_open),
    ):
        settings.load_from_local(mock_page)

    assert settings["primability_cutoff"] == "0.65"
    assert settings["dark_mode"] is True


def test_desktop_load_from_local_corrupted() -> None:
    """Test loading from a corrupted YAML file falls back to defaults."""
    settings = GUISettings()
    settings["primability_cutoff"] = "0.99"

    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False

    mock_path = MagicMock(spec=Path)
    mock_path.exists.return_value = True

    # Invalid YAML syntax
    m_open = mock_open(read_data="{invalid: yaml:")

    with (
        patch.object(settings, "_get_config_path", return_value=mock_path),
        patch("builtins.open", m_open),
    ):
        # Should not raise exception
        settings.load_from_local(mock_page)

    # Since it failed, the settings should still be whatever they were
    assert settings["primability_cutoff"] == "0.99"


def test_web_save_and_load_client_storage() -> None:
    """Test saving and loading settings on web using page.client_storage."""
    settings = GUISettings()
    settings["primability_cutoff"] = "0.85"
    settings["dark_mode"] = True

    import asyncio

    mock_page = MagicMock()
    mock_page.web = True
    mock_page.run_task = MagicMock(
        side_effect=lambda fn, *args, **kwargs: asyncio.run(fn(*args, **kwargs))
    )

    # Set up SharedPreferences mock
    mock_storage: dict[str, Any] = {}

    async def mock_set(key: str, value: Any) -> bool:
        mock_storage[key] = value
        return True

    async def mock_get(key: str) -> Any:
        return mock_storage.get(key)

    async def mock_contains(key: str) -> bool:
        return key in mock_storage

    mock_prefs = MagicMock()
    mock_prefs.set = MagicMock(side_effect=mock_set)
    mock_prefs.get = MagicMock(side_effect=mock_get)
    mock_prefs.contains_key = MagicMock(side_effect=mock_contains)

    with patch("flet.SharedPreferences", return_value=mock_prefs):
        # Test Save
        settings.save_to_local(mock_page)

        # Assert key exists in storage with prefix
        assert "amplifyp.settings.primability_cutoff" in mock_storage
        assert mock_storage["amplifyp.settings.primability_cutoff"] == "0.85"
        assert mock_storage["amplifyp.settings.dark_mode"] is True

        # Test Load
        new_settings = GUISettings()
        new_settings.load_from_local(mock_page)

        assert new_settings["primability_cutoff"] == "0.85"
        assert new_settings["dark_mode"] is True


def test_dark_mode_system_coercion() -> None:
    """Test that setting dark_mode to 'system' is coerced correctly."""
    settings = GUISettings()
    settings["dark_mode"] = False
    assert settings.get("dark_mode") is False
    settings["dark_mode"] = "system"
    assert settings.get("dark_mode") == "system"


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_web_save_snapshot_and_key_filtering() -> None:
    """Test that web saves snapshot values at call time and filter keys."""
    import asyncio

    settings = GUISettings()
    settings["primer_info_panel_height"] = 420

    mock_page = MagicMock()
    mock_page.web = True
    tasks: list[asyncio.Task[Any]] = []

    def _run_task(fn: Any, *args: Any, **kwargs: Any) -> asyncio.Task[Any]:
        t = asyncio.create_task(fn(*args, **kwargs))
        tasks.append(t)
        return t

    mock_page.run_task = MagicMock(side_effect=_run_task)

    mock_storage: dict[str, Any] = {}

    async def slow_set(key: str, value: Any) -> bool:
        await asyncio.sleep(0.01)
        mock_storage[key] = value
        return True

    mock_prefs = MagicMock()
    mock_prefs.set = MagicMock(side_effect=slow_set)

    with patch("flet.SharedPreferences", return_value=mock_prefs):
        t1 = settings.save_to_local(
            mock_page, keys=("primer_info_panel_height",)
        )
        # Mutate live settings after scheduling; snapshot must keep 420
        settings["primer_info_panel_height"] = 480
        await t1
        assert mock_storage == {
            "amplifyp.settings.primer_info_panel_height": 420
        }


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_web_startup_awaits_preferences_before_init() -> None:
    """Test web startup awaits preferences before logging and controller."""
    import asyncio

    from amplifyp.gui.app import main as app_main

    mock_page = MagicMock()
    mock_page.web = True
    tasks: list[asyncio.Task[Any]] = []

    def _run_task(fn: Any, *a: Any, **kw: Any) -> asyncio.Task[Any]:
        t = asyncio.create_task(fn(*a, **kw))
        tasks.append(t)
        return t

    mock_page.run_task = MagicMock(side_effect=_run_task)

    events: list[str] = []

    async def delayed_contains(key: str) -> bool:
        await asyncio.sleep(0.01)
        return key == "amplifyp.settings.dark_mode"

    async def delayed_get(key: str) -> Any:
        return True

    mock_prefs = MagicMock()
    mock_prefs.contains_key = MagicMock(side_effect=delayed_contains)
    mock_prefs.get = MagicMock(side_effect=delayed_get)

    mock_ctrl = MagicMock()
    mock_ctrl.initialise = MagicMock(
        side_effect=lambda **kw: events.append(f"ctrl_init:{kw}")
    )

    with (
        patch("flet.SharedPreferences", return_value=mock_prefs),
        patch(
            "amplifyp.gui.app.initialise_logging",
            side_effect=lambda **kw: events.append("logging"),
        ),
        patch("amplifyp.gui.app.GUIController", return_value=mock_ctrl) as cls,
    ):
        app_main(mock_page)
        assert events == []
        await asyncio.gather(*tasks)
        assert events == ["logging", "ctrl_init:{'load_settings': False}"]
        passed_settings = cls.call_args.kwargs["settings"]
        assert passed_settings["dark_mode"] is True


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_web_startup_recovers_from_preference_error(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test web startup logs preference errors and continues with defaults."""
    import asyncio
    import logging

    from amplifyp.gui.app import main as app_main

    mock_page = MagicMock()
    mock_page.web = True
    tasks: list[asyncio.Task[Any]] = []

    def _run_task(fn: Any, *a: Any, **kw: Any) -> asyncio.Task[Any]:
        t = asyncio.create_task(fn(*a, **kw))
        tasks.append(t)
        return t

    mock_page.run_task = MagicMock(side_effect=_run_task)

    events: list[str] = []

    async def failing_contains(key: str) -> bool:
        raise RuntimeError("client storage unavailable")

    mock_prefs = MagicMock()
    mock_prefs.contains_key = MagicMock(side_effect=failing_contains)

    mock_ctrl = MagicMock()
    mock_ctrl.initialise = MagicMock(
        side_effect=lambda **kw: events.append(f"ctrl_init:{kw}")
    )

    with (
        caplog.at_level(logging.WARNING, logger="amplifyp.gui.app"),
        patch("flet.SharedPreferences", return_value=mock_prefs),
        patch(
            "amplifyp.gui.app.initialise_logging",
            side_effect=lambda **kw: events.append("logging"),
        ),
        patch("amplifyp.gui.app.GUIController", return_value=mock_ctrl) as cls,
    ):
        app_main(mock_page)
        await asyncio.gather(*tasks)
        assert events == ["logging", "ctrl_init:{'load_settings': False}"]
        passed_settings = cls.call_args.kwargs["settings"]
        assert passed_settings["dark_mode"] == "system"
        assert "Failed to load web preferences, using defaults" in caplog.text


def test_info_panel_overflow_height_fallback() -> None:
    """Test that huge saved primer_info_panel_height falls back to 260."""
    from amplifyp.gui.views.input.primer.info_panel import PrimerInfoPanel

    settings = GUISettings()
    settings._settings["primer_info_panel_height"] = 10**400
    panel = PrimerInfoPanel(settings, "monospace")
    panel._apply_height_mode(True)
    assert panel.height == 260.0


def test_handle_resize_clamps_temporary_drag_height() -> None:
    """Test handle_resize clamps temporary drag height in auto mode."""
    from amplifyp.gui.views.input.layout import handle_resize

    mock_view = MagicMock()
    mock_view.app_page.width = 800
    mock_view.app_page.height = 600
    mock_view.right_fraction = 0.5
    mock_view.settings.get = MagicMock(return_value=False)

    # Temporary drag height 700 should be clamped to 600 - 200 = 400
    mock_view.primer_input.primer_info_panel.height = 700.0
    handle_resize(mock_view, MagicMock())
    assert mock_view.primer_input.primer_info_panel.height == 400.0

    # Content-sized mode (height=None) should stay None
    mock_view.primer_input.primer_info_panel.height = None
    handle_resize(mock_view, MagicMock())
    assert mock_view.primer_input.primer_info_panel.height is None

    # Fixed-height mode applies height mode on resize
    mock_view.settings.get = MagicMock(return_value=True)
    handle_resize(mock_view, MagicMock())
    mock_view.primer_input.primer_info_panel._apply_height_mode.assert_called_once_with(
        True
    )


def test_gui_settings_edge_cases_and_coercions() -> None:
    """Test GUISettings coercions, views, error fallbacks, and storage."""
    from typing import Any, cast

    from amplifyp.dna import Primer
    from amplifyp.gui.colours import GUIColours

    # 1. Initialisation with settings_dict
    s1 = GUISettings(
        settings_dict={"colour_deficient": "yes", "dark_mode": "true"}
    )
    assert s1["colour_deficient"] == "yes"
    assert GUIColours.colour_deficient_mode is True
    assert GUIColours.dark_mode is True

    # 2. Type coercions in __setitem__
    s = GUISettings()
    # dark_mode non-string and string permutations
    s["dark_mode"] = True
    assert s["dark_mode"] is True
    s["dark_mode"] = "system"
    assert cast(Any, s["dark_mode"]) == "system"
    s["dark_mode"] = "true"
    assert s["dark_mode"] is True
    s["dark_mode"] = "false"
    assert s["dark_mode"] is False

    # bool setting
    s["auto_reload_on_startup"] = "yes"
    assert s["auto_reload_on_startup"] is True
    s["auto_reload_on_startup"] = False
    assert s["auto_reload_on_startup"] is False

    # int setting
    s._settings["log_rotation_max_bytes"] = 100
    s["log_rotation_max_bytes"] = "200"
    assert s["log_rotation_max_bytes"] == 200
    s["log_rotation_max_bytes"] = "invalid"
    assert s["log_rotation_max_bytes"] == "invalid"

    # float setting
    s._settings["primability_cutoff"] = 0.8
    s["primability_cutoff"] = "0.9"
    assert s["primability_cutoff"] == 0.9
    s["primability_cutoff"] = "invalid"
    assert s["primability_cutoff"] == "invalid"

    # colour_deficient
    s["colour_deficient"] = "true"
    assert s["colour_deficient"] is True
    s["colour_deficient"] = False
    assert s["colour_deficient"] is False

    # setting keys not yet in _settings
    s._settings.pop("colour_deficient", None)
    s["colour_deficient"] = "yes"
    assert GUIColours.colour_deficient_mode is True

    s._settings.pop("dark_mode", None)
    s["dark_mode"] = "yes"
    assert GUIColours.dark_mode is True

    # 3. Keys, items, iteration, membership
    assert "dark_mode" in s.keys()
    assert len(s.items()) > 0
    assert "dark_mode" in list(iter(s))
    assert "dark_mode" in s

    # 4. Factory methods
    _ = s.get_replication_settings()
    _ = s.get_primer_dimer_settings()
    _ = s.get_tm_settings()
    # verify cached tm_settings branch
    _ = s.get_tm_settings()

    # 5. _safe_float and _safe_int fallback handlers
    assert s._safe_float("missing_key", 9.9) == 9.9
    s._settings["primability_cutoff"] = None
    assert s._safe_float("primability_cutoff", 9.9) == 9.9
    s._settings["primability_cutoff"] = "invalid"
    assert s._safe_float("primability_cutoff", 9.9) == 9.9
    s._settings["primability_cutoff"] = ""
    assert s._safe_float("primability_cutoff", 9.9) == 9.9

    assert s._safe_int("missing_key", 5) == 5
    s._settings["log_rotation_max_bytes"] = None
    assert s._safe_int("log_rotation_max_bytes", 5) == 5
    s._settings["log_rotation_max_bytes"] = "invalid"
    assert s._safe_int("log_rotation_max_bytes", 5) == 5
    s._settings["log_rotation_max_bytes"] = ""
    assert s._safe_int("log_rotation_max_bytes", 5) == 5

    # 6. calculate_primer_tm with Lander method and SantaLucia fallback
    s["tm_method"] = "Lander / Amplify 4"
    tm_val = s.calculate_primer_tm(Primer("ACGTACGT"))
    assert tm_val > 0

    s["tm_method"] = "SantaLucia 1998 / Owczarzy 2008 (Default)"
    tm_val_sl = s.calculate_primer_tm(Primer("ACGTACGT"))
    assert tm_val_sl > 0

    # Degenerate sequence falls back to Lander
    tm_val_fallback = s.calculate_primer_tm(Primer("ACGTACGN"))
    assert tm_val_fallback >= 0

    # 7. from_dict type coercions
    s._settings["log_rotation_max_bytes"] = 100
    s._settings["primability_cutoff"] = 0.8
    s.from_dict(
        {
            "dark_mode": "true",
            "log_rotation_max_bytes": "invalid",
            "primability_cutoff": "invalid",
            "auto_reload_on_startup": True,
        }
    )
    s.from_dict(
        {
            "dark_mode": "system",
            "log_rotation_max_bytes": "200",
            "primability_cutoff": "0.7",
            "auto_reload_on_startup": False,
        }
    )

    # 8. _get_config_path on Windows with and without APPDATA
    with (
        patch("sys.platform", "win32"),
        patch.dict("os.environ", {"APPDATA": r"C:\AppData"}, clear=True),
    ):
        path = s._get_config_path()
        assert str(path).startswith(r"C:\AppData")

    with (
        patch("sys.platform", "win32"),
        patch.dict("os.environ", {}, clear=True),
    ):
        path = s._get_config_path()
        assert "Roaming" in str(path)

    # 9. Local storage exception handling on desktop
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    with patch.object(Path, "mkdir", side_effect=OSError("save error")):
        s.save_to_local(mock_page)

    with (
        patch.object(Path, "exists", return_value=True),
        patch("builtins.open", side_effect=OSError("load error")),
    ):
        s.load_from_local(mock_page)


def test_sequence_separator_synchronisation() -> None:
    """Test sequence separator bidirectional synchronisation."""
    # 1. Initialise with sequence_separator
    settings_a = GUISettings({"sequence_separator": " "})
    assert settings_a["dimer_sequence_separator"] == " "
    assert settings_a["sequence_separator"] == " "

    # 2. Initialise with legacy dimer_sequence_separator
    settings_b = GUISettings({"dimer_sequence_separator": "-"})
    assert settings_b["sequence_separator"] == "-"
    assert settings_b["dimer_sequence_separator"] == "-"

    # 3. from_dict with legacy dimer_sequence_separator
    settings_c = GUISettings()
    settings_c.from_dict({"dimer_sequence_separator": " / "})
    assert settings_c["sequence_separator"] == " / "
    assert settings_c["dimer_sequence_separator"] == " / "


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_web_save_awaits_previous_task_and_no_run_task() -> None:
    """Test web save awaits failing prev task and handles missing run_task."""
    import asyncio
    from types import SimpleNamespace

    settings = GUISettings()
    page_without_run_task = SimpleNamespace(web=True)
    assert settings.save_to_local(page_without_run_task) is None

    async def failing_prev() -> None:
        await asyncio.sleep(0)
        raise RuntimeError("previous web save failed")

    settings._web_save_task = asyncio.create_task(failing_prev())

    mock_page = MagicMock()
    mock_page.web = True
    mock_page.run_task = MagicMock(
        side_effect=lambda fn, *a, **kw: asyncio.create_task(fn(*a, **kw))
    )
    mock_storage: dict[str, Any] = {}

    async def fake_set(key: str, val: Any) -> None:
        mock_storage[key] = val

    mock_prefs = MagicMock()
    mock_prefs.set = MagicMock(side_effect=fake_set)

    settings["dark_mode"] = True
    with patch("flet.SharedPreferences", return_value=mock_prefs):
        task = settings.save_to_local(mock_page, keys=("dark_mode",))
        assert task is not None
        await task

    assert mock_storage == {"amplifyp.settings.dark_mode": True}
