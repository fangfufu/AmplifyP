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

import asyncio
import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast
from unittest.mock import MagicMock, mock_open, patch

import flet as ft
import pytest
import yaml

from amplifyp.gui.settings import GUISettings


def _make_desktop_page() -> MagicMock:
    """Create a mock desktop Flet page."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.web = False
    return mock_page


def _load_desktop_yaml(settings: GUISettings, yaml_content: str) -> None:
    """Load YAML content into settings using a mocked local config file."""
    mock_page = _make_desktop_page()
    mock_path = MagicMock(spec=Path)
    mock_path.exists.return_value = True
    m_open = mock_open(read_data=yaml_content)

    with (
        patch.object(settings, "_get_config_path", return_value=mock_path),
        patch("builtins.open", m_open),
    ):
        settings.load_from_local(mock_page)


def _make_web_page(
    *, sync: bool = False
) -> tuple[MagicMock, list[asyncio.Task[Any]]]:
    """Create a mock web Flet page and a list collecting scheduled tasks."""
    mock_page = MagicMock()
    mock_page.web = True
    tasks: list[asyncio.Task[Any]] = []

    if sync:
        mock_page.run_task = MagicMock(
            side_effect=lambda fn, *args, **kwargs: asyncio.run(
                fn(*args, **kwargs)
            )
        )
    else:

        def _run_task(fn: Any, *args: Any, **kwargs: Any) -> asyncio.Task[Any]:
            task = asyncio.create_task(fn(*args, **kwargs))
            tasks.append(task)
            return task

        mock_page.run_task = MagicMock(side_effect=_run_task)

    return mock_page, tasks


def _make_mock_prefs(
    storage: dict[str, Any] | None = None,
    *,
    set_delay: float = 0.0,
) -> tuple[MagicMock, dict[str, Any]]:
    """Create a mock SharedPreferences backed by an in-memory dict."""
    mock_storage: dict[str, Any] = {} if storage is None else storage

    async def mock_set(key: str, value: Any) -> bool:
        if set_delay > 0:
            await asyncio.sleep(set_delay)
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
    return mock_prefs, mock_storage


async def _run_web_startup(
    mock_prefs: MagicMock,
) -> tuple[list[str], list[str], GUISettings]:
    """Run app_main on a mock web page and return events and settings."""
    from amplifyp.gui.app import main as app_main

    mock_page, tasks = _make_web_page()
    events: list[str] = []

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
        events_before_gather = list(events)
        await asyncio.gather(*tasks)
        passed_settings: GUISettings = cls.call_args.kwargs["settings"]

    return events_before_gather, events, passed_settings


@pytest.mark.parametrize(  # type: ignore[untyped-decorator]
    ("platform", "env", "expected_substrings", "expected_posix"),
    [
        (
            "win32",
            {"APPDATA": r"C:\Users\test\AppData\Roaming"},
            ("AppData", "AmplifyP"),
            None,
        ),
        ("win32", {}, ("Roaming", "AmplifyP"), None),
        (
            "darwin",
            {"HOME": "/Users/test"},
            (),
            "/Users/test/Library/Application Support/AmplifyP/settings.yaml",
        ),
        (
            "linux",
            {"XDG_CONFIG_HOME": "/custom/config"},
            (),
            "/custom/config/amplifyp/settings.yaml",
        ),
        (
            "linux",
            {"HOME": "/home/test"},
            (),
            "/home/test/.config/amplifyp/settings.yaml",
        ),
    ],
)
def test_get_config_path(
    platform: str,
    env: dict[str, str],
    expected_substrings: tuple[str, ...],
    expected_posix: str | None,
) -> None:
    """Test that OS-specific config paths are resolved correctly."""
    settings = GUISettings()
    with (
        patch("sys.platform", platform),
        patch.dict("os.environ", env, clear=True),
    ):
        path = settings._get_config_path()
        assert path.name == "settings.yaml"
        for substring in expected_substrings:
            assert substring in str(path)
        if expected_posix is not None:
            assert path.as_posix() == expected_posix


def test_desktop_save_to_local() -> None:
    """Test saving settings to local YAML file on desktop."""
    settings = GUISettings()
    settings["primability_cutoff"] = "0.75"

    mock_page = _make_desktop_page()
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
    _load_desktop_yaml(settings, "primability_cutoff: '0.65'\ndark_mode: true")

    assert settings["primability_cutoff"] == "0.65"
    assert settings["dark_mode"] is True


def test_desktop_load_from_local_corrupted() -> None:
    """Test loading from a corrupted YAML file falls back to defaults."""
    settings = GUISettings()
    settings["primability_cutoff"] = "0.99"

    # Invalid YAML syntax; should not raise exception
    _load_desktop_yaml(settings, "{invalid: yaml:")

    # Since it failed, the settings should still be whatever they were
    assert settings["primability_cutoff"] == "0.99"


def test_web_save_and_load_client_storage() -> None:
    """Test saving and loading settings on web using page.client_storage."""
    settings = GUISettings()
    settings["primability_cutoff"] = "0.85"
    settings["dark_mode"] = True

    mock_page, _ = _make_web_page(sync=True)
    mock_prefs, mock_storage = _make_mock_prefs()

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
    settings = GUISettings()
    settings["primer_info_panel_height"] = 420

    mock_page, _ = _make_web_page()
    mock_prefs, mock_storage = _make_mock_prefs(set_delay=0.01)

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

    async def delayed_contains(key: str) -> bool:
        await asyncio.sleep(0.01)
        return key == "amplifyp.settings.dark_mode"

    async def delayed_get(key: str) -> Any:
        return True

    mock_prefs = MagicMock()
    mock_prefs.contains_key = MagicMock(side_effect=delayed_contains)
    mock_prefs.get = MagicMock(side_effect=delayed_get)

    events_before_gather, events, passed_settings = await _run_web_startup(
        mock_prefs
    )
    assert events_before_gather == []
    assert events == ["logging", "ctrl_init:{'load_settings': False}"]
    assert passed_settings["dark_mode"] is True


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_web_startup_recovers_from_preference_error(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test web startup logs preference errors and continues with defaults."""

    async def failing_contains(key: str) -> bool:
        raise RuntimeError("client storage unavailable")

    mock_prefs = MagicMock()
    mock_prefs.contains_key = MagicMock(side_effect=failing_contains)

    with caplog.at_level(logging.WARNING, logger="amplifyp.gui.app"):
        _, events, passed_settings = await _run_web_startup(mock_prefs)

    assert events == ["logging", "ctrl_init:{'load_settings': False}"]
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
    from amplifyp.dna import Primer
    from amplifyp.gui.colours import GUIColours

    # 1. Initialisation with settings_dict
    s1 = GUISettings(
        settings_dict={"colour_deficient": "yes", "dark_mode": "true"}
    )
    assert s1["colour_deficient"] is True
    assert s1["dark_mode"] is True
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
    for method, key, default in (
        (s._safe_float, "primability_cutoff", 9.9),
        (s._safe_int, "log_rotation_max_bytes", 5),
    ):
        assert method("missing_key", default) == default
        for bad_val in (None, "invalid", ""):
            s._settings[key] = bad_val
            assert method(key, default) == default

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

    # 8. Local storage exception handling on desktop
    mock_page = _make_desktop_page()
    with patch.object(Path, "mkdir", side_effect=OSError("save error")):
        s.save_to_local(mock_page)

    with (
        patch.object(Path, "exists", return_value=True),
        patch("builtins.open", side_effect=OSError("load error")),
    ):
        s.load_from_local(mock_page)


def test_sequence_separator_synchronisation() -> None:
    """Test sequence separator bidirectional synchronisation."""
    for key, sep in (
        ("sequence_separator", " "),
        ("dimer_sequence_separator", "-"),
    ):
        settings = GUISettings({key: sep})
        assert settings["sequence_separator"] == sep
        assert settings["dimer_sequence_separator"] == sep

    settings_from_dict = GUISettings()
    settings_from_dict.from_dict({"dimer_sequence_separator": " / "})
    assert settings_from_dict["sequence_separator"] == " / "
    assert settings_from_dict["dimer_sequence_separator"] == " / "


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_web_save_awaits_previous_task_and_no_run_task() -> None:
    """Test web save awaits failing prev task and handles missing run_task."""
    settings = GUISettings()
    page_without_run_task = SimpleNamespace(web=True)
    assert settings.save_to_local(page_without_run_task) is None

    async def failing_prev() -> None:
        await asyncio.sleep(0)
        raise RuntimeError("previous web save failed")

    settings._web_save_task = asyncio.create_task(failing_prev())

    mock_page, _ = _make_web_page()
    mock_prefs, mock_storage = _make_mock_prefs()

    settings["dark_mode"] = True
    with patch("flet.SharedPreferences", return_value=mock_prefs):
        task = settings.save_to_local(mock_page, keys=("dark_mode",))
        assert task is not None
        await task

    assert mock_storage == {"amplifyp.settings.dark_mode": True}


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_web_save_prunes_stale_weights_and_fallback_load() -> None:
    """Test web save prunes stale weights and fallback load reads extras."""
    from amplifyp.gui.colours import GUIColours

    settings = GUISettings()
    mock_page, _ = _make_web_page()

    # 1. Fallback load_from_local_async (no get_keys) reads extra weight idx 15
    storage: dict[str, Any] = {
        f"amplifyp.settings.match_weight_{i}": str(i) for i in range(16)
    }
    mock_prefs, _ = _make_mock_prefs(storage=storage)
    del mock_prefs.get_keys

    with patch("flet.SharedPreferences", return_value=mock_prefs):
        await settings.load_from_local_async(mock_page)
    assert settings["match_weight_15"] == "15"

    # 2. Full web save with get_keys + remove prunes stale lengthwise weights
    storage_prune: dict[str, Any] = {
        "amplifyp.settings.match_weight_20": "99",
        "amplifyp.settings.run_weight_10": "88",
        "amplifyp.settings.other_key": "keep",
    }

    async def _set_key(k: str, v: Any) -> bool:
        storage_prune[k] = v
        return True

    async def _get_keys(prefix: str) -> list[str]:
        return [k for k in storage_prune if k.startswith(prefix)]

    async def _remove(k: str) -> bool:
        storage_prune.pop(k, None)
        return True

    prune_prefs = MagicMock()
    prune_prefs.set = _set_key
    prune_prefs.get_keys = _get_keys
    prune_prefs.remove = _remove

    with patch("flet.SharedPreferences", return_value=prune_prefs):
        task = settings.save_to_local(mock_page)
        assert task is not None
        await task

    assert "amplifyp.settings.match_weight_20" not in storage_prune
    assert "amplifyp.settings.run_weight_10" not in storage_prune
    assert storage_prune["amplifyp.settings.other_key"] == "keep"

    # 3. reset_to_defaults handles non-"system" default dark_mode cleanly
    prev_dark = GUIColours.dark_mode
    try:
        fake_defaults = SimpleNamespace(
            _settings={"dark_mode": True, "colour_deficient": False}
        )
        with patch(
            "amplifyp.gui.settings.GUISettings", return_value=fake_defaults
        ):
            settings.reset_to_defaults()
        assert GUIColours.dark_mode is True
    finally:
        GUIColours.dark_mode = prev_dark


def test_from_dict_validates_and_canonicalizes_weight_keys() -> None:
    """Validate weight keys and normalize numeric suffixes."""
    settings = GUISettings()
    original_default = settings["match_weight_default"]
    original_zero = settings["match_weight_0"]

    # 1. Invalid-only weight keys must leave existing weights unchanged
    settings.from_dict({"match_weight_typo": "99", "run_weight_bad": "42"})
    assert settings["match_weight_default"] == original_default
    assert settings["match_weight_0"] == original_zero
    assert "match_weight_typo" not in settings._settings
    assert "run_weight_bad" not in settings._settings

    # 2. Leading-zero numeric suffix is canonicalized to match_weight_0
    settings.from_dict({"match_weight_00": "77"})
    assert settings["match_weight_0"] == "77"
    assert "match_weight_00" not in settings._settings

    # 3. Exact canonical key takes precedence over non-canonical duplicate
    settings.from_dict({"match_weight_0": "55", "match_weight_00": "88"})
    assert settings["match_weight_0"] == "55"
    assert "match_weight_00" not in settings._settings

    # 4. Default-only replacement still clears non-default overrides
    settings.from_dict({"match_weight_default": "12"})
    assert settings["match_weight_default"] == "12"
    assert "match_weight_0" not in settings._settings
