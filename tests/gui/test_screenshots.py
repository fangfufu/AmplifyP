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

"""Tests for GUI screenshot export functionality."""

from pathlib import Path
from unittest.mock import AsyncMock, MagicMock

import flet as ft
import pytest

import main as main_module
from amplifyp.gui.utils.system import (
    capture_all_views_async,
    capture_view_screenshot_async,
)


def test_cli_screenshots_flags(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test CLI argument parsing with --screenshots flag."""
    mock_run = MagicMock()
    monkeypatch.setattr(ft, "run", mock_run)

    main_module.cli(
        [
            "--state",
            "test_state.yaml",
            "--screenshots",
            "--screenshots-dir",
            "docs/gui/images",
        ]
    )

    assert main_module.state_file == "test_state.yaml"
    assert main_module.export_screenshots is True
    assert main_module.screenshots_dir == "docs/gui/images"


def test_cli_window_dimensions(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test CLI argument parsing with window dimension flags."""
    mock_run = MagicMock()
    monkeypatch.setattr(ft, "run", mock_run)

    main_module.cli(["--window-width", "1920", "--window-height", "1080"])

    assert main_module.window_width == 1920
    assert main_module.window_height == 1080


def test_cli_screenshots_requires_state(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test CLI argument error when --screenshots lacks --state."""
    mock_run = MagicMock()
    monkeypatch.setattr(ft, "run", mock_run)

    with pytest.raises(SystemExit):
        main_module.cli(["--screenshots"])


@pytest.mark.asyncio  # type: ignore[untyped-decorator]
async def test_capture_view_screenshot_async(tmp_path: Path) -> None:
    """Test capture_view_screenshot_async writing PNG bytes."""
    mock_page = MagicMock(spec=ft.Page)
    fake_bytes = b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR"
    mock_page.take_screenshot = AsyncMock(return_value=fake_bytes)

    out_file = await capture_view_screenshot_async(
        mock_page, "input_view", screenshots_dir=tmp_path
    )

    assert out_file == tmp_path / "input_view.png"
    assert out_file.exists()
    assert out_file.read_bytes() == fake_bytes

    mock_page.take_screenshot = AsyncMock(return_value=None)
    with pytest.raises(RuntimeError):
        await capture_view_screenshot_async(
            mock_page, "input_view", screenshots_dir=tmp_path
        )


@pytest.mark.ci  # type: ignore[untyped-decorator]
def test_capture_all_views_async(tmp_path: Path) -> None:
    """Test capture_all_views_async switching views and writing PNG files."""
    import asyncio as _asyncio

    async def _run() -> None:
        mock_controller = MagicMock()
        mock_page = MagicMock(spec=ft.Page)
        fake_bytes = b"PNG_DATA"
        mock_page.take_screenshot = AsyncMock(return_value=fake_bytes)

        mock_controller.page = mock_page
        mock_controller.input_view = MagicMock()
        mock_controller.pcr_view = MagicMock()
        mock_controller.dimers_view = MagicMock()
        mock_controller.auto_close = False
        mock_controller.confirm_exit_async = AsyncMock()

        pcr_btn = MagicMock()
        pcr_btn.disabled = False
        mock_controller.pcr_button_ref.current = pcr_btn

        dimers_btn = MagicMock()
        dimers_btn.disabled = False
        mock_controller.dimers_button_ref.current = dimers_btn

        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(Path, "cwd", lambda: tmp_path)
            await capture_all_views_async(mock_controller)

        assert (tmp_path / "screenshots" / "input_view.png").exists()
        assert (tmp_path / "screenshots" / "pcr_view.png").exists()
        assert (tmp_path / "screenshots" / "primer_dimer_view.png").exists()

        assert mock_controller.switch_view.call_count == 3
        mock_controller.pcr_view.run_pcr.assert_called_once()
        mock_controller.dimers_view.run_analysis.assert_called_once()

    _asyncio.run(_run())
