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

"""Tests for amplifyp.version."""

from __future__ import annotations

import json
import subprocess
import sys
from importlib.metadata import PackageNotFoundError
from typing import Any
from unittest.mock import MagicMock, mock_open, patch

from amplifyp.version import (
    _get_sha,
    _parse_to_tuple,
    fetch_latest_release_version,
    get_full_sha,
    get_git_sha,
    get_version,
    is_newer_version,
    should_check_for_updates,
)


def test_get_sha_from_core_and_gui_git_sha_modules() -> None:
    mock_core_sha = MagicMock()
    mock_core_sha.GIT_SHA = "core123"
    mock_core_sha.GIT_FULL_SHA = "corefull123"
    with patch.dict(sys.modules, {"amplifyp.git_sha": mock_core_sha}):
        assert get_git_sha() == "core123"
        assert get_full_sha() == "corefull123"

    mock_gui_sha = MagicMock()
    mock_gui_sha.GIT_SHA = "gui1234"
    mock_gui_sha.GIT_FULL_SHA = "guifull1234"
    with patch.dict(
        sys.modules,
        {"amplifyp.git_sha": None, "amplifyp.gui.git_sha": mock_gui_sha},
    ):
        assert get_git_sha() == "gui1234"
        assert get_full_sha() == "guifull1234"


def test_get_sha_js_git_and_fallbacks() -> None:
    mock_js = MagicMock()
    mock_js.window.__APP_SHA__ = "jssha99"
    with patch.dict(
        sys.modules,
        {"amplifyp.git_sha": None, "amplifyp.gui.git_sha": None, "js": mock_js},
    ):
        assert get_git_sha() == "jssha99"

    mock_result = MagicMock(returncode=0, stdout="deadbee\n")
    with (
        patch.dict(
            sys.modules,
            {
                "amplifyp.git_sha": None,
                "amplifyp.gui.git_sha": None,
                "js": None,
            },
        ),
        patch("subprocess.run", return_value=mock_result),
    ):
        assert get_git_sha() == "deadbee"

    with (
        patch.dict(
            sys.modules,
            {
                "amplifyp.git_sha": None,
                "amplifyp.gui.git_sha": None,
                "js": None,
            },
        ),
        patch("subprocess.run", side_effect=subprocess.SubprocessError),
        patch("os.path.exists", return_value=False),
    ):
        assert _get_sha(full=False) == "unknown"

    # RuntimeError from subprocess.run (e.g. Pyodide/Emscripten) falls back
    # cleanly to reading .git/HEAD rather than propagating.
    files = {
        "HEAD": "ref: refs/heads/dev\n",
        "dev": "fedcba9876543210fedcba9876543210fedcba98\n",
    }

    def open_runtime_fallback(file: Any, *args: Any, **kwargs: Any) -> Any:
        filename = str(file).replace("\\", "/")
        for suffix, data in files.items():
            if filename.endswith(suffix):
                return mock_open(read_data=data).return_value
        raise FileNotFoundError(filename)

    def git_exists_runtime_fallback(path: str) -> bool:
        norm = str(path).replace("\\", "/")
        return norm.endswith((".git", ".git/HEAD", ".git/refs/heads/dev"))

    with (
        patch.dict(
            sys.modules,
            {
                "amplifyp.git_sha": None,
                "amplifyp.gui.git_sha": None,
                "js": None,
            },
        ),
        patch(
            "subprocess.run",
            side_effect=RuntimeError("can't start new thread"),
        ),
        patch("os.path.exists", side_effect=git_exists_runtime_fallback),
        patch("builtins.open", side_effect=open_runtime_fallback),
    ):
        assert _get_sha(full=False) == "fedcba9"
        assert _get_sha(full=True) == "fedcba9876543210fedcba9876543210fedcba98"


def test_get_sha_dot_git_and_dist_sha() -> None:
    files = {
        "HEAD": "ref: refs/heads/main\n",
        "main": "0123456789abcdef0123456789abcdef01234567\n",
    }

    def open_side_effect(file: Any, *args: Any, **kwargs: Any) -> Any:
        filename = str(file).replace("\\", "/")
        for suffix, data in files.items():
            if filename.endswith(suffix):
                return mock_open(read_data=data).return_value
        raise FileNotFoundError(filename)

    def dot_git_exists(path: str) -> bool:
        norm = str(path).replace("\\", "/")
        return norm.endswith((".git/HEAD", ".git/refs/heads/main"))

    with (
        patch.dict(
            sys.modules,
            {
                "amplifyp.git_sha": None,
                "amplifyp.gui.git_sha": None,
                "js": None,
            },
        ),
        patch("subprocess.run", side_effect=OSError),
        patch("os.path.exists", side_effect=dot_git_exists),
        patch("builtins.open", side_effect=open_side_effect),
    ):
        assert _get_sha(full=False) == "0123456"
        assert _get_sha(full=True) == "0123456789abcdef0123456789abcdef01234567"

    with (
        patch.dict(
            sys.modules,
            {
                "amplifyp.git_sha": None,
                "amplifyp.gui.git_sha": None,
                "js": None,
            },
        ),
        patch("subprocess.run", side_effect=OSError),
        patch("os.path.exists", side_effect=lambda p: ".git-sha" in str(p)),
        patch("builtins.open", mock_open(read_data="dist123\n")),
    ):
        assert _get_sha(full=False) == "dist123"

    full_hex = "a" * 40
    with (
        patch.dict(
            sys.modules,
            {
                "amplifyp.git_sha": None,
                "amplifyp.gui.git_sha": None,
                "js": None,
            },
        ),
        patch("subprocess.run", side_effect=OSError),
        patch("os.path.exists", side_effect=lambda p: ".git-sha" in str(p)),
        patch("builtins.open", mock_open(read_data=f"{full_hex}\n")),
    ):
        assert _get_sha(full=False) == "a" * 7
        assert _get_sha(full=True) == full_hex

    # packed-refs fallback when loose ref file does not exist
    packed_files = {
        ".git/HEAD": "ref: refs/heads/main\n",
        ".git/packed-refs": (
            "# pack-refs\n"
            "fedcba9876543210fedcba9876543210fedcba98 refs/heads/main\n"
        ),
    }

    def packed_open(file: Any, *args: Any, **kwargs: Any) -> Any:
        filename = str(file).replace("\\", "/")
        for suffix, data in packed_files.items():
            if filename.endswith(suffix):
                return mock_open(read_data=data).return_value
        raise FileNotFoundError(filename)

    def packed_exists(path: str) -> bool:
        norm = str(path).replace("\\", "/")
        return norm.endswith((".git/HEAD", ".git/packed-refs"))

    with (
        patch.dict(
            sys.modules,
            {
                "amplifyp.git_sha": None,
                "amplifyp.gui.git_sha": None,
                "js": None,
            },
        ),
        patch("subprocess.run", side_effect=OSError),
        patch("os.path.exists", side_effect=packed_exists),
        patch("builtins.open", side_effect=packed_open),
    ):
        assert _get_sha(full=False) == "fedcba9"
        assert _get_sha(full=True) == "fedcba9876543210fedcba9876543210fedcba98"


def test_get_version_uses_package_and_metadata_fallback() -> None:
    assert get_version(git_sha_fn=lambda: "abc1234").endswith("(abc1234)")

    import amplifyp

    orig_version = getattr(amplifyp, "__version__", None)
    try:
        if hasattr(amplifyp, "__version__"):
            del amplifyp.__version__
        with patch("importlib.metadata.version", return_value="1.2.3"):
            assert (
                get_version(git_sha_fn=lambda: "def5678") == "1.2.3 (def5678)"
            )
        with patch(
            "importlib.metadata.version", side_effect=PackageNotFoundError
        ):
            assert (
                get_version(git_sha_fn=lambda: "def5678") == "unknown (def5678)"
            )
    finally:
        if orig_version is not None:
            amplifyp.__version__ = orig_version


def test_version_comparison_and_update_schedule() -> None:
    assert _parse_to_tuple("v1.2.3") == (1, 2, 3, 0, 4, 0)
    assert _parse_to_tuple("1.27.0 (abc1234)") == (1, 27, 0, 0, 4, 0)
    assert is_newer_version("v1.2.0", "1.1.9 (abc)")
    assert not is_newer_version("v1.27.0", "1.27.0 (abc1234)")
    assert is_newer_version("v1.28.0", "1.27.0 (abc1234)")
    assert not is_newer_version("v1.0.0rc1", "1.0.0")
    assert should_check_for_updates("At Startup", 0.0, 0.0)
    assert not should_check_for_updates("Disabled", 0.0, 1e9)
    assert should_check_for_updates("Once per Day", 0.0, 86400.0)
    assert not should_check_for_updates("Once per Week", 0.0, 100.0)
    assert should_check_for_updates("Once per Month", 0.0, 2592000.0)
    assert not should_check_for_updates("Unknown", 0.0, 1e9)


def test_fetch_latest_release_version() -> None:
    payload = json.dumps({"tag_name": " v2.0.0 "}).encode("utf-8")
    mock_resp = MagicMock()
    mock_resp.read.return_value = payload
    mock_cm = MagicMock()
    mock_cm.__enter__.return_value = mock_resp

    with patch("urllib.request.urlopen", return_value=mock_cm):
        assert fetch_latest_release_version() == "v2.0.0"

    mock_resp.read.return_value = json.dumps({"tag_name": "   "}).encode(
        "utf-8"
    )
    with patch("urllib.request.urlopen", return_value=mock_cm):
        assert fetch_latest_release_version() is None

    with patch("urllib.request.urlopen", side_effect=OSError("offline")):
        assert fetch_latest_release_version() is None


def test_get_sha_uses_checkout_root_independent_of_cwd(tmp_path: Any) -> None:
    """_get_sha runs git in AmplifyP checkout and ignores unrelated repos."""
    import os

    import amplifyp.version as version_mod

    expected_root = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(version_mod.__file__)))
    )
    mock_result = MagicMock(returncode=0, stdout="1234abc\n")
    with (
        patch.dict(
            sys.modules,
            {
                "amplifyp.git_sha": None,
                "amplifyp.gui.git_sha": None,
                "js": None,
            },
        ),
        patch(
            "os.path.exists",
            side_effect=lambda p: str(p) == os.path.join(expected_root, ".git"),
        ),
        patch("subprocess.run", return_value=mock_result) as mock_run,
    ):
        assert _get_sha(full=False) == "1234abc"
        assert mock_run.call_args.kwargs["cwd"] == expected_root

    # When AmplifyP has no .git directory (e.g. packaged install inside another
    # repo's working directory), git must not be invoked against caller cwd.
    other_repo = tmp_path / "other_repo"
    (other_repo / ".git").mkdir(parents=True)
    orig_cwd = os.getcwd()
    try:
        os.chdir(other_repo)
        with (
            patch.dict(
                sys.modules,
                {
                    "amplifyp.git_sha": None,
                    "amplifyp.gui.git_sha": None,
                    "js": None,
                },
            ),
            patch("os.path.exists", return_value=False),
            patch("subprocess.run") as mock_run,
        ):
            assert _get_sha(full=False) == "unknown"
            mock_run.assert_not_called()
    finally:
        os.chdir(orig_cwd)
