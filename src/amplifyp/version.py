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

"""Application version, git SHA resolution, and update-check helpers."""

from __future__ import annotations

import importlib
import json
import logging
import os
import re
import subprocess
import urllib.request
from collections.abc import Callable
from importlib.metadata import PackageNotFoundError
from typing import Any, cast

logger = logging.getLogger(__name__)

__all__ = [
    "LATEST_RELEASE_API_URL",
    "_get_sha",
    "_parse_to_tuple",
    "fetch_latest_release_version",
    "get_full_sha",
    "get_git_sha",
    "get_version",
    "is_newer_version",
    "should_check_for_updates",
]

LATEST_RELEASE_API_URL = (
    "https://api.github.com/repos/fangfufu/AmplifyP/releases/latest"
)


def _get_sha(full: bool = False) -> str:
    """Retrieve the git commit SHA (either 40-char full or 7-char short)."""
    attr = "GIT_FULL_SHA" if full else "GIT_SHA"
    for mod_name in ("amplifyp.git_sha", "amplifyp.gui.git_sha"):
        try:
            sha_mod = importlib.import_module(mod_name)
            imported_sha = getattr(sha_mod, attr, None)
            if imported_sha and imported_sha != "unknown":
                return str(imported_sha)
        except ImportError:
            pass

    try:
        js_mod = importlib.import_module("js")
        if hasattr(js_mod, "window") and hasattr(js_mod.window, "__APP_SHA__"):
            sha = str(js_mod.window.__APP_SHA__)
            if sha and sha != "unknown":
                return sha
    except (ImportError, AttributeError):
        pass

    repo_root = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
    git_dir = os.path.join(repo_root, ".git")

    if os.path.exists(git_dir):
        try:
            cmd = (
                ["git", "rev-parse", "HEAD"]
                if full
                else ["git", "rev-parse", "--short", "HEAD"]
            )
            result = subprocess.run(  # noqa: S603
                cmd,
                cwd=repo_root,
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout.strip()
        except (OSError, RuntimeError, subprocess.SubprocessError):
            pass

    try:
        head_path = os.path.join(git_dir, "HEAD")
        if os.path.exists(head_path):
            with open(head_path, encoding="utf-8") as f:
                head_content = f.read().strip()
            if head_content.startswith("ref: "):
                ref_path = head_content.removeprefix("ref: ").strip()
                ref_file = os.path.join(git_dir, ref_path)
                if os.path.exists(ref_file):
                    with open(ref_file, encoding="utf-8") as f:
                        full_sha = f.read().strip()
                    return full_sha if full else full_sha[:7]
                packed_refs = os.path.join(git_dir, "packed-refs")
                if os.path.exists(packed_refs):
                    with open(packed_refs, encoding="utf-8") as f:
                        for line in f:
                            stripped = line.strip()
                            if not stripped or stripped.startswith(("#", "^")):
                                continue
                            parts = stripped.split(maxsplit=1)
                            if len(parts) == 2 and parts[1] == ref_path:
                                packed_sha = parts[0]
                                return packed_sha if full else packed_sha[:7]
            else:
                return head_content if full else head_content[:7]
    except OSError:
        pass

    try:
        dist_sha_path = os.path.join(
            os.path.dirname(__file__), "..", "..", ".git-sha"
        )
        dist_sha_path = os.path.normpath(dist_sha_path)
        if os.path.exists(dist_sha_path):
            with open(dist_sha_path, encoding="utf-8") as f:
                sha_val = f.read().strip()
            if (
                not full
                and len(sha_val) == 40
                and all(c in "0123456789abcdefABCDEF" for c in sha_val)
            ):
                return sha_val[:7]
            return sha_val
    except OSError:
        pass

    return "unknown"


def get_git_sha() -> str:
    """Get the short git commit SHA (7 chars), or 'unknown' if not available."""
    return _get_sha(full=False)


def get_full_sha() -> str:
    """Get the full git commit SHA (40 chars), or 'unknown' if not available."""
    return _get_sha(full=True)


def get_version(git_sha_fn: Callable[[], str] | None = None) -> str:
    """Return version string like 'v0.0.1 (abc1234f)' or 'v0.0.1 (unknown)'."""
    local_logger = logging.getLogger(__name__)
    try:
        from amplifyp import __version__ as pkg_version
    except ImportError:
        try:
            from importlib.metadata import version

            pkg_version = version("amplifyp")
        except PackageNotFoundError:
            local_logger.debug("amplifyp package version not found")
            pkg_version = "unknown"

    sha_fn = get_git_sha if git_sha_fn is None else git_sha_fn
    git_sha = sha_fn()
    return f"{pkg_version} ({git_sha})"


def fetch_latest_release_version() -> str | None:
    """Fetch the latest release version tag from GitHub API."""
    req = urllib.request.Request(
        LATEST_RELEASE_API_URL,
        headers={"User-Agent": "AmplifyP-Updater"},
    )

    try:
        with urllib.request.urlopen(req, timeout=5) as response:  # noqa: S310
            data = json.loads(response.read().decode("utf-8"))
            if isinstance(data, dict):
                payload = cast("dict[str, Any]", data)
                tag_name = payload.get("tag_name")
                if isinstance(tag_name, str):
                    stripped = tag_name.strip()
                    if stripped:
                        return stripped
    except Exception as e:
        logger.warning("Failed to check for updates: %s", e)
    return None


def _parse_to_tuple(v: str) -> tuple[int, ...]:
    """Parse a version string into a comparable tuple of integers."""
    v = re.sub(r"\s+\([^)]+\)$", "", v.strip())
    if v.startswith("v"):
        v = v[1:]

    match = re.match(r"^([\d.]+)(.*)$", v)
    if not match:
        return (0, 0, 0, 0, 0, 0)

    release_str, pre_str = match.groups()
    release_list = [int(x) for x in release_str.split(".") if x.isdigit()]

    while release_list and release_list[-1] == 0:
        release_list.pop()

    while len(release_list) < 4:
        release_list.append(0)

    if not pre_str or pre_str.startswith("+"):
        pre_type = 4
        pre_num = 0
    else:
        pre_str = pre_str.lower()
        pre_word_match = re.search(r"([a-z]+)", pre_str)
        if pre_word_match:
            pre_word = pre_word_match.group(1)
            if pre_word in ("a", "alpha"):
                pre_type = 1
            elif pre_word in ("b", "beta"):
                pre_type = 2
            elif pre_word in ("rc", "preview"):
                pre_type = 3
            elif pre_word in ("post", "rev", "patch"):
                pre_type = 5
            else:
                pre_type = 0
        else:
            pre_type = 0

        num_match = re.search(r"(\d+)", pre_str)
        pre_num = int(num_match.group(1)) if num_match else 0

    return (*release_list, pre_type, pre_num)


def is_newer_version(latest_tag: str, current_ver: str) -> bool:
    """Compare latest tag with current version to see if update is available."""
    latest_tuple = _parse_to_tuple(latest_tag)
    current_tuple = _parse_to_tuple(current_ver)

    return latest_tuple > current_tuple


def should_check_for_updates(
    frequency: str, last_check: float, current_time: float
) -> bool:
    """Determine if update checking is due based on configured frequency."""
    if frequency == "At Startup":
        return True
    if frequency == "Disabled":
        return False

    elapsed = current_time - last_check

    if frequency == "Once per Day":
        return elapsed >= 86400.0
    if frequency == "Once per Week":
        return elapsed >= 604800.0
    if frequency == "Once per Month":
        return elapsed >= 2592000.0

    return False
