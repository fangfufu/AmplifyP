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

"""Git commitment details and versioning utilities."""

import os
import subprocess
from importlib.metadata import PackageNotFoundError


def _get_sha(full: bool = False) -> str:
    """Retrieve the git commit SHA (either 40-char full or 7-char short)."""
    try:
        if full:
            from amplifyp.gui.git_sha import GIT_FULL_SHA as imported_sha
        else:
            from amplifyp.gui.git_sha import GIT_SHA as imported_sha

        if imported_sha and imported_sha != "unknown":
            return str(imported_sha)
    except ImportError:
        pass

    try:
        import js  # type: ignore[import-not-found, unused-ignore]

        if hasattr(js, "window") and hasattr(js.window, "__APP_SHA__"):
            sha = str(js.window.__APP_SHA__)
            if sha and sha != "unknown":
                return sha
    except (ImportError, AttributeError):
        pass

    try:
        cmd = (
            ["git", "rev-parse", "HEAD"]
            if full
            else ["git", "rev-parse", "--short", "HEAD"]
        )
        result = subprocess.run(  # noqa: S603
            cmd,
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass

    try:
        git_dir = os.path.join(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(
                        os.path.dirname(
                            os.path.dirname(os.path.abspath(__file__))
                        )
                    )
                )
            ),
            ".git",
        )
        head_path = os.path.join(git_dir, "HEAD")
        if os.path.exists(head_path):
            with open(head_path, encoding="utf-8") as f:
                head_content = f.read().strip()
            if head_content.startswith("ref: refs/heads/"):
                ref_path = head_content.replace("ref: refs/heads/", "")
                ref_file = os.path.join(git_dir, ref_path)
                if os.path.exists(ref_file):
                    with open(ref_file, encoding="utf-8") as f:
                        full_sha = f.read().strip()
                    return full_sha if full else full_sha[:7]
            else:
                return head_content if full else head_content[:7]
    except OSError:
        pass

    try:
        dist_sha_path = os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "..", ".git-sha"
        )
        dist_sha_path = os.path.normpath(dist_sha_path)
        if os.path.exists(dist_sha_path):
            with open(dist_sha_path, encoding="utf-8") as f:
                return f.read().strip()
    except OSError:
        pass

    return "unknown"


def get_git_sha() -> str:
    """Get the short git commit SHA (7 chars), or 'unknown' if not available."""
    return _get_sha(full=False)


def get_full_sha() -> str:
    """Get the full git commit SHA (40 chars), or 'unknown' if not available."""
    return _get_sha(full=True)


def get_version() -> str:
    """Return version string like 'v0.0.1 (abc1234f)' or 'v0.0.1 (unknown)'."""
    import logging

    logger = logging.getLogger(__name__)
    try:
        from amplifyp import __version__ as pkg_version
    except ImportError:
        try:
            from importlib.metadata import version

            pkg_version = version("amplifyp")
        except PackageNotFoundError:
            logger.debug("amplifyp package version not found")
            pkg_version = "unknown"

    git_sha = get_git_sha()
    return f"{pkg_version} ({git_sha})"
