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


def get_git_sha() -> str:
    """Get the short git commit SHA (7 chars), or 'unknown' if not available."""
    try:
        from amplifyp.gui.git_sha import GIT_SHA

        if GIT_SHA and GIT_SHA != "unknown":
            return str(GIT_SHA)
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
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except OSError:
        pass

    try:
        git_dir = os.path.join(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                )
            ),
            ".git",
        )
        head_path = os.path.join(git_dir, "HEAD")
        if os.path.exists(head_path):
            with open(head_path) as f:
                head_content = f.read().strip()
            if head_content.startswith("ref: refs/heads/"):
                ref_path = head_content.replace("ref: refs/heads/", "")
                ref_file = os.path.join(git_dir, ref_path)
                if os.path.exists(ref_file):
                    with open(ref_file) as f:
                        full_sha = f.read().strip()
                    return full_sha[:7]
            else:
                return head_content[:7]
    except OSError:
        pass

    try:
        dist_sha_path = os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "..", ".git-sha"
        )
        dist_sha_path = os.path.normpath(dist_sha_path)
        if os.path.exists(dist_sha_path):
            with open(dist_sha_path) as f:
                return f.read().strip()
    except OSError:
        pass

    return "unknown"


def get_full_sha() -> str:
    """Get the full git commit SHA (40 chars), or 'unknown' if not available."""
    try:
        from amplifyp.gui.git_sha import GIT_FULL_SHA

        if GIT_FULL_SHA and GIT_FULL_SHA != "unknown":
            return str(GIT_FULL_SHA)
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
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except OSError:
        pass

    try:
        git_dir = os.path.join(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                )
            ),
            ".git",
        )
        for ref in ("refs/heads/main", "refs/heads/master"):
            ref_file = os.path.join(git_dir, ref)
            if os.path.exists(ref_file):
                with open(ref_file) as f:
                    return f.read().strip()
    except OSError:
        pass

    try:
        dist_sha_path = os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "..", ".git-sha"
        )
        dist_sha_path = os.path.normpath(dist_sha_path)
        if os.path.exists(dist_sha_path):
            with open(dist_sha_path) as f:
                return f.read().strip()
    except OSError:
        pass

    return "unknown"


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
