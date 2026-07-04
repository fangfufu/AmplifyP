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

"""Version checking utilities for AmplifyP."""

import json
import logging
import urllib.request

logger = logging.getLogger(__name__)


def fetch_latest_release_version() -> str | None:
    """Fetch the latest release version tag from GitHub API."""
    url = "https://api.github.com/repos/fangfufu/AmplifyP/releases/latest"
    req = urllib.request.Request(  # noqa: S310
        url,
        headers={"User-Agent": "AmplifyP-Updater"},
    )
    try:
        with urllib.request.urlopen(req, timeout=5) as response:  # noqa: S310
            data = json.loads(response.read().decode())
            tag_name = data.get("tag_name")
            if isinstance(tag_name, str):
                return tag_name.strip()
    except Exception as e:
        logger.warning("Failed to check for updates: %s", e)
    return None


def is_newer_version(latest_tag: str, current_ver: str) -> bool:
    """Compare latest tag with current version to see if update is available."""

    def parse_to_tuple(v: str) -> tuple[int, ...]:
        if v.startswith("v"):
            v = v[1:]
        base = v.split("-")[0].split("+")[0]
        try:
            return tuple(int(x) for x in base.split(".") if x.isdigit())
        except ValueError:
            return (0,)

    latest_tuple = parse_to_tuple(latest_tag)
    current_tuple = parse_to_tuple(current_ver)

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
