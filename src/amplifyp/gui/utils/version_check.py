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
            data = json.loads(response.read().decode("utf-8"))
            if isinstance(data, dict):
                tag_name = data.get("tag_name")
                if isinstance(tag_name, str):
                    return tag_name.strip()
    except Exception as e:
        logger.warning("Failed to check for updates: %s", e)
    return None


def _parse_to_tuple(v: str) -> tuple[int, ...]:
    """Parse a version string into a comparable tuple of integers."""
    import re

    if v.startswith("v"):
        v = v[1:]

    match = re.match(r"^([\d.]+)(.*)$", v)
    if not match:
        return (0, 0, 0, 0, 0, 0)

    release_str, pre_str = match.groups()
    release_list = [int(x) for x in release_str.split(".") if x.isdigit()]

    while release_list and release_list[-1] == 0:
        release_list.pop()

    # Pad release list to a fixed length of 4 to ensure consistent
    # tuple comparison
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
