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

"""Update controller for background version checking."""

import asyncio
import logging
import time
from typing import Any

logger = logging.getLogger(__name__)


class UpdateManager:
    """Checks for newer application versions in the background."""

    def __init__(self, controller: Any) -> None:
        """Initialize UpdateManager with a reference to the main controller."""
        self.controller = controller

    def on_update_found(self, latest_version: str) -> None:
        """Update header version text when a new version is found."""
        ctrl = self.controller
        if hasattr(ctrl, "header") and ctrl.header:
            ctrl.header.set_update_available(latest_version)

    async def check_updates_async(self) -> None:
        """Run update checking asynchronously without blocking main thread."""
        from amplifyp import __version__ as current_version
        from amplifyp.gui2.utils.system import (
            fetch_latest_release_version,
            is_newer_version,
            should_check_for_updates,
        )

        ctrl = self.controller
        frequency = ctrl.settings.get(
            "version_checking_frequency", "Once per Month"
        )
        try:
            last_check = float(
                ctrl.settings.get("last_version_check_timestamp", 0.0)
            )
        except (TypeError, ValueError):
            last_check = 0.0
        current_time = float(time.time())

        if not should_check_for_updates(frequency, last_check, current_time):
            return

        loop = asyncio.get_running_loop()
        latest_tag = await loop.run_in_executor(
            None, fetch_latest_release_version
        )

        if latest_tag is not None:
            ctrl.settings["last_version_check_timestamp"] = current_time
            ctrl.settings.save_to_local()

            if is_newer_version(latest_tag, current_version):
                self.on_update_found(latest_tag)
