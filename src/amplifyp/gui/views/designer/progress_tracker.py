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

"""Shared progress bar tracker for the primer designer views.

Tracks a determinate/indeterminate progress bar plus a tick label and
flushes the page on Flet's event loop. The flush task runs on the loop
thread via ``page.run_task``, so progress ticks written by a background
analysis thread are reliably delivered to the client without waiting
for an unrelated client event (such as a window drag) to trigger a
render.
"""

from __future__ import annotations

import asyncio
import logging
import threading

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings

logger = logging.getLogger(__name__)

_FLUSH_INTERVAL_SECONDS = 0.05


class ProgressTracker:
    """Progress bar + label with event-loop-driven page flushing.

    ``show`` builds the loading body and schedules a flush task on the
    Flet event loop. ``update_progress`` (safe to call from a background
    analysis thread) writes the bar value and label text; the flush
    task picks the values up on its next cycle. ``hide`` stops the
    flush task and drops the control references.
    """

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings,
        hint_text: str,
    ) -> None:
        """Initialise the tracker.

        Args:
            page: Flet page used to schedule the flush task.
            settings: GUI settings for font sizes.
            hint_text: Sub-hint line shown under the bar and label.
        """
        self._page = page
        self._settings = settings
        self._hint_text = hint_text
        self._bar: ft.ProgressBar | None = None
        self._label: ft.Text | None = None
        self._body: ft.Container | None = None
        self._stop = threading.Event()
        self._stop.set()

    @property
    def bar(self) -> ft.ProgressBar | None:
        """The progress bar control while active, else ``None``."""
        return self._bar

    @property
    def label(self) -> ft.Text | None:
        """The progress label control while active, else ``None``."""
        return self._label

    @property
    def loading_body(self) -> ft.Container | None:
        """The loading body container while active, else ``None``."""
        return self._body

    @property
    def is_active(self) -> bool:
        """Whether the tracker has been shown and is not yet stopped."""
        return self._bar is not None and not self._stop.is_set()

    def show(self, total: int = 0) -> ft.Container:
        """Build the loading body and start the event-loop flush task.

        Args:
            total: Total steps. When 0 (unknown), an indeterminate bar
                is shown.

        Returns:
            The loading body container for the host view to insert.
        """
        # Stop any in-flight flush task from a previous show() call so a
        # re-run does not leak a looping task.
        self._stop.set()

        font_small = self._settings.get("font_size_small", 12)
        self._bar = ft.ProgressBar(
            value=0.0 if total > 0 else None,
            expand=True,
            color=GUIColours.PRIMARY,
            bgcolor=GUIColours.SURFACE_VARIANT,
            bar_height=8,
            border_radius=4,
        )
        self._label = ft.Text(
            f"0 / {total}" if total > 0 else "Analysing\u2026",
            italic=True,
            size=font_small,
            color=GUIColours.TEXT_ON_SURFACE,
        )
        self._body = ft.Container(
            content=ft.Column(
                [
                    ft.Row(
                        [
                            self._bar,
                            self._label,
                        ],
                        spacing=10,
                        vertical_alignment=ft.CrossAxisAlignment.CENTER,
                    ),
                    ft.Text(
                        self._hint_text,
                        size=font_small,
                        color=GUIColours.TEXT_ON_SURFACE,
                        opacity=0.6,
                    ),
                ],
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                alignment=ft.MainAxisAlignment.CENTER,
                spacing=8,
            ),
            expand=True,
            alignment=ft.Alignment(0, 0),
            padding=ft.Padding(24, 0, 24, 0),
        )

        stop = threading.Event()
        self._stop = stop
        try:
            self._page.run_task(self._flush_task, stop)
        except (RuntimeError, TypeError) as ex:
            logger.debug("Could not schedule progress flush task: %s", ex)
        return self._body

    async def _flush_task(self, stop: threading.Event) -> None:
        """Flush the page on the event loop until ``stop`` is set.

        Args:
            stop: Event that halts the flush loop.
        """
        while not stop.is_set():
            try:
                self._page.update()
            except RuntimeError:
                logger.debug(
                    "Progress flush stopped: page update failed",
                    exc_info=True,
                )
                break
            await asyncio.sleep(_FLUSH_INTERVAL_SECONDS)

    def update_progress(self, done: int, total: int) -> None:
        """Write current progress values; the flush task renders them.

        Safe to call from a background analysis thread. Stops the flush
        task when the final tick is received.

        Args:
            done: Number of steps completed so far.
            total: Total number of steps.
        """
        if self._bar is None or self._label is None:
            return
        fraction = done / total if total > 0 else 0.0
        self._bar.value = fraction
        self._label.value = f"{done} / {total} ({round(fraction * 100)}%)"
        if done == total:
            self._stop.set()

    def hide(self) -> None:
        """Stop the flush task and drop the control references."""
        self._stop.set()
        self._bar = None
        self._label = None
        self._body = None
