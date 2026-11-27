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

"""Shared cooperative asyncio helpers for primer designers."""

from __future__ import annotations

import asyncio
import math
import time
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import threading
    from collections.abc import Iterator

ASYNC_YIELD_SLEEP_SECONDS = 0.001
DEFAULT_ASYNC_YIELD_INTERVAL = 0.02


def validate_yield_interval(yield_interval: float) -> None:
    """Validate that ``yield_interval`` is finite and non-negative."""
    if (
        isinstance(yield_interval, bool)
        or not math.isfinite(yield_interval)
        or yield_interval < 0
    ):
        raise ValueError("yield_interval must be finite and non-negative.")


async def run_cooperative_analysis(
    steps: Iterator[tuple[int, int]],
    *,
    cancel_event: threading.Event | None,
    yield_interval: float,
) -> None:
    """Consume a ``(done, total)`` step iterator, yielding periodically.

    Yields once before starting work (unless already cancelled) and again
    between steps whenever ``yield_interval`` seconds have elapsed so
    single-threaded runtimes (such as Pyodide / WebAssembly) can flush UI
    updates and handle abort events without spawning OS threads.
    """
    if cancel_event is None or not cancel_event.is_set():
        await asyncio.sleep(ASYNC_YIELD_SLEEP_SECONDS)
    last_yield = time.monotonic()
    for done, total in steps:
        if done < total and (cancel_event is None or not cancel_event.is_set()):
            now = time.monotonic()
            if now - last_yield >= yield_interval:
                await asyncio.sleep(ASYNC_YIELD_SLEEP_SECONDS)
                last_yield = time.monotonic()
