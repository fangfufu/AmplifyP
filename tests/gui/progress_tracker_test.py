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

"""Tests for the shared ProgressTracker (designer 1D/2D progress bar).

Regression focus: the flush must be scheduled on the Flet event loop via
``page.run_task`` (thread-safe) rather than a raw ``threading.Thread``
calling ``page.update()`` from a foreign thread, which Flet 0.8x does not
reliably deliver to the client.
"""

import asyncio
import threading
from unittest.mock import MagicMock, patch

import flet as ft

from amplifyp.gui.settings import GUISettings
from amplifyp.gui.views.designer.progress_tracker import ProgressTracker


async def _instant_sleep(_delay: float) -> None:
    """Awaitable stand-in for ``asyncio.sleep`` (no real delay)."""
    return None


def _make_tracker(page: ft.Page) -> ProgressTracker:
    return ProgressTracker(page=page, settings=GUISettings(), hint_text="hint")


def test_tracker_initial_state_inactive() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)
    assert tracker.bar is None
    assert tracker.label is None
    assert tracker.loading_body is None
    assert tracker.is_active is False


def test_show_determinate_builds_bar_and_schedules_flush() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)

    body = tracker.show(total=10)

    assert isinstance(body, ft.Container)
    assert body is tracker.loading_body
    assert tracker.bar is not None
    assert tracker.bar.value == 0.0
    assert tracker.label is not None
    assert tracker.label.value == "0 / 10"
    assert tracker.is_active is True
    # Flush must be scheduled on the event loop via run_task (not a thread).
    page.run_task.assert_called_once()
    scheduled = page.run_task.call_args[0]
    assert scheduled[0] == tracker._flush_task  # bound-method equality
    assert isinstance(scheduled[1], threading.Event)


def test_show_indeterminate_when_total_zero() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)

    tracker.show(total=0)

    assert tracker.bar is not None
    assert tracker.bar.value is None  # indeterminate
    assert tracker.label is not None
    assert "Analysing" in (tracker.label.value or "")


def test_show_reschedules_on_reentry() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)

    tracker.show(total=5)
    tracker.show(total=8)

    assert page.run_task.call_count == 2
    label = tracker.label
    assert label is not None
    assert label.value == "0 / 8"


def test_show_reentry_stops_prior_flush_task() -> None:
    """A second show() must stop the previous flush task (no leak)."""
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)

    tracker.show(total=5)
    first_stop = page.run_task.call_args[0][1]
    assert first_stop.is_set() is False

    tracker.show(total=8)

    # The first flush task's stop event is now set, so its loop exits.
    assert first_stop.is_set() is True
    # The current stop is a fresh, unset event.
    second_stop = page.run_task.call_args[0][1]
    assert second_stop is not first_stop
    assert second_stop.is_set() is False


def test_update_progress_writes_values_without_stopping() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)
    tracker.show(total=6)

    tracker.update_progress(3, 6)

    bar = tracker.bar
    label = tracker.label
    assert bar is not None
    assert label is not None
    assert abs((bar.value or 0.0) - 0.5) < 1e-9
    assert label.value == "3 / 6 (50%)"
    assert tracker.is_active is True  # not the final tick


def test_update_progress_rounds_percentage() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)
    tracker.show(total=3)

    tracker.update_progress(1, 3)

    bar = tracker.bar
    label = tracker.label
    assert bar is not None
    assert label is not None
    assert abs((bar.value or 0.0) - (1 / 3)) < 1e-9
    assert label.value == "1 / 3 (33%)"


def test_update_progress_guards_zero_total() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)
    tracker.show(total=0)

    tracker.update_progress(1, 0)  # must not raise ZeroDivisionError

    bar = tracker.bar
    label = tracker.label
    assert bar is not None
    assert label is not None
    assert (bar.value or 0.0) == 0.0
    assert label.value == "1 / 0 (0%)"


def test_final_tick_stops_flush() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)
    tracker.show(total=4)
    assert tracker.is_active is True

    tracker.update_progress(4, 4)

    assert tracker.is_active is False


def test_update_progress_noop_when_hidden() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)
    tracker.hide()

    tracker.update_progress(1, 4)  # should not raise

    assert tracker.bar is None
    assert tracker.label is None
    assert tracker.is_active is False


def test_hide_stops_and_clears() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)
    tracker.show(total=4)
    assert tracker.is_active is True

    tracker.hide()

    assert tracker.bar is None
    assert tracker.label is None
    assert tracker.loading_body is None
    assert tracker.is_active is False


def test_flush_task_repeats_until_stopped() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)
    stop = threading.Event()
    calls = []

    def fake_update() -> None:
        calls.append(1)
        if len(calls) >= 3:
            stop.set()

    page.update = fake_update
    with patch.object(asyncio, "sleep", new=_instant_sleep):
        asyncio.run(tracker._flush_task(stop))

    assert len(calls) >= 3
    assert stop.is_set()


def test_flush_task_stops_on_runtime_error() -> None:
    page = MagicMock(spec=ft.Page)
    page.update.side_effect = RuntimeError("boom")
    tracker = _make_tracker(page)
    stop = threading.Event()  # never set

    with patch.object(asyncio, "sleep", new=_instant_sleep):
        asyncio.run(tracker._flush_task(stop))

    assert page.update.call_count == 1
    assert stop.is_set() is False  # loop broke, did not wait on the event


def test_flush_task_exits_immediately_when_pre_stopped() -> None:
    page = MagicMock(spec=ft.Page)
    tracker = _make_tracker(page)
    stop = threading.Event()
    stop.set()

    with patch.object(asyncio, "sleep", new=_instant_sleep):
        asyncio.run(tracker._flush_task(stop))

    assert page.update.call_count == 0
