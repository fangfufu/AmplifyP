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

"""Base view component for primer designer views."""

from __future__ import annotations

import asyncio
import logging
import threading
from collections.abc import Callable, Coroutine
from typing import Any

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils import data_helpers
from amplifyp.gui.utils import system as system_utils
from amplifyp.gui.utils.gui_helpers import NotificationHelper, show_error_dialog
from amplifyp.io import safe_dump_yaml, safe_load_yaml

logger = logging.getLogger(__name__)


class BaseDesignerView(ft.Row):  # type: ignore[misc]
    """Base view providing split layout, resizers, and card manager."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
        on_add_primer: Callable[[str, str], tuple[bool, str]] | None = None,
    ) -> None:
        """Initialise BaseDesignerView.

        Args:
            page: Flet page instance.
            input_data: Application GUIInput instance.
            settings: Application GUISettings instance.
            on_add_primer: Optional callback to add a primer to Input view.
        """
        super().__init__(expand=True, spacing=0)
        self.app_page = page
        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()
        self.on_add_primer = on_add_primer
        self._active_cards: list[ft.Card] = []
        self._background_tasks: set[asyncio.Task[Any]] = set()
        self._analysis_running: bool = False
        self._cancel_event: threading.Event | None = None

        # Left panel horizontal and vertical divider resizers
        self.main_h_divider = ft.GestureDetector(
            on_pan_update=self._on_h_pan_update,
            content=ft.Container(
                width=6,
                bgcolor=GUIColours.DIVIDER_GREY,
                border_radius=3,
                margin=ft.Margin.symmetric(horizontal=4),
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
        )

        self.left_v_divider = ft.GestureDetector(
            on_pan_update=self._on_v_pan_update,
            content=ft.Container(
                height=6,
                bgcolor=GUIColours.DIVIDER_GREY,
                border_radius=3,
                margin=ft.Margin.symmetric(vertical=4),
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_UP_DOWN,
        )

        # Right-hand panel components
        self.right_cards_list = ft.ListView(
            expand=True, spacing=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.clear_cards_button = ft.TextButton(
            "Clear Cards",
            icon=ft.Icons.DELETE_SWEEP,
            tooltip="Clear All Cards",
            on_click=self._clear_all_cards,
            visible=False,
        )

        self.right_title = ft.Text(
            "Detail Cards",
            weight=ft.FontWeight.BOLD,
            size=self.settings.get("font_size_header", 18),
        )

        self.right_header = ft.Row(
            [
                self.right_title,
                self.clear_cards_button,
            ],
            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
        )

        self.right_container = ft.Container(
            content=ft.Column(
                [
                    self.right_header,
                    ft.Container(content=self.right_cards_list, expand=True),
                ],
                spacing=8,
            ),
            expand=True,
            padding=10,
        )

        # Placeholder for left container to be initialized by subclasses
        self.left_container = ft.Container(expand=True, padding=5)
        self.top_left_container = ft.Container()
        self.bottom_left_container = ft.Container()

    def _on_h_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle horizontal resizing between left and right panels."""
        delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
        if self.left_container.width is None:
            self.left_container.expand = None
            self.right_container.expand = True
            page_w = (
                self.app_page.width
                if hasattr(self.app_page, "width")
                and isinstance(self.app_page.width, (int, float))
                else 800.0
            )
            self.left_container.width = max(
                250.0, float(page_w) * 0.5 + delta_x
            )
        else:
            current_w = float(self.left_container.width)
            self.left_container.width = max(250.0, current_w + delta_x)
        try:
            if self.app_page:
                self.update()
        except RuntimeError:
            pass

    def _on_v_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle vertical resizing of top-left container."""
        delta_y = getattr(e.local_delta, "y", 0.0) if e.local_delta else 0.0
        if self.top_left_container.height is None:
            page_h = (
                self.app_page.height
                if hasattr(self.app_page, "height")
                and isinstance(self.app_page.height, (int, float))
                else 600.0
            )
            self.top_left_container.expand = None
            self.bottom_left_container.expand = True
            self.top_left_container.height = float(page_h) * 0.5
        current_h = float(self.top_left_container.height or 240.0)
        self.top_left_container.height = max(110.0, current_h + delta_y)
        try:
            if self.app_page:
                self.update()
        except RuntimeError:
            pass

    def _show_notification(self, message: str) -> None:
        """Show a notification message in the snackbar.

        Args:
            message: Text message to display.
        """
        if not hasattr(self, "_notification_helper"):
            self._notification_helper = NotificationHelper(self.app_page)
        self._notification_helper.show_message(message)

    def _handle_add_primer(self, name: str, seq: str) -> None:
        """Add a primer from a designer card to the Input tab's primer list."""
        if not data_helpers.clean_sequence(seq):
            show_error_dialog(
                self.app_page,
                "Invalid Primer",
                "The primer sequence contains no valid nucleotides.",
            )
            return

        if self.on_add_primer:
            added, final_name = self.on_add_primer(name, seq)
        else:
            added, final_name = self.input_data.add_primer(name, seq)

        if final_name:
            if added:
                self._show_notification(
                    f"Added '{final_name}' to Input primer list."
                )
            else:
                self._show_notification(
                    f"Primer '{final_name}' is already in the Input primer "
                    "list."
                )

    def _on_background_task_done(self, task: asyncio.Task[Any]) -> None:
        """Discard a completed fallback task and log any unhandled error."""
        self._background_tasks.discard(task)
        if task.cancelled():
            return
        exc = task.exception()
        if exc is not None:
            logger.error(
                "Unhandled exception in background task %s",
                task.get_name(),
                exc_info=(type(exc), exc, exc.__traceback__),
            )

    def _schedule_on_event_loop(
        self, func: Callable[..., Coroutine[Any, Any, Any]], *args: Any
    ) -> bool:
        """Schedule an async function on the Flet page event loop.

        Marshals UI mutations made by background analysis threads onto the
        loop that owns the controls. Prefers the page-owned loop via
        ``page.run_task``. If the page cannot schedule tasks but a loop is
        already running on the current thread, the task is created there.
        When the page is detached and no loop exists, the callback is
        discarded: UI mutations must never run on a worker thread or a
        temporary event loop.

        Args:
            func: The async function to execute on the event loop.
            *args: Positional arguments passed to ``func``.

        Returns:
            True if the callback was scheduled, False if it was discarded.
        """
        page = self.app_page
        if page:
            try:
                page.run_task(func, *args)
                return True
            except (RuntimeError, TypeError, AttributeError):
                pass
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            logger.warning(
                "Discarded %s: no page event loop available",
                getattr(func, "__qualname__", repr(func)),
            )
            return False
        task = loop.create_task(func(*args))
        self._background_tasks.add(task)
        task.add_done_callback(self._on_background_task_done)
        return True

    def _start_analysis_worker(
        self,
        sync_worker: Callable[[], None],
        async_worker: Callable[[], Coroutine[Any, Any, None]],
        on_discarded: Callable[[], None] | None = None,
    ) -> None:
        """Run analysis in a background thread or cooperatively on the loop.

        Uses a daemon thread on desktop/server runtimes. On single-threaded
        WebAssembly/Pyodide runtimes (or if starting a thread fails with
        ``RuntimeError``), schedules ``async_worker`` on the Flet event loop.

        Args:
            sync_worker: Synchronous callable executed in a background thread.
            async_worker: Coroutine factory scheduled on the Flet event loop
                when threads are unavailable.
            on_discarded: Optional cleanup callback invoked if scheduling
                ``async_worker`` fails because no event loop is available.
        """
        if system_utils.can_use_background_threads():
            try:
                threading.Thread(target=sync_worker, daemon=True).start()
                return
            except RuntimeError as ex:
                logger.info(
                    "Background thread unavailable (%s); running analysis on "
                    "event loop",
                    ex,
                )
        scheduled = self._schedule_on_event_loop(async_worker)
        if scheduled is False and on_discarded is not None:
            on_discarded()

    def _get_analyse_button(self) -> Any:
        """Return the subclass form's analyse button control, if present."""
        form = getattr(self, "form", None)
        return getattr(form, "analyse_button", None)

    def _is_stale_run(self, cancel_event: threading.Event | None) -> bool:
        """Return True if ``cancel_event`` no longer matches the active run."""
        return (
            cancel_event is not None and self._cancel_event is not cancel_event
        )

    def _guard_not_running(self, action_message: str) -> bool:
        """Notify and return True if an analysis is currently running."""
        if self._analysis_running:
            self._show_notification(action_message)
            return True
        return False

    def _set_button_abort_mode(self, abort_mode: bool) -> None:
        """Toggle the analyse button between Analyse and Abort appearance."""
        button = self._get_analyse_button()
        if button is None:
            return
        if abort_mode:
            button.content = "Abort"
            button.icon = ft.Icons.STOP
            button.tooltip = "Stop analysis and keep results so far"
        else:
            button.content = "Analyse"
            button.icon = ft.Icons.PLAY_ARROW
            button.tooltip = "Run Primer Truncation Analysis"

    def _restore_loading_ui_if_active(self) -> None:
        """Restore the results panel if a loading progress bar is still shown.

        Subclasses override this to tear down their specific loading indicator.
        """

    def _finish_analysis_ui(
        self, cancel_event: threading.Event | None = None
    ) -> None:
        """Restore loading UI, button mode, and page state after analysis."""
        if self._is_stale_run(cancel_event):
            return
        self._restore_loading_ui_if_active()
        self._analysis_running = False
        button = self._get_analyse_button()
        if button is not None:
            button.disabled = False
        self._set_button_abort_mode(False)
        self._cancel_event = None
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    async def _on_analysis_finished(
        self, cancel_event: threading.Event | None = None
    ) -> None:
        """Restore the analyse button and flush the page after analysis."""
        self._finish_analysis_ui(cancel_event)

    def _bring_card_to_top_or_add(
        self, card_id: str, create_card_fn: Callable[[], ft.Card]
    ) -> None:
        """Bring existing card to the top of the list or create and prepend it.

        Args:
            card_id: Unique identifier string for the card.
            create_card_fn: Callable factory creating a new card if absent.
        """
        for ctrl in list(self.right_cards_list.controls):
            if getattr(ctrl, "_card_id", None) == card_id:
                self.right_cards_list.controls.remove(ctrl)
                self.right_cards_list.controls.insert(0, ctrl)
                if ctrl in self._active_cards:
                    self._active_cards.remove(ctrl)
                    self._active_cards.insert(0, ctrl)
                self._update_cards_header_visibility()
                try:
                    if self.app_page:
                        self.app_page.update()
                except RuntimeError:
                    pass
                return

        new_card = create_card_fn()
        self._active_cards.insert(0, new_card)
        self.right_cards_list.controls.insert(0, new_card)
        self._update_cards_header_visibility()
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    def _dismiss_card(self, card: ft.Card) -> None:
        """Remove a card from the right detail panel.

        Args:
            card: Card control to dismiss.
        """
        if card in self._active_cards:
            self._active_cards.remove(card)
        if card in self.right_cards_list.controls:
            self.right_cards_list.controls.remove(card)
        self._update_cards_header_visibility()
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    def _clear_all_cards(self, e: Any = None) -> None:
        """Clear all right-hand panel cards."""
        self._active_cards.clear()
        self.right_cards_list.controls.clear()
        self._update_cards_header_visibility()
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    def _update_cards_header_visibility(self) -> None:
        """Update visibility of clear button based on card count."""
        has_cards = len(self.right_cards_list.controls) > 0
        self.clear_cards_button.visible = has_cards

    def _update_panel_colours(self) -> None:
        """Refresh panel borders and split dividers from active theme."""
        self.top_left_container.border = ft.Border.all(
            1, GUIColours.OUTLINE_VARIANT
        )
        self.bottom_left_container.border = ft.Border.all(
            1, GUIColours.OUTLINE_VARIANT
        )
        if isinstance(self.main_h_divider.content, ft.Container):
            self.main_h_divider.content.bgcolor = GUIColours.DIVIDER_GREY
        if isinstance(self.left_v_divider.content, ft.Container):
            self.left_v_divider.content.bgcolor = GUIColours.DIVIDER_GREY

    async def _save_parameters_yaml(
        self,
        dialog_title: str,
        file_name: str,
        params: dict[str, Any],
    ) -> None:
        """Save a dictionary of parameters to a YAML file via file picker.

        Args:
            dialog_title: Dialog window title.
            file_name: Default file name.
            params: Dictionary of parameters to serialize.
        """
        yaml_str = safe_dump_yaml(params, sort_keys=False)
        await data_helpers.save_and_write_file(
            page=self.app_page,
            dialog_title=dialog_title,
            file_name=file_name,
            allowed_extensions=["yaml", "yml"],
            content=yaml_str,
            show_notification=self._show_notification,
            success_message_desktop="Parameters saved successfully.",
            success_message_web="Parameters ready for download!",
        )

    async def _load_parameters_yaml(
        self,
        dialog_title: str,
    ) -> dict[str, Any] | None:
        """Prompt user for a YAML parameter file and parse it.

        Args:
            dialog_title: Dialog window title.

        Returns:
            Parsed dictionary of parameters, or None if cancelled/invalid.
        """
        content = await data_helpers.pick_and_read_file(
            page=self.app_page,
            dialog_title=dialog_title,
            allowed_extensions=["yaml", "yml"],
            show_notification=self._show_notification,
        )
        if content is None:
            return None

        try:
            params = safe_load_yaml(content)
            if not isinstance(params, dict):
                self._show_notification("Error: Invalid parameter file format.")
                return None
            return params
        except Exception as ex:
            self._show_notification(f"Error parsing parameter file: {ex}")
            return None
