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

from collections.abc import Callable
from typing import Any

import flet as ft
import yaml

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils import data_helpers
from amplifyp.gui.utils.gui_helpers import NotificationHelper


class BaseDesignerView(ft.Row):  # type: ignore[misc]
    """Base view providing split layout, resizers, and card manager."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
    ) -> None:
        """Initialise BaseDesignerView.

        Args:
            page: Flet page instance.
            input_data: Application GUIInput instance.
            settings: Application GUISettings instance.
        """
        super().__init__(expand=True, spacing=0)
        self.app_page = page
        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()
        self._active_cards: list[ft.Card] = []

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
        yaml_str = yaml.safe_dump(params, sort_keys=False)
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
            params = yaml.safe_load(content)
            if not isinstance(params, dict):
                self._show_notification("Error: Invalid parameter file format.")
                return None
            return params
        except Exception as ex:
            self._show_notification(f"Error parsing parameter file: {ex}")
            return None
