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

"""2D Primer Designer View for the Flet application."""

import logging
from typing import cast

import flet as ft

from amplifyp.dimer import PrimerDimerGenerator
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.gui_helpers import show_error_dialog
from amplifyp.gui.views.designer_2d.designer_2d_form import Designer2DForm
from amplifyp.gui.views.designer_2d.dismissible_2d_card import Dismissible2DCard
from amplifyp.gui.views.designer_2d.grid_2d_results import Grid2DResultsView
from amplifyp.primer_designer_2d import PrimerDesigner2D, PrimerDimers2D

logger = logging.getLogger(__name__)


class Designer2DView(ft.Row):  # type: ignore[misc]
    """2D Primer Designer view with resizable panels and cards."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
    ) -> None:
        """Initialise the Designer2DView."""
        super().__init__(expand=True, spacing=0)
        self.app_page = page
        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()
        self._cached_designer: PrimerDesigner2D | None = None
        self._active_cards: list[Dismissible2DCard] = []

        # Form component for 2D input controls and parameters
        self.form = Designer2DForm(
            settings=self.settings,
            on_submit_callback=self._run_designer_event,
        )

        # Top-left container (50% default vertical height)
        self.top_left_container = ft.Container(
            content=self.form,
            expand=1,
            padding=10,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            border_radius=5,
        )

        # Vertical divider resizer for left panel (top vs bottom)
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

        # Grid view component for bottom-left container
        self.results_grid = Grid2DResultsView(
            settings=self.settings,
            on_select_step_callback=self._on_grid_step_selected,
        )

        # Bottom-left container (50% default vertical height)
        self.bottom_left_container = ft.Container(
            content=self.results_grid,
            expand=1,
            padding=10,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            border_radius=5,
        )

        # Left main panel column (50% top, 50% bottom)
        self.left_panel_column = ft.Column(
            [
                self.top_left_container,
                self.left_v_divider,
                self.bottom_left_container,
            ],
            spacing=0,
            expand=True,
        )
        self.left_container = ft.Container(
            content=self.left_panel_column,
            expand=1,  # 50% default viewport width
            padding=5,
        )

        # Horizontal divider resizer between left and right panel
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

        # Right panel for dismissible cards
        self.right_cards_list = ft.ListView(
            expand=True, spacing=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.clear_cards_button = ft.TextButton(
            "Clear Cards",
            icon=ft.Icons.DELETE_SWEEP,
            tooltip="Clear All 2D Pair Cards",
            on_click=self._clear_all_cards,
            visible=False,
        )
        self.right_header = ft.Row(
            [
                ft.Text(
                    "2D Primer Pair Detail Cards",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_header", 18),
                ),
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
            expand=1,  # 50% default viewport width
            padding=10,
        )

        # Assembly into main Row controls
        self.controls = [
            self.left_container,
            self.main_h_divider,
            self.right_container,
        ]

    def _on_h_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle horizontal resizing between left and right panels."""
        delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
        if self.left_container.width is None:
            page_w = (
                self.app_page.width
                if hasattr(self.app_page, "width")
                and isinstance(self.app_page.width, (int, float))
                else 800.0
            )
            base_w = float(page_w) * 0.5
            self.left_container.expand = None
            self.right_container.expand = 1
            self.left_container.width = base_w
        current_w = float(self.left_container.width or 300.0)
        self.left_container.width = max(250.0, current_w + delta_x)
        try:
            if self.app_page:
                self.update()
        except RuntimeError:
            pass

    def _on_v_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle vertical resizing between top-left and bottom-left panels."""
        delta_y = getattr(e.local_delta, "y", 0.0) if e.local_delta else 0.0
        if self.top_left_container.height is None:
            page_h = (
                self.app_page.height
                if hasattr(self.app_page, "height")
                and isinstance(self.app_page.height, (int, float))
                else 600.0
            )
            base_h = float(page_h) * 0.5
            self.top_left_container.expand = None
            self.bottom_left_container.expand = True
            self.top_left_container.height = base_h
        current_h = float(self.top_left_container.height or 300.0)
        self.top_left_container.height = max(150.0, current_h + delta_y)
        try:
            if self.app_page:
                self.update()
        except RuntimeError:
            pass

    def _run_designer_event(self) -> None:
        """Run 2D primer truncation analysis based on form inputs."""
        try:
            (
                fwd_dna,
                fwd_min_len,
                rev_dna,
                rev_min_len,
                threshold,
                max_overlap,
                filter_metric,
            ) = self.form.validate_and_get_params()
        except ValueError:
            return

        pd_settings = self.settings.get_primer_dimer_settings()
        generator = PrimerDimerGenerator(settings=pd_settings)

        try:
            designer = PrimerDesigner2D(
                fwd_dna=fwd_dna,
                fwd_min_length=fwd_min_len,
                rev_dna=rev_dna,
                rev_min_length=rev_min_len,
                generator=generator,
                threshold=threshold,
                max_overlap=max_overlap,
                filter_metric=filter_metric,
            )
            self._cached_designer = designer
            self.results_grid.update_grid(designer)
        except Exception as ex:
            logger.exception("Failed to run 2D primer designer")
            show_error_dialog(
                self.app_page,
                "Analysis Error",
                f"Error performing 2D primer design: {ex}",
            )

    def _on_grid_step_selected(self, step: PrimerDimers2D) -> None:
        """Handle selection of a step from the 2D grid results.

        Args:
            step: Selected PrimerDimers2D instance.
        """
        fwd_len = len(step.fwd_fwd.primer_1.seq)
        rev_len = len(step.rev_rev.primer_1.seq)
        card_id = f"card_2d_{fwd_len}_{rev_len}"

        # Check if card already exists in list
        for existing in self._active_cards:
            if existing._card_id == card_id:
                return

        card = Dismissible2DCard(
            card_id=card_id,
            step=step,
            settings=self.settings,
            dismiss_callback=self._dismiss_card,
            font_family=self.settings.get("font_family", "Roboto Mono"),
        )

        self._active_cards.insert(0, card)
        self.right_cards_list.controls.insert(0, card)
        self.clear_cards_button.visible = True

        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    def _dismiss_card(self, card: ft.Card) -> None:
        """Dismiss a card from the right panel.

        Args:
            card: Card control to dismiss.
        """
        card_2d = cast(Dismissible2DCard, card)
        if card_2d in self._active_cards:
            self._active_cards.remove(card_2d)
        if card_2d in self.right_cards_list.controls:
            self.right_cards_list.controls.remove(card_2d)

        self.clear_cards_button.visible = len(self._active_cards) > 0

        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    def _clear_all_cards(
        self, e: ft.Event[ft.TextButton] | None = None
    ) -> None:
        """Clear all active 2D detail cards."""
        self._active_cards.clear()
        self.right_cards_list.controls.clear()
        self.clear_cards_button.visible = False
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass
