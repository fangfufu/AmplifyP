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

"""1D Primer Designer View for the Flet application."""

import logging
import traceback

import flet as ft

from amplifyp.dimer import PrimerDimer, PrimerDimerGenerator
from amplifyp.dna import DNA
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.gui_helpers import show_error_dialog
from amplifyp.gui.views.designer_1d.designer_1d_form import Designer1DForm
from amplifyp.gui.views.designer_1d.dismissible_self_dimer_card import (
    DismissibleSelfDimerCard,
)
from amplifyp.gui.views.designer_1d.primer_item_card import PrimerItemCard
from amplifyp.gui.views.designer_1d.quality_bar_chart import QualityBarChart
from amplifyp.primer_designer_1d import PrimerDesigner1D

logger = logging.getLogger(__name__)


class PrimerDesignerView(ft.Row):  # type: ignore[misc]
    """1D Primer Designer view with resizable panels and cards."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
    ) -> None:
        """Initialise the PrimerDesignerView."""
        super().__init__(expand=True, spacing=0)
        self.app_page = page
        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()
        self._cached_designer: PrimerDesigner1D | None = None

        # Form component for input controls and parameters
        self.form = Designer1DForm(
            settings=self.settings,
            on_submit_callback=self._run_designer_event,
        )

        # Container for top-left controls
        self.top_left_container = ft.Container(
            content=self.form,
            height=240,  # 1/3rd default vertical space
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

        # Bottom-left output list of primers
        self.primer_list = ft.ListView(
            expand=True, spacing=6, scroll=ft.ScrollMode.ALWAYS
        )
        self.bottom_left_container = ft.Container(
            content=ft.Column(
                [
                    ft.Text(
                        "Generated Primers",
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_subheader", 16),
                    ),
                    ft.Container(content=self.primer_list, expand=True),
                ],
                spacing=6,
            ),
            expand=True,
            padding=10,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            border_radius=5,
        )

        # Left main panel container
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
            expand=True,
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

        # Top-right panel: Vertical Quality Bar Chart
        self.chart_content_container = ft.Container(
            expand=True,
        )
        self.top_right_chart_container = ft.Container(
            content=ft.Column(
                [
                    ft.Container(
                        content=ft.Text(
                            "Self-Dimer Quality by Primer Size (bp)",
                            weight=ft.FontWeight.BOLD,
                            size=self.settings.get("font_size_subheader", 16),
                        ),
                        padding=ft.Padding(10, 10, 10, 0),
                    ),
                    self.chart_content_container,
                ],
                spacing=4,
            ),
            height=240,  # 1/3rd default vertical space
            padding=0,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            border_radius=5,
        )
        self.chart_content_container.content = self._build_chart([])

        # Vertical divider resizer between top-right and bottom-right
        self.right_v_divider = ft.GestureDetector(
            on_pan_update=self._on_right_v_pan_update,
            content=ft.Container(
                height=6,
                bgcolor=GUIColours.DIVIDER_GREY,
                border_radius=3,
                margin=ft.Margin.symmetric(vertical=4),
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_UP_DOWN,
        )

        # Bottom-right panel for observing primer self-dimer cards
        self.right_cards_list = ft.ListView(
            expand=True, spacing=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.clear_cards_button = ft.TextButton(
            "Clear Cards",
            icon=ft.Icons.DELETE_SWEEP,
            tooltip="Clear All Self-Dimer Cards",
            on_click=self._clear_all_cards,
            visible=False,
        )
        self.right_header = ft.Row(
            [
                ft.Text(
                    "Self-Dimer Cards",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_header", 18),
                ),
                self.clear_cards_button,
            ],
            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
        )
        self.bottom_right_container = ft.Container(
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

        self.right_panel_column = ft.Column(
            [
                self.top_right_chart_container,
                self.right_v_divider,
                self.bottom_right_container,
            ],
            spacing=0,
            expand=True,
        )
        self.right_container = ft.Container(
            content=self.right_panel_column,
            expand=True,
            padding=5,
        )

        # Assembly into main Row controls
        self.controls = [
            self.left_container,
            self.main_h_divider,
            self.right_container,
        ]

    # --- Field property aliases to preserve backwards compatibility ---
    @property
    def dna_input(self) -> ft.TextField:
        """Get the DNA sequence input field."""
        return self.form.dna_input

    @property
    def min_len_input(self) -> ft.TextField:
        """Get the minimum length input field."""
        return self.form.min_len_input

    @property
    def mode_dropdown(self) -> ft.Dropdown:
        """Get the direction mode dropdown."""
        return self.form.mode_dropdown

    @property
    def max_quality_input(self) -> ft.TextField:
        """Get the max quality input field."""
        return self.form.max_quality_input

    @property
    def max_overlap_input(self) -> ft.TextField:
        """Get the max overlap input field."""
        return self.form.max_overlap_input

    @property
    def analyse_button(self) -> ft.FilledButton:
        """Get the analyse button control."""
        return self.form.analyse_button

    @property
    def error_text(self) -> ft.Text:
        """Get the error display text control."""
        return self.form.error_text

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
        """Handle vertical resizing of top-left control panel."""
        delta_y = getattr(e.local_delta, "y", 0.0) if e.local_delta else 0.0
        current_h = float(self.top_left_container.height or 160.0)
        self.top_left_container.height = max(110.0, current_h + delta_y)
        try:
            if self.app_page:
                self.update()
        except RuntimeError:
            pass

    def _on_right_v_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle vertical resizing of top-right chart panel."""
        delta_y = getattr(e.local_delta, "y", 0.0) if e.local_delta else 0.0
        current_h = float(self.top_right_chart_container.height or 140.0)
        self.top_right_chart_container.height = max(70.0, current_h + delta_y)
        if self._cached_designer and self._cached_designer.all_dimers:
            self.chart_content_container.content = self._build_chart(
                list(self._cached_designer.all_dimers)
            )
        try:
            if self.app_page:
                self.update()
        except RuntimeError:
            pass

    def _build_chart(self, dimers: list[PrimerDimer]) -> ft.Control:
        """Helper to build quality chart with container height."""
        container_h = float(self.top_right_chart_container.height or 140.0)
        return QualityBarChart.build_chart(
            dimers=dimers,
            container_height=container_h,
            on_primer_selected=self._on_primer_selected,
        )

    def _run_designer_event(self, e: ft.ControlEvent | None = None) -> None:
        """Event handler wrapper for running analysis."""
        self.run_designer()

    def run_designer(self) -> bool:
        """Validate inputs, run 1D primer design analysis, and update UI."""
        params = self.form.validate_and_get_params()
        if params is None:
            return False

        clean_seq, min_length, mode, threshold, max_overlap = params
        self.primer_list.controls.clear()

        try:
            dna_obj = DNA(clean_seq)
            pd_settings = self.settings.get_primer_dimer_settings()
            generator = PrimerDimerGenerator(settings=pd_settings)
            designer = PrimerDesigner1D(
                dna=dna_obj,
                min_length=min_length,
                mode=mode,
                generator=generator,
                threshold=threshold,
                max_overlap=max_overlap,
            )
            self._cached_designer = designer

            # Update top-right quality bar chart
            self.chart_content_container.content = self._build_chart(
                list(designer.all_dimers)
            )

            for step_idx, dimer in enumerate(designer.all_dimers):
                item_card = PrimerItemCard(
                    dimer=dimer,
                    step_index=step_idx,
                    mode=mode,
                    settings=self.settings,
                    on_select_callback=self._on_primer_selected,
                )
                self.primer_list.controls.append(item_card)

        except (ValueError, RuntimeError, OSError) as ex:
            logger.exception("1D Primer Design failed: %s", ex)
            self.form.show_error(f"Error: {ex}")
            show_error_dialog(
                self.app_page,
                "Error running Primer Designer",
                f"{ex}\n{traceback.format_exc()}",
            )
            self.app_page.update()
            return False

        self.app_page.update()
        return True

    def _on_primer_selected(self, dimer: PrimerDimer, step_index: int) -> None:
        """Handle primer selection: add or raise self-dimer card."""
        card_id = f"1d_dimer_{dimer.primer_1.seq}_{step_index}"

        # Check if card already exists -> bring to top
        for ctrl in self.right_cards_list.controls:
            if getattr(ctrl, "_card_id", None) == card_id:
                self.right_cards_list.controls.remove(ctrl)
                self.right_cards_list.controls.insert(0, ctrl)
                self._update_cards_header_visibility()
                self.app_page.update()
                return

        # Create new dismissible self-dimer card
        font_family = self.settings.get("font_family", "Roboto Mono")
        card = DismissibleSelfDimerCard(
            card_id=card_id,
            dimer=dimer,
            settings=self.settings,
            dismiss_callback=self._dismiss_card,
            font_family=font_family,
            step_index=step_index,
        )
        self.right_cards_list.controls.insert(0, card)
        self._update_cards_header_visibility()
        self.app_page.update()

    def _dismiss_card(self, card: ft.Card) -> None:
        """Callback to remove a dismissed self-dimer card."""
        if card in self.right_cards_list.controls:
            self.right_cards_list.controls.remove(card)
            self._update_cards_header_visibility()
            self.app_page.update()

    def _update_cards_header_visibility(self) -> None:
        """Update visibility of clear button based on card count."""
        has_cards = len(self.right_cards_list.controls) > 0
        self.clear_cards_button.visible = has_cards

    def _clear_all_cards(
        self, e: ft.Event[ft.TextButton] | None = None
    ) -> None:
        """Clear all right-hand panel self-dimer cards."""
        self.right_cards_list.controls.clear()
        self._update_cards_header_visibility()
        self.app_page.update()
