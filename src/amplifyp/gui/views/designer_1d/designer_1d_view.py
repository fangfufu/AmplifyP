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

from __future__ import annotations

import logging
import traceback

import flet as ft

from amplifyp.dimer import PrimerDimer, PrimerDimerGenerator
from amplifyp.dna import DNA
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.data_helpers import clean_sequence
from amplifyp.gui.utils.gui_helpers import show_error_dialog
from amplifyp.gui.views.designer.designer_view_base import BaseDesignerView
from amplifyp.gui.views.designer_1d.designer_1d_form import Designer1DForm
from amplifyp.gui.views.designer_1d.dismissible_self_dimer_card import (
    DismissibleSelfDimerCard,
)
from amplifyp.gui.views.designer_1d.primer_item_card import PrimerItemCard
from amplifyp.gui.views.designer_1d.quality_bar_chart import QualityBarChart
from amplifyp.primer_designer_1d import PrimerDesigner1D

logger = logging.getLogger(__name__)


class PrimerDesignerView(BaseDesignerView):
    """1D Primer Designer view with resizable panels and cards."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
    ) -> None:
        """Initialise the PrimerDesignerView."""
        super().__init__(page=page, input_data=input_data, settings=settings)
        self._cached_designer: PrimerDesigner1D | None = None

        # Form component for input controls and parameters
        self.form = Designer1DForm(
            settings=self.settings,
            on_submit_callback=self._run_designer_event,
            on_save_callback=self._save_designer_1d_click,
            on_load_callback=self._load_designer_1d_click,
            on_clear_all_callback=self._clear_all,
        )

        # Top-left container for form
        self.top_left_container = ft.Container(
            content=self.form,
            height=240,
            padding=10,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            border_radius=5,
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

        # Top-right panel: Vertical Quality Bar Chart
        self.chart_content_container = ft.Container(expand=True)
        self.top_right_chart_container = ft.Container(
            content=ft.Column(
                [
                    ft.Container(
                        content=ft.Text(
                            "Self-Dimer Quality by Primer Size (nt)",
                            weight=ft.FontWeight.BOLD,
                            size=self.settings.get("font_size_subheader", 16),
                        ),
                        padding=ft.Padding(10, 10, 10, 0),
                    ),
                    self.chart_content_container,
                ],
                spacing=4,
            ),
            height=240,
            padding=0,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            border_radius=5,
        )
        self.chart_content_container.content = self._build_chart([])

        # Right-side vertical divider resizer
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

        # Customise right-hand panel header
        self.right_title.value = "Self-Dimer Cards"
        self.clear_cards_button.tooltip = "Clear All Self-Dimer Cards"

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

    # --- Property accessors for backward compatibility and test access ---
    @property
    def dna_input(self) -> ft.TextField:
        """Get the DNA sequence input field."""
        return self.form.dna_input

    @property
    def length_display(self) -> ft.TextField:
        """Get the length display field."""
        return self.form.length_display

    @property
    def min_len_input(self) -> ft.TextField:
        """Get the minimum length input field."""
        return self.form.min_len_input

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
    def clear_all_button(self) -> ft.FilledTonalButton:
        """Get the clear all button control."""
        return self.form.clear_all_button

    @property
    def error_text(self) -> ft.Text:
        """Get the error display text control."""
        return self.form.error_text

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

        def _factory() -> DismissibleSelfDimerCard:
            font_family = self.settings.get("font_family", "Roboto Mono")
            return DismissibleSelfDimerCard(
                card_id=card_id,
                dimer=dimer,
                settings=self.settings,
                dismiss_callback=self._dismiss_card,
                font_family=font_family,
                step_index=step_index,
            )

        self._bring_card_to_top_or_add(card_id, _factory)

    def _clear_all(self, e: ft.ControlEvent | None = None) -> None:
        """Clear all inputs, parameters, error messages, and results."""
        self.form.dna_input.value = ""
        self.form.length_display.value = "0"
        self.form.min_len_input.value = ""
        self.form.max_quality_input.value = ""
        self.form.max_overlap_input.value = ""
        self.form.clear_errors()
        self.primer_list.controls.clear()
        self._cached_designer = None
        self.chart_content_container.content = self._build_chart([])
        self._clear_all_cards()
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    async def _save_designer_1d_click(self, e: ft.ControlEvent) -> None:
        """Save Designer 1D parameters to a YAML file."""
        params = {
            "dna": (self.form.dna_input.value or ""),
            "min_length": (self.form.min_len_input.value or ""),
            "max_quality": (self.form.max_quality_input.value or ""),
            "max_overlap": (self.form.max_overlap_input.value or ""),
        }
        await self._save_parameters_yaml(
            dialog_title="Save Designer 1D Parameters",
            file_name="designer_1d_parameters.yaml",
            params=params,
        )

    async def _load_designer_1d_click(self, e: ft.ControlEvent) -> None:
        """Load Designer 1D parameters from a YAML file."""
        params = await self._load_parameters_yaml(
            dialog_title="Load Designer 1D Parameters"
        )
        if params is None:
            return

        dna_val = params.get("dna")
        dna_str = str(dna_val) if dna_val is not None else ""
        self.form.dna_input.value = dna_str
        self.form.length_display.value = str(len(clean_sequence(dna_str)))

        min_len_val = params.get("min_length")
        self.form.min_len_input.value = (
            str(min_len_val) if min_len_val is not None else ""
        )

        max_q_val = params.get("max_quality")
        self.form.max_quality_input.value = (
            str(max_q_val) if max_q_val is not None else ""
        )
        max_ov_val = params.get("max_overlap")
        self.form.max_overlap_input.value = (
            str(max_ov_val) if max_ov_val is not None else ""
        )

        self.form.clear_errors()
        self.app_page.update()

        self._show_notification("Parameters loaded successfully.")
        self.run_designer()
