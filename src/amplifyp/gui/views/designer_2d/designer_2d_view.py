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

from __future__ import annotations

import logging

import flet as ft

from amplifyp.dimer import PrimerDimerGenerator
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.gui_helpers import show_error_dialog
from amplifyp.gui.views.designer.designer_view_base import BaseDesignerView
from amplifyp.gui.views.designer_2d.designer_2d_form import Designer2DForm
from amplifyp.gui.views.designer_2d.dismissible_2d_card import Dismissible2DCard
from amplifyp.gui.views.designer_2d.grid_2d_results import Grid2DResultsView
from amplifyp.primer_designer_2d import PrimerDesigner2D, PrimerDimers2D

logger = logging.getLogger(__name__)


class Designer2DView(BaseDesignerView):
    """2D Primer Designer view with resizable panels and cards."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
    ) -> None:
        """Initialise the Designer2DView."""
        super().__init__(page=page, input_data=input_data, settings=settings)
        self._cached_designer: PrimerDesigner2D | None = None

        # Form component for 2D input controls and parameters
        self.form = Designer2DForm(
            settings=self.settings,
            on_submit_callback=self._run_designer_event,
            on_save_callback=self._save_designer_2d_click,
            on_load_callback=self._load_designer_2d_click,
            on_clear_all_callback=self._clear_all,
        )

        # Top-left container (50% default vertical height)
        self.top_left_container = ft.Container(
            content=self.form,
            expand=1,
            padding=10,
            border=ft.Border.all(1, GUIColours.OUTLINE_VARIANT),
            border_radius=5,
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
            expand=1,
            padding=5,
        )

        # Customise right-hand panel header
        self.right_title.value = "2D Primer Pair Detail Cards"
        self.clear_cards_button.tooltip = "Clear All 2D Pair Cards"

        self.right_container = ft.Container(
            content=ft.Column(
                [
                    self.right_header,
                    ft.Container(content=self.right_cards_list, expand=True),
                ],
                spacing=8,
            ),
            expand=1,
            padding=10,
        )

        # Assembly into main Row controls
        self.controls = [
            self.left_container,
            self.main_h_divider,
            self.right_container,
        ]

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
            self._clear_all_cards()
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

        def _factory() -> Dismissible2DCard:
            return Dismissible2DCard(
                card_id=card_id,
                step=step,
                settings=self.settings,
                dismiss_callback=self._dismiss_card,
                font_family=self.settings.get("font_family", "Roboto Mono"),
            )

        self._bring_card_to_top_or_add(card_id, _factory)

    def update_ui(self) -> None:
        """Update 2D designer view components when settings change."""
        if self._cached_designer is not None:
            self.results_grid.update_grid(self._cached_designer)
        for card in self._active_cards:
            if isinstance(card, Dismissible2DCard):
                card.on_settings_change()
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    def _clear_all(self, e: ft.ControlEvent | None = None) -> None:
        """Clear inputs, parameters, error messages, grid, and cards."""
        self.form.fwd_dna_input.value = ""
        self.form.fwd_min_len_input.value = ""
        self.form.rev_dna_input.value = ""
        self.form.rev_min_len_input.value = ""
        self.form.max_quality_input.value = ""
        self.form.max_overlap_input.value = ""
        self.form.clear_errors()
        self._cached_designer = None
        self.results_grid.clear_grid()
        self._clear_all_cards()
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    async def _save_designer_2d_click(self, e: ft.ControlEvent) -> None:
        """Save Designer 2D parameters to a YAML file."""
        params = {
            "fwd_dna": (self.form.fwd_dna_input.value or ""),
            "fwd_min_length": (self.form.fwd_min_len_input.value or ""),
            "rev_dna": (self.form.rev_dna_input.value or ""),
            "rev_min_length": (self.form.rev_min_len_input.value or ""),
            "max_quality": (self.form.max_quality_input.value or ""),
            "max_overlap": (self.form.max_overlap_input.value or ""),
        }
        await self._save_parameters_yaml(
            dialog_title="Save Designer 2D Parameters",
            file_name="designer_2d_parameters.yaml",
            params=params,
        )

    async def _load_designer_2d_click(self, e: ft.ControlEvent) -> None:
        """Load Designer 2D parameters from a YAML file."""
        params = await self._load_parameters_yaml(
            dialog_title="Load Designer 2D Parameters"
        )
        if params is None:
            return

        self.form.fwd_dna_input.value = str(params.get("fwd_dna", ""))
        fwd_min_val = params.get("fwd_min_length")
        self.form.fwd_min_len_input.value = (
            str(fwd_min_val) if fwd_min_val is not None else ""
        )
        self.form.rev_dna_input.value = str(params.get("rev_dna", ""))
        rev_min_val = params.get("rev_min_length")
        self.form.rev_min_len_input.value = (
            str(rev_min_val) if rev_min_val is not None else ""
        )

        # Support both new max_quality and legacy quality_filter
        max_q_val = params.get("max_quality", params.get("quality_filter"))
        self.form.max_quality_input.value = (
            str(max_q_val) if max_q_val is not None else ""
        )

        # Support both new max_overlap and legacy overlap_filter
        max_ov_val = params.get("max_overlap", params.get("overlap_filter"))
        self.form.max_overlap_input.value = (
            str(max_ov_val) if max_ov_val is not None else ""
        )

        self.form.clear_errors()
        self.app_page.update()

        self._show_notification("Parameters loaded successfully.")
        self._run_designer_event()
