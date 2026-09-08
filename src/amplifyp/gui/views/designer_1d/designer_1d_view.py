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
import threading
import traceback
from collections.abc import Callable

import flet as ft

from amplifyp.dimer import PrimerDimer, PrimerDimerGenerator
from amplifyp.dna import DNA, DNADirection, DNAType
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.data_helpers import clean_sequence
from amplifyp.gui.utils.gui_helpers import BorderedCheckbox, show_error_dialog
from amplifyp.gui.views.designer.designer_view_base import BaseDesignerView
from amplifyp.gui.views.designer.progress_tracker import ProgressTracker
from amplifyp.gui.views.designer_1d.designer_1d_form import Designer1DForm
from amplifyp.gui.views.designer_1d.dismissible_self_dimer_card import (
    DismissibleSelfDimerCard,
)
from amplifyp.gui.views.designer_1d.primer_item_card import PrimerItemCard
from amplifyp.gui.views.designer_1d.quality_bar_chart import QualityBarChart
from amplifyp.primer_designer_1d import PrimerDesigner1D
from amplifyp.repliconf import Repliconf

logger = logging.getLogger(__name__)


class PrimerDesignerView(BaseDesignerView):
    """1D Primer Designer view with resizable panels and cards."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
        on_run_pcr: Callable[[str, str], None] | None = None,
    ) -> None:
        """Initialise the PrimerDesignerView."""
        super().__init__(page=page, input_data=input_data, settings=settings)
        self.on_run_pcr = on_run_pcr
        self._cached_designer: PrimerDesigner1D | None = None
        self._analysis_running = False
        self._cancel_event: threading.Event | None = None

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

        # Bottom-left output: progress bar (while loading) + primer list
        self.progress_tracker = ProgressTracker(
            page=page,
            settings=self.settings,
            hint_text="Analysing primer truncations\u2026",
        )
        self.primer_list = ft.ListView(
            expand=True, spacing=6, scroll=ft.ScrollMode.ALWAYS
        )
        self._primer_list_header = ft.Text(
            "Generated Primers",
            weight=ft.FontWeight.BOLD,
            size=self.settings.get("font_size_subheader", 16),
        )
        self._primer_list_body = ft.Container(
            content=self.primer_list, expand=True
        )
        self.bottom_left_container = ft.Container(
            content=ft.Column(
                [
                    self._primer_list_header,
                    self._primer_list_body,
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
    def check_template_checkbox(self) -> BorderedCheckbox | ft.Checkbox:
        """Get the check against template checkbox control."""
        return self.form.filter_dna_checkbox

    @property
    def filter_dna_checkbox(self) -> BorderedCheckbox | ft.Checkbox:
        """Get the check against template checkbox control."""
        return self.form.filter_dna_checkbox

    @property
    def check_dna_checkbox(self) -> BorderedCheckbox | ft.Checkbox:
        """Get the check against template checkbox control."""
        return self.form.filter_dna_checkbox

    @property
    def max_binding_sites_input(self) -> ft.TextField:
        """Get the max binding sites input field."""
        return self.form.max_binding_sites_input

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
        """Event handler wrapper for running analysis or aborting it.

        While an analysis is running, the button acts as an abort button and
        the click requests a stop instead of starting a new run.
        """
        if self._analysis_running and self._cancel_event is not None:
            self._cancel_event.set()
            return
        self._start_designer()

    def _set_button_abort_mode(self, abort_mode: bool) -> None:
        """Toggle the analyse button between Analyse and Abort appearance."""
        button = self.form.analyse_button
        if abort_mode:
            button.content = "Abort"
            button.icon = ft.Icons.STOP
            button.tooltip = "Stop analysis and keep results so far"
        else:
            button.content = "Analyse"
            button.icon = ft.Icons.PLAY_ARROW
            button.tooltip = "Run Primer Truncation Analysis"

    def show_loading(self, total: int = 0) -> None:
        """Replace primer list with a progress bar while analysis runs.

        The tracker schedules a flush task on the Flet event loop that
        calls ``app_page.update()`` every 50 ms so the bar animates
        smoothly regardless of analysis speed.

        Args:
            total: Total truncation steps. When 0, bar is indeterminate.
        """
        loading_body = self.progress_tracker.show(total)
        col = self.bottom_left_container.content
        if isinstance(col, ft.Column):
            col.controls = [self._primer_list_header, loading_body]
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    def update_progress(self, done: int, total: int) -> None:
        """Write current progress values; the flush task renders them.

        The bar value and label are updated on every call. The flush
        task runs on the Flet event loop at ~20 fps, so the analysis
        thread is never blocked by rendering. Stops the flush task when
        the final tick is received.

        Args:
            done: Number of truncation steps completed so far.
            total: Total number of truncation steps.
        """
        self.progress_tracker.update_progress(done, total)

    def _restore_primer_list(self) -> None:
        """Restore bottom-left panel to show the primer list."""
        self.progress_tracker.hide()
        col = self.bottom_left_container.content
        if isinstance(col, ft.Column):
            col.controls = [
                self._primer_list_header,
                self._primer_list_body,
            ]

    def _compute_origin_counts(
        self,
        designer: PrimerDesigner1D,
        template_dna: DNA | None,
    ) -> list[int | None]:
        """Count template binding sites per dimer.

        Runs in the analysis worker thread; performs no UI mutation.

        Args:
            designer: The completed 1D primer designer.
            template_dna: Template DNA for origin counting, or None.

        Returns:
            Per-dimer binding-site counts, or None entries when no template.
        """
        if template_dna is None:
            return [None] * len(designer.all_dimers)

        counts: list[int | None] = []
        for dimer in designer.all_dimers:
            repliconf = Repliconf(template_dna, dimer.primer_1)
            repliconf.search()
            counts.append(
                len(repliconf.origin_db.fwd) + len(repliconf.origin_db.rev)
            )
        return counts

    def _update_chart_and_primer_list(
        self,
        designer: PrimerDesigner1D,
        origin_counts: list[int | None],
        mode: DNADirection,
    ) -> None:
        """Update the quality chart and populate the primer list with cards.

        Must run on the Flet event loop.

        Args:
            designer: The completed 1D primer designer.
            origin_counts: Per-dimer binding-site counts, or None entries.
            mode: The truncation mode used for the design.
        """
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
                on_run_pcr_callback=self._handle_run_pcr,
                origin_count=origin_counts[step_idx],
            )
            self.primer_list.controls.append(item_card)

    async def _on_analysis_success(
        self,
        designer: PrimerDesigner1D,
        origin_counts: list[int | None],
        mode: DNADirection,
    ) -> None:
        """Populate the results UI on the event loop after analysis."""
        self._restore_primer_list()
        self._update_chart_and_primer_list(designer, origin_counts, mode)
        if designer.aborted:
            self._show_notification(
                "Analysis aborted — showing primers analysed so far."
            )

    async def _on_analysis_error(self, ex: Exception, tb: str) -> None:
        """Show the analysis failure UI on the event loop."""
        self.form.show_error(f"Error: {ex}")
        show_error_dialog(
            self.app_page,
            "Error running Primer Designer",
            f"{ex}\n{tb}",
        )
        self._restore_primer_list()

    async def _on_analysis_finished(
        self, cancel_event: threading.Event | None = None
    ) -> None:
        """Restore the analyse button and flush the page after analysis."""
        if cancel_event is not None and self._cancel_event is not cancel_event:
            return
        self.form.analyse_button.disabled = False
        self._set_button_abort_mode(False)
        self._cancel_event = None
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

    def _start_designer(self) -> None:
        """Validate inputs, show progress bar, and run analysis in a thread."""
        if self._analysis_running:
            return
        params = self.form.validate_and_get_params()
        if params is None:
            return

        (
            clean_seq,
            min_length,
            mode,
            threshold,
            max_overlap,
            filter_dna_enabled,
            max_binding_sites,
        ) = params
        self.primer_list.controls.clear()

        template_dna: DNA | None = None
        if filter_dna_enabled:
            clean_tpl = clean_sequence(self.input_data.template)
            if not clean_tpl:
                self.form.show_error(
                    "Template DNA sequence is required when check against "
                    "template is enabled. Please enter a template in the "
                    "Input view."
                )
                try:
                    if self.app_page:
                        self.app_page.update()
                except RuntimeError:
                    pass
                return
            t_type = (
                DNAType.CIRCULAR
                if self.input_data.template_circular
                else DNAType.LINEAR
            )
            template_dna = DNA(clean_tpl, dna_type=t_type)

        dna_obj = DNA(clean_seq)
        # Total truncation steps is deterministic before threading.
        total_steps = len(dna_obj.seq) - min_length + 1

        self.show_loading(total=total_steps)
        self._analysis_running = True
        cancel_event = threading.Event()
        self._cancel_event = cancel_event
        self._set_button_abort_mode(True)
        try:
            if self.app_page:
                self.app_page.update()
        except RuntimeError:
            pass

        pd_settings = self.settings.get_primer_dimer_settings()
        generator = PrimerDimerGenerator(settings=pd_settings)

        def _on_progress(done: int, total: int) -> None:
            """Forward every progress tick to the progress bar."""
            self.update_progress(done, total)

        def _run_analysis() -> None:
            """Execute 1D analysis in a background thread.

            Computation (PrimerDesigner1D and Repliconf origin counting)
            stays in the worker thread; UI mutations are marshalled onto
            the Flet event loop.
            """
            try:
                designer = PrimerDesigner1D(
                    dna=dna_obj,
                    min_length=min_length,
                    mode=mode,
                    generator=generator,
                    threshold=threshold,
                    max_overlap=max_overlap,
                    template=template_dna,
                    max_origin_count=max_binding_sites,
                    on_progress=_on_progress,
                    cancel_event=cancel_event,
                )
                self._cached_designer = designer
                origin_counts = self._compute_origin_counts(
                    designer, template_dna
                )
                self._schedule_on_event_loop(
                    self._on_analysis_success, designer, origin_counts, mode
                )
            except (ValueError, RuntimeError, OSError) as ex:
                logger.exception("1D Primer Design failed: %s", ex)
                self._schedule_on_event_loop(
                    self._on_analysis_error, ex, traceback.format_exc()
                )
            finally:
                self._analysis_running = False
                self._schedule_on_event_loop(
                    self._on_analysis_finished, cancel_event
                )

        threading.Thread(target=_run_analysis, daemon=True).start()

    def run_designer(self) -> bool:
        """Validate inputs, run 1D primer design analysis, and update UI.

        Returns:
            True on success, False on validation failure or analysis error.

        .. deprecated::
            Use :meth:`_start_designer` for new callers. This synchronous
            wrapper is retained for backwards compatibility with existing tests.
        """
        params = self.form.validate_and_get_params()
        if params is None:
            return False

        (
            clean_seq,
            min_length,
            mode,
            threshold,
            max_overlap,
            filter_dna_enabled,
            max_binding_sites,
        ) = params
        self.primer_list.controls.clear()

        try:
            template_dna: DNA | None = None
            if filter_dna_enabled:
                clean_tpl = clean_sequence(self.input_data.template)
                if not clean_tpl:
                    self.form.show_error(
                        "Template DNA sequence is required when check against "
                        "template is enabled. Please enter a template in the "
                        "Input view."
                    )
                    self.app_page.update()
                    return False
                t_type = (
                    DNAType.CIRCULAR
                    if self.input_data.template_circular
                    else DNAType.LINEAR
                )
                template_dna = DNA(clean_tpl, dna_type=t_type)

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
                template=template_dna,
                max_origin_count=max_binding_sites,
            )
            self._cached_designer = designer
            self._update_chart_and_primer_list(
                designer,
                self._compute_origin_counts(designer, template_dna),
                mode,
            )

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

    def _handle_run_pcr(self, primer_seq: str, primer_name: str) -> None:
        """Handle running PCR using template with a primer from card."""
        if self.on_run_pcr:
            self.on_run_pcr(primer_seq, primer_name)
        else:
            clean_tpl = clean_sequence(self.input_data.template)
            if not clean_tpl:
                show_error_dialog(
                    self.app_page,
                    "Template Required",
                    "Please enter a DNA template in the Input view before "
                    "running PCR.",
                )

    def _on_primer_selected(self, dimer: PrimerDimer, step_index: int) -> None:
        """Handle primer selection: add or raise self-dimer card."""
        card_id = f"1d_dimer_{dimer.primer_1.seq}_{step_index}"

        origin_count: int | None = None
        if self.form.filter_dna_checkbox.value:
            clean_tpl = clean_sequence(self.input_data.template)
            if clean_tpl:
                t_type = (
                    DNAType.CIRCULAR
                    if self.input_data.template_circular
                    else DNAType.LINEAR
                )
                tpl_dna = DNA(clean_tpl, dna_type=t_type)
                repliconf = Repliconf(tpl_dna, dimer.primer_1)
                repliconf.search()
                origin_count = len(repliconf.origin_db.fwd) + len(
                    repliconf.origin_db.rev
                )

        def _factory() -> DismissibleSelfDimerCard:
            font_family = self.settings.get("font_family", "Roboto Mono")
            return DismissibleSelfDimerCard(
                card_id=card_id,
                dimer=dimer,
                settings=self.settings,
                dismiss_callback=self._dismiss_card,
                font_family=font_family,
                step_index=step_index,
                on_run_pcr_callback=self._handle_run_pcr,
                origin_count=origin_count,
            )

        self._bring_card_to_top_or_add(card_id, _factory)

    def _clear_all(self, e: ft.ControlEvent | None = None) -> None:
        """Clear all inputs, parameters, error messages, and results."""
        self.form.dna_input.value = ""
        self.form.length_display.value = "0"
        self.form.min_len_input.value = ""
        self.form.max_quality_input.value = ""
        self.form.max_overlap_input.value = ""
        self.form.filter_dna_checkbox.value = False
        self.form.max_binding_sites_input.value = ""
        self.form.max_binding_sites_input.disabled = True
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
            "filter_dna": bool(self.form.filter_dna_checkbox.value),
            "max_binding_sites": (
                self.form.max_binding_sites_input.value or ""
            ),
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

        filter_dna_val = params.get("filter_dna")
        self.form.filter_dna_checkbox.value = bool(filter_dna_val)
        self.form.max_binding_sites_input.disabled = (
            not self.form.filter_dna_checkbox.value
        )
        max_sites_val = params.get("max_binding_sites")
        self.form.max_binding_sites_input.value = (
            str(max_sites_val) if max_sites_val is not None else ""
        )

        self.form.clear_errors()
        self.app_page.update()

        self._show_notification("Parameters loaded successfully.")
        self.run_designer()
