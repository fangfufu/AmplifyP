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
from collections.abc import Sequence

import flet as ft

from amplifyp.dimer import PrimerDimer, PrimerDimerGenerator
from amplifyp.dna import DNA, DNADirection
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.utils.data_helpers import clean_sequence
from amplifyp.gui.utils.gui_helpers import show_error_dialog
from amplifyp.gui.views.designer.dismissible_self_dimer_card import (
    DismissibleSelfDimerCard,
)
from amplifyp.primer_designer import PrimerDesigner1D

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

        # Input fields (Top-left control panel)
        self.dna_input = ft.TextField(
            label="DNA Sequence",
            hint_text="e.g. ATGCGTACGT...",
            expand=True,
            multiline=False,
            autofocus=True,
            on_submit=self._run_designer_event,
        )
        self.min_len_input = ft.TextField(
            label="Min Length (bp)",
            value="18",
            expand=True,
            on_submit=self._run_designer_event,
        )
        self.mode_dropdown = ft.Dropdown(
            label="Direction",
            options=[
                ft.dropdown.Option("FWD", "Forward"),
                ft.dropdown.Option("REV", "Reverse"),
            ],
            value="FWD",
            expand=True,
        )
        pd_settings = self.settings.get_primer_dimer_settings()

        self.max_quality_input = ft.TextField(
            label="Max Quality",
            hint_text="Unconstrained if empty",
            value=f"{pd_settings.threshold:.1f}",
            expand=True,
            on_submit=self._run_designer_event,
        )
        self.max_overlap_input = ft.TextField(
            label="Max Overlap (bp)",
            hint_text="Unconstrained if empty",
            value=str(pd_settings.min_overlap),
            expand=True,
            on_submit=self._run_designer_event,
        )
        self.analyse_button = ft.FilledButton(
            "Analyse",
            icon=ft.Icons.PLAY_ARROW,
            tooltip="Run 1D Primer Truncation Analysis",
            on_click=self._run_designer_event,
        )
        self.error_text = ft.Text(
            "", color=GUIColours.ERROR_RED, visible=False, size=12
        )

        # Container for top-left controls
        self.top_left_controls = ft.Column(
            [
                ft.Text(
                    "1D Truncation Parameters",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_subheader", 16),
                ),
                ft.Row(
                    [
                        self.dna_input,
                    ]
                ),
                ft.Row(
                    [
                        self.min_len_input,
                        self.mode_dropdown,
                    ],
                    alignment=ft.MainAxisAlignment.START,
                    vertical_alignment=ft.CrossAxisAlignment.CENTER,
                    spacing=14,
                ),
                ft.Row(
                    [
                        self.max_quality_input,
                        self.max_overlap_input,
                        ft.Container(
                            content=self.analyse_button,
                            margin=ft.Margin.only(left=24),
                        ),
                    ],
                    alignment=ft.MainAxisAlignment.START,
                    vertical_alignment=ft.CrossAxisAlignment.CENTER,
                    spacing=14,
                ),
                self.error_text,
            ],
            spacing=8,
        )
        self.top_left_container = ft.Container(
            content=self.top_left_controls,
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
            expand=1,
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
            content=self._build_quality_bar_chart([]),
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
            expand=3,
            padding=5,
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
            base_w = float(page_w) * 0.25
            self.left_container.expand = None
            self.right_container.expand = 3
            self.left_container.width = base_w
        current_w = float(self.left_container.width or 200.0)
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
            self.chart_content_container.content = (
                self._build_quality_bar_chart(self._cached_designer.all_dimers)
            )
        try:
            if self.app_page:
                self.update()
        except RuntimeError:
            pass

    def _build_quality_bar_chart(
        self, dimers: Sequence[PrimerDimer]
    ) -> ft.Control:
        """Build vertical quality bar chart for generated self-dimers."""
        if not dimers:
            return ft.Container(
                content=ft.Text(
                    "No analysis results yet. Enter sequence and "
                    "click Analyse.",
                    size=12,
                    italic=True,
                    color=GUIColours.TEXT_ON_SURFACE,
                ),
                alignment=ft.Alignment(0, 0),
                expand=True,
            )

        max_quality = max((d.quality for d in dimers), default=1.0)
        max_quality = max(max_quality, 1.0)

        container_h = float(self.top_right_chart_container.height or 140.0)
        max_bar_h = max(10.0, container_h - 98.0)
        bar_stack_height = max_bar_h + 20.0

        bar_controls: list[ft.Control] = []

        for step_idx, dimer in enumerate(dimers):
            primer_seq = dimer.primer_1.seq
            primer_len = len(primer_seq)
            q_val = dimer.quality

            bar_h = max(6.0, (q_val / max_quality) * max_bar_h)

            bar_with_score = ft.Column(
                [
                    ft.Text(
                        f"{q_val:.1f}",
                        size=11,
                        weight=ft.FontWeight.BOLD,
                    ),
                    ft.Container(
                        width=22,
                        height=bar_h,
                        bgcolor=GUIColours.PRIMARY,
                        border_radius=4,
                        tooltip=(
                            f"Step {step_idx + 1}: {primer_len} bp\n"
                            f"Quality: {q_val:.1f}"
                        ),
                    ),
                ],
                height=bar_stack_height,
                alignment=ft.MainAxisAlignment.END,
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=2,
            )

            bar_column = ft.Column(
                [
                    bar_with_score,
                    ft.Text(f"{primer_len}bp", size=11),
                ],
                alignment=ft.MainAxisAlignment.END,
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=2,
            )

            bar_item = ft.Container(
                content=bar_column,
                padding=ft.Padding(6, 4, 6, 26),
                ink=True,
                on_click=lambda _ev, d=dimer, idx=step_idx: (
                    self._on_primer_selected(d, idx)
                ),
            )
            bar_controls.append(bar_item)

        return ft.Row(
            bar_controls,
            scroll=ft.ScrollMode.ALWAYS,
            alignment=ft.MainAxisAlignment.START,
            vertical_alignment=ft.CrossAxisAlignment.END,
            spacing=4,
        )

    def _run_designer_event(self, e: ft.ControlEvent) -> None:
        """Event handler wrapper for running analysis."""
        self.run_designer()

    def run_designer(self) -> bool:
        """Validate inputs, run 1D primer design analysis, and update UI."""
        self.error_text.visible = False
        self.error_text.value = ""
        self.primer_list.controls.clear()

        dna_raw = self.dna_input.value or ""
        clean_seq = clean_sequence(dna_raw)
        if not clean_seq:
            self._show_error("Please enter a valid DNA sequence.")
            self.app_page.update()
            return False

        min_len_raw = (self.min_len_input.value or "").strip()
        if not min_len_raw.isdigit():
            self._show_error("Minimum length must be a positive integer.")
            self.app_page.update()
            return False

        min_length = int(min_len_raw)
        if min_length <= 0:
            self._show_error("Minimum length must be greater than 0.")
            self.app_page.update()
            return False

        if min_length > len(clean_seq):
            self._show_error(
                f"Minimum length ({min_length}) cannot exceed sequence length "
                f"({len(clean_seq)})."
            )
            self.app_page.update()
            return False

        mode_val = self.mode_dropdown.value or "FWD"
        mode = DNADirection.FWD if mode_val == "FWD" else DNADirection.REV

        max_q_raw = (self.max_quality_input.value or "").strip()
        threshold: float | None = None
        if max_q_raw:
            try:
                threshold = float(max_q_raw)
            except ValueError:
                self._show_error("Max Quality must be a valid number.")
                self.app_page.update()
                return False

        max_overlap_raw = (self.max_overlap_input.value or "").strip()
        max_overlap: int | None = None
        if max_overlap_raw:
            if not max_overlap_raw.isdigit():
                self._show_error("Max Overlap must be a non-negative integer.")
                self.app_page.update()
                return False
            max_overlap = int(max_overlap_raw)

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

            font_size_default = self.settings.get("font_size_default", 14)
            font_size_header = self.settings.get("font_size_header", 18)

            # Update top-right quality bar chart
            self.chart_content_container.content = (
                self._build_quality_bar_chart(designer.all_dimers)
            )

            is_reverse = mode == DNADirection.REV
            seq_align = ft.TextAlign.RIGHT if is_reverse else ft.TextAlign.LEFT

            for step_idx, dimer in enumerate(designer.all_dimers):
                seq = dimer.primer_1.seq
                length = len(seq)
                overlap_str = f"Overlap: {dimer.overlap} bp"

                item_card = ft.Card(
                    content=ft.Container(
                        content=ft.Row(
                            [
                                ft.Column(
                                    [
                                        ft.Text(
                                            f"{length} bp",
                                            weight=ft.FontWeight.BOLD,
                                            size=font_size_header,
                                            text_align=ft.TextAlign.LEFT,
                                        ),
                                        ft.Container(
                                            content=ft.Text(
                                                seq,
                                                size=font_size_default,
                                                font_family="Roboto Mono",
                                                text_align=seq_align,
                                            ),
                                            alignment=(
                                                ft.Alignment(1, 0)
                                                if is_reverse
                                                else ft.Alignment(-1, 0)
                                            ),
                                        ),
                                    ],
                                    spacing=2,
                                    expand=True,
                                    horizontal_alignment=ft.CrossAxisAlignment.START,
                                ),
                                ft.Column(
                                    [
                                        ft.Container(
                                            content=ft.Text(
                                                f"Quality: {dimer.quality:.1f}",
                                                weight=ft.FontWeight.BOLD,
                                                size=font_size_default,
                                                color=GUIColours.DIAGRAM_BLACK,
                                            ),
                                            bgcolor=GUIColours.SELECTED_ROW_BG,
                                            padding=ft.Padding(6, 3, 6, 3),
                                            border_radius=4,
                                        ),
                                        ft.Container(
                                            content=ft.Text(
                                                overlap_str,
                                                weight=ft.FontWeight.BOLD,
                                                size=font_size_default,
                                                color=GUIColours.DIAGRAM_BLACK,
                                            ),
                                            bgcolor=GUIColours.SELECTED_ROW_BG,
                                            padding=ft.Padding(6, 3, 6, 3),
                                            border_radius=4,
                                        ),
                                    ],
                                    spacing=4,
                                    horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                                    alignment=ft.MainAxisAlignment.CENTER,
                                ),
                            ],
                            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                            vertical_alignment=ft.CrossAxisAlignment.CENTER,
                        ),
                        padding=10,
                        ink=True,
                        on_click=lambda _ev, d=dimer, idx=step_idx: (
                            self._on_primer_selected(d, idx)
                        ),
                    )
                )
                self.primer_list.controls.append(item_card)

        except (ValueError, RuntimeError, OSError) as ex:
            logger.exception("1D Primer Design failed: %s", ex)
            self._show_error(f"Error: {ex}")
            show_error_dialog(
                self.app_page,
                "Error running Primer Designer",
                f"{ex}\n{traceback.format_exc()}",
            )
            self.app_page.update()
            return False

        self.app_page.update()
        return True

    def _show_error(self, message: str) -> None:
        """Display validation error message."""
        self.error_text.value = message
        self.error_text.visible = True

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

    def _clear_all_cards(self, e: ft.ControlEvent) -> None:
        """Clear all right-hand panel self-dimer cards."""
        self.right_cards_list.controls.clear()
        self._update_cards_header_visibility()
        self.app_page.update()
