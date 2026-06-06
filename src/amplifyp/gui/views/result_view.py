# Copyright (C) 2026 Fufu Fang
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

"""Result View for the Flet application."""

import traceback
from typing import Any

import flet as ft
import flet.canvas as cv

from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.pcr import PCR


class ResultView(ft.Column):  # type: ignore[misc]
    """Result view for rendering PCR custom execution targets."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
    ) -> None:
        """Initialize the ResultView."""
        super().__init__(expand=True)
        self.app_page = page

        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()

        self.result_list = ft.ListView(
            expand=True, spacing=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.diagram_canvas = cv.Canvas(
            expand=True,
            shapes=[],
        )
        self.diagram_stack = ft.Stack(
            controls=[self.diagram_canvas],
        )
        self.diagram_scrollable = ft.Column(
            controls=[self.diagram_stack],
            scroll=ft.ScrollMode.ALWAYS,
            horizontal_alignment=ft.CrossAxisAlignment.STRETCH,
        )
        self.diagram_container = ft.Container(
            content=self.diagram_scrollable,
            visible=False,
            border=ft.Border.all(1, GUIColors.OUTLINE),
            border_radius=5,
            padding=10,
            height=300,
        )
        self.divider = ft.GestureDetector(
            on_pan_update=self._on_pan_update,
            content=ft.Container(
                height=5,
                bgcolor=GUIColors.DIVIDER_GREY,
                border_radius=5,
                margin=ft.Margin.symmetric(vertical=5),
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_UP_DOWN,
            visible=False,
        )

        self.clear_button = ft.TextButton(
            "Clear",
            icon=ft.Icons.DELETE_SWEEP,
            tooltip="Clear All Cards",
            on_click=self._clear_all_cards,
            visible=False,
        )
        self.cards_header = ft.Row(
            [
                ft.Text(
                    "Details",
                    weight=ft.FontWeight.BOLD,
                    size=self.settings.get("font_size_header", 18),
                ),
                self.clear_button,
            ],
            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            visible=False,
        )

        self.controls = [
            self.diagram_container,
            self.divider,
            self.cards_header,
            ft.Container(content=self.result_list, expand=True),
        ]
        self.app_page.on_resize = self._handle_resize

    def _handle_resize(self, e: ft.ControlEvent) -> None:
        """Handle window resizing by redrawing the PCR diagram."""
        if self.diagram_container.visible:
            self.run_pcr(keep_cards=True)

    def _on_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle vertical resizing of the diagram container."""
        delta_y = getattr(e.local_delta, "y", 0.0) if e.local_delta else 0.0
        self.diagram_container.height = max(
            150.0, float(self.diagram_container.height or 300.0) + delta_y
        )
        self.app_page.update()

    def run_pcr(self, keep_cards: bool = False) -> bool:
        """Execute the PCR simulation and update the UI."""
        saved_cards = self._reset_pcr_ui(keep_cards)
        success = True
        try:
            pcr = self._execute_pcr_simulation()
            num_amplicons = len(pcr.amplicons)

            if num_amplicons == 0:
                self.result_list.controls.append(
                    ft.Text("No amplicons found.", selectable=True)
                )
            else:
                self.diagram_container.visible = True
                self.divider.visible = True

                target_length = len(pcr.template)
                v_target, h_margin, c_width, t_width = (
                    self._calculate_canvas_dimensions(
                        target_length, num_amplicons
                    )
                )

                if target_length > 0:
                    self._draw_template_baseline(
                        v_target, h_margin, c_width, t_width, target_length
                    )

                fwd_bindings, rev_bindings = self._collect_primer_bindings(pcr)

                self._draw_primers(
                    fwd_bindings,
                    rev_bindings,
                    target_length,
                    t_width,
                    h_margin,
                    v_target,
                )

                self._draw_amplicons(
                    pcr,
                    target_length,
                    t_width,
                    h_margin,
                    v_target,
                    c_width,
                )

        except Exception as ex:
            self.result_list.controls.append(
                ft.Text(
                    f"Error: {ex}\n{traceback.format_exc()}",
                    color=GUIColors.ERROR_RED,
                )
            )
            from amplifyp.gui.util import show_error_dialog

            show_error_dialog(self.app_page, "Error running PCR", str(ex))
            success = False

        if keep_cards:
            self.result_list.controls.extend(saved_cards)
            self._update_cards_header_visibility()

        self.app_page.update()
        return success

    def _reset_pcr_ui(self, keep_cards: bool) -> list[Any]:
        """Reset the result view UI controls and canvas shapes."""
        saved_cards = list(self.result_list.controls) if keep_cards else []
        self.result_list.controls.clear()
        self._update_cards_header_visibility()
        self.diagram_canvas.shapes.clear()
        self.diagram_stack.controls.clear()
        self.diagram_stack.controls.append(self.diagram_canvas)
        self.diagram_container.visible = False
        self.divider.visible = False
        return saved_cards

    def _execute_pcr_simulation(self) -> PCR:
        """Clean sequences, build DNA and PCR objects, and run simulation."""
        from amplifyp.gui.util import clean_sequence

        clean_template = clean_sequence(self.input_data.template)
        t_type = (
            DNAType.CIRCULAR
            if self.input_data.template_circular
            else DNAType.LINEAR
        )
        template_dna = DNA(clean_template, dna_type=t_type)

        rep_settings = self.settings.get_replication_settings()
        pcr = PCR(template_dna, settings=rep_settings)

        primers = self.input_data.get_active_primers()
        for p in primers:
            name = p["name"]
            seq = clean_sequence(p["seq"])
            pcr.add_primer(Primer(sequence=seq, name=name))
        pcr.predict_amplicons()
        return pcr

    def _calculate_canvas_dimensions(
        self, target_length: int, num_amplicons: int
    ) -> tuple[float, float, float, float]:
        """Calculate drawing coordinates and set canvas/stack heights."""
        v_target = 100.0  # Y position of target baseline
        h_margin = 20.0  # X padding
        c_width = (
            max(600.0, self.app_page.width - 80.0)
            if self.app_page.width
            else 800.0
        )
        self.diagram_canvas.width = c_width
        t_width = c_width - (2.0 * h_margin)

        # Calculate vertical size based on number of amplicons
        v_frag_start = v_target + 40
        v_frag_step = 35
        canvas_height = v_frag_start + num_amplicons * v_frag_step + 30.0
        self.diagram_canvas.height = canvas_height
        self.diagram_stack.height = canvas_height

        return v_target, h_margin, c_width, t_width

    def _draw_template_baseline(
        self,
        v_target: float,
        h_margin: float,
        c_width: float,
        t_width: float,
        target_length: int,
    ) -> None:
        """Draw boundary lines, baseline, ticks, and texts."""
        # Draw vertical boundary lines at start and end of template
        boundary_paint = ft.Paint(
            color=GUIColors.DIAGRAM_BLACK,
            style=ft.PaintingStyle.STROKE,
            stroke_width=1.0,
        )
        self.diagram_canvas.shapes.append(
            cv.Path(
                [
                    cv.Path.MoveTo(h_margin, v_target),
                    cv.Path.LineTo(h_margin, v_target - 65),
                ],
                paint=boundary_paint,
            )
        )
        self.diagram_canvas.shapes.append(
            cv.Path(
                [
                    cv.Path.MoveTo(c_width - h_margin, v_target),
                    cv.Path.LineTo(c_width - h_margin, v_target - 65),
                ],
                paint=boundary_paint,
            )
        )

        # Draw 1 and target_length text
        self.diagram_canvas.shapes.append(
            cv.Text(
                h_margin,
                v_target - 85,
                "1",
                style=ft.TextStyle(
                    size=self.settings.get("font_size_map_baseline", 16),
                    weight=ft.FontWeight.BOLD,
                    color=GUIColors.DIAGRAM_BLACK,
                ),
            )
        )
        self.diagram_canvas.shapes.append(
            cv.Text(
                c_width - h_margin,
                v_target - 85,
                str(target_length),
                style=ft.TextStyle(
                    size=self.settings.get("font_size_map_baseline", 16),
                    weight=ft.FontWeight.BOLD,
                    color=GUIColors.DIAGRAM_BLACK,
                ),
                alignment=ft.Alignment(1.0, -1.0),
            )
        )

        # Target Baseline
        self.diagram_canvas.shapes.append(
            cv.Path(
                [
                    cv.Path.MoveTo(h_margin, v_target),
                    cv.Path.LineTo(c_width - h_margin, v_target),
                ],
                paint=ft.Paint(
                    color=GUIColors.DIAGRAM_BLACK,
                    style=ft.PaintingStyle.STROKE,
                    stroke_width=2.5,
                ),
            )
        )

        # Ticks dynamically scaled
        tick_interval = 100
        if target_length > 10000:
            tick_interval = 1000
        elif target_length > 5000:
            tick_interval = 500

        tick_paint = ft.Paint(
            color=GUIColors.DIAGRAM_BLACK,
            style=ft.PaintingStyle.STROKE,
            stroke_width=1.0,
        )
        for bp in range(0, target_length + 1, tick_interval):
            x_pos = h_margin + (bp / target_length * t_width)
            self.diagram_canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_pos, v_target - 3),
                        cv.Path.LineTo(x_pos, v_target + 3),
                    ],
                    paint=tick_paint,
                )
            )

    def _collect_primer_bindings(
        self, pcr: PCR
    ) -> tuple[
        dict[int, tuple[str, float, Any, Any]],
        dict[int, tuple[str, float, Any, Any]],
    ]:
        """Collect and group unique forward and reverse primer binding sites."""
        fwd_bindings = {}
        rev_bindings = {}
        for amp in pcr.amplicons:
            fwd_conf = next(
                (
                    c
                    for c in pcr.amplicon_generator.repliconfs
                    if c.primer is amp.fwd_origin
                ),
                None,
            )
            rev_conf = next(
                (
                    c
                    for c in pcr.amplicon_generator.repliconfs
                    if c.primer is amp.rev_origin
                ),
                None,
            )
            if fwd_conf is None or rev_conf is None:
                continue
            fwd_origin_point = fwd_conf.origin(amp.start)
            rev_origin_point = rev_conf.origin(amp.end)
            if fwd_origin_point is None or rev_origin_point is None:
                continue
            fwd_quality = fwd_origin_point.quality
            rev_quality = rev_origin_point.quality

            # Scale triangle size S based on quality score
            fwd_s = 6.0 + (max(0.1, min(1.0, fwd_quality)) * 10.0)
            rev_s = 6.0 + (max(0.1, min(1.0, rev_quality)) * 10.0)

            fwd_bindings[amp.start.index] = (
                amp.fwd_origin.name,
                fwd_s,
                fwd_conf,
                amp.start,
            )
            rev_bindings[amp.end.index] = (
                amp.rev_origin.name,
                rev_s,
                rev_conf,
                amp.end,
            )
        return fwd_bindings, rev_bindings

    def _draw_primers(
        self,
        fwd_bindings: dict[int, tuple[str, float, Any, Any]],
        rev_bindings: dict[int, tuple[str, float, Any, Any]],
        target_length: int,
        t_width: float,
        h_margin: float,
        v_target: float,
    ) -> None:
        """Draw forward and reverse primers, labels, and detectors."""
        # Draw FWD primers (blue, float above baseline, pointing down)
        for start_idx, (name, S, fwd_conf, fwd_var) in fwd_bindings.items():
            x_pos = (
                h_margin + (start_idx / target_length * t_width)
                if target_length
                else h_margin
            )
            # Connector line from baseline to floating tip
            self.diagram_canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_pos, v_target),
                        cv.Path.LineTo(x_pos, v_target - 25),
                    ],
                    paint=ft.Paint(
                        color=GUIColors.FWD_PRIMER,
                        style=ft.PaintingStyle.STROKE,
                        stroke_width=1.0,
                    ),
                )
            )
            # Down-pointing triangle of size S:
            # Tip: (x_pos, v_target - 25)
            # Top-Left: (x_pos - S/2, v_target - 25 - S)
            # Top-Right: (x_pos + S/2, v_target - 25 - S)
            self.diagram_canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_pos, v_target - 25),
                        cv.Path.LineTo(x_pos - S / 2, v_target - 25 - S),
                        cv.Path.LineTo(x_pos + S / 2, v_target - 25 - S),
                        cv.Path.Close(),
                    ],
                    paint=ft.Paint(
                        color=GUIColors.FWD_PRIMER,
                        style=ft.PaintingStyle.FILL,
                    ),
                )
            )
            # Add rotated label in the stack
            self.diagram_stack.controls.append(
                ft.Text(
                    name,
                    color=GUIColors.FWD_PRIMER,
                    size=self.settings.get("font_size_map_primer", 13),
                    weight=ft.FontWeight.BOLD,
                    left=x_pos - 15,
                    top=v_target - 25 - S - 38,
                    rotate=ft.Rotate(-1.5708),
                )
            )

            def fwd_tap(
                _: Any,
                n: str = name,
                i: int = start_idx,
                c: Any = fwd_conf,
                v: Any = fwd_var,
            ) -> None:
                self._show_context_map(n, i, c, v)

            # Click overlay
            self.diagram_stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=fwd_tap,
                    content=ft.Container(
                        bgcolor=GUIColors.TRANSPARENT,
                        width=20,
                        height=25 + S,
                    ),
                    left=x_pos - 10,
                    top=v_target - 25 - S,
                )
            )

        # Draw REV primers (red, float below baseline, pointing up)
        for end_idx, (name, S, rev_conf, rev_var) in rev_bindings.items():
            x_pos = (
                h_margin + (end_idx / target_length * t_width)
                if target_length
                else h_margin
            )
            # Connector line from baseline to floating tip
            self.diagram_canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_pos, v_target),
                        cv.Path.LineTo(x_pos, v_target + 25),
                    ],
                    paint=ft.Paint(
                        color=GUIColors.REV_PRIMER,
                        style=ft.PaintingStyle.STROKE,
                        stroke_width=1.0,
                    ),
                )
            )
            # Up-pointing triangle of size S:
            # Tip: (x_pos, v_target + 25)
            # Bottom-Left: (x_pos - S/2, v_target + 25 + S)
            # Bottom-Right: (x_pos + S/2, v_target + 25 + S)
            self.diagram_canvas.shapes.append(
                cv.Path(
                    [
                        cv.Path.MoveTo(x_pos, v_target + 25),
                        cv.Path.LineTo(x_pos - S / 2, v_target + 25 + S),
                        cv.Path.LineTo(x_pos + S / 2, v_target + 25 + S),
                        cv.Path.Close(),
                    ],
                    paint=ft.Paint(
                        color=GUIColors.REV_PRIMER,
                        style=ft.PaintingStyle.FILL,
                    ),
                )
            )
            # Add rotated label in the stack
            self.diagram_stack.controls.append(
                ft.Text(
                    name,
                    color=GUIColors.REV_LABEL,
                    size=self.settings.get("font_size_map_primer", 13),
                    weight=ft.FontWeight.BOLD,
                    left=x_pos - 15,
                    top=v_target + 25 + S + 10,
                    rotate=ft.Rotate(-1.5708),
                )
            )

            def rev_tap(
                _: Any,
                n: str = name,
                i: int = end_idx,
                c: Any = rev_conf,
                v: Any = rev_var,
            ) -> None:
                self._show_context_map(n, i, c, v)

            # Click overlay
            self.diagram_stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=rev_tap,
                    content=ft.Container(
                        bgcolor=GUIColors.TRANSPARENT,
                        width=20,
                        height=25 + S,
                    ),
                    left=x_pos - 10,
                    top=v_target,
                )
            )

    def _draw_amplicons(
        self,
        pcr: PCR,
        target_length: int,
        t_width: float,
        h_margin: float,
        v_target: float,
        c_width: float,
    ) -> None:
        """Draw amplicons, text labels, and detectors."""
        v_frag_start = v_target + 40
        v_frag_step = 35

        for idx, amp in enumerate(pcr.amplicons):
            start_idx = amp.start.index
            end_idx = amp.end.index

            # X coords
            x_start = (
                h_margin + (start_idx / target_length * t_width)
                if target_length
                else h_margin
            )
            x_end = (
                h_margin + (end_idx / target_length * t_width)
                if target_length
                else (c_width - h_margin)
            )

            y_pos = v_frag_start + (idx * v_frag_step)

            # Dynamic thickness based on the quality score
            if amp.q_score < 300:
                bar_height = 8.0
            elif amp.q_score < 700:
                bar_height = 5.5
            elif amp.q_score < 1500:
                bar_height = 3.5
            elif amp.q_score < 4000:
                bar_height = 2.0
            else:
                bar_height = 1.0

            # Amplicon Bar (Black, filling path)
            amp_paint = ft.Paint(
                color=GUIColors.DIAGRAM_BLACK,
                style=ft.PaintingStyle.FILL,
            )
            if amp.circular:
                self.diagram_canvas.shapes.append(
                    cv.Path(
                        [
                            cv.Path.MoveTo(h_margin, y_pos),
                            cv.Path.LineTo(c_width - h_margin, y_pos),
                            cv.Path.LineTo(
                                c_width - h_margin,
                                y_pos + bar_height,
                            ),
                            cv.Path.LineTo(h_margin, y_pos + bar_height),
                            cv.Path.Close(),
                        ],
                        paint=amp_paint,
                    )
                )
                label_x = h_margin + (t_width / 2.0)
            else:
                self.diagram_canvas.shapes.append(
                    cv.Path(
                        [
                            cv.Path.MoveTo(x_start, y_pos),
                            cv.Path.LineTo(x_end, y_pos),
                            cv.Path.LineTo(x_end, y_pos + bar_height),
                            cv.Path.LineTo(x_start, y_pos + bar_height),
                            cv.Path.Close(),
                        ],
                        paint=amp_paint,
                    )
                )
                label_x = x_start + ((x_end - x_start) / 2.0)

            # Amplicon Length Text (just number, black, centered)
            self.diagram_canvas.shapes.append(
                cv.Text(
                    label_x,
                    y_pos + bar_height + 5,
                    str(len(amp.product)),
                    style=ft.TextStyle(
                        size=self.settings.get("font_size_map_amplicon", 13),
                        color=GUIColors.DIAGRAM_BLACK,
                    ),
                    alignment=ft.Alignment(0.0, -1.0),
                )
            )

            # Click overlay for the amplicon
            amp_width = max(
                10.0,
                (x_end - x_start)
                if not amp.circular
                else (c_width - 2 * h_margin),
            )
            amp_left = x_start if not amp.circular else h_margin
            self.diagram_stack.controls.append(
                ft.GestureDetector(
                    mouse_cursor=ft.MouseCursor.CLICK,
                    on_tap=lambda _, a=amp: self._show_amplicon_dialog(a),
                    content=ft.Container(
                        bgcolor=GUIColors.TRANSPARENT,
                        width=amp_width,
                        height=20 + bar_height,
                    ),
                    left=amp_left,
                    top=y_pos - 3,
                )
            )

    def _show_context_map(
        self,
        primer_name: str,
        padded_idx: int,
        conf: Any,
        var: Any,
    ) -> None:
        """Create and show context map card below the overview map."""
        card_id = f"context_{primer_name}_{padded_idx}"
        for ctrl in self.result_list.controls:
            if getattr(ctrl, "_card_id", None) == card_id:
                self.result_list.controls.remove(ctrl)
                self.result_list.controls.insert(0, ctrl)
                self._update_cards_header_visibility()
                self.app_page.update()
                return

        context_card = self._create_context_card(
            primer_name, padded_idx, conf, var
        )
        self.result_list.controls.insert(0, context_card)
        self._update_cards_header_visibility()
        self.app_page.update()

    def _get_template_substring(
        self, template: DNA, start: int, length: int
    ) -> str:
        """Get substring of template DNA, supporting circular wrapping."""
        N_len = len(template)
        if template.type == DNAType.CIRCULAR:
            return "".join(
                template.seq[i % N_len] for i in range(start, start + length)
            )
        else:
            return "".join(
                template.seq[i] if 0 <= i < N_len else "-"
                for i in range(start, start + length)
            )

    def _format_context_lines(
        self,
        primer_name: str,
        padded_idx: int,
        conf: Any,
        origin: Any,
        L: int,
        N: int,
        direction: DNADirection,
    ) -> tuple[str, str, str]:
        """Format context alignment lines for binding site context map."""
        if direction == DNADirection.FWD:
            start_genomic = (padded_idx - L) % N
            primer_display_seq = conf.primer.seq
            primer_label = f"{primer_name} (Forward)"
            primer_ends = ("5'", "3'")
            strength_display = origin.binding_strength_str[::-1]
        else:
            start_genomic = padded_idx % N
            primer_display_seq = conf.primer.seq[::-1]
            primer_label = f"{primer_name} (Reverse)"
            primer_ends = ("3'", "5'")
            strength_display = origin.binding_strength_str

        end_genomic = start_genomic + L
        start_num_str = str(((start_genomic) % N) + 1)
        end_num_str = str(((end_genomic - 1) % N) + 1)

        # Construct primer line:
        primer_line = (
            f"{primer_label:<29}"
            f"{primer_ends[0]}-{primer_display_seq}-{primer_ends[1]}"
        )

        # Construct strength line:
        bonds_line = f"{' ' * 12}{' ' * 20}{strength_display}"

        # Construct arrows:
        arrow_line = f"{' ' * 12}{' ' * 20}↓{' ' * (L - 2)}↓"

        # Construct numbers:
        num_line_chars = [" "] * (12 + 20 + L + 20)
        col1 = 12 + 20 - len(start_num_str) // 2
        for idx, char in enumerate(start_num_str):
            num_line_chars[col1 + idx] = char
        col2 = (12 + 20 + L - 1) - len(end_num_str) // 2
        for idx, char in enumerate(end_num_str):
            num_line_chars[col2 + idx] = char
        top_line = "".join(num_line_chars).rstrip()

        # Combined top line:
        top_line = f"{top_line}\n{arrow_line}"

        upstream_seq = self._get_template_substring(
            conf.template, start_genomic - 20, 20
        )
        binding_seq = self._get_template_substring(
            conf.template, start_genomic, L
        )
        downstream_seq = self._get_template_substring(
            conf.template, start_genomic + L, 20
        )

        context_line_prefix = "Context  "
        bottom_line = (
            f"{bonds_line}\n"
            f"{context_line_prefix}5'-{upstream_seq}"
            f"{binding_seq}{downstream_seq}-3'"
        )

        return top_line, primer_line, bottom_line

    def _create_context_card(
        self,
        primer_name: str,
        padded_idx: int,
        conf: Any,
        var: Any,
    ) -> ft.Card:
        """Create an ft.Card showing visually aligned replication context map.

        Includes a close button.
        """
        card_ref = ft.Ref[ft.Card]()

        def remove_card(e: Any) -> None:
            if card_ref.current in self.result_list.controls:
                self.result_list.controls.remove(card_ref.current)
                self._update_cards_header_visibility()
                self.app_page.update()

        origin = conf.origin(var)
        if origin is None:
            return ft.Card(
                content=ft.Text("Error: Replication origin not found")
            )
        primer_type = (
            "Forward Primer"
            if var.direction == DNADirection.FWD
            else "Reverse Primer"
        )
        L = len(conf.primer)
        N = len(conf.template)

        from amplifyp.gui.util import create_overlapped_sequence_view

        top_line, primer_line, bottom_line = self._format_context_lines(
            primer_name=primer_name,
            padded_idx=padded_idx,
            conf=conf,
            origin=origin,
            L=L,
            N=N,
            direction=var.direction,
        )

        font_family = self.settings.get("font_family", "Roboto Mono")
        font_size = self.settings.get("font_size_default", 14)
        diagram_text = create_overlapped_sequence_view(
            top_line=top_line,
            mid_line=primer_line,
            bottom_line=bottom_line,
            font_family=font_family,
            font_size=font_size,
        )

        card = ft.Card(
            ref=card_ref,
            content=ft.Container(
                padding=10,
                content=ft.Column(
                    [
                        ft.Row(
                            [
                                ft.Text(
                                    f"Context Map - {primer_name} "
                                    f"({primer_type}) Binding Site",
                                    weight=ft.FontWeight.BOLD,
                                    size=self.settings.get(
                                        "font_size_subheader", 16
                                    ),
                                ),
                                ft.IconButton(
                                    icon=ft.Icons.CLOSE,
                                    icon_size=18,
                                    tooltip="Dismiss",
                                    on_click=remove_card,
                                ),
                            ],
                            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                        ),
                        ft.Container(
                            content=ft.Row(
                                [diagram_text],
                                scroll=ft.ScrollMode.ALWAYS,
                            ),
                            padding=12,
                            border_radius=6,
                            border=ft.Border.all(1, GUIColors.OUTLINE_VARIANT),
                        ),
                    ],
                    tight=True,
                ),
            ),
        )
        card._card_id = f"context_{primer_name}_{padded_idx}"
        return card

    def _create_amplicon_card(self, amp: Any) -> ft.Card:
        """Create a Flet Card for displaying details of a single amplicon.

        Includes a close button.
        """
        card_ref = ft.Ref[ft.Card]()

        def remove_card(e: Any) -> None:
            if card_ref.current in self.result_list.controls:
                self.result_list.controls.remove(card_ref.current)
                self._update_cards_header_visibility()
                self.app_page.update()

        full_seq = str(amp.product.seq)
        fwd_len = len(amp.fwd_origin.seq)
        rev_len = len(amp.rev_origin.seq)
        mid_len = len(amp.product) - (fwd_len + rev_len)

        if len(full_seq) >= fwd_len + rev_len:
            fwd_part = full_seq[:fwd_len]
            mid_part = full_seq[fwd_len:-rev_len]
            rev_part = full_seq[-rev_len:]
        else:
            fwd_part = full_seq
            mid_part = ""
            rev_part = ""

        font_family = self.settings.get("font_family", "Roboto Mono")
        sequence_text = ft.Text(
            spans=[
                ft.TextSpan(
                    fwd_part,
                    style=ft.TextStyle(
                        color=GUIColors.FWD_PRIMER,
                        weight=ft.FontWeight.BOLD,
                    ),
                ),
                ft.TextSpan(
                    mid_part,
                    style=ft.TextStyle(
                        color=GUIColors.TEXT_ON_SURFACE,
                    ),
                ),
                ft.TextSpan(
                    rev_part,
                    style=ft.TextStyle(
                        color=GUIColors.REV_LABEL,
                        weight=ft.FontWeight.BOLD,
                    ),
                ),
            ],
            font_family=font_family,
            size=self.settings.get("font_size_body", 13),
            selectable=True,
        )

        card = ft.Card(
            ref=card_ref,
            content=ft.Container(
                padding=10,
                content=ft.Column(
                    [
                        ft.Row(
                            [
                                ft.Text(
                                    f"Amplicon: {len(amp.product)} bp",
                                    weight=ft.FontWeight.BOLD,
                                    size=self.settings.get(
                                        "font_size_subheader", 16
                                    ),
                                    selectable=True,
                                ),
                                ft.IconButton(
                                    icon=ft.Icons.CLOSE,
                                    icon_size=18,
                                    tooltip="Dismiss",
                                    on_click=remove_card,
                                ),
                            ],
                            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                        ),
                        ft.Text(
                            spans=[
                                ft.TextSpan("▶(primer: "),
                                ft.TextSpan(
                                    amp.fwd_origin.name,
                                    style=ft.TextStyle(
                                        color=GUIColors.FWD_PRIMER,
                                        weight=ft.FontWeight.BOLD,
                                    ),
                                ),
                                ft.TextSpan(f") — {mid_len} bp — (primer: "),
                                ft.TextSpan(
                                    amp.rev_origin.name,
                                    style=ft.TextStyle(
                                        color=GUIColors.REV_LABEL,
                                        weight=ft.FontWeight.BOLD,
                                    ),
                                ),
                                ft.TextSpan(
                                    f")◀      Q = {amp.q_score:.1f} "
                                    f"({amp.q_score_report_str(verbose=True)})"
                                ),
                            ],
                            selectable=True,
                            size=self.settings.get("font_size_body", 13),
                        ),
                        ft.Text(
                            "Amplified Sequence:",
                            weight=ft.FontWeight.BOLD,
                            selectable=True,
                        ),
                        ft.Container(
                            content=sequence_text,
                            padding=12,
                            border_radius=6,
                            border=ft.Border.all(1, GUIColors.OUTLINE_VARIANT),
                        ),
                    ]
                ),
            ),
        )
        card_id = (
            f"amplicon_{amp.fwd_origin.name}_{amp.rev_origin.name}_"
            f"{amp.start.index}_{amp.end.index}"
        )
        card._card_id = card_id
        return card

    def _show_amplicon_dialog(self, amp: Any) -> None:
        """Show details card of the selected amplicon below the overview map."""
        card_id = (
            f"amplicon_{amp.fwd_origin.name}_{amp.rev_origin.name}_"
            f"{amp.start.index}_{amp.end.index}"
        )
        for ctrl in self.result_list.controls:
            if getattr(ctrl, "_card_id", None) == card_id:
                self.result_list.controls.remove(ctrl)
                self.result_list.controls.insert(0, ctrl)
                self._update_cards_header_visibility()
                self.app_page.update()
                return

        amplicon_card = self._create_amplicon_card(amp)
        self.result_list.controls.insert(0, amplicon_card)
        self._update_cards_header_visibility()
        self.app_page.update()

    def _update_cards_header_visibility(self) -> None:
        """Toggle header visibility based on list content."""
        has_cards = len(self.result_list.controls) > 0
        self.cards_header.visible = has_cards
        self.clear_button.visible = has_cards

    def _clear_all_cards(self, e: Any) -> None:
        """Clear all detail cards below the overview map."""
        self.result_list.controls.clear()
        self._update_cards_header_visibility()
        self.app_page.update()
