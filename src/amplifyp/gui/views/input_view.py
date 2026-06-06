# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Input view for DNA template and primers."""

import threading
from typing import Any

import flet as ft

from amplifyp.gui.settings import GUIColors, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.util import clean_sequence, format_sequence


class InputView(ft.Row):  # type: ignore[misc]
    """Input view for DNA template and primers."""

    def __init__(
        self,
        page: ft.Page,
        input_data: GUIInput | None = None,
        settings: GUISettings | None = None,
        on_change: Any | None = None,
        on_stop_editing: Any | None = None,
    ) -> None:
        """Initialize the InputView."""
        super().__init__(
            expand=True, vertical_alignment=ft.CrossAxisAlignment.STRETCH
        )
        self.app_page = page
        self.input_data = input_data if input_data is not None else GUIInput()
        self.settings = settings if settings is not None else GUISettings()
        self.on_change = on_change
        self.on_stop_editing_callback = on_stop_editing
        self._focus_timer: threading.Timer | None = None
        self.focused_primer_index: int | None = None

        font_family = self.settings.get("font_family", "Roboto Mono")

        self.template_sequence = ft.TextField(
            multiline=True,
            max_lines=None,
            border=ft.InputBorder.NONE,
            text_align=ft.TextAlign.LEFT,
            hint_text="Enter DNA sequence here...",
            content_padding=0,
            on_change=self.on_change_handler,
            on_focus=self.handle_field_focus,
            on_blur=self.handle_field_blur,
            on_submit=self.handle_field_submit,
            text_style=ft.TextStyle(font_family=font_family),
        )
        self.template_circular = ft.Checkbox(
            label="Circular Template",
            value=False,
            on_change=self.on_change_handler,
        )
        self.circular_container = ft.Container(
            content=self.template_circular,
            border=ft.Border.all(1, GUIColors.OUTLINE),
            border_radius=5,
            padding=ft.Padding(10, 0, 10, 0),
            height=32,
            alignment=ft.Alignment(0, 0),
        )

        self.name_column_width = 150.0
        self.primers_list = ft.ListView(
            expand=True, spacing=0, padding=0, scroll=ft.ScrollMode.ALWAYS
        )
        self.primers_header = ft.Row(
            [
                ft.Container(
                    content=ft.Text(
                        "Active",
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_small", 12),
                    ),
                    width=55,
                    alignment=ft.Alignment(0, 0),
                ),
                ft.Container(
                    width=4,
                    bgcolor=GUIColors.DIVIDER_GREY,
                    margin=0,
                    height=36,
                ),
                ft.Container(
                    content=ft.Text(
                        "Name",
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_small", 12),
                    ),
                    width=self.name_column_width,
                    padding=ft.Padding(5, 0, 0, 0),
                    alignment=ft.Alignment(-1, 0),
                ),
                ft.GestureDetector(
                    on_pan_update=self.on_primer_divider_pan,
                    content=ft.Container(
                        width=4,
                        bgcolor=GUIColors.DIVIDER_GREY,
                        margin=0,
                        height=36,
                    ),
                    mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
                ),
                ft.Container(
                    content=ft.Text(
                        "Sequence",
                        weight=ft.FontWeight.BOLD,
                        size=self.settings.get("font_size_small", 12),
                    ),
                    expand=True,
                    padding=ft.Padding(5, 0, 0, 0),
                    alignment=ft.Alignment(-1, 0),
                ),
            ],
            alignment=ft.MainAxisAlignment.START,
            height=36,
            spacing=0,
        )
        self.primers_header_container = ft.Container(
            content=self.primers_header,
            padding=0,
            border=ft.Border(bottom=ft.BorderSide(4, GUIColors.DIVIDER_GREY)),
            height=40,
        )

        self.clear_primers_button = ft.OutlinedButton(
            "Clear All",
            icon=ft.Icons.DELETE_OUTLINE,
            tooltip="Clear All Primers",
            on_click=self.clear_primers,
            height=32,
        )

        # Primer Info Panel UI Components
        self.info_header = ft.Container(
            content=ft.Text(
                "Primer: -",
                weight=ft.FontWeight.BOLD,
                size=self.settings.get("font_size_default", 14),
                color=GUIColors.DIAGRAM_BLACK,
            ),
            bgcolor=GUIColors.INFO_HEADER_BG,
            padding=ft.Padding(10, 5, 10, 5),
            alignment=ft.Alignment(-1, 0),
        )

        self.info_seq_text = ft.Text(
            "",
            font_family=font_family,
            size=self.settings.get("font_size_body", 13),
        )
        self.info_tm_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
        )
        self.info_pairs_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
        )
        self.info_redundancy_text = ft.Text(
            "",
            size=self.settings.get("font_size_body", 13),
        )
        self.info_dimer_text = ft.Text(
            "",
            color=GUIColors.ERROR_RED,
            size=self.settings.get("font_size_body", 13),
        )

        self.primer_info_panel = ft.Container(
            content=ft.Column(
                [
                    self.info_header,
                    ft.Container(
                        content=ft.Column(
                            [
                                self.info_seq_text,
                                self.info_tm_text,
                                self.info_pairs_text,
                                self.info_redundancy_text,
                                self.info_dimer_text,
                            ],
                            spacing=3,
                        ),
                        padding=ft.Padding(10, 5, 10, 10),
                    ),
                ],
                spacing=0,
            ),
            border=ft.Border.all(1, GUIColors.OUTLINE),
            border_radius=5,
            visible=False,
        )

        self.top_container = ft.Container(
            content=ft.Column(
                [
                    ft.Row(
                        [
                            ft.Text(
                                "Template Sequence", weight=ft.FontWeight.BOLD
                            ),
                            self.circular_container,
                        ],
                        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                        height=40,
                    ),
                    ft.Container(
                        content=ft.ListView(
                            [self.template_sequence],
                            expand=True,
                            scroll=ft.ScrollMode.ALWAYS,
                        ),
                        expand=True,
                        border=ft.Border.all(1, GUIColors.OUTLINE),
                        border_radius=5,
                        padding=0,
                    ),
                ],
                expand=True,
                spacing=5,
            ),
            expand=5,
        )

        self.right_fraction = 0.5

        self.bottom_container = ft.Container(
            content=ft.Column(
                [
                    ft.Row(
                        [
                            ft.Text("Primers", weight=ft.FontWeight.BOLD),
                            self.clear_primers_button,
                        ],
                        alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
                        height=40,
                    ),
                    ft.Container(
                        content=ft.Column(
                            [
                                self.primers_header_container,
                                self.primers_list,
                            ],
                            expand=True,
                            spacing=0,
                        ),
                        expand=True,
                        border=ft.Border.all(1, GUIColors.OUTLINE),
                        border_radius=5,
                        padding=0,
                    ),
                    self.primer_info_panel,
                ],
                expand=True,
                spacing=5,
            ),
            expand=5,
        )

        self.divider = ft.GestureDetector(
            on_pan_update=self.on_pan_update,
            content=ft.Container(
                width=5,
                bgcolor=GUIColors.DIVIDER_GREY,
                border_radius=5,
                margin=ft.Margin.symmetric(horizontal=5),
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
        )

        self.controls = [
            self.top_container,
            self.divider,
            self.bottom_container,
        ]

        # Sync initial UI state
        self.update_ui()

    def handle_field_focus(self, e: ft.ControlEvent) -> None:
        """Handle focus on input fields to cancel auto-trigger timer."""
        if self._focus_timer is not None:
            self._focus_timer.cancel()
            self._focus_timer = None

        if e.control.data is not None:
            self.focused_primer_index = e.control.data
            self.update_primer_info_panel()

    def handle_field_blur(self, e: ft.ControlEvent) -> None:
        """Handle blur on input fields to trigger results page after a delay."""
        if self._focus_timer is not None:
            self._focus_timer.cancel()
            self._focus_timer = None

        def timer_callback() -> None:
            if not self.page:
                return
            self.sync_to_state()
            if self.on_stop_editing_callback:
                self.on_stop_editing_callback()

        self._focus_timer = threading.Timer(0.15, timer_callback)
        self._focus_timer.daemon = True
        try:
            self._focus_timer.start()
        except RuntimeError:
            # Pyodide (WebAssembly) may not support starting OS threads when
            # SharedArrayBuffer / cross-origin isolation is unavailable.
            # Fall back to a direct synchronous invocation of the callback.
            self._focus_timer = None
            timer_callback()

    def handle_field_submit(self, e: ft.ControlEvent) -> None:
        """Handle submission (Enter key) to immediately trigger results."""
        if self._focus_timer is not None:
            self._focus_timer.cancel()
            self._focus_timer = None
        self.sync_to_state()
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback()

    def will_unmount(self) -> None:
        """Clean up when the view is unmounted."""
        if self._focus_timer is not None:
            self._focus_timer.cancel()
            self._focus_timer = None

    def _extract_primer_data_from_ui(self) -> list[dict[str, Any]]:
        """Extract primer data from UI controls."""
        ui_primers = []
        for container in self.primers_list.controls:
            if not isinstance(container, ft.Container):
                continue
            row = container.content
            if not isinstance(row, ft.Row) or len(row.controls) < 5:
                continue

            checkbox_control = row.controls[0]
            checkbox = (
                checkbox_control.content
                if isinstance(checkbox_control, ft.Container)
                else checkbox_control
            )
            name_tf = row.controls[2]
            seq_tf = row.controls[4]

            name_val = str(name_tf.value or "").strip()
            seq_val = clean_sequence(str(seq_tf.value or ""))
            is_active = bool(checkbox.value)

            ui_primers.append(
                {
                    "name": name_val,
                    "seq": seq_val,
                    "active": is_active,
                    "container": container,
                    "checkbox": checkbox,
                }
            )
        return ui_primers

    def _apply_activation_rules(
        self, ui_primers: list[dict[str, Any]]
    ) -> tuple[list[dict[str, Any]], bool]:
        """Apply auto-activation and auto-inactivation and deletion rules.

        Returns (filtered_ui_primers, rebuild_needed).
        """
        filtered_primers = []
        should_rebuild = False

        for i, p in enumerate(ui_primers):
            name_val = p["name"]
            seq_val = p["seq"]
            is_active = p["active"]
            checkbox = p["checkbox"]

            # Deletion rule: if both are empty, delete it
            if not name_val and not seq_val:
                should_rebuild = True
                continue

            prev_primer = (
                self.input_data.primers[i]
                if i < len(self.input_data.primers)
                else None
            )
            prev_name = (
                str(prev_primer.get("name", "")).strip() if prev_primer else ""
            )
            prev_seq = (
                clean_sequence(str(prev_primer.get("seq", "")))
                if prev_primer
                else ""
            )

            # Auto-activation rule: if it was previously incomplete
            # and both fields are now filled
            if (not prev_name or not prev_seq) and (name_val and seq_val):
                is_active = True
                checkbox.value = True

            # Auto-inactivation rule: if either is empty, set active to False
            if not name_val or not seq_val:
                if is_active:
                    is_active = False
                    checkbox.value = False
                    should_rebuild = True

            p["active"] = is_active
            filtered_primers.append(p)

        return filtered_primers, should_rebuild

    def _detect_and_mark_duplicates(
        self, ui_primers: list[dict[str, Any]]
    ) -> tuple[list[dict[str, Any]], bool]:
        """Detect duplicates and mark container colors.

        Returns (primers_list_for_state, rebuild_needed).
        """
        names_count: dict[str, int] = {}
        seqs_count: dict[str, int] = {}
        for p in ui_primers:
            n_lower = p["name"].lower()
            s_lower = p["seq"].lower()
            if n_lower:
                names_count[n_lower] = names_count.get(n_lower, 0) + 1
            if s_lower:
                seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1

        primers = []
        should_rebuild = False
        for p in ui_primers:
            container = p["container"]
            n_lower = p["name"].lower()
            s_lower = p["seq"].lower()

            is_dup = False
            if n_lower and names_count.get(n_lower, 0) > 1:
                is_dup = True
            if s_lower and seqs_count.get(s_lower, 0) > 1:
                is_dup = True

            new_color = GUIColors.DUPLICATE_BG if is_dup else None
            if container.bgcolor != new_color:
                container.bgcolor = new_color

            primers.append(
                {
                    "name": p["name"],
                    "seq": p["seq"],
                    "active": p["active"],
                }
            )
        return primers, should_rebuild

    def sync_to_state(self) -> None:
        """Sync current UI controls back to the central state."""
        self.template_sequence.value = self.template_sequence.value or ""
        self.input_data.template = clean_sequence(
            str(self.template_sequence.value)
        )
        self.input_data.template_circular = bool(self.template_circular.value)

        ui_primers = self._extract_primer_data_from_ui()
        filtered_ui_primers, should_rebuild = self._apply_activation_rules(
            ui_primers
        )
        primers, dup_rebuild = self._detect_and_mark_duplicates(
            filtered_ui_primers
        )
        should_rebuild = should_rebuild or dup_rebuild

        # If the number of non-empty primers changed, we rebuild
        if len(primers) != len(
            [
                p
                for p in self.input_data.primers
                if p.get("name", "").strip() or p.get("seq", "").strip()
            ]
        ):
            should_rebuild = True

        self.input_data.primers = primers

        if should_rebuild:
            self.update_ui()

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        font_family = self.settings.get("font_family", "Roboto Mono")
        self.template_sequence.text_style = ft.TextStyle(
            font_family=font_family
        )
        self.template_sequence.value = self.input_data.template
        self.template_circular.value = self.input_data.template_circular

        self.primers_list.controls.clear()

        # Filter out all empty primers from state first
        clean_primers = [
            p
            for p in self.input_data.primers
            if str(p.get("name", "")).strip()
            or clean_sequence(str(p.get("seq", ""))).strip()
        ]

        # Always append one trailing empty row
        clean_primers.append({"name": "", "seq": "", "active": False})
        self.input_data.primers = clean_primers

        # Check for duplicates to highlight them on load/update
        names_count: dict[str, int] = {}
        seqs_count: dict[str, int] = {}
        for p in self.input_data.primers:
            n_lower = str(p.get("name", "")).strip().lower()
            s_lower = clean_sequence(str(p.get("seq", ""))).lower()
            if n_lower:
                names_count[n_lower] = names_count.get(n_lower, 0) + 1
            if s_lower:
                seqs_count[s_lower] = seqs_count.get(s_lower, 0) + 1

        for idx, p in enumerate(self.input_data.primers):
            name_val = p["name"]
            seq_val = p["seq"]
            is_active = p.get("active", True)

            n_lower = name_val.strip().lower()
            s_lower = seq_val.lower()
            is_dup = False
            if n_lower and names_count.get(n_lower, 0) > 1:
                is_dup = True
            if s_lower and seqs_count.get(s_lower, 0) > 1:
                is_dup = True

            checkbox = ft.Checkbox(
                value=is_active,
                on_change=self.on_change_handler,
            )
            checkbox_container = ft.Container(
                content=checkbox,
                width=55,
                alignment=ft.Alignment(0, 0),
            )
            name_edit = ft.TextField(
                value=name_val,
                hint_text="New Primer Name",
                dense=True,
                content_padding=ft.Padding(5, 0, 0, 0),
                height=30,
                width=self.name_column_width,
                border=ft.InputBorder.NONE,
                data=idx,
                on_focus=self.handle_field_focus,
                on_blur=self.handle_field_blur,
                on_submit=self.handle_field_submit,
            )
            seq_edit = ft.TextField(
                value=seq_val,
                hint_text="New Primer Sequence",
                dense=True,
                expand=True,
                content_padding=ft.Padding(5, 0, 0, 0),
                height=30,
                border=ft.InputBorder.NONE,
                text_style=ft.TextStyle(font_family=font_family),
                data=idx,
                on_focus=self.handle_field_focus,
                on_blur=self.handle_field_blur,
                on_submit=self.handle_field_submit,
            )
            divider = ft.GestureDetector(
                on_pan_update=self.on_primer_divider_pan,
                content=ft.Container(
                    width=4,
                    bgcolor=GUIColors.DIVIDER_GREY,
                    margin=0,
                    height=30,
                ),
                mouse_cursor=ft.MouseCursor.RESIZE_LEFT_RIGHT,
            )

            active_divider = ft.Container(
                width=4,
                bgcolor=GUIColors.DIVIDER_GREY,
                margin=0,
                height=30,
            )

            row = ft.Row(
                [
                    checkbox_container,
                    active_divider,
                    name_edit,
                    divider,
                    seq_edit,
                ],
                expand=True,
                spacing=0,
            )

            row_container = ft.Container(
                content=row,
                bgcolor=GUIColors.DUPLICATE_BG if is_dup else None,
                padding=0,
                height=30,
            )
            self.primers_list.controls.append(row_container)

        self.update_primer_info_panel()

    def on_primer_divider_pan(self, e: ft.DragUpdateEvent) -> None:
        """Handle dragging the vertical divider between Name and Sequence."""
        delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
        self.name_column_width = max(
            80.0, min(300.0, self.name_column_width + delta_x)
        )
        # Update the width of the Name header control
        self.primers_header.controls[2].width = self.name_column_width
        # Update the width of all Name TextFields in the list controls
        for container in self.primers_list.controls:
            if isinstance(container, ft.Container):
                row = container.content
                if isinstance(row, ft.Row) and len(row.controls) >= 3:
                    name_tf = row.controls[2]
                    if isinstance(name_tf, ft.TextField):
                        name_tf.width = self.name_column_width
        self.app_page.update()

    def on_change_handler(self, e: ft.ControlEvent) -> None:
        """Handle change in input fields."""
        self.sync_to_state()
        self.update_primer_info_panel()
        if self.on_change:
            self.on_change(e)

    def clear_primers(self, e: ft.ControlEvent) -> None:
        """Clear all primers."""
        self.input_data.primers = []
        self.focused_primer_index = None
        self.update_ui()
        if self.on_change:
            self.on_change(e)
        if self.on_stop_editing_callback:
            self.on_stop_editing_callback()

    def on_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle resizing the bottom (right) container via the divider."""
        page_width = self.app_page.width
        if isinstance(page_width, (int, float)) and page_width > 0:
            delta_x = getattr(e.local_delta, "x", 0.0) if e.local_delta else 0.0
            # Calculate current pixel width of the right container
            current_width = page_width * self.right_fraction
            new_width = max(200.0, current_width - delta_x)
            # Ensure the left container stays at least 200px wide
            new_width = min(new_width, page_width - 200.0)

            # Recalculate fractions and set relative expand weights
            self.right_fraction = new_width / page_width
            self.bottom_container.expand = int(self.right_fraction * 1000)
            self.top_container.expand = int((1.0 - self.right_fraction) * 1000)
            self.app_page.update()

    def get_template(self) -> str:
        """Get the current template sequence."""
        return self.input_data.template

    def is_circular(self) -> bool:
        """Check if the template is circular."""
        return self.input_data.template_circular

    def get_primers(self) -> list[dict[str, Any]]:
        """Get the list of active primers."""
        return self.input_data.get_active_primers()

    def get_all_primers_state(self) -> list[dict[str, Any]]:
        """Get all primers (active and inactive) for serialization."""
        primers: list[dict[str, Any]] = []
        for p in self.input_data.primers:
            if (
                not str(p.get("name", "")).strip()
                and not clean_sequence(str(p.get("seq", ""))).strip()
            ):
                continue
            primers.append(
                {
                    "name": p["name"],
                    "seq": format_sequence(p["seq"]),
                    "active": p.get("active", True),
                }
            )
        return primers

    def get_state(self) -> dict[str, Any]:
        """Get the current input data state for serialization."""
        self.sync_to_state()
        return self.input_data.to_dict()

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current input data from deserialized data."""
        self.input_data.from_dict(state)
        self.update_ui()
        self.app_page.update()

    def update_primer_info_panel(self) -> None:
        """Update the primer information panel based on the focused primer."""
        if self.focused_primer_index is None:
            self.primer_info_panel.visible = False
            self.app_page.update()
            return

        try:
            primers = self.input_data.primers
            if (
                self.focused_primer_index < 0
                or self.focused_primer_index >= len(primers)
            ):
                self.primer_info_panel.visible = False
                self.app_page.update()
                return

            p_data = primers[self.focused_primer_index]
            name_val = p_data.get("name", "").strip()
            seq_val = clean_sequence(p_data.get("seq", ""))

            if not seq_val:
                self.primer_info_panel.visible = False
                self.app_page.update()
                return

            from amplifyp.dimer import PrimerDimerGenerator
            from amplifyp.dna import Primer
            from amplifyp.melting import (
                calculate_tm_lander_amplify4,
                calculate_tm_santalucia_1998_owczarzy_2008,
            )

            primer_obj = Primer(sequence=seq_val, name=name_val)

            # Header
            header_text = (
                f"Primer: {name_val}" if name_val else f"Primer: {seq_val}"
            )
            self.info_header.content.value = header_text

            # Sequence details
            self.info_seq_text.value = (
                f"{len(primer_obj)} bp:   {primer_obj.seq.upper()}"
            )

            # Melting temperature
            tm_method = self.settings.get(
                "tm_method", "SantaLucia 1998 / Owczarzy 2008"
            )
            if tm_method == "Lander / Amplify 4":
                from amplifyp.settings import TMSettings

                tm_settings = TMSettings(
                    dna_conc=self.settings._safe_float("tm_dna_conc", 50.0),
                    monovalent_salt_conc=self.settings._safe_float(
                        "tm_mono_salt", 50.0
                    ),
                )
                tm = calculate_tm_lander_amplify4(primer_obj, tm_settings)
            else:
                from amplifyp.settings import TMSettings

                tm_settings = TMSettings(
                    dna_conc=self.settings._safe_float("tm_dna_conc", 50.0),
                    dnap_conc=self.settings._safe_float("tm_dnap_conc", 0.0),
                    monovalent_salt_conc=self.settings._safe_float(
                        "tm_mono_salt", 50.0
                    ),
                    divalent_salt_conc=self.settings._safe_float(
                        "tm_div_salt", 0.0
                    ),
                    dnTP_conc=self.settings._safe_float("tm_dNTP_conc", 0.0),
                )
                tm = calculate_tm_santalucia_1998_owczarzy_2008(
                    primer_obj, tm_settings
                )

            self.info_tm_text.value = f"Tm = {tm:.2f}°C"

            # Base counts / percentage
            self.info_pairs_text.value = (
                f"{primer_obj.count_at()} AT Pairs, "
                f"{primer_obj.count_cg()} GC Pairs, "
                f"{primer_obj.ratio_at() * 100:.1f}% AT"
            )

            # Redundancy
            if primer_obj.redundant_base_count == 0:
                self.info_redundancy_text.value = "No redundant bases."
            else:
                self.info_redundancy_text.value = (
                    f"{primer_obj.redundant_base_count} redundant bases "
                    f"(redundancy fold = {primer_obj.redundancy_fold})."
                )

            # Dimer potential check
            pd_settings = self.settings.get_primer_dimer_settings()
            generator = PrimerDimerGenerator(settings=pd_settings)

            # Check self-dimer and cross-dimer potential
            # against all active primers
            active_primers = self.input_data.get_active_primers()
            max_dimer = None

            # Always check self-dimer
            res_self = generator.generate_primer_dimer(primer_obj, primer_obj)
            if (
                res_self.quality > pd_settings.threshold
                and res_self.overlap > pd_settings.min_overlap
            ):
                max_dimer = res_self

            for ap in active_primers:
                ap_name = ap.get("name", "").strip()
                ap_seq = clean_sequence(ap.get("seq", ""))
                if not ap_seq:
                    continue
                try:
                    ap_obj = Primer(sequence=ap_seq, name=ap_name)
                    res = generator.generate_primer_dimer(primer_obj, ap_obj)
                    if (
                        res.quality > pd_settings.threshold
                        and res.overlap > pd_settings.min_overlap
                    ):
                        if max_dimer is None or res.quality > max_dimer.quality:
                            max_dimer = res
                except ValueError:
                    continue

            if max_dimer is not None:
                self.info_dimer_text.value = (
                    "Potential Primer Dimer with quality = "
                    f"{max_dimer.quality:.1f} "
                    f"and overlap = {max_dimer.overlap}"
                )
                self.info_dimer_text.visible = True
            else:
                self.info_dimer_text.value = ""
                self.info_dimer_text.visible = False

            self.primer_info_panel.visible = True

        except Exception:
            self.primer_info_panel.visible = False

        self.app_page.update()
