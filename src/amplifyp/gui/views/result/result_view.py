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

from amplifyp.dna import DNA, DNAType, Primer
from amplifyp.gui.settings import MAX_AMPLICONS_RENDER, GUIColors, GUISettings
from amplifyp.gui.user_data import GUIInput
from amplifyp.pcr import PCR

from .amplicon_drawing import AmpliconDetailCard
from .diagram_panel import ResultDiagramPanel
from .primer_drawing import ReplicationContextCard


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
        self._cached_pcr: PCR | None = None
        self._cached_state_key: tuple[dict[str, Any], dict[str, Any]] | None = (
            None
        )

        self.result_list = ft.ListView(
            expand=True, spacing=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.diagram_panel = ResultDiagramPanel(
            page=self.app_page,
            settings=self.settings,
            on_primer_click=self._show_context_map,
            on_amplicon_click=self._show_amplicon_dialog,
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
            self.diagram_panel,
            self.cards_header,
            ft.Container(content=self.result_list, expand=True),
        ]
        self.app_page.on_resize = self._handle_resize

    @property
    def diagram_stack(self) -> ft.Stack:
        """Shortcut property to the diagram stack for test compatibility."""
        return self.diagram_panel.diagram_stack

    @property
    def diagram_container(self) -> ft.Container:
        """Shortcut property to the diagram container for test compatibility."""
        return self.diagram_panel.diagram_container

    @property
    def divider(self) -> ft.GestureDetector:
        """Shortcut property to the divider resizer for test compatibility."""
        return self.diagram_panel.divider

    def _handle_resize(self, e: ft.ControlEvent) -> None:
        """Handle window resizing by redrawing the PCR diagram."""
        if self.diagram_panel.diagram_container.visible:
            self.run_pcr(keep_cards=True)

    def run_pcr(self, keep_cards: bool = False) -> bool:
        """Execute the PCR simulation and update the UI."""
        saved_cards = self._reset_pcr_ui(keep_cards)
        success = True
        try:
            current_state_key = (
                self.input_data.to_dict(),
                self.settings.to_dict(),
            )
            if (
                self._cached_state_key == current_state_key
                and self._cached_pcr is not None
            ):
                pcr = self._cached_pcr
            else:
                pcr = self._execute_pcr_simulation()
                self._cached_pcr = pcr
                self._cached_state_key = current_state_key
            num_amplicons = len(pcr.amplicons)

            if num_amplicons == 0:
                self.result_list.controls.append(
                    ft.Text("No amplicons found.", selectable=True)
                )
            else:
                self.diagram_panel.render_diagram(pcr)
                if num_amplicons > MAX_AMPLICONS_RENDER:
                    self.result_list.controls.append(
                        ft.Container(
                            content=ft.Text(
                                f"Warning: {num_amplicons} amplicons "
                                "found. Only the top "
                                f"{MAX_AMPLICONS_RENDER} (sorted by "
                                "quality score) are displayed to "
                                "prevent UI freeze.",
                                color=GUIColors.ERROR_RED,
                                weight=ft.FontWeight.BOLD,
                            ),
                            padding=10,
                        )
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
        self.diagram_panel.reset_ui()
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

        def dismiss(card: ft.Card) -> None:
            if card in self.result_list.controls:
                self.result_list.controls.remove(card)
                self._update_cards_header_visibility()
                self.app_page.update()

        context_card = ReplicationContextCard(
            primer_name=primer_name,
            padded_idx=padded_idx,
            conf=conf,
            var=var,
            settings=self.settings,
            dismiss_callback=dismiss,
        )
        self.result_list.controls.insert(0, context_card)
        self._update_cards_header_visibility()
        self.app_page.update()

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

        def dismiss(card: ft.Card) -> None:
            if card in self.result_list.controls:
                self.result_list.controls.remove(card)
                self._update_cards_header_visibility()
                self.app_page.update()

        amplicon_card = AmpliconDetailCard(
            amp=amp,
            settings=self.settings,
            dismiss_callback=dismiss,
        )
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
