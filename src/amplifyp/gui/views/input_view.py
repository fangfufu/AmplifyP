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

"""Input View for the Flet application."""

from typing import Any

import flet as ft

from amplifyp.gui.state import GUIState
from amplifyp.gui.util import clean_sequence, format_sequence


class InputView(ft.Column):  # type: ignore[misc]
    """Input view for DNA template and primers."""

    def __init__(
        self,
        page: ft.Page,
        state: GUIState | None = None,
        on_change: Any | None = None,
    ) -> None:
        """Initialize the InputView."""
        super().__init__(expand=True)
        self.app_page = page
        self.state = state if state is not None else GUIState()
        self.on_change = on_change

        self.template_sequence = ft.TextField(
            multiline=True,
            expand=True,
            label="Template Sequence",
            text_align=ft.TextAlign.LEFT,
            hint_text="Enter DNA sequence here...",
            on_change=self.on_change_handler,
            text_style=ft.TextStyle(font_family="Roboto Mono"),
        )
        self.template_circular = ft.Checkbox(
            label="Circular Template",
            value=False,
            on_change=self.on_change_handler,
        )

        self.primers_list = ft.ListView(expand=True, spacing=2)
        self.primer_name_input = ft.TextField(label="Primer Name", expand=1)
        self.primer_seq_input = ft.TextField(
            label="Primer Sequence",
            expand=3,
            text_style=ft.TextStyle(font_family="Roboto Mono"),
        )

        self.top_container = ft.Container(
            content=ft.Column(
                [self.template_sequence, self.template_circular], expand=True
            ),
            expand=True,
        )

        self.bottom_container = ft.Container(
            content=ft.Column(
                [
                    ft.Text("Primers", weight=ft.FontWeight.BOLD),
                    ft.Container(
                        content=self.primers_list,
                        expand=True,
                        border=ft.Border.all(1, ft.Colors.OUTLINE),
                        border_radius=5,
                        padding=5,
                    ),
                    ft.Row(
                        [
                            self.primer_name_input,
                            self.primer_seq_input,
                            ft.FilledButton(
                                "Add", on_click=self.add_primer_clicked
                            ),
                        ]
                    ),
                ],
                expand=True,
            ),
            height=300,
        )

        self.divider = ft.GestureDetector(
            on_pan_update=self.on_pan_update,
            content=ft.Container(
                height=5,
                bgcolor=ft.Colors.GREY_400,
                border_radius=5,
                margin=ft.Margin.symmetric(vertical=5),
            ),
            mouse_cursor=ft.MouseCursor.RESIZE_UP_DOWN,
        )

        self.controls = [
            self.top_container,
            self.divider,
            self.bottom_container,
        ]

        # Sync initial UI state
        self.update_ui()

    def sync_to_state(self) -> None:
        """Sync current UI controls back to the central state."""
        self.state.template = clean_sequence(
            str(self.template_sequence.value or "")
        )
        self.state.template_circular = bool(self.template_circular.value)

        primers = []
        for control in self.primers_list.controls:
            if isinstance(control, ft.ListTile) and isinstance(
                control.data, dict
            ):
                is_active = True
                if isinstance(control.leading, ft.Checkbox):
                    is_active = bool(control.leading.value)
                primers.append(
                    {
                        "name": str(control.data.get("name", "")),
                        "seq": clean_sequence(str(control.data.get("seq", ""))),
                        "active": is_active,
                    }
                )
        self.state.primers = primers

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central state."""
        self.template_sequence.value = self.state.template
        self.template_circular.value = self.state.template_circular

        self.primers_list.controls.clear()
        for p in self.state.primers:
            name_val = p["name"]
            seq_val = p["seq"]
            is_active = p.get("active", True)
            self.primers_list.controls.append(
                ft.ListTile(
                    dense=True,
                    content_padding=ft.Padding(0, 0, 0, 0),
                    leading=ft.Checkbox(
                        value=is_active, on_change=self.on_change_handler
                    ),
                    title=ft.Text(
                        f"{name_val}, {seq_val}",
                        font_family="Roboto Mono",
                    ),
                    data={"name": name_val, "seq": seq_val},
                    trailing=ft.IconButton(
                        ft.Icons.DELETE,
                        on_click=lambda e, n=name_val, s=seq_val: (
                            self.remove_primer(e, n, s)
                        ),
                    ),
                )
            )

    def on_change_handler(self, e: ft.ControlEvent) -> None:
        """Handle change in input fields."""
        self.sync_to_state()
        if self.on_change:
            self.on_change(e)

    def add_primer_clicked(self, e: ft.ControlEvent) -> None:
        """Handle adding a new primer."""
        self.sync_to_state()  # Sync first to preserve any un-synced edits
        if self.primer_name_input.value and self.primer_seq_input.value:
            name_val = str(self.primer_name_input.value)
            seq_val = clean_sequence(str(self.primer_seq_input.value))

            # Add to state and then update UI
            self.state.primers.append(
                {"name": name_val, "seq": seq_val, "active": True}
            )
            self.update_ui()

            self.primer_name_input.value = ""
            self.primer_seq_input.value = ""
            self.app_page.update()

            if self.on_change:
                self.on_change(e)

    def remove_primer(self, e: ft.ControlEvent, name: str, seq: str) -> None:
        """Handle removing a primer."""
        self.sync_to_state()  # Sync first to preserve any un-synced edits
        # Update state directly
        self.state.primers = [
            p
            for p in self.state.primers
            if not (p["name"] == name and p["seq"] == seq)
        ]
        self.update_ui()
        self.app_page.update()

        if self.on_change:
            self.on_change(e)

    def on_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle resizing the bottom container via the divider."""
        delta_y = getattr(e.local_delta, "y", 0.0) if e.local_delta else 0.0
        self.bottom_container.height = max(
            100.0, float(self.bottom_container.height or 300.0) - delta_y
        )
        self.app_page.update()

    def get_template(self) -> str:
        """Get the current template sequence."""
        return self.state.template

    def is_circular(self) -> bool:
        """Check if the template is circular."""
        return self.state.template_circular

    def get_primers(self) -> list[dict[str, Any]]:
        """Get the list of active primers."""
        return self.state.get_active_primers()

    def get_all_primers_state(self) -> list[dict[str, Any]]:
        """Get all primers (active and inactive) for serialization."""
        primers: list[dict[str, Any]] = []
        for p in self.state.primers:
            primers.append(
                {
                    "name": p["name"],
                    "seq": format_sequence(p["seq"]),
                    "active": p.get("active", True),
                }
            )
        return primers

    def get_state(self) -> dict[str, Any]:
        """Get the current state for serialization."""
        self.sync_to_state()
        return self.state.to_dict()

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current state from deserialized data."""
        self.state.from_dict(state)
        self.update_ui()
        self.app_page.update()
