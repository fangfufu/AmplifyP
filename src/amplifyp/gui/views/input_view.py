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


class InputView(ft.Column):  # type: ignore[misc]
    """Input view for DNA template and primers."""

    def __init__(self, page: ft.Page) -> None:
        """Initialize the InputView."""
        super().__init__(expand=True)
        self.app_page = page

        self.template_sequence = ft.TextField(
            multiline=True,
            expand=True,
            label="Template Sequence",
            text_align=ft.TextAlign.LEFT,
            hint_text="Enter DNA sequence here...",
        )
        self.template_circular = ft.Checkbox(
            label="Circular Template", value=False
        )

        self.primers_list = ft.ListView(expand=True, spacing=2)
        self.primer_name_input = ft.TextField(label="Primer Name", expand=1)
        self.primer_seq_input = ft.TextField(label="Primer Sequence", expand=3)

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

    def add_primer_clicked(self, e: ft.ControlEvent) -> None:
        """Handle adding a new primer."""
        if self.primer_name_input.value and self.primer_seq_input.value:
            name_val = str(self.primer_name_input.value)
            seq_val = str(self.primer_seq_input.value)
            self.primers_list.controls.append(
                ft.ListTile(
                    dense=True,
                    content_padding=ft.padding.all(0),
                    leading=ft.Checkbox(value=True),
                    title=ft.Text(f"{name_val}, {seq_val}"),
                    data={"name": name_val, "seq": seq_val},
                    trailing=ft.IconButton(
                        ft.Icons.DELETE,
                        on_click=lambda e, n=name_val, s=seq_val: (
                            self.remove_primer(e, n, s)
                        ),
                    ),
                )
            )
            self.primer_name_input.value = ""
            self.primer_seq_input.value = ""
            self.app_page.update()

    def remove_primer(self, e: ft.ControlEvent, name: str, seq: str) -> None:
        """Handle removing a primer."""
        for control in self.primers_list.controls[:]:
            if (
                isinstance(control, ft.ListTile)
                and isinstance(control.data, dict)
                and control.data.get("name") == name
                and control.data.get("seq") == seq
            ):
                self.primers_list.controls.remove(control)
        self.app_page.update()

    def on_pan_update(self, e: ft.DragUpdateEvent) -> None:
        """Handle resizing the bottom container via the divider."""
        delta_y = getattr(e.local_delta, "y", 0.0) if e.local_delta else 0.0
        self.bottom_container.height = max(
            100.0, float(self.bottom_container.height or 300.0) - delta_y
        )
        self.app_page.update()

    def get_template(self) -> str:
        """Get the current template sequence."""
        return str(self.template_sequence.value or "")

    def is_circular(self) -> bool:
        """Check if the template is circular."""
        return bool(self.template_circular.value)

    def get_primers(self) -> list[dict[str, Any]]:
        """Get the list of active primers."""
        primers: list[dict[str, Any]] = []
        for control in self.primers_list.controls:
            if isinstance(control, ft.ListTile):
                # If checkbox exists and is unchecked, skip it
                if (
                    isinstance(control.leading, ft.Checkbox)
                    and not control.leading.value
                ):
                    continue
                if not isinstance(control.data, dict):
                    continue
                primers.append(
                    {
                        "name": str(control.data.get("name", "")),
                        "seq": str(control.data.get("seq", "")),
                    }
                )
        return primers

    def get_all_primers_state(self) -> list[dict[str, Any]]:
        """Get all primers (active and inactive) for serialization."""

        def format_seq(seq_val: str | None) -> str:
            if not seq_val:
                return ""
            # Remove literal escaped whitespace characters
            clean_seq = (
                str(seq_val)
                .replace("\\n", "")
                .replace("\\t", "")
                .replace("\\r", "")
            )
            clean_seq = "".join(clean_seq.split()).upper()
            return "\n".join(
                [clean_seq[i : i + 80] for i in range(0, len(clean_seq), 80)]
            )

        primers: list[dict[str, Any]] = []
        for c in self.primers_list.controls:
            if isinstance(c, ft.ListTile) and isinstance(c.data, dict):
                is_active = True
                if isinstance(c.leading, ft.Checkbox):
                    is_active = bool(c.leading.value)
                primers.append(
                    {
                        "name": str(c.data.get("name", "")),
                        "seq": format_seq(str(c.data.get("seq", ""))),
                        "active": is_active,
                    }
                )
        return primers

    def get_state(self) -> dict[str, Any]:
        """Get the current state for serialization."""

        def format_seq(seq_val: str | None) -> str:
            if not seq_val:
                return ""
            # Remove literal escaped whitespace characters
            clean_seq = (
                str(seq_val)
                .replace("\\n", "")
                .replace("\\t", "")
                .replace("\\r", "")
            )
            clean_seq = "".join(clean_seq.split()).upper()
            return "\n".join(
                [clean_seq[i : i + 80] for i in range(0, len(clean_seq), 80)]
            )

        return {
            "template": format_seq(self.template_sequence.value),
            "template_circular": self.template_circular.value,
            "primers": self.get_all_primers_state(),
        }

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current state from deserialized data."""
        if "template" in state:
            self.template_sequence.value = state["template"].replace("\n", "")

        if "template_circular" in state:
            self.template_circular.value = state["template_circular"]

        if "primers" in state:
            self.primers_list.controls.clear()
            for p in state["primers"]:
                p_seq = p.get("seq", "")
                seq_str = (
                    "".join(p_seq) if isinstance(p_seq, list) else p_seq
                ).replace("\n", "")
                is_active = p.get("active", True)
                self.primers_list.controls.append(
                    ft.ListTile(
                        dense=True,
                        content_padding=ft.padding.all(0),
                        leading=ft.Checkbox(value=is_active),
                        title=ft.Text(f"{p['name']}, {seq_str}"),
                        data={"name": p["name"], "seq": seq_str},
                        trailing=ft.IconButton(
                            ft.Icons.DELETE,
                            on_click=lambda e, n=p["name"], s=seq_str: (
                                self.remove_primer(e, n, s)
                            ),
                        ),
                    )
                )
        self.app_page.update()
