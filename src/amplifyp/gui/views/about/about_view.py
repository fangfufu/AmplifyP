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

"""About View for the Flet application."""

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings


class AboutView(ft.ListView):  # type: ignore[misc]
    """About view for showing application version info."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings | None = None,
    ) -> None:
        """Initialise the AboutView."""
        super().__init__(
            expand=True, spacing=20, padding=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.app_page = page
        self.settings = settings if settings is not None else GUISettings()

        from amplifyp.gui.util import get_full_sha, get_version

        full_sha = get_full_sha()
        app_version = get_version()
        font_size_default = self.settings.get("font_size_default", 14)

        # ponytail: simple ListView containing icon, version, and commit SHA
        self.controls = [
            ft.Text(
                "About AmplifyP",
                weight=ft.FontWeight.BOLD,
                size=24,
            ),
            ft.Divider(),
            ft.Row(
                [
                    ft.Image(
                        src="/images/favicon.png",
                        height=64,
                        fit=ft.BoxFit.CONTAIN,
                    ),
                    ft.Column(
                        [
                            ft.Text(
                                "AmplifyP",
                                size=20,
                                weight=ft.FontWeight.BOLD,
                            ),
                            ft.Text(
                                "A GUI application for PCR replication "
                                "simulation and primer-dimer analysis.",
                                size=font_size_default,
                            ),
                        ],
                        spacing=5,
                    ),
                ],
                spacing=20,
            ),
            ft.Divider(),
            ft.Container(
                content=ft.Column(
                    [
                        ft.Row(
                            [
                                ft.Text(
                                    "Version: ",
                                    weight=ft.FontWeight.BOLD,
                                    size=font_size_default,
                                ),
                                ft.Text(
                                    app_version,
                                    selectable=True,
                                    size=font_size_default,
                                ),
                            ],
                            alignment=ft.MainAxisAlignment.START,
                            vertical_alignment=ft.CrossAxisAlignment.CENTER,
                        ),
                        ft.Row(
                            [
                                ft.Text(
                                    "Repository: ",
                                    weight=ft.FontWeight.BOLD,
                                    size=font_size_default,
                                ),
                                ft.Text(
                                    spans=[
                                        ft.TextSpan(
                                            "github.com/fangfufu/AmplifyP",
                                            url="https://github.com/fangfufu/AmplifyP",
                                            style=ft.TextStyle(
                                                decoration=ft.TextDecoration.UNDERLINE,
                                                color=GUIColours.LINK_BLUE,
                                            ),
                                        )
                                    ],
                                    size=font_size_default,
                                ),
                            ],
                            alignment=ft.MainAxisAlignment.START,
                            vertical_alignment=ft.CrossAxisAlignment.CENTER,
                        ),
                        ft.Row(
                            [
                                ft.Text(
                                    "Full Git SHA: ",
                                    weight=ft.FontWeight.BOLD,
                                    size=font_size_default,
                                ),
                                ft.Text(
                                    full_sha,
                                    selectable=True,
                                    size=font_size_default,
                                ),
                            ],
                            alignment=ft.MainAxisAlignment.START,
                            vertical_alignment=ft.CrossAxisAlignment.CENTER,
                        ),
                    ],
                    spacing=10,
                ),
                padding=ft.Padding(10, 5, 10, 5),
            ),
            ft.Divider(),
            ft.Text(
                "Attribution",
                weight=ft.FontWeight.BOLD,
                size=18,
            ),
            ft.Text(
                "• Amplify4: This project is built upon the logic "
                "and methodology of the "
                "original Amplify4 software by William Engels. We "
                "preserve the simulation "
                "models and algorithms of the original while "
                "offering a modern, robust, "
                "and accessible cross-platform implementation.",
                size=font_size_default,
            ),
            ft.Text(
                "• Roboto Mono Font: Licenced under the Apache "
                "License, Version 2.0. "
                "Copyright 2015 The Roboto Mono Project Authors.",
                size=font_size_default,
            ),
            ft.Divider(),
            ft.Text(
                "Licence",
                weight=ft.FontWeight.BOLD,
                size=18,
            ),
            ft.Text(
                "This project is licenced under the GNU General Public "
                "License v3.0 (GPL-3.0).",
                size=font_size_default,
            ),
        ]
