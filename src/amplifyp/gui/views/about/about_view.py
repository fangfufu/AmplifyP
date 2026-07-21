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

        from amplifyp.gui.utils.system import get_full_sha, get_version

        full_sha = get_full_sha()
        app_version = get_version()
        font_size_default = self.settings.get("font_size_default", 14)

        self.controls = [
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
                        ft.Row(
                            [
                                ft.Text(
                                    "Licence: ",
                                    weight=ft.FontWeight.BOLD,
                                    size=font_size_default,
                                ),
                                ft.Text(
                                    spans=[
                                        ft.TextSpan(
                                            "This project is licenced "
                                            "under the "
                                        ),
                                        ft.TextSpan(
                                            "GNU General Public "
                                            "License v3.0 (GPL-3.0)",
                                            url="https://github.com/"
                                            "fangfufu/AmplifyP/blob/"
                                            "stable/LICENSE",
                                            style=ft.TextStyle(
                                                decoration=ft.TextDecoration.UNDERLINE,
                                                color=GUIColours.LINK_BLUE,
                                            ),
                                        ),
                                        ft.TextSpan("."),
                                    ],
                                    size=font_size_default,
                                ),
                            ]
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
                spans=[
                    ft.TextSpan("• "),
                    ft.TextSpan(
                        "Amplify 4",
                        url="https://engels.genetics.wisc.edu/amplify/",
                        style=ft.TextStyle(
                            decoration=ft.TextDecoration.UNDERLINE,
                            color=GUIColours.LINK_BLUE,
                        ),
                    ),
                    ft.TextSpan(
                        ": This project is based on the original "
                        "software by William Engels."
                    ),
                ],
                size=font_size_default,
            ),
            ft.Text(
                spans=[
                    ft.TextSpan("• "),
                    ft.TextSpan(
                        "Roboto Mono Font",
                        url="https://fonts.google.com/specimen/Roboto+Mono",
                        style=ft.TextStyle(
                            decoration=ft.TextDecoration.UNDERLINE,
                            color=GUIColours.LINK_BLUE,
                        ),
                    ),
                    ft.TextSpan(
                        ": Licenced under the SIL Open Font License, "
                        "Version 1.1. Copyright 2015 The Roboto Mono "
                        "Project Authors."
                    ),
                ],
                size=font_size_default,
            ),
        ]
