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

"""Centralized GUI color definitions and dynamic color shifting."""

from typing import cast

import flet as ft


class _GUIColorsMeta(type):
    """Metaclass for GUIColors to support dynamic color shifting."""

    @property
    def color_deficient_mode(cls) -> bool:
        """Get color deficient mode status."""
        return cls._color_deficient_mode

    @color_deficient_mode.setter
    def color_deficient_mode(cls, value: bool) -> None:
        """Set color deficient mode status."""
        cls._color_deficient_mode = value

    @property
    def dark_mode(cls) -> bool:
        """Get dark mode status."""
        return cls._dark_mode

    @dark_mode.setter
    def dark_mode(cls, value: bool) -> None:
        """Set dark mode status."""
        cls._dark_mode = value

    @property
    def TEXT_ON_SURFACE(cls) -> str:
        """Get standard text color on surface."""
        return cast(str, ft.Colors.ON_SURFACE)

    @property
    def SUCCESS_GREEN(cls) -> str:
        """Get color-blind friendly success color."""
        # Blue is highly distinguishable for most common
        # red-green color blindness.
        return cast(
            str,
            ft.Colors.BLUE_400
            if cls._color_deficient_mode
            else ft.Colors.GREEN_400,
        )

    @property
    def SURFACE(cls) -> str:
        """Get standard surface color."""
        return cast(str, ft.Colors.SURFACE)

    @property
    def OUTLINE_VARIANT(cls) -> str:
        """Get outline variant color."""
        return cast(str, ft.Colors.OUTLINE_VARIANT)

    @property
    def OUTLINE(cls) -> str:
        """Get outline color."""
        return cast(str, ft.Colors.OUTLINE)

    @property
    def ERROR_RED(cls) -> str:
        """Get color-blind friendly error color."""
        # Use Orange/Vermilion instead of Red in color deficient mode.
        return cast(
            str,
            ft.Colors.ORANGE_800
            if cls._color_deficient_mode
            else ft.Colors.RED,
        )

    @property
    def DIVIDER_GREY(cls) -> str:
        """Get divider grey color."""
        return cast(str, ft.Colors.GREY_400)

    @property
    def MUTED_GREY(cls) -> str:
        """Get muted/greyed out text color."""
        return cast(str, ft.Colors.GREY_500)

    @property
    def DUPLICATE_BG(cls) -> str:
        """Get duplicate warning background color."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.ORANGE_900
                if cls._color_deficient_mode
                else ft.Colors.RED_900,
            )
        return cast(
            str,
            ft.Colors.ORANGE_100
            if cls._color_deficient_mode
            else ft.Colors.RED_100,
        )

    @property
    def SELECTED_ROW_BG(cls) -> str:
        """Get selected/focused row background color."""
        return cast(
            str,
            ft.Colors.BLUE_900 if cls._dark_mode else ft.Colors.BLUE_50,
        )

    @property
    def DIAGRAM_BLACK(cls) -> str:
        """Get diagram black/white color depending on dark mode."""
        return cast(
            str,
            ft.Colors.WHITE if cls._dark_mode else ft.Colors.BLACK,
        )

    @property
    def INFO_HEADER_BG(cls) -> str:
        """Get background color for info header."""
        return cast(
            str,
            ft.Colors.GREY_800 if cls._dark_mode else ft.Colors.GREY_200,
        )

    @property
    def FWD_PRIMER(cls) -> str:
        """Get forward primer color."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.BLUE_300
                if cls._color_deficient_mode
                else ft.Colors.BLUE_400,
            )
        # Sky blue / clear blue for forward primer in color deficient mode
        return cast(
            str,
            ft.Colors.BLUE_600
            if cls._color_deficient_mode
            else ft.Colors.BLUE_800,
        )

    @property
    def REV_PRIMER(cls) -> str:
        """Get reverse primer color."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.ORANGE_300
                if cls._color_deficient_mode
                else ft.Colors.RED_ACCENT_200,
            )
        # Vermilion / orange-red for reverse primer in color deficient mode
        # (instead of red-accent)
        return cast(
            str,
            ft.Colors.ORANGE_700
            if cls._color_deficient_mode
            else ft.Colors.RED_ACCENT_700,
        )

    @property
    def REV_LABEL(cls) -> str:
        """Get reverse primer label color."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.ORANGE_200
                if cls._color_deficient_mode
                else ft.Colors.RED_300,
            )
        return cast(
            str,
            ft.Colors.ORANGE_900
            if cls._color_deficient_mode
            else ft.Colors.RED_800,
        )

    @property
    def WHITE(cls) -> str:
        """Get standard white color."""
        return cast(str, ft.Colors.WHITE)

    @property
    def TRANSPARENT(cls) -> str:
        """Get transparent color."""
        return cast(str, ft.Colors.TRANSPARENT)


class GUIColors(metaclass=_GUIColorsMeta):
    """Centralized semantic color constants for the GUI."""

    _color_deficient_mode = False
    _dark_mode = False
