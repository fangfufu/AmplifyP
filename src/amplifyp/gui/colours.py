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

"""Centralised GUI colour definitions and dynamic colour shifting."""

from typing import ClassVar, cast

import flet as ft


class _GUIColoursMeta(type):
    """Metaclass for GUIColours to support dynamic colour shifting."""

    @property
    def colour_deficient_mode(cls) -> bool:
        """Get colour deficient mode status."""
        return cls._colour_deficient_mode

    @colour_deficient_mode.setter
    def colour_deficient_mode(cls, value: bool) -> None:
        """Set colour deficient mode status."""
        cls._colour_deficient_mode = value

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
        """Get standard text colour on surface."""
        return cast(str, ft.Colors.ON_SURFACE)

    @property
    def SUCCESS_GREEN(cls) -> str:
        """Get colour-blind friendly success colour."""
        # Blue is highly distinguishable for most common
        # red-green colour blindness.
        return cast(
            str,
            ft.Colors.BLUE_400
            if cls._colour_deficient_mode
            else ft.Colors.GREEN_400,
        )

    @property
    def SURFACE(cls) -> str:
        """Get standard surface colour."""
        return cast(str, ft.Colors.SURFACE)

    @property
    def OUTLINE_VARIANT(cls) -> str:
        """Get outline variant colour."""
        return cast(str, ft.Colors.OUTLINE_VARIANT)

    @property
    def OUTLINE(cls) -> str:
        """Get outline colour."""
        return cast(str, ft.Colors.OUTLINE)

    @property
    def ERROR_RED(cls) -> str:
        """Get colour-blind friendly error colour."""
        # Use Orange/Vermilion instead of Red in colour deficient mode.
        return cast(
            str,
            ft.Colors.ORANGE_800
            if cls._colour_deficient_mode
            else ft.Colors.RED,
        )

    @property
    def DIVIDER_GREY(cls) -> str:
        """Get divider grey colour."""
        return cast(str, ft.Colors.GREY_400)

    @property
    def MUTED_GREY(cls) -> str:
        """Get muted/greyed out text colour."""
        return cast(str, ft.Colors.GREY_500)

    @property
    def DUPLICATE_BG(cls) -> str:
        """Get duplicate warning background colour."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.ORANGE_900
                if cls._colour_deficient_mode
                else ft.Colors.RED_900,
            )
        return cast(
            str,
            ft.Colors.ORANGE_100
            if cls._colour_deficient_mode
            else ft.Colors.RED_100,
        )

    @property
    def SELECTED_ROW_BG(cls) -> str:
        """Get selected/focused row background colour."""
        return cast(
            str,
            ft.Colors.BLUE_900 if cls._dark_mode else ft.Colors.BLUE_50,
        )

    @property
    def DIAGRAM_BLACK(cls) -> str:
        """Get diagram black/white colour depending on dark mode."""
        return cast(str, ft.Colors.ON_SURFACE)

    @property
    def INFO_HEADER_BG(cls) -> str:
        """Get background colour for info header."""
        return cast(
            str,
            ft.Colors.GREY_800 if cls._dark_mode else ft.Colors.GREY_200,
        )

    @property
    def FWD_PRIMER(cls) -> str:
        """Get forward primer colour."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.BLUE_300
                if cls._colour_deficient_mode
                else ft.Colors.BLUE_400,
            )
        # Sky blue / clear blue for forward primer in colour deficient mode
        return cast(
            str,
            ft.Colors.BLUE_600
            if cls._colour_deficient_mode
            else ft.Colors.BLUE_800,
        )

    @property
    def REV_PRIMER(cls) -> str:
        """Get reverse primer colour."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.ORANGE_300
                if cls._colour_deficient_mode
                else ft.Colors.RED_ACCENT_200,
            )
        # Vermilion / orange-red for reverse primer in colour deficient mode
        # (instead of red-accent)
        return cast(
            str,
            ft.Colors.ORANGE_700
            if cls._colour_deficient_mode
            else ft.Colors.RED_ACCENT_700,
        )

    @property
    def REV_LABEL(cls) -> str:
        """Get reverse primer label colour."""
        if cls._dark_mode:
            return cast(
                str,
                ft.Colors.ORANGE_200
                if cls._colour_deficient_mode
                else ft.Colors.RED_300,
            )
        return cast(
            str,
            ft.Colors.ORANGE_900
            if cls._colour_deficient_mode
            else ft.Colors.RED_800,
        )

    @property
    def LINK_BLUE(cls) -> str:
        """Get standard hyperlink blue colour."""
        return cast(
            str,
            ft.Colors.BLUE_400 if cls._dark_mode else ft.Colors.BLUE_700,
        )

    @property
    def WHITE(cls) -> str:
        """Get standard white colour."""
        return cast(str, ft.Colors.WHITE)

    @property
    def TRANSPARENT(cls) -> str:
        """Get transparent colour."""
        return cast(str, ft.Colors.TRANSPARENT)

    @property
    def GRADIENT_GREY(cls) -> str:
        """Get grey colour for temperature gradient."""
        return cast(str, ft.Colors.GREY_500)

    @property
    def GRADIENT_GREEN(cls) -> str:
        """Get green/blue gradient colour for temperature."""
        return cast(
            str,
            ft.Colors.BLUE_600
            if cls._colour_deficient_mode
            else ft.Colors.GREEN_600,
        )

    @property
    def GRADIENT_YELLOW(cls) -> str:
        """Get yellow/orange gradient colour for temperature."""
        return cast(
            str,
            ft.Colors.ORANGE_600
            if cls._colour_deficient_mode
            else ft.Colors.YELLOW_600,
        )

    @property
    def GRADIENT_RED(cls) -> str:
        """Get red colour for temperature gradient."""
        return cast(str, ft.Colors.RED_700)

    @property
    def GRADIENT_BLUE(cls) -> str:
        """Get blue colour for Cool-Warm temperature gradient."""
        return cast(str, ft.Colors.BLUE_700)

    @property
    def GRADIENT_MIDPOINT(cls) -> str:
        """Get midpoint colour for Cool-Warm temperature gradient."""
        return cast(str, ft.Colors.WHITE if cls._dark_mode else ft.Colors.BLACK)


class GUIColours(metaclass=_GUIColoursMeta):
    """Centralised semantic colour constants for the GUI."""

    _colour_deficient_mode = False
    _dark_mode = False


class ColourInterpolator:
    """Helper to convert color names to hex codes and interpolate colors."""

    MAP: ClassVar[dict[str, str]] = {
        "grey500": "#9e9e9e",
        "blue600": "#1e88e5",
        "green600": "#388e3c",
        "orange600": "#f57c00",
        "yellow600": "#fdd835",
        "red700": "#d32f2f",
        "blue700": "#1565c0",
        "white": "#ffffff",
        "black": "#000000",
    }

    @classmethod
    def get_hex(cls, color_str: str) -> str:
        """Get hex representation of a color name or return hex string."""
        if color_str.startswith("#"):
            return color_str
        norm_name = color_str.lower().replace("_", "")
        return cls.MAP.get(norm_name, color_str)

    @classmethod
    def to_rgb(cls, color_str: str) -> tuple[int, int, int]:
        """Convert a colour name or hex string to an (R, G, B) tuple."""
        hex_clean = cls.get_hex(color_str).lstrip("#")
        if len(hex_clean) == 3:
            hex_clean = "".join(c * 2 for c in hex_clean)
        if len(hex_clean) != 6:
            raise ValueError(f"Invalid hex colour: {color_str}")
        return (
            int(hex_clean[0:2], 16),
            int(hex_clean[2:4], 16),
            int(hex_clean[4:6], 16),
        )

    @staticmethod
    def to_hex(r: int, g: int, b: int) -> str:
        """Convert an (R, G, B) tuple to a hex colour string."""
        return f"#{r:02x}{g:02x}{b:02x}"

    @classmethod
    def interpolate(cls, color_1: str, color_2: str, t: float) -> str:
        """Interpolate between two colors by factor t (0 to 1)."""
        r1, g1, b1 = cls.to_rgb(color_1)
        r2, g2, b2 = cls.to_rgb(color_2)
        r = int(r1 + (r2 - r1) * t)
        g = int(g1 + (g2 - g1) * t)
        b = int(b1 + (b2 - b1) * t)
        return cls.to_hex(r, g, b)


def tm_colour(tm_value: float | None, scheme: str) -> str | None:
    """Return a Flet colour or hex string for the given Tm.

    Returns None when scheme is "None" or unrecognised.

    Args:
        tm_value: The melting temperature in degrees Celsius, or None.
        scheme: One of "None", "Cool-Warm", or "Traffic Light".

    Returns:
        A Flet Colors enum value/hex string, or None.
    """
    if tm_value is None or scheme == "None":
        return None

    if scheme == "Traffic Light":
        if tm_value <= 0:
            return ColourInterpolator.get_hex(GUIColours.GRADIENT_GREY)
        if tm_value < 55:
            t = tm_value / 55.0
            return ColourInterpolator.interpolate(
                GUIColours.GRADIENT_GREY, GUIColours.GRADIENT_GREEN, t
            )
        if tm_value < 65:
            t = (tm_value - 55.0) / 10.0
            return ColourInterpolator.interpolate(
                GUIColours.GRADIENT_GREEN, GUIColours.GRADIENT_YELLOW, t
            )
        if tm_value < 75:
            t = (tm_value - 65.0) / 10.0
            return ColourInterpolator.interpolate(
                GUIColours.GRADIENT_YELLOW, GUIColours.GRADIENT_RED, t
            )
        return ColourInterpolator.get_hex(GUIColours.GRADIENT_RED)

    if scheme == "Cool-Warm":
        if tm_value <= 35:
            return ColourInterpolator.get_hex(GUIColours.GRADIENT_BLUE)
        if tm_value < 55:
            t = (tm_value - 35.0) / 20.0
            return ColourInterpolator.interpolate(
                GUIColours.GRADIENT_BLUE, GUIColours.GRADIENT_MIDPOINT, t
            )
        if tm_value < 75:
            t = (tm_value - 55.0) / 20.0
            return ColourInterpolator.interpolate(
                GUIColours.GRADIENT_MIDPOINT, GUIColours.GRADIENT_RED, t
            )
        return ColourInterpolator.get_hex(GUIColours.GRADIENT_RED)

    return None
