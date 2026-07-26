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

"""Centralised GUI colour definitions and dynamic colour shifting."""

from __future__ import annotations

from typing import ClassVar


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
        """Set dark mode mode status."""
        cls._dark_mode = value

    @property
    def TEXT_ON_SURFACE(cls) -> str:
        """Get standard text colour on surface."""
        return "#000000" if not cls._dark_mode else "#ffffff"

    @property
    def SUCCESS_GREEN(cls) -> str:
        """Get colour-blind friendly success colour."""
        return "#1e88e5" if cls.colour_deficient_mode else "#388e3c"

    @property
    def SURFACE(cls) -> str:
        """Get standard surface colour."""
        return "#ffffff" if not cls._dark_mode else "#121212"

    @property
    def SURFACE_VARIANT(cls) -> str:
        """Get standard surface variant colour."""
        return "#f5f5f5" if not cls._dark_mode else "#1e1e1e"

    @property
    def PRIMARY(cls) -> str:
        """Get primary theme colour."""
        return "#6750a4" if not cls._dark_mode else "#d0bcff"

    @property
    def OUTLINE_VARIANT(cls) -> str:
        """Get outline variant colour."""
        return "#9e9e9e" if cls._dark_mode else "#e0e0e0"

    @property
    def OUTLINE(cls) -> str:
        """Get outline colour."""
        return "#757575" if cls._dark_mode else "#9e9e9e"

    @property
    def ERROR_RED(cls) -> str:
        """Get colour-blind friendly error colour."""
        return "#e65100" if cls.colour_deficient_mode else "#d32f2f"

    @property
    def DIVIDER_GREY(cls) -> str:
        """Get divider grey colour."""
        return "#bdbdbd"

    @property
    def UPDATE_AVAILABLE_COLOUR(cls) -> str:
        """Get colour-blind friendly update available notification colour."""
        return "#29b6f6" if cls.colour_deficient_mode else "#69f0ae"

    @property
    def MUTED_GREY(cls) -> str:
        """Get muted/greyed out text colour."""
        return "#757575"

    @property
    def DUPLICATE_BG(cls) -> str:
        """Get duplicate warning background colour."""
        if cls._dark_mode:
            return "#e65100" if cls.colour_deficient_mode else "#b71c1c"
        return "#fff3e0" if cls.colour_deficient_mode else "#ffebee"

    @property
    def FOCUSED_DUPLICATE_BG(cls) -> str:
        """Get background colour for a duplicate warning row that is focused."""
        if cls._dark_mode:
            return "#e64a19" if cls.colour_deficient_mode else "#c62828"
        return "#ffccbc" if cls.colour_deficient_mode else "#ffcdd2"

    @property
    def SELECTED_ROW_BG(cls) -> str:
        """Get selected/focused row background colour."""
        return "#1a237e" if cls._dark_mode else "#d0ebfc"

    @property
    def DIAGRAM_BLACK(cls) -> str:
        """Get diagram black/white colour depending on dark mode."""
        return "#ffffff" if cls._dark_mode else "#000000"

    @property
    def DIAGRAM_BG(cls) -> str:
        """Get diagram background colour depending on dark mode."""
        return "#212121" if cls._dark_mode else "#f5f5f5"

    @property
    def INFO_HEADER_BG(cls) -> str:
        """Get background colour for info header."""
        return "#1e1e1e" if cls._dark_mode else "#eeeeee"

    @property
    def GUTTER_BG(cls) -> str:
        """Get background colour for template display gutter."""
        return "#1e1e1e" if cls._dark_mode else "#eeeeee"

    @property
    def FWD_PRIMER(cls) -> str:
        """Get forward primer colour."""
        if cls._dark_mode:
            return "#64b5f6" if cls.colour_deficient_mode else "#42a5f5"
        return "#1565c0" if cls.colour_deficient_mode else "#0d47a1"

    @property
    def REV_PRIMER(cls) -> str:
        """Get reverse primer colour."""
        if cls._dark_mode:
            return "#ffb74d" if cls.colour_deficient_mode else "#ff8a65"
        return "#e65100" if cls.colour_deficient_mode else "#c62828"

    @property
    def REV_LABEL(cls) -> str:
        """Get reverse primer label colour."""
        if cls._dark_mode:
            return "#ffcc80" if cls.colour_deficient_mode else "#ff8a65"
        return "#bf360c" if cls.colour_deficient_mode else "#880e4f"

    @property
    def LINK_BLUE(cls) -> str:
        """Get standard hyperlink blue colour."""
        return "#42a5f5" if cls._dark_mode else "#1565c0"

    @property
    def PURPLE(cls) -> str:
        """Get standard purple colour."""
        return "#ce93d8" if cls._dark_mode else "#7b1fa2"

    @property
    def WHITE(cls) -> str:
        """Get standard white colour."""
        return "#ffffff"

    @property
    def TRANSPARENT(cls) -> str:
        """Get transparent colour."""
        return "transparent"

    @property
    def GRADIENT_GREY(cls) -> str:
        """Get grey colour for temperature gradient."""
        return "#9e9e9e"

    @property
    def GRADIENT_GREEN(cls) -> str:
        """Get green/blue gradient colour for temperature."""
        return "#1e88e5" if cls.colour_deficient_mode else "#388e3c"

    @property
    def GRADIENT_YELLOW(cls) -> str:
        """Get yellow/orange gradient colour for temperature."""
        return "#f57c00" if cls.colour_deficient_mode else "#fdd835"

    @property
    def GRADIENT_RED(cls) -> str:
        """Get red colour for temperature gradient."""
        return "#d32f2f"

    @property
    def GRADIENT_BLUE(cls) -> str:
        """Get blue colour for Cool-Warm temperature gradient."""
        return "#1565c0"

    @property
    def GRADIENT_MIDPOINT(cls) -> str:
        """Get midpoint colour for Cool-Warm temperature gradient."""
        return "#000000" if cls._dark_mode else "#ffffff"

    @property
    def GRID_2D_BLUE(cls) -> str:
        """Get soft desaturated blue colour for 2D grid results colour map."""
        return "#2d4a68" if cls._dark_mode else "#d0e1f9"

    @property
    def GRID_2D_ORANGE(cls) -> str:
        """Get soft desaturated orange colour for 2D grid results colour map."""
        return "#68482d" if cls._dark_mode else "#fdebd0"

    @property
    def GRID_2D_GREEN(cls) -> str:
        """Get soft desaturated green colour for 2D grid results colour map."""
        if cls.colour_deficient_mode:
            return cls.GRID_2D_BLUE
        return "#2d683b" if cls._dark_mode else "#d0f9d0"

    @property
    def GRID_2D_YELLOW(cls) -> str:
        """Get soft desaturated yellow colour for 2D grid results colour map."""
        if cls.colour_deficient_mode:
            return cls.GRID_2D_ORANGE
        return "#68632d" if cls._dark_mode else "#fdf6d0"

    @property
    def GRID_2D_RED(cls) -> str:
        """Get soft desaturated red colour for 2D grid results colour map."""
        return "#682d2d" if cls._dark_mode else "#f9d0d0"


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
    """Return a hex colour string for the given Tm.

    Returns None when scheme is "None" or unrecognised.

    Args:
        tm_value: The melting temperature in degrees Celsius, or None.
        scheme: One of "None", "Cool-Warm", or "Traffic Light".

    Returns:
        A hex colour string, or None.
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


def get_text_contrast_colour(bg_hex: str) -> str:
    """Return high contrast text colour (dark or white) for a background hex.

    Args:
        bg_hex: Background hex colour string.

    Returns:
        Hex colour string for high contrast text.
    """
    try:
        r, g, b = ColourInterpolator.to_rgb(bg_hex)
        luminance = 0.299 * r + 0.587 * g + 0.114 * b
        return "#ffffff" if luminance < 128 else "#111827"
    except ValueError:
        return GUIColours.TEXT_ON_SURFACE


def designer_2d_colour(
    value: float | None, min_val: float, max_val: float, scheme: str
) -> str | None:
    """Return a hex colour string for a 2D grid value using specified scheme.

    Inverts the scale so that lower quality scores (which represent better
    primer pairs with lower dimerisation risk) receive hotter / higher
    intensity colours.

    Args:
        value: The numerical quality/metric score to colour.
        min_val: Minimum metric value across the grid.
        max_val: Maximum metric value across the grid.
        scheme: One of "None", "Cool-Warm", "Traffic Light", or "Blue-Orange".

    Returns:
        Hex colour string, or None if scheme is disabled or invalid.
    """
    if value is None or scheme == "None":
        return None

    if max_val <= min_val:
        t = 0.5
    else:
        t = 1.0 - max(0.0, min(1.0, (value - min_val) / (max_val - min_val)))

    if scheme == "Cool-Warm":
        return ColourInterpolator.interpolate(
            GUIColours.GRID_2D_BLUE, GUIColours.GRID_2D_RED, t
        )

    if scheme == "Traffic Light":
        if t < 0.5:
            norm_t = t * 2.0
            return ColourInterpolator.interpolate(
                GUIColours.GRID_2D_RED, GUIColours.GRID_2D_YELLOW, norm_t
            )
        norm_t = (t - 0.5) * 2.0
        return ColourInterpolator.interpolate(
            GUIColours.GRID_2D_YELLOW, GUIColours.GRID_2D_GREEN, norm_t
        )

    if scheme == "Blue-Orange":
        return ColourInterpolator.interpolate(
            GUIColours.GRID_2D_BLUE, GUIColours.GRID_2D_ORANGE, t
        )

    return None
