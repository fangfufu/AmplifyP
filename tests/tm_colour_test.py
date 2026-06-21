"""Unit tests for the tm_colour() colour-mapping helper."""

from collections.abc import Generator

import flet as ft
import pytest

from amplifyp.gui.colours import GUIColours, tm_colour


@pytest.fixture(autouse=True)  # type: ignore[untyped-decorator]
def reset_colour_modes() -> Generator[None, None, None]:
    """Reset colour-deficient and dark mode between tests."""
    GUIColours.colour_deficient_mode = False
    GUIColours.dark_mode = False
    yield
    GUIColours.colour_deficient_mode = False
    GUIColours.dark_mode = False


class TestTmColourNone:
    """Tests for the "None" scheme."""

    def test_none_scheme_returns_none_for_any_tm(self) -> None:
        """Scheme "None" always returns None regardless of temperature."""
        assert tm_colour(30.0, "None") is None
        assert tm_colour(55.0, "None") is None
        assert tm_colour(80.0, "None") is None

    def test_none_tm_value_returns_none(self) -> None:
        """None tm_value always returns None regardless of scheme."""
        assert tm_colour(None, "Traffic Light") is None
        assert tm_colour(None, "Cool-Warm") is None

    def test_unrecognised_scheme_returns_none(self) -> None:
        """Unrecognised scheme strings return None."""
        assert tm_colour(55.0, "Rainbow") is None
        assert tm_colour(55.0, "") is None


class TestTmColourTrafficLight:
    """Tests for the "Traffic Light" scheme."""

    def test_low_tm_gives_grey(self) -> None:
        """Tm <= 0 deg C gives grey in normal mode."""
        assert tm_colour(-5.0, "Traffic Light") == "#9e9e9e"
        assert tm_colour(0.0, "Traffic Light") == "#9e9e9e"

    def test_green_midpoint(self) -> None:
        """Tm = 55 deg C gives green in normal mode."""
        assert tm_colour(55.0, "Traffic Light") == "#388e3c"

    def test_yellow_midpoint(self) -> None:
        """Tm = 65 deg C gives yellow in normal mode."""
        assert tm_colour(65.0, "Traffic Light") == "#fdd835"

    def test_red_midpoint(self) -> None:
        """Tm >= 75 deg C gives red in normal mode."""
        assert tm_colour(75.0, "Traffic Light") == "#d32f2f"
        assert tm_colour(80.0, "Traffic Light") == "#d32f2f"

    def test_interpolation_grey_to_green(self) -> None:
        """Tm between 0 and 55 interpolates between grey and green."""
        # Halfway point at 27.5 deg C (t = 0.5)
        # Grey: #9e9e9e (158, 158, 158)
        # Green: #388e3c (56, 142, 60)
        # R = 158 + (56 - 158) * 0.5 = 107
        # G = 158 + (142 - 158) * 0.5 = 150
        # B = 158 + (60 - 158) * 0.5 = 109
        # Hex: #6b966d
        assert tm_colour(27.5, "Traffic Light") == "#6b966d"

    def test_interpolation_green_to_yellow(self) -> None:
        """Tm between 55 and 65 interpolates between green and yellow."""
        # Halfway point at 60 deg C (t = 0.5)
        # Green: #388e3c (56, 142, 60)
        # Yellow: #fdd835 (253, 216, 53)
        # R = 56 + (253 - 56) * 0.5 = 154
        # G = 142 + (216 - 142) * 0.5 = 179
        # B = 60 + (53 - 60) * 0.5 = 56
        # Hex: #9ab338
        assert tm_colour(60.0, "Traffic Light") == "#9ab338"

    def test_interpolation_yellow_to_red(self) -> None:
        """Tm between 65 and 75 interpolates between yellow and red."""
        # Halfway point at 70 deg C (t = 0.5)
        # Yellow: #fdd835 (253, 216, 53)
        # Red: #d32f2f (211, 47, 47)
        # R = 253 + (211 - 253) * 0.5 = 232
        # G = 216 + (47 - 216) * 0.5 = 131
        # B = 53 + (47 - 53) * 0.5 = 50
        # Hex: #e88332
        assert tm_colour(70.0, "Traffic Light") == "#e88332"

    def test_colour_deficient_grey_boundary(self) -> None:
        """Colour-deficient mode: Tm <= 0 deg C still gives grey."""
        GUIColours.colour_deficient_mode = True
        assert tm_colour(-5.0, "Traffic Light") == "#9e9e9e"

    def test_colour_deficient_blue_midpoint(self) -> None:
        """Colour-deficient mode: Tm = 55 deg C gives blue."""
        GUIColours.colour_deficient_mode = True
        assert tm_colour(55.0, "Traffic Light") == "#1e88e5"

    def test_colour_deficient_orange_midpoint(self) -> None:
        """Colour-deficient mode: Tm = 65 deg C gives orange."""
        GUIColours.colour_deficient_mode = True
        assert tm_colour(65.0, "Traffic Light") == "#f57c00"

    def test_colour_deficient_red_midpoint(self) -> None:
        """Colour-deficient mode: Tm >= 75 deg C gives red."""
        GUIColours.colour_deficient_mode = True
        assert tm_colour(75.0, "Traffic Light") == "#d32f2f"

    def test_colour_deficient_interpolation_blue_to_orange(self) -> None:
        """Colour-deficient mode: interpolates between blue and orange."""
        # Halfway point at 60 deg C (t = 0.5)
        # Blue: #1e88e5 (30, 136, 229)
        # Orange: #f57c00 (245, 124, 0)
        # R = 30 + (245 - 30) * 0.5 = 137
        # G = 136 + (124 - 136) * 0.5 = 130
        # B = 229 + (0 - 229) * 0.5 = 114
        # Hex: #898272
        GUIColours.colour_deficient_mode = True
        assert tm_colour(60.0, "Traffic Light") == "#898272"


class TestTmColourCoolWarm:
    """Tests for the "Cool-Warm" scheme."""

    def test_very_cold_gives_blue_700(self) -> None:
        """Tm well below 45 deg C gives the coldest blue."""
        assert tm_colour(30.0, "Cool-Warm") is ft.Colors.BLUE_700

    def test_boundary_45_gives_blue_500(self) -> None:
        """Tm exactly at 45 deg C is in the 45-50 band (BLUE_500)."""
        assert tm_colour(45.0, "Cool-Warm") is ft.Colors.BLUE_500

    def test_mid_low_gives_blue_300(self) -> None:
        """Tm in 50-55 band gives BLUE_300."""
        assert tm_colour(52.0, "Cool-Warm") is ft.Colors.BLUE_300

    def test_near_midpoint_gives_on_surface_variant(self) -> None:
        """Tm in the 58-62 midpoint band gives ON_SURFACE_VARIANT."""
        assert tm_colour(60.0, "Cool-Warm") is ft.Colors.ON_SURFACE_VARIANT

    def test_warm_gives_red_300(self) -> None:
        """Tm in the 65-70 band gives RED_300."""
        assert tm_colour(67.0, "Cool-Warm") is ft.Colors.RED_300

    def test_very_hot_gives_red_700(self) -> None:
        """Tm >= 75 deg C gives the hottest red."""
        assert tm_colour(80.0, "Cool-Warm") is ft.Colors.RED_700
