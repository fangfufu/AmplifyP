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

    def test_high_tm_gives_green(self) -> None:
        """Tm >= 60 deg C gives green in normal mode."""
        colour = tm_colour(65.0, "Traffic Light")
        assert colour is ft.Colors.GREEN_600

    def test_mid_tm_gives_amber(self) -> None:
        """50 <= Tm < 60 deg C gives amber in normal mode."""
        colour = tm_colour(55.0, "Traffic Light")
        assert colour is ft.Colors.ORANGE_600

    def test_low_tm_gives_red(self) -> None:
        """Tm < 50 deg C gives red in normal mode."""
        colour = tm_colour(45.0, "Traffic Light")
        assert colour is ft.Colors.RED_700

    def test_boundary_60(self) -> None:
        """Exactly 60 deg C is green (>= 60)."""
        assert tm_colour(60.0, "Traffic Light") is ft.Colors.GREEN_600

    def test_boundary_50(self) -> None:
        """Exactly 50 deg C is amber (>= 50, < 60)."""
        assert tm_colour(50.0, "Traffic Light") is ft.Colors.ORANGE_600

    def test_colour_deficient_high_tm_gives_blue(self) -> None:
        """Colour-deficient mode: Tm >= 60 deg C gives blue not green."""
        GUIColours.colour_deficient_mode = True
        colour = tm_colour(65.0, "Traffic Light")
        assert colour is ft.Colors.BLUE_600

    def test_colour_deficient_low_tm_gives_red(self) -> None:
        """Colour-deficient mode: Tm < 50 deg C still gives red."""
        GUIColours.colour_deficient_mode = True
        colour = tm_colour(40.0, "Traffic Light")
        assert colour is ft.Colors.RED_700

    def test_colour_deficient_mid_tm_gives_orange(self) -> None:
        """Colour-deficient mode: mid-range Tm gives orange."""
        GUIColours.colour_deficient_mode = True
        colour = tm_colour(55.0, "Traffic Light")
        assert colour is ft.Colors.ORANGE_600


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
