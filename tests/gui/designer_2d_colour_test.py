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

"""Unit tests for the designer_2d_colour() colour-mapping helper."""

from collections.abc import Generator

import pytest

from amplifyp.gui.colours import (
    GUIColours,
    designer_2d_colour,
    get_text_contrast_colour,
)


@pytest.fixture(autouse=True)  # type: ignore[untyped-decorator]
def reset_colour_modes() -> Generator[None, None, None]:
    """Reset colour-deficient and dark mode between tests."""
    GUIColours.colour_deficient_mode = False
    GUIColours.dark_mode = False
    yield
    GUIColours.colour_deficient_mode = False
    GUIColours.dark_mode = False


class TestDesigner2DColourNone:
    """Tests for the 'None' scheme."""

    def test_none_scheme_returns_none_for_any_value(self) -> None:
        """Scheme 'None' always returns None regardless of score."""
        assert designer_2d_colour(10.0, 0.0, 100.0, "None") is None
        assert designer_2d_colour(50.0, 0.0, 100.0, "None") is None

    def test_none_value_returns_none(self) -> None:
        """None value always returns None regardless of scheme."""
        assert designer_2d_colour(None, 0.0, 100.0, "Cool-Warm") is None
        assert designer_2d_colour(None, 0.0, 100.0, "Traffic Light") is None

    def test_unrecognised_scheme_returns_none(self) -> None:
        """Unrecognised scheme strings return None."""
        assert designer_2d_colour(50.0, 0.0, 100.0, "Rainbow") is None
        assert designer_2d_colour(50.0, 0.0, 100.0, "") is None


class TestDesigner2DColourCoolWarm:
    """Tests for the 'Cool-Warm' scheme."""

    def test_min_value_returns_red(self) -> None:
        """Minimum score (best quality) returns hotter red hex."""
        col = designer_2d_colour(0.0, 0.0, 100.0, "Cool-Warm")
        assert col == "#f9d0d0"  # GUIColours.GRID_2D_RED

    def test_max_value_returns_blue(self) -> None:
        """Maximum score (worst quality) returns cool blue hex."""
        col = designer_2d_colour(100.0, 0.0, 100.0, "Cool-Warm")
        assert col == "#d0e1f9"  # GUIColours.GRID_2D_BLUE

    def test_equal_min_max_returns_midpoint(self) -> None:
        """Equal min and max value defaults to factor t = 0.5."""
        col = designer_2d_colour(50.0, 50.0, 50.0, "Cool-Warm")
        assert col is not None
        assert isinstance(col, str)


class TestDesigner2DColourTrafficLight:
    """Tests for the 'Traffic Light' scheme."""

    def test_traffic_light_min_max_and_midpoint(self) -> None:
        """Traffic light returns green at min, yellow at mid, red at max."""
        col_min = designer_2d_colour(0.0, 0.0, 100.0, "Traffic Light")
        col_mid = designer_2d_colour(50.0, 0.0, 100.0, "Traffic Light")
        col_max = designer_2d_colour(100.0, 0.0, 100.0, "Traffic Light")

        assert col_min == "#d0f9d0"  # GRID_2D_GREEN
        assert col_mid == "#fdf6d0"  # GRID_2D_YELLOW
        assert col_max == "#f9d0d0"  # GRID_2D_RED

    def test_traffic_light_colour_deficient_mode(self) -> None:
        """Traffic light adapts green/yellow in colour-deficient mode."""
        GUIColours.colour_deficient_mode = True
        col_min = designer_2d_colour(0.0, 0.0, 100.0, "Traffic Light")
        col_mid = designer_2d_colour(50.0, 0.0, 100.0, "Traffic Light")

        assert col_min == "#d0e1f9"  # GRID_2D_BLUE
        assert col_mid == "#fdebd0"  # GRID_2D_ORANGE


class TestDesigner2DColourBlueOrange:
    """Tests for the 'Blue-Orange' scheme."""

    def test_blue_orange_min_and_max(self) -> None:
        """Blue-Orange returns orange at min (best) and blue at max (worst)."""
        col_min = designer_2d_colour(0.0, 0.0, 100.0, "Blue-Orange")
        col_max = designer_2d_colour(100.0, 0.0, 100.0, "Blue-Orange")

        assert col_min == "#fdebd0"
        assert col_max == "#d0e1f9"


class TestGetTextContrastColour:
    """Tests for the get_text_contrast_colour helper."""

    def test_light_background_returns_dark_text(self) -> None:
        """Light background returns dark text for high contrast."""
        assert get_text_contrast_colour("#d0e1f9") == "#111827"
        assert get_text_contrast_colour("#ffffff") == "#111827"

    def test_dark_background_returns_white_text(self) -> None:
        """Dark background returns white text for high contrast."""
        assert get_text_contrast_colour("#2d4a68") == "#ffffff"
        assert get_text_contrast_colour("#000000") == "#ffffff"
