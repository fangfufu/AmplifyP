# conftest.py
from collections.abc import Generator

import pytest

from amplifyp.gui.colours import GUIColours


@pytest.fixture(autouse=True)  # type: ignore[untyped-decorator]
def reset_global_colour_modes() -> Generator[None, None, None]:
    """Reset colour modes between all tests to prevent leakage."""
    GUIColours.colour_deficient_mode = False
    GUIColours.dark_mode = False
    yield
    GUIColours.colour_deficient_mode = False
    GUIColours.dark_mode = False
