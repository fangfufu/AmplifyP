# conftest.py
from collections.abc import Generator
from typing import Any

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


@pytest.fixture(scope="session")  # type: ignore[untyped-decorator]
def browser_context_args(
    browser_context_args: dict[str, Any],
) -> dict[str, Any]:
    """Block service workers to prevent caching/stale code issues.

    This ensures E2E tests always use the latest build.
    """
    return {
        **browser_context_args,
        "service_workers": "block",
    }
