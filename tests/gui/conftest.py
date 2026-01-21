"""Pytest configuration and fixtures for GUI tests."""

import copy
from collections.abc import Generator

import pytest

from amplifyp import gui
from amplifyp.settings import DEFAULT_SETTINGS


@pytest.fixture(autouse=True)  # type: ignore[untyped-decorator]
def reset_gui_settings() -> Generator[None, None, None]:
    """Fixture to reset gui.settings to defaults after each test.

    This ensures that any modifications to the global gui.settings object
    during a test do not leak into other tests.
    """
    yield
    gui.settings = copy.deepcopy(DEFAULT_SETTINGS)
