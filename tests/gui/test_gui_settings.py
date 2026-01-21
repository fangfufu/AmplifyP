"""Test GUI settings using Amplify 4's settings module."""

from amplifyp import gui
from amplifyp.settings import DEFAULT_SETTINGS, Settings


def test_gui_settings_exists() -> None:
    """Verify that gui.settings exists and is of the correct type."""
    assert hasattr(gui, "settings"), "gui module should have a 'settings' attribute"
    assert isinstance(gui.settings, Settings), "gui.settings should be an instance of Settings"


def test_gui_settings_initialization() -> None:
    """Verify that gui.settings is initialized with correct values."""
    # Direct comparison fails because helper classes don't implement __eq__ and deepcopy
    # creates new instances, so we check individual fields or values.
    # In the original script, we checked primability_cutoff and stability_cutoff.
    assert gui.settings.primability_cutoff == DEFAULT_SETTINGS.primability_cutoff
    assert gui.settings.stability_cutoff == DEFAULT_SETTINGS.stability_cutoff


def test_gui_settings_modification() -> None:
    """Verify that gui.settings can be modified and persists."""
    gui.settings.primability_cutoff = 0.5
    assert gui.settings.primability_cutoff == 0.5
