"""Test GUI settings using Amplify 4's settings module."""

from amplifyp import gui
from amplifyp.settings import DEFAULT_SETTINGS, Settings


def test_gui_settings_exists() -> None:
    """Verify that gui.settings exists and is of the correct type."""
    assert hasattr(gui, "settings"), "gui module should have a 'settings' attribute"
    assert isinstance(gui.settings, Settings), (
        "gui.settings should be an instance of Settings"
    )


def test_gui_settings_initialization() -> None:
    """Verify that gui.settings is initialized with correct values."""
    # Direct comparison fails because helper classes don't implement __eq__ and deepcopy
    # creates new instances, so we check individual fields or values.
    # In the original script, we checked primability_cutoff and stability_cutoff.
    assert gui.settings.primability_cutoff == DEFAULT_SETTINGS.primability_cutoff
    assert gui.settings.stability_cutoff == DEFAULT_SETTINGS.stability_cutoff


def test_gui_settings_modification() -> None:
    """Verify that gui.settings can be modified and persists."""
    # We should probably reset settings after this test or use a fixture if we
    # wanted full isolation, but since this runs in a separate process or we want
    # to test state persistence (as implied by "persist"), we'll follow the script's
    # logic. However, purely modifying global state in tests is risky. For now,
    # I will mirror the script's logic to pass the "verify it persists" requirement
    # within the scope of the test execution.

    # Save original value to restore it later if needed (good practice)
    original_cutoff = gui.settings.primability_cutoff

    try:
        gui.settings.primability_cutoff = 0.5
        assert gui.settings.primability_cutoff == 0.5
    finally:
        # Restore state to avoid polluting other tests if they run in the same process
        gui.settings.primability_cutoff = original_cutoff
