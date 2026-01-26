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

"""Tests for GUI simulation functionality."""

import os
import sys
import threading
import unittest
from unittest.mock import MagicMock, patch

from tests.gui.mocks import create_mock_ctk

# Ensure src is in path
sys.path.insert(0, os.path.abspath("src"))


# Create mock classes for Tkinter
class MockTk:
    """Mock Tk class."""

    def __init__(self) -> None:
        """Initialize MockTk."""

    def title(self, t: str) -> None:
        """Set title."""

    def geometry(self, g: str) -> None:
        """Set geometry."""

    def configure(self, **kwargs: object) -> None:
        """Configure widget."""

    def mainloop(self) -> None:
        """Start mainloop."""


class MockFrame:
    """Mock Frame class."""

    def __init__(self, master: object = None) -> None:
        """Initialize MockFrame."""
        self.master = master

    def pack(self, **kwargs: object) -> None:
        """Pack widget."""

    def config(self, **kwargs: object) -> None:
        """Configure widget."""

    def after(self, ms: int, func: object = None, *args: object) -> None:
        """Schedule callback."""
        # Verify after is called

    def quit(self) -> None:
        """Quit application."""


class MockToplevel:
    """Mock Toplevel class."""

    def __init__(self, master: object = None) -> None:
        """Initialize MockToplevel."""

    def title(self, text: str) -> None:
        """Set title."""

    def geometry(self, geom: str) -> None:
        """Set geometry."""

    def destroy(self) -> None:
        """Destroy widget."""


# Prepare the mocks
mock_tk = MagicMock()
mock_tk.Tk = MockTk
mock_tk.Toplevel = MockToplevel
mock_tk.Misc = object  # For type hinting
mock_tk.END = "end"

mock_ttk = MagicMock()
mock_ttk.Frame = MockFrame

# IMPORTANT: Link ttk to tk so 'from tkinter import ttk' works as expected
mock_tk.ttk = mock_ttk


class TestGUISimulation(unittest.TestCase):
    """Test GUI Simulation."""

    def setUp(self) -> None:
        """Set up test fixtures."""
        from amplifyp.dna import Primer

        self.sys_modules_backup = sys.modules.copy()

        # Patch tkinter modules
        self.patcher = patch.dict(
            sys.modules,
            {
                "tkinter": mock_tk,
                "tkinter.ttk": mock_ttk,
                "customtkinter": create_mock_ctk(),
            },
        )
        self.patcher.start()

        # Reload amplifyp.gui to use mocked tkinter
        # Reload amplifyp.gui to use mocked tkinter
        for mod in ["amplifyp.gui", "amplifyp.gui.gui"]:
            if mod in sys.modules:
                del sys.modules[mod]

        # We need to import inside the patched environment
        from amplifyp.gui import AmplifyPApp

        self.AmplifyPApp = AmplifyPApp

        # Instantiate app
        # mock_ctk is the MagicMock returned by create_mock_ctk()
        # The classes are attributes of it.
        # But wait, self.AmplifyPApp was imported AFTER patch, so it uses the
        # mocks in sys.modules['customtkinter'] which is the object provided by
        # create_mock_ctk().
        # So create_mock_ctk().CTk should be the class.
        # But I need to instantiate it.
        # simpler: use the one from sys.modules
        mock_ctk = sys.modules["customtkinter"]
        self.root = mock_ctk.CTk()
        self.app = self.AmplifyPApp(self.root)

        # Setup some mock widgets that simulate_pcr uses
        self.app.template_text = MagicMock()
        self.app.template_text.get.return_value = "GATC" * 10

        self.app.primability_var = MagicMock()
        self.app.primability_var.get.return_value = 0.5

        self.app.stability_var = MagicMock()
        self.app.stability_var.get.return_value = 0.5

        self.app.simulate_btn = MagicMock()

        self.app.tree = MagicMock()
        self.app.tree.get_children.return_value = []

        self.app.primers_data = [Primer("GATC", "P1")]

        # Mock self.after to capture callbacks
        self.app.after = MagicMock()
        self.app.config = MagicMock()

    def tearDown(self) -> None:
        """Clean up."""
        self.patcher.stop()

        # Restore sys.modules to original state
        sys.modules.clear()
        sys.modules.update(self.sys_modules_backup)

    @patch("amplifyp.gui.gui.messagebox")
    @patch("amplifyp.gui.gui.AmpliconGenerator")
    @patch("amplifyp.gui.gui.Repliconf")
    def test_simulate_pcr_starts_thread(
        self,
        MockRepliconf: MagicMock,
        MockAmpliconGenerator: MagicMock,
        MockMessagebox: MagicMock,
    ) -> None:
        """Test that simulating PCR starts a new thread."""
        # Setup mocks
        mock_gen = MockAmpliconGenerator.return_value
        mock_gen.get_amplicons.return_value = []

        # Run simulate_pcr
        self.app.simulate_pcr()

        # Check if button disabled
        self.app.simulate_btn.configure.assert_called_with(state="disabled")
        self.app.configure.assert_called_with(cursor="watch")

        # Check if thread started
        self.assertTrue(hasattr(self.app, "simulation_thread"))
        self.assertIsInstance(self.app.simulation_thread, threading.Thread)
        # We don't check is_alive() because thread might finish instantly due to
        # mocking or small input
        # or small input

        # Check if after called to schedule check
        self.app.after.assert_called()

        # Wait for thread to finish
        self.app.simulation_thread.join(timeout=1)

        # Manually trigger check_simulation completion logic
        # We assume the thread finished and populated results (or error)

        # Trigger completion
        self.app.on_simulation_complete()

        # Check if on_simulation_complete logic ran
        # Button enabled
        self.app.simulate_btn.configure.assert_called_with(state="normal")
        self.app.configure.assert_called_with(cursor="")

        # Check results are populated (meaning thread ran)
        # If AmpliconGenerator was real, it produced results. If mock,
        # it produced return_value.
        # Since we found out it's running real one, it should have results.
        self.assertIsNotNone(self.app.simulation_result)


if __name__ == "__main__":
    unittest.main()
