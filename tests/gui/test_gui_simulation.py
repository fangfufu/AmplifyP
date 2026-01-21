import unittest
from unittest.mock import MagicMock, patch
import sys
import os
import threading
import time

# Ensure src is in path
sys.path.insert(0, os.path.abspath("src"))

# Create mock classes for Tkinter
class MockTk:
    def __init__(self):
        pass
    def title(self, t): pass
    def geometry(self, g): pass
    def configure(self, **kwargs): pass
    def mainloop(self): pass

class MockFrame:
    def __init__(self, master=None):
        self.master = master
    def pack(self, **kwargs):
        pass
    def config(self, **kwargs):
        pass
    def after(self, ms, func=None, *args):
        # Verify after is called
        pass
    def quit(self):
        pass

class MockToplevel:
    def __init__(self, master=None):
        pass
    def title(self, text):
        pass
    def geometry(self, geom):
        pass
    def destroy(self):
        pass

# Prepare the mocks
mock_tk = MagicMock()
mock_tk.Tk = MockTk
mock_tk.Toplevel = MockToplevel
mock_tk.Misc = object # For type hinting
mock_tk.END = "end"

mock_ttk = MagicMock()
mock_ttk.Frame = MockFrame

# IMPORTANT: Link ttk to tk so 'from tkinter import ttk' works as expected
mock_tk.ttk = mock_ttk

with patch.dict(sys.modules, {'tkinter': mock_tk, 'tkinter.ttk': mock_ttk}):
    import tkinter as tk
    from tkinter import ttk
    from amplifyp.gui import AmplifyPApp
    from amplifyp.dna import Primer

class TestGUISimulation(unittest.TestCase):
    def setUp(self):
        # Instantiate app
        self.root = MockTk()
        self.app = AmplifyPApp(self.root)

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

    @patch('amplifyp.gui.messagebox')
    @patch('amplifyp.gui.AmpliconGenerator')
    @patch('amplifyp.gui.Repliconf')
    def test_simulate_pcr_starts_thread(self, MockRepliconf, MockAmpliconGenerator, MockMessagebox):
        # Setup mocks
        mock_gen = MockAmpliconGenerator.return_value
        mock_gen.get_amplicons.return_value = []

        # Run simulate_pcr
        self.app.simulate_pcr()

        # Check if button disabled
        self.app.simulate_btn.config.assert_called_with(state="disabled")
        self.app.config.assert_called_with(cursor="watch")

        # Check if thread started
        self.assertTrue(hasattr(self.app, 'simulation_thread'))
        self.assertIsInstance(self.app.simulation_thread, threading.Thread)
        # We don't check is_alive() because thread might finish instantly due to mocking/small input

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
        self.app.simulate_btn.config.assert_called_with(state="normal")
        self.app.config.assert_called_with(cursor="")

        # Check results are populated (meaning thread ran)
        # If AmpliconGenerator was real, it produced results. If mock, it produced return_value.
        # Since we found out it's running real one, it should have results.
        self.assertIsNotNone(self.app.simulation_result)

if __name__ == '__main__':
    unittest.main()
