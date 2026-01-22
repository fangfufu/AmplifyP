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

"""Test GUI Save/Load state functionality."""

import json
import os
import sys
import unittest
from typing import cast
from unittest.mock import MagicMock, patch

from tests.gui.mocks import create_mock_ctk, create_mock_tk


class TestGUIState(unittest.TestCase):
    """Test Case for GUI Save/Load state."""

    def setUp(self) -> None:
        """Set up the test environment."""
        self.sys_modules_backup = sys.modules.copy()

        mock_tk, mock_ttk = create_mock_tk()
        mock_ctk = create_mock_ctk()

        self.patcher = patch.dict(
            sys.modules,
            {
                "tkinter": mock_tk,
                "tkinter.ttk": mock_ttk,
                "customtkinter": mock_ctk,
            },
        )
        self.patcher.start()

        # Reload amplifyp.gui to use mocked modules
        if "amplifyp.gui" in sys.modules:
            del sys.modules["amplifyp.gui"]

        from amplifyp.gui import AmplifyPApp

        self.AmplifyPApp = AmplifyPApp

        # Instantiate app
        self.root = mock_ctk.CTk()
        self.app = self.AmplifyPApp(self.root)

        # We don't need update() because mocks don't loop

    def tearDown(self) -> None:
        """Clean up after tests."""
        self.patcher.stop()
        sys.modules.clear()
        sys.modules.update(self.sys_modules_backup)

    def test_save_state(self) -> None:
        """Test saving state to a JSON file."""
        # 1. Setup App State
        # template_text is MockCTkTextbox (MockCTkWidget)
        self.app.template_text.insert("1.0", "ATGC")

        from amplifyp.dna import Primer

        p1 = Primer("AAAA", "P1")
        self.app.primers_data.append(p1)

        self.app.primability_var.set(0.8)
        self.app.stability_var.set(0.2)

        # 2. Mock File Dialog
        # filedialog is imported in gui.py from tkinter.
        # Since we mocked tkinter, it is mock_tk.filedialog
        # But gui.py: 'from tkinter import filedialog' -> mock_tk.filedialog
        import tkinter
        import tkinter.filedialog
        import tkinter.messagebox

        mock_asksaveas = cast(MagicMock, tkinter.filedialog.asksaveasfilename)
        test_file = "test_save_state.json"
        mock_asksaveas.return_value = test_file

        # 3. Call Save
        try:
            self.app.save_state()

            # 4. Verify File Content
            with open(test_file, encoding="utf-8") as f:
                data = json.load(f)

            self.assertEqual(data["version"], 1)
            self.assertEqual(data["template"], "ATGC")
            self.assertEqual(len(data["primers"]), 1)
            self.assertEqual(data["primers"][0]["name"], "P1")
            self.assertEqual(data["primers"][0]["sequence"], "AAAA")
            self.assertEqual(data["settings"]["primability_cutoff"], 0.8)
            self.assertEqual(data["settings"]["stability_cutoff"], 0.2)

            # Check messagebox
            mock_showinfo = cast(MagicMock, tkinter.messagebox.showinfo)
            mock_showinfo.assert_called_with(
                "Success", "State saved successfully."
            )

        finally:
            if os.path.exists(test_file):
                os.remove(test_file)

    def test_load_state(self) -> None:
        """Test loading state from a JSON file."""
        # 1. Create Test File
        test_file = "test_load_state.json"
        state = {
            "version": 1,
            "template": "GGCC",
            "primers": [{"name": "P2", "sequence": "TTTT"}],
            "settings": {"primability_cutoff": 0.9, "stability_cutoff": 0.1},
        }
        with open(test_file, "w", encoding="utf-8") as f:
            json.dump(state, f)

        import tkinter
        import tkinter.filedialog
        import tkinter.messagebox

        mock_askopen = cast(MagicMock, tkinter.filedialog.askopenfilename)
        mock_askopen.return_value = test_file

        try:
            # 2. Call Load
            self.app.load_state()

            # 3. Verify App State
            # MockCTkTextbox.get returns self._text
            self.assertEqual(
                self.app.template_text.get("1.0", "end-1c"), "GGCC"
            )
            self.assertEqual(len(self.app.primers_data), 1)
            self.assertEqual(self.app.primers_data[0].name, "P2")
            self.assertEqual(self.app.primers_data[0].seq, "TTTT")

            # MockDoubleVar setup
            self.assertEqual(self.app.primability_var.get(), 0.9)
            self.assertEqual(self.app.stability_var.get(), 0.1)

            # Check Listbox contains the item
            # MockListbox.get(0)
            self.assertEqual(self.app.primers_list.get(0), "P2: TTTT")

            mock_showinfo = cast(MagicMock, tkinter.messagebox.showinfo)
            mock_showinfo.assert_called_with(
                "Success", "State loaded successfully."
            )

        finally:
            if os.path.exists(test_file):
                os.remove(test_file)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
