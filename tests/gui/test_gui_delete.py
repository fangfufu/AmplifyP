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

"""Test GUI Delete Primer functionality."""

import tkinter as tk
import unittest
from unittest.mock import MagicMock, patch

from amplifyp.gui import AmplifyPApp, Primer


class TestGUIDelete(unittest.TestCase):
    """Test Case for GUI Delete Primer."""

    root: tk.Tk

    @classmethod
    def setUpClass(cls) -> None:
        """Create root window once for the class."""
        cls.root = tk.Tk()

    @classmethod
    def tearDownClass(cls) -> None:
        """Destroy root window after all tests."""
        cls.root.destroy()

    def setUp(self) -> None:
        """Set up the test environment."""
        self.app = AmplifyPApp(self.root)
        # Prevent the mainloop from running if we ever call it (we won't)
        self.app.update()

    def tearDown(self) -> None:
        """Clean up after tests."""
        self.app.destroy()

    @patch("amplifyp.gui.messagebox.showwarning")
    def test_delete_primer(self, mock_showwarning: MagicMock) -> None:
        """Test deleting a primer."""
        # 1. Add Primers
        p1 = Primer("AAAA", "P1")
        p2 = Primer("TTTT", "P2")
        self.app.primers_data.append(p1)
        self.app.primers_data.append(p2)
        self.app.primers_list.insert(tk.END, f"{p1.name}: {p1.seq}")
        self.app.primers_list.insert(tk.END, f"{p2.name}: {p2.seq}")

        # 2. Select P1 (index 0)
        self.app.primers_list.selection_set(0)

        # 3. Delete
        self.app.delete_primer()

        # 4. Verify
        self.assertEqual(len(self.app.primers_data), 1)
        self.assertEqual(self.app.primers_data[0].name, "P2")
        # Listbox should have 1 item
        self.assertEqual(self.app.primers_list.size(), 1)
        self.assertEqual(self.app.primers_list.get(0), "P2: TTTT")

    @patch("amplifyp.gui.messagebox.showwarning")
    def test_delete_no_selection(self, mock_showwarning: MagicMock) -> None:
        """Test delete with no selection."""
        # 1. Add Primer
        p1 = Primer("AAAA", "P1")
        self.app.primers_data.append(p1)
        self.app.primers_list.insert(tk.END, "P1: AAAA")

        # 2. Ensure clear selection
        self.app.primers_list.selection_clear(0, tk.END)

        # 3. Try Delete
        self.app.delete_primer()

        # 4. Verify Warning
        mock_showwarning.assert_called_with(
            "Warning", "Please select a primer to delete."
        )

        # Verify nothing deleted
        self.assertEqual(len(self.app.primers_data), 1)
        self.assertEqual(self.app.primers_list.size(), 1)


if __name__ == "__main__":
    unittest.main()
