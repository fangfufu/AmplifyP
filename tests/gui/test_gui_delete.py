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

import sys
import unittest
from unittest.mock import patch

from tests.gui.mocks import create_mock_ctk, create_mock_tk


class TestGUIDelete(unittest.TestCase):
    """Test Case for GUI Delete Primer."""

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

        from amplifyp.gui import AmplifyPApp, Primer

        self.AmplifyPApp = AmplifyPApp
        self.Primer = Primer

        # Instantiate app
        self.root = mock_ctk.CTk()
        self.app = self.AmplifyPApp(self.root)

    def tearDown(self) -> None:
        """Clean up after tests."""
        self.patcher.stop()
        sys.modules.clear()
        sys.modules.update(self.sys_modules_backup)

    def test_delete_primer(self) -> None:
        """Test deleting a primer."""
        # 1. Add Primers
        p1 = self.Primer("AAAA", "P1")
        p2 = self.Primer("TTTT", "P2")
        self.app.primers_data.append(p1)
        self.app.primers_data.append(p2)

        # mock listbox insert
        self.app.primers_list.insert("end", f"{p1.name}: {p1.seq}")
        self.app.primers_list.insert("end", f"{p2.name}: {p2.seq}")

        # 2. Select P1 (index 0)
        # MockListbox needs helper for manual selection setup in test
        # self.app.primers_list.selection_set(0) -> adds to _selection
        self.app.primers_list.selection_set(0)

        # 3. Delete
        self.app.delete_primer()

        # 4. Verify
        self.assertEqual(len(self.app.primers_data), 1)
        self.assertEqual(self.app.primers_data[0].name, "P2")
        # Listbox should have 1 item
        self.assertEqual(self.app.primers_list.size(), 1)
        self.assertEqual(self.app.primers_list.get(0), f"{p2.name}: {p2.seq}")

    def test_delete_no_selection(self) -> None:
        """Test delete with no selection."""
        # 1. Add Primer
        p1 = self.Primer("AAAA", "P1")
        self.app.primers_data.append(p1)
        self.app.primers_list.insert("end", "P1: AAAA")

        # 2. Ensure clear selection
        # Default mock is empty anyway
        self.app.primers_list.selection_clear(0, "end")

        # 3. Try Delete
        self.app.delete_primer()

        # 4. Verify Warning
        import tkinter
        import tkinter.messagebox
        from typing import cast
        from unittest.mock import MagicMock

        mock_showwarning = cast(MagicMock, tkinter.messagebox.showwarning)
        mock_showwarning.assert_called_with(
            "Warning", "Please select a primer to delete."
        )

        # Verify nothing deleted
        self.assertEqual(len(self.app.primers_data), 1)
        self.assertEqual(self.app.primers_list.size(), 1)


if __name__ == "__main__":
    unittest.main()
