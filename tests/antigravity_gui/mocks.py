from unittest.mock import MagicMock


class MockCTk:
    """Mock CTk class."""

    def __init__(self) -> None:
        """Initialize MockCTk."""

    def mainloop(self) -> None:
        """Start mainloop."""

    def set_appearance_mode(self, mode: str) -> None:
        """Set appearance mode."""

    def set_default_color_theme(self, theme: str) -> None:
        """Set default color theme."""

    def title(self, text: str) -> None:
        """Set title."""

    def geometry(self, geom: str) -> None:
        """Set geometry."""

    def configure(self, **kwargs: object) -> None:
        """Configure."""


class MockCTkFrame:
    """Mock CTkFrame class."""

    def __init__(self, master: object = None, **kwargs: object) -> None:
        """Initialize MockCTkFrame."""
        self.master = master
        self.pack = MagicMock()
        self.config = MagicMock()
        self.configure = MagicMock()
        self.after = MagicMock()
        self.destroy = MagicMock()
        self.grid_columnconfigure = MagicMock()
        self.grid_rowconfigure = MagicMock()

    def quit(self) -> None:
        """Quit application."""


class MockCTkToplevel:
    """Mock CTkToplevel class."""

    def __init__(self, master: object = None) -> None:
        """Initialize MockCTkToplevel."""

    def title(self, text: str) -> None:
        """Set title."""

    def geometry(self, geom: str) -> None:
        """Set geometry."""

    def destroy(self) -> None:
        """Destroy widget."""


class MockCTkWidget:
    """Mock generic CTk widget."""

    def __init__(self, master: object = None, **kwargs: object) -> None:
        """Initialize."""
        self.master = master
        self._text = ""
        self.pack = MagicMock()
        self.grid = MagicMock()
        self.configure = MagicMock()
        self.bind = MagicMock()
        self.event_generate = MagicMock()

    def get(self, *args: object) -> object:
        """Get value."""
        # Handle arguments for Textbox (index1, index2)
        if len(args) >= 2:
            return self._text
        return self._text

    def insert(self, index: object, text: str, tags: object = None) -> None:
        """Insert."""
        # Very simple append logic for verification
        self._text += text

    def delete(self, *args: object) -> None:
        """Delete."""
        self._text = ""


# --- Tkinter Mocks ---


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

    def update(self) -> None:
        """Update."""

    def destroy(self) -> None:
        """Destroy."""


class MockFrame:
    """Mock Frame class."""

    def __init__(self, master: object = None) -> None:
        """Initialize MockFrame."""
        self.master = master

    def pack(self, **kwargs: object) -> None:
        """Pack widget."""

    def config(self, **kwargs: object) -> None:
        """Configure widget."""


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


class MockStringVar:
    """Mock StringVar."""

    def __init__(self, value: str = "") -> None:
        self._value = value

    def get(self) -> str:
        return self._value

    def set(self, value: str) -> None:
        self._value = value


class MockDoubleVar:
    """Mock DoubleVar."""

    def __init__(self, value: float = 0.0) -> None:
        self._value = value

    def get(self) -> float:
        return self._value

    def set(self, value: float) -> None:
        self._value = value


class MockListbox:
    """Mock Listbox."""

    def __init__(self, master: object = None, **kwargs: object) -> None:
        self._items: list[str] = []
        self._selection: set[int] = set()

    def pack(self, **kwargs: object) -> None:
        pass

    def insert(self, index: object, *elements: str) -> None:
        # Simplification: always append for "end"
        for e in elements:
            self._items.append(e)

    def delete(self, first: int, last: object = None) -> None:
        if last is None:
            if 0 <= first < len(self._items):
                self._items.pop(first)
        elif last == "end" or last == "END":  # tk.END is string "end"
            del self._items[first:]
        else:
            # integer?
            pass

    def size(self) -> int:
        return len(self._items)

    def get(self, first: int, last: object = None) -> str | tuple[str, ...]:
        if last is None:
            return self._items[first]
        return tuple(self._items[first : len(self._items)])  # Simplified logic

    def curselection(self) -> tuple[int, ...]:
        return tuple(sorted(self._selection))

    def selection_set(self, first: int, last: object = None) -> None:
        self._selection.add(first)

    def selection_clear(self, first: int, last: object = None) -> None:
        self._selection.clear()

    def nearest(self, y: int) -> int:
        return 0

    def activate(self, index: int) -> None:
        pass

    def bind(
        self,
        sequence: str | None = None,
        func: object = None,
        add: bool | None = None,
    ) -> None:
        pass


class MockMenu:
    """Mock Menu."""

    def __init__(self, master: object = None, tearoff: int = 0) -> None:
        pass

    def add_command(self, **kwargs: object) -> None:
        pass

    def add_separator(self) -> None:
        pass

    def add_cascade(self, **kwargs: object) -> None:
        pass

    def tk_popup(self, x: int, y: int) -> None:
        pass


# Helper to create the mock dicts
def create_mock_tk() -> tuple[MagicMock, MagicMock]:
    mock_tk = MagicMock()
    mock_tk.Tk = MockTk
    mock_tk.Toplevel = MockToplevel
    mock_tk.Frame = MockFrame
    mock_tk.StringVar = MockStringVar
    mock_tk.DoubleVar = MockDoubleVar
    mock_tk.Listbox = MockListbox
    mock_tk.Menu = MockMenu
    mock_tk.END = "end"
    mock_tk.Event = MagicMock

    mock_ttk = MagicMock()
    mock_ttk.Frame = MockFrame
    mock_ttk.LabelFrame = MockFrame
    mock_ttk.Scrollbar = MagicMock()
    mock_ttk.Treeview = MagicMock()
    mock_ttk.Label = MagicMock()
    mock_ttk.Entry = MagicMock()
    mock_ttk.Button = MagicMock()

    mock_tk.filedialog = MagicMock()
    mock_tk.messagebox = MagicMock()

    # Link ttk
    mock_tk.ttk = mock_ttk

    return mock_tk, mock_ttk


def create_mock_ctk() -> MagicMock:
    mock_ctk = MagicMock()
    mock_ctk.CTk = MockCTk
    mock_ctk.CTkFrame = MockCTkFrame
    mock_ctk.CTkToplevel = MockCTkToplevel
    mock_ctk.CTkButton = MockCTkWidget
    mock_ctk.CTkEntry = MockCTkWidget
    mock_ctk.CTkLabel = MockCTkWidget
    mock_ctk.CTkTextbox = MockCTkWidget
    mock_ctk.set_appearance_mode = MagicMock()
    mock_ctk.set_default_color_theme = MagicMock()
    return mock_ctk
