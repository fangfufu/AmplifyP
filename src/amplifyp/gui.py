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

"""GUI application module for AmplifyP."""

__all__ = ["AmplifyPApp", "Primer", "PrimerStatsDialog", "settings"]

import copy
import json
import os
import sys
import threading
import tkinter as tk
import traceback
from abc import ABC, abstractmethod
from tkinter import filedialog, messagebox, ttk
from typing import Any

import customtkinter as ctk

from amplifyp.amplicon import Amplicon, AmpliconGenerator
from amplifyp.dna import (
    DNA,
    DNAType,
    Primer,
)
from amplifyp.repliconf import Repliconf
from amplifyp.settings import (
    GLOBAL_REPLICATION_SETTINGS,
    ReplicationSettings,
)

# Set theme
ctk.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
ctk.set_default_color_theme(
    "blue"
)  # Themes: "blue" (standard), "green", "dark-blue"

# Ensure we can import amplifyp when running as a script
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(
        0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
    )

settings: ReplicationSettings = copy.deepcopy(GLOBAL_REPLICATION_SETTINGS)


class PrimerStatsDialog(ctk.CTkToplevel):  # type: ignore[misc]
    """A dialog window for displaying analysis statistics for a single primer.

    This window calculates and lists all potential binding sites for a given
    primer on the template DNA, showing detailed scores for primability,
    stability, and quality.

    Attributes:
        primer (Primer): The primer being analyzed.
        template_seq (str): The template DNA sequence string.
    """

    def __init__(
        self, parent: tk.Misc, primer: Primer, template_seq: str
    ) -> None:
        """Initialize the PrimerStatsDialog.

        Args:
            parent (tk.Misc): The parent widget.
            primer (Primer): The primer object to analyze.
            template_seq (str): The raw template DNA sequence.
        """
        super().__init__(parent)
        self.title(f"Primer Properties: {primer.name}")
        self.geometry("800x400")

        self.primer = primer
        self.template_seq = template_seq

        self.create_widgets()
        self.analyze()

    def create_widgets(self) -> None:
        """Initialize and layout the GUI widgets (Treeview, Scrollbar)."""
        columns = ("pos", "strand", "primability", "stability", "quality")
        self.tree = ttk.Treeview(self, columns=columns, show="headings")
        self.tree.heading("pos", text="Position")
        self.tree.heading("strand", text="Strand")
        self.tree.heading("primability", text="Primability")
        self.tree.heading("stability", text="Stability")
        self.tree.heading("quality", text="Quality")

        # Add scrollbar
        # Note: CTk doesn't have a native Treeview, so we keep using
        # ttk.Treeview but we can try to style it or just let it be.
        scrollbar = ttk.Scrollbar(
            self, orient="vertical", command=self.tree.yview
        )
        self.tree.configure(yscrollcommand=scrollbar.set)

        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def analyze(self) -> None:
        """Perform the analysis and populate the Treeview with results."""
        try:
            template = DNA(self.template_seq, DNAType.LINEAR, "Template")
            rc = Repliconf(template, self.primer, settings)
            rc.search()  # This populates origin_idx

            # Helper to add rows
            def add_rows(indices: list[Any]) -> None:
                for i in indices:
                    origin = rc.origin(i)
                    self.tree.insert(
                        "",
                        tk.END,
                        values=(
                            i.index,
                            i.direction,
                            f"{origin.primability:.4f}",
                            f"{origin.stability:.4f}",
                            f"{origin.quality:.4f}",
                        ),
                    )

            add_rows(rc.origin_db.fwd)
            add_rows(rc.origin_db.rev)

        except Exception as e:
            traceback.print_exc()
            messagebox.showerror(
                "Analysis Error",
                f"Failed to analyze primer: {e}\nSee console for details.",
            )
            self.destroy()


# --- Plotter Classes (Moved from gui2/plotter.py) ---
class PlotterThing(ABC):
    """Abstract base class for things that can be plotted."""

    def __init__(self) -> None:
        """Initialize the PlotterThing."""
        self.bounds = (0, 0, 0, 0)  # x, y, width, height

    @abstractmethod
    def plot(self, canvas: tk.Canvas) -> None:
        """Draw the object on the canvas.

        Args:
            canvas (tk.Canvas): The canvas to draw on.
        """


class StringThing(PlotterThing):
    """Draws a text string."""

    def __init__(
        self, text: str, x: float, y: float, anchor: str = "nw", **kwargs: Any
    ) -> None:
        """Initialize a StringThing.

        Args:
            text (str): The text to draw.
            x (float): X coordinate.
            y (float): Y coordinate.
            anchor (str, optional): The anchor point. Defaults to "nw".
            kwargs (Any): Additional arguments for canvas.create_text.
        """
        super().__init__()
        self.text = text
        self.x = x
        self.y = y
        self.anchor = anchor
        self.kwargs = kwargs  # e.g. font, fill

    def plot(self, canvas: tk.Canvas) -> None:
        """Draw the text string."""
        canvas.create_text(  # type: ignore[call-overload]
            self.x, self.y, text=self.text, anchor=self.anchor, **self.kwargs
        )


class VStringThing(PlotterThing):
    """Draws a vertical text string (rotated 90 degrees counter-clockwise)."""

    def __init__(
        self, text: str, x: float, y: float, anchor: str = "ne", **kwargs: Any
    ) -> None:
        """Initialize a VStringThing.

        Args:
            text (str): The text to draw.
            x (float): X coordinate.
            y (float): Y coordinate.
            anchor (str, optional): The anchor point. Defaults to "ne".
            kwargs (Any): Additional arguments.
        """
        super().__init__()
        self.text = text
        self.x = x
        self.y = y
        self.anchor = anchor
        self.kwargs = kwargs

    def plot(self, canvas: tk.Canvas) -> None:
        """Draw the vertical text string."""
        # Tkinter canvas doesn't support text rotation natively easily.
        # We can use the 'angle' parameter if supported (Tk 8.6+).
        canvas.create_text(  # type: ignore[call-overload]
            self.x,
            self.y,
            text=self.text,
            angle=90,
            anchor=self.anchor,
            **self.kwargs,
        )


class VStringRectThing(PlotterThing):
    """Draws a vertical text string within a rectangle (conceptually)."""

    def __init__(
        self, text: str, x: float, y: float, height: float, **kwargs: Any
    ) -> None:
        """Initialize a VStringRectThing.

        Args:
            text (str): The text to draw.
            x (float): X coordinate.
            y (float): Y coordinate.
            height (float): height of the rectangle (unused in simplified plot).
            kwargs (Any): Additional arguments.
        """
        super().__init__()
        self.text = text
        self.x = x
        self.y = y
        self.height = height
        self.kwargs = kwargs

    def plot(self, canvas: tk.Canvas) -> None:
        """Draw the rotated text."""
        canvas.create_text(
            self.x,
            self.y,
            text=self.text,
            angle=90,
            # Anchor at bottom-left (which corresponds to top-left after
            # rotation) - check this
            anchor="sw",
            **self.kwargs,
        )


class LineThing(PlotterThing):
    """Draws a line or path."""

    def __init__(self, coords: list[float], **kwargs: Any) -> None:
        """Initialize a LineThing.

        Args:
            coords (list[float]): List of coordinates [x1, y1, x2, y2, ...].
            kwargs (Any): Additional arguments for canvas.create_line.
        """
        super().__init__()
        self.coords = coords
        self.kwargs = kwargs  # width, fill, etc.

    def plot(self, canvas: tk.Canvas) -> None:
        """Draw the line."""
        canvas.create_line(*self.coords, **self.kwargs)


class RectThing(PlotterThing):
    """Draws a rectangle."""

    def __init__(
        self, x1: float, y1: float, x2: float, y2: float, **kwargs: Any
    ) -> None:
        """Initialize a RectThing.

        Args:
            x1 (float): Left x.
            y1 (float): Top y.
            x2 (float): Right x.
            y2 (float): Bottom y.
            kwargs (Any): Additional arguments for canvas.create_rectangle.
        """
        super().__init__()
        self.coords = (x1, y1, x2, y2)
        self.kwargs = kwargs  # outline, fill, width, tags

    def plot(self, canvas: tk.Canvas) -> None:
        """Draw the rectangle."""
        canvas.create_rectangle(*self.coords, **self.kwargs)


class ArrowThing(PlotterThing):
    """Draws a simplistic arrow head."""

    def __init__(
        self, points: list[float], fill: str = "black", **kwargs: Any
    ) -> None:
        """Initialize an ArrowThing.

        Args:
            points (list[float]): List of polygon points.
            fill (str, optional): Fill color. Defaults to "black".
            kwargs (Any): Additional arguments for canvas.create_polygon.
        """
        super().__init__()
        self.points = points
        self.fill = fill
        self.kwargs = kwargs

    def plot(self, canvas: tk.Canvas) -> None:
        """Draw the arrow (polygon)."""
        canvas.create_polygon(
            *self.points, fill=self.fill, outline="", **self.kwargs
        )


class GroupThing(PlotterThing):
    """A group of plotters."""

    def __init__(self, items: list[PlotterThing]) -> None:
        """Initialize a GroupThing.

        Args:
            items (list[PlotterThing]): List of plotter items.
        """
        super().__init__()
        self.items = items

    def plot(self, canvas: tk.Canvas) -> None:
        """Draw all items in the group."""
        for item in self.items:
            item.plot(canvas)


# --- MapView Class (Moved from gui2/map_view.py) ---
class MapView(tk.Canvas):
    """A Canvas widget for visualizing DNA amplicons."""

    def __init__(self, master: Any, **kwargs: Any) -> None:
        """Initialize the MapView.

        Args:
            master (Any): The parent widget.
            kwargs (Any): Additional keyword arguments for the Canvas.
        """
        super().__init__(master, bg="white", **kwargs)
        self.plotters: list[PlotterThing] = []
        self.target_length = 0
        self.width = 0.0

        # Drawing constants
        self.h_margin = 15.0  # Horizontal margin
        self.v_margin_top = 5.0
        self.v_tick_top = 35.0  # Top of ticks
        self.v_baseline = 105.0  # Target baseline Y position
        self.v_arrow_offset = 17.0  # Arrow offset from baseline
        self.v_primer_name = 37.0  # Primer name offset
        self.v_frag_start = 130.0  # Start of fragment bars
        self.v_frag_spacing = 30.0  # Spacing between fragments
        self.tick_height_up = 5.0
        self.tick_height_down = 5.0

        self.bind("<Configure>", self.on_resize)

    def on_resize(self, event: tk.Event) -> None:
        """Handle canvas resize events."""
        self.width = float(event.width)
        # Re-draw if we have data
        if self.plotters:
            self.draw()

    def set_data(
        self, template: DNA, amplicons: list[Amplicon], primers: list[Any]
    ) -> None:
        """Set the data to be visualized and prepare plotters."""
        self.template = template
        self.amplicons = amplicons
        self.primers = primers  # These are "used primers" conceptually for now
        self.target_length = len(template)
        self.prepare_plotters()
        self.draw()

    def prepare_plotters(self) -> None:
        """Create the PlotterThing objects based on the data."""
        self.plotters = []
        if self.target_length == 0:
            return

        # We need self.width to calculate positions, but it might not be set yet
        # if not mapped.
        # on_resize will trigger redraw.
        if self.width <= 1.0:
            self.width = float(self.winfo_width())
            if self.width <= 1.0:
                self.width = 800.0  # Default fallback

        twidth = self.width - 2 * self.h_margin
        points_per_base = twidth / self.target_length

        def basex(base_idx: int) -> float:
            return self.h_margin + base_idx * points_per_base

        # 1. Base Numbers
        self.plotters.append(
            StringThing(
                "1", self.h_margin, self.v_margin_top, font=("Arial", 10)
            )
        )
        self.plotters.append(
            StringThing(
                str(self.target_length),
                self.width - self.h_margin,
                self.v_margin_top,
                anchor="ne",
                font=("Arial", 10),
            )
        )

        # 2. Baseline
        self.plotters.append(
            LineThing(
                [
                    self.h_margin,
                    self.v_baseline,
                    self.width - self.h_margin,
                    self.v_baseline,
                ],
                width=2,
                fill="black",
            )
        )

        # 3. Ticks
        for tick_base in range(100, self.target_length, 100):
            x = basex(tick_base)
            is_large = tick_base % 1000 == 0
            factor = 1.5 if is_large else 1.0

            # Draw tick
            self.plotters.append(
                LineThing(
                    [
                        x,
                        self.v_baseline - (self.tick_height_up * factor),
                        x,
                        self.v_baseline + (self.tick_height_down * factor),
                    ],
                    width=1,
                    fill="black",
                )
            )

        # 4. Amplicons / Fragments
        current_y = self.v_frag_start
        for amp in self.amplicons:
            # Linear logic mostly for now
            start_x = basex(amp.start.index)
            end_x = basex(amp.end.index)

            # If circular and wraps, visual representation is tricky.
            # Simplified: just draw from start to end if linear-ish,
            # or two parts.
            # Using simple bar for now.

            if amp.circular:
                # TODO: visual for circular wrap
                pass

            width = end_x - start_x
            if width < 0:  # Wrap around or circular
                # Just draw nothing for now or handle later
                pass

            # Bar
            self.plotters.append(
                RectThing(
                    start_x,
                    current_y,
                    end_x,
                    current_y + 10,
                    fill="black",
                    outline="",
                )
            )

            # Size label centered
            center_x = start_x + width / 2
            self.plotters.append(
                StringThing(
                    str(len(amp.product)) + " bp",
                    center_x,
                    current_y + 12,
                    anchor="n",
                    font=("Arial", 9),
                )
            )

            current_y += self.v_frag_spacing

        # Update scroll region
        total_height = current_y + 50
        self.configure(scrollregion=(0, 0, self.width, total_height))

    def draw(self) -> None:
        """Clear canvas and execute plotters."""
        self.delete("all")
        # Re-calculate positions if width changed
        self.prepare_plotters()

        for p in self.plotters:
            p.plot(self)


class AmplifyPApp(ctk.CTkFrame):  # type: ignore[misc]
    """The main application frame for the AmplifyP GUI.

    This class manages the main window layout, including inputs for template DNA
    and primers, settings configuration, and the display of simulation results.
    """

    def __init__(self, master: tk.Tk | None = None) -> None:
        """Initialize the AmplifyP application.

        Args:
            master (tk.Tk, optional): The root window. Defaults to None.
        """
        super().__init__(master)
        if master:
            master.title("AmplifyP")
            master.geometry("1000x800")
            # Create a menu
            self.create_menu()

        self.pack(fill="both", expand=True)
        self.create_widgets()

    def create_menu(self) -> None:
        """Create and attach the main application menu bar."""
        # Menus are still attached to the root tk/ctk object
        if not isinstance(self.master, (tk.Tk, ctk.CTk)):
            return

        menubar = tk.Menu(self.master)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Save State", command=self.save_state)
        file_menu.add_command(label="Load State", command=self.load_state)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.quit)

        menubar.add_cascade(label="File", menu=file_menu)
        self.master.configure(menu=menubar)

    def create_context_menu(
        self, widget: tk.Text | ttk.Entry | ctk.CTkEntry
    ) -> None:
        """Attach a standard Cut/Copy/Paste context menu to a widget.

        Args:
            widget (tk.Text | ttk.Entry | ctk.CTkEntry): The widget to attach
                the menu to.
        """
        menu = tk.Menu(widget, tearoff=0)
        menu.add_command(
            label="Cut", command=lambda: widget.event_generate("<<Cut>>")
        )
        menu.add_command(
            label="Copy", command=lambda: widget.event_generate("<<Copy>>")
        )
        menu.add_command(
            label="Paste", command=lambda: widget.event_generate("<<Paste>>")
        )

        def show_menu(event: tk.Event) -> None:
            menu.tk_popup(event.x_root, event.y_root)

        widget.bind("<Button-3>", show_menu)

    def _create_labeled_entry(
        self,
        parent: ctk.CTkFrame,
        label: str,
        variable: tk.Variable,
        width: int = 140,
        **pack_kwargs: Any,
    ) -> ctk.CTkEntry:
        """Create and pack a labeled entry widget with context menu.

        Args:
            parent: The parent widget.
            label: Text for the label.
            variable: The variable to bind.
            width: Width of the entry.
            pack_kwargs: Additional arguments for packing the entry widget (e.g.
                fill, expand).

        Returns:
            The created and packed entry widget.
        """
        ctk.CTkLabel(parent, text=label).pack(side="left", padx=5)
        entry = ctk.CTkEntry(parent, textvariable=variable, width=width)

        pack_args = {"side": "left", "padx": 5}
        pack_args.update(pack_kwargs)

        entry.pack(**pack_args)
        self.create_context_menu(entry)
        return entry

    def create_widgets(self) -> None:
        """Initialize and arrange all widgets in the main frame."""
        # --- Template DNA Section ---
        template_frame = ctk.CTkFrame(self)
        template_frame.pack(fill="x", padx=10, pady=5)

        # Label for the frame since CTkFrame doesn't have 'text' argument like
        # LabelFrame
        ctk.CTkLabel(
            template_frame, text="Template DNA", font=("", 14, "bold")
        ).pack(anchor="w", padx=10, pady=(5, 0))

        self.template_text = ctk.CTkTextbox(template_frame, height=100)
        self.template_text.pack(fill="x", padx=5, pady=5)
        self.create_context_menu(self.template_text)

        # --- Primers Section ---
        primers_frame = ctk.CTkFrame(self)
        primers_frame.pack(fill="both", expand=True, padx=10, pady=5)

        ctk.CTkLabel(primers_frame, text="Primers", font=("", 14, "bold")).pack(
            anchor="w", padx=10, pady=(5, 0)
        )

        input_frame = ctk.CTkFrame(primers_frame, fg_color="transparent")
        input_frame.pack(fill="x", padx=5, pady=5)

        self.primer_name_var = tk.StringVar()
        self._create_labeled_entry(
            input_frame, "Name:", self.primer_name_var, width=120
        )

        self.primer_seq_var = tk.StringVar()
        self._create_labeled_entry(
            input_frame,
            "Sequence:",
            self.primer_seq_var,
            width=140,
            fill="x",
            expand=True,
        )

        ctk.CTkButton(
            input_frame, text="Add", command=self.add_primer, width=60
        ).pack(side="left", padx=5)
        ctk.CTkButton(
            input_frame,
            text="Delete",
            command=self.delete_primer,
            width=60,
            fg_color="red",
            hover_color="darkred",
        ).pack(side="left", padx=5)
        ctk.CTkButton(
            input_frame, text="Analyze", command=self.analyze_primer, width=80
        ).pack(side="left", padx=5)

        # Listbox for keys - CTk doesn't have Listbox. Use standard tk.Listbox
        # or a ScrollableFrame. Keeping tk.Listbox for now for simplicity,
        # maybe wrap in a frame.
        self.primers_list = tk.Listbox(primers_frame, height=5)
        self.primers_list.pack(fill="both", expand=True, padx=5, pady=5)
        self.primers_data: list[Primer] = []

        # Context menu for listbox
        self.list_menu = tk.Menu(self.primers_list, tearoff=0)
        self.list_menu.add_command(label="Delete", command=self.delete_primer)

        def show_list_menu(event: tk.Event) -> None:
            # Select the item under mouse if not already selected
            index = self.primers_list.nearest(event.y)  # type: ignore[no-untyped-call]
            if index not in self.primers_list.curselection():  # type: ignore[no-untyped-call]
                self.primers_list.selection_clear(0, tk.END)
                self.primers_list.selection_set(index)
                self.primers_list.activate(index)
            self.list_menu.tk_popup(event.x_root, event.y_root)

        self.primers_list.bind("<Button-3>", show_list_menu)

        # --- Settings and Action ---
        settings_frame = ctk.CTkFrame(self)
        settings_frame.pack(fill="x", padx=10, pady=5)

        self.primability_var = tk.DoubleVar(value=settings.primability_cutoff)
        self.primability_entry = self._create_labeled_entry(
            settings_frame,
            "Primability Cutoff:",
            self.primability_var,
            width=60,
        )

        self.stability_var = tk.DoubleVar(value=settings.stability_cutoff)
        self.stability_entry = self._create_labeled_entry(
            settings_frame,
            "Stability Cutoff:",
            self.stability_var,
            width=60,
        )

        self.simulate_btn = ctk.CTkButton(
            settings_frame, text="Simulate PCR", command=self.simulate_pcr
        )
        self.simulate_btn.pack(side="right", padx=10)

        # --- Results Section ---
        self.tree = ttk.Treeview(
            self,
            columns=("seq", "len", "start", "end", "fwd", "rev", "q_score"),
            show="headings",
        )
        self.tree.heading("seq", text="Sequence")
        self.tree.heading("len", text="Length")
        self.tree.heading("start", text="Start")
        self.tree.heading("end", text="End")
        self.tree.heading("fwd", text="Fwd Origin")
        self.tree.heading("rev", text="Rev Origin")
        self.tree.heading("q_score", text="Q-Score")
        self.tree.pack(fill="both", expand=True, padx=10, pady=10)

        # --- Map View ---
        # Map frame to hold map and scrollbar
        self.map_frame = ctk.CTkFrame(self)
        self.map_frame.pack(fill="both", expand=True, padx=10, pady=5)
        self.map_frame.grid_columnconfigure(0, weight=1)
        self.map_frame.grid_rowconfigure(0, weight=1)

        self.map_view = MapView(self.map_frame, highlightthickness=0)
        self.map_view.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

        self.map_scrollbar = ctk.CTkScrollbar(
            self.map_frame, command=self.map_view.yview
        )
        self.map_scrollbar.grid(row=0, column=1, sticky="ns")

        self.map_view.configure(yscrollcommand=self.map_scrollbar.set)

    def add_primer(self) -> None:
        """Parse inputs and add a new primer to the list."""
        name = self.primer_name_var.get().strip()
        seq = self.primer_seq_var.get().strip()

        if not seq:
            messagebox.showerror("Error", "Primer sequence is required.")
            return

        try:
            # Basic validation by creating the object
            p = Primer(seq, name if name else None)
            self.primers_data.append(p)
            display_text = f"{p.name}: {p.seq}"
            self.primers_list.insert(tk.END, display_text)
            self.primer_name_var.set("")
            self.primer_seq_var.set("")
        except Exception as e:
            messagebox.showerror("Error", f"Invalid primer: {e}")

    def delete_primer(self) -> None:
        """Remove the currently selected primer from the list."""
        selection = self.primers_list.curselection()  # type: ignore[no-untyped-call]
        if not selection:
            messagebox.showwarning(
                "Warning", "Please select a primer to delete."
            )
            return

        index = selection[0]

        # Remove from data
        try:
            del self.primers_data[index]
        except IndexError:
            messagebox.showerror(
                "Error", "Could not delete primer: Index out of range."
            )
            return

        # Remove from listbox
        self.primers_list.delete(index)

    def analyze_primer(self) -> None:
        """Open the PrimerStatsDialog for the selected primer."""
        # Get selected primer
        selection = self.primers_list.curselection()  # type: ignore[no-untyped-call]
        if not selection:
            messagebox.showwarning(
                "Warning", "Please select a primer from the list."
            )
            return

        index = selection[0]
        primer = self.primers_data[index]

        # Get Template
        template_str = self.template_text.get("1.0", tk.END).strip()
        if not template_str:
            messagebox.showerror("Error", "Template DNA is required.")
            return

        PrimerStatsDialog(self, primer, template_str)

    def simulate_pcr(self) -> None:
        """Execute the PCR simulation and display potential amplicons."""
        # 1. Get Template
        template_str = self.template_text.get("1.0", tk.END).strip()
        if not template_str:
            messagebox.showerror("Error", "Template DNA is required.")
            return

        try:
            template = DNA(template_str, DNAType.LINEAR, "Main Template")
        except Exception as e:
            messagebox.showerror("Error", f"Invalid Template: {e}")
            return

        # 2. Get Settings
        try:
            # Update global settings
            settings.primability_cutoff = self.primability_var.get()
            settings.stability_cutoff = self.stability_var.get()
        except ValueError:
            messagebox.showerror(
                "Error", "Invalid numeric values for settings."
            )
            return

        # UI Updates
        self.simulate_btn.configure(state="disabled")
        self.configure(cursor="watch")

        # Prepare data for thread
        sim_settings = copy.deepcopy(settings)
        sim_primers = list(self.primers_data)

        self.simulation_result = None
        self.simulation_error = None

        def worker() -> None:
            try:
                # 3. Build Generator
                generator = AmpliconGenerator(template)

                # 4. Add Repliconfs (Combinations of Primer + Template +
                # Settings)
                for primer in sim_primers:
                    try:
                        rc = Repliconf(template, primer, sim_settings)
                        rc.search()  # Important: Must search for origins first!
                        generator.add_repliconf(rc)
                    except Exception as e:
                        print(f"Failed to process primer {primer.name}: {e}")

                # 5. Generate
                self.simulation_result = generator.get_amplicons()
            except Exception as e:
                self.simulation_error = e

        self.simulation_thread = threading.Thread(target=worker, daemon=True)
        self.simulation_thread.start()
        self.after(100, self.check_simulation)

    def check_simulation(self) -> None:
        """Check if the simulation thread is still running."""
        if self.simulation_thread.is_alive():
            self.after(100, self.check_simulation)
        else:
            self.on_simulation_complete()

    def on_simulation_complete(self) -> None:
        """Handle completion of the simulation."""
        self.simulate_btn.configure(state="normal")
        self.configure(cursor="")

        if self.simulation_error:
            messagebox.showerror(
                "Error", f"Simulation failed: {self.simulation_error}"
            )
            return

        # 6. Display Results
        for item in self.tree.get_children():
            self.tree.delete(item)

        if self.simulation_result:
            for amp in self.simulation_result:
                self.tree.insert(
                    "",
                    tk.END,
                    values=(
                        amp.product.seq,
                        len(amp.product),
                        amp.start,
                        amp.end,
                        amp.fwd_origin.name,
                        amp.rev_origin.name,
                        f"{amp.q_score:.2f}",
                    ),
                )

            # Update map view
            try:
                # Reconstruct DNA object from template text to be safe.
                # We need it here for the map context.
                # simulation_result contains amplicons referring to product,
                # but map view needs the full template object.
                template_str = self.template_text.get("1.0", tk.END).strip()
                template = DNA(template_str, DNAType.LINEAR, "Main Template")
                self.map_view.set_data(
                    template, self.simulation_result, self.primers_data
                )
            except Exception as e:
                print(f"Error updating map view: {e}")

    def save_state(self) -> None:
        """Serialize the current application state to a JSON file."""
        file_path = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        if not file_path:
            return

        state = {
            "version": 1,
            "template": self.template_text.get("1.0", tk.END).strip(),
            "primers": [
                {"name": p.name, "sequence": p.seq} for p in self.primers_data
            ],
            "settings": {
                "primability_cutoff": self.primability_var.get(),
                "stability_cutoff": self.stability_var.get(),
            },
        }

        try:
            with open(file_path, "w") as f:
                json.dump(state, f, indent=2)
            messagebox.showinfo("Success", "State saved successfully.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save state: {e}")

    def load_state(self) -> None:
        """Deserialize application state from a JSON file."""
        file_path = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if not file_path:
            return

        try:
            with open(file_path) as f:
                state = json.load(f)

            # Validate version if needed in future
            # if state.get("version") != 1: ...

            # Restore Template
            self.template_text.delete("1.0", tk.END)
            self.template_text.insert(tk.END, state.get("template", ""))

            # Restore Settings
            settings_data = state.get("settings", {})
            if "primability_cutoff" in settings_data:
                self.primability_var.set(settings_data["primability_cutoff"])
            if "stability_cutoff" in settings_data:
                self.stability_var.set(settings_data["stability_cutoff"])

            # Restore Primers
            self.primers_data = []
            self.primers_list.delete(0, tk.END)
            for p_data in state.get("primers", []):
                p = Primer(p_data["sequence"], p_data.get("name"))
                self.primers_data.append(p)
                self.primers_list.insert(tk.END, f"{p.name}: {p.seq}")

            messagebox.showinfo("Success", "State loaded successfully.")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to load state: {e}")


if __name__ == "__main__":
    root = ctk.CTk()
    app = AmplifyPApp(root)
    root.mainloop()
