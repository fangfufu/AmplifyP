"""Calculates PCR stats."""

import copy
import json
import os
import sys
import threading
import tkinter as tk
import traceback
from tkinter import filedialog, messagebox, ttk

# Ensure we can import amplifyp when running as a script
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))

from amplifyp.amplicon import AmpliconGenerator
from amplifyp.dna import (
    DNA,
    DNADirection,
    DNAType,
    Primer,
)
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS

settings = copy.deepcopy(DEFAULT_SETTINGS)


class PrimerStatsDialog(tk.Toplevel):
    """Dialog to display primer properties."""

    def __init__(self, parent: tk.Misc, primer: Primer, template_seq: str) -> None:
        """Initialize the dialog."""
        super().__init__(parent)
        self.title(f"Primer Properties: {primer.name}")
        self.geometry("800x400")

        self.primer = primer
        self.template_seq = template_seq

        self.create_widgets()
        self.analyze()

    def create_widgets(self) -> None:
        """Create and arrange widgets."""
        columns = ("pos", "strand", "primability", "stability", "quality")
        self.tree = ttk.Treeview(self, columns=columns, show="headings")
        self.tree.heading("pos", text="Position")
        self.tree.heading("strand", text="Strand")
        self.tree.heading("primability", text="Primability")
        self.tree.heading("stability", text="Stability")
        self.tree.heading("quality", text="Quality")

        # Add scrollbar
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)

        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def analyze(self) -> None:
        """Calculate and display stats."""
        try:
            template = DNA(self.template_seq, DNAType.LINEAR, "Template")
            rc = Repliconf(template, self.primer, settings)
            rc.search()  # This populates origin_idx

            # Helper to add rows
            def add_rows(indices: list[int], direction: DNADirection) -> None:
                strand_str = "Forward" if direction == DNADirection.FWD else "Reverse"
                for i in indices:
                    origin = rc.origin(direction, i)
                    self.tree.insert(
                        "",
                        tk.END,
                        values=(
                            i,
                            strand_str,
                            f"{origin.primability:.4f}",
                            f"{origin.stability:.4f}",
                            f"{origin.quality:.4f}",
                        ),
                    )

            add_rows(rc.origin_db.fwd, DNADirection.FWD)
            add_rows(rc.origin_db.rev, DNADirection.REV)

        except Exception as e:
            traceback.print_exc()
            messagebox.showerror(
                "Analysis Error",
                f"Failed to analyze primer: {e}\nSee console for details.",
            )
            self.destroy()


class AmplifyPApp(ttk.Frame):
    """Main application class for AmplifyP GUI."""

    def __init__(self, master: tk.Tk | None = None) -> None:
        """Initialize the application."""
        super().__init__(master)
        if master:
            master.title("AmplifyP")
            master.geometry("800x600")
            # Create a menu
            self.create_menu()

        self.pack(fill="both", expand=True)
        self.create_widgets()

    def create_menu(self) -> None:
        """Create application menu."""
        if not isinstance(self.master, tk.Tk):
            return

        menubar = tk.Menu(self.master)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Save State", command=self.save_state)
        file_menu.add_command(label="Load State", command=self.load_state)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.quit)

        menubar.add_cascade(label="File", menu=file_menu)
        self.master.configure(menu=menubar)

    def create_context_menu(self, widget: tk.Text | ttk.Entry) -> None:
        """Create a right-click context menu for a widget."""
        menu = tk.Menu(widget, tearoff=0)
        menu.add_command(label="Cut", command=lambda: widget.event_generate("<<Cut>>"))
        menu.add_command(
            label="Copy", command=lambda: widget.event_generate("<<Copy>>")
        )
        menu.add_command(
            label="Paste", command=lambda: widget.event_generate("<<Paste>>")
        )

        def show_menu(event: tk.Event) -> None:
            menu.tk_popup(event.x_root, event.y_root)

        widget.bind("<Button-3>", show_menu)

    def create_widgets(self) -> None:
        """Create and arrange widgets."""
        # --- Template DNA Section ---
        template_frame = ttk.LabelFrame(self, text="Template DNA")
        template_frame.pack(fill="x", padx=10, pady=5)

        self.template_text = tk.Text(template_frame, height=5)
        self.template_text.pack(fill="x", padx=5, pady=5)
        self.create_context_menu(self.template_text)

        # --- Primers Section ---
        primers_frame = ttk.LabelFrame(self, text="Primers")
        primers_frame.pack(fill="both", expand=True, padx=10, pady=5)

        input_frame = ttk.Frame(primers_frame)
        input_frame.pack(fill="x", padx=5, pady=5)

        ttk.Label(input_frame, text="Name:").pack(side="left")
        self.primer_name_var = tk.StringVar()
        name_entry = ttk.Entry(input_frame, textvariable=self.primer_name_var, width=15)
        name_entry.pack(side="left", padx=5)
        self.create_context_menu(name_entry)

        ttk.Label(input_frame, text="Sequence:").pack(side="left")
        self.primer_seq_var = tk.StringVar()
        seq_entry = ttk.Entry(input_frame, textvariable=self.primer_seq_var)
        seq_entry.pack(side="left", fill="x", expand=True, padx=5)
        self.create_context_menu(seq_entry)

        ttk.Button(input_frame, text="Add", command=self.add_primer).pack(side="left")
        ttk.Button(input_frame, text="Delete", command=self.delete_primer).pack(
            side="left", padx=5
        )
        ttk.Button(input_frame, text="Analyze", command=self.analyze_primer).pack(
            side="left"
        )

        # Listbox for keys
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
        settings_frame = ttk.Frame(self)
        settings_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(settings_frame, text="Primability Cutoff:").pack(side="left")
        self.primability_var = tk.DoubleVar(value=settings.primability_cutoff)
        self.primability_entry = ttk.Entry(
            settings_frame, textvariable=self.primability_var, width=10
        )
        self.primability_entry.pack(side="left", padx=5)
        self.create_context_menu(self.primability_entry)

        ttk.Label(settings_frame, text="Stability Cutoff:").pack(side="left")
        self.stability_var = tk.DoubleVar(value=settings.stability_cutoff)
        self.stability_entry = ttk.Entry(
            settings_frame, textvariable=self.stability_var, width=10
        )
        self.stability_entry.pack(side="left", padx=5)
        self.create_context_menu(self.stability_entry)

        self.simulate_btn = ttk.Button(
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

    def add_primer(self) -> None:
        """Add primer to the list."""
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
        """Delete selected primer."""
        selection = self.primers_list.curselection()  # type: ignore[no-untyped-call]
        if not selection:
            messagebox.showwarning("Warning", "Please select a primer to delete.")
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
        """Analyze selected primer against template."""
        # Get selected primer
        selection = self.primers_list.curselection()  # type: ignore[no-untyped-call]
        if not selection:
            messagebox.showwarning("Warning", "Please select a primer from the list.")
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
        """Run the simulation."""
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
            messagebox.showerror("Error", "Invalid numeric values for settings.")
            return

        # UI Updates
        self.simulate_btn.config(state="disabled")
        self.config(cursor="watch")

        # Prepare data for thread
        sim_settings = copy.deepcopy(settings)
        sim_primers = list(self.primers_data)

        self.simulation_result = None
        self.simulation_error = None

        def worker() -> None:
            try:
                # 3. Build Generator
                generator = AmpliconGenerator(template)

                # 4. Add Repliconfs (Combinations of Primer + Template + Settings)
                for primer in sim_primers:
                    try:
                        rc = Repliconf(template, primer, sim_settings)
                        rc.search()  # Important: Must search for origins first!
                        generator.add(rc)
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
        self.simulate_btn.config(state="normal")
        self.config(cursor="")

        if self.simulation_error:
            messagebox.showerror("Error", f"Simulation failed: {self.simulation_error}")
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

    def save_state(self) -> None:
        """Save the current state to a JSON file."""
        file_path = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        if not file_path:
            return

        state = {
            "version": 1,
            "template": self.template_text.get("1.0", tk.END).strip(),
            "primers": [{"name": p.name, "sequence": p.seq} for p in self.primers_data],
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
        """Load state from a JSON file."""
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
    root = tk.Tk()
    app = AmplifyPApp(root)
    root.mainloop()
