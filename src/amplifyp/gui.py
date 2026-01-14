# -*- coding: utf-8 -*-
"""Calculates PCR stats."""

import os
import sys
import traceback
import tkinter as tk
from tkinter import messagebox, ttk
from typing import List

# Ensure we can import amplifyp when running as a script
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))

from amplifyp.amplicon import AmpliconGenerator
from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS
import copy

settings = copy.deepcopy(DEFAULT_SETTINGS)


class PrimerStatsDialog(tk.Toplevel):
    """Dialog to display primer properties."""

    def __init__(self, parent: tk.Tk, primer: Primer, template_seq: str) -> None:
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
            # Use strictness 0 to find all matches (or very low matches)
            # Use strictness 0 to find all matches (or very low matches)
            # Use a copy of the global settings to avoid modifying them
            local_settings = copy.deepcopy(settings)
            local_settings.primability_cutoff = 0.0
            local_settings.stability_cutoff = 0.0

            rc = Repliconf(template, self.primer, local_settings)
            rc.search()  # This populates origin_idx

            # Helper to add rows
            def add_rows(indices: List[int], direction: DNADirection) -> None:
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

            add_rows(rc.origin_idx.fwd, DNADirection.FWD)
            add_rows(rc.origin_idx.rev, DNADirection.REV)

        except Exception as e:
            traceback.print_exc()
            messagebox.showerror(
                "Analysis Error",
                f"Failed to analyze primer: {e}\nSee console for details.",
            )
            self.destroy()


class AmplifyPApp(tk.Tk):
    """Main application class for AmplifyP GUI."""

    def __init__(self) -> None:
        """Initialize the application."""
        super().__init__()
        self.title("AmplifyP")
        self.geometry("800x600")

        self.create_widgets()

    def create_widgets(self) -> None:
        """Create and arrange widgets."""
        # --- Template DNA Section ---
        template_frame = ttk.LabelFrame(self, text="Template DNA")
        template_frame.pack(fill="x", padx=10, pady=5)

        self.template_text = tk.Text(template_frame, height=5)
        self.template_text.pack(fill="x", padx=5, pady=5)

        # --- Primers Section ---
        primers_frame = ttk.LabelFrame(self, text="Primers")
        primers_frame.pack(fill="both", expand=True, padx=10, pady=5)

        input_frame = ttk.Frame(primers_frame)
        input_frame.pack(fill="x", padx=5, pady=5)

        ttk.Label(input_frame, text="Name:").pack(side="left")
        self.primer_name_var = tk.StringVar()
        ttk.Entry(input_frame, textvariable=self.primer_name_var, width=15).pack(
            side="left", padx=5
        )

        ttk.Label(input_frame, text="Sequence:").pack(side="left")
        self.primer_seq_var = tk.StringVar()
        ttk.Entry(input_frame, textvariable=self.primer_seq_var).pack(
            side="left", fill="x", expand=True, padx=5
        )

        ttk.Button(input_frame, text="Add", command=self.add_primer).pack(side="left")
        ttk.Button(input_frame, text="Analyze", command=self.analyze_primer).pack(
            side="left", padx=5
        )

        # Listbox for keys
        self.primers_list = tk.Listbox(primers_frame, height=5)
        self.primers_list.pack(fill="both", expand=True, padx=5, pady=5)
        self.primers_data: List[Primer] = []

        # --- Settings and Action ---
        settings_frame = ttk.Frame(self)
        settings_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(settings_frame, text="Primability Cutoff:").pack(side="left")
        self.primability_var = tk.DoubleVar(value=settings.primability_cutoff)
        ttk.Entry(settings_frame, textvariable=self.primability_var, width=10).pack(
            side="left", padx=5
        )

        ttk.Label(settings_frame, text="Stability Cutoff:").pack(side="left")
        self.stability_var = tk.DoubleVar(value=settings.stability_cutoff)
        ttk.Entry(settings_frame, textvariable=self.stability_var, width=10).pack(
            side="left", padx=5
        )

        ttk.Button(settings_frame, text="Simulate PCR", command=self.simulate_pcr).pack(
            side="right", padx=10
        )

        # --- Results Section ---
        self.tree = ttk.Treeview(
            self,
            columns=("seq", "len", "start", "end", "fwd", "rev"),
            show="headings",
        )
        self.tree.heading("seq", text="Sequence")
        self.tree.heading("len", text="Length")
        self.tree.heading("start", text="Start")
        self.tree.heading("end", text="End")
        self.tree.heading("fwd", text="Fwd Origin")
        self.tree.heading("rev", text="Rev Origin")
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

        # 3. Build Generator
        generator = AmpliconGenerator(template)

        # 4. Add Repliconfs (Combinations of Primer + Template + Settings)
        for primer in self.primers_data:
            try:
                rc = Repliconf(template, primer, settings)
                rc.search()  # Important: Must search for origins first!
                generator.add(rc)
            except Exception as e:
                print(f"Failed to process primer {primer.name}: {e}")

        # 5. Generate
        try:
            amplicons = generator.generate_amplicons()
        except Exception as e:
            messagebox.showerror("Error", f"Simulation failed: {e}")
            return

        # 6. Display Results
        for item in self.tree.get_children():
            self.tree.delete(item)

        for amp in amplicons:
            self.tree.insert(
                "",
                tk.END,
                values=(
                    amp.sequence.seq,
                    len(amp.sequence),
                    "N/A",
                    "N/A",
                    amp.fwd_origin.name,
                    amp.rev_origin.name,
                ),
            )


if __name__ == "__main__":
    app = AmplifyPApp()
    app.mainloop()
