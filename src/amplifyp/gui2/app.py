"""The main AmplifyP GUI application (v2)."""

import threading
import tkinter as tk

import customtkinter as ctk

from amplifyp.amplicon import Amplicon, AmpliconGenerator
from amplifyp.dna import DNA, DNAType, Primer
from amplifyp.gui2.map_view import MapView
from amplifyp.repliconf import Repliconf
from amplifyp.settings import GLOBAL_REPLICATION_SETTINGS

ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")


class AmplifyPApp(ctk.CTk):  # type: ignore[misc]
    """Main application window."""

    def __init__(self) -> None:
        """Initialize the AmplifyPApp."""
        super().__init__()

        self.title("AmplifyP")
        self.geometry("1000x800")

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=0)  # Input area
        self.grid_rowconfigure(1, weight=1)  # Map view

        self.create_widgets()

        # State
        self.primers: list[Primer] = []
        self.template: DNA | None = None
        self.amplicons: list[Amplicon] = []

    def create_widgets(self) -> None:
        """Create UI elements."""
        # --- Top Frame (Inputs) ---
        self.top_frame = ctk.CTkFrame(self)
        self.top_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=10)

        # Template Input
        ctk.CTkLabel(self.top_frame, text="Template DNA:").grid(
            row=0, column=0, sticky="w", padx=5
        )
        self.template_text = ctk.CTkTextbox(
            self.top_frame, height=60, width=600
        )
        self.template_text.grid(row=0, column=1, rowspan=2, padx=5, pady=5)

        # Controls
        self.btn_simulate = ctk.CTkButton(
            self.top_frame, text="Simulate PCR", command=self.run_simulation
        )
        self.btn_simulate.grid(row=0, column=2, padx=10, pady=5)

        # Primers Input (Minimal for now)
        ctk.CTkLabel(self.top_frame, text="Add Primer (Seq):").grid(
            row=2, column=0, sticky="w", padx=5
        )
        self.primer_entry = ctk.CTkEntry(self.top_frame, width=400)
        self.primer_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")
        self.btn_add_primer = ctk.CTkButton(
            self.top_frame, text="Add", command=self.add_primer, width=60
        )
        self.btn_add_primer.grid(
            row=2, column=1, padx=(410, 5), pady=5, sticky="w"
        )

        self.lbl_primer_count = ctk.CTkLabel(self.top_frame, text="Primers: 0")
        self.lbl_primer_count.grid(row=2, column=2, padx=10)

        # --- Middle Frame (Map) ---
        # We need to wrap MapView (tk.Canvas) in a CTk frame or place it
        # directly.
        # Placing directly in grid.
        self.map_frame = ctk.CTkFrame(self)
        self.map_frame.grid(
            row=1, column=0, sticky="nsew", padx=10, pady=(0, 10)
        )
        self.map_frame.grid_columnconfigure(0, weight=1)
        self.map_frame.grid_rowconfigure(0, weight=1)

        self.map_view = MapView(self.map_frame, highlightthickness=0)
        self.map_view.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

        # Scrollbar for map
        self.scrollbar = ctk.CTkScrollbar(
            self.map_frame, command=self.map_view.yview
        )
        self.scrollbar.grid(row=0, column=1, sticky="ns")
        self.map_view.configure(yscrollcommand=self.scrollbar.set)

        # --- Bottom Frame (Output/Log) ---
        self.log_text = ctk.CTkTextbox(self, height=100)
        self.log_text.grid(row=2, column=0, sticky="ew", padx=10, pady=(0, 10))

    def log(self, message: str) -> None:
        """Log a message to the bottom text area."""
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)

    def add_primer(self) -> None:
        """Add primer from entry."""
        seq = self.primer_entry.get().strip()
        if not seq:
            return

        try:
            name = f"P{len(self.primers) + 1}"
            p = Primer(seq, name=name)
            self.primers.append(p)
            self.lbl_primer_count.configure(
                text=f"Primers: {len(self.primers)}"
            )
            self.primer_entry.delete(0, tk.END)
            self.log(f"Added primer: {name} ({len(seq)} bp)")
        except Exception as e:
            self.log(f"Error adding primer: {e}")

    def run_simulation(self) -> None:
        """Run the PCR simulation."""
        t_seq = self.template_text.get("1.0", tk.END).strip()
        if not t_seq:
            self.log("Error: No template sequence.")
            return

        try:
            self.template = DNA(t_seq, DNAType.LINEAR, "Template")
        except Exception as e:
            self.log(f"Error parsing template: {e}")
            return

        self.btn_simulate.configure(state="disabled")
        self.log("Starting simulation...")

        def worker() -> None:
            if self.template is None:
                return
            try:
                gen = AmpliconGenerator(self.template)
                for p in self.primers:
                    rc = Repliconf(
                        self.template, p, GLOBAL_REPLICATION_SETTINGS
                    )
                    rc.search()
                    gen.add_repliconf(rc)

                self.amplicons = gen.get_amplicons()
                self.after(0, self.on_simulation_complete)
            except Exception as e:
                err_msg = str(e)
                self.after(0, lambda: self.log(f"Simulation failed: {err_msg}"))
                self.after(
                    0, lambda: self.btn_simulate.configure(state="normal")
                )

        threading.Thread(target=worker, daemon=True).start()

    def on_simulation_complete(self) -> None:
        """Handle results."""
        self.log(f"Simulation complete. Found {len(self.amplicons)} amplicons.")
        self.map_view.set_data(self.template, self.amplicons, self.primers)  # type: ignore[arg-type]
        self.btn_simulate.configure(state="normal")


if __name__ == "__main__":
    app = AmplifyPApp()
    app.mainloop()
