"""The graphical MapView for displaying PCR results."""

import tkinter as tk
from typing import Any

from amplifyp.amplicon import Amplicon
from amplifyp.dna import DNA
from amplifyp.gui2.plotter import (
    LineThing,
    PlotterThing,
    RectThing,
    StringThing,
)


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

            # Check for wrapping (start > end implies wrapping in circular DNA)
            if end_x < start_x:
                # Circular wrapping visualization
                # Segment 1: start_x to end of template
                right_limit = basex(self.target_length)
                self.plotters.append(
                    RectThing(
                        start_x,
                        current_y,
                        right_limit,
                        current_y + 10,
                        fill="black",
                        outline="",
                    )
                )

                # Segment 2: start of template to end_x
                left_limit = basex(0)
                self.plotters.append(
                    RectThing(
                        left_limit,
                        current_y,
                        end_x,
                        current_y + 10,
                        fill="black",
                        outline="",
                    )
                )

                # Label placement on the longer segment
                len_seg1 = right_limit - start_x
                len_seg2 = end_x - left_limit

                if len_seg1 >= len_seg2:
                    center_x = start_x + len_seg1 / 2
                else:
                    center_x = left_limit + len_seg2 / 2

            else:
                # Linear visualization
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
                width = end_x - start_x
                center_x = start_x + width / 2

            # Size label centered
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
