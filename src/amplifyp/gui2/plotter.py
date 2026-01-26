"""Helper classes for drawing on the MapView canvas."""

import tkinter as tk
from abc import ABC, abstractmethod
from typing import Any


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
