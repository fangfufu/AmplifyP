# Copyright (C) 2026 AmplifyP Contributors
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

"""Directed index classes for AmplifyP."""

from __future__ import annotations

from dataclasses import dataclass, field

from .dna import DNADirection


@dataclass(slots=True, frozen=True)
class DirIdx:
    """A class representing a directed index.

    This class encapsulates a specific position (index) on a DNA strand,
    along with the direction (forward or reverse) of that strand, and
    optional cached primability and stability scores.

    Attributes:
        direction (DNADirection): The direction of the DNA strand.
        index (int): The integer index position.
        primability (float | None): Cached primability score.
        stability (float | None): Cached stability score.
    """

    direction: DNADirection
    index: int
    primability: float | None = field(default=None, compare=False, repr=False)
    stability: float | None = field(default=None, compare=False, repr=False)

    def __hash__(self) -> int:
        """Compute the hash of the DirIdx.

        Returns:
            int: The hash based on direction and index.
        """
        return hash((self.direction, self.index))

    def __int__(self) -> int:
        """Return the index value as an integer.

        Returns:
            int: The index.
        """
        return self.index

    def __index__(self) -> int:
        """Return the index value for slicing/indexing usage.

        Returns:
            int: The index.
        """
        return self.index

    def __add__(self, other: object) -> DirIdx:
        """Add an integer or another DirIdx to this index.

        Args:
            other (int | DirIdx): The value to add.

        Returns:
            DirIdx: A new DirIdx with the updated index value. The direction
                remains unchanged.
        """
        if isinstance(other, int):
            return DirIdx(self.direction, self.index + other)
        if not isinstance(other, DirIdx):
            return NotImplemented
        return DirIdx(self.direction, self.index + other.index)

    def __sub__(self, other: object) -> DirIdx:
        """Subtract an integer or another DirIdx from this index.

        Args:
            other (int | DirIdx): The value to subtract.

        Returns:
            DirIdx: A new DirIdx with the updated index value. The direction
                remains unchanged.
        """
        if isinstance(other, int):
            return DirIdx(self.direction, self.index - other)
        if not isinstance(other, DirIdx):
            return NotImplemented
        return DirIdx(self.direction, self.index - other.index)

    def __eq__(self, other: object) -> bool:
        """Check equality with an integer or another DirIdx.

        Args:
            other (object): The object to compare.

        Returns:
            bool: True if equal, False otherwise. If comparing with int, checks
                index only.
        """
        if isinstance(other, int):
            return self.index == other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.direction == other.direction and self.index == other.index

    def __lt__(self, other: object) -> bool:
        """Check if index is less than another.

        Args:
            other (object): Value to compare.

        Returns:
            bool: True if less than.
        """
        if isinstance(other, int):
            return self.index < other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index < other.index

    def __gt__(self, other: object) -> bool:
        """Check if index is greater than another.

        Args:
            other (object): Value to compare.

        Returns:
            bool: True if greater than.
        """
        if isinstance(other, int):
            return self.index > other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index > other.index

    def __le__(self, other: object) -> bool:
        """Check if index is less than or equal to another.

        Args:
            other (object): Value to compare.

        Returns:
            bool: True if less than or equal.
        """
        if isinstance(other, int):
            return self.index <= other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index <= other.index

    def __ge__(self, other: object) -> bool:
        """Check if index is greater than or equal to another.

        Args:
            other (object): Value to compare.

        Returns:
            bool: True if greater than or equal.
        """
        if isinstance(other, int):
            return self.index >= other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index >= other.index

    def __str__(self) -> str:
        """Return the string representation of the index.

        Returns:
            str: The index value as a string.
        """
        return f"{self.index}"


@dataclass(slots=True)
class DirIdxDb:
    """A database for storing valid replication origin locations.

    This container holds lists of indices where valid replication origins were
    found, segregated by strand direction (forward and reverse). It also tracks
    whether a search operation has been performed.

    Attributes:
        fwd (list[DirIdx]): A list of locations (DirIdx) for origins on the
            forward strand.
        rev (list[DirIdx]): A list of locations (DirIdx) for origins on the
            reverse strand.
        searched (bool): Indicates if the search has been executed. Defaults to
            False.
    """

    fwd: list[DirIdx]
    rev: list[DirIdx]
    searched: bool = False

    def clear(self) -> None:
        """Clear the database.

        Removes stored indices and resets search flag.
        """
        self.fwd.clear()
        self.rev.clear()
        self.searched = False

    def __getitem__(self, key: tuple[DNADirection, int]) -> DirIdx:
        """Retrieve a stored DirIdx by direction and list index.

        Args:
            key (tuple[DNADirection, int]): A tuple specifying the direction
                (FWD/REV) and the integer index within that list.

        Returns:
            DirIdx: The requested directed index.
        """
        direction, i = key
        if direction == DNADirection.FWD:
            return self.fwd[i]
        else:
            return self.rev[i]
