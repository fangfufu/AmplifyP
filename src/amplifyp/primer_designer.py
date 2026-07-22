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

"""Primer design module for AmplifyP."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .dimer import PrimerDimer, PrimerDimerGenerator
from .dna import DNA, DNADirection, Primer

if TYPE_CHECKING:
    from collections.abc import Iterator

DEFAULT_PRIMER_DIMER_GENERATOR = PrimerDimerGenerator()


class PrimerDesigner1D:
    """Performs 1D primer design by truncating a DNA sequence.

    Calculates self-dimers for each truncation step from the initial length
    down to target length `n`.

    In Forward mode (`DNADirection.FWD`), self-dimers are generated while
    truncating bases from the 3' end. In Reverse/Backward mode
    (`DNADirection.REV`), self-dimers are generated while truncating bases
    from the 5' end.
    """

    __slots__ = ("_dimers", "_dna", "_generator", "_min_length", "_mode")

    def __init__(
        self,
        dna: DNA,
        min_length: int,
        mode: DNADirection = DNADirection.FWD,
        generator: PrimerDimerGenerator = DEFAULT_PRIMER_DIMER_GENERATOR,
    ) -> None:
        """Initialises a new PrimerDesigner1D object and runs the analysis.

        Args:
            dna (DNA): The input DNA object.
            min_length (int): The minimum primer sequence length to reach via
                truncation.
            mode (DNADirection, optional): The truncation mode.
                `DNADirection.FWD` for 3' end truncation; `DNADirection.REV`
                for 5' end truncation. Defaults to `DNADirection.FWD`.
            generator (PrimerDimerGenerator, optional): Custom primer dimer
                generator to use for self-dimer analysis. Defaults to
                `PrimerDimerGenerator()`.

        Raises:
            ValueError: If minimum length is non-positive or greater than
                sequence length.
        """
        if min_length <= 0:
            raise ValueError("Target length n must be greater than 0.")
        if len(dna.seq_upper) < min_length:
            msg = (
                f"Target length n ({min_length}) cannot exceed initial "
                f"sequence length ({len(dna.seq_upper)})."
            )
            raise ValueError(msg)

        self._dna: DNA = dna
        self._min_length: int = min_length
        self._mode: DNADirection = mode
        self._generator: PrimerDimerGenerator = generator
        self._dimers: list[PrimerDimer] = []

        self._analyse(dna.seq_upper)

    @property
    def dna(self) -> DNA:
        """The input DNA object."""
        return self._dna

    @property
    def min_length(self) -> int:
        """The minimum primer sequence length."""
        return self._min_length

    @property
    def mode(self) -> DNADirection:
        """The truncation mode used for primer design."""
        return self._mode

    @property
    def generator(self) -> PrimerDimerGenerator:
        """The primer dimer generator used for self-dimer analysis."""
        return self._generator

    @property
    def results(self) -> tuple[PrimerDimer, ...]:
        """The stored self-dimers from initial sequence down to length n."""
        return tuple(self._dimers)

    @property
    def best_primer(self) -> int:
        """Return the index of the PrimerDimer with the lowest quality score.

        A lower self-dimer quality score indicates lower dimer formation
        potential, representing the optimal primer choice.

        Returns:
            int: The 0-indexed step position of the self-dimer with the lowest
                quality score.

        Raises:
            RuntimeError: If no results have been generated.
        """
        if not self._dimers:
            raise RuntimeError("No analysis steps recorded.")
        return min(
            range(len(self._dimers)), key=lambda i: self._dimers[i].quality
        )

    def filter_by_quality(
        self, max_quality: float, sort: bool = True
    ) -> tuple[int, ...]:
        """Filter results by quality score.

        Returns 0-indexed step positions for self-dimers with quality
        score <= `max_quality`.

        Args:
            max_quality (float): The maximum quality score threshold.
            sort (bool, optional): If True, sort step positions by quality
                score in ascending order. Defaults to True.

        Returns:
            tuple[int, ...]: Step indices matching the criteria.
        """
        matching_indices = [
            i for i, d in enumerate(self._dimers) if d.quality <= max_quality
        ]
        if sort:
            matching_indices.sort(key=lambda i: self._dimers[i].quality)
        return tuple(matching_indices)

    def sorted_by_quality(self, reverse: bool = False) -> tuple[int, ...]:
        """Return self-dimer result step indices sorted by quality score.

        Args:
            reverse (bool, optional): If True, sort in descending order of
                quality score (worst first). Defaults to False (best first).

        Returns:
            tuple[int, ...]: Step indices sorted by quality score.
        """
        return tuple(
            sorted(
                range(len(self._dimers)),
                key=lambda i: self._dimers[i].quality,
                reverse=reverse,
            )
        )

    def __len__(self) -> int:
        """Return the number of design truncation steps."""
        return len(self._dimers)

    def __getitem__(self, index: int) -> PrimerDimer:
        """Get the self-dimer at 0-indexed step position `index`.

        Args:
            index (int): The 0-indexed step position.

        Returns:
            PrimerDimer: The self-dimer object at the given index.
        """
        return self._dimers[index]

    def __iter__(self) -> Iterator[PrimerDimer]:
        """Iterate over stored self-dimers."""
        return iter(self._dimers)

    def __contains__(self, item: object) -> bool:
        """Check if a PrimerDimer is contained in the results."""
        return item in self._dimers

    def __repr__(self) -> str:
        """Return an unambiguous string representation of PrimerDesigner1D."""
        return (
            f"PrimerDesigner1D(dna={self._dna!r}, "
            f"min_length={self._min_length}, mode={self._mode!r})"
        )

    def __str__(self) -> str:
        """Return a user-friendly string representation of PrimerDesigner1D."""
        return (
            f"PrimerDesigner1D({len(self._dimers)} steps, "
            f"min_length={self._min_length}, mode={self._mode.name})"
        )

    def get_dimer(self, index: int) -> PrimerDimer:
        """Query stored self-dimer by 0-indexed step position `index`.

        Args:
            index (int): The 0-indexed step position.

        Returns:
            PrimerDimer: The self-dimer corresponding to the step at `index`.
        """
        return self[index]

    def _analyse(self, initial_seq: str) -> None:
        """Perform the 1D truncation analysis on the initial sequence.

        This method iteratively truncates the sequence from either the 5' or 3'
        end (depending on the mode) and generates a self-dimer for each
        resulting sequence until the minimum length is reached.

        Args:
            initial_seq (str): The initial sequence string to analyze.
        """
        self._dimers.clear()
        current_seq = initial_seq

        while len(current_seq) >= self._min_length:
            primer = Primer(current_seq)
            dimer = self._generator.generate_primer_dimer(primer, primer)
            self._dimers.append(dimer)

            if len(current_seq) == self._min_length:
                break

            if self._mode == DNADirection.FWD:
                # Chop off a base from the 3' end (rightmost base)
                current_seq = current_seq[:-1]
            else:
                # Chop off a base from the 5' end (leftmost base)
                current_seq = current_seq[1:]
