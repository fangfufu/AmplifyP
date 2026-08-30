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
from .repliconf import Repliconf

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

    __slots__ = (
        "_dimers",
        "_dna",
        "_generator",
        "_max_origin_count",
        "_max_overlap",
        "_min_length",
        "_mode",
        "_template",
        "_threshold",
    )

    def __init__(
        self,
        dna: DNA,
        min_length: int,
        mode: DNADirection = DNADirection.FWD,
        generator: PrimerDimerGenerator = DEFAULT_PRIMER_DIMER_GENERATOR,
        threshold: float | None = None,
        max_overlap: int | None = None,
        template: DNA | None = None,
        max_origin_count: int | None = None,
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
            threshold (float | None, optional): Upper bound for quality score
                filter. Defaults to None.
            max_overlap (int | None, optional): Upper bound for overlap length
                filter. Defaults to None.
            template (DNA | None, optional): Optional template DNA sequence
                used to calculate replication origins. Defaults to None.
            max_origin_count (int | None, optional): Upper bound for number of
                replication origins filter. Defaults to None.

        Raises:
            ValueError: If minimum length is non-positive or greater than
                sequence length, if max_origin_count is specified without
                template, or if max_origin_count is negative.
        """
        if min_length <= 0:
            raise ValueError("Target length n must be greater than 0.")
        if len(dna.seq_upper) < min_length:
            msg = (
                f"Target length n ({min_length}) cannot exceed initial "
                f"sequence length ({len(dna.seq_upper)})."
            )
            raise ValueError(msg)
        if max_origin_count is not None:
            if template is None:
                msg = (
                    "Template DNA must be provided when max_origin_count "
                    "filter is specified."
                )
                raise ValueError(msg)
            if max_origin_count < 0:
                raise ValueError("max_origin_count must be non-negative.")

        self._dna: DNA = dna
        self._min_length: int = min_length
        self._mode: DNADirection = mode
        self._generator: PrimerDimerGenerator = generator
        self._threshold: float | None = threshold
        self._max_overlap: int | None = max_overlap
        self._template: DNA | None = template
        self._max_origin_count: int | None = max_origin_count
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
    def threshold(self) -> float | None:
        """The maximum quality score cutoff filter, if set."""
        return self._threshold

    @property
    def max_overlap(self) -> int | None:
        """The maximum overlap length cutoff filter, if set."""
        return self._max_overlap

    @property
    def template(self) -> DNA | None:
        """The optional template DNA object, if set."""
        return self._template

    @property
    def max_origin_count(self) -> int | None:
        """The maximum replication origin count cutoff filter, if set."""
        return self._max_origin_count

    @property
    def all_dimers(self) -> tuple[PrimerDimer, ...]:
        """The stored self-dimers from initial sequence down to length n."""
        return tuple(self._dimers)

    @property
    def best_score(self) -> tuple[int, float]:
        """Return the (index, quality) pair of the optimal self-dimer.

        A lower self-dimer quality score indicates lower dimer formation
        potential, representing the optimal primer choice.

        Returns:
            tuple[int, float]: A (step_index, quality_score) tuple for the
                self-dimer with the lowest quality score.

        Raises:
            RuntimeError: If no results have been generated.
        """
        if not self._dimers:
            raise RuntimeError("No analysis steps recorded.")
        return self.quality_score(sorted=True)[0]

    def quality_score(
        self, sorted: bool = False
    ) -> tuple[tuple[int, float], ...]:
        """Return (index, quality) pairs for all self-dimer results.

        Args:
            sorted (bool, optional): If True, sort pairs by quality score in
                ascending order (best quality first). Defaults to False.

        Returns:
            tuple[tuple[int, float], ...]: Sequence of (step_index, quality)
                tuples.
        """
        pairs = [(i, d.quality) for i, d in enumerate(self._dimers)]
        if sorted:
            pairs.sort(key=lambda item: item[1])
        return tuple(pairs)

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
            origin_count: int | None = None
            if self._template is not None:
                repliconf = Repliconf(self._template, primer)
                repliconf.search()
                origin_count = len(repliconf.origin_db.fwd) + len(
                    repliconf.origin_db.rev
                )

            if (
                (self._threshold is None or dimer.quality <= self._threshold)
                and (
                    self._max_overlap is None
                    or dimer.overlap <= self._max_overlap
                )
                and (
                    self._max_origin_count is None
                    or (
                        origin_count is not None
                        and origin_count <= self._max_origin_count
                    )
                )
            ):
                self._dimers.append(dimer)

            if len(current_seq) == self._min_length:
                break

            if self._mode == DNADirection.FWD:
                # Chop off a base from the 3' end (rightmost base)
                current_seq = current_seq[:-1]
            else:
                # Chop off a base from the 5' end (leftmost base)
                current_seq = current_seq[1:]
