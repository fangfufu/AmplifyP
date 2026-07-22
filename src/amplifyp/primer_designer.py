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

__all__ = ["PrimerDesigner1D"]


class PrimerDesigner1D:
    """Performs 1D primer design by truncating a DNA sequence.

    Calculates self-dimers for each truncation step from the initial length
    down to target length `n`.

    In Forward mode (`DNADirection.FWD`), self-dimers are generated while
    truncating bases from the 3' end. In Reverse/Backward mode
    (`DNADirection.REV`), self-dimers are generated while truncating bases
    from the 5' end.
    """

    __slots__ = ("_dimers", "_dna", "_generator", "_mode", "_n")

    def __init__(
        self,
        dna: DNA,
        n: int | None = None,
        mode: DNADirection = DNADirection.FWD,
        *,
        min_length: int | None = None,
        target_length: int | None = None,
        generator: PrimerDimerGenerator | None = None,
    ) -> None:
        """Initialises a new PrimerDesigner1D object and runs the analysis.

        Args:
            dna (DNA): The input DNA object.
            n (int, optional): The minimum primer sequence length to reach via
                truncation.
            mode (DNADirection, optional): The truncation mode.
                `DNADirection.FWD` for 3' end truncation; `DNADirection.REV`
                for 5' end truncation. Defaults to `DNADirection.FWD`.
            min_length (int, optional): Explicit alias for minimum primer length
                `n`.
            target_length (int, optional): Alias for `n`.
            generator (PrimerDimerGenerator, optional): Custom primer dimer
                generator to use for self-dimer analysis.

        Raises:
            TypeError: If `dna` is not a DNA object, `mode` is not a
                `DNADirection` instance, or `generator` is invalid.
            ValueError: If minimum length is non-positive or greater than
                sequence length.
        """
        if not isinstance(dna, DNA):  # pyright: ignore[reportUnnecessaryIsInstance]
            raise TypeError("dna must be a DNA object.")
        if not isinstance(mode, DNADirection):  # pyright: ignore[reportUnnecessaryIsInstance]
            raise TypeError("mode must be a DNADirection instance.")
        if generator is not None and not isinstance(  # pyright: ignore[reportUnnecessaryIsInstance]
            generator, PrimerDimerGenerator
        ):
            raise TypeError(
                "generator must be a PrimerDimerGenerator instance."
            )

        length_param = min_length if min_length is not None else target_length
        if length_param is not None and n is not None and length_param != n:
            raise ValueError(
                "Cannot specify conflicting values for n and min_length."
            )

        eff_n = length_param if length_param is not None else n
        if eff_n is None:
            raise ValueError(
                "Minimum primer length (n or min_length) must be provided."
            )

        clean_seq = dna.seq_upper

        if eff_n <= 0:
            raise ValueError("Target length n must be greater than 0.")
        if len(clean_seq) < eff_n:
            msg = (
                f"Target length n ({eff_n}) cannot exceed initial "
                f"sequence length ({len(clean_seq)})."
            )
            raise ValueError(msg)

        self._dna: DNA = dna
        self._n: int = eff_n
        self._mode: DNADirection = mode
        self._generator: PrimerDimerGenerator = (
            generator if generator is not None else PrimerDimerGenerator()
        )
        self._dimers: list[PrimerDimer] = []

        self._analyse(clean_seq)

    @property
    def dna(self) -> DNA:
        """The input DNA object."""
        return self._dna

    @property
    def n(self) -> int:
        """The target length n."""
        return self._n

    @property
    def min_length(self) -> int:
        """The minimum primer sequence length."""
        return self._n

    @property
    def target_length(self) -> int:
        """The minimum primer sequence length (alias for min_length)."""
        return self._n

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
    def best_dimer(self) -> PrimerDimer:
        """Return the PrimerDimer with the lowest quality score among steps.

        A lower self-dimer quality score indicates lower dimer formation
        potential, representing the optimal primer choice.

        Returns:
            PrimerDimer: The self-dimer object with the lowest quality score.

        Raises:
            RuntimeError: If no results have been generated.
        """
        if not self._dimers:
            raise RuntimeError("No analysis steps recorded.")
        return min(self._dimers, key=lambda d: d.quality)

    @property
    def worst_dimer(self) -> PrimerDimer:
        """Return the PrimerDimer with the highest quality score among steps.

        Returns:
            PrimerDimer: The self-dimer object with the highest quality score.

        Raises:
            RuntimeError: If no results have been generated.
        """
        if not self._dimers:
            raise RuntimeError("No analysis steps recorded.")
        return max(self._dimers, key=lambda d: d.quality)

    def filter_by_quality(self, max_quality: float) -> tuple[PrimerDimer, ...]:
        """Filter results by quality score.

        Returns self-dimers with quality score <= `max_quality`.

        Args:
            max_quality (float): The maximum quality score threshold.

        Returns:
            tuple[PrimerDimer, ...]: Self-dimers matching the criteria.
        """
        return tuple(d for d in self._dimers if d.quality <= max_quality)

    def sorted_by_quality(
        self, reverse: bool = False
    ) -> tuple[PrimerDimer, ...]:
        """Return self-dimer results sorted by quality score.

        Args:
            reverse (bool, optional): If True, sort in descending order of
                quality score (worst first). Defaults to False (best first).

        Returns:
            tuple[PrimerDimer, ...]: Self-dimers sorted by quality score.
        """
        return tuple(
            sorted(self._dimers, key=lambda d: d.quality, reverse=reverse)
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
            f"PrimerDesigner1D(dna={self._dna!r}, n={self._n}, "
            f"mode={self._mode!r})"
        )

    def __str__(self) -> str:
        """Return a user-friendly string representation of PrimerDesigner1D."""
        return (
            f"PrimerDesigner1D({len(self._dimers)} steps, "
            f"target length={self._n}, mode={self._mode.name})"
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
        self._dimers.clear()
        current_seq = initial_seq

        while len(current_seq) >= self._n:
            primer = Primer(current_seq)
            dimer = self._generator.generate_primer_dimer(primer, primer)
            self._dimers.append(dimer)

            if len(current_seq) == self._n:
                break

            if self._mode == DNADirection.FWD:
                # Chop off a base from the 3' end (rightmost base)
                current_seq = current_seq[:-1]
            else:
                # Chop off a base from the 5' end (leftmost base)
                current_seq = current_seq[1:]
