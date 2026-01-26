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

"""Primer dimer analysis module for AmplifyP."""

from __future__ import annotations

import itertools
from dataclasses import dataclass

from .dna import Primer
from .settings import (
    DEFAULT_PRIMER_DIMER_SETTINGS,
    PrimerDimerSettings,
)


@dataclass(slots=True, frozen=True)
class PrimerDimer:
    """Represents a primer dimer.

    Attributes:
        primer_1 (Primer): The first primer.
        primer_2 (Primer): The second primer.
        quality (float): The calculated dimer quality score.
        overlap (int): The length of the overlap region.
        p1_pos (int): The starting position of the alignment on p2 (0-indexed).
    """

    primer_1: Primer
    primer_2: Primer
    quality: float
    overlap: int
    p1_pos: int


class PrimerDimerGenerator:
    """Generates and analyses primer dimers directly from primers."""

    def __init__(
        self, settings: PrimerDimerSettings = DEFAULT_PRIMER_DIMER_SETTINGS
    ):
        """Initialize the PrimerDimerGenerator.

        Args:
            settings (PrimerDimerSettings): Settings to use for dimer
                generation.
        """
        self.settings = settings
        self.primers: list[Primer] = []
        self.primer_dimers: list[PrimerDimer] = []
        self.__analysed: bool = False

    def add(self, primer: Primer) -> None:
        """Add a primer to the generator.

        Args:
            primer (Primer): The primer to add.
        """
        self.primers.append(primer)

    def remove(self, primer: Primer) -> None:
        """Remove a primer from the generator.

        Args:
            primer (Primer): The primer to remove.
        """
        self.primers.remove(primer)

    def clear(self) -> None:
        """Clear all primers and results from the generator."""
        self.primers.clear()
        self.primer_dimers.clear()
        self.__analysed = False

    @property
    def analysed(self) -> bool:
        """Return whether the primers have been analysed."""
        return self.__analysed

    def generate_primer_dimer(self, p1: Primer, p2: Primer) -> PrimerDimer:
        """Calculate the dimer potential (quality) between two primers.

        The algorithm mimics the Amplify4 implementation. It aligns the 3' end
        of the shorter primer (p1) with different positions on the longer
        primer (p2) and scores the antiparallel overlap.

        Args:
            p1 (Primer): The first primer.
            p2 (Primer): The second primer.

        Returns:
            PrimerDimer: The best dimer alignment found.
        """
        # Ensure p1 is the shorter primer (or equal)
        if len(p1) < len(p2):
            short_p, long_p = p1, p2
        else:
            short_p, long_p = p2, p1

        n1 = len(short_p)
        n2 = len(long_p)

        seq1 = short_p.seq.upper()
        seq2 = long_p.seq.upper()

        best_quality: float = float("-inf")
        best_pos: int = 0

        # Iterate through all possible starting positions on p2
        # leftEnd in Swift implementation.
        for left_end in range(n2):
            q: float = 0.0
            # Overlap length
            # Swift: for rightEnd in leftEnd..<min(leftEnd + n1 , n2)
            # Length of segment = (min_end - left_end)

            # We iterate offset from 0 to overlap_len
            current_overlap = min(n1, n2 - left_end)

            for offset in range(current_overlap):
                # Align p1's 3' end (index n1-1) with p2 at (left_end)
                # as offset increases, we move 5' on p1 and 3' on p2

                # P1 index: n1 - 1 - offset
                c1 = seq1[n1 - 1 - offset]

                # P2 index: left_end + offset
                c2 = seq2[left_end + offset]

                # Sum weights
                q += self.settings.weights[c1, c2]

            if q >= best_quality:
                best_quality = q
                best_pos = left_end

        overlap_len = min(n1, n2 - best_pos)

        return PrimerDimer(
            primer_1=short_p,
            primer_2=long_p,
            overlap=overlap_len,
            quality=best_quality,
            p1_pos=best_pos,
        )

    def analyse_primers(self) -> None:
        """Analyse all pairs of primers for primer dimers."""
        self.primer_dimers.clear()
        for p1, p2 in itertools.combinations_with_replacement(self.primers, 2):
            res = self.generate_primer_dimer(p1, p2)
            if (
                res.quality > self.settings.threshold
                and res.overlap > self.settings.min_overlap
            ):
                self.primer_dimers.append(res)
        self.primer_dimers.sort(key=lambda x: x.quality, reverse=True)
        self.__analysed = True
