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

"""Primer dimer analysis module for AmplifyP."""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

from .dna import Primer
from .settings import (
    GLOBAL_PRIMER_DIMER_SETTINGS,
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
    settings: PrimerDimerSettings = field(
        default_factory=lambda: GLOBAL_PRIMER_DIMER_SETTINGS
    )

    @property
    def binding_strength_str(self) -> str:
        """The binding strength of each base-pair as a string."""
        p1 = self.primer_1
        p2 = self.primer_2
        if len(p1) < len(p2):
            short_p, long_p = p1, p2
        else:
            short_p, long_p = p2, p1

        seq1 = short_p.seq.upper()
        seq2 = long_p.seq.upper()
        n1 = len(short_p)

        bonds: list[str] = []
        for offset in range(self.overlap):
            c1 = seq1[n1 - 1 - offset]
            c2 = seq2[self.p1_pos + offset]
            try:
                score = self.settings.weights[c1, c2]
            except KeyError:
                score = 0.0

            if score < 0:
                symbol = " "
            elif score < self.settings.symbol_threshold:
                symbol = ":"
            else:
                symbol = "|"
            bonds.append(symbol)
        return "".join(bonds)


class PrimerDimerGenerator:
    """Generates and analyses primer dimers directly from primers."""

    def __init__(
        self, settings: PrimerDimerSettings = GLOBAL_PRIMER_DIMER_SETTINGS
    ):
        """Initialise the PrimerDimerGenerator.

        Args:
            settings (PrimerDimerSettings): Settings to use for dimer
                generation.
        """
        self.settings = settings
        self.primers: list[Primer] = []
        self.primer_dimers: list[PrimerDimer] = []
        self.__analysed: bool = False

    def add_primer(self, primer: Primer) -> None:
        """Add a primer to the generator.

        Args:
            primer (Primer): The primer to add.
        """
        self.primers.append(primer)

    def remove_primer(self, primer: Primer) -> None:
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

    def _calculate_dimer_stats(
        self, s1: str, s2: str, n1: int, n2: int
    ) -> tuple[float, int, int]:
        return calculate_dimer_stats(s1, s2, n1, n2, settings=self.settings)

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
        return generate_primer_dimer(p1, p2, settings=self.settings)

    def analyse_primers(self) -> None:
        """Analyse all pairs of primers for primer dimers.

        This method iterates through all unique pairs of primers added to the
        generator, calculates their dimer potential, and stores those that
        meet the specified quality and overlap thresholds.
        """
        self.primer_dimers.clear()

        primer_props = [(p, len(p), p.seq_upper) for p in self.primers]
        memo: dict[tuple[str, str], tuple[float, int, int]] = {}

        threshold = self.settings.threshold
        min_overlap = self.settings.min_overlap

        for (p1, l1, seq1), (
            p2,
            l2,
            seq2,
        ) in itertools.combinations_with_replacement(primer_props, 2):
            if l1 < l2:
                short_p, long_p, n1, n2 = p1, p2, l1, l2
                s1, s2 = seq1, seq2
            else:
                short_p, long_p, n1, n2 = p2, p1, l2, l1
                s1, s2 = seq2, seq1

            memo_key = (s1, s2)
            if memo_key in memo:
                best_quality, best_pos, overlap_len = memo[memo_key]
            else:
                best_quality, best_pos, overlap_len = calculate_dimer_stats(
                    s1, s2, n1, n2, settings=self.settings
                )
                memo[memo_key] = (best_quality, best_pos, overlap_len)

            if best_quality > threshold and overlap_len > min_overlap:
                res = PrimerDimer(
                    primer_1=short_p,
                    primer_2=long_p,
                    overlap=overlap_len,
                    quality=best_quality,
                    p1_pos=best_pos,
                    settings=self.settings,
                )
                self.primer_dimers.append(res)

        self.primer_dimers.sort(key=lambda x: x.quality, reverse=True)
        self.__analysed = True


def calculate_dimer_stats(
    s1: str,
    s2: str,
    n1: int,
    n2: int,
    settings: PrimerDimerSettings = GLOBAL_PRIMER_DIMER_SETTINGS,
) -> tuple[float, int, int]:
    """Calculate the best dimer alignment quality, position, and overlap length.

    Args:
        s1 (str): First sequence (shorter or equal).
        s2 (str): Second sequence (longer or equal).
        n1 (int): Length of first sequence.
        n2 (int): Length of second sequence.
        settings (PrimerDimerSettings): Settings to use.

    Returns:
        tuple[float, int, int]: Best quality score, best position, and overlap
            length.
    """
    best_quality: float = float("-inf")
    best_pos: int = 0
    weights = settings.weights

    for left_end in range(n2):
        q: float = 0.0
        current_overlap = min(n1, n2 - left_end)

        for offset in range(current_overlap):
            c1 = s1[n1 - 1 - offset]
            c2 = s2[left_end + offset]

            try:
                q += weights[c1, c2]
            except KeyError:
                pass

        if q >= best_quality:
            best_quality = q
            best_pos = left_end

    overlap_len = min(n1, n2 - best_pos)
    return best_quality, best_pos, overlap_len


def generate_primer_dimer(
    p1: Primer,
    p2: Primer,
    settings: PrimerDimerSettings = GLOBAL_PRIMER_DIMER_SETTINGS,
) -> PrimerDimer:
    """Calculate the dimer potential (quality) between two primers.

    The algorithm mimics the Amplify4 implementation. It aligns the 3' end
    of the shorter primer (p1) with different positions on the longer
    primer (p2) and scores the antiparallel overlap.

    Args:
        p1 (Primer): The first primer.
        p2 (Primer): The second primer.
        settings (PrimerDimerSettings): Settings to use.

    Returns:
        PrimerDimer: The best dimer alignment found.
    """
    if len(p1) < len(p2):
        short_p, long_p = p1, p2
    else:
        short_p, long_p = p2, p1

    n1 = len(short_p)
    n2 = len(long_p)

    seq1 = short_p.seq_upper
    seq2 = long_p.seq_upper

    best_quality, best_pos, overlap_len = calculate_dimer_stats(
        seq1, seq2, n1, n2, settings=settings
    )

    return PrimerDimer(
        primer_1=short_p,
        primer_2=long_p,
        overlap=overlap_len,
        quality=best_quality,
        p1_pos=best_pos,
        settings=settings,
    )
