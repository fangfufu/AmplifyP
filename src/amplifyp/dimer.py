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
class DimerResult:
    """Result of a primer dimer analysis between two primers.

    Attributes:
        p1 (Primer): The shorter primer (or first if equal length).
        p2 (Primer): The longer primer (or second if equal length).
        overlap (int): The length of the overlap region.
        quality (float): The calculated dimer quality score.
        serious (bool): Whether the dimer is considered serious based on
            overlap and quality thresholds.
        p1_pos (int): The starting position of the alignment on p2 (0-indexed).
    """

    p1: Primer
    p2: Primer
    overlap: int
    quality: float
    serious: bool
    p1_pos: int


def calculate_dimer(
    p1: Primer,
    p2: Primer,
    settings: PrimerDimerSettings = DEFAULT_PRIMER_DIMER_SETTINGS,
) -> DimerResult:
    """Calculate the dimer potential (quality) between two primers.

    The algorithm mimics the Amplify4 implementation. It aligns the 3' end of
    the shorter primer (p1) with different positions on the longer primer (p2)
    and scores the antiparallel overlap.

    Args:
        p1 (Primer): The first primer.
        p2 (Primer): The second primer.
        settings (PrimerDimerSettings, optional): The settings object containing
            weights, overlap cutoff, and threshold cutoff. Uses
            DEFAULT_PRIMER_DIMER_SETTINGS if not provided.

    Returns:
        DimerResult: The best dimer alignment found.
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
            q += settings.weights[c1, c2]

        if q >= best_quality:
            best_quality = q
            best_pos = left_end

    overlap_len = min(n1, n2 - best_pos)

    is_serious = (best_quality > settings.threshold) and (
        overlap_len >= settings.overlap
    )

    return DimerResult(
        p1=short_p,
        p2=long_p,
        overlap=overlap_len,
        quality=best_quality,
        serious=is_serious,
        p1_pos=best_pos,
    )


def analyze_group(
    primers: list[Primer],
    settings: PrimerDimerSettings = DEFAULT_PRIMER_DIMER_SETTINGS,
) -> list[DimerResult]:
    """Analyze a group of primers for all potential dimers.

    Includes self-dimers.

    Args:
        primers (list[Primer]): List of primers to analyze.
        settings (PrimerDimerSettings, optional): Settings object.

    Returns:
        list[DimerResult]: A list of all significant (serious) dimer results,
            or all results if filtered?
            Usually users only care about "serious" dimers.
            But the User Request said "analyse a group of primer at once".
            Let's return only the serious ones to keep it clean, or maybe
            sort by quality?
            Let's return all serious ones.
    """
    results: list[DimerResult] = []

    # Analyze all combinations with replacement (includes self-dimers)
    # Swift uses `combinations` of some sort? No, the user just said analyze a
    # group.
    # Usually we compare every primer with every other primer AND itself.
    # So `itertools.combinations_with_replacement` covers (A, A), (A, B),
    # (B, B).
    # Order matters? transform(A, B) vs transform(B, A).
    # calculate_dimer normalizes by length.
    # If len(A) == len(B), calculate_dimer(A, B) swaps to (B, A) effectively
    # (assuming p2=A in my logic if p1 passed as B).
    # Wait, my logic: if len(p1) < len(p2): short=p1. else: short=p2.
    # If equal: short=p2 (second arg).
    # So calculate_dimer(A, B) -> short=B, long=A.
    # calculate_dimer(B, A) -> short=A, long=B.
    # This might produce slightly different result if the scoring is asymmetric?
    # Or if the alignment finds a different max?
    # Ideally the dimer potential is symmetric physically, but the algorithm
    # scans "short along long".
    # If equal length, scanning A along B vs B along A should be symmetric
    # IF the weights are symmetric.
    # `DEFAULT_PRIMER_DIMER_WEIGHTS` matrix...
    # `[-20, -20, -20, 30...`
    # It looks symmetric for base pairs (A-T vs T-A).
    # So we probably only need `combinations_with_replacement`.

    for p1, p2 in itertools.combinations_with_replacement(primers, 2):
        res = calculate_dimer(
            p1,
            p2,
            settings=settings,
        )
        if res.serious:
            results.append(res)

    # Sort results by quality descending
    results.sort(key=lambda x: x.quality, reverse=True)

    return results
