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

"""2D primer design module for AmplifyP."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import TYPE_CHECKING

from .amplicon import AmpliconGenerator
from .dimer import PrimerDimer, PrimerDimerGenerator
from .dna import DNA, Primer
from .repliconf import Repliconf

if TYPE_CHECKING:
    from collections.abc import Iterator

DEFAULT_PRIMER_DIMER_GENERATOR = PrimerDimerGenerator()


class FilterMetric(StrEnum):
    """Enumeration of filter metrics for 2D primer design.

    Attributes:
        MAX: Filter steps based on the maximum quality or overlap among all 4
            dimers.
        MEAN: Filter steps based on the mean quality or overlap across all 4
            dimers.
    """

    MAX = "max"
    MEAN = "mean"


@dataclass(slots=True, frozen=True)
class PrimerDimers2D:
    """Container class storing 4 primer dimers for forward/reverse primers.

    Attributes:
        fwd_fwd (PrimerDimer): Forward self-dimer.
        rev_rev (PrimerDimer): Reverse self-dimer.
        fwd_rev (PrimerDimer): Cross-dimer with forward 3' end aligned against
            reverse primer.
        rev_fwd (PrimerDimer): Cross-dimer with reverse 3' end aligned against
            forward primer.
    """

    fwd_fwd: PrimerDimer
    rev_rev: PrimerDimer
    fwd_rev: PrimerDimer
    rev_fwd: PrimerDimer

    @property
    def max_quality(self) -> float:
        """The maximum quality score among the 4 primer dimers."""
        return max(
            self.fwd_fwd.quality,
            self.rev_rev.quality,
            self.fwd_rev.quality,
            self.rev_fwd.quality,
        )

    @property
    def mean_quality(self) -> float:
        """The mean quality score across the 4 primer dimers."""
        return (
            self.fwd_fwd.quality
            + self.rev_rev.quality
            + self.fwd_rev.quality
            + self.rev_fwd.quality
        ) / 4.0

    @property
    def max_overlap(self) -> int:
        """The maximum overlap length among the 4 primer dimers."""
        return max(
            self.fwd_fwd.overlap,
            self.rev_rev.overlap,
            self.fwd_rev.overlap,
            self.rev_fwd.overlap,
        )

    @property
    def mean_overlap(self) -> float:
        """The mean overlap length across the 4 primer dimers."""
        return (
            self.fwd_fwd.overlap
            + self.rev_rev.overlap
            + self.fwd_rev.overlap
            + self.rev_fwd.overlap
        ) / 4.0


class PrimerDesigner2D:
    """Performs 2D primer design by truncating DNA sequences.

    Forward primers are truncated from the 3' end (Forward mode), while reverse
    primers are truncated from the 5' end (Reverse mode).

    Calculates 4 primer dimers (forward self-dimer, reverse self-dimer,
    forward-reverse cross-dimer, reverse-forward cross-dimer) for each
    truncation combination from initial lengths down to target minimum lengths
    `fwd_min_length` and `rev_min_length`.
    """

    __slots__ = (
        "_filter_metric",
        "_fwd_dna",
        "_fwd_min_length",
        "_generator",
        "_max_amplicon_count",
        "_max_overlap",
        "_rev_dna",
        "_rev_min_length",
        "_steps",
        "_template",
        "_threshold",
    )

    def __init__(
        self,
        fwd_dna: DNA,
        fwd_min_length: int,
        rev_dna: DNA,
        rev_min_length: int,
        generator: PrimerDimerGenerator = DEFAULT_PRIMER_DIMER_GENERATOR,
        threshold: float | None = None,
        max_overlap: int | None = None,
        filter_metric: FilterMetric = FilterMetric.MAX,
        template: DNA | None = None,
        max_amplicon_count: int | None = None,
    ) -> None:
        """Initialises a new PrimerDesigner2D object and runs the analysis.

        Args:
            fwd_dna (DNA): The input forward DNA object.
            fwd_min_length (int): Minimum sequence length for forward primer.
            rev_dna (DNA): The input reverse DNA object.
            rev_min_length (int): Minimum sequence length for reverse primer.
            generator (PrimerDimerGenerator, optional): Custom primer dimer
                generator to use for analysis. Defaults to
                `DEFAULT_PRIMER_DIMER_GENERATOR`.
            threshold (float | None, optional): Upper bound for quality score
                filter. Defaults to None.
            max_overlap (int | None, optional): Upper bound for overlap length
                filter. Defaults to None.
            filter_metric (FilterMetric, optional): Metric to use for filtering
                (MAX or MEAN). Defaults to `FilterMetric.MAX`.
            template (DNA | None, optional): Optional template DNA sequence
                used to calculate predicted amplicons. Defaults to None.
            max_amplicon_count (int | None, optional): Upper bound for number of
                predicted amplicons filter. Defaults to None.

        Raises:
            ValueError: If minimum lengths are non-positive or exceed sequence
                lengths, if max_amplicon_count is specified without template, or
                if max_amplicon_count is negative.
        """
        if fwd_min_length <= 0:
            raise ValueError("Forward target length n must be greater than 0.")
        if len(fwd_dna.seq_upper) < fwd_min_length:
            msg = (
                f"Forward target length n ({fwd_min_length}) cannot exceed "
                f"initial sequence length ({len(fwd_dna.seq_upper)})."
            )
            raise ValueError(msg)

        if rev_min_length <= 0:
            raise ValueError("Reverse target length n must be greater than 0.")
        if len(rev_dna.seq_upper) < rev_min_length:
            msg = (
                f"Reverse target length n ({rev_min_length}) cannot exceed "
                f"initial sequence length ({len(rev_dna.seq_upper)})."
            )
            raise ValueError(msg)

        if max_amplicon_count is not None:
            if template is None:
                msg = (
                    "Template DNA must be provided when max_amplicon_count "
                    "filter is specified."
                )
                raise ValueError(msg)
            if max_amplicon_count < 0:
                raise ValueError("max_amplicon_count must be non-negative.")

        self._fwd_dna: DNA = fwd_dna
        self._fwd_min_length: int = fwd_min_length
        self._rev_dna: DNA = rev_dna
        self._rev_min_length: int = rev_min_length
        self._generator: PrimerDimerGenerator = generator
        self._threshold: float | None = threshold
        self._max_overlap: int | None = max_overlap
        self._filter_metric: FilterMetric = FilterMetric(filter_metric)
        self._template: DNA | None = template
        self._max_amplicon_count: int | None = max_amplicon_count
        self._steps: list[PrimerDimers2D] = []

        self._analyse()

    @property
    def fwd_dna(self) -> DNA:
        """The input forward DNA object."""
        return self._fwd_dna

    @property
    def fwd_min_length(self) -> int:
        """The minimum sequence length for forward primer."""
        return self._fwd_min_length

    @property
    def rev_dna(self) -> DNA:
        """The input reverse DNA object."""
        return self._rev_dna

    @property
    def rev_min_length(self) -> int:
        """The minimum sequence length for reverse primer."""
        return self._rev_min_length

    @property
    def generator(self) -> PrimerDimerGenerator:
        """The primer dimer generator used for analysis."""
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
    def filter_metric(self) -> FilterMetric:
        """The filter metric used for evaluation (MAX or MEAN)."""
        return self._filter_metric

    @property
    def template(self) -> DNA | None:
        """The optional template DNA object, if set."""
        return self._template

    @property
    def max_amplicon_count(self) -> int | None:
        """The maximum amplicon count cutoff filter, if set."""
        return self._max_amplicon_count

    @property
    def all_steps(self) -> tuple[PrimerDimers2D, ...]:
        """The stored 2D steps from sequence truncation."""
        return tuple(self._steps)

    @property
    def best_score(self) -> tuple[int, float]:
        """Return the (index, quality) pair of the optimal 2D step.

        Returns:
            tuple[int, float]: A (step_index, quality_score) tuple for the step
                with the lowest quality score based on active `filter_metric`.

        Raises:
            RuntimeError: If no results have been generated.
        """
        if not self._steps:
            raise RuntimeError("No analysis steps recorded.")
        return self.quality_score(sorted=True)[0]

    def quality_score(
        self, sorted: bool = False
    ) -> tuple[tuple[int, float], ...]:
        """Return (index, quality) pairs for all 2D step results.

        Args:
            sorted (bool, optional): If True, sort pairs by quality score in
                ascending order (best quality first). Defaults to False.

        Returns:
            tuple[tuple[int, float], ...]: Sequence of (step_index, quality)
                tuples.
        """
        if self._filter_metric == FilterMetric.MAX:
            pairs = [(i, s.max_quality) for i, s in enumerate(self._steps)]
        else:
            pairs = [(i, s.mean_quality) for i, s in enumerate(self._steps)]
        if sorted:
            pairs.sort(key=lambda item: item[1])
        return tuple(pairs)

    def __len__(self) -> int:
        """Return the number of stored design steps."""
        return len(self._steps)

    def __getitem__(self, index: int) -> PrimerDimers2D:
        """Get the 2D step at 0-indexed position `index`.

        Args:
            index (int): The 0-indexed step position.

        Returns:
            PrimerDimers2D: The step object at the given index.
        """
        return self._steps[index]

    def __iter__(self) -> Iterator[PrimerDimers2D]:
        """Iterate over stored 2D steps."""
        return iter(self._steps)

    def __contains__(self, item: object) -> bool:
        """Check if a PrimerDimers2D is contained in the results."""
        return item in self._steps

    def __repr__(self) -> str:
        """Return an unambiguous string representation of PrimerDesigner2D."""
        return (
            f"PrimerDesigner2D(fwd_dna={self._fwd_dna!r}, "
            f"fwd_min_length={self._fwd_min_length}, "
            f"rev_dna={self._rev_dna!r}, "
            f"rev_min_length={self._rev_min_length})"
        )

    def __str__(self) -> str:
        """Return a user-friendly string representation of PrimerDesigner2D."""
        return (
            f"PrimerDesigner2D({len(self._steps)} steps, "
            f"fwd_min_length={self._fwd_min_length}, "
            f"rev_min_length={self._rev_min_length})"
        )

    def get_step(self, index: int) -> PrimerDimers2D:
        """Query stored 2D step by 0-indexed position `index`.

        Args:
            index (int): The 0-indexed position.

        Returns:
            PrimerDimers2D: The 2D step corresponding to `index`.
        """
        return self[index]

    def _analyse(self) -> None:
        """Perform 2D truncation analysis on forward and reverse sequences.

        Forward sequence is truncated from 3' end, while reverse sequence is
        truncated from 5' end.
        """
        self._steps.clear()

        # Forward sequence truncation (3' end chopping)
        fwd_seqs: list[str] = []
        seq = self._fwd_dna.seq_upper
        while len(seq) >= self._fwd_min_length:
            fwd_seqs.append(seq)
            if len(seq) == self._fwd_min_length:
                break
            seq = seq[:-1]

        # Reverse sequence truncation (5' end chopping)
        rev_seqs: list[str] = []
        seq = self._rev_dna.seq_upper
        while len(seq) >= self._rev_min_length:
            rev_seqs.append(seq)
            if len(seq) == self._rev_min_length:
                break
            seq = seq[1:]

        for fwd_seq in fwd_seqs:
            fwd_p = Primer(fwd_seq)
            for rev_seq in rev_seqs:
                rev_p = Primer(rev_seq)

                if self._template is not None:
                    amp_gen = AmpliconGenerator(self._template)
                    fwd_conf = Repliconf(self._template, fwd_p)
                    amp_gen.add_repliconf(fwd_conf)
                    if fwd_seq != rev_seq:
                        rev_conf = Repliconf(self._template, rev_p)
                        amp_gen.add_repliconf(rev_conf)
                    amplicon_count = len(amp_gen.get_amplicons())
                    if (
                        self._max_amplicon_count is not None
                        and amplicon_count > self._max_amplicon_count
                    ):
                        continue

                d_ff = self._generator.generate_primer_dimer(
                    fwd_p, fwd_p, reorder=False
                )
                d_rr = self._generator.generate_primer_dimer(
                    rev_p, rev_p, reorder=False
                )
                d_fr = self._generator.generate_primer_dimer(
                    fwd_p, rev_p, reorder=False
                )
                d_rf = self._generator.generate_primer_dimer(
                    rev_p, fwd_p, reorder=False
                )

                step = PrimerDimers2D(
                    fwd_fwd=d_ff,
                    rev_rev=d_rr,
                    fwd_rev=d_fr,
                    rev_fwd=d_rf,
                )

                if self._filter_metric == FilterMetric.MAX:
                    q_val = step.max_quality
                    o_val = float(step.max_overlap)
                else:
                    q_val = step.mean_quality
                    o_val = step.mean_overlap

                if self._threshold is not None and q_val > self._threshold:
                    continue
                if self._max_overlap is not None and o_val > self._max_overlap:
                    continue

                self._steps.append(step)
