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

from .dimer import PrimerDimer, PrimerDimerGenerator
from .dna import DNA, DNADirection, DNAType, Primer
from .repliconf import Repliconf
from .settings import GLOBAL_REPLICATION_SETTINGS, ReplicationSettings

if TYPE_CHECKING:
    import threading
    from collections.abc import Callable, Iterator

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
    amplicon_count: int | None = None

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
        "_aborted",
        "_cancel_event",
        "_filter_metric",
        "_fwd_dna",
        "_fwd_min_length",
        "_generator",
        "_max_amplicon_count",
        "_max_overlap",
        "_on_progress",
        "_replication_settings",
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
        on_progress: Callable[[int, int], None] | None = None,
        cancel_event: threading.Event | None = None,
        replication_settings: ReplicationSettings = GLOBAL_REPLICATION_SETTINGS,
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
            on_progress (Callable[[int, int], None] | None, optional): Optional
                callback invoked after each (fwd, rev) pair is processed.
                 Receives ``(done, total)`` where *done* counts completed pairs
                 and *total* is the full combination count. Defaults to None.
            cancel_event (threading.Event | None, optional): Optional event
                used to abort the analysis. When set, the analysis stops
                after the current (fwd, rev) pair is processed and the steps
                analysed so far are retained. Defaults to None.
            replication_settings (ReplicationSettings, optional): Custom
                replication settings used for template amplicon searches.
                Defaults to `GLOBAL_REPLICATION_SETTINGS`.

        Raises:
            ValueError: If minimum lengths are non-positive or exceed sequence
                lengths, if max_amplicon_count is specified without template, or
                if max_amplicon_count is negative.
        """
        fwd_seqs = fwd_dna.generate_truncations(
            fwd_min_length,
            DNADirection.FWD,
            prefix="Forward target length n",
        )
        rev_seqs = rev_dna.generate_truncations(
            rev_min_length,
            DNADirection.REV,
            prefix="Reverse target length n",
        )

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
        self._on_progress: Callable[[int, int], None] | None = on_progress
        self._cancel_event: threading.Event | None = cancel_event
        self._replication_settings: ReplicationSettings = replication_settings
        self._aborted: bool = False
        self._steps: list[PrimerDimers2D] = []

        self._analyse((fwd_seqs, rev_seqs))

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
    def replication_settings(self) -> ReplicationSettings:
        """The replication settings used for template amplicon searches."""
        return self._replication_settings

    @property
    def aborted(self) -> bool:
        """Whether the analysis was stopped before completion."""
        return self._aborted

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

    @staticmethod
    def _count_linear_amplicons(starts: list[int], ends: list[int]) -> int:
        """Count ordered (start, end) pairs with start < end."""
        n_ends = len(ends)
        if not starts or not n_ends:
            return 0
        count = 0
        j = 0
        for s in starts:
            while j < n_ends and ends[j] <= s:
                j += 1
            if j == n_ends:
                break
            count += n_ends - j
        return count

    def _fails_max_dimer(self, dimer: PrimerDimer) -> bool:
        """Return True if a dimer alone exceeds threshold or max_overlap."""
        return (
            self._threshold is not None and dimer.quality > self._threshold
        ) or (
            self._max_overlap is not None and dimer.overlap > self._max_overlap
        )

    def _get_origin_indices(
        self,
        template: DNA,
        primer: Primer,
        cache: dict[str, tuple[list[int], list[int], int]],
        is_circular: bool,
    ) -> tuple[list[int], list[int], int]:
        """Return sorted fwd indices, rev indices, and self-amplicon count."""
        key = primer.seq_upper
        cached = cache.get(key)
        if cached is not None:
            return cached
        conf = Repliconf(
            template,
            primer,
            settings=self._replication_settings,
        )
        conf.search()
        fwd_idx = [o.index for o in conf.origin_db.fwd]
        rev_idx = [o.index for o in conf.origin_db.rev]
        self_count = (
            len(fwd_idx) * len(rev_idx)
            if is_circular
            else self._count_linear_amplicons(fwd_idx, rev_idx)
        )
        res = (fwd_idx, rev_idx, self_count)
        cache[key] = res
        return res

    def _evaluate_pair(
        self,
        template: DNA | None,
        fwd_p: Primer,
        d_ff: PrimerDimer,
        fwd_origins: tuple[list[int], list[int], int] | None,
        rev_p: Primer,
        d_rr: PrimerDimer,
        origin_cache: dict[str, tuple[list[int], list[int], int]],
        is_max_metric: bool,
        is_circular: bool,
    ) -> PrimerDimers2D | None:
        """Evaluate one (fwd, rev) truncation pair against all filters."""
        amplicon_count: int | None = None
        if template is not None and fwd_origins is not None:
            fwd_f, fwd_r, fwd_self = fwd_origins
            if fwd_p.seq_upper == rev_p.seq_upper:
                amplicon_count = fwd_self
            else:
                rev_f, rev_r, rev_self = self._get_origin_indices(
                    template, rev_p, origin_cache, is_circular
                )
                if is_circular:
                    amplicon_count = (len(fwd_f) + len(rev_f)) * (
                        len(fwd_r) + len(rev_r)
                    )
                else:
                    amplicon_count = (
                        fwd_self
                        + rev_self
                        + self._count_linear_amplicons(fwd_f, rev_r)
                        + self._count_linear_amplicons(rev_f, fwd_r)
                    )
            if amplicon_count < 1:
                return None
            if (
                self._max_amplicon_count is not None
                and amplicon_count > self._max_amplicon_count
            ):
                return None

        gen = self._generator.generate_primer_dimer
        d_fr = gen(fwd_p, rev_p, reorder=False)
        if is_max_metric and self._fails_max_dimer(d_fr):
            return None

        d_rf = gen(rev_p, fwd_p, reorder=False)
        if is_max_metric and self._fails_max_dimer(d_rf):
            return None

        step = PrimerDimers2D(
            fwd_fwd=d_ff,
            rev_rev=d_rr,
            fwd_rev=d_fr,
            rev_fwd=d_rf,
            amplicon_count=amplicon_count,
        )

        if not is_max_metric and (
            (
                self._threshold is not None
                and step.mean_quality > self._threshold
            )
            or (
                self._max_overlap is not None
                and step.mean_overlap > self._max_overlap
            )
        ):
            return None

        return step

    def _analyse(
        self,
        initial_truncations: tuple[list[str], list[str]] | None = None,
    ) -> None:
        """Perform 2D truncation analysis on forward and reverse sequences.

        Forward sequence is truncated from 3' end, while reverse sequence is
        truncated from 5' end.
        """
        self._steps.clear()
        if initial_truncations is not None:
            fwd_seqs, rev_seqs = initial_truncations
        else:
            fwd_seqs = self._fwd_dna.generate_truncations(
                self._fwd_min_length, DNADirection.FWD
            )
            rev_seqs = self._rev_dna.generate_truncations(
                self._rev_min_length, DNADirection.REV
            )

        total = len(fwd_seqs) * len(rev_seqs)
        done = 0

        is_max_metric = self._filter_metric == FilterMetric.MAX
        template = self._template
        is_circular = template is not None and template.type == DNAType.CIRCULAR
        gen = self._generator.generate_primer_dimer
        origin_cache: dict[str, tuple[list[int], list[int], int]] = {}
        rev_dimer_cache: list[tuple[Primer, PrimerDimer, bool] | None] = [
            None
        ] * len(rev_seqs)

        for fwd_seq in fwd_seqs:
            if self._cancel_event is not None and self._cancel_event.is_set():
                self._aborted = True
                break
            fwd_p = Primer(fwd_seq)
            d_ff = gen(fwd_p, fwd_p, reorder=False)
            fwd_fails_max = is_max_metric and self._fails_max_dimer(d_ff)

            fwd_origins: tuple[list[int], list[int], int] | None = None
            if template is not None and not fwd_fails_max:
                fwd_origins = self._get_origin_indices(
                    template, fwd_p, origin_cache, is_circular
                )

            for r_idx, rev_seq in enumerate(rev_seqs):
                if (
                    self._cancel_event is not None
                    and self._cancel_event.is_set()
                ):
                    self._aborted = True
                    break

                if not fwd_fails_max:
                    rev_entry = rev_dimer_cache[r_idx]
                    if rev_entry is None:
                        rev_p = Primer(rev_seq)
                        d_rr = gen(rev_p, rev_p, reorder=False)
                        rev_fails_max = is_max_metric and self._fails_max_dimer(
                            d_rr
                        )
                        rev_entry = (rev_p, d_rr, rev_fails_max)
                        rev_dimer_cache[r_idx] = rev_entry
                    else:
                        rev_p, d_rr, rev_fails_max = rev_entry

                    if not rev_fails_max:
                        step = self._evaluate_pair(
                            template,
                            fwd_p,
                            d_ff,
                            fwd_origins,
                            rev_p,
                            d_rr,
                            origin_cache,
                            is_max_metric,
                            is_circular,
                        )
                        if step is not None:
                            self._steps.append(step)

                done += 1
                if self._on_progress is not None:
                    self._on_progress(done, total)
