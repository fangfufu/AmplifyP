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

"""Amplicon-related classes for AmplifyP."""

from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import ClassVar

from .dna import DNA, DNADirection, DNAType, Primer
from .errors import (
    DuplicateRepliconfError,
    InvalidAmpliconRangeError,
    InvalidEndDirectionError,
    InvalidIndexOrderError,
    InvalidStartDirectionError,
    TemplateMismatchError,
)
from .repliconf import DirIdx, Repliconf


class AmpliconQuality(StrEnum):
    """Quality tier classifications for an amplicon based on q_score."""

    GOOD = "good"
    OKAY = "okay"
    MODERATE = "moderate"
    WEAK = "weak"
    VERY_WEAK = "very weak"


class AmpliconRanking(StrEnum):
    """Supported sorting modes for amplicon lists."""

    POSITION_THEN_LENGTH = "Position, then length"
    QUALITY_SCORE = "Quality score"
    POSITION_THEN_QUALITY = "Position, then quality"


@dataclass(slots=True, frozen=True)
class Amplicon:
    """A class representing a DNA amplification product (Amplicon).

    An amplicon is defined by the DNA sequence resulting from a PCR reaction,
    bounded by the forward and reverse primers. It stores metadata about
    the primers involved and the genomic location of the product.

    Attributes:
        product (DNA): The full DNA sequence of the amplicon.
        fwd_origin (Primer): The primer that initiated the forward strand
            synthesis.
        rev_origin (Primer): The primer that initiated the reverse strand
            synthesis.
        start (DirIdx): The starting index of the amplicon on the template DNA.
        end (DirIdx): The ending index of the amplicon on the template DNA.
        q_score (float): A calculated quality score for the amplicon.
        circular (bool): Whether the amplicon is circular.
    """

    _Q_SCORE_THRESHOLDS: ClassVar = (
        (300, AmpliconQuality.GOOD, " amplification"),
        (700, AmpliconQuality.OKAY, " amplification"),
        (1500, AmpliconQuality.MODERATE, " amplification"),
        (
            4000,
            AmpliconQuality.WEAK,
            " amplification — might be visible on an agarose gel",
        ),
        (
            float("inf"),
            AmpliconQuality.VERY_WEAK,
            " amplification — probably not visible on an agarose gel",
        ),
    )

    product: DNA
    fwd_origin: Primer
    rev_origin: Primer
    start: DirIdx
    end: DirIdx
    q_score: float
    circular: bool

    def __post_init__(self) -> None:
        """Validate the consistency of the amplicon indices.

        Raises:
            ValueError: If `start` is not forward, `end` is not reverse, or
                if `end` index precedes or equals `start` index on linear DNA.
        """
        if self.start.direction != DNADirection.FWD:
            raise InvalidStartDirectionError()
        if self.end.direction != DNADirection.REV:
            raise InvalidEndDirectionError()
        if (self.start.index >= self.end.index) and not self.circular:
            raise InvalidIndexOrderError()

    @staticmethod
    def quality_tier_for_score(q_score: float) -> AmpliconQuality:
        """Return the quality tier corresponding to a numeric q_score."""
        for limit, tier, _ in Amplicon._Q_SCORE_THRESHOLDS:
            if q_score < limit or limit == float("inf"):
                return tier
        return AmpliconQuality.VERY_WEAK  # pragma: no cover

    @staticmethod
    def split_sequence(
        full_seq: str, fwd_len: int, rev_len: int
    ) -> tuple[str, str, str]:
        """Split a product sequence into `(fwd_primer, middle, rev_primer)`."""
        if fwd_len < 0 or rev_len < 0:
            raise ValueError("Primer lengths must be non-negative.")
        if len(full_seq) < fwd_len + rev_len:
            return full_seq, "", ""
        fwd_part = full_seq[:fwd_len]
        if rev_len == 0:
            return fwd_part, full_seq[fwd_len:], ""
        return fwd_part, full_seq[fwd_len:-rev_len], full_seq[-rev_len:]

    @property
    def quality_tier(self) -> AmpliconQuality:
        """Return the quality tier corresponding to this amplicon's q_score."""
        return self.quality_tier_for_score(self.q_score)

    @property
    def inter_primer_length(self) -> int:
        """Return template base count between forward and reverse primers.

        May be negative if the forward and reverse primer footprints overlap.
        """
        return len(self.product) - (
            len(self.fwd_origin.seq) + len(self.rev_origin.seq)
        )

    @property
    def sequence_parts(self) -> tuple[str, str, str]:
        """Split amplicon product into `(fwd_primer, middle, rev_primer)`."""
        return self.split_sequence(
            str(self.product.seq),
            len(self.fwd_origin.seq),
            len(self.rev_origin.seq),
        )

    def q_score_report_str(self, verbose: bool = False) -> str:
        """Generate textual report based on the amplicon q-score.

        Uses a set of predefined thresholds to categorise the quality score
        into descriptive strings ranging from "good" to "very weak".

        Args:
            verbose (bool): Whether to return verbose description of the
                q-score.

        Returns:
            str: The textual report describing the amplification quality.
        """
        for limit, tier, verbose_suffix in self._Q_SCORE_THRESHOLDS:
            if self.q_score < limit or limit == float("inf"):
                suffix = verbose_suffix if verbose else ""
                circ = " (Circular)" if self.circular else ""
                return f"{tier.value}{suffix}{circ}"

        return ""  # pragma: no cover


def sort_amplicons(
    amplicons: Sequence[Amplicon],
    ranking_mode: AmpliconRanking | str = AmpliconRanking.POSITION_THEN_LENGTH,
) -> list[Amplicon]:
    """Sort amplicons according to the chosen ranking mode.

    Args:
        amplicons: Sequence of amplicons to sort.
        ranking_mode: Mode to sort by. Options:
            - 'Position, then length': Primary 5' start pos, secondary length,
              tertiary quality score.
            - 'Quality score': Primary quality score (ascending), secondary
              start pos, tertiary length.
            - 'Position, then quality': Primary 5' start pos, secondary
              quality score, tertiary length.

    Returns:
        Sorted list of amplicons.
    """
    resolved_mode = AmpliconRanking(ranking_mode)
    if resolved_mode == AmpliconRanking.POSITION_THEN_LENGTH:
        return sorted(
            amplicons,
            key=lambda a: (a.start.index, len(a.product), a.q_score),
        )
    if resolved_mode == AmpliconRanking.POSITION_THEN_QUALITY:
        return sorted(
            amplicons,
            key=lambda a: (a.start.index, a.q_score, len(a.product)),
        )
    # "Quality score": lower q_score is better quality
    return sorted(
        amplicons,
        key=lambda a: (a.q_score, a.start.index, len(a.product)),
    )


class AmpliconGenerator:
    """A generator for predicting PCR amplicons from a template.

    This class manages a list of replication configurations (Repliconfs)
    associated with a single DNA template. It simulates the PCR process by
    finding valid combinations of forward and reverse priming sites that could
    produce an amplicon.

    Attributes:
        template (DNA): The DNA sequence used as the template for amplification.
        repliconfs (list[Repliconf]): A list of `Repliconf` objects, each
            representing a primer and its binding properties on the template.
    """

    def __init__(self, template: DNA) -> None:
        """Initialise an AmpliconGenerator.

        Args:
            template (DNA): The template DNA sequence.
        """
        self.template = template
        self.repliconfs: list[Repliconf] = []
        self._repliconf_set: set[Repliconf] = set()

    def add_repliconf(self, repliconf: Repliconf) -> None:
        """Add a replication configuration to the generator.

        Args:
            repliconf (Repliconf): The configuration to add.

        Raises:
            ValueError: If the `repliconf` uses a different template than the
                generator.
            DuplicateRepliconfError: If the `repliconf` has already been added.
        """
        if self.template != repliconf.template:
            raise TemplateMismatchError()

        if repliconf in self._repliconf_set:
            raise DuplicateRepliconfError()

        self.repliconfs.append(repliconf)
        self._repliconf_set.add(repliconf)

    def remove_repliconf(self, repliconf: Repliconf) -> None:
        """Remove a replication configuration from the generator.

        Args:
            repliconf (Repliconf): The configuration to remove.

        Raises:
            ValueError: If the `repliconf` is not present in the generator.
        """
        self.repliconfs.remove(repliconf)
        self._repliconf_set.remove(repliconf)

    def get_amplicon_quality_score(
        self,
        fwd_conf: Repliconf,
        rev_conf: Repliconf,
        start: DirIdx,
        end: DirIdx,
        length: int,
    ) -> float:
        """Calculate a quality score for a potential amplicon.

        The score is derived from the length of the amplicon and the quality
        scores of the forward and reverse priming events.

        Args:
            fwd_conf (Repliconf): The configuration providing the forward
                primer.
            rev_conf (Repliconf): The configuration providing the reverse
                primer.
            start (DirIdx): The start index of the forward primer.
            end (DirIdx): The end index of the reverse primer.
            length (int): The length of the amplicon.

        Returns:
            float: The calculated quality score.
        """
        fwd_quality = fwd_conf.origin(start).quality
        rev_quality = rev_conf.origin(end).quality
        denom = (fwd_quality * rev_quality) ** 2
        if denom == 0:
            return float("inf")
        return length / denom

    def _construct_amplicon_sequence(
        self,
        fwd_conf: Repliconf,
        rev_conf: Repliconf,
        start: DirIdx,
        end: DirIdx,
    ) -> tuple[DNA | None, bool]:
        """Construct the amplicon sequence from start and end points.

        Handles both linear and circular DNA templates. For linear DNA, ensures
        start < end. For circular DNA, handles wrapping around the origin if
        end < start.

        Args:
            fwd_conf (Repliconf): Configuration for the forward primer.
            rev_conf (Repliconf): Configuration for the reverse primer.
            start (DirIdx): Start index.
            end (DirIdx): End index.

        Returns:
            tuple[DNA | None, bool]: A tuple containing the constructed DNA
                sequence (or None if the range is invalid for linear DNA) and
                a boolean indicating if the product is circular.

        Raises:
            InvalidAmpliconRangeError: If the start index equals the end
                index on a linear DNA template.
        """
        if start < end:
            # Linear DNA template branch
            seq = (
                fwd_conf.primer
                + self.template[start:end]
                + rev_conf.rev_comp_primer
            )
            return seq, False

        if self.template.type == DNAType.CIRCULAR:
            # Circular DNA handling
            seq = (
                fwd_conf.primer
                + self.template[start:]
                + self.template[:end]
                + rev_conf.rev_comp_primer
            )
            return seq, True

        if (start > end) and (self.template.type == DNAType.LINEAR):
            # Not possible on linear DNA
            return None, False

        raise InvalidAmpliconRangeError()

    def get_amplicons(self) -> list[Amplicon]:
        """Generate all possible amplicons based on added configurations.

        Triggers the search for origins in all added `Repliconf` objects (if
        not already searched). Then, collects all forward and reverse origins
        and iterates through their Cartesian product to identify valid
        amplicon combinations.

        Returns:
            list[Amplicon]: A list of all generated `Amplicon` objects sorted
                by discovery order.
        """
        amplicons: list[Amplicon] = []

        # Ensure all repliconfs are searched and collect origins
        all_fwd_origins: list[tuple[Repliconf, DirIdx]] = []
        all_rev_origins: list[tuple[Repliconf, DirIdx]] = []

        for repliconf in self.repliconfs:
            if not repliconf.searched:
                repliconf.search()
            # Collect (conf, index) pairs
            for start in repliconf.origin_db.fwd:
                all_fwd_origins.append((repliconf, start))
            for end in repliconf.origin_db.rev:
                all_rev_origins.append((repliconf, end))

        # Iterate over combinations
        is_circular = self.template.type == DNAType.CIRCULAR
        for fwd_conf, start in all_fwd_origins:
            for rev_conf, end in all_rev_origins:
                if not is_circular and start.index >= end.index:
                    continue

                seq, circular = self._construct_amplicon_sequence(
                    fwd_conf, rev_conf, start, end
                )

                if seq is None:  # pragma: no cover
                    continue

                q_score = self.get_amplicon_quality_score(
                    fwd_conf,
                    rev_conf,
                    start,
                    end,
                    (len(seq) - len(fwd_conf.primer) - len(rev_conf.primer)),
                )
                amplicons.append(
                    Amplicon(
                        seq,
                        fwd_conf.primer,
                        rev_conf.primer,
                        start,
                        end,
                        q_score,
                        circular,
                    )
                )
        return amplicons
