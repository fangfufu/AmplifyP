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

"""Amplicon-related classes for AmplifyP."""

from dataclasses import dataclass

from amplifyp.dna import DNADirection, DNAType
from amplifyp.repliconf import DirIdx

from .dna import DNA, Primer
from .errors import DuplicateRepliconfError
from .repliconf import Repliconf


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
                if `end` index precedes `start` index.
        """
        if self.start.direction != DNADirection.FWD:
            raise ValueError("Start direction must be forward.")
        if self.end.direction != DNADirection.REV:
            raise ValueError("End direction must be reverse.")
        if (self.start.index > self.end.index) and not self.circular:
            raise ValueError(
                "End index must be greater than start index for linear DNA."
            )

    def q_score_report_str(self, verbose: bool = False) -> str:
        """Generate textual report based on the amplicon q-score.

        Uses a set of predefined thresholds to categorize the quality score
        into descriptive strings ranging from "good" to "very weak".

        Args:
            verbose (bool): Whether to return verbose description of the
                q-score.

        Returns:
            str: The textual report describing the amplification quality.
        """
        thresholds = [
            (300, "good", " amplification"),
            (700, "okay", " amplification"),
            (1500, "moderate", " amplification"),
            (
                4000,
                "weak",
                " amplification — might be visible on an agarose gel",
            ),
        ]

        desc, verbose_suffix = (
            "very weak",
            (" amplification — probably not visible on an agarose gel"),
        )

        for limit, text, v_text in thresholds:
            if self.q_score < limit:
                desc, verbose_suffix = text, v_text
                break

        result = desc + (verbose_suffix if verbose else "")
        if self.circular:
            result += " (Circular)"
        return result


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
        """Initialize an AmpliconGenerator.

        Args:
            template (DNA): The template DNA sequence.
        """
        self.template = template
        self.repliconfs: list[Repliconf] = []

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
            raise ValueError(
                "The Repliconf contains a different template to the "
                "AmpliconGenerator."
            )

        if repliconf not in self.repliconfs:
            self.repliconfs.append(repliconf)
        else:
            raise DuplicateRepliconfError(
                "The Repliconf is already in the AmpliconGenerator."
            )

    def remove_repliconf(self, repliconf: Repliconf) -> None:
        """Remove a replication configuration from the generator.

        Args:
            repliconf (Repliconf): The configuration to remove.

        Raises:
            ValueError: If the `repliconf` is not present in the generator.
        """
        self.repliconfs.remove(repliconf)

    def get_amplicon_quality_score(
        self,
        fwd_conf: Repliconf,
        rev_conf: Repliconf,
        start: DirIdx,
        end: DirIdx,
        len: int,
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
            len (int): The length of the amplicon.

        Returns:
            float: The calculated quality score.
        """
        fwd_quality = fwd_conf.origin(start).quality
        rev_quality = rev_conf.origin(end).quality
        return len / (fwd_quality * rev_quality) ** 2

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
            NotImplementedError: If the start index is greater than the end
                index on a linear DNA template (which should be unreachable
                logic if called correctly).
        """
        seq = None
        circular = False
        if start < end:
            # Linear DNA template branch
            seq = (
                fwd_conf.primer
                + self.template[start:end]
                + rev_conf.primer.reverse_complement()
            )
        elif self.template.type == DNAType.CIRCULAR:
            # Circular DNA handling
            seq = (
                fwd_conf.primer
                + self.template[start:]
                + self.template[:end]
                + rev_conf.primer.reverse_complement()
            )
            circular = True
        elif (start > end) and (self.template.type == DNAType.LINEAR):
            # Not possible on linear DNA
            pass
        else:
            raise NotImplementedError(
                "Attempted to search for an amplicon with the start index "
                "bigger than the end index on a linear DNA template."
            )
        return seq, circular

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
        all_fwd_origins = []
        all_rev_origins = []

        for repliconf in self.repliconfs:
            if not repliconf.searched:
                repliconf.search()
            # Collect (conf, index) pairs
            for start in repliconf.origin_db.fwd:
                all_fwd_origins.append((repliconf, start))
            for end in repliconf.origin_db.rev:
                all_rev_origins.append((repliconf, end))

        # Iterate over combinations
        for fwd_conf, start in all_fwd_origins:
            for rev_conf, end in all_rev_origins:
                seq, circular = self._construct_amplicon_sequence(
                    fwd_conf, rev_conf, start, end
                )

                if seq is None:
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
