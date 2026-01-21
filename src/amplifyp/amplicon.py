"""Amplicon-related classes for AmplifyP."""

from dataclasses import dataclass

from amplifyp.dna import DNADirection
from amplifyp.repliconf import DirIdx

from .dna import DNA, Primer
from .errors import DuplicateRepliconfError
from .repliconf import Repliconf


@dataclass(slots=True, frozen=True)
class Amplicon:
    """A class representing an amplicon.

    An amplicon is a piece of DNA that is the product of the amplification
    events. It is defined by the sequence that has been amplified, and the
    forward and reverse replication origins.

    Attributes:
        product (DNA): The amplification product.
        fwd_origin (Primer): The forward replication origin.
        rev_origin (Primer): The reverse replication origin.
        start (int): The start index of the amplicon relative to the template.
        end (int): The end index of the amplicon relative to the template.
    """

    product: DNA
    fwd_origin: Primer
    rev_origin: Primer
    start: DirIdx
    end: DirIdx
    q_score: float

    def __post_init__(self) -> None:
        """Validate the amplicon indices."""
        if self.start.direction != DNADirection.FWD:
            raise ValueError("Start direction must be forward.")
        if self.end.direction != DNADirection.REV:
            raise ValueError("End direction must be reverse.")
        if self.end < self.start:
            raise ValueError("End index must be after start index.")


class AmpliconGenerator:
    """A class for generating amplicons from a template DNA sequence.

    This class takes a template DNA and a list of replication configurations
    (Repliconf) to generate amplicons. It ensures that all replication
    configurations use the same template DNA.

    Attributes:
        template (DNA): The template DNA sequence.
        repliconfs (List[Repliconf]): A list of replication configurations.
    """

    def __init__(self, template: DNA) -> None:
        """Constructs an AmpliconGenerator object.

        Args:
            template (DNA): The template DNA sequence.
        """
        self.template = template
        self.repliconfs: list[Repliconf] = []

    def add(self, repliconf: Repliconf) -> None:
        """Adds a replication configuration to the AmpliconGenerator.

        Args:
            repliconf (Repliconf): The replication configuration to add.

        Raises:
            ValueError: If the Repliconf contains a different template than
                        the one in the AmpliconGenerator.
        """
        if self.template != repliconf.template:
            raise ValueError(
                "The Repliconf contains a different template to the AmpliconGenerator."
            )

        if repliconf not in self.repliconfs:
            self.repliconfs.append(repliconf)
        else:
            raise DuplicateRepliconfError(
                "The Repliconf is already in the AmpliconGenerator."
            )

    def remove(self, repliconf: Repliconf) -> None:
        """Removes a replication configuration from the AmpliconGenerator.

        Args:
            repliconf (Repliconf): The replication configuration to remove.

        Raises:
            ValueError: If the Repliconf is not in the AmpliconGenerator.
        """
        self.repliconfs.remove(repliconf)

    def get_amplicon_quality_score(
        self,
        fwd_conf: Repliconf,
        rev_conf: Repliconf,
        start_db_idx: int,
        end_db_idx: int,
        start: DirIdx,
        end: DirIdx,
    ) -> float:
        """Get the quality score of the amplicons."""
        fwd_quality = fwd_conf.origin_from_db(DNADirection.FWD, start_db_idx).quality
        rev_quality = rev_conf.origin_from_db(DNADirection.REV, end_db_idx).quality
        return (int(end - start)) / (fwd_quality * rev_quality) ** 2

    def get_amplicons(self) -> list[Amplicon]:
        """Get amplicons from the added replication configurations.

        Returns:
            List[Amplicon]: A list of generated Amplicons.
        """
        amplicons: list[Amplicon] = []

        for repliconf in self.repliconfs:
            if not repliconf.searched:
                repliconf.search()

        for fwd_conf in self.repliconfs:
            for start_db_idx, start in enumerate(fwd_conf.target_start):
                for rev_conf in self.repliconfs:
                    for end_db_idx, end in enumerate(rev_conf.target_end):
                        if start < end:
                            # Generate amplicon sequence from template slice.
                            # Python slicing handles ends gracefully.
                            seq = (
                                fwd_conf.primer
                                + self.template[start:end]
                                + rev_conf.primer.reverse_complement()
                            )
                            q_score = self.get_amplicon_quality_score(
                                fwd_conf, rev_conf, start_db_idx, end_db_idx, start, end
                            )
                            amplicons.append(
                                Amplicon(
                                    seq,
                                    fwd_conf.primer,
                                    rev_conf.primer,
                                    start,
                                    end,
                                    q_score,
                                )
                            )
        return amplicons
