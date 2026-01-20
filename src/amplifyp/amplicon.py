"""Amplicon-related classes for AmplifyP."""

from dataclasses import dataclass

from .dna import DNA, Primer
from .errors import DuplicateRepliconfError
from .repliconf import Repliconf


@dataclass(slots=True)
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
    start: int
    end: int

    def __post_init__(self) -> None:
        """Post-initialization checks for the Amplicon."""


class AmpliconGenerator:
    """A class for generating amplicons from a template DNA sequence.

    This class takes a template DNA and a list of replication configurations
    (Repliconf) to generate amplicons. It ensures that all replication
    configurations use the same template DNA.

    Attributes:
        template (DNA): The template DNA sequence.
        repliconfs (List[Repliconf]): A list of replication configurations.
    """

    # pylint: disable=too-few-public-methods

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
            for start in fwd_conf.amplicon_start:
                for rev_conf in self.repliconfs:
                    for end in rev_conf.amplicon_end:
                        if start < end:
                            # Generate amplicon sequence from template slice.
                            # Python slicing handles ends gracefully.
                            seq = self.template[start:end]
                            amplicons.append(
                                Amplicon(
                                    seq, fwd_conf.primer, rev_conf.primer, start, end
                                )
                            )
        return amplicons
