# -*- coding: utf-8 -*-
"""Amplicon-related classes for AmplifyP."""

from dataclasses import dataclass
from typing import List

from .dna import DNA, Primer
from .repliconf import Repliconf


@dataclass
class Amplicon:
    """
    A class representing an amplicon.

    An amplicon is a piece of DNA or RNA that is the source and/or product of
    amplification or replication events. It is defined by the sequence that has
    been amplified, and the forward and reverse replication origins.

    Attributes:
        sequence (DNA): The DNA sequence of the amplicon.
        fwd_origin (Primer): The forward replication origin.
        rev_origin (Primer): The reverse replication origin.
    """

    sequence: DNA
    fwd_origin: Primer
    rev_origin: Primer

    def __post_init__(self) -> None:
        """Post-initialization checks for the Amplicon."""


class AmpliconGenerator:
    """
    A class for generating amplicons from a template DNA sequence.

    This class takes a template DNA and a list of replication configurations
    (Repliconf) to generate amplicons. It ensures that all replication
    configurations use the same template DNA.

    Attributes:
        template (DNA): The template DNA sequence.
        repliconfs (List[Repliconf]): A list of replication configurations.
    """

    # pylint: disable=too-few-public-methods

    def __init__(self, template: DNA) -> None:
        """
        Constructs an AmpliconGenerator object.

        Args:
            template (DNA): The template DNA sequence.
        """
        self.template = template
        self.repliconfs: List[Repliconf] = []

    def add(self, repliconf: Repliconf) -> None:
        """
        Adds a replication configuration to the AmpliconGenerator.

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

        if not self.repliconfs:
            self.repliconfs.append(repliconf)
