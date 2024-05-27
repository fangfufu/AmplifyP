# -*- coding: utf-8 -*-
"""Amplify P - amplicon related."""

from dataclasses import dataclass
from typing import List

from .dna import DNA, Primer
from .repliconf import Repliconf


@dataclass
class Amplicon:
    """
    A class representing the amplicons

    An amplicon is defined by the sequence that had been amplified. It starts
    with the forward replication origin, it ends with the reverse replication
    origin.
    """

    sequence: DNA
    fwd_origin: Primer
    rev_origin: Primer

    def __post_init__(self) -> None:
        pass


class AmpliconGenerator:
    """
    A class for generating the amplicons
    """

    # pylint: disable=too-few-public-methods

    def __init__(self, template: DNA) -> None:
        """Construct an AmpliconGenerator"""
        self.template = template
        self.repliconfs: List[Repliconf] = []

    def add(self, repliconf: Repliconf) -> None:
        """Add a Repliconf to the AmpliconGenerator"""
        if self.template != repliconf.template:
            raise ValueError(
                "The Repliconf contains a different template to the AmpliconGenerator."
            )

        if not self.repliconfs:
            self.repliconfs.append(repliconf)
