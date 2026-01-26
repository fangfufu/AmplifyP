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

"""This module contains the PCR class, which represents a PCR reaction."""

from amplifyp.amplicon import Amplicon, AmpliconGenerator
from amplifyp.dna import DNA, Primer
from amplifyp.errors import DuplicatedPrimerError, PrimerNotFoundError
from amplifyp.repliconf import Repliconf
from amplifyp.settings import GLOBAL_REPLICATION_SETTINGS, ReplicationSettings


class PCR:
    """A class representing a PCR reaction.

    Attributes:
        template (DNA): The template DNA.
        primers (list[Primer]): The primers used in the reaction.
        repliconfs (list[Replicon]): The replication configurations.
        amplicon_generator (AmpliconGenerator): The amplicon generator.
        amplicons (list[Amplicon]): The predicted amplicons.
    """

    def __init__(
        self,
        template: DNA,
        settings: ReplicationSettings = GLOBAL_REPLICATION_SETTINGS,
    ):
        """Initialize a PCR reaction."""
        self.template = template
        self.settings = settings
        self.__primers: list[Primer] = []
        self.amplicon_generator = AmpliconGenerator(self.template)
        self.__amplicons: list[Amplicon] = []

    def add_primer(self, primer: Primer) -> None:
        """Add a primer to the PCR reaction.

        Args:
            primer (Primer): The primer to add.

        Raises:
            DuplicatedPrimerError: If the primer is already added.
        """
        if primer in self.__primers:
            raise DuplicatedPrimerError(f"Primer {primer} already added")
        self.__primers.append(primer)
        repliconf = Repliconf(self.template, primer)
        self.amplicon_generator.add_repliconf(repliconf)

    def remove_primer(self, primer: Primer) -> None:
        """Remove a primer from the PCR reaction.

        Args:
            primer (Primer): The primer to remove.

        Raises:
            PrimerNotFoundError: If the primer is not found.
        """
        if primer not in self.__primers:
            raise PrimerNotFoundError(f"Primer {primer} not found")
        self.__primers.remove(primer)
        repliconf = Repliconf(self.template, primer)
        self.amplicon_generator.remove_repliconf(repliconf)

    def add_primers(self, primers: list[Primer]) -> None:
        """Add multiple primers to the PCR reaction.

        Args:
            primers (list[Primer]): The primers to add.
        """
        for primer in primers:
            self.add_primer(primer)

    @property
    def primers(self) -> list[Primer]:
        """Get the primers used in the PCR reaction."""
        return self.__primers.copy()

    def predict_amplicons(self) -> int:
        """Predict the amplicons for the PCR reaction."""
        self.__amplicons = self.amplicon_generator.get_amplicons()
        return len(self.__amplicons)

    @property
    def amplicons(self) -> list[Amplicon]:
        """Get the predicted amplicons."""
        return self.__amplicons.copy()
