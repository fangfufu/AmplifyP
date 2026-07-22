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

"""PCR reaction simulation module.

This module provides the `PCR` class, which represents a complete Polymerase
Chain Reaction (PCR) simulation. It manages templates, primers, and the
internal amplicon generator to predict all possible amplification products.
"""

from collections.abc import Iterable

from amplifyp.amplicon import Amplicon, AmpliconGenerator
from amplifyp.dna import DNA, Primer
from amplifyp.errors import (
    DuplicatedNameError,
    DuplicatedSequenceError,
    PrimerNotFoundError,
)
from amplifyp.repliconf import Repliconf
from amplifyp.settings import GLOBAL_REPLICATION_SETTINGS, ReplicationSettings


class PCR:
    """A class representing a PCR reaction.

    Attributes:
        template (DNA): The template DNA.
        settings (ReplicationSettings): The replication settings used.
        amplicon_generator (AmpliconGenerator): The internal amplicon generator.
        primers (list[Primer]): The primers configured for the reaction.
        amplicons (list[Amplicon]): The predicted amplicons.
    """

    def __init__(
        self,
        template: DNA,
        settings: ReplicationSettings = GLOBAL_REPLICATION_SETTINGS,
    ) -> None:
        """Initialise a PCR reaction.

        Args:
            template (DNA): The template DNA for the PCR reaction.
            settings (ReplicationSettings, optional): The replication settings.
                Defaults to GLOBAL_REPLICATION_SETTINGS.
        """
        self.template = template
        self.settings = settings
        self.amplicon_generator = AmpliconGenerator(self.template)

        self._primers: list[Primer] = []
        self._primer_repliconfs: dict[Primer, Repliconf] = {}
        self._amplicons: list[Amplicon] = []

    def add_primer(self, primer: Primer) -> None:
        """Add a primer to the PCR reaction.

        Args:
            primer (Primer): The primer to add.

        Raises:
            DuplicatedNameError: If a primer with the same name is
                already added.
            DuplicatedSequenceError: If a primer with the same sequence
                is already added.
        """
        seq_upper = primer.seq.upper()
        for p in self._primers:
            if p.name == primer.name:
                raise DuplicatedNameError(primer.name)
            if p.seq.upper() == seq_upper:
                raise DuplicatedSequenceError(primer.seq)

        repliconf = Repliconf(self.template, primer, self.settings)
        self.amplicon_generator.add_repliconf(repliconf)

        self._primers.append(primer)
        self._primer_repliconfs[primer] = repliconf
        self._amplicons.clear()

    def remove_primer(self, primer: Primer) -> None:
        """Remove a primer from the PCR reaction.

        Args:
            primer (Primer): The primer to remove.

        Raises:
            PrimerNotFoundError: If the primer is not found.
        """
        if primer not in self._primer_repliconfs:
            raise PrimerNotFoundError(primer)

        repliconf = self._primer_repliconfs.pop(primer)
        self._primers.remove(primer)
        self.amplicon_generator.remove_repliconf(repliconf)
        self._amplicons.clear()

    def add_primers(self, primers: Iterable[Primer]) -> None:
        """Add multiple primers to the PCR reaction.

        Args:
            primers (Iterable[Primer]): The primers to add.
        """
        for primer in primers:
            self.add_primer(primer)

    @property
    def primers(self) -> list[Primer]:
        """Get the primers used in the PCR reaction."""
        return self._primers.copy()

    def predict_amplicons(self) -> int:
        """Predict the amplicons for the PCR reaction.

        This updates the internal `amplicons` list by generating
        possible amplicons based on the current template and primers.

        Returns:
            int: The number of amplicons predicted.
        """
        self._amplicons = self.amplicon_generator.get_amplicons()
        return len(self._amplicons)

    @property
    def amplicons(self) -> list[Amplicon]:
        """Get the predicted amplicons.

        Returns:
            list[Amplicon]: A list of predicted amplicon objects.
        """
        return self._amplicons.copy()

    def __len__(self) -> int:
        """Return the number of primers currently in the PCR reaction."""
        return len(self._primers)

    def __contains__(self, primer: Primer) -> bool:
        """Check if a primer is part of the PCR reaction."""
        return primer in self._primer_repliconfs

    def __repr__(self) -> str:
        """Return a string representation of the PCR reaction."""
        return (
            f"PCR(template={self.template!r}, primers={len(self._primers)}, "
            f"amplicons={len(self._amplicons)})"
        )
