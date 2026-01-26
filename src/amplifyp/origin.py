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

"""Replication origin-related classes for AmplifyP."""

from dataclasses import dataclass, field
from math import trunc

from .dna import DNA, Primer
from .settings import (
    GLOBAL_REPLICATION_SETTINGS,
    BasePairWeightsTbl,
    LengthWiseWeightTbl,
    ReplicationSettings,
)


@dataclass(frozen=True, slots=True)
class ReplicationOrigin:
    """A class representing a potential replication origin.

    A replication origin is a site on the template DNA where the primer binds
    and replication initiates. This class calculates the binding properties
    of a primer to a specific target sequence, including primability, stability,
    and an overall quality score.

    Attributes:
        target (str): The target DNA sequence as a string, aligned with the
            primer. It must be in the 3'-5' orientation relative to the
            primer's binding.
        primer (str): The primer sequence as a string, typically in 3'-5'
            orientation for calculation purposes.
        settings (Settings): The configuration object containing scoring tables
            and cutoff thresholds.
    """

    target: str
    primer: str
    settings: ReplicationSettings = field(
        default_factory=lambda: GLOBAL_REPLICATION_SETTINGS
    )

    def __post_init__(self) -> None:
        """Validate that the target and primer have equal lengths.

        Raises:
            ValueError: If the lengths of `target` and `primer` do not match.
        """
        if len(self.target) != len(self.primer):
            raise ValueError(
                "The target has to have the same length as the primer."
            )

    @property
    def primability(self) -> float:
        """Calculate the primability score of the replication origin.

        Primability estimates the likelihood of the primer binding to the
        target.
        It is a weighted average of base-pairing scores, where weights are
        determined by the position of the base pair (length-wise match
        weights).

        Returns:
            float: The primability score, ranging from 0.0 to 1.0.
        """
        m: LengthWiseWeightTbl = self.settings.match_weight
        S: BasePairWeightsTbl = self.settings.base_pair_scores
        numerator: float = 0
        denominator: float = 0
        for k, (i, j) in enumerate(zip(self.primer, self.target, strict=False)):
            numerator += m[k] * S[i, j]
            denominator += m[k] * S.row_max(i)
        score = numerator / denominator
        return score

    @property
    def stability(self) -> float:
        """Calculate the stability score of the replication origin.

        Stability measures the thermodynamic or structural strength of the
        primer-target duplex. It considers consecutive runs of matching bases,
        applying weights based on the length of these runs.

        Note:
            The formula used here is a direct translation from the Amplify 4
            source code, which differs from the formula described in the
            Amplify 4 README.

        Returns:
            float: The stability score, ranging from 0.0 to 1.0.
        """
        r = self.settings.run_weights
        S = self.settings.base_pair_scores
        numerator: float = 0
        denominator: float = 0
        this_run_len: float = 0
        this_run_score: float = 0
        for i, j in zip(self.primer, self.target, strict=False):
            denominator += S.row_max(i)
            if S[i, j] > 0:
                this_run_len += 1
                this_run_score += S[i, j]
            else:
                # N.B. that each run group is scored using the same run score!
                # We have to allow a running length of 0 here.
                numerator += r[int(max(0, this_run_len - 1))] * this_run_score
                this_run_len = 0
                this_run_score = 0
        # Allows for finishing during a run:
        numerator += r[int(max(0, this_run_len - 1))] * this_run_score
        # We multiply the denominator by the largest score that this primer
        # can obtain.
        score = numerator / (denominator * r[int(max(0, len(self.primer) - 1))])
        return score

    @property
    def quality(self) -> float:
        """Calculate the overall quality score of the replication origin.

        The quality score combines primability and stability into a single
        metric.
        It averages the two scores after adjusting for their respective
        cutoffs.

        Returns:
            float: The quality score. Can be negative if scores are below
                cutoffs.
        """
        cutoffs = (
            self.settings.primability_cutoff + self.settings.stability_cutoff
        )
        if not self.settings.amplify4_compatibility_mode:
            return (self.primability + self.stability - cutoffs) / (2 - cutoffs)
        else:
            primability = trunc(self.primability * 100) / 100
            stability = trunc(self.stability * 100) / 100
            return (primability + stability - cutoffs) / (2 - cutoffs)


class Amplify4RevOrigin(ReplicationOrigin):
    """A helper class for creating an Amplify4-style reverse replication origin.

    This class facilitates the creation of a `ReplicationOrigin` for the reverse
    strand by automatically complementing the target sequence and using default
    settings. It is primarily used for testing and ensuring compatibility with
    Amplify 4 logic.
    """

    def __init__(self, target: str, primer: str) -> None:
        """Initialize an Amplify4RevOrigin object.

        Args:
            target (str): The target DNA sequence. The complement will be used.
            primer (str): The primer sequence.
        """
        super().__init__(
            target=DNA(target).complement().seq,
            primer=primer,
            settings=ReplicationSettings(amplify4_compatibility_mode=True),
        )


class Amplify4FwdOrigin(ReplicationOrigin):
    """A helper class for creating an Amplify4-style forward replication origin.

    This class facilitates the creation of a `ReplicationOrigin` for the
    forward strand by automatically reversing the target and primer sequences
    (as required by the internal calculation logic) and using default settings.

    It is primarily used for testing and ensuring compatibility with Amplify 4
    logic.
    """

    def __init__(self, target: str, primer: str) -> None:
        """Initialize an Amplify4FwdOrigin object.

        Args:
            target (str): The target DNA sequence.
            primer (str): The primer sequence.
        """
        super().__init__(
            target=DNA(target).reverse().seq,
            primer=Primer(primer).reverse().seq,
            settings=ReplicationSettings(amplify4_compatibility_mode=True),
        )
