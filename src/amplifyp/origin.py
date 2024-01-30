# -*- coding: utf-8 -*-
"""Amplify P - replication origin related."""

from dataclasses import dataclass
from .dna import DNA, Primer
from .settings import (
    Settings,
    LengthWiseWeightTbl,
    BasePairWeightsTbl,
)


@dataclass(frozen=True, slots=True)
class ReplicationOrigin:
    """A class representing the origin of replication.

    Attributes:
        target (str): The target DNA sequence as a string, in 3'-5' orientation.
        primer (str): The primer sequence as a string, in 3'-5' orientation.
        replication_config (ReplisomeConfig): The configuration for the replisome.
            Defaults to DEFAULT_REPLISOME_CONFIG.

    Raises:
        ValueError: If the length of the target is not equal to the length of the primer.
    """

    target: str
    primer: str
    settings: Settings

    def __post_init__(self) -> None:
        """Validates that the length of the target and primer are equal."""
        if len(self.target) != len(self.primer):
            raise ValueError("The target has to have the same length as the primer.")

    @property
    def primability(self) -> float:
        """Returns the primability of the origin.

        Returns:
            float: The primability of the origin.
        """
        m: LengthWiseWeightTbl = self.settings.match_weight
        S: BasePairWeightsTbl = (  # pylint: disable=invalid-name
            self.settings.base_pair_scores
        )
        numerator: float = 0
        denominator: float = 0
        for k, (i, j) in enumerate(zip(self.primer, self.target)):
            numerator += m[k] * S[i, j]
            denominator += m[k] * S.row_max(i)
        score = numerator / denominator
        return score

    @property
    def stability(self) -> float:
        """Returns the stability of the origin.

        Please note that the formula in Amplify 4's README is incorrect. This is
        direct translation from the source code.

        Returns:
            float: The stability of the origin.
        """
        r = self.settings.run_weights  # pylint: disable=invalid-name
        S = self.settings.base_pair_scores  # pylint: disable=invalid-name
        numerator: float = 0
        denominator: float = 0
        this_run_len: float = 0
        this_run_score: float = 0
        for i, j in zip(self.primer, self.target):
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
        """Returns the quality of the origin.

        Returns:
            float: The quality of the origin.
        """
        cutoffs = self.settings.primability_cutoff + self.settings.stability_cutoff
        return (self.primability + self.stability - cutoffs) / (2 - cutoffs)


class Amplify4RevOrigin(ReplicationOrigin):
    """
    Amplify4-style reverse replication origin.

    This class simplifies the process of constructing an Amplify4 style reverse
    replication origin. This is mainly used in testing.
    """

    def __init__(self, target: str, primer: str) -> None:
        """Constructor."""
        super().__init__(
            target=DNA(target).complement().sequence, primer=primer, settings=Settings()
        )


class Amplify4FwdOrigin(ReplicationOrigin):
    """
    Amplify4-style forward replication origin.

    This class simplifies the process of constructing an Amplify4 style forward
    replication origin. This is mainly used in testing.
    """

    def __init__(self, target: str, primer: str) -> None:
        """Constructor."""
        super().__init__(
            target=DNA(target).reverse().sequence,
            primer=Primer(primer).reverse().sequence,
            settings=Settings(),
        )
