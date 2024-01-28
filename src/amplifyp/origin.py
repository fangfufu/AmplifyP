# -*- coding: utf-8 -*-
"""Amplify P - Replication related."""

from dataclasses import dataclass

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

    def stability(self) -> float:
        """Returns the stability of the origin.

        Returns:
            float: The stability of the origin.
        """
        r = self.settings.run_weight
        S = self.settings.base_pair_scores  # pylint: disable=invalid-name
        numerator: float = 0
        denominator: float = 0
        Rn: float = 0  # pylint: disable=invalid-name
        for k, (i, j) in enumerate(zip(self.primer, self.target)):
            numerator += r[k] * S[i, j]
            denominator += S.row_max(i)
            if r[k] > Rn:
                Rn = r[k]  # pylint: disable=invalid-name
            if S[i, j] <= 0:
                break
        score = numerator / (Rn * denominator)
        return score

    def quality(self) -> float:
        """Returns the quality of the origin.

        Returns:
            float: The quality of the origin.
        """
        cutoffs = self.settings.primability_cutoff + self.settings.stability_cutoff
        return (self.primability() + self.stability() - cutoffs) / (2 - cutoffs)
