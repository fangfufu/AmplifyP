"""Replication origin-related classes for AmplifyP."""

from dataclasses import dataclass

from .dna import DNA, Primer
from .settings import BasePairWeightsTbl, LengthWiseWeightTbl, Settings


@dataclass(frozen=True, slots=True)
class ReplicationOrigin:
    """A class representing a replication origin.

    A replication origin is a specific sequence in a genome at which replication
    is initiated. This class calculates various properties of a replication
    origin, such as its primability, stability, and quality, based on the
    target and primer sequences and the provided settings.

    Attributes:
        target (str): The target DNA sequence as a string, in 3'-5'
          orientation.
        primer (str): The primer sequence as a string, in 3'-5' orientation.
        settings (Settings): The configuration for the replisome.

    Raises:
        ValueError: If the length of the target is not equal to the length of
                    the primer.
    """

    target: str
    primer: str
    settings: Settings

    def __post_init__(self) -> None:
        """Validate that the length of the target and primer are equal."""
        if len(self.target) != len(self.primer):
            raise ValueError("The target has to have the same length as the primer.")

    @property
    def primability(self) -> float:
        """Calculate the primability of the replication origin.

        Primability is a measure of how well a primer can bind to a target DNA
        sequence. It is calculated based on the base pairing scores and match
        weights defined in the settings.

        Returns:
            float: The primability of the origin, a value between 0 and 1.
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
        """Calculate the stability of the replication origin.

        Stability is a measure of the strength of the binding between the primer
        and the target DNA sequence. It is calculated based on the base pairing
        scores and run weights defined in the settings.

        Please note that the formula in Amplify 4's README is incorrect. This
        is a direct translation from the source code.

        Returns:
            float: The stability of the origin, a value between 0 and 1.
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
        """Calculate the quality of the replication origin.

        The quality is a combined measure of primability and stability. It is
        calculated as the average of the primability and stability scores,
        adjusted by their respective cutoffs.

        Returns:
            float: The quality of the origin, a value between 0 and 1.
        """
        cutoffs = self.settings.primability_cutoff + self.settings.stability_cutoff
        return (self.primability + self.stability - cutoffs) / (2 - cutoffs)


class Amplify4RevOrigin(ReplicationOrigin):
    """A class representing an Amplify4-style reverse replication origin.

    This class simplifies the process of constructing an Amplify4-style reverse
    replication origin. It is mainly used in testing and for compatibility with
    the Amplify4 software.
    """

    def __init__(self, target: str, primer: str) -> None:
        """Construct an Amplify4RevOrigin object.

        Args:
            target (str): The target DNA sequence.
            primer (str): The primer sequence.
        """
        super().__init__(
            target=DNA(target).complement().seq, primer=primer, settings=Settings()
        )


class Amplify4FwdOrigin(ReplicationOrigin):
    """A class representing an Amplify4-style forward replication origin.

    This class simplifies the process of constructing an Amplify4-style forward
    replication origin. It is mainly used in testing and for compatibility with
    the Amplify4 software.
    """

    def __init__(self, target: str, primer: str) -> None:
        """Construct an Amplify4FwdOrigin object.

        Args:
            target (str): The target DNA sequence.
            primer (str): The primer sequence.
        """
        super().__init__(
            target=DNA(target).reverse().seq,
            primer=Primer(primer).reverse().seq,
            settings=Settings(),
        )
