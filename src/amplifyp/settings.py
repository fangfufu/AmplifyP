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

"""Settings-related classes and constants for AmplifyP.

This module contains the configuration classes and default values used for
calculating primability, stability, and other properties of replication
origins. The default values are based on the Amplify4 software.

For more information on the original Amplify4 software, see:
https://github.com/wrengels/Amplify4
"""

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Final

from .dna import Nucleotides


class LengthWiseWeightTbl:
    """A class representing a length-wise weight table.

    This class is used to score sequences based on their length, such as
    runs of identical nucleotides. It allows defining a default score (weight)
    that applies generally, along with specific overrides for certain lengths.

    Attributes:
        default_weight (float): The default weight returned for lengths not
            found in the overrides.
        overrides (dict[int, float]): A dictionary mapping specific lengths
            (int) to their corresponding weights (float).
    """

    def __init__(
        self,
        default_weight: float = 0,
        overrides: dict[int, float] | None = None,
    ) -> None:
        """Initialize a LengthWiseWeightTbl object.

        Args:
            default_weight (float, optional): The default weight to use when a
                length is not specified in overrides. Defaults to 0.
            overrides (dict[int, float], optional): A dictionary mapping
                specific lengths to their weights. Defaults to None.
        """
        if overrides is None:
            overrides = {}
        self.__default_weight = default_weight
        self.__overrides = overrides

    def __copy__(self) -> "LengthWiseWeightTbl":
        """Return a deep copy of this object.

        Returns:
            LengthWiseWeightTbl: A new object with the same weights.
        """
        return deepcopy(self)

    def copy(self) -> "LengthWiseWeightTbl":
        """Return a deep copy of this object.

        Returns:
            LengthWiseWeightTbl: A new object with the same weights.
        """
        return self.__copy__()

    def __getitem__(self, key: int) -> float:
        """Return the weight for a given length.

        Args:
            key (int): The length to look up.

        Returns:
            float: The weight associated with the length. Returns the specific
            override if present, otherwise returns the default weight.
        """
        if key in self.__overrides:
            return self.__overrides[key]
        return self.__default_weight

    def __setitem__(self, key: int, value: float) -> None:
        """Set the weight for a specific length.

        Args:
            key (int): The length to set the weight for.
            value (float): The weight value.
        """
        self.__overrides[key] = value


class BasePairWeightsTbl:
    """A class representing a nucleotide pairwise weight table.

    This table stores scores (weights) for every possible pair of nucleotides.
    It is typically used to score the interaction (binding strength) between
    a primer and a target DNA sequence.

    Attributes:
        row (str): The string of nucleotides corresponding to the table rows.
        col (str): The string of nucleotides corresponding to the table columns.
    """

    def __init__(self, row: str, col: str, weight: list[list[float]]) -> None:
        """Construct a BasePairWeightsTbl object.

        Args:
            row (str): A string representing the row labels (nucleotides) of
                the table. Gap characters are allowed but not stored in the
                table.
            col (str): A string representing the column labels (nucleotides) of
                the table.
            weight (list[list[float]]): A 2D list (matrix) of weights.
                `weight[i][j]` corresponds to the pair (`row[i]`, `col[j]`).

        Raises:
            ValueError: If the dimensions of the provided weight matrix do not
                match the lengths of the row and column strings (excluding
                gaps).
        """
        self.__row = row.upper()
        self.__col = col.upper()
        self.__weight: dict[tuple[str, str], float] = {}
        self.__row_max: dict[str, float] = {}

        # Optimized lookups
        self.__matrix: list[list[float]] = [
            [0.0] * len(self.__col) for _ in range(len(self.__row))
        ]
        self.__row_map: list[int] = [-1] * 128
        self.__col_map: list[int] = [-1] * 128

        # Populate row map
        for i, char in enumerate(self.__row):
            code = ord(char)
            if 0 <= code < 128:
                self.__row_map[code] = i
                # Handle lowercase
                lower_code = ord(char.lower())
                if 0 <= lower_code < 128:
                    self.__row_map[lower_code] = i

        # Populate col map
        for i, char in enumerate(self.__col):
            code = ord(char)
            if 0 <= code < 128:
                self.__col_map[code] = i
                # Handle lowercase
                lower_code = ord(char.lower())
                if 0 <= lower_code < 128:
                    self.__col_map[lower_code] = i

        # Expected row and column lengths
        exp_row_len = len(row) if Nucleotides.GAP not in row else len(row) - 1
        exp_col_len = len(col) if Nucleotides.GAP not in col else len(col) - 1

        if len(weight) != exp_row_len:
            raise ValueError(
                "BasePairWeightsTbl: row length mismatch at initialisation."
            )

        for i, row_val in enumerate(self.__row):
            if row_val != Nucleotides.GAP:
                # We never put the gap symbol in the table, hence the -1.
                if len(weight[i]) != exp_col_len:
                    raise ValueError(
                        "BasePairWeightsTbl: column length mismatch at "
                        "initialisation."
                    )
                self.__row_max[row_val] = max(weight[i])
            for j, col_val in enumerate(self.__col):
                val: float | int
                if Nucleotides.GAP in [row_val, col_val]:
                    val = 0
                else:
                    val = weight[i][j]
                self.__weight[row_val, col_val] = val
                self.__matrix[i][j] = val

    def row(self) -> str:
        """Return the string of row nucleotides.

        Returns:
        Returns:
            str: The row nucleotides, excluding the last character if it was a
                gap placeholder.
        """
        return self.__row[:-1]

    def column(self) -> str:
        """Return the string of column nucleotides.

        Returns:
        Returns:
            str: The column nucleotides, excluding the last character if it was
                a gap placeholder.
        """
        return self.__col[:-1]

    def row_max(self, row: str) -> float:
        """Return the maximum weight available for a given row nucleotide.

        Args:
            row (str): The nucleotide character for the row.

        Returns:
            float: The maximum weight value in that row.
        """
        row = row.upper()
        return self.__row_max[row]

    def __getitem__(self, key: tuple[str, str]) -> float:
        """Retrieve the weight for a specific nucleotide pair.

        Args:
            key (tuple[str, str]): A tuple containing two nucleotide characters.
                The lookup is case-insensitive.

        Returns:
            float: The weight associated with the pair.
        """
        i, j = key
        try:
            r = self.__row_map[ord(i)]
            c = self.__col_map[ord(j)]
            if r != -1 and c != -1:
                return self.__matrix[r][c]
        except (TypeError, IndexError):
            pass

        # Fallback to dictionary lookup
        i = i.upper()
        j = j.upper()
        return self.__weight[i, j]

    def __setitem__(self, key: tuple[str, str], value: float) -> None:
        """Set the weight for a specific nucleotide pair.

        Args:
            key (tuple[str, str]): A tuple containing two nucleotide characters.
                The key is case-insensitive.
            value (float): The new weight value.
        """
        i, j = key
        # We use upper case keys to be consistent with getitem expectation
        i_upper = i.upper()
        j_upper = j.upper()
        self.__weight[i_upper, j_upper] = value

        try:
            r = self.__row_map[ord(i)]
            c = self.__col_map[ord(j)]
            if r != -1 and c != -1:
                self.__matrix[r][c] = value
        except (TypeError, IndexError):
            pass

    def __len__(self) -> int:
        """Return the total number of entries in the weight table.

        Returns:
            int: The product of the number of row and column nucleotides.
        """
        return len(self.row()) * len(self.column())

    def __str__(self) -> str:
        """Return a string representation of the weight table.

        Returns:
        Returns:
            str: A string representation of the internal dictionary mapping
                pairs to weights.
        """
        return str(self.__weight)

    def __copy__(self) -> "BasePairWeightsTbl":
        """Return a deep copy of this object.

        Returns:
            BasePairWeightsTbl: A new object with the same weights.
        """
        return deepcopy(self)

    def copy(self) -> "BasePairWeightsTbl":
        """Return a deep copy of this object.

        Returns:
            BasePairWeightsTbl: A new object with the same weights.
        """
        return self.__copy__()


DEFAULT_MATCH_WEIGHTS: Final[LengthWiseWeightTbl] = LengthWiseWeightTbl(
    default_weight=1,
    overrides={
        0: 30,
        1: 20,
        2: 10,
        3: 10,
        4: 9,
        5: 9,
        6: 8,
        7: 7,
        8: 6,
        9: 5,
        10: 5,
        11: 4,
        12: 3,
        13: 2,
        14: 1,
    },
)
DEFAULT_RUN_WEIGHTS: Final[LengthWiseWeightTbl] = LengthWiseWeightTbl(
    default_weight=186,
    overrides={0: 100, 1: 150, 2: 175, 3: 182, 4: 186},
)

DEFAULT_BASE_PAIR_WEIGHTS: Final[BasePairWeightsTbl] = BasePairWeightsTbl(
    row=Nucleotides.PRIMER,
    col=Nucleotides.TEMPLATE,
    weight=[
        [100, 0, 0, 0, 30],
        [0, 100, 0, 0, 30],
        [0, 0, 100, 0, 30],
        [0, 0, 0, 100, 30],
        [0, 70, 0, 70, 30],
        [70, 70, 0, 0, 30],
        [0, 70, 70, 0, 30],
        [70, 0, 0, 70, 30],
        [0, 0, 70, 70, 30],
        [70, 0, 70, 0, 30],
        [50, 50, 0, 50, 30],
        [0, 50, 50, 50, 30],
        [50, 50, 50, 0, 30],
        [50, 0, 50, 50, 30],
        [30, 30, 30, 30, 30],
    ],
)


DEFAULT_PRIMABILITY_CUTOFF: Final[float] = 0.8


DEFAULT_STABILITY_CUTOFF: Final[float] = 0.4


@dataclass(slots=True)
class ReplicationSettings:
    """A configuration class for replication settings.

    This class aggregates all parameters required for analyzing replication
    origins. It includes scoring tables for base pairing and run lengths, as
    well as cutoff thresholds for filtering results.

    Attributes:
        base_pair_scores (BasePairWeightsTbl): Table of weights for nucleotide
            pairs. Defaults to `DEFAULT_BASE_PAIR_WEIGHTS`.
        match_weight (LengthWiseWeightTbl): Weights based on the length of
            matching segments. Defaults to `DEFAULT_MATCH_WEIGHTS`.
        run_weights (LengthWiseWeightTbl): Weights based on the length of
            consecutive runs. Defaults to `DEFAULT_RUN_WEIGHTS`.
        primability_cutoff (float): The minimum score required for a site to be
            considered a valid priming site (primability). Defaults to
            `DEFAULT_PRIMABILITY_CUTOFF`.
        stability_cutoff (float): The minimum stability score required.
            Defaults to `DEFAULT_STABILITY_CUTOFF`.
        amplify4_compatibility_mode (bool): For scoring, whether to use Amplify4
            compatibility mode. Defaults to `False`.
    """

    base_pair_scores: BasePairWeightsTbl = field(
        default_factory=lambda: DEFAULT_BASE_PAIR_WEIGHTS.copy()
    )
    match_weight: LengthWiseWeightTbl = field(
        default_factory=lambda: DEFAULT_MATCH_WEIGHTS.copy()
    )
    run_weights: LengthWiseWeightTbl = field(
        default_factory=lambda: DEFAULT_RUN_WEIGHTS.copy()
    )
    primability_cutoff: float = DEFAULT_PRIMABILITY_CUTOFF
    stability_cutoff: float = DEFAULT_STABILITY_CUTOFF
    amplify4_compatibility_mode: bool = False


GLOBAL_REPLICATION_SETTINGS: ReplicationSettings = ReplicationSettings()


@dataclass(slots=True)
class TMSettings:
    """Configuration for melting temperature calculations.

    Attributes:
        dna_conc (float): Total strand concentration in nM. Defaults to 50.
        dnap_conc (float): DNA Polymerase concentration. Not currently used in
            standard Tm. Defaults to 0.
        monovalent_salt_conc (float): Concentration of monovalent cations (Na+,
            K+, Tris+) in mM. Defaults to 50.
        divalent_salt_conc (float): Concentration of divalent cations (Mg++) in
            mM. Defaults to 1.5.
        dnTP_conc (float): Concentration of dNTPs in mM. Defaults to 0.
    """

    dna_conc: float = 50.0
    dnap_conc: float = 0.0
    monovalent_salt_conc: float = 50.0
    divalent_salt_conc: float = 1.5
    dnTP_conc: float = 0.0


GLOBAL_TM_SETTINGS: TMSettings = TMSettings()

DEFAULT_AMPLIFY4_TM_ENTHALPY: Final[BasePairWeightsTbl] = BasePairWeightsTbl(
    row=Nucleotides.SINGLE + Nucleotides.WILDCARD,
    col=Nucleotides.SINGLE + Nucleotides.WILDCARD,
    weight=[
        [110, 78, 58, 119, 94],
        [56, 91, 60, 58, 81],
        [65, 86, 91, 78, 78],
        [111, 65, 56, 110, 105],
        [65, 86, 60, 78, 78],
    ],
)

DEFAULT_AMPLIFY4_TM_ENTROPY: Final[BasePairWeightsTbl] = BasePairWeightsTbl(
    row=Nucleotides.SINGLE + Nucleotides.WILDCARD,
    col=Nucleotides.SINGLE + Nucleotides.WILDCARD,
    weight=[
        [266, 208, 129, 278, 220],
        [135, 240, 169, 129, 168],
        [173, 239, 240, 208, 215],
        [267, 173, 135, 266, 210],
        [173, 239, 169, 208, 215],
    ],
)


@dataclass(slots=True)
class Amplify4TMSettings:
    """Configuration specific to the Amplify4 melting temperature algorithm.

    Attributes:
        dna_conc (float): Total strand concentration in nM. Defaults to 50.
        monovalent_salt_conc (float): Concentration of monovalent cations (Na+,
            K+, Tris+) in mM. Defaults to 50.
        enthalpy (list[list[int]]): Enthalpy values (5x5 matrix).
            Defaults to DEFAULT_AMPLIFY4_ENTHALPY.
        entropy (list[list[int]]): Entropy values (5x5 matrix).
            Defaults to DEFAULT_AMPLIFY4_ENTROPY.
    """

    dna_conc: float = 50.0
    monovalent_salt_conc: float = 50.0
    enthalpy: BasePairWeightsTbl = field(
        default_factory=lambda: DEFAULT_AMPLIFY4_TM_ENTHALPY.copy()
    )
    entropy: BasePairWeightsTbl = field(
        default_factory=lambda: DEFAULT_AMPLIFY4_TM_ENTROPY.copy()
    )


GLOBAL_AMPLIFY4_TM_SETTINGS: Amplify4TMSettings = Amplify4TMSettings()

DEFAULT_PRIMER_DIMER_WEIGHTS: Final[BasePairWeightsTbl] = BasePairWeightsTbl(
    row=Nucleotides.PRIMER,
    col=Nucleotides.PRIMER,
    weight=[
        [-20, -20, -20, 30, 5, -20, -20, 5, 5, -20, -3, -3, -20, -3, -8],
        [-20, -20, 20, -20, -20, -20, 0, -20, 0, 0, -20, -7, -7, -7, -10],
        [-20, 20, -20, -20, 0, 0, 0, -20, -20, -20, -7, -7, -7, -20, -10],
        [30, -20, -20, -20, -20, 5, -20, 5, -20, 5, -3, -20, -3, -3, -8],
        [5, -20, 0, -20, -20, -8, -10, -8, -10, 3, -12, -13, -5, -5, -9],
        [-20, -20, 0, 5, -8, -20, -10, -8, 3, -10, -12, -5, -13, -5, -9],
        [-20, 0, 0, -20, -10, -10, 0, -20, -10, -10, -13, -7, -7, -13, -10],
        [5, -20, -20, 5, -8, -8, -20, 5, -8, -8, -3, -12, -12, -3, -8],
        [5, 0, -20, -20, -10, 3, -10, -8, -20, -8, -5, -13, -5, -12, -9],
        [-20, 0, -20, 5, 3, -10, -10, -8, -8, -20, -5, -5, -13, -12, -9],
        [-3, -20, -7, -3, -12, -12, -13, -3, -5, -5, -9, -10, -10, -4, -8],
        [-3, -7, -7, -20, -13, -5, -7, -12, -13, -5, -10, -11, -6, -10, -9],
        [-20, -7, -7, -3, -5, -13, -7, -12, -5, -13, -10, -6, -11, -10, -9],
        [-3, -7, -20, -3, -5, -5, -13, -3, -12, -12, -4, -10, -10, -9, -8],
        [-8, -10, -10, -8, -9, -9, -10, -8, -9, -9, -8, -9, -9, -8, -9],
    ],
)

DEFAULT_PRIMER_DIMER_OVERLAP: Final[int] = 3
DEFAULT_PRIMER_DIMER_THRESHOLD: Final[float] = 60.0


@dataclass(slots=True, frozen=True)
class PrimerDimerSettings:
    """Settings for primer dimer analysis.

    Attributes:
        weights (BasePairWeightsTbl): Scoring table for base pair interactions.
        min_overlap (int): Minimum overlap length for a primer dimer.
        threshold (float): Minimum quality score for a primer dimer.
    """

    weights: BasePairWeightsTbl = field(
        default_factory=lambda: DEFAULT_PRIMER_DIMER_WEIGHTS.copy()
    )
    min_overlap: int = DEFAULT_PRIMER_DIMER_OVERLAP
    threshold: float = DEFAULT_PRIMER_DIMER_THRESHOLD


GLOBAL_PRIMER_DIMER_SETTINGS: PrimerDimerSettings = PrimerDimerSettings()
