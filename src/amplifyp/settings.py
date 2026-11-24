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

"""Settings-related classes and constants for AmplifyP.

This module contains the configuration classes and default values used for
calculating primability, stability, and other properties of replication
origins. The default values are based on the Amplify4 software.

For more information on the original Amplify4 software, see:
https://github.com/wrengels/Amplify4
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, Final, Self, cast

from .dna import Nucleotides
from .errors import ColumnLengthMismatchError, RowLengthMismatchError


def _coerce_float(value: Any, default: float) -> float:
    """Convert value to a finite float, falling back to default on failure."""
    if value is None or value == "" or isinstance(value, bool):
        return default
    try:
        parsed = float(value)
    except (ValueError, TypeError):
        return default
    return parsed if math.isfinite(parsed) else default


def _coerce_int(value: Any, default: int) -> int:
    """Convert value to int, falling back to default on failure."""
    if value is None or value == "" or isinstance(value, bool):
        return default
    try:
        return int(value)
    except (ValueError, TypeError):
        return default


def coerce_bool(value: Any, default: bool) -> bool:
    """Convert bool, int, or string flag values into a boolean."""
    if value is None or value == "":
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        if isinstance(value, float) and not math.isfinite(value):
            return default
        return value != 0
    if isinstance(value, str):
        norm = value.strip().lower()
        if norm in ("true", "1", "yes", "on"):
            return True
        if norm in ("false", "0", "no", "off"):
            return False
    return default


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
        """Initialise a LengthWiseWeightTbl object.

        Args:
            default_weight (float, optional): The default weight to use when a
                length is not specified in overrides. Defaults to 0.
            overrides (dict[int, float], optional): A dictionary mapping
                specific lengths to their weights. Defaults to None.
        """
        if overrides is None:
            overrides = {}
        self._default_weight = default_weight
        self._overrides = overrides

    @property
    def default_weight(self) -> float:
        """Return the default weight."""
        return self._default_weight

    @property
    def overrides(self) -> dict[int, float]:
        """Return the dictionary of length overrides."""
        return self._overrides

    def __copy__(self) -> LengthWiseWeightTbl:
        """Return a deep copy of this object.

        Returns:
            LengthWiseWeightTbl: A new object with the same weights.
        """
        return LengthWiseWeightTbl(
            default_weight=self._default_weight,
            overrides=dict(self._overrides),
        )

    def copy(self) -> LengthWiseWeightTbl:
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
        return self._overrides.get(key, self._default_weight)

    def __setitem__(self, key: int, value: float) -> None:
        """Set the weight for a specific length.

        Args:
            key (int): The length to set the weight for.
            value (float): The weight value.
        """
        self._overrides[key] = value

    def __eq__(self, other: object) -> bool:
        """Check equality between two LengthWiseWeightTbl instances."""
        if not isinstance(other, LengthWiseWeightTbl):
            return NotImplemented
        return (
            self._default_weight == other._default_weight
            and self._overrides == other._overrides
        )

    def __repr__(self) -> str:
        """Return a string representation of the LengthWiseWeightTbl object."""
        return (
            f"LengthWiseWeightTbl(default_weight={self._default_weight!r}, "
            f"overrides={self._overrides!r})"
        )

    def to_dict(self, prefix: str) -> dict[str, Any]:
        """Serialise weights into a flat dictionary keyed by ``prefix``."""
        result: dict[str, Any] = {f"{prefix}_default": self._default_weight}
        for idx in sorted(self._overrides):
            result[f"{prefix}_{idx}"] = self._overrides[idx]
        return result

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        prefix: str,
        default_tbl: LengthWiseWeightTbl,
    ) -> LengthWiseWeightTbl:
        """Construct a LengthWiseWeightTbl from a flat dictionary."""
        default_key = f"{prefix}_default"
        prefix_underscore = f"{prefix}_"
        has_flat_keys = default_key in data or any(
            isinstance(cast("object", k), str)
            and k.startswith(prefix_underscore)
            and k[len(prefix_underscore) :].isascii()
            and k[len(prefix_underscore) :].isdigit()
            for k in data
        )
        if not has_flat_keys:
            return default_tbl.copy()

        default_weight = _coerce_float(
            data.get(default_key), default_tbl.default_weight
        )
        overrides: dict[int, float] = (
            {} if default_key in data else dict(default_tbl.overrides)
        )
        for key, val in data.items():
            if (
                isinstance(cast("object", key), str)
                and key.startswith(prefix_underscore)
                and key != default_key
            ):
                suffix = key[len(prefix_underscore) :]
                if suffix.isascii() and suffix.isdigit():
                    idx = int(suffix)
                    fallback_weight = (
                        default_tbl.default_weight
                        if default_key in data
                        else default_tbl[idx]
                    )
                    overrides[idx] = _coerce_float(val, fallback_weight)
        return cls(default_weight=default_weight, overrides=overrides)


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
        self._row = row.upper()
        self._col = col.upper()
        self._weight: dict[tuple[str, str], float] = {}
        self._row_max: dict[str, float] = {}

        # Optimised lookups (supporting 256 ASCII / Latin-1 characters)
        self._matrix: list[list[float]] = [
            [0.0] * len(self._col) for _ in range(len(self._row))
        ]
        self._row_map: list[int] = [-1] * 256
        self._col_map: list[int] = [-1] * 256

        # Populate row map
        for i, char in enumerate(self._row):
            code = ord(char)
            if 0 <= code < 256:
                self._row_map[code] = i
                lower_code = ord(char.lower())
                if 0 <= lower_code < 256:
                    self._row_map[lower_code] = i

        # Populate col map
        for i, char in enumerate(self._col):
            code = ord(char)
            if 0 <= code < 256:
                self._col_map[code] = i
                lower_code = ord(char.lower())
                if 0 <= lower_code < 256:
                    self._col_map[lower_code] = i

        # Expected row and column lengths
        exp_row_len = len(row) if Nucleotides.GAP not in row else len(row) - 1
        exp_col_len = len(col) if Nucleotides.GAP not in col else len(col) - 1

        if len(weight) != exp_row_len:
            raise RowLengthMismatchError()

        gap = Nucleotides.GAP

        for i, row_val in enumerate(self._row):
            if row_val != gap:
                # We never put the gap symbol in the table, hence the -1.
                if len(weight[i]) != exp_col_len:
                    raise ColumnLengthMismatchError()
                self._row_max[row_val] = max(weight[i])
            for j, col_val in enumerate(self._col):
                val: float | int
                if row_val == gap or col_val == gap:
                    val = 0
                else:
                    val = weight[i][j]
                self._weight[row_val, col_val] = val
                self._matrix[i][j] = val

    def row(self) -> str:
        """Return the string of row nucleotides.

        Returns:
            str: The row nucleotides, excluding the last character if it was a
                gap placeholder.
        """
        if self._row.endswith(Nucleotides.GAP):
            return self._row[:-1]
        return self._row

    def column(self) -> str:
        """Return the string of column nucleotides.

        Returns:
            str: The column nucleotides, excluding the last character if it was
                a gap placeholder.
        """
        if self._col.endswith(Nucleotides.GAP):
            return self._col[:-1]
        return self._col

    def row_max(self, row: str) -> float:
        """Return the maximum weight available for a given row nucleotide.

        Args:
            row (str): The nucleotide character for the row.

        Returns:
            float: The maximum weight value in that row.
        """
        row = row.upper()
        return self._row_max[row]

    def __getitem__(self, key: tuple[str, str]) -> float:
        """Retrieve the weight for a specific nucleotide pair.

        Args:
            key (tuple[str, str]): A tuple containing two nucleotide characters.
                The lookup is case-insensitive.

        Returns:
            float: The weight associated with the pair.
        """
        i, j = key
        code_i = ord(i)
        code_j = ord(j)
        if 0 <= code_i < 256 and 0 <= code_j < 256:
            r = self._row_map[code_i]
            c = self._col_map[code_j]
            if r != -1 and c != -1:
                return self._matrix[r][c]

        i_u, j_u = i.upper(), j.upper()
        return self._weight[(i_u, j_u)]

    def __setitem__(self, key: tuple[str, str], value: float) -> None:
        """Set the weight for a specific nucleotide pair.

        Args:
            key (tuple[str, str]): A tuple containing two nucleotide characters.
                The key is case-insensitive.
            value (float): The new weight value.
        """
        i, j = key
        i_upper = i.upper()
        j_upper = j.upper()
        self._weight[i_upper, j_upper] = value

        code_i = ord(i)
        code_j = ord(j)
        if 0 <= code_i < 256 and 0 <= code_j < 256:
            r = self._row_map[code_i]
            c = self._col_map[code_j]
            if r != -1 and c != -1:
                self._matrix[r][c] = value
                self._row_max[i_upper] = max(self._matrix[r])

    def __len__(self) -> int:
        """Return the total number of entries in the weight table.

        Returns:
            int: The product of the number of row and column nucleotides.
        """
        return len(self.row()) * len(self.column())

    def __eq__(self, other: object) -> bool:
        """Check equality between two BasePairWeightsTbl instances."""
        if not isinstance(other, BasePairWeightsTbl):
            return NotImplemented
        return (
            self._row == other._row
            and self._col == other._col
            and self._weight == other._weight
        )

    def __str__(self) -> str:
        """Return a string representation of the weight table.

        Returns:
            str: A string representation of the internal dictionary mapping
                pairs to weights.
        """
        return str(self._weight)

    def __repr__(self) -> str:
        """Return a string representation of the BasePairWeightsTbl object."""
        return f"BasePairWeightsTbl(row={self.row()!r}, col={self.column()!r})"

    def __copy__(self) -> BasePairWeightsTbl:
        """Return a deep copy of this object.

        Returns:
            BasePairWeightsTbl: A new object with the same weights.
        """
        new_tbl = BasePairWeightsTbl.__new__(BasePairWeightsTbl)
        new_tbl._row = self._row
        new_tbl._col = self._col
        new_tbl._weight = dict(self._weight)
        new_tbl._row_max = dict(self._row_max)
        new_tbl._matrix = [list(r) for r in self._matrix]
        new_tbl._row_map = list(self._row_map)
        new_tbl._col_map = list(self._col_map)
        return new_tbl

    def copy(self) -> BasePairWeightsTbl:
        """Return a deep copy of this object.

        Returns:
            BasePairWeightsTbl: A new object with the same weights.
        """
        return self.__copy__()

    def to_dict(self, prefix: str) -> dict[str, Any]:
        """Serialise weights into a flat dictionary keyed by ``prefix``."""
        result: dict[str, Any] = {}
        for r_char in self.row():
            for c_char in self.column():
                result[f"{prefix}_{r_char}_{c_char}"] = self[r_char, c_char]
        return result

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        prefix: str,
        default_tbl: BasePairWeightsTbl,
    ) -> BasePairWeightsTbl:
        """Construct a BasePairWeightsTbl from a flat dictionary."""
        prefix_underscore = f"{prefix}_"
        has_prefix_keys = any(
            isinstance(cast("object", k), str)
            and k.startswith(prefix_underscore)
            for k in data
        )
        if not has_prefix_keys:
            return default_tbl.copy()

        weights: list[list[float]] = []
        for r_char in default_tbl.row():
            row_vals: list[float] = []
            for c_char in default_tbl.column():
                key = f"{prefix_underscore}{r_char}_{c_char}"
                try:
                    cell_default = default_tbl[r_char, c_char]
                except KeyError:
                    cell_default = 0.0
                row_vals.append(_coerce_float(data.get(key), cell_default))
            weights.append(row_vals)
        return cls(
            row=default_tbl._row,
            col=default_tbl._col,
            weight=weights,
        )


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

DEFAULT_PRIMER_DIMER_SYMBOL_THRESHOLD: Final[float] = 10.0


@dataclass(slots=True)
class ReplicationSettings:
    """A configuration class for replication settings.

    This class aggregates all parameters required for analysing replication
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
        primer_dimer_scores (BasePairWeightsTbl): Table of weights for primer
            dimer. Defaults to `DEFAULT_PRIMER_DIMER_WEIGHTS`.
    """

    base_pair_scores: BasePairWeightsTbl = field(
        default_factory=DEFAULT_BASE_PAIR_WEIGHTS.copy
    )
    primer_dimer_scores: BasePairWeightsTbl = field(
        default_factory=DEFAULT_PRIMER_DIMER_WEIGHTS.copy
    )
    match_weight: LengthWiseWeightTbl = field(
        default_factory=DEFAULT_MATCH_WEIGHTS.copy
    )
    run_weights: LengthWiseWeightTbl = field(
        default_factory=DEFAULT_RUN_WEIGHTS.copy
    )
    primability_cutoff: float = DEFAULT_PRIMABILITY_CUTOFF
    stability_cutoff: float = DEFAULT_STABILITY_CUTOFF
    amplify4_compatibility_mode: bool = False

    def to_dict(self) -> dict[str, Any]:
        """Serialise replication settings into a flat dictionary."""
        data: dict[str, Any] = {
            "primability_cutoff": self.primability_cutoff,
            "stability_cutoff": self.stability_cutoff,
            "amplify4_compatibility_mode": self.amplify4_compatibility_mode,
        }
        data.update(self.base_pair_scores.to_dict("bp_score"))
        if self.primer_dimer_scores != DEFAULT_PRIMER_DIMER_WEIGHTS:
            data.update(self.primer_dimer_scores.to_dict("rep_pd_score"))
        data.update(self.match_weight.to_dict("match_weight"))
        data.update(self.run_weights.to_dict("run_weight"))
        return data

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Construct ReplicationSettings from a dictionary."""
        compat_val = data.get(
            "amplify4_compatibility_mode", data.get("amp4_compat", False)
        )
        return cls(
            base_pair_scores=BasePairWeightsTbl.from_dict(
                data, "bp_score", DEFAULT_BASE_PAIR_WEIGHTS
            ),
            primer_dimer_scores=BasePairWeightsTbl.from_dict(
                data, "rep_pd_score", DEFAULT_PRIMER_DIMER_WEIGHTS
            ),
            match_weight=LengthWiseWeightTbl.from_dict(
                data, "match_weight", DEFAULT_MATCH_WEIGHTS
            ),
            run_weights=LengthWiseWeightTbl.from_dict(
                data, "run_weight", DEFAULT_RUN_WEIGHTS
            ),
            primability_cutoff=_coerce_float(
                data.get("primability_cutoff"), DEFAULT_PRIMABILITY_CUTOFF
            ),
            stability_cutoff=_coerce_float(
                data.get("stability_cutoff"), DEFAULT_STABILITY_CUTOFF
            ),
            amplify4_compatibility_mode=coerce_bool(compat_val, False),
        )


GLOBAL_REPLICATION_SETTINGS: ReplicationSettings = ReplicationSettings()


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
        dntp_conc (float): Concentration of dNTPs in mM. Defaults to 0.
        enthalpy (BasePairWeightsTbl): Enthalpy values (5x5 matrix).
            Defaults to DEFAULT_AMPLIFY4_TM_ENTHALPY.
        entropy (BasePairWeightsTbl): Entropy values (5x5 matrix).
            Defaults to DEFAULT_AMPLIFY4_TM_ENTROPY.
    """

    dna_conc: float = 50.0
    dnap_conc: float = 0.0
    monovalent_salt_conc: float = 50.0
    divalent_salt_conc: float = 1.5
    dntp_conc: float = 0.0
    enthalpy: BasePairWeightsTbl = field(
        default_factory=DEFAULT_AMPLIFY4_TM_ENTHALPY.copy
    )
    entropy: BasePairWeightsTbl = field(
        default_factory=DEFAULT_AMPLIFY4_TM_ENTROPY.copy
    )

    def to_dict(self) -> dict[str, Any]:
        """Serialise thermodynamic concentration and matrix settings."""
        data: dict[str, Any] = {
            "tm_dna_conc": self.dna_conc,
            "tm_dnap_conc": self.dnap_conc,
            "tm_mono_salt": self.monovalent_salt_conc,
            "tm_div_salt": self.divalent_salt_conc,
            "tm_dNTP_conc": self.dntp_conc,
        }
        if self.enthalpy != DEFAULT_AMPLIFY4_TM_ENTHALPY:
            data.update(self.enthalpy.to_dict("tm_enthalpy"))
        if self.entropy != DEFAULT_AMPLIFY4_TM_ENTROPY:
            data.update(self.entropy.to_dict("tm_entropy"))
        return data

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Construct TMSettings from a dictionary."""
        defaults = cls()
        return cls(
            dna_conc=_coerce_float(
                data.get("tm_dna_conc", data.get("dna_conc")),
                defaults.dna_conc,
            ),
            dnap_conc=_coerce_float(
                data.get("tm_dnap_conc", data.get("dnap_conc")),
                defaults.dnap_conc,
            ),
            monovalent_salt_conc=_coerce_float(
                data.get("tm_mono_salt", data.get("monovalent_salt_conc")),
                defaults.monovalent_salt_conc,
            ),
            divalent_salt_conc=_coerce_float(
                data.get("tm_div_salt", data.get("divalent_salt_conc")),
                defaults.divalent_salt_conc,
            ),
            dntp_conc=_coerce_float(
                data.get("tm_dNTP_conc", data.get("dntp_conc")),
                defaults.dntp_conc,
            ),
            enthalpy=BasePairWeightsTbl.from_dict(
                data, "tm_enthalpy", DEFAULT_AMPLIFY4_TM_ENTHALPY
            ),
            entropy=BasePairWeightsTbl.from_dict(
                data, "tm_entropy", DEFAULT_AMPLIFY4_TM_ENTROPY
            ),
        )


GLOBAL_TM_SETTINGS: TMSettings = TMSettings()


DEFAULT_PRIMER_DIMER_OVERLAP: Final[int] = 3
DEFAULT_PRIMER_DIMER_THRESHOLD: Final[float] = 60.0


@dataclass(slots=True, frozen=True)
class PrimerDimerSettings:
    """Settings for primer dimer analysis.

    Attributes:
        weights (BasePairWeightsTbl): Scoring table for base pair interactions.
        min_overlap (int): Minimum overlap length for a primer dimer.
        threshold (float): Minimum quality score for a primer dimer.
        symbol_threshold (float): Threshold score for ':' bond symbol
            in primer dimers. Defaults to
            `DEFAULT_PRIMER_DIMER_SYMBOL_THRESHOLD`.
    """

    weights: BasePairWeightsTbl = field(
        default_factory=DEFAULT_PRIMER_DIMER_WEIGHTS.copy
    )
    min_overlap: int = DEFAULT_PRIMER_DIMER_OVERLAP
    threshold: float = DEFAULT_PRIMER_DIMER_THRESHOLD
    symbol_threshold: float = DEFAULT_PRIMER_DIMER_SYMBOL_THRESHOLD

    def to_dict(self) -> dict[str, Any]:
        """Serialise primer dimer settings into a flat dictionary."""
        data: dict[str, Any] = {
            "pd_min_overlap": self.min_overlap,
            "pd_threshold": self.threshold,
            "symbol_threshold": self.symbol_threshold,
        }
        data.update(self.weights.to_dict("pd_score"))
        return data

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Construct PrimerDimerSettings from a dictionary."""
        return cls(
            weights=BasePairWeightsTbl.from_dict(
                data, "pd_score", DEFAULT_PRIMER_DIMER_WEIGHTS
            ),
            min_overlap=_coerce_int(
                data.get("pd_min_overlap", data.get("min_overlap")),
                DEFAULT_PRIMER_DIMER_OVERLAP,
            ),
            threshold=_coerce_float(
                data.get("pd_threshold", data.get("threshold")),
                DEFAULT_PRIMER_DIMER_THRESHOLD,
            ),
            symbol_threshold=_coerce_float(
                data.get("symbol_threshold"),
                DEFAULT_PRIMER_DIMER_SYMBOL_THRESHOLD,
            ),
        )


GLOBAL_PRIMER_DIMER_SETTINGS: PrimerDimerSettings = PrimerDimerSettings()
