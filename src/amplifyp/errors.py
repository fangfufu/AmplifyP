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

"""Custom exceptions for AmplifyP."""


class AmplifyPError(Exception):
    """Base exception for all AmplifyP errors."""


class DuplicateRepliconfError(AmplifyPError, ValueError):
    """Exception raised when attempting to add a duplicate Repliconf.

    This error occurs when a `Repliconf` object is added to an
    `AmpliconGenerator` that already contains an identical configuration (same
    primer and template).
    """

    def __init__(
        self,
        repliconf: object | None = None,
        message: str = "The Repliconf is already in the AmpliconGenerator.",
    ) -> None:
        """Initialise the duplicate repliconf error.

        Args:
            repliconf: The duplicate Repliconf object, if available.
            message: The error message describing the duplicate repliconf.
        """
        super().__init__(message)
        self.repliconf = repliconf


class PrimerNotFoundError(AmplifyPError, ValueError):
    """Exception raised when a primer is not found."""

    def __init__(self, primer: object) -> None:
        """Initialise the primer not found error.

        Args:
            primer: The primer object that was not found.
        """
        super().__init__(f"Primer {primer} not found")
        self.primer = primer


class DuplicatedNameError(AmplifyPError, ValueError):
    """Exception raised when a primer name is already added."""

    def __init__(self, name: str) -> None:
        """Initialise the duplicated name error.

        Args:
            name: The primer name that was already added.
        """
        super().__init__(f"Primer name '{name}' already added")
        self.name = name


class DuplicatedSequenceError(AmplifyPError, ValueError):
    """Exception raised when a primer sequence is already added."""

    def __init__(self, sequence: str) -> None:
        """Initialise the duplicated sequence error.

        Args:
            sequence: The primer sequence that was already added.
        """
        super().__init__(f"Primer sequence '{sequence}' already added")
        self.sequence = sequence


class InvalidDNATypeError(AmplifyPError, TypeError):
    """Exception raised when a DNA type is not valid."""

    def __init__(self, message: str = "Invalid DNA type.") -> None:
        """Initialise the invalid DNA type error.

        Args:
            message: The error message describing the invalid DNA type.
        """
        super().__init__(message)


class InvalidDNASequenceError(AmplifyPError, ValueError):
    """Exception raised when a DNA sequence contains invalid characters."""

    def __init__(self, invalid_chars: set[str] | list[str]) -> None:
        """Initialise the invalid DNA sequence error.

        Args:
            invalid_chars: The set or list of invalid characters found in the
                DNA sequence.
        """
        sorted_chars = ", ".join(sorted(invalid_chars))
        super().__init__(
            f"The DNA sequence contains invalid characters: {sorted_chars}"
        )
        self.invalid_chars = invalid_chars


class ReplicationOriginLengthError(AmplifyPError, ValueError):
    """Exception raised when target and primer have mismatched lengths.

    This occurs when the length of target DNA is different from the primer DNA.
    """

    def __init__(
        self,
        target_len: int | None = None,
        primer_len: int | None = None,
        message: str = "The target has to have the same length as the primer.",
    ) -> None:
        """Initialise the replication origin length error.

        Args:
            target_len: Length of the target DNA, if available.
            primer_len: Length of the primer DNA, if available.
            message: The error message describing the length mismatch.
        """
        if target_len is not None and primer_len is not None:
            message = (
                f"{message} (Target length: {target_len}, "
                f"Primer length: {primer_len})"
            )
        super().__init__(message)
        self.target_len = target_len
        self.primer_len = primer_len


class BasePairTableDimensionError(AmplifyPError, ValueError):
    """Base exception for weight table dimension errors."""


class RowLengthMismatchError(BasePairTableDimensionError):
    """Exception raised when a weight table's row length does not match.

    The row length must match the expected size at initialisation.
    """

    def __init__(
        self,
        message: str = (
            "BasePairWeightsTbl: row length mismatch at initialisation."
        ),
    ) -> None:
        """Initialise the row length mismatch error.

        Args:
            message: The error message describing the row length mismatch.
        """
        super().__init__(message)


class ColumnLengthMismatchError(BasePairTableDimensionError):
    """Exception raised when a weight table's column length does not match.

    The column length must match the expected size at initialisation.
    """

    def __init__(
        self,
        message: str = (
            "BasePairWeightsTbl: column length mismatch at initialisation."
        ),
    ) -> None:
        """Initialise the column length mismatch error.

        Args:
            message: The error message describing the column length mismatch.
        """
        super().__init__(message)


class AmpliconDirectionError(AmplifyPError, ValueError):
    """Base exception for amplicon direction errors."""


class InvalidStartDirectionError(AmpliconDirectionError):
    """Exception raised when an amplicon start direction is not FWD."""

    def __init__(
        self, message: str = "Start direction must be forward."
    ) -> None:
        """Initialise the invalid start direction error.

        Args:
            message: The error message describing the invalid start direction.
        """
        super().__init__(message)


class InvalidEndDirectionError(AmpliconDirectionError):
    """Exception raised when an amplicon end direction is not REV."""

    def __init__(self, message: str = "End direction must be reverse.") -> None:
        """Initialise the invalid end direction error.

        Args:
            message: The error message describing the invalid end direction.
        """
        super().__init__(message)


class InvalidIndexOrderError(AmplifyPError, ValueError):
    """Exception raised when the end index precedes the start index.

    This is only checked on linear DNA templates.
    """

    def __init__(
        self,
        message: str = (
            "End index must be greater than start index for linear DNA."
        ),
    ) -> None:
        """Initialise the invalid index order error.

        Args:
            message: The error message describing the invalid index order.
        """
        super().__init__(message)


class TemplateMismatchError(AmplifyPError, ValueError):
    """Exception raised when a Repliconf template differs.

    This happens when the template of the Repliconf does not match the
    template of the AmpliconGenerator.
    """

    def __init__(
        self,
        message: str = (
            "The Repliconf contains a different template to the "
            "AmpliconGenerator."
        ),
    ) -> None:
        """Initialise the template mismatch error.

        Args:
            message: The error message describing the template mismatch.
        """
        super().__init__(message)


class InvalidAmpliconRangeError(AmplifyPError, NotImplementedError):
    """Exception raised when amplicon start index is greater than end index.

    This error is raised when attempting to search for an amplicon on a linear
    DNA template where start index is greater than end index.
    """

    def __init__(
        self,
        message: str = (
            "Attempted to search for an amplicon with the start index "
            "bigger than the end index on a linear DNA template."
        ),
    ) -> None:
        """Initialise the invalid amplicon range error.

        Args:
            message: The error message describing the invalid amplicon range.
        """
        super().__init__(message)


class InsufficientThermodynamicDataError(AmplifyPError, ValueError):
    """Exception raised when a sequence lacks thermodynamic data."""

    def __init__(
        self,
        message: str = "The sequence lacks thermodynamic data.",
    ) -> None:
        """Initialise the insufficient thermodynamic data error.

        Args:
            message: The error message describing the issue.
        """
        super().__init__(message)
