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

"""Unit tests for custom exceptions in amplifyp.errors."""

from amplifyp.errors import (
    AmpliconDirectionError,
    AmplifyPError,
    BasePairTableDimensionError,
    ColumnLengthMismatchError,
    DuplicatedNameError,
    DuplicatedSequenceError,
    DuplicateRepliconfError,
    InsufficientThermodynamicDataError,
    InvalidAmpliconRangeError,
    InvalidDNASequenceError,
    InvalidDNATypeError,
    InvalidEndDirectionError,
    InvalidIndexOrderError,
    InvalidStartDirectionError,
    PrimerNotFoundError,
    ReplicationOriginLengthError,
    RowLengthMismatchError,
    TemplateMismatchError,
)


def test_base_exception_inheritance() -> None:
    """Test that all custom exceptions inherit from AmplifyPError."""
    exceptions = [
        DuplicateRepliconfError(),
        PrimerNotFoundError("test_primer"),
        DuplicatedNameError("name"),
        DuplicatedSequenceError("ACGT"),
        InvalidDNATypeError(),
        InvalidDNASequenceError({"X"}),
        ReplicationOriginLengthError(),
        RowLengthMismatchError(),
        ColumnLengthMismatchError(),
        InvalidStartDirectionError(),
        InvalidEndDirectionError(),
        InvalidIndexOrderError(),
        TemplateMismatchError(),
        InvalidAmpliconRangeError(),
        InsufficientThermodynamicDataError(),
    ]
    for exc in exceptions:
        assert isinstance(exc, AmplifyPError)
        assert isinstance(exc, Exception)


def test_exception_attributes() -> None:
    """Test that context attributes are stored properly on exceptions."""
    p_err = PrimerNotFoundError("primer_1")
    assert p_err.primer == "primer_1"

    n_err = DuplicatedNameError("name_1")
    assert n_err.name == "name_1"

    s_err = DuplicatedSequenceError("ACGT")
    assert s_err.sequence == "ACGT"

    seq_err = InvalidDNASequenceError({"Z", "X"})
    assert seq_err.invalid_chars == {"Z", "X"}

    rep_obj = object()
    dup_err = DuplicateRepliconfError(repliconf=rep_obj)
    assert dup_err.repliconf is rep_obj

    orig_err = ReplicationOriginLengthError(target_len=10, primer_len=15)
    assert orig_err.target_len == 10
    assert orig_err.primer_len == 15
    assert "(Target length: 10, Primer length: 15)" in str(orig_err)


def test_base_pair_table_dimension_error_hierarchy() -> None:
    """Test BasePairTableDimensionError hierarchy."""
    row_err = RowLengthMismatchError()
    col_err = ColumnLengthMismatchError()

    assert isinstance(row_err, BasePairTableDimensionError)
    assert isinstance(col_err, BasePairTableDimensionError)
    assert isinstance(row_err, AmplifyPError)
    assert isinstance(row_err, ValueError)


def test_amplicon_direction_error_hierarchy() -> None:
    """Test AmpliconDirectionError hierarchy."""
    start_err = InvalidStartDirectionError()
    end_err = InvalidEndDirectionError()

    assert isinstance(start_err, AmpliconDirectionError)
    assert isinstance(end_err, AmpliconDirectionError)
    assert isinstance(start_err, AmplifyPError)
    assert isinstance(start_err, ValueError)
