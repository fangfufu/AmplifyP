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

import pytest

from amplifyp.errors import DuplicatedPrimerError, PrimerNotFoundError
from amplifyp.pcr import PCR
from tests.examples.amplify4_examples import (
    amplify4_circular_example,
    amplify4_linear_example,
    fragment_177bp,
    fragment_190bp,
    fragment_220bp,
    fragment_566bp,
    fragment_660bp,
    fragment_1250bp,
    fragment_1280bp,
    fragment_1626bp,
    fragment_2033bp,
    fragment_2516bp,
    primer_11bp,
    primer_1701,
    primer_10289,
    primer_10290,
)


def test_pcr_circular() -> None:
    """Test PCR on a circular template."""
    pcr = PCR(amplify4_circular_example)
    pcr.add_primer(primer_11bp)
    pcr.add_primer(primer_1701)
    pcr.add_primer(primer_10289)
    pcr.add_primer(primer_10290)
    pcr.predict_amplicons()
    observed = [a.product for a in pcr.amplicons]
    observed.sort(key=len)

    expected = [
        fragment_177bp,
        fragment_190bp,
        fragment_220bp,
        fragment_566bp,
        fragment_660bp,
        fragment_1250bp,
        fragment_1280bp,
        fragment_1626bp,
        fragment_2033bp,
        fragment_2516bp,
    ]
    expected.sort(key=len)

    assert observed == expected


def test_pcr_linear() -> None:
    """Test PCR on a linear template."""
    pcr = PCR(amplify4_linear_example)
    pcr.add_primer(primer_11bp)
    pcr.add_primer(primer_1701)
    pcr.add_primer(primer_10289)
    pcr.add_primer(primer_10290)
    pcr.predict_amplicons()
    observed = [a.product for a in pcr.amplicons]
    observed.sort(key=len)

    expected = [
        fragment_177bp,
        fragment_190bp,
        fragment_220bp,
        fragment_566bp,
        fragment_660bp,
        fragment_2033bp,
        fragment_2516bp,
    ]
    expected.sort(key=len)

    assert observed == expected


def test_pcr_add_primer_duplicates() -> None:
    """Test that adding a duplicate primer raises DuplicatedPrimerError."""
    pcr = PCR(amplify4_linear_example)
    pcr.add_primer(primer_11bp)

    with pytest.raises(DuplicatedPrimerError):
        pcr.add_primer(primer_11bp)


def test_pcr_remove_primer() -> None:
    """Test removing a primer."""
    pcr = PCR(amplify4_linear_example)
    pcr.add_primer(primer_11bp)
    assert primer_11bp in pcr.primers

    pcr.remove_primer(primer_11bp)
    assert primer_11bp not in pcr.primers

    # Remove non-existent
    with pytest.raises(PrimerNotFoundError):
        pcr.remove_primer(primer_11bp)


def test_pcr_add_primers() -> None:
    """Test adding multiple primers at once."""
    pcr = PCR(amplify4_linear_example)
    primers = [primer_11bp, primer_1701]
    pcr.add_primers(primers)

    assert len(pcr.primers) == 2
    assert primer_11bp in pcr.primers
    assert primer_1701 in pcr.primers
