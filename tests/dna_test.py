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

"""Tests related to the DNA class."""

import pytest

from amplifyp.dna import DNA, DNADirection, DNAType, Primer


def test_dna_init() -> None:
    """Test the initialization of a DNA object with a given sequence."""
    dna = DNA("ATCG")
    assert dna.seq == "ATCG"
    assert dna.type == DNAType.LINEAR
    assert dna.name == "ATCG"
    assert dna.direction == DNADirection.FWD


def test_dna_name_setter() -> None:
    """Test that the name of a DNA object can be set correctly."""
    dna = DNA("ATCG")
    dna.name = "test"
    assert dna.name == "test"


def test_dna_lower() -> None:
    """Test that the `lower` method of the `DNA` class."""
    dna = DNA("ATCG")
    assert dna.lower().seq == "atcg"


def test_dna_upper() -> None:
    """Test the `upper` method of the `DNA` class."""
    dna = DNA("atcg")
    assert dna.upper().seq == "ATCG"


def test_dna_complement() -> None:
    """Test the complement method of the DNA class."""
    dna = DNA("ATCG")
    assert dna.complement().seq == "TAGC"


def test_dna_reverse() -> None:
    """Test the reverse method of the DNA class."""
    dna = DNA("ATCG")
    assert dna.reverse().seq == "GCTA"


def test_dna_eq() -> None:
    """Test that the equality operator works for DNA objects."""
    dna1 = DNA("ATCG")
    dna2 = DNA("ATCG")
    assert dna1 == dna2
    assert dna1 != ""


def test_dna_len() -> None:
    """Test the length of a DNA sequence."""
    dna = DNA("ATCG")
    assert len(dna) == 4


def test_dna_getitem() -> None:
    """Test the __getitem__ method of the DNA class."""
    dna = DNA("ATCG")
    assert dna[1:3].seq == "TC"


def test_dna_str() -> None:
    """Test the string representation of the DNA class."""
    dna = DNA("ATCG")
    assert str(dna) == "DNA: ATCG, LINEAR, FWD"


def test_primer_init() -> None:
    """Test the initialization of a Primer object."""
    primer = Primer("ATCG")
    assert primer.seq == "ATCG"
    assert primer.type == DNAType.PRIMER
    assert primer.direction == DNADirection.FWD
    assert primer.name == "ATCG"


def test_dna_invalid_type() -> None:
    """Test DNA initialised with invalid type."""
    with pytest.raises(TypeError):
        DNA("A", 4)  # type: ignore[arg-type]


def test_dna_invalid_char() -> None:
    """Test DNA initialised with invalid characters."""
    with pytest.raises(ValueError):
        DNA("L")


def test_dna_hash() -> None:
    """Test that the DNA hash function generates a hash."""
    assert hash(DNA("A"))


def test_dna_count_at() -> None:
    """Test the count_at method."""
    dna = DNA("ATCGW", dna_type=DNAType.PRIMER)
    assert dna.count_at() == 3  # A, T, W
    assert (
        DNA("atcgw", dna_type=DNAType.PRIMER).count_at() == 3
    )  # Case insensitive
    assert DNA("CG", dna_type=DNAType.PRIMER).count_at() == 0


def test_dna_count_cg() -> None:
    """Test the count_cg method."""
    dna = DNA("ATCGS", dna_type=DNAType.PRIMER)
    assert dna.count_cg() == 3  # C, G, S
    assert (
        DNA("atcgs", dna_type=DNAType.PRIMER).count_cg() == 3
    )  # Case insensitive
    assert DNA("AT", dna_type=DNAType.PRIMER).count_cg() == 0


def test_dna_ratio_at() -> None:
    """Test the ratio_at method."""
    dna = DNA("AAAA")
    assert dna.ratio_at() == pytest.approx(1.0)
    dna = DNA("ATCG")
    assert dna.ratio_at() == pytest.approx(0.5)
    dna = DNA("")
    assert dna.ratio_at() == pytest.approx(0.0)


def test_dna_ratio_cg() -> None:
    """Test the ratio_cg method."""
    dna = DNA("CCCC")
    assert dna.ratio_cg() == pytest.approx(1.0)
    dna = DNA("ATCG")
    assert dna.ratio_cg() == pytest.approx(0.5)
    dna = DNA("")
    assert dna.ratio_cg() == pytest.approx(0.0)

def test_dna_reverse_complement_optimization() -> None:
    """Test reverse_complement correctness and consistency with reverse().complement()."""
    # Test 1: Standard sequence
    seq = "ATGC"
    dna = DNA(seq)
    rc = dna.reverse_complement()

    # Expected: Reverse of ATGC is CGTA. Complement of CGTA is GCAT.
    assert rc.seq == "GCAT"
    # Should be bool True (FWD equivalent) as per legacy behavior of reverse().complement()
    assert rc.direction is True
    assert rc.type == dna.type

    # Test 2: Verify consistency with reverse().complement() logic manually
    # reverse() -> NOT direction. complement() -> NOT (NOT direction) -> direction.
    assert rc.direction == bool(dna.direction)

    # Test 3: Ambiguous nucleotides
    # R (G/A) -> Y (C/T)
    # Y (C/T) -> R (G/A)
    # N -> N
    # Sequence: RYN
    # Reverse: NYR
    # Complement: N (N->N), Y (Y->R), R (R->Y) -> NRY
    seq_ambig = "RYN"
    dna_ambig = DNA(seq_ambig, dna_type=DNAType.PRIMER)
    rc_ambig = dna_ambig.reverse_complement()
    assert rc_ambig.seq == "NRY"

    # Test 4: Lowercase handling
    # "atgc" -> reverse "cgta" -> complement "gcat"
    dna_lower = DNA("atgc")
    assert dna_lower.reverse_complement().seq == "gcat"
