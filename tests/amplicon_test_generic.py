"""Generic tests related to AmpliconGenerator."""

import pytest

from amplifyp.amplicon import AmpliconGenerator
from amplifyp.dna import DNA, Primer
from amplifyp.errors import DuplicateRepliconfError
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS


def test_duplicate_repliconf_error() -> None:
    """Test that adding a duplicate Repliconf raises DuplicateRepliconfError."""
    dna = DNA("ATGC" * 10, name="Template")
    primer = Primer("ATGC", "Primer1")
    repliconf = Repliconf(dna, primer, DEFAULT_SETTINGS)
    generator = AmpliconGenerator(dna)

    generator.add(repliconf)

    with pytest.raises(DuplicateRepliconfError):
        generator.add(repliconf)


def test_remove_repliconf() -> None:
    """Test that a Repliconf can be removed from the AmpliconGenerator."""
    dna = DNA("ATGC" * 10, name="Template")
    primer = Primer("ATGC", "Primer1")
    repliconf = Repliconf(dna, primer, DEFAULT_SETTINGS)
    generator = AmpliconGenerator(dna)

    generator.add(repliconf)
    assert repliconf in generator.repliconfs

    generator.remove(repliconf)
    assert repliconf not in generator.repliconfs


def test_remove_repliconf_not_found() -> None:
    """Test that removing a non-existent Repliconf raises ValueError."""
    dna = DNA("ATGC" * 10, name="Template")
    primer = Primer("ATGC", "Primer1")
    repliconf = Repliconf(dna, primer, DEFAULT_SETTINGS)
    generator = AmpliconGenerator(dna)

    with pytest.raises(ValueError):
        generator.remove(repliconf)
