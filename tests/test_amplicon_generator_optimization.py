import pytest
from amplifyp.dna import DNA, Primer
from amplifyp.amplicon import AmpliconGenerator
from amplifyp.repliconf import Repliconf
from amplifyp.errors import DuplicateRepliconfError

def test_internal_consistency():
    """Verify that repliconfs list and internal set are kept in sync."""
    dna = DNA("ATGC" * 10, name="Template")
    primer1 = Primer("ATGC", "Primer1")
    primer2 = Primer("GCAT", "Primer2")
    repliconf1 = Repliconf(dna, primer1)
    repliconf2 = Repliconf(dna, primer2)
    generator = AmpliconGenerator(dna)

    # Initial state
    assert len(generator.repliconfs) == 0
    assert len(generator._repliconfs_set) == 0

    # Add first repliconf
    generator.add_repliconf(repliconf1)
    assert len(generator.repliconfs) == 1
    assert len(generator._repliconfs_set) == 1
    assert repliconf1 in generator.repliconfs
    assert repliconf1 in generator._repliconfs_set

    # Add second repliconf
    generator.add_repliconf(repliconf2)
    assert len(generator.repliconfs) == 2
    assert len(generator._repliconfs_set) == 2
    assert repliconf2 in generator.repliconfs
    assert repliconf2 in generator._repliconfs_set

    # Try duplicate
    with pytest.raises(DuplicateRepliconfError):
        generator.add_repliconf(repliconf1)

    # State should remain unchanged
    assert len(generator.repliconfs) == 2
    assert len(generator._repliconfs_set) == 2

    # Remove repliconf1
    generator.remove_repliconf(repliconf1)
    assert len(generator.repliconfs) == 1
    assert len(generator._repliconfs_set) == 1
    assert repliconf1 not in generator.repliconfs
    assert repliconf1 not in generator._repliconfs_set
    assert repliconf2 in generator.repliconfs
    assert repliconf2 in generator._repliconfs_set

    # Remove repliconf2
    generator.remove_repliconf(repliconf2)
    assert len(generator.repliconfs) == 0
    assert len(generator._repliconfs_set) == 0

    # Try remove non-existent
    with pytest.raises(ValueError):
        generator.remove_repliconf(repliconf1)
