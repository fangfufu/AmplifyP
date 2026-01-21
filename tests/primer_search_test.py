"""Search potential replication origin using Amplify 4's built-in examples."""

from amplifyp.dna import DNADirection
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS

from .examples.amplify4_examples import (
    amplify4_linear_example,
    primer_2223,
    primer_10290,
    primer_20049,
)


def test_linear_search_primer_10290() -> None:
    """Test: search for primer 10290 in the Amplify 4 example template."""
    repliconf = Repliconf(amplify4_linear_example, primer_10290, DEFAULT_SETTINGS)
    repliconf.search()
    origin_idx = repliconf.origin_db

    origin_1 = repliconf.origin(origin_idx.fwd[0])
    assert origin_1.primability == 1.0
    assert origin_1.stability == 1.0
    assert origin_1.quality == 1.0

    origin_2 = repliconf.origin(origin_idx.fwd[1])
    assert origin_2.primability == 0.8059701492537313
    assert origin_2.stability == 0.4717741935483871
    assert origin_2.quality == 0.09718042850264788


def test_linear_search_primer_20049() -> None:
    """Test: search for primer 20049 in the Amplify 4 example template."""
    repliconf = Repliconf(amplify4_linear_example, primer_20049, DEFAULT_SETTINGS)
    repliconf.search()
    origin_index = repliconf.origin_db

    origin_1 = repliconf.origin(origin_index.fwd[0])
    assert origin_1.primability == 0.8188405797101449
    assert origin_1.stability == 0.4838709677419355
    assert origin_1.quality == 0.12838943431510044

    origin_2 = repliconf.origin(origin_index.rev[0])
    assert origin_2.primability == 0.855072463768116
    assert origin_2.stability == 0.9005376344086021
    assert origin_2.quality == 0.6945126227208976


def test_linear_search_primer_2223() -> None:
    """Test: search for primer 2223 in the Amplify 4 example template."""
    repliconf = Repliconf(amplify4_linear_example, primer_2223, DEFAULT_SETTINGS)
    repliconf.search()
    origin_index = repliconf.origin_db

    origin_1 = repliconf.origin(DNADirection.FWD, origin_index.fwd[0])
    assert origin_1.primability == 1.0
    assert origin_1.stability == 1.0
    assert origin_1.quality == 1.0

    origin_2 = repliconf.origin(origin_index.rev[0])
    assert origin_2.primability == 0.8208955223880597
    assert origin_2.stability == 0.4306451612903226
    assert origin_2.quality == 0.06442585459797757

    origin_3 = repliconf.origin(origin_index.rev[1])
    assert origin_3.primability == 0.8208955223880597
    assert origin_3.stability == 0.5419354838709678
    assert origin_3.quality == 0.2035387578237841
