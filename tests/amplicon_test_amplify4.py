"""Amplify 4-specific tests for AmpliconGenerator."""

from amplifyp.amplicon import AmpliconGenerator
from amplifyp.dna import DNADirection
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS
from tests.examples.amplify4_examples import (
    amplify4_linear_example,
    fragment_177bp,
    fragment_220bp,
    fragment_566bp,
    fragment_660bp,
    fragment_2033bp,
    fragment_2516bp,
    primer_11bp,
    primer_1701,
    primer_10289,
    primer_10290,
)

repliconf_11bp = Repliconf(amplify4_linear_example, primer_11bp, DEFAULT_SETTINGS)
repliconf_1701 = Repliconf(amplify4_linear_example, primer_1701, DEFAULT_SETTINGS)
repliconf_10289 = Repliconf(amplify4_linear_example, primer_10289, DEFAULT_SETTINGS)
repliconf_10290 = Repliconf(amplify4_linear_example, primer_10290, DEFAULT_SETTINGS)


def test_amplify4_11bp_1701_10289_10290() -> None:
    """Test amplicon generation using Amplify 4 examples."""
    amplicon_generator = AmpliconGenerator(amplify4_linear_example)
    amplicon_generator.add(repliconf_11bp)
    amplicon_generator.add(repliconf_1701)
    amplicon_generator.add(repliconf_10289)
    amplicon_generator.add(repliconf_10290)

    amplicons = amplicon_generator.get_amplicons()
    encountered_660 = False
    encountered_177 = False
    encountered_220 = False
    for amplicon in amplicons:
        if len(amplicon.product.seq) == 660:
            assert amplicon.product.seq.upper() == fragment_660bp.seq.upper()
            assert amplicon.q_score == 629
            encountered_660 = True
        elif len(amplicon.product.seq) == 177:
            assert amplicon.product.seq.upper() == fragment_177bp.seq.upper()
            assert amplicon.q_score == 146
            encountered_177 = True
        elif len(amplicon.product.seq) == 220:
            assert amplicon.product.seq.upper() == fragment_220bp.seq.upper()
            assert amplicon.q_score == 189
            encountered_220 = True
    assert encountered_660
    assert encountered_177
    assert encountered_220


def test_amplify4_1701_10289_10290() -> None:
    """Test amplicon generation using Amplify 4 examples."""
    amplicon_generator = AmpliconGenerator(amplify4_linear_example)
    amplicon_generator.add(repliconf_1701)
    amplicon_generator.add(repliconf_10289)
    amplicon_generator.add(repliconf_10290)

    amplicons = amplicon_generator.get_amplicons()
    encountered_2033 = False
    encountered_2516 = False
    encountered_566 = False
    for amplicon in amplicons:
        if len(amplicon.product.seq) == 2033:
            assert amplicon.product.seq.upper() == fragment_2033bp.seq.upper()
            assert amplicon.q_score == 1993
            encountered_2033 = True
        elif len(amplicon.product.seq) == 2516:
            assert amplicon.product.seq.upper() == fragment_2516bp.seq.upper()
            assert amplicon.q_score == 2476
            encountered_2516 = True
        elif len(amplicon.product.seq) == 566:
            assert amplicon.product.seq.upper() == fragment_566bp.seq.upper()
            assert amplicon.q_score == 55696.52822187758
            encountered_566 = True
    assert encountered_2033
    assert encountered_2516
    assert encountered_566


def test_10290_10289_566bp_manual() -> None:
    """Test 10290_10289_566bp manually.

    This calculates the 566bp amplicon generated from the built-in template using
    the built-in primers with the label 10290 and 10289. The target sequence
    length is 526bp. The whole amplicon is 566bp.
    """
    repliconf_10290.search()
    repliconf_10289.search()
    assert repliconf_10290.origin_db.fwd == [710, 2177]
    origin_fwd = repliconf_10290.origin_from_db(DNADirection.FWD, 1)
    assert origin_fwd.primability == 0.8059701492537313
    assert origin_fwd.stability == 0.4717741935483871
    assert origin_fwd.quality == 0.09718042850264788
    origin_rev = repliconf_10289.origin_from_db(DNADirection.REV, 0)
    assert origin_rev.primability == 1.0
    assert origin_rev.stability == 1.0
    assert origin_rev.quality == 1.0
    fwd_quality = origin_fwd.quality
    rev_quality = origin_rev.quality
    assert 526 / (fwd_quality * rev_quality) ** 2 == 55696.52822187758
