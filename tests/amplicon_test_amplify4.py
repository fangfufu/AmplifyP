"""Amplify 4-specific tests for AmpliconGenerator."""

from amplifyp.amplicon import AmpliconGenerator
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
            # assert amplicon.q_score == 68702
            encountered_566 = True
    assert encountered_2033
    assert encountered_2516
    assert encountered_566
