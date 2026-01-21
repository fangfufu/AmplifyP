"""Amplify 4-specific tests for AmpliconGenerator."""

import logging

import tests.examples.amplify4_examples as ex
from amplifyp.amplicon import AmpliconGenerator
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS


def test_amplify4_examples() -> None:
    """Test amplicon generation using Amplify 4 examples."""
    # Create the repliconfs
    repliconf_11bp = Repliconf(
        ex.amplify4_linear_example, ex.primer_11bp, DEFAULT_SETTINGS
    )
    repliconf_1701 = Repliconf(
        ex.amplify4_linear_example, ex.primer_1701, DEFAULT_SETTINGS
    )
    repliconf_10289 = Repliconf(
        ex.amplify4_linear_example, ex.primer_10289, DEFAULT_SETTINGS
    )
    repliconf_10290 = Repliconf(
        ex.amplify4_linear_example, ex.primer_10290, DEFAULT_SETTINGS
    )

    amplicon_generator = AmpliconGenerator(ex.amplify4_linear_example)
    amplicon_generator.add(repliconf_11bp)
    amplicon_generator.add(repliconf_1701)
    amplicon_generator.add(repliconf_10289)
    amplicon_generator.add(repliconf_10290)

    amplicons = amplicon_generator.get_amplicons()
    for amplicon in amplicons:
        logging.info("Forward primer: %s", amplicon.fwd_origin)
        logging.info("Reverse primer: %s", amplicon.rev_origin)
        logging.info("Amplicon length: %s", len(amplicon.product.seq))
        logging.info("Amplicon quality: %s", amplicon.q_score)
        if len(amplicon.product.seq) == 660:
            assert amplicon.product.seq.upper() == ex.fragment_660bp.seq.upper()
            logging.info("Amplicon 660bp is correct")
        elif len(amplicon.product.seq) == 177:
            assert amplicon.product.seq.upper() == ex.fragment_177bp.seq.upper()
            logging.info("Amplicon 177bp is correct")
        elif len(amplicon.product.seq) == 220:
            assert amplicon.product.seq.upper() == ex.fragment_220bp.seq.upper()
            logging.info("Amplicon 220bp is correct")
        logging.info("Amplicon: %s", amplicon.product.seq)
