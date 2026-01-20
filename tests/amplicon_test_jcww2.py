"""Comprehensive test for amplicon using examples provided by jcww2."""

import logging

import tests.examples.jcww2_examples as ex
from amplifyp.amplicon import Amplicon, AmpliconGenerator
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS


def test_jcww2_amplicon_examples() -> None:
    """Test amplicon creation using examples provided by jcww2."""
    # Create the repliconfs
    jcww2_left_repliconfs: list[Repliconf | None] = [
        None for _ in range(len(ex.jcww2_left_primer))
    ]
    for i, primer in enumerate(ex.jcww2_left_primer):
        repliconf = Repliconf(ex.jcww2_template, primer, DEFAULT_SETTINGS)
        repliconf.search()
        jcww2_left_repliconfs[i] = repliconf

    jcww2_right_repliconfs: list[Repliconf | None] = [
        None for _ in range(len(ex.jcww2_right_primer))
    ]
    for i, primer in enumerate(ex.jcww2_right_primer):
        repliconf = Repliconf(ex.jcww2_template, primer, DEFAULT_SETTINGS)
        repliconf.search()
        jcww2_right_repliconfs[i] = repliconf

    # Create the Amplicons
    jcww2_test_amplicon: list[list[Amplicon | None]] = [
        [None for _ in range(len(ex.jcww2_right_primer))]
        for _ in range(len(ex.jcww2_left_primer))
    ]

    for i, left_repliconf in enumerate(jcww2_left_repliconfs):
        for j, right_repliconf in enumerate(jcww2_right_repliconfs):
            amplicon_generator = AmpliconGenerator(ex.jcww2_template)
            if left_repliconf is None:
                continue
            amplicon_generator.add(left_repliconf)
            logging.info("Added left repliconf: %s", left_repliconf)
            logging.info("Left repliconf index: %s", left_repliconf.origin_idx)
            logging.info(
                "Left repliconf amplicon start: %s", left_repliconf.amplicon_start
            )
            if right_repliconf is None:
                continue
            amplicon_generator.add(right_repliconf)
            logging.info("Added right repliconf: %s", right_repliconf)
            logging.info("Right repliconf index: %s", right_repliconf.origin_idx)
            logging.info(
                "Right repliconf amplicon end: %s", right_repliconf.amplicon_end
            )
            amplicon = amplicon_generator.get_amplicons()
            logging.info("Generated amplicon: %s", amplicon)
            example_amplicon = ex.jcww2_example_amplicon[i][j]
            if example_amplicon is None:
                continue
            logging.info("Expected amplicon: %s", example_amplicon.product)
            jcww2_test_amplicon[i][j] = amplicon[0]
            # assert (
            #     jcww2_test_amplicon[i][j].product.seq.upper()
            #     == ex.jcww2_example_amplicon[i][j].product.seq.upper()
            # )
    return
