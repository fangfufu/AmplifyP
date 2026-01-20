"""Comprehensive test for amplicon using examples provided by jcww2."""

import tests.examples.jcww2_examples as ex
from amplifyp.amplicon import Amplicon
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
            amplicon = Amplicon(left_repliconf, right_repliconf, DEFAULT_SETTINGS)
            jcww2_test_amplicon[i][j] = amplicon

    return
