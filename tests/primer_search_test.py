# -*- coding: utf-8 -*-
"""Test primer search using Amplify 4 built-in examples."""

import logging

from amplifyp.dna import DNADirection
from amplifyp.replication import Repliconf
from amplifyp.settings import Settings
from .amplify4_examples_test import (
    amplify4_linear_example,
    primer_20049,
    primer_10290,
    primer_2223,
)


def test_linear_search_primer_10290() -> None:
    """Test: search for primer 10290 in the Amplify 4 example template."""
    repliconf = Repliconf(amplify4_linear_example, primer_10290, Settings)
    repliconf.search()
    origin_id = repliconf.origin_id
    logging.info(origin_id)

    origin_1 = repliconf.origin(DNADirection.FWD, origin_id[DNADirection.FWD][0])
    assert origin_1.primability == 1.0
    assert origin_1.stability == 1.0
    assert origin_1.quality == 1.0

    origin_2 = repliconf.origin(DNADirection.FWD, origin_id[DNADirection.FWD][1])
    assert origin_2.primability == 0.8059701492537313
    assert origin_2.stability == 0.4717741935483871
    assert origin_2.quality == 0.09718042850264788


def test_linear_search_primer_20049() -> None:
    """Test: search for primer 20049 in the Amplify 4 example template."""
    repliconf = Repliconf(amplify4_linear_example, primer_20049, Settings)
    repliconf.search()
    origin_id = repliconf.origin_id
    logging.info(origin_id)

    origin_1 = repliconf.origin(DNADirection.FWD, origin_id[DNADirection.FWD][0])
    logging.info(origin_1)
    assert origin_1.primability == 0.8188405797101449
    assert origin_1.stability == 0.4838709677419355
    assert origin_1.quality == 0.12838943431510044

    origin_2 = repliconf.origin(DNADirection.REV, origin_id[DNADirection.REV][0])
    logging.info(origin_2)
    assert origin_2.primability == 0.855072463768116
    assert origin_2.stability == 0.9005376344086021
    assert origin_2.quality == 0.6945126227208976


def test_linear_search_primer_2223() -> None:
    """Test: search for primer 2223 in the Amplify 4 example template."""
    repliconf = Repliconf(amplify4_linear_example, primer_2223, Settings)
    repliconf.search()
    origin_id = repliconf.origin_id
    logging.info(origin_id)

    origin_1 = repliconf.origin(DNADirection.FWD, origin_id[DNADirection.FWD][0])
    logging.info(origin_1)
    assert origin_1.primability == 1.0
    assert origin_1.stability == 1.0
    assert origin_1.quality == 1.0

    origin_2 = repliconf.origin(DNADirection.REV, origin_id[DNADirection.REV][0])
    logging.info(origin_2)
    assert origin_2.primability == 0.8208955223880597
    assert origin_2.stability == 0.4306451612903226
    assert origin_2.quality == 0.06442585459797757

    origin_3 = repliconf.origin(DNADirection.REV, origin_id[DNADirection.REV][1])
    logging.info(origin_2)
    assert origin_3.primability == 0.8208955223880597
    assert origin_3.stability == 0.5419354838709678
    assert origin_3.quality == 0.2035387578237841
