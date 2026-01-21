"""Comprehensive tests for Repliconf using jcww2 examples."""

import tests.examples.jcww2_examples as ex

from amplifyp.dna import DNADirection
from amplifyp.repliconf import Repliconf
from amplifyp.settings import DEFAULT_SETTINGS


def test_repliconf_primer_jcww2_l1() -> None:
    """Test repliconf with primer_jcww2_l1 (GFPxVP249TPGF)."""
    template = ex.jcww2_template
    primer = ex.jcww2_left_primer[0]
    settings = DEFAULT_SETTINGS

    repliconf = Repliconf(template, primer, settings)
    repliconf.search()

    # Assertions based on generated output
    assert len(repliconf.origin_db.fwd) == 1
    assert len(repliconf.origin_db.rev) == 0
    assert repliconf.origin_db.fwd == [4482]

    origin = repliconf.origin(DNADirection.FWD, 4482)
    assert origin.primability == 1.0
    assert origin.stability == 1.0
    assert origin.quality == 1.0


def test_repliconf_primer_jcww2_l2() -> None:
    """Test repliconf with primer_jcww2_l2 (VP249GTd87_GF)."""
    template = ex.jcww2_template
    primer = ex.jcww2_left_primer[1]
    settings = DEFAULT_SETTINGS

    repliconf = Repliconf(template, primer, settings)
    repliconf.search()

    # Assertions based on generated output
    assert len(repliconf.origin_db.fwd) == 1
    assert len(repliconf.origin_db.rev) == 0
    assert repliconf.origin_db.fwd == [4743]

    origin = repliconf.origin(DNADirection.FWD, 4743)
    assert origin.primability == 0.9137931034482759
    assert origin.stability == 0.7281362007168459
    assert origin.quality == 0.552411630206402


def test_repliconf_primer_jcww2_l3() -> None:
    """Test repliconf with primer_jcww2_l3 (VP249GTd68_GF)."""
    template = ex.jcww2_template
    primer = ex.jcww2_left_primer[2]
    settings = DEFAULT_SETTINGS

    repliconf = Repliconf(template, primer, settings)
    repliconf.search()

    # Assertions based on generated output
    assert len(repliconf.origin_db.fwd) == 1
    assert len(repliconf.origin_db.rev) == 0
    assert repliconf.origin_db.fwd == [4686]

    origin = repliconf.origin(DNADirection.FWD, 4686)
    assert origin.primability == 0.9022988505747126
    assert origin.stability == 0.6935483870967742
    assert origin.quality == 0.4948090470893586


def test_repliconf_primer_jcww2_l4() -> None:
    """Test repliconf with primer_jcww2_l4 (VP249GTd41_GF)."""
    template = ex.jcww2_template
    primer = ex.jcww2_left_primer[3]
    settings = DEFAULT_SETTINGS

    repliconf = Repliconf(template, primer, settings)
    repliconf.search()

    # Assertions based on generated output
    assert len(repliconf.origin_db.fwd) == 1
    assert len(repliconf.origin_db.rev) == 0
    assert repliconf.origin_db.fwd == [4604]

    origin = repliconf.origin(DNADirection.FWD, 4604)
    assert origin.primability == 0.9421965317919075
    assert origin.stability == 0.8002551485328959
    assert origin.quality == 0.6780646004060044


def test_repliconf_primer_jcww2_l5() -> None:
    """Test repliconf with primer_jcww2_l5 (VP249GTd32_GF)."""
    template = ex.jcww2_template
    primer = ex.jcww2_left_primer[4]
    settings = DEFAULT_SETTINGS

    repliconf = Repliconf(template, primer, settings)
    repliconf.search()

    # Assertions based on generated output
    assert len(repliconf.origin_db.fwd) == 1
    assert len(repliconf.origin_db.rev) == 0
    assert repliconf.origin_db.fwd == [4583]

    origin = repliconf.origin(DNADirection.FWD, 4583)
    assert origin.primability == 0.9329608938547486
    assert origin.stability == 0.8039702233250621
    assert origin.quality == 0.671163896474763


def test_repliconf_primer_jcww2_l6() -> None:
    """Test repliconf with primer_jcww2_l6 (VP249GTd21_GF)."""
    template = ex.jcww2_template
    primer = ex.jcww2_left_primer[5]
    settings = DEFAULT_SETTINGS

    repliconf = Repliconf(template, primer, settings)
    repliconf.search()

    # Assertions based on generated output
    assert len(repliconf.origin_db.fwd) == 1
    assert len(repliconf.origin_db.rev) == 0
    assert repliconf.origin_db.fwd == [4545]

    origin = repliconf.origin(DNADirection.FWD, 4545)
    assert origin.primability == 0.9022988505747126
    assert origin.stability == 0.7012544802867383
    assert origin.quality == 0.5044416635768135


def test_repliconf_primer_jcww2_r1() -> None:
    """Test repliconf with primer_jcww2_r1 (TetRSTR)."""
    template = ex.jcww2_template
    primer = ex.jcww2_right_primer[0]
    settings = DEFAULT_SETTINGS

    repliconf = Repliconf(template, primer, settings)
    repliconf.search()

    # Assertions based on generated output
    assert len(repliconf.origin_db.fwd) == 0
    assert len(repliconf.origin_db.rev) == 1
    assert repliconf.origin_db.rev == [1616]

    origin = repliconf.origin(DNADirection.REV, 1616)
    assert origin.primability == 0.9937106918238994
    assert origin.stability == 0.9777777777777777
    assert origin.quality == 0.9643605870020963


def test_repliconf_primer_jcww2_r2() -> None:
    """Test repliconf with primer_jcww2_r2 (VP249GTd12C_GR)."""
    template = ex.jcww2_template
    primer = ex.jcww2_right_primer[1]
    settings = DEFAULT_SETTINGS

    repliconf = Repliconf(template, primer, settings)
    repliconf.search()

    # Assertions based on generated output
    assert len(repliconf.origin_db.fwd) == 0
    assert len(repliconf.origin_db.rev) == 1
    assert repliconf.origin_db.rev == [5071]

    origin = repliconf.origin(DNADirection.REV, 5071)
    assert origin.primability == 0.9080459770114943
    assert origin.stability == 0.7025089605734767
    assert origin.quality == 0.5131936719812137


def test_repliconf_primer_jcww2_r3() -> None:
    """Test repliconf with primer_jcww2_r3 (VP249GTd24C_GR)."""
    template = ex.jcww2_template
    primer = ex.jcww2_right_primer[2]
    settings = DEFAULT_SETTINGS

    repliconf = Repliconf(template, primer, settings)
    repliconf.search()

    # Assertions based on generated output
    assert len(repliconf.origin_db.fwd) == 0
    assert len(repliconf.origin_db.rev) == 1
    assert repliconf.origin_db.rev == [5031]

    origin = repliconf.origin(DNADirection.REV, 5031)
    assert origin.primability == 0.9157303370786517
    assert origin.stability == 0.7367271505376344
    assert origin.quality == 0.5655718595203575
