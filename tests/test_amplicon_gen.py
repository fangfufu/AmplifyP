"""Verify amplicon generation."""
import logging
from amplifyp.dna import DNA, Primer
from amplifyp.settings import Settings
from amplifyp.repliconf import Repliconf
from amplifyp.amplicon import AmpliconGenerator

def test_verify_amplicon_gen() -> None:
    """Run verification for amplicon generation."""
    # Template: AAAAA + (GT)*10 + GGGGG (30 bp)
    seq = "AAAAA" + "GT" * 10 + "GGGGG"
    template = DNA(seq, name="CleanTemplate")

    primer_fwd = Primer("AAAAA", "PF")
    primer_rev = Primer("CCCCC", "PR")

    settings = Settings()
    settings.primability_cutoff = 0.0
    settings.stability_cutoff = -100.0

    rc_fwd = Repliconf(template, primer_fwd, settings)
    rc_fwd.search()

    rc_rev = Repliconf(template, primer_rev, settings)
    rc_rev.search()

    gen = AmpliconGenerator(template)
    gen.add(rc_fwd)
    gen.add(rc_rev)

    amplicons = gen.generate_amplicons()
    logging.info(f"Generated {len(amplicons)} amplicons.")

    longest_amp = None
    max_len = 0

    for amp in amplicons:
        if amp.fwd_origin == primer_fwd and amp.rev_origin == primer_rev:
            if len(amp.sequence) > max_len:
                max_len = len(amp.sequence)
                longest_amp = amp

    assert longest_amp is not None, "No valid amplicon found"
    
    logging.info(f"Longest Amplicon Len: {max_len}")
    logging.info(f"Seq: {longest_amp.sequence.seq}")

    # Expected 30
    assert max_len == 30, f"Expected length 30, got {max_len}"
    assert longest_amp.sequence.seq == seq, "Sequence mismatch"
