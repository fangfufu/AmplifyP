"""Verify amplicon generation."""

import logging

# pylint: disable=duplicate-code
from amplifyp.dna import DNA, Primer
from amplifyp.settings import Settings
from amplifyp.repliconf import Repliconf
from amplifyp.amplicon import AmpliconGenerator

logging.basicConfig(level=logging.ERROR)


def verify() -> None:
    """Run verification."""
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
    print(f"Generated {len(amplicons)} amplicons.")

    longest_amp = None
    max_len = 0

    for amp in amplicons:
        if amp.fwd_origin == primer_fwd and amp.rev_origin == primer_rev:
            if len(amp.sequence) > max_len:
                max_len = len(amp.sequence)
                longest_amp = amp

    if longest_amp:
        print(f"Longest Amplicon Len: {max_len}")
        print(f"Seq: {longest_amp.sequence.seq}")
        # Expected 30
        if max_len != 30:
            raise AssertionError(f"Expected length 30, got {max_len}")
        if longest_amp.sequence.seq != seq:
            raise AssertionError("Sequence mismatch")
        print("Verification passed! Full length amplicon found.")
    else:
        raise AssertionError("No valid amplicon found")


if __name__ == "__main__":
    verify()
