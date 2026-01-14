"""Verify primer analysis logic."""

import os
import sys

# Add src to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

# pylint: disable=wrong-import-position

from amplifyp.dna import DNA, DNADirection, DNAType, Primer
from amplifyp.repliconf import Repliconf
from amplifyp.settings import Settings


def verify() -> None:
    """Run verification."""
    # 1. Setup
    template_seq = "GATCGATCGATCGATCGATC" * 5  # 100bp
    primer_seq = "GATC"  # match

    template = DNA(template_seq, DNAType.LINEAR, "Template")
    primer = Primer(primer_seq, "Primer1")

    print(f"Template: {template_seq}")
    print(f"Primer: {primer_seq}")

    # 2. Settings (0.0 cutoffs)
    settings = Settings()
    settings.primability_cutoff = 0.0
    settings.stability_cutoff = 0.0

    # 3. Search
    print("Searching...")
    rc = Repliconf(template, primer, settings)
    rc.search()

    # 4. Verify results
    print(f"Found {len(rc.origin_idx.fwd)} forward origins")
    print(f"Found {len(rc.origin_idx.rev)} reverse origins")

    found_any = False

    for direction in [DNADirection.FWD, DNADirection.REV]:
        indices = (
            rc.origin_idx.fwd if direction == DNADirection.FWD else rc.origin_idx.rev
        )
        d_str = "FWD" if direction == DNADirection.FWD else "REV"

        for i in indices:
            origin = rc.origin(direction, i)
            print(
                f"[{d_str} at {i}] Primability: {origin.primability:.4f}, "
                f"Stability: {origin.stability:.4f}, Quality: {origin.quality:.4f}"
            )
            found_any = True

            # Basic sanity check
            if origin.primability < 0 or origin.stability < 0:
                print("FAIL: Negative scores")
                sys.exit(1)

    if not found_any:
        print(
            "FAIL: No origins found with 0.0 cutoff (should find matches for "
            "GATC in GATC...)"
        )
        sys.exit(1)

    print("SUCCESS: Logic verified.")


if __name__ == "__main__":
    verify()
