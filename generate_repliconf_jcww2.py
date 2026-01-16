#!/usr/bin/env python3
"""Script to generate Repliconf objects using examples from jcww2_examples.py.

This script demonstrates how to create Repliconf instances with the provided
DNA templates and primers.
"""

import logging
import os
import sys

# Ensure project root is in sys.path
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.append(project_root)

# Also add src to path if project root has src dir layout
src_path = os.path.join(project_root, "src")
if src_path not in sys.path:
    sys.path.append(src_path)

# Add examples directory to path to import jcww2_examples
examples_path = os.path.join(project_root, "tests", "examples")
if examples_path not in sys.path:
    sys.path.append(examples_path)

import jcww2_examples as ex  # noqa: E402

from amplifyp.dna import DNADirection  # noqa: E402
from amplifyp.repliconf import Repliconf  # noqa: E402
from amplifyp.settings import DEFAULT_SETTINGS  # noqa: E402


def generate_repliconfs() -> None:
    """Generates and prints Repliconf objects for jcww2 examples."""
    # Configure logging to see debug output if needed
    logging.basicConfig(level=logging.ERROR, format="%(message)s")

    template = ex.jcww2_linear_example
    settings = DEFAULT_SETTINGS

    # List of primers from the example file
    primers = [
        ex.primer_jcww2_l1,
        ex.primer_jcww2_l2,
        ex.primer_jcww2_l3,
        ex.primer_jcww2_l4,
        ex.primer_jcww2_l5,
        ex.primer_jcww2_l6,
        ex.primer_jcww2_r1,
        ex.primer_jcww2_r2,
        ex.primer_jcww2_r3,
    ]

    print(f"Template Name: {template.name}")
    print(f"Template Length: {len(template)}")
    print("-" * 40)

    for i, primer in enumerate(primers):
        print(f"Processing Primer {i + 1}: {primer.name}")

        # Create Repliconf
        repliconf = Repliconf(template, primer, settings)

        # Search for origins
        repliconf.search()

        # Output results
        print(f"  Primer Sequence: {primer.seq}")
        print(f"  Forward Origins Found: {len(repliconf.origin_idx.fwd)}")
        print(f"  Reverse Origins Found: {len(repliconf.origin_idx.rev)}")

        if repliconf.origin_idx.fwd:
            print(f"  Forward Origin Indices: {repliconf.origin_idx.fwd}")
            for idx in repliconf.origin_idx.fwd:
                origin = repliconf.origin(DNADirection.FWD, idx)
                print(f"    Index {idx}:")
                print(f"      Primability: {origin.primability}")
                print(f"      Stability:   {origin.stability}")
                print(f"      Quality:     {origin.quality}")

        if repliconf.origin_idx.rev:
            print(f"  Reverse Origin Indices: {repliconf.origin_idx.rev}")
            for idx in repliconf.origin_idx.rev:
                origin = repliconf.origin(DNADirection.REV, idx)
                print(f"    Index {idx}:")
                print(f"      Primability: {origin.primability}")
                print(f"      Stability:   {origin.stability}")
                print(f"      Quality:     {origin.quality}")
        print("-" * 40)


if __name__ == "__main__":
    generate_repliconfs()
