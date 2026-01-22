# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Test case clarifying primer binding logic.

This test serves as documentation and verification of how primers bind to
DNA templates in AmplifyP. It specifically addresses confusion about
directionality (Sense vs Antisense) and how searches are performed.
"""

from amplifyp.dna import DNA, Primer
from amplifyp.repliconf import Repliconf
from amplifyp.settings import Settings


def test_primer_binding_logic_clarification() -> None:
    """Verify that forward and reverse primers bind to the correct strands.

    In PCR and AmplifyP:
    1.  **Forward Primer**: Identical to the Sense (+) strand (5'->3').
        - It binds to the Antisense (-) strand (3'->5').
        - In AmplifyP, `search` with `DNADirection.FWD` looks for sequences
          identical to the primer on the Sense strand.

    2.  **Reverse Primer**: Reverse Complement of the Sense (+) strand.
        - It binds to the Sense (+) strand (5'->3').
        - In AmplifyP, `search` with `DNADirection.REV` looks for sequences
          identical to the primer on the Antisense strand (which corresponds
          to the Complement of the Sense strand).

    This test sets up a scenario where we know exactly where the primers should bind
    and asserts that `Repliconf` finds them in the expected `DNADirection` buckets.
    """
    # Template: 5'-AAAA ATGC AAAA-3' (Sense)
    # Sense sequence: "AAAAATGCAAAA"
    # Antisense sequence (3'->5'): "TTTTTACGTTTT"
    template_str = "AAAAATGCAAAA"
    template = DNA(template_str, name="TestTemplate")
    # Use strict cutoff to avoid partial matches (like 3/4 bases) confusing the test results
    # Note: Cutoff is exclusive (>), so we use 0.99 to allow 1.0 matches.
    settings = Settings(primability_cutoff=0.99)

    # --- Scenario 1: Forward Primer ---
    # Primer Sequence: 5'-ATGC-3'
    # This is identical to the middle part of the Sense strand.
    # Ideally binds to Antisense "TACG" (3'-5').
    # AmplifyP should find this in the FWD direction (Identity to Sense).
    primer_fwd = Primer("ATGC", "ForwardPrimer")
    repliconf_fwd = Repliconf(template, primer_fwd, settings)
    repliconf_fwd.search()

    # Assert we found exactly one binding site in FWD bucket
    assert len(repliconf_fwd.origin_db.fwd) == 1
    assert len(repliconf_fwd.origin_db.rev) == 0

    # Check the found origin
    origin_fwd = repliconf_fwd.origin(repliconf_fwd.origin_db.fwd[0])
    # The target found on the template should be "CGTA" (which is ATGC reversed)
    # because FWD search compares Primer(reversed) to Sense(reversed).
    # "ATGC" reversed is "CGTA".
    assert origin_fwd.target == "CGTA"
    assert origin_fwd.primability == 1.0

    # --- Scenario 2: Reverse Primer ---
    # Primer Sequence: 5'-GCAT-3'
    # This is the Reverse Complement of "ATGC" (Sense).
    # Ideally binds to Sense "ATGC" (5'-3').
    # AmplifyP should find this in the REV bucket.
    # REV search checks Identity to Antisense.
    # Antisense has "TACG" (reading 3'-5').
    # Primer "GCAT".
    # Wait. Antisense "TACG" corresponds to Sense "ATGC".
    # If Primer is "GCAT". Primer(rev) is "TACG".
    # Match "TACG" vs "TACG" (Antisense).
    primer_rev = Primer("GCAT", "ReversePrimer")
    repliconf_rev = Repliconf(template, primer_rev, settings)
    repliconf_rev.search()

    # Assert we found exactly one binding site in REV bucket
    assert len(repliconf_rev.origin_db.fwd) == 0
    assert len(repliconf_rev.origin_db.rev) == 1

    # Check the found origin
    origin_rev = repliconf_rev.origin(repliconf_rev.origin_db.rev[0])
    # The target found on the Antisense template should be "TACG".
    # This corresponds to the Antisense sequence segment complementary to "ATGC".
    assert origin_rev.target == "TACG"
    assert origin_rev.primability == 1.0
