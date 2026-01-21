"""Debug script for Repliconf."""

import logging

from amplifyp.dna import DNA, DNADirection, Primer
from amplifyp.repliconf import Repliconf
from amplifyp.settings import Settings

# Setup logging to see debug output from Repliconf
logging.basicConfig(level=logging.DEBUG)

TEMPLATE_STR = "ATGC" * 10  # 40 bp
# ATGCATGC...
# 01234567

template = DNA(TEMPLATE_STR, name="TestTemplate")
primer_fwd = Primer("ATGC", "PF")  # ATGC -> binds to GCAT on template
# Template is ATGC... so GCAT is rev comp of ATGC.
# Wait.
# Fwd Primer "ATGC"
# Template 5'-...ATGC...-3'
# Primer binds to 3'-...TACG...-5' strand?
# PCR primers:
# Forward primer is identical to the Sense (+) strand (5'->3').
# It binds to the Antisense (-) strand (3'->5').
# Extension goes 5'->3' on the new strand, which matches the Sense strand.
# So Fwd Primer "ATGC" should match "ATGC" in the template sequence if template is
# Sense strand.
# "ATGC" in template -> FWD match.

primer_rev = Primer("GCAT", "PR")
# Reverse primer is Reverse Complement of the Sense strand.
# Sense: ...GCAT...
# Rev Primer: ATGC (if seq is ATGC)
# Wait, let's just see what Repliconf does.

settings = Settings()
settings.primability_cutoff = 0.0
settings.stability_cutoff = -100.0

print("-" * 20)
print("Forward Primer Test")
rc = Repliconf(template, primer_fwd, settings)
rc.search()

print(f"Template Seq (FWD padded): {rc.template_seq[DNADirection.FWD]}")
print(f"Origins FWD: {rc.origin_db.fwd}")
# amplicon_start calculation: [x - len(primer) for x in fwd]
print(f"Amplicon Start: {rc.amplicon_start}")

print("-" * 20)
print("Reverse Primer Test")
rc2 = Repliconf(template, primer_rev, settings)
rc2.search()
print(f"Origins REV: {rc2.origin_db.rev}")
# amplicon_end calculation: [x + len(primer) for x in rev]
print(f"Amplicon End: {rc2.amplicon_end}")
