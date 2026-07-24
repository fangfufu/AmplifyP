# AmplifyP Python API Guide

This guide is for bioinformaticians, developers, and researchers who want to use
**AmplifyP** programmatically as a Python library.

The `amplifyp` library provides a fully typed, high-performance API to run
Polymerase Chain Reaction (PCR) simulations, analyse potential primer dimer
risks, and compute thermodynamic melting properties of DNA.

## 1. Defining DNA and Primers

The core of the library revolves around the `DNA` and `Primer` classes located
in `amplifyp.dna`.

### The `DNA` Class

The `DNA` class represents a template sequence. It supports different topologies
(linear or circular) and orientations.

```python
from amplifyp.dna import DNA, DNAType, DNADirection

# Define a linear DNA template
template_seq = "CATGATGAAATAACATAAGGTGGTCCCGTCGAAAGCCGAAGCGCAGAGCT"
linear_template = DNA(template_seq, DNAType.LINEAR, name="MyLinearTemplate")

# Define a circular DNA template (e.g. plasmids)
circular_template = DNA(template_seq, DNAType.CIRCULAR, name="MyPlasmid")

# Basic properties
print(f"Name: {linear_template.name}")
print(f"Sequence Length: {len(linear_template)} bp")
print(f"Topology: {linear_template.type.name}")  # LINEAR or CIRCULAR
```

### The `Primer` Class

The `Primer` class inherits from `DNA` but enforces a direction of 5' to 3' and
represents primer sequences, supporting IUPAC degenerate base codes.

```python
from amplifyp.dna import Primer

# Define primers (can include ambiguous codes like R, Y, M, S, W, K, V, H, D, B, N)
fwd_primer = Primer("CGACTGGGCAAAGGAAATCC", name="Fwd1701")
deg_primer = Primer("CGACTGGGYAAAGGAMATCC", name="DegeneratePrimer")

# Calculate redundancy fold (e.g. how many unique oligonucleotides are represented)
print(f"Redundancy fold: {deg_primer.redundancy_fold}")
# Count degenerate/redundant bases
print(f"Redundant base count: {deg_primer.redundant_base_count}")
```

## 2. PCR Simulation Engine

The `PCR` class in `amplifyp.pcr` orchestrates the simulation. It handles adding
primers, running the matching algorithm against the template, and predicting
potential amplicons.

```python
from amplifyp.dna import DNA, Primer, DNAType
from amplifyp.pcr import PCR

# 1. Define your DNA template and primers
template = DNA(
    "CATGATGAAATAACATAAGGTGGTCCCGTCGAAAGCCGAAGCGCAGAGCTGTCATACTCGAAGGCAGCAGCGAC"
    "CTTCATCTCGTCGAAAGCGAGTACGCAAAGCTTGTCGGCGTCATCAACTCCATCACTGTCCATTAGGTCTATGA"
    "CCACATCCAAACATCCTCTTTTTATGTCCACATCTGATAACCATCTGTACAAAGTCGTACGACTGGGCAAAGGA"
    "AATCCTTTTTTGTACAGATGGTTATACGCTCGAGGGCCTGCGGTGTGGAGACAAATAGCTGTAGAAATGTCGTC"
    "GGAATTGAACGTAGCTCTTTGTCCACCATTCTTCAGTATCCGTATCTGCGTGTCCGTGAAGATTTTGCGTAGAG"
    "ACTCCTCCAACTGTTGAGACTCCCTCAGCTGCTGCTCTAAACGACGCATTTCGTACTCCAAAGTACGAATTTTT"
    "TCCCTCAAGCTCTTATTTTCATTAAACAATGAACAGGACCTAACGCACAGTCACGTTATTGTTTACATAAATGA",
    DNAType.LINEAR,
    name="TestTemplate",
)

primer_fwd = Primer("CGACTGGGCAAAGGAAATCC", name="1701")
primer_rev = Primer("GTGGGTATCACAAATTTGGG", name="10290")

# 2. Initialise the PCR simulation
pcr = PCR(template)

# 3. Add primers to the reaction
pcr.add_primer(primer_fwd)
pcr.add_primer(primer_rev)

# 4. Predict amplicons (returns the count of predicted products)
num_amplicons = pcr.predict_amplicons()
print(f"Predicted {num_amplicons} amplicon(s).")

# 5. Access and print the amplicons
for amp in pcr.amplicons:
    print(f"Product Sequence: {amp.product.seq}")
    print(f"Product Length: {len(amp.product)} bp")
    print(f"Left Primer: {amp.origin_start.primer.name}")
    print(f"Right Primer: {amp.origin_end.primer.name}")
    print(f"Quality Score (Q): {amp.q_score:.2f}")
    print(
        f"Primability/Stability (Left): {amp.origin_start.primability:.1f}% / {amp.origin_start.stability:.1f}%"
    )
    print(
        f"Primability/Stability (Right): {amp.origin_end.primability:.1f}% / {amp.origin_end.stability:.1f}%"
    )
```

## 3. Primer Dimer Analysis

The `PrimerDimerGenerator` class in `amplifyp.dimer` identifies and scores the
potential of primers to hybridise with each other (self-dimer or cross-dimer
formation) and act as templates for DNA synthesis.

```python
from amplifyp.dna import Primer
from amplifyp.dimer import PrimerDimerGenerator

primer1 = Primer("CGACTGGGCAAAGGAAATCC", name="1701")
primer2 = Primer("GTGGGTATCACAAATTTGGG", name="10290")

# Initialise generator
dimer_gen = PrimerDimerGenerator()
dimer_gen.add_primer(primer1)
dimer_gen.add_primer(primer2)

# Run dimer analysis
dimer_gen.analyse_primers()

# Retrieve results (sorted by binding quality/strength)
for dimer in dimer_gen.primer_dimers:
    print(f"Dimer between: {dimer.primer_1.name} and {dimer.primer_2.name}")
    print(f"Overlap length: {dimer.overlap} bp")
    print(f"Quality Score: {dimer.quality:.1f}")
    print(f"Alignment Visual representation:")
    # Print the visual base pairing alignment
    # e.g., '|| |||' indicating match strength symbols
    print(f"  5' {dimer.primer_1.seq} 3'")
    print(f"     {dimer.binding_strength_str}")
    print(f"  3' {dimer.primer_2.seq[::-1]} 5'")
    print("-" * 30)
```

## 4. Thermodynamic Melting Temperature (Tm) Calculations

AmplifyP provides two models for computing primer melting temperatures in
`amplifyp.melting`:

1. **SantaLucia (1998)** with **Owczarzy (2008)** salt corrections (Default
   modern model).
2. **Lander-Amplify4** model (Legacy matrix-based model from original Amplify4).

```python
from amplifyp.dna import Primer
from amplifyp.melting import (
    calculate_tm_santalucia_1998_owczarzy_2008,
    calculate_tm_lander_amplify4,
)

primer = Primer("CGACTGGGCAAAGGAAATCC", name="FwdPrimer")

# 1. Calculate modern Tm (degrees Celsius)
tm_modern = calculate_tm_santalucia_1998_owczarzy_2008(primer)
print(f"Modern Tm: {tm_modern:.2f}°C")

# 2. Calculate legacy Amplify4 Tm (degrees Celsius)
tm_legacy = calculate_tm_lander_amplify4(primer)
print(f"Legacy Tm: {tm_legacy:.2f}°C")
```

## 5. Customising Simulation and Scoring Settings

The scoring thresholds, pairwise scores, and physical constants (like salt and
primer concentrations) are fully customisable.

```python
from amplifyp.settings import (
    ReplicationSettings,
    TMSettings,
    PrimerDimerSettings,
)

# 1. Customise PCR matching thresholds (stability and primability cutoffs)
custom_replication_settings = ReplicationSettings(
    primability_cutoff=85.0,  # Require 85% primability minimum
    stability_cutoff=80.0,  # Require 80% stability minimum
)

# 2. Customise thermodynamic parameters (concentration in M, e.g. 200 nM primers, 50 mM salt)
custom_tm_settings = TMSettings(
    oligo_concentration=2.0e-7,
    monovalent_concentration=0.05,
    divalent_concentration=0.0015,
)

# 3. Customise primer dimer thresholds
custom_dimer_settings = PrimerDimerSettings(
    threshold=40.0,  # Minimum dimer score to flag
    min_overlap=8,  # Minimum overlapping bp
)
```

You can pass these custom settings objects to class constructors:

- `PCR(template, settings=custom_replication_settings)`
- `PrimerDimerGenerator(settings=custom_dimer_settings)`
- `calculate_tm_santalucia_1998_owczarzy_2008(primer, settings=custom_tm_settings)`

## 6. Primer Design (1D and 2D Truncation Analysis)

AmplifyP includes primer design functionality to systematically evaluate dimer
formation across sequence truncations.

### 1D Primer Design (`PrimerDesigner1D`)

`PrimerDesigner1D` in `amplifyp.primer_designer_1d` performs 1D truncation
analysis by iteratively shortening a DNA sequence from either the 3' end
(`DNADirection.FWD`) or the 5' end (`DNADirection.REV`) down to a target minimum
length (`min_length`). At each step, self-dimer formation potential is
evaluated.

```python
from amplifyp.dna import DNA, DNADirection
from amplifyp.primer_designer_1d import PrimerDesigner1D

dna_template = DNA("CGACTGGGCAAAGGAAATCCGTGA", name="CandidatePrimer")

# Initialise 1D designer (truncating 3' end down to length 18)
designer_1d = PrimerDesigner1D(
    dna=dna_template,
    min_length=18,
    mode=DNADirection.FWD,
)

# Best truncation step (lowest quality score indicates minimal self-dimer potential)
best_step_idx, best_quality = designer_1d.best_score
best_dimer = designer_1d[best_step_idx]

print(
    f"Optimal length: {len(best_dimer.primer_1)} bp (Quality: {best_quality:.2f})"
)
print(f"Optimal sequence: {best_dimer.primer_1.seq}")

# Iterate over all truncation steps
for idx, dimer in enumerate(designer_1d.all_dimers):
    print(
        f"Step {idx}: length {len(dimer.primer_1)} bp, quality score = {dimer.quality:.2f}"
    )
```

### 2D Primer Design (`PrimerDesigner2D`)

`PrimerDesigner2D` in `amplifyp.primer_designer_2d` performs 2D truncation
analysis on candidate forward and reverse primer pairs. The forward sequence is
truncated from the 3' end while the reverse sequence is truncated from the 5'
end.

For each pair of truncated lengths down to `fwd_min_length` and
`rev_min_length`, 4 primer dimer configurations are calculated:

- Forward self-dimer (`fwd_fwd`)
- Reverse self-dimer (`rev_rev`)
- Forward-reverse cross-dimer (`fwd_rev`)
- Reverse-forward cross-dimer (`rev_fwd`)

Steps can be filtered using quality threshold limits and an evaluation metric
(`FilterMetric.MAX` or `FilterMetric.MEAN`).

```python
from amplifyp.dna import DNA
from amplifyp.primer_designer_2d import FilterMetric, PrimerDesigner2D

fwd_template = DNA("CGACTGGGCAAAGGAAATCCGTGA", name="FwdCandidate")
rev_template = DNA("GTGGGTATCACAAATTTGGGGACC", name="RevCandidate")

# Initialise 2D designer (truncating forward and reverse candidates)
designer_2d = PrimerDesigner2D(
    fwd_dna=fwd_template,
    fwd_min_length=18,
    rev_dna=rev_template,
    rev_min_length=18,
    filter_metric=FilterMetric.MAX,
)

# Identify best combined truncation step
best_step_idx, best_quality = designer_2d.best_score
best_step = designer_2d[best_step_idx]

print(
    f"Optimal step index: {best_step_idx} (Max Quality Score: {best_quality:.2f})"
)
print(f"Forward primer sequence: {best_step.fwd_fwd.primer_1.seq}")
print(f"Reverse primer sequence: {best_step.rev_rev.primer_1.seq}")
print(f"Cross-dimer quality (Fwd x Rev): {best_step.fwd_rev.quality:.2f}")
```
