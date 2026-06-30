# AmplifyP Python API Guide

This guide is for bioinformaticians, developers, and researchers who want to use
**AmplifyP** programmatically as a Python library.

The `amplifyp` library provides a fully typed, high-performance API to run
Polymerase Chain Reaction (PCR) simulations, analyse potential primer dimer
risks, and compute thermodynamic melting properties of DNA.

______________________________________________________________________

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

______________________________________________________________________

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
    name="TestTemplate"
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
    print(f"Primability/Stability (Left): {amp.origin_start.primability:.1f}% / {amp.origin_start.stability:.1f}%")
    print(f"Primability/Stability (Right): {amp.origin_end.primability:.1f}% / {amp.origin_end.stability:.1f}%")
```

______________________________________________________________________

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

______________________________________________________________________

## 4. Thermodynamic Melting Temperature (Tm) Calculations

AmplifyP provides two models for computing primer melting temperatures in
`amplifyp.melting`:

1. **SantaLucia (1998)** with **Owczarzy (2008)** salt corrections (Default
   modern model).
1. **Lander-Amplify4** model (Legacy matrix-based model from original Amplify4).

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

______________________________________________________________________

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
    stability_cutoff=80.0,    # Require 80% stability minimum
)

# 2. Customise thermodynamic parameters (concentration in M, e.g. 200 nM primers, 50 mM salt)
custom_tm_settings = TMSettings(
    oligo_concentration=2.0e-7,
    monovalent_concentration=0.05,
    divalent_concentration=0.0015,
)

# 3. Customise primer dimer thresholds
custom_dimer_settings = PrimerDimerSettings(
    threshold=40.0,      # Minimum dimer score to flag
    min_overlap=8,       # Minimum overlapping bp
)
```

You can pass these custom settings objects to class constructors:

- `PCR(template, settings=custom_replication_settings)`
- `PrimerDimerGenerator(settings=custom_dimer_settings)`
- `calculate_tm_santalucia_1998_owczarzy_2008(primer, settings=custom_tm_settings)`
