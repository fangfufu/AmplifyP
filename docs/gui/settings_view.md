# AmplifyP GUI — Settings & Preferences View

You can customise the matching and dimer algorithms. To open settings, click the
gear icon or settings button in the GUI:

- **General**: Customise default directories, fonts, and state persistence.
- **PCR Match Settings**:
  - **Maximum Effective Primer**: Defaults to 30 bp. If a primer is longer, only
    the 30 bases from the 3' end are evaluated for matches (though the entire
    primer is amplified).
  - **Cutoffs**: Set the minimum **Primability Cutoff** and **Stability Cutoff**
    percentages. Matches below these thresholds are ignored.
  - **Pairwise Weights**: Customise base pairing scores and weighting factors
    for runs of matches.
- **Tm Calculations**: Configure monovalent and divalent salt concentrations,
  primer concentrations, and thermodynamic parameters (SantaLucia vs.
  Lander-Amplify4).
- **Primer Dimers**: Customise weights and minimum scoring thresholds for dimer
  detection.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
