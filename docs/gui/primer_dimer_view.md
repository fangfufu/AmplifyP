# AmplifyP GUI — Primer Dimers View

Select the **Primer Dimers** tab to run a primer-primer hybridisation analysis.
AmplifyP evaluates all combinations of active primers (both self-dimers and
cross-dimers) to warn you about potential secondary structures.

![Primer Dimers View](images/primer_dimer_view.png)

- **Dimer Alignments**: Shows the optimal antiparallel alignment between the two
  primers (usually aligning the 3' end of the shorter primer with the longer
  primer).
- **Match Quality**: Dimers are sorted by quality score. A high score suggests a
  high risk of primer-dimer formation, which can poison your PCR reaction.
- **Visualisation**: Identifies complementary base pairs using `|` (perfect
  match) and `:` (weak/ambiguous match).

______________________________________________________________________

[Return to GUI Manual Index](README.md)
