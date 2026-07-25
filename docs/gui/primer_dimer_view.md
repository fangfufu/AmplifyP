# AmplifyP GUI — Primer Dimers View

The **Primer Dimers View** evaluates potential secondary structure interactions
and hybridisation risks between all active primers configured in the
[Input View](input_view.md). It predicts both self-dimers (a primer binding to
another molecule of itself) and cross-dimers (two different primers binding to
each other) that could consume reagents or inhibit target amplification in PCR
reactions.

![Primer Dimers View](images/primer_dimer_view.png)

## Dimer Hybridisation Analysis

- **Pairwise Evaluation**: The algorithm checks all unique pairs of active
  primers ($\\frac{N(N+1)}{2}$ combinations for $N$ active primers, including
  self-dimer pairs).
- **Antiparallel 3' Alignment**: For each primer pair, the 3' end of the shorter
  primer (`primer_1`) is aligned against all candidate positions along the
  longer primer (`primer_2`). If both primers are of equal length and distinct,
  both directional orientations are evaluated and the higher-scoring alignment
  is selected.
- **Scoring & Filtering Criteria**:
  - **Quality Score**: Each base pair interaction in the overlapping region is
    scored using a pairwise nucleotide weighting matrix (default settings based
    on Amplify4).
  - **Minimum Overlap (`pd_min_overlap`)**: Only alignments meeting or exceeding
    the minimum overlap length (default: 3 bp) are retained.
  - **Quality Threshold (`pd_threshold`)**: Only dimers with a total quality
    score exceeding the minimum threshold (default: 60.0) are reported.
  - Both filtering parameters can be customised in [Settings](settings_view.md).
- **Automatic Result Caching**: Analysis results are cached automatically based
  on active primer sequences and dimer settings. Re-analysis is triggered only
  when primer data or dimer settings are modified.
- **Sorting**: Identified primer dimers are sorted in descending order of
  quality score, presenting the highest-risk interactions at the top of the
  list.

## Results Display & Performance Controls

- **Scrollable List**: Dimer results are presented in a continuous vertical
  scrollable list (`ListView`).
- **No Dimers Detected**: If no primer pairs pass the overlap and quality score
  thresholds, a centred notice is displayed:
  > *"No primer dimers detected above threshold."*
- **Safe Rendering Limit**: To prevent UI freezes when evaluating large primer
  sets with many secondary structure matches, rendering is capped at the top 100
  strongest binding dimers (`MAX_DIMERS_RENDER`). If more than 100 dimers are
  detected, a red warning notification is displayed at the top of the list:
  > [!WARNING]
  > *X primer dimers detected. Only the top 100 strongest binding dimers are
  > displayed to prevent UI freeze.*

## Primer Dimer Alignment Cards

Each detected dimer interaction is displayed in an individual card (`DimerCard`)
comprising a summary header and an antiparallel sequence alignment diagram:

- **Card Header**:
  - **Interaction Title**: Displays primer names as `{Primer Name} self-dimer`
    for self-dimer matches, or `{Primer 1} vs {Primer 2}` for cross-dimer
    matches. Title text is selectable.
  - **Quality Score Badge**: Displays `Quality: {score}` (rounded to 1 decimal
    place) in a highlighted badge. High quality scores indicate strong binding
    potential and high dimerisation risk.
  - **Overlap Length Badge**: Displays `Overlap: {overlap} bp` in a highlighted
    badge.
- **Visual Alignment Diagram**:
  - Rendered in a monospace font (`Roboto Mono` or configured font family)
    inside a bordered container with horizontal scrolling support.
  - Displays a 5-line visually aligned stack:
    1. **Top Primer Name**: Name of the longer primer (`primer_2`) rendered in
       bold purple.
    2. **Top Sequence ($5' \\to 3'$)**: Sequence of the longer primer displayed
       in standard $5' \\to 3'$ direction: `5'-{seq2}-3'`.
    3. **Binding Interface Line**: Visual bond symbols indicating interaction
       strength at each overlapping base position:
       - `|` (Vertical bar): Strong base pair match (pairwise score $\\ge
         10.0$).
       - `:` (Colon): Weak or ambiguous base pair match ($0.0 \\le \\text{score}
         < 10.0$).
       - ` ` (Space): Mismatch or repulsive interaction ($\\text{score} < 0.0$).
    4. **Bottom Sequence ($3' \\to 5'$)**: Sequence of the shorter primer
       displayed in antiparallel $3' \\to 5'$ orientation
       (`3'-{seq1[::-1]}-5'`), horizontally indented to match its alignment
       position.
    5. **Bottom Primer Name**: Name of the shorter primer (`primer_1`) rendered
       in bold purple, right-padded under its 5' end.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
