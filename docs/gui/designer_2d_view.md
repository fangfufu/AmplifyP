# AmplifyP GUI — Designer 2D View

The **Designer 2D View** performs paired two-sequence 2D primer truncation
analysis. It evaluates candidate forward and reverse DNA sequences across
progressive length truncations down to specified target minimum lengths,
evaluating all pair combinations and computing self-dimerisation and
cross-dimerisation quality scores and overlap lengths.

For each pair combination step, 3 primer dimer alignments are evaluated by
default (a 4th is optional, see below):

1. **Forward Self-Dimer**: Self-dimerisation potential of the forward primer.
2. **Reverse Self-Dimer**: Self-dimerisation potential of the reverse primer.
3. **Forward-Reverse Cross-Dimer**: Cross-dimerisation alignment with the
   forward primer 3' end against the reverse primer.
4. **Reverse-Forward Cross-Dimer** (optional): Cross-dimerisation alignment with
   the reverse primer 3' end against the forward primer. Only evaluated and
   displayed when the **Show Reverse-Forward cross-dimer** checkbox in GUI
   Settings (Designer 2D tile) is enabled. It is disabled by default.

The view features a two-panel split layout with interactive resizers, providing
user control over split panel widths and heights.

## Layout & Panel Resizing

- **Two-Panel Split**:
  - **Left Column**: Contains the **2D Truncation Parameters** form (top-left)
    and the **2D Truncation Results Grid** matrix (bottom-left).
  - **Right Panel**: Contains the **2D Primer Pair Detail Cards** panel
    displaying dismissible detail cards for selected pair combinations.
- **Interactive Resizers**:
  - **Main Horizontal Divider**: Drag the vertical splitter bar (minimum 250 px
    left width) to resize the left column relative to the right cards panel.
  - **Left Vertical Divider**: Drag the horizontal splitter bar (minimum 110 px
    top height, default 380 px) to adjust the height of the top-left parameters
    form relative to the bottom-left results grid matrix.

## 2D Truncation Parameters (Top-Left Panel)

The **2D Truncation Parameters** form configures the forward and reverse
candidate DNA sequences, minimum length constraints, and dimer filtering
thresholds:

- **Forward Candidate Primer Sequence**: Text field for inputting the forward
  candidate sequence. Raw sequence input is automatically cleaned to remove
  invalid characters. Must contain at least one valid base.
- **Length (nt)**: Unmodifiable text field next to the forward sequence field
  displaying the current cleaned nucleotide length of the forward candidate
  primer sequence.
- **Fwd Min Length (nt)**: Minimum forward primer length. Required field (no
  default). The forward sequence is truncated from the 3' end base-by-base down
  to this minimum length. Must be a positive integer greater than 0 and cannot
  exceed the forward sequence length.
- **Reverse Candidate Primer Sequence**: Text field for inputting the reverse
  candidate sequence. Raw sequence input is automatically cleaned. Must contain
  at least one valid base.
- **Rev Comp**: Outlined button directly following the reverse sequence text
  field to immediately reverse-complement the reverse candidate sequence in
  place.
- **Length (nt)**: Unmodifiable text field next to the reverse sequence field
  displaying the current cleaned nucleotide length of the reverse candidate
  primer sequence.
- **Rev Min Length (nt)**: Minimum reverse primer length. Required field (no
  default). The reverse sequence is truncated from the 5' end base-by-base down
  to this minimum length. Must be a positive integer greater than 0 and cannot
  exceed the reverse sequence length.
- **Max Quality**: Upper bound quality cutoff. Non-negative integer (no
  default). Leave blank for unconstrained quality filtering. Pair combinations
  with quality scores exceeding this cutoff are excluded.
- **Max Overlap (bp)**: Upper bound overlap length cutoff in base pairs.
  Non-negative integer (no default). Leave blank for unconstrained overlap
  filtering. Pair combinations with overlap lengths exceeding this cutoff are
  excluded.
- **Check against template**: Bordered tickbox enabling evaluation of candidate
  primer pairs against the template DNA sequence from the Input view. When
  enabled, predicted amplicons are computed for each candidate pair and at least
  1 amplicon is enforced (pairs generating 0 amplicons are excluded). If no
  template sequence is present in the Input view, an error notification is
  displayed.
- **Max Amplicons**: Maximum allowed predicted amplicons on the template
  sequence. Enabled only when **Check against template** is selected. Leave
  blank for unconstrained amplicon evaluation (all pairs matching quality and
  overlap thresholds and having at least 1 amplicon are displayed alongside
  their predicted amplicon count). Must be a positive integer greater than 0;
  candidate pairs exceeding this limit are excluded.
- **Analyse Button**: Positioned on the right side of the bottom row. Triggers
  validation and runs the 2D primer truncation analysis. Can also be executed by
  pressing Enter inside any text input field.
- **Save / Load / Clear All Parameters**:
  - **Save Button**: Saves the current form parameters to a YAML file using a
    file save dialog.
  - **Load Button**: Opens a file picker dialog to import parameters from a
    `.yaml` or `.yml` file. Automatically populates input fields, clears
    previous errors, and executes analysis.
  - **Clear All Button**: Clears all form parameters and the analysis results
    (grid and cards).

## 2D Truncation Results Grid (Bottom-Left Panel)

The **2D Truncation Results Grid** presents all valid forward-reverse primer
truncation combinations in a colour-coded matrix grid. While an analysis is
running, the grid is replaced by a determinate **progress bar** labelled
*"Analysing primer combinations…"* that fills as each forward-reverse
combination is evaluated. The **Analyse** button changes to an **Abort** button
during analysis; clicking **Abort** cancels the running analysis and keeps the
results analysed so far, showing a notification: *"Analysis aborted — showing
results analysed so far."*

- **Empty / No Match State**: Displays an italicised placeholder message prior
  to running an analysis, or an error message if no truncation combinations
  match active quality and overlap filters.
- **Grid Axes**:
  - **Columns**: Forward primer lengths sorted in descending order (e.g.
    `24 nt`, `23 nt`, ...). Top-left header origin cell reads `Rev \ Fwd`.
  - **Rows**: Reverse primer lengths sorted in descending order.
- **Cell Representation & Quality Score**:
  - Each cell displays the maximum dimer quality score across the evaluated
    dimer alignments for that pair, rounded to the nearest integer (e.g. `42`).
  - **Colour Mapping & Text Contrast**: Cell background colours are assigned
    dynamically based on quality score. Text contrast colour automatically
    switches between dark and light text for optimal legibility.
  - **Colour Schemes**: Configurable via GUI Settings (under Designer 2D
    Settings): `None`, `Cool-Warm`, `Traffic Light`, `Blue-Orange` (default), or
    `Greyscale`.
- **Optimal Pair Highlighting**:
  - Cell(s) with the minimum quality score (optimal low dimerisation risk) are
    highlighted with a green border and tagged with a star
    `★ Best Quality (Lowest Score)`.
- **Interactive Selection**:
  - Clicking a grid cell highlights it with a primary border and opens or raises
    its detailed pair card in the right panel.
- **Tooltips**:
  - Hovering over a cell displays the `★ Best Quality (Lowest Score)` marker if
    applicable, the Forward and Reverse lengths (nt), `Max Quality`,
    `Max Overlap` (bp), and `Amplicons: {count}` when template evaluation was
    performed.
- **Header Legend Badges**:
  - Displays summary badges above the grid for the evaluation metric
    (`Metric: Max Quality`), best quality score (e.g. `★ Best Quality: 42`), and
    active colour scheme (e.g. `Colour Map: Blue-Orange (120 - 42)`).
- **Scrollbar**: Supported with a top horizontal scrollbar for wide matrices.

## 2D Primer Pair Detail Cards (Right Panel)

The right panel displays detailed pair alignment cards in a vertical scrollable
list:

- **Header Controls**:
  - **Panel Header**: Displays title **2D Primer Pair Detail Cards**.
  - **Clear Cards Button**: Appears when one or more cards are open. Clicking
    **Clear Cards** dismisses all open cards simultaneously.
- **Card Selection & Positioning**:
  - Selecting a cell in the results grid creates a dismissible detail card
    positioned at the top of the right panel list.
  - If a card for that pair combination already exists in the list, it is
    automatically raised to the top of the list.
- **Card Contents**:
  - **Card Header**: Displays title
    `2D Primer Pair (Forward: {fwd_len} nt, Reverse: {rev_len} nt)` alongside a
    close/dismiss button.
  - **Title Metric Badges**: Highlighted summary badges for
    `Max Quality: {score}`, `Mean Quality: {score}`,
    `Max Overlap: {overlap} bp`, and `Mean Overlap: {overlap} bp`. If template
    evaluation was performed, an `Amplicons: {count}` badge is also displayed.
  - **Run PCR Button**: Directly simulates PCR with the candidate forward and
    reverse primer pair against the template in the PCR view.
  - **Primer Details Container**: Displays read-only sequence fields for both
    primers alongside individual copy buttons and thermodynamic/composition
    badges:
    - **Forward Primer & Reverse Primer Sequence Fields**: Monospace text
      fields.
    - **Melting Temperature ($T_m$)**: `Tm: {value}°C` calculated using
      configured thermodynamic settings (shown as `Tm: N/A` when it cannot be
      calculated).
    - **% AT Content**: `% AT: {percentage}%`.
  - **Dimer Alignment Subcontainers**: Rendered in distinct bordered boxes for
    each evaluated dimer alignment:
    1. **Forward Self-Dimer (Fwd-Fwd)**
    2. **Reverse Self-Dimer (Rev-Rev)**
    3. **Forward-Reverse Cross-Dimer (Fwd-Rev)**
    4. **Reverse-Forward Cross-Dimer (Rev-Fwd)** — only when the **Show
       Reverse-Forward cross-dimer** setting is enabled.
    - Each subcontainer displays label header, metric badges
      (`Quality: {score}`, `Overlap: {length} bp`), and an antiparallel sequence
      alignment diagram rendered in a monospace font.
- **Individual Dismissal**: Clicking the close button on an individual card
  removes it from the panel.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
