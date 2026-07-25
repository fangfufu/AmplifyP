# AmplifyP GUI — Designer 2D View

The **Designer 2D View** performs paired two-sequence 2D primer truncation
analysis. It evaluates candidate forward and reverse DNA sequences across
progressive length truncations down to specified target minimum lengths,
evaluating all pair combinations and computing self-dimerisation and
cross-dimerisation quality scores and overlap lengths.

For each pair combination step, 4 primer dimer alignments are evaluated:

1. **Forward Self-Dimer**: Self-dimerisation potential of the forward primer.
2. **Reverse Self-Dimer**: Self-dimerisation potential of the reverse primer.
3. **Forward-Reverse Cross-Dimer**: Cross-dimerisation alignment with the
   forward primer 3' end against the reverse primer.
4. **Reverse-Forward Cross-Dimer**: Cross-dimerisation alignment with the
   reverse primer 3' end against the forward primer.

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
  - **Left Vertical Divider**: Drag the horizontal splitter bar (minimum 150 px
    top height) to adjust the height of the top-left parameters form relative to
    the bottom-left results grid matrix.

## 2D Truncation Parameters (Top-Left Panel)

The **2D Truncation Parameters** form configures the forward and reverse
candidate DNA sequences, minimum length constraints, dimer filtering thresholds,
and filter evaluation metric:

- **Forward Candidate Primer Sequence**: Text field for inputting the forward
  candidate sequence. Raw sequence input is automatically cleaned to remove
  invalid characters. Must contain at least one valid base.
- **Fwd Min Len (bp)**: Minimum forward primer length integer (default: `18`).
  The forward sequence is truncated from the 3' end base-by-base down to this
  minimum length. Must be a positive integer greater than 0 and cannot exceed
  the forward sequence length.
- **Reverse Candidate Primer Sequence**: Text field for inputting the reverse
  candidate sequence. Raw sequence input is automatically cleaned. Must contain
  at least one valid base.
- **Rev Min Len (bp)**: Minimum reverse primer length integer (default: `18`).
  The reverse sequence is truncated from the 5' end base-by-base down to this
  minimum length. Must be a positive integer greater than 0 and cannot exceed
  the reverse sequence length.
- **Quality Filter**: Maximum quality score threshold cutoff. Defaults to the
  configured dimer quality threshold (e.g. `60.0`). Leave blank for
  unconstrained quality filtering. Pair combinations with quality scores
  exceeding this cutoff are excluded.
- **Overlap Filter (bp)**: Maximum overlap length constraint integer. Leave
  blank for unconstrained overlap filtering. Pair combinations with overlap
  lengths exceeding this constraint are excluded.
- **Metric**: Dropdown selector to choose how quality scores and overlap lengths
  are evaluated across the 4 dimer alignments for filtering and grid
  representation:
  - **Max** (default): Evaluates pair filtering and grid values based on the
    maximum quality score / overlap length among all 4 dimer alignments.
  - **Mean**: Evaluates pair filtering and grid values based on the mean
    (average) quality score / overlap length across all 4 dimer alignments.
- **Analyse Button**: Triggers validation and runs the 2D primer truncation
  analysis. Can also be executed by pressing Enter inside any text input field.
- **Save / Load Parameters**:
  - **Save Button**: Saves the current form parameters to a YAML file using a
    file save dialog.
  - **Load Button**: Opens a file picker dialog to import parameters from a
    `.yaml` or `.yml` file. Automatically populates input fields, clears
    previous errors, and executes analysis.

## 2D Truncation Results Grid (Bottom-Left Panel)

The **2D Truncation Results Grid** presents all valid forward-reverse primer
truncation combinations in a colour-coded matrix grid:

- **Empty / No Match State**: Displays an italicised placeholder message prior
  to running an analysis, or an error message if no truncation combinations
  match active quality and overlap filters.
- **Grid Axes**:
  - **Columns**: Forward primer lengths sorted in descending order (e.g.
    `24 bp`, `23 bp`, ...). Top-left header origin cell reads `Rev \ Fwd`.
  - **Rows**: Reverse primer lengths sorted in descending order.
- **Cell Representation & Quality Score**:
  - Each cell displays the evaluated quality score formatted to 1 decimal place
    (e.g. `42.0`), calculated using the selected metric (**Max** or **Mean**).
  - **Colour Mapping & Text Contrast**: Cell background colours are assigned
    dynamically based on quality score. Text contrast colour automatically
    switches between dark and light text for optimal legibility.
  - **Colour Schemes**: Configurable via GUI Settings (under Designer 2D
    Settings): `None`, `Cool-Warm`, `Traffic Light`, or `Blue-Orange`.
- **Optimal Pair Highlighting**:
  - Cell(s) with the minimum quality score (optimal low dimerisation risk) are
    highlighted with a green border and tagged with a star
    `★ Best Quality (Lowest Score)`.
- **Interactive Selection**:
  - Clicking a grid cell highlights it with a primary border and opens or raises
    its detailed pair card in the right panel.
- **Tooltips**:
  - Hovering over a cell displays details including Forward length, Reverse
    length, Max/Mean Quality, Max/Mean Overlap (bp), and Best Quality indicator
    if applicable.
- **Header Legend Badges**:
  - Displays summary badges above the grid for active metric
    (`Metric: Max/Mean Quality`), best quality score (e.g.
    `★ Best Quality: 42.0`), and active colour scheme (e.g.
    `Colour Map: Blue-Orange (120.0 - 42.0)`).
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
    `2D Primer Pair (Forward: {fwd_len} bp, Reverse: {rev_len} bp)` alongside a
    close/dismiss button.
  - **Title Metric Badges**: Highlighted summary badges for
    `Max Quality: {score}`, `Mean Quality: {score}`,
    `Max Overlap: {overlap} bp`, and `Mean Overlap: {overlap} bp`.
  - **Primer Details Container**: Displays read-only sequence fields for both
    primers alongside individual copy buttons and thermodynamic/composition
    badges:
    - **Forward Primer & Reverse Primer Sequence Fields**: Monospace text
      fields.
    - **Melting Temperature ($T_m$)**: `Tm: {value}°C` calculated using
      configured thermodynamic settings.
    - **AT Pairs**: `AT Pairs: {count}`
    - **GC Pairs**: `GC Pairs: {count}`
    - **% AT Content**: `% AT: {percentage}%`
  - **4 Dimer Alignment Subcontainers**: Rendered in distinct bordered boxes for
    each of the 4 dimer alignments:
    1. **Forward Self-Dimer (Fwd-Fwd)**
    2. **Reverse Self-Dimer (Rev-Rev)**
    3. **Forward-Reverse Cross-Dimer (Fwd-Rev)**
    4. **Reverse-Forward Cross-Dimer (Rev-Fwd)**
    - Each subcontainer displays label header, metric badges
      (`Quality: {score}`, `Overlap: {length} bp`), and an antiparallel sequence
      alignment diagram rendered in a monospace font.
- **Individual Dismissal**: Clicking the close button on an individual card
  removes it from the panel.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
