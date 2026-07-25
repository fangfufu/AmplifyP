# AmplifyP GUI — Designer 2D View

The **Designer 2D View** (`Designer2DView`) performs paired two-sequence 2D
primer truncation analysis. It evaluates candidate forward and reverse DNA
sequences across progressive length truncations down to target minimum lengths
(`fwd_min_len`, `rev_min_len`), evaluating all pair combinations and computing
self-dimerisation and cross-dimerisation quality scores and overlap lengths.

For each pair combination step, 4 primer dimer alignments are evaluated:

1. **Forward Self-Dimer (`fwd_fwd`)**: Self-dimerisation potential of the
   forward primer.
2. **Reverse Self-Dimer (`rev_rev`)**: Self-dimerisation potential of the
   reverse primer.
3. **Forward-Reverse Cross-Dimer (`fwd_rev`)**: Cross-dimerisation alignment
   with the forward primer 3' end against the reverse primer.
4. **Reverse-Forward Cross-Dimer (`rev_fwd`)**: Cross-dimerisation alignment
   with the reverse primer 3' end against the forward primer.

The view features a two-panel split layout with interactive resizers
(`GestureDetector`), providing user control over split panel widths and heights.

## Layout & Panel Resizing

- **Two-Panel Split**:
  - **Left Column**: Contains the **2D Truncation Parameters** form (top-left)
    and the **2D Truncation Results Grid** matrix (bottom-left).
  - **Right Panel**: Contains the **2D Primer Pair Detail Cards** panel
    displaying dismissible detail cards for selected pair combinations.
- **Interactive Resizers**:
  - **Main Horizontal Divider**: Drag the vertical splitter bar
    (`main_h_divider`, minimum 250 px left width) to resize the left column
    relative to the right cards panel.
  - **Left Vertical Divider**: Drag the horizontal splitter bar
    (`left_v_divider`, minimum 150 px top height) to adjust the height of the
    top-left parameters form relative to the bottom-left results grid matrix.

## 2D Truncation Parameters (Top-Left Panel)

The **2D Truncation Parameters** form (`Designer2DForm`) configures the forward
and reverse candidate DNA sequences, minimum length constraints, dimer filtering
thresholds, and filter evaluation metric:

- **Forward Candidate Primer Sequence** (`fwd_dna_input`): Text field for
  inputting the forward candidate sequence. Raw sequence input is sanitised via
  `clean_sequence` to remove invalid characters. Must contain at least one valid
  base.
- **Fwd Min Len (bp)** (`fwd_min_len_input`): Minimum forward primer length
  integer (default: `18`, width 160). The forward sequence is truncated from the
  3' end base-by-base down to this minimum length. Must be a positive integer
  greater than 0 and cannot exceed the forward sequence length.
- **Reverse Candidate Primer Sequence** (`rev_dna_input`): Text field for
  inputting the reverse candidate sequence. Raw sequence input is sanitised via
  `clean_sequence`. Must contain at least one valid base.
- **Rev Min Len (bp)** (`rev_min_len_input`): Minimum reverse primer length
  integer (default: `18`, width 160). The reverse sequence is truncated from the
  5' end base-by-base down to this minimum length. Must be a positive integer
  greater than 0 and cannot exceed the reverse sequence length.
- **Quality Filter** (`quality_filter_input`): Maximum quality score threshold
  cutoff float. Defaults to the configured dimer quality threshold (`threshold`,
  e.g. `60.0`). Leave blank for unconstrained quality filtering. Pair
  combinations with quality scores exceeding this cutoff are excluded.
- **Overlap Filter (bp)** (`overlap_filter_input`): Maximum overlap length
  constraint integer. Leave blank for unconstrained overlap filtering. Pair
  combinations with overlap lengths exceeding this constraint are excluded.
- **Metric** (`filter_metric_dropdown`): Dropdown selector to choose how quality
  scores and overlap lengths are evaluated across the 4 dimer alignments for
  filtering and grid representation:
  - **Max (`MAX`)** (default): Evaluates pair filtering and grid values based on
    the maximum quality score / overlap length among all 4 dimer alignments.
  - **Mean (`MEAN`)**: Evaluates pair filtering and grid values based on the
    mean (average) quality score / overlap length across all 4 dimer alignments.
- **Analyse Button** (`analyse_button`): Triggers validation and runs the 2D
  primer truncation analysis (`_run_designer_event`). Can also be executed by
  pressing Enter inside any text input field.
- **Save / Load Parameters**:
  - **Save Button**: Saves the current form parameters (`fwd_dna`,
    `fwd_min_length`, `rev_dna`, `rev_min_length`, `quality_filter`,
    `overlap_filter`, `filter_metric`) to a YAML file
    (`designer_2d_parameters.yaml`) using a file save dialog
    (`save_and_write_file`).
  - **Load Button**: Opens a file picker dialog (`pick_and_read_file`) to import
    parameters from a `.yaml` or `.yml` file. Automatically populates input
    fields, clears previous errors, and executes analysis.

## 2D Truncation Results Grid (Bottom-Left Panel)

The **2D Truncation Results Grid** (`Grid2DResultsView`) presents all valid
forward-reverse primer truncation combinations in a colour-coded matrix grid:

- **Empty / No Match State**: Displays an italicised placeholder message prior
  to running an analysis, or an error message if no truncation combinations
  match active quality and overlap filters.
- **Grid Axes**:
  - **Columns**: Forward primer lengths (`fwd_lengths`) sorted in descending
    order (e.g. `24 bp`, `23 bp`, ...). Top-left header origin cell reads
    `Rev \ Fwd`.
  - **Rows**: Reverse primer lengths (`rev_lengths`) sorted in descending order.
- **Cell Representation & Quality Score**:
  - Each cell displays the evaluated quality score formatted to 1 decimal place
    (`{q_val:.1f}`), calculated using the selected metric (**Max** or **Mean**).
  - **Colour Mapping & Text Contrast**: Cell background colours are assigned
    dynamically via `designer_2d_colour` based on quality score. Text contrast
    colour (`get_text_contrast_colour`) automatically switches between dark and
    light text for optimal legibility.
  - **Colour Schemes**: Configurable via GUI Settings
    (`designer_2d_colour_scheme` setting in `Designer2DTile`): `None`,
    `Cool-Warm`, `Traffic Light`, or `Blue-Orange`.
- **Optimal Pair Highlighting**:
  - Cell(s) with the minimum quality score (optimal low dimerisation risk) are
    highlighted with a 2 px green border (`GUIColours.SUCCESS_GREEN`) and tagged
    with a star `★ Best Quality (Lowest Score)`.
- **Interactive Selection**:
  - Clicking a grid cell highlights it with a 3 px primary border
    (`GUIColours.PRIMARY`) and opens or raises its detailed pair card in the
    right panel.
- **Tooltips**:
  - Hovering over a cell displays details including Forward length, Reverse
    length, Max/Mean Quality, Max/Mean Overlap (bp), and Best Quality indicator
    if applicable.
- **Header Legend Badges**:
  - Displays summary badges above the grid for active metric
    (`Metric: Max/Mean Quality`), best quality score
    (`★ Best Quality: {min_q:.1f}`), and active colour scheme
    (`Colour Map: {scheme} ({max_q:.1f} - {min_q:.1f})`).
- **Scrollbar**: Supported with a top horizontal scrollbar
  (`ScrollbarOrientation.TOP`) for wide matrices.

## 2D Primer Pair Detail Cards (Right Panel)

The right panel displays detailed pair alignment cards (`Dismissible2DCard`) in
a vertical scrollable list (`ListView`):

- **Header Controls**:
  - **Panel Header**: Displays title **2D Primer Pair Detail Cards**.
  - **Clear Cards Button**: Appears (`visible=True`) when one or more cards are
    open. Clicking **Clear Cards** (`_clear_all_cards`) dismisses all open cards
    simultaneously.
- **Card Selection & Positioning**:
  - Selecting a cell in the results grid creates a dismissible detail card
    positioned at the top of the right panel list.
  - If a card for that pair combination already exists (identified by
    `card_2d_{fwd_len}_{rev_len}`), it is automatically raised to the top of the
    list.
- **Card Contents**:
  - **Card Header**: Displays title
    `2D Primer Pair (Forward: {fwd_len} bp, Reverse: {rev_len} bp)` alongside a
    close/dismiss button (`x`).
  - **Title Metric Badges**: Highlighted summary badges for
    `Max Quality: {score}`, `Mean Quality: {score}`,
    `Max Overlap: {overlap} bp`, and `Mean Overlap: {overlap} bp`.
  - **Primer Details Container**: Displays read-only sequence fields for both
    primers alongside individual copy buttons (`ft.Clipboard().set()`) and
    thermodynamic/composition badges:
    - **Forward Primer & Reverse Primer Sequence Fields**: Monospace text fields
      (`Roboto Mono`).
    - **Melting Temperature ($T_m$)**: `Tm: {val}°C` calculated using configured
      thermodynamic settings.
    - **AT Pairs**: `AT Pairs: {count}`
    - **GC Pairs**: `GC Pairs: {count}`
    - **% AT Content**: `% AT: {pct:.1f}%`
  - **4 Dimer Alignment Subcontainers**: Rendered in distinct bordered boxes for
    each of the 4 dimer alignments:
    1. **Forward Self-Dimer (Fwd-Fwd)**
    2. **Reverse Self-Dimer (Rev-Rev)**
    3. **Forward-Reverse Cross-Dimer (Fwd-Rev)**
    4. **Reverse-Forward Cross-Dimer (Rev-Fwd)**
    - Each subcontainer displays label header, metric badges
      (`Quality: {score}`, `Overlap: {length} bp`), and an antiparallel sequence
      alignment diagram (`create_overlapped_sequence_view`) rendered in
      monospace font (`Roboto Mono`).
- **Individual Dismissal**: Clicking the `x` button on an individual card
  removes it from the panel.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
