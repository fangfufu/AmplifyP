# AmplifyP GUI — Designer 1D View

The **Designer 1D View** performs single-sequence 1D primer truncation analysis.
It evaluates a candidate DNA sequence across progressive length truncations down
to a specified minimum length, calculating self-dimerisation potential,
alignment quality, and binding statistics at each truncation step.

The view features a two-column layout with resizable panel dividers, allowing
user control over split panel widths and heights.

## Layout & Panel Resizing

- **Two-Column Split**:
  - **Left Column**: Contains the **1D Truncation Parameters** form (top) and
    the **Generated Primers** list (bottom).
  - **Right Column**: Contains the **Self-Dimer Quality by Primer Size (nt)**
    chart (top) and the **Self-Dimer Cards** panel (bottom).
- **Interactive Resizers**:
  - **Main Horizontal Divider**: Drag the vertical splitter bar (minimum 250 px
    left width) to resize the two side-by-side main panels.
  - **Left Vertical Divider**: Drag the horizontal splitter bar (minimum 110 px
    top height, default 240 px) to adjust the height of the top-left parameters
    form relative to the bottom-left primer list.
  - **Right Vertical Divider**: Drag the horizontal splitter bar (minimum 70 px
    top height, default 240 px) to adjust chart height relative to the
    bottom-right cards container. Resizing dynamically re-scales and rebuilds
    the quality bar chart.

## 1D Truncation Parameters (Top-Left Panel)

The **1D Truncation Parameters** form configures the input DNA sequence, length
constraints, and dimer filtering parameters. Truncation always proceeds in the
forward direction, removing bases from the 3' end and maintaining the 5'
terminus:

- **Candidate Primer Sequence**: Enter the nucleotide sequence to analyse. Raw
  sequence input is automatically cleaned to filter out non-nucleotide
  characters. Must contain at least one valid base.
- **Length (nt)**: Read-only counter next to the sequence field showing the
  cleaned sequence length in nucleotides. Updates as you type.
- **Min Length (nt)**: Minimum primer length. Required field (no default); leave
  the "e.g. 18" hint in mind when filling it in. The analysis generates
  candidate primers base-by-base from full sequence length down to this minimum
  value. Must be a positive integer greater than 0 and cannot exceed the input
  sequence length.
- **Max Quality**: Maximum quality score threshold. Integer only (no default).
  Leave blank for unconstrained quality filtering.
- **Max Overlap (bp)**: Maximum overlap length constraint. Non-negative integer
  (no default). Leave blank for unconstrained overlap filtering.
- **Check against template**: Checkbox to enable checking candidate primers
  against the template DNA sequence defined in the Input view.
- **Max Binding Sites**: Maximum number of template binding sites (replication
  origins) allowed for each candidate primer. Enabled when "Check against
  template" is checked. Leave blank ("Unconstrained if empty") to evaluate and
  display binding site counts without filtering out primers. When specified,
  must be a positive integer greater than 0.
- **Analyse Button**: Positioned in the bottom parameters row after Max Binding
  Sites. Triggers validation and runs the 1D primer truncation analysis. Can
  also be executed by pressing Enter inside any input field.
- **Save / Load / Clear All Parameters**:
  - **Save Button**: Saves the current form parameters (including sequence
    filter and max binding sites) to a YAML file using a file save dialog.
  - **Load Button**: Opens a file picker dialog to import parameters from a
    `.yaml` or `.yml` file. Automatically populates input fields, clears
    previous errors, and executes analysis.
  - **Clear All Button**: Clears all form parameters and the analysis results
    (primer list, chart, and cards).

## Generated Primers List (Bottom-Left Panel)

The bottom-left panel displays all candidate primer lengths generated during
truncation analysis in a vertical scrollable list. While an analysis is running,
the list is replaced by a determinate **progress bar** labelled *"Analysing
primer truncations…"* that fills as each truncation length is evaluated.

- **Progress & Abort**: During analysis, the **Analyse** button changes to an
  **Abort** button ("Stop analysis and keep results so far"). Clicking **Abort**
  cancels the running analysis and keeps the primers analysed so far, showing a
  notification: *"Analysis aborted — showing primers analysed so far."*
- **Primer Item Cards**: Each card represents a candidate primer step:
  - **Length Header**: Displays primer length in nucleotides (e.g. `20 nt`).
  - **Sequence Field**: Read-only monospace text field, left-aligned.
  - **Metrics Badges**: Displays highlighted badges for **Quality**
    (`Quality: {score}`, rounded to the nearest integer), **Overlap**
    (`Overlap: {length} bp`), and optionally **Binding Sites**
    (`Sites: {count}`) when "Check against template" is enabled.
  - **PCR Button**: A dedicated button on each card to run the PCR simulation
    using the template DNA sequence and this candidate primer directly.
- **Interactive Selection**: Clicking anywhere on a primer item card selects
  that primer step and opens or brings to top its detailed self-dimer card in
  the right-hand panel.

## Self-Dimer Quality Bar Chart (Top-Right Panel)

The top-right panel renders an interactive vertical quality bar chart inside a
horizontally scrollable container:

- **Empty State**: Displays *"No analysis results yet. Enter sequence and click
  Analyse."* prior to running an analysis.
- **Bar Display**: Each generated primer truncation step is represented by a
  vertical bar:
  - **Score**: Self-dimer quality score value printed above each bar (rounded to
    the nearest integer).
  - **Bar Height**: Scaled relative to the maximum quality score among all
    generated primers and current panel container height.
  - **X-Axis Label**: Displays primer length in nucleotides (e.g. `20 nt`) below
    each bar.
  - **Tooltip**: Hovering over a bar displays details including truncation step
    number, primer length, and quality score.
- **Interactive Selection**: Clicking a chart bar selects the corresponding
  primer and displays or raises its self-dimer card in the bottom-right panel.

## Self-Dimer Cards Panel (Bottom-Right Panel)

The bottom-right panel displays detailed self-dimer alignment cards in a
vertical scrollable list:

- **Header Controls**:
  - **Panel Header**: Displays title **Self-Dimer Cards**.
  - **Clear Cards Button**: Appears when one or more cards are present. Clicking
    **Clear Cards** dismisses all open cards simultaneously.
- **Card Selection & Positioning**:
  - Selecting a primer from the generated primers list or quality bar chart
    creates a dismissible card positioned at the top of the cards list.
  - If a card for that specific primer step already exists in the list, it is
    automatically raised to the top of the stack.
- **Card Contents**:
  - **Card Header**: Displays title `Self-dimer ({length} nt)` alongside a
    close/dismiss button and a **Run PCR** button to simulate PCR using the
    template with this primer.
  - **Thermodynamic & Base Composition Badges**:
    - **Quality Score**: `Quality: {score}`
    - **Overlap**: `Overlap: {overlap} bp`
    - **Binding Sites**: `Sites: {count}` (if "Check against template" is
      enabled)
    - **Melting Temperature ($T_m$)**: `Tm: {value}°C` calculated using
      configured thermodynamic settings (shown as `Tm: N/A` when it cannot be
      calculated).
    - **% AT Content**: `% AT: {percentage}%`.
  - **Antiparallel Alignment Diagram**: Rendered in a monospace font displaying
    a 5-line structural alignment diagram:
    1. Top primer name
    2. Top sequence ($5' \\to 3'$)
    3. Bond interaction line (`|` for strong matches $\\ge 10.0$, `:` for weak
       matches $0.0 \\le \\text{score} < 10.0$, space for mismatches)
    4. Bottom sequence ($3' \\to 5'$)
    5. Bottom primer name
- **Individual Dismissal**: Clicking the close button on an individual card
  removes it from the panel.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
