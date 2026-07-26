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
  - **Right Column**: Contains the **Self-Dimer Quality by Primer Size (bp)**
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

The **1D Truncation Parameters** form configures the input DNA sequence,
truncation direction, length constraints, and dimer filtering parameters:

- **Candidate Primer Sequence**: Enter the nucleotide sequence to analyse. Raw
  sequence input is automatically cleaned to filter out non-nucleotide
  characters. Must contain at least one valid base.
- **Min Length (bp)**: Minimum primer length integer (default: `18`). The
  analysis generates candidate primers base-by-base from full sequence length
  down to this minimum value. Must be a positive integer greater than 0 and
  cannot exceed the input sequence length.
- **Primer Direction**: Dropdown selector for truncation orientation:
  - **Forward**: Truncates sequence from the 3' end (maintains the 5' terminus).
  - **Reverse**: Truncates sequence from the 5' end (maintains the 3' terminus).
- **Max Quality**: Maximum quality score threshold. Defaults to the configured
  primer dimer quality threshold (e.g. `60.0`). Leave blank for unconstrained
  quality filtering.
- **Max Overlap (bp)**: Maximum overlap length constraint. Defaults to the
  minimum dimer overlap setting (e.g. `3`). Leave blank for unconstrained
  overlap filtering.
- **Analyse Button**: Triggers validation and runs the 1D primer truncation
  analysis. Can also be executed by pressing Enter inside any input field.
- **Save / Load Parameters**:
  - **Save Button**: Saves the current form parameters to a YAML file using a
    file save dialog.
  - **Load Button**: Opens a file picker dialog to import parameters from a
    `.yaml` or `.yml` file. Automatically populates input fields, clears
    previous errors, and executes analysis.

## Generated Primers List (Bottom-Left Panel)

The bottom-left panel displays all candidate primer lengths generated during
truncation analysis in a vertical scrollable list:

- **Primer Item Cards**: Each card represents a candidate primer step:
  - **Length Header**: Displays primer length in base pairs (e.g. `20 bp`).
  - **Sequence Field**: Read-only monospace text field. Text is left-aligned for
    Forward direction primers and right-aligned for Reverse direction primers.
  - **Metrics Badges**: Displays highlighted badges for **Quality**
    (`Quality: {score}`) and **Overlap** (`Overlap: {length} bp`).
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
  - **X-Axis Label**: Displays primer length in base pairs (e.g. `20 bp`) below
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
  - **Card Header**: Displays title `Self-dimer ({length} bp)` alongside a
    close/dismiss button.
  - **Thermodynamic & Base Composition Badges**:
    - **Quality Score**: `Quality: {score}`
    - **Overlap**: `Overlap: {overlap} bp`
    - **Melting Temperature ($T_m$)**: `Tm: {value}°C` calculated using
      configured thermodynamic settings.
    - **AT Pairs**: `AT Pairs: {count}`
    - **GC Pairs**: `GC Pairs: {count}`
    - **% AT Content**: `% AT: {percentage}%`
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
