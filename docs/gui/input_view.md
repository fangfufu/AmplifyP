# AmplifyP GUI — Input View

The **Input View** (also referred to as the Substrate View) is the main
configuration interface in AmplifyP. It allows you to specify the target DNA
template sequence, configure DNA topology, and manage the list of primers used
for PCR simulation and primer dimer analysis.

The window is split horizontally into a **Template Sequence Panel** on the left
and a **Primer List Panel** on the right, separated by an interactive, draggable
vertical divider.

![Input View](images/input_view.png)

## Template Sequence Panel (Left Panel)

The left panel provides controls for entering, formatting, and manipulating the
target DNA template sequence.

- **Sequence Input & Auto-Cleaning**:
  - Multiline monospace text field for entering or editing sequence data.
  - Automatic filtering: Non-biological characters (spaces, numbers, line
    breaks, punctuation) are automatically cleaned out. Valid base characters
    include standard nucleotides (`A`, `T`, `C`, `G`, `U`) and IUPAC
    degenerate/ambiguous bases (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`,
    `V`, `N`).
- **Dynamic Line Number Gutter**:
  - The left-hand gutter displays the starting base index for each row.
  - Line numbers update dynamically as sequence text or wrapping length changes.
  - Gutter width scales automatically based on the number of digits in the
    template sequence length.
- **Bases Per Line & Wrap Control**:
  - **Wrap Selector**: A popup menu in the status bar allows choosing the number
    of bases displayed per line (preset options from 10 to 100 in steps of 10,
    or **Auto**).
  - **Auto Wrapping**: When set to **Auto**, the wrap length is calculated
    dynamically based on available panel width.
  - **Horizontal Scrollbar**: Automatic text wrapping at window borders is
    disabled in favour of horizontal scrolling to maintain strict visual
    alignment.
- **Topology Toggle**:
  - **Circular Checkbox**: Toggle the template topology between **Linear** and
    **Circular**. Circular topology is essential for simulating plasmid PCR,
    allowing predicted amplicons to wrap around origin base 1.
- **Sequence Casing Controls**:
  - **Upper**: Converts currently selected bases in the template sequence to
    upper case.
  - **Lower**: Converts currently selected bases in the template sequence to
    lower case.
  - *(Note: If no text is selected when clicking upper or lower case buttons, a
    notification prompt appears).*
- **Template File Operations**:
  - **Load**: Opens a file picker to import a template sequence from a plain
    text file (`.txt`).
  - **Save**: Saves the current template sequence to a plain text file (`.txt`).
  - **Clear**: Clears the template sequence text field.
- **Status Bar**:
  - Located at the bottom of the template box, displaying sequence metadata:
    - **Total Bases**: Shows total cleaned sequence length (e.g.,
      `Total Bases: 1500`) when the text box is not focused.
    - **Insertion Point**: Displays current cursor position (e.g.,
      `Insertion Point After Base: 450`) when focused without selection.
    - **Selected Range**: Displays base index range (e.g.,
      `Selected Bases: 120 - 250`) when text is highlighted.
    - **Bases per Line Menu**: Shows current wrap setting with popup menu
      trigger.

## Primer List Panel (Right Panel)

The right panel manages the list of primers used for PCR simulation and dimer
checking.

### Interactive Primer Table

The table comprises the following columns:

- **Drag Handle**: Visual drag icon (`⋮⋮`) for live drag-and-drop row
  reordering.
- **Active Checkbox**: Checkbox to include or exclude individual primers from
  simulation.
  - **Header Checkbox**: Tri-state checkbox in the table header to toggle all
    primers active or inactive at once.
- **Name Column**: Multiline editable text field for primer identification. The
  column width can be resized by dragging the vertical divider in the header.
- **Sequence Column**: Monospace text field for primer sequence entry ($5' \\to
  3'$). Sequences are automatically cleaned and validated.
- **Melting Temperature ($T_m$) Column** *(Optional)*:
  - Enabled when `show_primer_temperature` is set in Settings.
  - Displays calculated $T_m$ (°C) in real-time for each sequence.
  - Can be colour-coded dynamically based on the active $T_m$ colour scheme.

### Row Selection & Reordering

- **Single-Click Selection**: Click any row to highlight it and focus its name
  field.
- **Double-Click Range Selection**: Double-click a row to set an anchor
  (`click_a`), then double-click a second row (`click_b`) to highlight the range
  of primers between them.
- **Live Drag-and-Drop**: Click and drag the drag handle (`⋮⋮`) to move single
  rows or contiguous highlighted blocks up or down in real time.
- **Header Reordering Controls**:
  - **Add Primer (+)**: Inserts a new empty primer row below the currently
    selected or focused row.
  - **Delete Primer (-)**: Deletes all currently highlighted primer rows.
  - **Move Up (↑)** / **Move Down (↓)**: Moves selected primer row(s) up or down
    by one position.
- **Auto-Append Row**: Entering a valid sequence into the final row
  automatically appends a new empty row to the list.

### Toolbar Actions

Located at the top of the Primer List Panel:

- **Load**: Imports primers from CSV (`.csv`), TSV (`.tsv`), or text (`.txt`)
  files.
- **Save**: Exports all filled primers to a TSV file (`.tsv`).
- **Copy**: Copies highlighted (or currently focused) primers to the system
  clipboard in TSV format (`Name\tSequence`).
- **Paste**: Pastes primers from the system clipboard starting at the selected
  row position or end of the list. Supports TSV, CSV, or raw sequence lines.
- **Rev Comp**: Calculates and replaces the sequence of all highlighted primers
  with their reverse complement sequence.
- **Delete**: Deletes highlighted primers.
- **Clear**: Resets the primer list to a single empty row.

### Primer Information & Analysis Panel

Clicking or focusing a primer displays a detailed analysis card (positionable at
the top or bottom of the list via Settings):

- **Header & Sequence**: Displays primer name, length in base pairs, and full
  sequence ($5' \\to 3'$).
- **Melting Temperature**: Displays exact $T_m$ (°C).
- **Base Composition**: Displays count of AT pairs, GC pairs, and overall AT
  percentage.
- **Degeneracy & Redundancy**: Displays redundant IUPAC base count and
  calculated redundancy fold.
- **Self-Dimer Analysis**: Performs real-time self-dimer prediction. If
  self-dimer binding exceeds threshold settings, an interactive visual dimer
  alignment card is rendered within the info panel.

### Validation & Error Warnings

- Real-time validation flags illegal base characters, missing required fields
  (for active primers), duplicate primer names, and duplicate primer sequences.
- Cells with errors display explicit red error messages.
- **Error Banner**: If active primers contain invalid format or duplicate
  entries, a prominent red warning banner is displayed at the bottom of the
  panel:
  > *"PCR and Primer Dimer views are disabled because one or more selected
  > primers are invalid, or have duplicated names/sequences."*

## Layout & Panel Resizing

- **Central Divider**: Drag the vertical divider bar between the left and right
  panels to adjust relative panel widths (minimum 200 px width per panel).
- **Column Resizing**: Drag the vertical divider between the Name and Sequence
  header columns to resize column proportions. Column ratios scale
  proportionally on window resize.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
