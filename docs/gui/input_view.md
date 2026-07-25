# AmplifyP GUI — Input View

The **Input View** (also referred to as the Substrate View) is the main window
where you configure your experiment. It is divided into a target sequence panel
on the left and a primer list panel on the right.

![Input View](images/input_view.png)

## Target Sequence Panel (Left)

- **Sequence Input**: Paste or type your target DNA sequence. Valid base
  characters are standard nucleotides (`A`, `T`, `C`, `G`) and
  wildcard/ambiguous bases (`N`). Non-nucleotide characters (such as spaces or
  numbers) are automatically filtered out.
- **Topology Toggle**: Toggle the DNA topology between **Linear** and
  **Circular** (useful for simulating plasmid PCR).
- **Line Numbers**: The left column displays the base index at the start of each
  row, which adjusts dynamically if you resize the window.
- **Base Style Actions**: Choose to convert sequences to uppercase or lowercase,
  or invert (flip/reverse-complement) the entire sequence.

## Primer List Panel (Right)

- **Active Primers**: Use the checkbox in the **Active** column to check primers
  you wish to include in the PCR simulation. Only checked primers will be used.
- **Manage Primers**:
  - Click the **+** button to add a new primer row.
  - Click the **-** button to delete selected primers.
  - Double-click or click inside cells to edit the **Sequence**, **Name**, or
    **Notes**.
- **Actions**:
  - **Flip Sequence (🔄)**: Reverses and complements the sequence of the selected
    primer.
  - **Save / Load States**: Save all template sequences, primers, and settings
    to a local YAML file (`.yaml`), or reload a saved configuration to resume
    your work.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
