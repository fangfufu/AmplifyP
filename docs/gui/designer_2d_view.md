# AmplifyP GUI — Designer 2D View

The **Designer 2D** tab performs 2D primer truncation analysis. It takes
separate forward and reverse DNA sequences and truncates each from the 3' end
down to a minimum length, evaluating all combinations of forward-reverse primer
pairs.

## Truncation Parameters (Top-Left Panel)

- **Forward Candidate Primer Sequence**: The forward primer target sequence.
- **Fwd Min Len (bp)**: Minimum length for forward primer truncation (default
  18).
- **Reverse Candidate Primer Sequence**: The reverse primer target sequence.
- **Rev Min Len (bp)**: Minimum length for reverse primer truncation (default
  18).
- **Quality Filter**: Optional minimum quality threshold for primer pairs. Pairs
  below this score are excluded. Leave empty for no filtering.
- **Overlap Filter (bp)**: Optional maximum overlap constraint for primer pairs.
  Leave empty for no filtering.
- **Metric**: Choose **Max** (use the maximum dimer quality across all steps) or
  **Mean** (use the average) to rank and filter results.

## Results Grid (Bottom-Left Panel)

A colour-coded grid showing all forward-reverse primer pair combinations. Each
cell represents a pair at a specific truncation step, with quality scores mapped
to colours for quick visual identification. Click on a cell to view the detailed
pair analysis in the right panel.

## 2D Primer Pair Detail Cards (Right Panel)

Selecting a grid cell opens a detail card showing the forward and reverse primer
sequences, their alignment, quality scores, and binding statistics. Cards can be
dismissed individually or cleared with the **Clear Cards** button.

## Save / Load

Parameters (both sequences, minimum lengths, filters, metric) can be saved to a
YAML file and reloaded to resume work.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
