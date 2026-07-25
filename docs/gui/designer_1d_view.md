# AmplifyP GUI — Designer 1D View

The **Designer 1D** tab performs 1D primer truncation analysis. It takes a
single DNA sequence and progressively truncates it from the 3' end down to a
minimum length, evaluating self-dimer potential at each truncation step.

## Truncation Parameters (Top-Left Panel)

- **Candidate Primer Sequence**: Enter the target DNA sequence to analyse. Valid
  nucleotide characters (`A`, `T`, `C`, `G`) are accepted; invalid characters
  are filtered out.
- **Min Length (bp)**: The minimum primer length for truncation (default 18).
  The analysis generates one candidate primer per base from the full length down
  to this minimum.
- **Primer Direction**: Choose **Forward** (truncation from the 3' end) or
  **Reverse** (truncation from the 5' end).
- **Max Quality**: Optional quality threshold. Self-dimers with quality scores
  above this value are excluded from results. Leave empty for no filtering.
- **Max Overlap (bp)**: Optional maximum overlap constraint. Self-dimers with
  overlap greater than this value are excluded. Leave empty for no filtering.

## Generated Primers (Bottom-Left Panel)

Each truncation step produces a candidate primer. The list shows all generated
primers with their sequence, length, and quality metrics. Click on a primer to
display its self-dimer analysis card in the right panel.

## Self-Dimer Quality Chart (Top-Right Panel)

A bar chart displaying the self-dimer quality score for each primer by size
(bp). The chart is resizable — drag the horizontal divider to adjust height.
Click on a bar to select the corresponding primer and view its self-dimer card.

## Self-Dimer Cards (Bottom-Right Panel)

Clicking a primer or chart bar opens a detailed self-dimer card showing the
alignment, match quality, and binding statistics. Cards can be dismissed
individually or cleared all at once with the **Clear Cards** button.

## Save / Load

Parameters (DNA sequence, min length, direction, filters) can be saved to a YAML
file and reloaded to resume work.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
