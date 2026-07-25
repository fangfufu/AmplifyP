# AmplifyP GUI — PCR View

Click the **PCR** tab to run the simulation. The screen is split into a visual
lane/map at the top and a detailed textual analysis panel at the bottom.

![PCR View](images/pcr_view.png)

## The PCR Map (Top Panel)

- **Primer Match Arrows**: Red (forward) and blue (reverse) arrows show where
  primers bind to the target DNA. The size of the arrow indicates the strength
  of the match. Hovering over an arrow displays basic statistics in the tooltip.
- **Predicted Amplicons**: Horizontal lines below the arrows represent the
  predicted PCR products. The label displays the length of the fragment in base
  pairs (including the primers).

## Textual Analysis (Bottom Panel)

Click on any primer match arrow or amplicon line in the map to see comprehensive
details in the output log:

- **Primer Match Details**: Alignments show exactly where the primer binds to
  the template. Base pairings are marked with a vertical bar (`|`) for exact
  matches and a colon (`:`) for ambiguous matches. It also shows:
  - **Primability**: A percentage score weighting matches heavier near the
    critical 3' end.
  - **Stability**: A percentage score reflecting the overall binding stability
    across the primer length.
  - **Quality**: A normalised score between 0.0 and 1.0 based on the sum of
    primability and stability.
- **Amplicon Details**: Selecting an amplicon displays its complete sequence
  (with primer regions highlighted), primer names, and the quality score ($Q$).
  - **Amplicon Q**: Short fragments amplify better. The quality score is defined
    as: $$Q = \\frac{\\text{Length of fragment in bp}}{(\\text{Left Match
    Quality} \\times \\text{Right Match Quality})^2}$$ For perfect primer
    matches, $Q$ is equal to the amplicon length. Poorer primer matches result
    in a larger $Q$ value (lower amplification efficiency).

______________________________________________________________________

[Return to GUI Manual Index](README.md)
