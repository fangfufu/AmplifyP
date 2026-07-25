# AmplifyP GUI User Manual

This user manual provides detailed guidance for using the graphical user
interface (GUI) of **AmplifyP**.

AmplifyP is a cross-platform application for simulating polymerase chain
reactions (PCR), predicting amplicons, evaluating primer dimers, and designing
candidate primers. It is a modern Python rewrite of William Engels's
**Amplify4** Mac program, supporting Linux, Windows, macOS, and web browsers.

<table>
  <tr>
    <td width="33%" align="center">
      <b>Template & Primer Management</b>
    </td>
    <td width="33%" align="center">
      <b>Primer Binding Site & Amplicon</b>
    </td>
    <td width="33%" align="center">
      <b>Primer Dimer Analysis</b>
    </td>
  </tr>
  <tr>
    <td width="33%" align="center">
      <img src="images/input_view.png" alt="Input View" width="100%">
    </td>
    <td width="33%" align="center">
      <img src="images/pcr_view.png" alt="PCR Results View" width="100%">
    </td>
    <td width="33%" align="center">
      <img src="images/primer_dimer_view.png" alt="Primer Dimer Analysis" width="100%">
    </td>
  </tr>
</table>

## Documentation by View

The user manual is divided into dedicated guides for each view in the interface:

1. **[Input View](input_view.md)**: Manage DNA template sequences (supporting
   linear and circular topologies), sequence validation and cleaning, and primer
   lists (adding, editing, toggling active status, computing reverse
   complements, and importing/exporting FASTA or CSV files).
2. **[PCR View](pcr_view.md)**: View interactive sequence maps displaying primer
   binding sites with directionality and match strength indicators, predicted
   amplicon fragment bars with quality scores ($Q$), circular wraparound
   amplicons, detailed breakdown cards, and export options.
3. **[Primer Dimers View](primer_dimer_view.md)**: Analyse potential
   self-dimerization and cross-dimerization risks across active primers, inspect
   base-pairing alignment structures, evaluate free energy stability ($\\Delta
   G$), and filter or sort candidate dimer pairs.
4. **[Designer 1D View](designer_1d_view.md)**: Perform single-sequence 1D
   primer truncation analysis across a template sequence to evaluate binding
   quality scores at varying lengths, identify optimal primer lengths, and
   detect internal priming sites.
5. **[Designer 2D View](designer_2d_view.md)**: Perform pairwise 2D primer
   truncation matrix analysis for forward and reverse primer pairs across length
   combinations to maximise amplicon yield while minimising dimer risks.
6. **[Settings & Preferences](settings_view.md)**: Customise algorithm cutoffs,
   thermodynamic parameters ($T_m$ calculation methods, salt concentrations,
   annealing temperatures), application themes, and operational preferences.

## Command Line Options

AmplifyP can be launched from the command line with options to configure window
dimensions, load pre-saved state files, run in web browser mode, or
automatically export view screenshots.

When running from source or using the command line launcher, pass arguments as
follows:

```bash
python src/main.py [OPTIONS]
```

or with the compiled binary executable:

```bash
amplifyp [OPTIONS]
```

### Available Arguments

| Option                     | Short | Description                                                                        |
| :------------------------- | :---- | :--------------------------------------------------------------------------------- |
| `--state <PATH>`           | `-f`  | Path to a YAML state file to load automatically on startup.                        |
| `--auto-close`             |       | Automatically quit the application after rendering completes (requires `--state`). |
| `--screenshots`            | `-s`  | Save PNG screenshots of views on load (requires `--state`).                        |
| `--screenshots-dir <DIR>`  |       | Target directory path where saved PNG screenshots are stored.                      |
| `--window-width <PIXELS>`  |       | Set initial application window width in pixels.                                    |
| `--window-height <PIXELS>` |       | Set initial application window height in pixels.                                   |
| `--web`                    |       | Launch application in web browser mode rather than desktop window mode.            |
| `--help`                   | `-h`  | Display command help message and exit.                                             |

### Usage Examples

- **Set application window size**:

  ```bash
  python src/main.py --window-width 1280 --window-height 800
  ```

- **Launch in web browser mode**:

  ```bash
  python src/main.py --web
  ```

- **Load a pre-saved YAML state file**:

  ```bash
  python src/main.py --state saved_session.yaml
  ```

- **Automated screenshot export and auto-close**:

  ```bash
  python src/main.py --state saved_session.yaml --screenshots --screenshots-dir ./screenshots --auto-close
  ```
