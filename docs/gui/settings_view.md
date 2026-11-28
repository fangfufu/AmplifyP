# AmplifyP GUI — Settings & Preferences View

The **Settings & Preferences View** provides comprehensive options for
customising algorithm parameters, visual appearance, primer list behaviour,
melting temperature ($T_m$) calculations, primer dimer threshold scoring,
designer view display options, logging diagnostics, and local setting
persistence.

Settings are organised into seven collapsible expansion tiles:

1. **General**
2. **Input**
3. **PCR**
4. **Primer Dimers**
5. **Designer 1D**
6. **Designer 2D**
7. **Diagnostics**

A **Reset to Default** button at the bottom opens a confirmation dialogue to
restore all configuration values to their factory defaults.

## 1. General

Configures application appearance, template auto-loading on startup, software
update checks, settings backup file operations, and factory reset.

- **Appearance**:
  - **Font Family**: Monospace font dropdown selector (`Roboto Mono` [default],
    `Courier New`, `Consolas`, `monospace`). Applied across sequence displays
    and alignment cards.
  - **Colour Scheme**: Theme dropdown selector (`Light`, `Dark`,
    `Light (Colour Deficient Friendly)`, `Dark (Colour Deficient Friendly)`,
    `System` [default], `System (Colour Deficient Friendly)`). Syncs dark mode
    and colour deficient mode settings.
  - **5'/3' Sequence Separator**: Sequence boundary spacing selector
    (`Space (' ')`, `Dash ('-')` [default]). Configures whether sequence ends in
    dimer alignment cards and PCR context cards format with whitespace
    (`5' ... 3'`) or hyphens (`5'-...-3'`).
- **Autosave**:
  - **Automatically reload last template and primers on startup**: Bordered
    checkbox (default: checked). When enabled, automatically saves the active
    session state and restores the last opened sequence and primer list upon
    launching.
- **Updates**:
  - **Version Checking Frequency**: Frequency for checking new release tags on
    GitHub (`At Startup`, `Once per Day` [default], `Once per Week`,
    `Once per Month`, `Disabled`).
  - **Check for Updates Button**: Manually triggers an asynchronous version
    check against the current application version. A status line next to the
    button shows *"Checking for updates..."* while the check runs, then the
    outcome in green (*"AmplifyP is up to date."* or *"New version {tag} is
    available!"*) or in red (*"Could not check for updates. Please check your
    network connection."*) when the check fails. When a newer release tag is
    found, the version text in the application header also changes to
    *"{version} (Update {tag} available!)"* and becomes a clickable link to the
    GitHub releases page.
- **Backup Settings**:
  - **Save Settings**: Opens a file dialogue to export all current settings into
    a YAML configuration file.
  - **Load Settings**: Opens a file dialogue to import settings from a YAML
    configuration file and immediately apply theme and parameter changes across
    the UI.
- **Reset Settings**:
  - **Reset Settings**: Opens a confirmation dialogue and, upon confirmation,
    resets all preferences across all tiles to their default values, saves the
    updated state to local storage, reconfigures logging, and refreshes the UI.

## 2. Input

Controls primer list layout, temperature displays, automatic activation
behaviour, and melting temperature ($T_m$) calculations within the
[Input View](input_view.md).

- **Primer List**:
  - **Primer Info Panel Position**: Dropdown selecting whether the primer
    details card is docked at the `Top` or `Bottom` [default] of the primer list
    panel.
  - **Fixed height primer info box**: Bordered checkbox (default: unchecked).
    Enforces a fixed height for the primer info panel to prevent layout shifting
    during primer selection.
  - **Auto-activate new valid primer**: Bordered checkbox (default: unchecked).
    Automatically checks the active checkbox for newly added primers once valid
    sequences are entered.
  - **Show primer temperature column**: Bordered checkbox (default: unchecked).
    Displays an extra $T_m$ column in the primer table showing live melting
    temperature calculations.
  - **Primer temperature colour scheme**: Dropdown palette selector (`None`
    [default], `Cool-Warm`, `Traffic Light`) for background colour-coding of
    $T_m$ values in the primer list.
- **Primer Melting Temperature (Tm)**:
  - **Tm Calculation Method**: Dropdown selector for the algorithm:
    - **SantaLucia 1998 / Owczarzy 2008 (Default)**: Nearest-neighbour
      thermodynamics with salt corrections for monovalent and divalent cations.
    - **Lander / Amplify 4**: Empirical formula from the original Amplify4
      program.
    - *(Note: If SantaLucia thermodynamics encounters unsupported non-standard
      or degenerate bases, it automatically falls back to Lander / Amplify 4).*
  - **Reaction Concentrations**:
    - **DNA Conc (nM)**: Total primer strand concentration in nM (default:
      `50.0`).
    - **DNA Pol Conc**: Polymerase concentration (default: `0.0`).
    - **Monovalent Salt Conc (mM)**: Concentration of monovalent cations
      ($Na^+$, $K^+$, $Tris^+$) in mM (default: `50.0`).
    - **Divalent Salt Conc (mM)**: Concentration of divalent cations ($Mg^{2+}$)
      in mM (default: `1.5`).
    - **dNTP Conc (mM)**: Concentration of dNTPs in mM (default: `0.0`).

## 3. PCR

Configures overview map layout, pairwise nucleotide scoring weights, lengthwise
match/run weights, primability/stability thresholds, and compatibility modes for
PCR replication calculations.

- **Overview Map**:
  - **Horizontal divider behaviour**: Dropdown controlling how the PCR overview
    map height is sized (`Auto-resizing, max 1/4 of the viewport`,
    `Auto-resizing, max 1/3 of the viewport`,
    `Auto-resizing, max 1/2 of the viewport`,
    `Auto-resizing, max 2/3 of the viewport`,
    `Auto-resizing, max 3/4 of the viewport` [default],
    `Save proportionality across session`).
  - **Primer label orientation**: Dropdown controlling primer name orientation
    on the PCR overview map (`Horizontal` [default], `45°`, `Vertical`).
- **Replication Origin**:
  - **Base Pair Weights**: A $15 \\times 5$ interactive table of pairwise
    weights for primer nucleotides (`G`, `A`, `T`, `C`, IUPAC degenerate codes
    `M`, `R`, `W`, `S`, `Y`, `K`, `V`, `H`, `D`, `B`, `N`) vs template
    nucleotides (`G`, `A`, `T`, `C`, `N`), with a **Reset** button to restore
    defaults.
  - **Match Weights**: An interactive table of positional match weights indexed
    from the primer 3' end (`3' Position` 1–15 plus `Default`), with add (`+`),
    remove (`-`), and **Reset** controls.
  - **Run Weights**: An interactive table of consecutive match run weights
    (`Run Length` 1–5 plus `Default`), with add (`+`), remove (`-`), and
    **Reset** controls.
  - **Primability Cutoff**: Minimum fractional primability score required for a
    match site to be recognised (default: `0.8`).
  - **Stability Cutoff**: Minimum fractional 3' stability score required for a
    match site to be recognised (default: `0.4`).
  - **Amplify 4 Compatibility Mode**: Bordered checkbox (default: unchecked).
    Enforces legacy Amplify4 scoring logic and weighting constraints.
  - **Improved Replication Origin Visualisation**: Bordered checkbox (default:
    checked). Enables enhanced visual match diagrams and complementary strand
    display in PCR view cards.
- **Amplicon**:
  - **Amplicon vertical ranking**: Dropdown controlling the vertical ordering of
    amplicon fragment bars in the PCR View diagram (`Position, then length`
    [default], `Quality score`, `Position, then quality`).
  - **Relative Q score**: Bordered checkbox (default: checked). When enabled,
    scales each amplicon's $Q$ score by dividing by the smallest positive finite
    $Q$ score across all predicted amplicons and multiplying by 100.

## 4. Primer Dimers

Configures scoring and thresholds for evaluating primer self-dimers and
cross-dimers in the Primer Dimers View and info panels.

- **Primer Dimer Weights**:
  - A $15 \\times 15$ interactive table defining pairwise scoring weights for
    all primer base vs primer base interactions (including IUPAC ambiguous
    bases), with a **Reset** button to restore default weights.
  - Default values are derived from Amplify4 thermodynamic binding weights.
- **Parameters**:
  - **Min Overlap**: Minimum overlap length (in base pairs) required for a
    primer dimer alignment to be evaluated (default: `3`).
  - **Threshold**: Minimum overall quality score threshold for a dimer
    interaction to be reported in the [Primer Dimers View](primer_dimer_view.md)
    (default: `60.0`).

## 5. Designer 1D

Configures display options for the self-dimer cards in the
[Designer 1D View](designer_1d_view.md).

- **Show melting temperature (Tm) on cards**: Bordered checkbox (default:
  checked). Displays a $T_m$ badge on each self-dimer card header.
- **Show % AT on cards**: Bordered checkbox (default: unchecked). Displays a %
  AT content badge on each self-dimer card header.

## 6. Designer 2D

Configures display options for pair-wise 2D primer truncation analysis in
Designer 2D View.

- **2D Results Grid colour map scheme**: Dropdown palette selector (`None`,
  `Cool-Warm`, `Traffic Light`, `Blue-Orange` [default], `Greyscale`). Controls
  the heat map colour scheme for displaying amplicon scores across 2D primer
  length combination grids.
- **Show Reverse-Forward cross-dimer**: Bordered checkbox (default: unchecked).
  When enabled, Designer 2D also evaluates and displays the Reverse-Forward
  cross-dimer alignment (reverse primer 3' end against the forward primer) in
  addition to the three default dimer alignments.

## 7. Diagnostics

Configures application logging, log file output, log level thresholds, and log
file rotation parameters.

- **Console Output**: Bordered checkbox (default: checked). Toggles logging
  messages to standard console output.
- **AmplifyP Log Level**: Dropdown setting the log level threshold for
  application logic (`DEBUG`, `INFO` [default], `WARNING`, `ERROR`, `CRITICAL`).
- **Flet Log Level**: Dropdown setting log level threshold for the underlying UI
  framework (`INFO` [default], `WARNING`, `ERROR`, `CRITICAL`).
- **File Logging**: Bordered checkbox (default: checked on desktop platforms,
  hidden on web, where *"File logging is not available in web mode."* is shown
  instead). Enables writing log messages to a disk file. The log rotation
  controls below are hidden while file logging is disabled.
- **Log File Path**: Dropdown selector with options:
  - `(Default)`: Writes to OS-specific application log directory (e.g.
    `~/.config/amplifyp/app.log` or `%APPDATA%\AmplifyP\app.log`).
  - `Select folder`: Launches a directory picker dialogue to select a custom
    destination folder for `app.log`.
- **Log Rotation**: Bordered checkbox (default: checked on desktop). Enables
  rotating log files when size limits are reached.
- **Max Log Size (MB)**: Dropdown selector for max file size before rotation
  (`1`, `5` [default], `10`, `20`, `50` MB).

## Reset & Local Persistence

- **Reset to Default Button**: Located at the bottom of the Settings View.
  Clicking **Reset to Default** opens a confirmation dialogue and, upon
  confirmation, resets all configuration settings across all tiles to factory
  defaults, saves the updated state to local storage, reconfigures logging, and
  refreshes the application UI.
- **Storage Locations**:
  - **Desktop Platforms (Linux, macOS, Windows)**: Settings are persisted in
    YAML format at the OS-specific config path (`settings.yaml` under
    `$XDG_CONFIG_HOME/amplifyp/`, `~/Library/Application Support/AmplifyP/`, or
    `%APPDATA%\AmplifyP\`).
  - **Web Platform**: Settings are saved directly to the browser's local
    storage.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
