# AmplifyP GUI — Settings & Preferences View

The **Settings & Preferences View** provides comprehensive options for
customising algorithm parameters, visual appearance, primer list behavior,
melting temperature ($T_m$) calculations, primer dimer threshold scoring, 2D
designer colour mapping, logging diagnostics, and local setting persistence.

Settings are organised into seven collapsible expansion tiles:

1. **General Settings**
2. **Primer List Settings**
3. **Primer Melting Temperature (Tm) Settings**
4. **PCR Settings**
5. **Primer Dimer Settings**
6. **Designer 2D Settings**
7. **Diagnostics Settings**

A **Reset to Default** button at the bottom restores all configuration values to
their factory defaults.

## 1. General Settings

Configures application appearance, template auto-loading on startup, software
update checks, and settings backup file operations.

- **Appearance**:
  - **Font Family**: Monospace font dropdown selector (`Roboto Mono` [default],
    `Courier New`, `Consolas`, `monospace`). Applied across sequence displays
    and alignment cards.
  - **Colour Scheme**: Theme dropdown selector (`Light`, `Dark`,
    `Light (Colour Deficient Friendly)`, `Dark (Colour Deficient Friendly)`,
    `System` [default], `System (Colour Deficient Friendly)`). Syncs dark mode
    and colour deficient mode settings.
- **Autosave / Reload**:
  - **Automatically reload last template and primers on startup**: Bordered
    checkbox (default: checked). When enabled, automatically saves the active
    session state and restores the last opened sequence and primer list upon
    launching.
- **Updates**:
  - **Version Checking Frequency**: Frequency for checking new release tags on
    GitHub (`At Startup`, `Once per Day`, `Once per Week`, `Once per Month`
    [default], `Disabled`).
  - **Check for Updates Button**: Manually triggers an asynchronous version
    check against the current application version. Reports whether AmplifyP is
    up to date or displays a notification when a new release tag is available.
- **Backup Settings**:
  - **Save Settings**: Opens a file dialog to export all current settings into a
    YAML configuration file.
  - **Load Settings**: Opens a file dialog to import settings from a YAML
    configuration file and immediately apply theme and parameter changes across
    the UI.

## 2. Primer List Settings

Controls table layout, duplicate warnings, temperature displays, and automatic
activation behavior within the [Input View](input_view.md) primer list.

- **Primer Info Panel Position**: Dropdown selecting whether the primer details
  card is docked at the `Top` or `Bottom` [default] of the primer list panel.
- **Fixed height primer info box**: Bordered checkbox (default: unchecked).
  Enforces a fixed height for the primer info panel to prevent layout shifting
  during primer selection.
- **Ignore inactive primers when checking for duplicate names**: Bordered
  checkbox (default: checked). Excludes inactive/deselected primers from
  triggering duplicate name warnings.
- **Ignore inactive primers when checking for duplicate sequences**: Bordered
  checkbox (default: checked). Excludes inactive/deselected primers from
  triggering duplicate sequence warnings.
- **Auto-activate new valid primer**: Bordered checkbox (default: unchecked).
  Automatically checks the active checkbox for newly added primers once valid
  sequences are entered.
- **Show primer temperature column**: Bordered checkbox (default: unchecked).
  Displays an extra $T_m$ column in the primer table showing live melting
  temperature calculations.
- **Primer temperature colour scheme**: Dropdown palette selector (`None`
  [default], `Cool-Warm`, `Traffic Light`) for background colour-coding of $T_m$
  values in the primer list.

## 3. Primer Melting Temperature (Tm) Settings

Configures thermodynamic conditions and calculation methods used for computing
primer melting temperatures ($T_m$).

- **Tm Calculation Method**: Dropdown selector for the algorithm:
  - **SantaLucia 1998 / Owczarzy 2008 (Default)**: Nearest-neighbor
    thermodynamics with salt corrections for monovalent and divalent cations.
  - **Lander / Amplify 4**: Empirical formula from the original Amplify4
    program.
  - *(Note: If SantaLucia thermodynamics encounters unsupported non-standard or
    degenerate bases, it automatically falls back to Lander / Amplify 4).*
- **Reaction Concentrations**:
  - **DNA Conc (nM)**: Total primer strand concentration in nM (default:
    `50.0`).
  - **DNA Pol Conc**: Polymerase concentration (default: `0.0`).
  - **Monovalent Salt Conc (mM)**: Concentration of monovalent cations ($Na^+$,
    $K^+$, $Tris^+$) in mM (default: `50.0`).
  - **Divalent Salt Conc (mM)**: Concentration of divalent cations ($Mg^{2+}$)
    in mM (default: `1.5`).
  - **dNTP Conc (mM)**: Concentration of dNTPs in mM (default: `0.0`).

## 4. PCR Settings

Configures pairwise nucleotide scoring weights, primability, stability
thresholds, and compatibility modes for PCR replication calculations.

- **Base Pair Weights Matrix**:
  - A $15 \\times 4$ interactive table of pairwise weights for primer
    nucleotides (`A`, `T`, `C`, `G`, IUPAC degenerate codes `R`, `Y`, `S`, `W`,
    `K`, `M`, `B`, `D`, `H`, `V`, `N`) vs template nucleotides (`A`, `T`, `C`,
    `G`).
  - Controls how base-pairing quality is scored during primer binding site
    search.
- **Parameters**:
  - **Primability Cutoff**: Minimum fractional primability score required for a
    match site to be recognized (default: `0.8`).
  - **Stability Cutoff**: Minimum fractional 3' stability score required for a
    match site to be recognized (default: `0.4`).
  - **Amplify4 Compatibility Mode**: Bordered checkbox (default: unchecked).
    Enforces legacy Amplify4 scoring logic and weighting constraints.
  - **Improved Primer Binding Site Visualisation**: Bordered checkbox (default:
    checked). Enables enhanced visual match diagrams and bond strength symbols
    in PCR view cards.

## 5. Primer Dimer Settings

Configures scoring and thresholds for evaluating primer self-dimers and
cross-dimers in the Primer Dimers View and info panels.

- **Primer Dimer Weights Matrix**:
  - A $15 \\times 15$ interactive table defining pairwise scoring weights for
    all primer base vs primer base interactions (including IUPAC ambiguous
    bases).
  - Default values are derived from Amplify4 thermodynamic binding weights.
- **Parameters**:
  - **Min Overlap**: Minimum overlap length (in base pairs) required for a
    primer dimer alignment to be evaluated (default: `3`).
  - **Threshold**: Minimum overall quality score threshold for a dimer
    interaction to be reported in the [Primer Dimers View](primer_dimer_view.md)
    (default: `60.0`).

## 6. Designer 2D Settings

Configures display options for pair-wise 2D primer truncation analysis in
Designer 2D View.

- **2D Results Grid colour map scheme**: Dropdown palette selector (`None`,
  `Cool-Warm`, `Traffic Light`, `Blue-Orange` [default]). Controls the heat map
  colour scheme for displaying amplicon scores across 2D primer length
  combination grids.

## 7. Diagnostics Settings

Configures application logging, log file output, log level thresholds, and log
file rotation parameters.

- **Console Output**: Bordered checkbox (default: checked). Toggles logging
  messages to standard console output.
- **AmplifyP Log Level**: Dropdown setting the log level threshold for
  application logic (`DEBUG`, `INFO` [default], `WARNING`, `ERROR`, `CRITICAL`).
- **Flet Log Level**: Dropdown setting log level threshold for the underlying UI
  framework (`INFO` [default], `WARNING`, `ERROR`, `CRITICAL`).
- **File Logging**: Bordered checkbox (default: checked on desktop platforms,
  hidden/disabled on web). Enables writing log messages to a disk file.
- **Log File Path**: Dropdown selector with options:
  - `(Default)`: Writes to OS-specific application log directory (e.g.
    `~/.config/amplifyp/app.log` or `%APPDATA%\AmplifyP\app.log`).
  - `Select folder`: Launches a directory picker dialog to select a custom
    destination folder for `app.log`.
- **Log Rotation**: Bordered checkbox (default: checked on desktop). Enables
  rotating log files when size limits are reached.
- **Max Log Size (MB)**: Dropdown selector for max file size before rotation
  (`1`, `5` [default], `10`, `20`, `50` MB).

## Reset & Local Persistence

- **Reset to Default Button**: Located at the bottom of the Settings View.
  Clicking **Reset to Default** resets all configuration settings across all
  tiles to factory defaults, saves the updated state to local storage,
  reconfigures logging, and refreshes the application UI.
- **Storage Locations**:
  - **Desktop Platforms (Linux, macOS, Windows)**: Settings are persisted in
    YAML format at the OS-specific config path (`settings.yaml` under
    `$XDG_CONFIG_HOME/amplifyp/`, `~/Library/Application Support/AmplifyP/`, or
    `%APPDATA%\AmplifyP\`).
  - **Web Platform**: Settings are saved directly to the browser's local
    storage.

______________________________________________________________________

[Return to GUI Manual Index](README.md)
