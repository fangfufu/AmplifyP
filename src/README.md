# AmplifyP Developers' and Contributors' Guide

This document is the comprehensive guide on how to program, modify, and test the
**AmplifyP** application and its core library.

## 1. Development Installation

AmplifyP requires Python 3.12 or higher.

To set up a local development environment, clone the repository, activate your
virtual environment, and install all runtime and optional dependencies:

```bash
# 1. Create and source virtual environment
python -m venv .venv
source .venv/bin/activate

# 2. Install package in editable mode with development/e2e dependencies
pip install -e ".[dev,e2e]"
```

## 2. Codebase Architecture

The project directory is structured as follows:

```
AmplifyP/
├── docs/                       # Markdown user and API documentation
│   ├── api_guide.md            # Python library programmatic usage guide
│   ├── gui/
│   │   ├── images/             # Visual assets (screenshots)
│   │   └── manual.md           # GUI application user manual
│   ├── history.md              # Project history and background
│   └── windows_setup.md        # Windows development and installer guide
├── pyproject.toml              # Project metadata, dependencies, and tools configuration
├── README.md                   # Repository overview and quick start guide
├── src/
│   ├── README.md               # This document (Development Guide)
│   ├── main.py                 # Flet GUI application command-line entry point
│   └── amplifyp/               # Core Python library
│       ├── __init__.py         # Package entry and version definition
│       ├── amplicon.py         # Amplicon prediction and product sequence generation
│       ├── dimer.py            # Primer self-dimer and cross-dimer analysis
│       ├── dir_idx.py          # Directional index calculation utilities
│       ├── dna.py              # DNA and Primer classes, IUB/IUPAC codes
│       ├── errors.py           # Custom exception types
│       ├── melting.py          # Thermodynamic melting temperature calculations (Tm)
│       ├── origin.py           # Origin matching, stability/primability scoring algorithms
│       ├── pcr.py              # PCR reaction simulation engine
│       ├── primer_designer_1d.py # 1D primer design and analysis engine
│       ├── primer_designer_2d.py # 2D primer pair design and matrix analysis engine
│       ├── py.typed            # PEP 561 marker for inline type annotations
│       ├── repliconf.py        # Replication configurations, padding, and caching
│       ├── settings.py         # Customisable simulation and biophysical constants
│       └── gui/                # Flet-based cross-platform GUI modules
└── tests/                      # Unit, integration, GUI, and end-to-end tests
    └── gui/                    # Flet GUI tests
```

## 3. Running the GUI with Hot-Reload

When developing the user interface, you can run the application using the **Flet
CLI** to enable **hot-reload**. This automatically refreshes the interface in
real-time as you edit and save the source files.

### Desktop Application with Hot-Reload

Runs the application as a standalone desktop window:

```bash
flet run -r src/main.py
```

*Note: The `-r` (or `--recursive`) flag is required so Flet watches changes
inside the `src/amplifyp` subdirectory.*

### Web Application with Hot-Reload

Runs the application as a local web app in your browser:

```bash
flet run -w -r src/main.py
```

To run it on a specific port (e.g. port `34521`):

```bash
flet run -w -r -p 34521 src/main.py
```

## 4. Command Line Arguments

`src/main.py` provides several command-line arguments for launching the
application, pre-loading state files, exporting screenshots, and configuring
window dimensions or web mode.

```bash
python src/main.py [options]
```

| Argument                   | Type    | Description                                                        |
| :------------------------- | :------ | :----------------------------------------------------------------- |
| `-h`, `--help`             | Flag    | Display the help message and exit.                                 |
| `-f`, `--state <path>`     | String  | Path to a YAML state file to load on startup.                      |
| `--auto-close`             | Flag    | Automatically quit after rendering completes (requires `--state`). |
| `-s`, `--screenshots`      | Flag    | Save PNG screenshots of views (requires `--state`).                |
| `--screenshots-dir <dir>`  | String  | Target directory for saved PNG screenshots.                        |
| `--window-width <pixels>`  | Integer | Set application window width in pixels.                            |
| `--window-height <pixels>` | Integer | Set application window height in pixels.                           |
| `--web`                    | Flag    | Launch in web browser mode.                                        |

## 5. Pyodide static builds (Client-Side Web)

To build the static web version that runs entirely client-side using Pyodide in
the browser (similar to what is hosted on GitHub Pages):

1. Compile the app using the static build script:

   ```bash
   ./build_web.sh
   ```

   *This compiles the Flet application and puts the output in the `dist`
   directory.*

2. Serve the static site locally to test it:

   ```bash
   python -m http.server 23455 -d dist
   ```

3. Open your browser and navigate to `http://localhost:23455`.

## 6. Code Quality & Verification

To maintain code quality, AmplifyP uses several code style, quality, and testing
checks.

### Git Hook Checks (`prek`)

Always run the `prek` checks to verify all code meets style, formatting, and
type-checking rules before making a commit. Run:

```bash
prek run --all-files
```

This runs:

- **`ruff`** for linting and code formatting checks.
- **`mypy`** for strict static type-checking.
- **`vulture`** for detecting unused code.
- **`yamlfmt`** for formatting configuration YAML files.
- **`typos`** for identifying spelling errors.
- **`mdformat`** for consistent markdown formatting.

### Static Type Checking with Pyright

Pyright runs automatically in CI when pushing to the `dev` branch. To reduce
local execution delays, it is not included in the default `prek` pre-commit
hooks.

You can run `pyright` manually before pushing:

```bash
pyright
```

### Running the Test Suite

We use `pytest` for all unit, integration, and browser tests.

#### 1. Unit and Integration Tests

To run the fast, non-GUI tests:

```bash
pytest
```

#### 2. End-to-End (E2E) Browser Tests

End-to-end (E2E) tests are not run automatically by `prek` local hooks to save
time. You must run them manually.

To run the Playwright browser automation tests that verify the Flet web server's
behaviour:

```bash
# Ensure playwright browser binaries are installed
playwright install

# Run the E2E tests
pytest -m e2e
```

*(Note: On headless CI environments like GitHub Actions runners, these tests are
executed under `xvfb`)*

## 7. Automated Release Process

The release process is automated via GitHub Actions using
[release-please](https://github.com/googleapis/release-please). Version bumps in
`pyproject.toml` and `src/amplifyp/__init__.py` are handled automatically based
on conventional commits.
