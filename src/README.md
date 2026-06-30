# AmplifyP Developer & Contributor Guide

This document is the comprehensive guide on how to program, modify, and test the
**AmplifyP** application and its core library.

______________________________________________________________________

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

______________________________________________________________________

## 2. Codebase Architecture

The project directory is structured as follows:

```
AmplifyP/
├── docs/                       # Markdown user and API documentation
│   ├── images/                 # Visual assets (screenshots)
│   ├── gui_guide.md            # GUI application user manual
│   └── api_guide.md            # Python library programmatic usage guide
├── pyproject.toml              # Project metadata, dependencies, and tools configuration
├── README.md                   # Repository overview and quick start guide
├── src/
│   ├── README.md               # This document (Development Guide)
│   ├── main.py                 # Flet GUI application command-line entry point
│   └── amplifyp/               # Core Python library
│       ├── __init__.py         # Package entry and version definition
│       ├── dna.py              # DNA and Primer classes, IUB/IUPAC codes
│       ├── pcr.py              # PCR reaction simulation engine
│       ├── repliconf.py        # Replication configurations, padding, and caching
│       ├── origin.py           # Origin matching, stability/primability scoring algorithms
│       ├── amplicon.py         # Amplicon prediction and product sequence generation
│       ├── dimer.py            # Primer self-dimer and cross-dimer analysis
│       ├── melting.py          # Thermodynamic melting temperature calculations (Tm)
│       ├── settings.py         # Customisable simulation and biophysical constants
│       ├── errors.py           # Custom exception types
│       └── gui/                # Flet-based cross-platform GUI modules
└── tests/                      # Unit, integration, and Playwright end-to-end tests
```

______________________________________________________________________

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

______________________________________________________________________

## 4. Pyodide static builds (Client-Side Web)

To build the static web version that runs entirely client-side using Pyodide in
the browser (similar to what is hosted on GitHub Pages):

1. Compile the app using the static build script:

   ```bash
   ./build_static.sh
   ```

   *This compiles the Flet application and puts the output in the `dist`
   directory.*

1. Serve the static site locally to test it:

   ```bash
   python -m http.server 23455 -d dist
   ```

1. Open your browser and navigate to `http://localhost:23455`.

______________________________________________________________________

## 5. Code Quality & Verification

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

### Running the Test Suite

We use `pytest` for all unit, integration, and browser tests.

#### 1. Unit and Integration Tests

To run the fast, non-GUI tests:

```bash
pytest
```

#### 2. End-to-End (E2E) Browser Tests

To run the Playwright browser automation tests that verify the Flet web server's
behaviour:

```bash
# Ensure playwright browser binaries are installed
playwright install

# Run the slow E2E tests
pytest --run-slow -m e2e
```

______________________________________________________________________

## 6. Version Release Checklist

When preparing a new release of AmplifyP, you must update the version string in
the following files:

1. **`pyproject.toml`**: Update the `version` field under `[project]`
   ```toml
   [project]
   version = "1.0.0"
   ```
1. **`src/amplifyp/__init__.py`**: Update the `__version__` variable
   ```python
   __version__ = "1.0.0"
   ```
