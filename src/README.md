# AmplifyP Development Guide

This document contains information for developers working on AmplifyP.

## Development Installation

AmplifyP requires Python 3.12 or higher.

To install development dependencies (testing and E2E automation):

```bash
pip install -e ".[dev,e2e]"
```

## Running the App with Hot-Reload

Using the **Flet CLI** (`flet run`) enables **hot-reload**, which automatically
updates the app interface in real-time as you modify the source files.

### Launching the App

- **Desktop Application with Hot-Reload (via Flet CLI)**:
  ```bash
  flet run -r src/main.py
  ```
- **Dynamic Web Application with Hot-Reload (via Flet CLI)**: To run the app as
  a dynamic local website in your web browser with hot-reloading:
  ```bash
  flet run -w -r src/main.py
  ```
  To specify a custom port (e.g., port `43425`):
  ```bash
  flet run -w -r -p 43425 src/main.py
  ```

### How Hot-Reload Works

When launching the application using the Flet CLI (`flet run`):

- **File Watcher**: The CLI monitors the files in the directory containing
  `main.py` for any changes.
- **Live Updates**: Saving changes to your code automatically updates the active
  application window or browser tab without needing to manually restart the
  process.
- **Recursive Watching**: Because the project's core source files are located in
  the `src/amplifyp` sub-directory, you must include the `-r` (or `--recursive`)
  flag to ensure changes to files in sub-directories are also detected.

## Local Pyodide Build & Test

To build and test the static web app locally using Pyodide:

```bash
./build_static.sh
python -m http.server 23455 -d dist
```

Then navigate to <http://localhost:23455> in your web browser.

## Development & Testing

Run the test suite using `pytest` to verify correctness:

```bash
# Run unit and integration tests
pytest
```

## Updating the Version Number

When releasing a new version of AmplifyP, you must update the version string in
two locations:

1. **`pyproject.toml`**: Update the `version` field under `[project]` (e.g.,
   `version = "v0.6.7"`).
1. **`src/amplifyp/__init__.py`**: Update the `__version__` variable (e.g.,
   `__version__ = "v0.6.7"`).
