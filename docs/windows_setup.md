# Windows Development Environment Setup Guide

This guide provides step-by-step instructions for setting up the developer and
contributor environment for **AmplifyP** on Windows. It covers installing
prerequisites, setting up Python virtual environments, running the application
with hot-reload, executing unit and end-to-end (E2E) tests, and compiling
release packages.

______________________________________________________________________

## 1. Prerequisites

Before setting up the project, ensure the following software is installed on
your Windows machine:

### A. Python (3.12 or Higher)

AmplifyP requires Python 3.12 or higher.

- **Installation via Winget**: Open PowerShell or Command Prompt and run:
  ```powershell
  winget install Python.Python.3.12
  ```
- **Manual Installation**: Download and run the installer from the official
  [Python Website](https://www.python.org/downloads/).
  > [!IMPORTANT]
  > During the Python installation process, ensure the checkbox for **"Add
  > Python to PATH"** (or similar) is ticked. This is necessary to access Python
  > from the command line.

### B. Git

Required for cloning the repository and managing version control.

- **Installation via Winget**:
  ```powershell
  winget install Git.Git
  ```
- **Manual Installation**: Download and install Git from
  [Git for Windows](https://gitforwindows.org/).

### C. Inno Setup (Optional - for compiling the installer)

Required only if you wish to build the standalone installer executable (`.exe`).

- **Installation via Winget**:
  ```powershell
  winget install JRSoftware.InnoSetup
  ```
- **Installation via Chocolatey**:
  ```powershell
  choco install innosetup -y
  ```

______________________________________________________________________

## 2. Cloning the Repository

Open PowerShell, Command Prompt, or Git Bash, and navigate to the directory
where you want to store the project. Run:

```powershell
git clone https://github.com/fangfufu/AmplifyP.git
cd AmplifyP
```

______________________________________________________________________

## 3. Creating and Activating the Virtual Environment

Set up a Python virtual environment under `.venv` at the root of the repository.
Follow the instructions for your preferred terminal shell:

### PowerShell

```powershell
# Create the virtual environment
python -m venv .venv

# Activate the virtual environment
.venv\Scripts\Activate.ps1
```

> [!NOTE]
> If you receive an execution policy error (e.g.
> `Script execution is disabled on this system`), run the following command to
> allow running scripts in the current process:
>
> ```powershell
> Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope Process
> ```

### Command Prompt (cmd.exe)

```cmd
:: Create the virtual environment
python -m venv .venv

:: Activate the virtual environment
.venv\Scripts\activate.bat
```

### Git Bash

```bash
# Create the virtual environment
python -m venv .venv

# Activate the virtual environment
source .venv/Scripts/activate
```

______________________________________________________________________

## 4. Installing Dependencies

With the virtual environment activated, install the package in editable mode
along with development and end-to-end testing dependencies:

```powershell
pip install -e ".[dev,e2e]"
```

______________________________________________________________________

## 5. Running the Application

### Running the Python Script Directly

To launch the graphical user interface (GUI) from source:

```powershell
python src/main.py
```

To launch the web version locally at port `34521`:

```powershell
python src/main.py --web
```

### Running with Hot-Reload (Flet CLI)

To enable hot-reload (so that UI changes are applied automatically when files
are saved), run the application using the Flet CLI:

- **Desktop application**:
  ```powershell
  flet run -r src/main.py
  ```
- **Web application**:
  ```powershell
  flet run -w -r src/main.py
  ```
- **Web application on a specific port** (e.g. port `34521`):
  ```powershell
  flet run -w -r -p 34521 src/main.py
  ```

______________________________________________________________________

## 6. Verification and Code Quality

To verify your changes and maintain high code quality, use the following tools:

### A. Pre-commit Verification (`prek`)

Before committing your work, run the quality checks using `prek`. This ensures
that your code complies with formatting, linting, and type-checking rules:

```powershell
prek run --all-files
```

This runs:

- **`ruff`** for linting and code formatting checks.
- **`mypy`** for strict static type-checking.
- **`vulture`** for detecting unused code.
- **`yamlfmt`** for formatting configuration YAML files.
- **`typos`** for identifying spelling errors.
- **`mdformat`** for consistent markdown formatting.

### B. Running Tests

We use `pytest` for unit, integration, and browser tests.

#### 1. Unit and Integration Tests

To run the fast, non-GUI tests:

```powershell
pytest
```

#### 2. End-to-End (E2E) Browser Tests

To run the Playwright browser automation tests that verify the web interface's
behaviour:

```powershell
# Ensure playwright browser binaries are installed
playwright install

# Run the slow E2E tests
pytest --run-slow -m e2e
```

______________________________________________________________________

## 7. Packaging and Building (Windows Releases)

To package the application and generate release assets on Windows, you can run
the `build_windows.ps1` script located at the root of the repository. This
script automates generating the Git SHA, building the Flet binary, archiving the
build as a ZIP file, and compiling the setup installer using Inno Setup.

Open PowerShell and run:

```powershell
.\build_windows.ps1
```

### Output Artefacts

After the script runs successfully, the following files and directories will be
created:

- **`build\AmplifyP`**: The compiled Flet Windows application folder.
- **`amplifyp-windows-<version>.zip`**: A ZIP archive of the application
  directory.
- **`amplifyp-windows-setup-<version>.exe`**: The standalone Windows
  installation wizard (only if Inno Setup's `iscc` compiler is installed and
  available in the system PATH).
