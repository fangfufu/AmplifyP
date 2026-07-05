# AmplifyP Documentation

Welcome to the documentation directory for **AmplifyP**, a modern Python rewrite
of William Engels's classic **Amplify4** program for simulating Polymerase Chain
Reactions (PCR) and analysing primer bindings.

This directory contains guides and instructions tailored for both users
(biologists and researchers) and developers/contributors.

______________________________________________________________________

## Table of Contents

1. [AmplifyP GUI User Guide](gui_guide.md) — For biologists and researchers
   using the graphical interface.
1. [AmplifyP Python API Guide](api_guide.md) — For developers and
   bioinformaticians using AmplifyP programmatically as a library.
1. [Windows Development Environment Setup Guide](windows_setup.md) — For setting
   up the development, testing, and packaging environment on Windows.

______________________________________________________________________

## 1. [GUI User Guide](gui_guide.md)

The **[GUI User Guide](gui_guide.md)** is a comprehensive manual for the desktop
and web-based graphical user interface. It explains how to interact with the
three primary workspaces:

- **Substrate View**: How to input and manage your target DNA sequence
  (topology, case styling, direction) and configure your active and custom
  primers.
- **PCR Results View**: How to interpret the visual PCR simulation map (binding
  site marks, product lines) and review predicted amplicon details.
- **Primer Dimers View**: How to analyse self-dimer and cross-dimer
  configurations, including scoring and binding alignments.
- **Project Management**: Saving and loading simulation states using standard
  YAML configuration files.

Refer to the [GUI User Guide](gui_guide.md) to get started with the visual tool.

______________________________________________________________________

## 2. [Python API Guide](api_guide.md)

The **[Python API Guide](api_guide.md)** details how to use `amplifyp` as a
programmatic library within Python scripts or bioinformatics pipelines. It
covers:

- **DNA & Primer Classes** (`amplifyp.dna`): Creating DNA and Primer objects,
  handling IUPAC degenerate bases, and calculating redundancy folds.
- **Simulating PCR** (`amplifyp.pcr`): Running PCR simulations programmatically
  to find binding sites and predict amplicons.
- **Biophysical Scoring Engine** (`amplifyp.origin`): How matching, primability,
  and stability scores are calculated using biophysical weight tables.
- **Melting Calculations** (`amplifyp.melting`): Computing thermodynamic melting
  temperatures ($T_m$) using nearest-neighbour models and salt corrections.
- **Primer Dimer Analysis** (`amplifyp.dimer`): Querying self-dimers and
  cross-dimers between primers to identify potential dimerisation risks.

Refer to the [Python API Guide](api_guide.md) for full code examples and API
details.

______________________________________________________________________

## 3. [Windows Setup Guide](windows_setup.md)

The **[Windows Development Environment Setup Guide](windows_setup.md)** is a
step-by-step developer setup guide for Windows. It includes:

- **Prerequisites**: Installing Python (3.12+), Git, and JRSoftware Inno Setup.
- **Environment Setup**: Cloning the repository, configuring the `.venv` virtual
  environment, and installing development dependencies.
- **Running & Debugging**: Running the Flet GUI application with hot-reload in
  desktop or web mode.
- **Testing**: Executing unit and integration tests with `pytest`, and running
  Playwright end-to-end (E2E) tests.
- **Packaging**: Building standalone executables using PyInstaller and
  generating the setup installer with Inno Setup.

Refer to the [Windows Setup Guide](windows_setup.md) for platform-specific
development workflows.

______________________________________________________________________

## Directory Assets

- **[images/](images/)**: Visual assets and screenshots used throughout the user
  and developer documentation guides.
