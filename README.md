[![AmplifyP post-commit](https://github.com/fangfufu/AmplifyP/actions/workflows/workflow.yml/badge.svg)](https://github.com/fangfufu/AmplifyP/actions/workflows/workflow.yml)
[![CodeQL](https://github.com/fangfufu/AmplifyP/actions/workflows/codeql.yml/badge.svg)](https://github.com/fangfufu/AmplifyP/actions/workflows/codeql.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/fangfufu/AmplifyP/stable.svg)](https://results.pre-commit.ci/latest/github/fangfufu/AmplifyP/stable)
[![codecov](https://codecov.io/gh/fangfufu/AmplifyP/graph/badge.svg?token=UNEJRZSVPJ)](https://codecov.io/gh/fangfufu/AmplifyP)
[![CodeFactor](https://www.codefactor.io/repository/github/fangfufu/amplifyp/badge)](https://www.codefactor.io/repository/github/fangfufu/amplifyp)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f3bffb3752794a728b3722120ca267fa)](https://app.codacy.com/gh/fangfufu/AmplifyP/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# AmplifyP

AmplifyP is a modern, high-performance Python rewrite of William Engels's
classic [Amplify4](https://github.com/wrengels/Amplify4) tool for simulating
**Polymerase Chain Reaction (PCR)**, with a web version available on
[GitHub Pages](https://fangfufu.github.io/AmplifyP/). It allows molecular
biologists and researchers to predict DNA amplification products (amplicons)
from a template sequence and primer set, accurately factoring in primability,
stability, and melting properties of primer binding sites.

______________________________________________________________________

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Web Version (Static)](#web-version-static)
  - [Local Graphical User Interface (GUI)](#local-graphical-user-interface-gui)
  - [Python API](#python-api)
    - [Basic PCR Simulation](#basic-pcr-simulation)
    - [Primer Dimer Analysis](#primer-dimer-analysis)
    - [Thermodynamic Melting Calculations](#thermodynamic-melting-calculations)
- [Project Structure](#project-structure)
- [Development](#development)
- [Attribution](#attribution)
- [License](#license)

______________________________________________________________________

## Features

- **PCR Simulation**: Predict potential amplicons, product sequences, and
  lengths using rigorous primer-template binding models. Supports both linear
  and circular DNA templates.
- **Cross-Platform GUI**: Built with [Flet](https://flet.dev/) (a modern
  Flutter-based UI framework) providing a beautiful, fully responsive app for
  desktop and web browsers.
- **Scoring Engine**: Calculates primability and stability scores using highly
  customisable length-wise and pairwise weight tables.
- **State Serialisation**: Save and load templates, primers, cutoffs, and
  replication configurations seamlessly using YAML files.
- **Programmatic Python API**: A developer-friendly, fully typed API to run
  simulations, dimer analyses, or thermodynamic melting calculations within your
  pipelines.
- **Primer Dimer Analysis**: Detect and score primer-dimer formation risks
  between primer pairs, including self-dimers and cross-dimers.
- **Melting Temperature Calculations**: Compute thermodynamic melting properties
  of primer-template and primer-primer hybrids using nearest-neighbor models.

______________________________________________________________________

## Installation

AmplifyP requires Python 3.12 or higher.

First, clone the repository and set up your Python virtual environment under
`.venv` at the root of the repository:

```bash
git clone https://github.com/fangfufu/AmplifyP.git
cd AmplifyP

# Set up and activate virtual environment
python -m venv .venv
source .venv/bin/activate
```

Install the package and its runtime dependencies:

```bash
pip install .
```

For developer setup and running tests, please refer to the
[Development Guide](src/README.md).

______________________________________________________________________

## Usage

### Web Version (Static)

Run AmplifyP entirely in your browser without any local installation! Visit the
live static web app:
**[https://fangfufu.github.io/AmplifyP/](https://fangfufu.github.io/AmplifyP/)**

To build and test the static web app locally, or to run the application with
hot-reload for development, see the [Development Guide](src/README.md).

### Local Graphical User Interface (GUI)

You can launch the Flet GUI application locally as a standalone desktop app.

#### 1. Launching the App

- **Desktop Application (via Python)**:
  ```bash
  python src/main.py
  ```

#### 2. Features and Workflows

The GUI offers intuitive workflows to:

1. **Input template DNA sequences** and configure linear/circular topologies.
1. **Manage multiple primers** with individual sequences and custom names.
1. **Customise thresholds** for primability and stability cutoffs.
1. **Simulate PCR** to visualise products, lanes, and details of binding sites.
1. **Analyze primer dimers** to identify potential primer-dimer formation risks.
1. **Save/Load projects** to resume work easily via YAML configuration files.

### Python API

Integrate AmplifyP into your custom bioinformatics workflows using the Python
API.

#### Basic PCR Simulation

```python
from amplifyp.dna import DNA, Primer, DNAType
from amplifyp.repliconf import Repliconf
from amplifyp.amplicon import AmpliconGenerator
from amplifyp.settings import GLOBAL_REPLICATION_SETTINGS

# 1. Define your DNA template and primers
template_seq = "AGCT..."  # Replace with your actual sequence
template = DNA(template_seq, DNAType.LINEAR, name="MyTemplate")

primer_fwd = Primer("AGCT...", name="FwdPrimer")
primer_rev = Primer("TCGA...", name="RevPrimer")

# 2. Initialise the Amplicon Generator
generator = AmpliconGenerator(template)

# 3. Create Replication Configurations for each primer
# This step calculates potential binding sites on the template
conf_fwd = Repliconf(template, primer_fwd, GLOBAL_REPLICATION_SETTINGS)
conf_rev = Repliconf(template, primer_rev, GLOBAL_REPLICATION_SETTINGS)

# 4. Search for origins (binding sites)
conf_fwd.search()
conf_rev.search()

# 5. Add configurations to the generator
generator.add(conf_fwd)
generator.add(conf_rev)

# 6. Generate amplicons
amplicons = generator.get_amplicons()

for amp in amplicons:
    print(f"Product: {amp.product.name}")
    print(f"Sequence: {amp.product.seq}")
    print(f"Length: {len(amp.product)}")
    print(f"Quality Score: {amp.q_score}")
    print("-" * 20)
```

#### Primer Dimer Analysis

```python
from amplifyp.dna import DNA, Primer, DNAType
from amplifyp.repliconf import Repliconf
from amplifyp.dimer import DimerGenerator
from amplifyp.settings import GLOBAL_REPLICATION_SETTINGS

template = DNA("AGCT...", DNAType.LINEAR, name="Template")
primer_a = Primer("AGCT...", name="PrimerA")
primer_b = Primer("TCGA...", name="PrimerB")

# Create replication configurations
conf_a = Repliconf(template, primer_a, GLOBAL_REPLICATION_SETTINGS)
conf_b = Repliconf(template, primer_b, GLOBAL_REPLICATION_SETTINGS)
conf_a.search()
conf_b.search()

# Generate dimer analysis
dimer_gen = DimerGenerator()
dimer_gen.add(conf_a)
dimer_gen.add(conf_b)

dimers = dimer_gen.get_dimers()
for dimer in dimers:
    print(f"Dimer: {dimer.primer1.name}-{dimer.primer2.name}")
    print(f"Score: {dimer.score}")
    print(f"Product: {dimer.product.seq}")
```

#### Thermodynamic Melting Calculations

```python
from amplifyp.dna import DNA, Primer
from amplifyp.melting import MeltingCalculator

primer = Primer("AGCTTCGAA", name="TestPrimer")
template = DNA("TTTCTAGCTTCGAAAGGG", name="Template")

calculator = MeltingCalculator()

# Calculate melting temperature of primer-template hybrid
tm = calculator.calculate_tm(primer, template)
print(f"Melting Temperature: {tm:.2f}°C")
```

______________________________________________________________________

## Project Structure

```
AmplifyP/
├── pyproject.toml              # Package configuration and dependencies
├── README.md                   # This file
├── src/
│   ├── main.py                 # GUI application entry point
│   ├── amplifyp/               # Core package
│   │   ├── amplicon.py         # Amplicon generation logic
│   │   ├── dimer.py            # Primer dimer analysis
│   │   ├── dna.py              # DNA/Primer sequence classes
│   │   ├── errors.py           # Custom exception types
│   │   ├── gui/                # Flet-based GUI components
│   │   ├── melting.py          # Thermodynamic melting calculations
│   │   ├── origin.py           # Primer binding site detection
│   │   ├── pcr.py              # PCR simulation engine
│   │   ├── repliconf.py        # Replication configuration
│   │   └── settings.py         # Global replication settings
│   └── README.md               # Development guide
├── tests/                      # Unit and integration tests
└── scripts/                    # Utility scripts
```

______________________________________________________________________

## Development

For development setup, running the app with hot-reload, building the static web
version, and running tests, see the [Development Guide](src/README.md).

______________________________________________________________________

## Attribution

- **Amplify4**: This project is built upon the logic and methodology of the
  original [Amplify4](https://github.com/wrengels/Amplify4) software by William
  Engels. We preserve the simulation models and algorithms of the original while
  offering a modern, robust, and accessible cross-platform implementation.
- **Roboto Mono Font**: Licensed under the Apache License, Version 2.0.
  Copyright 2015 The Roboto Mono Project Authors.

______________________________________________________________________

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE)
file for details.
