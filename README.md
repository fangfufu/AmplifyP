[![AmplifyP post-commit](https://github.com/fangfufu/AmplifyP/actions/workflows/workflow.yml/badge.svg)](https://github.com/fangfufu/AmplifyP/actions/workflows/workflow.yml)
[![CodeQL](https://github.com/fangfufu/AmplifyP/actions/workflows/github-code-scanning/codeql/badge.svg)](https://github.com/fangfufu/AmplifyP/actions/workflows/github-code-scanning/codeql)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/fangfufu/AmplifyP/master.svg)](https://results.pre-commit.ci/latest/github/fangfufu/AmplifyP/master)
[![codecov](https://codecov.io/gh/fangfufu/AmplifyP/graph/badge.svg?token=UNEJRZSVPJ)](https://codecov.io/gh/fangfufu/AmplifyP)
[![CodeFactor](https://www.codefactor.io/repository/github/fangfufu/amplifyp/badge)](https://www.codefactor.io/repository/github/fangfufu/amplifyp)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f3bffb3752794a728b3722120ca267fa)](https://app.codacy.com/gh/fangfufu/AmplifyP/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

# AmplifyP

AmplifyP is a Python rewrite of William Engels's
[Amplify4](https://github.com/wrengels/Amplify4), a tool for simulating
Polymerase Chain Reaction (PCR). It allows users to predict amplification
products (amplicons) from a given DNA template and a set of primers, taking into
 account the primability and stability of primer binding sites.

## Features

- **PCR Simulation**: Predict potential amplicons based on primer binding
properties.
- **Scoring System**: Calculates primability and stability scores for primer
binding sites using customizable weight tables.
- **GUI Application**: Includes a Tkinter-based graphical user interface for
easy interaction.
- **Python API**: Provides a flexible API for programmatic access to simulation
tools.

## Installation

To use AmplifyP, clone the repository and ensure you have a compatible Python
environment (Python 3.12+).

You probably want to set up your Python virtual environment first:
```
python -m venv venv
source venv/bin/activate
```

Then install AmplifyP:
```bash
git clone https://github.com/fangfufu/AmplifyP.git
cd AmplifyP
pip install .
```

## Usage

### Graphical User Interface

To launch the GUI:

```bash
python -m amplifyp.gui
```

The GUI allows you to:
1. Input a template DNA sequence.
2. Add multiple primers.
3. Configure primability and stability cutoffs.
4. Simulate PCR to view potential amplicons.
5. Analyze individual primers to see all potential binding sites.

### Python API

You can also use AmplifyP as a library in your Python scripts.

```python
from amplifyp.dna import DNA, Primer, DNAType
from amplifyp.repliconf import Repliconf
from amplifyp.amplicon import AmpliconGenerator
from amplifyp.settings import DEFAULT_SETTINGS

# 1. Define your DNA template and primers
template_seq = "AGCT..."  # Replace with your actual sequence
template = DNA(template_seq, DNAType.LINEAR, name="MyTemplate")

primer_fwd = Primer("AGCT...", name="FwdPrimer")
primer_rev = Primer("TCGA...", name="RevPrimer")

# 2. Initialize the Amplicon Generator
generator = AmpliconGenerator(template)

# 3. Create Replication Configurations for each primer
# This step calculates potential binding sites on the template
conf_fwd = Repliconf(template, primer_fwd, DEFAULT_SETTINGS)
conf_rev = Repliconf(template, primer_rev, DEFAULT_SETTINGS)

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

## Development

### Running Tests

To run the unit tests, install `pytest` and run it from the root of the
repository:

```bash
pip install pytest pytest-cov
pytest
```

### Todo

- [ ] Verification of circular target amplification

## Attribution

This project is based on the logic and methodology of
[Amplify4](https://github.com/wrengels/Amplify4) by William Engels. We aim to
preserve the accuracy of the original simulations while providing a modern
Python implementation.
