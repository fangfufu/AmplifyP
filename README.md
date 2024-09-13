[![AmplifyP post-commit](https://github.com/fangfufu/AmplifyP/actions/workflows/workflow.yml/badge.svg)](https://github.com/fangfufu/AmplifyP/actions/workflows/workflow.yml)
[![CodeQL](https://github.com/fangfufu/AmplifyP/actions/workflows/github-code-scanning/codeql/badge.svg)](https://github.com/fangfufu/AmplifyP/actions/workflows/github-code-scanning/codeql)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/fangfufu/AmplifyP/master.svg)](https://results.pre-commit.ci/latest/github/fangfufu/AmplifyP/master)
[![codecov](https://codecov.io/gh/fangfufu/AmplifyP/graph/badge.svg?token=UNEJRZSVPJ)](https://codecov.io/gh/fangfufu/AmplifyP)
[![CodeFactor](https://www.codefactor.io/repository/github/fangfufu/amplifyp/badge)](https://www.codefactor.io/repository/github/fangfufu/amplifyp)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f3bffb3752794a728b3722120ca267fa)](https://app.codacy.com/gh/fangfufu/AmplifyP/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

# AmplifyP
Python rewrite of  William Engels's
[Amplify4](https://github.com/wrengels/Amplify4) by two people who went to
Pembroke College Cambridge, hence AmplifyP.

This is pretty much a work in progress.

## How does it all work?
These are all notes for Fufu to understand.

### Replication configuration
A replication configuration (**Repliconf**) consists of a long strand
of DNA which acts as the **template**, and a **primer**.

The **search()** function enumerates all the possible replication origins, and
calculate the primability and stability score. For a potential replication
origin, if the primability and stability scores are above the threshold, the
index (location) of this potential replication origin is considered as valid.
It is recorded in **origin_idx.fwd** and **origin_idx.rev**. Note that during
the search, the template sequence is searched in both forward and reverse
direction. The process of search in the reverse direction involves reversing
the DNA template. (I will add more information about the reverse compliment
later.)

## Running the unit tests
Please invoke ``pytest`` at the root of the repository.
