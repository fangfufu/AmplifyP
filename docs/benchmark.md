# Benchmark Suite

AmplifyP includes a standalone benchmark and profiling script at
[`scripts/benchmark.py`](../scripts/benchmark.py) for measuring the performance
of core bioinformatics calculations, YAML state serialisation, settings caching,
and headless Flet GUI operations.

## What It Does

`scripts/benchmark.py` constructs deterministic synthetic DNA sequences and runs
isolated, repeatable timing loops using `time.perf_counter()`. For each
workload, it records both the **best** (minimum) and **mean** elapsed wall-clock
time in milliseconds across a configurable number of repetitions.

Key characteristics:

- **Deterministic synthetic inputs**: Helper `_make_dna()` uses a fixed-seed
  linear congruential generator (`seed` values `1`, `2`, `42`, `100+i`, `200+i`,
  `300+i`) to generate repeatable A/C/G/T templates, primer pools, and
  degenerate nucleotide sequences (`M, R, W, S, Y, K, N`) without external
  fixture files.
- **No GUI window required**: GUI benchmarks run headlessly against a mocked
  `ft.Page` (`width=1200`, `height=800`, `web=False`), isolating Python widget
  construction, state-key evaluation, layout solving, and editor updates.
- **Side-by-side baselines**: The `io` suite compares pure-Python
  `yaml.safe_load` / `yaml.safe_dump` directly against `parse_state_yaml()`,
  `serialise_state()`, and `ProjectState.to_yaml()` (which use PyYAML
  `CSafeLoader` / `CSafeDumper` when available).

## Benchmark Suites

### `core` Suite

| Workload                                          | Description                                                                                               |
| ------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| `PCR: 2 primers on 5 kb circular plasmid`         | End-to-end `PCR.predict_amplicons()` with two 20-nt primers on a 5 kb circular template                   |
| `PCR: 20 multiplex primers on 30 kb linear`       | Multiplex `PCR.predict_amplicons()` with 20 primers (22-nt) on a 30 kb linear template                    |
| `Repliconf.search: 21 primers on 30 kb linear`    | Isolated `Repliconf.search()` across 20 specific primers plus 1 degenerate primer (`ACGTMRWSYKVHDBNACGT`) |
| `Designer 1D: 23 steps (no template)`             | `PrimerDesigner1D` truncating a 40-nt primer down to 18 nt without template search                        |
| `Designer 1D: 23 steps + 5 kb circular`           | `PrimerDesigner1D` with origin counting on a 5 kb circular template                                       |
| `Designer 1D: 23 steps + 30 kb linear`            | `PrimerDesigner1D` with origin counting on a 30 kb linear template                                        |
| `Designer 2D: 23x23 pairs (no template)`          | `PrimerDesigner2D` evaluating 529 primer pairs without a template                                         |
| `Designer 2D: 23x23 pairs + 5 kb linear`          | `PrimerDesigner2D` with linear amplicon counting on a 5 kb linear template                                |
| `Designer 2D: 23x23 pairs + 5 kb circular`        | `PrimerDesigner2D` with circular amplicon counting on a 5 kb circular template                            |
| `calculate_dimer_stats: 10,000 pairs`             | Direct pairwise dimer alignment over 100×100 primers containing degenerate bases                          |
| `PrimerDimerGenerator: 100 primers (5,050 pairs)` | All-vs-all dimer discovery across 100 unique 25-nt primers                                                |
| `calculate_tm: 1,000 calls (SantaLucia/Owczarzy)` | 1,000 nearest-neighbour melting temperature calculations                                                  |
| `clean_sequence: 60 kb wrapped FASTA (100x)`      | Sequence whitespace/escape stripping on 60 kb line-wrapped input (100 calls)                              |
| `DNA ops: slice, add, reverse_complement (500x)`  | Isolated `DNA` slicing, concatenation (`+`), and `Primer.reverse_complement()` calls                      |

### `io` Suite

| Workload                                              | Description                                                                 |
| ----------------------------------------------------- | --------------------------------------------------------------------------- |
| `yaml.safe_dump (SafeDumper, 30 kb + 50 primers)`     | Pure-Python YAML serialisation baseline for a 30 kb template and 50 primers |
| `serialise_state() (CSafeDumper, 30 kb + 50 primers)` | Direct `_StateDumper` YAML serialisation on the pre-built state dictionary  |
| `ProjectState.to_yaml() (30 kb + 50 primers)`         | End-to-end `ProjectState.to_dict()` + `_StateDumper` serialisation          |
| `yaml.safe_load (SafeLoader, 30 kb + 50 primers)`     | Pure-Python YAML load baseline followed by `ProjectState.from_dict()`       |
| `parse_state_yaml() (30 kb + 50 primers)`             | Production YAML load using `CSafeLoader` when available                     |
| `GUISettings.get_replication_settings() (100x)`       | Cached retrieval of `ReplicationSettings` from `GUISettings` (100 calls)    |
| `GUISettings.get_primer_dimer_settings() (100x)`      | Cached retrieval of `PrimerDimerSettings` from `GUISettings` (100 calls)    |

### `gui` Suite

| Workload                                             | Description                                                                           |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------- |
| `SettingsView.__init__ (359 controls)`               | Constructing the full settings form and score matrix tables                           |
| `Grid2DResultsView._on_cell_click (100x, 169 cells)` | 100 consecutive cell selections on a 13×13 2D designer matrix                         |
| `PCRView.run_pcr (cold 30 kb)`                       | Initial `PCRView.run_pcr()` simulation and diagram layout/drawing on a 30 kb template |
| `PCRView.run_pcr (20x resize cache hit)`             | 20 repeated diagram redraws with `_cached_pcr` cache hit (simulating window resizes)  |
| `PCRView.open_all_cards (multiple amplicons)`        | Batch opening all predicted amplicon detail cards via `open_all_cards()`              |
| `DimerView.run_analysis (cold)`                      | Initial `DimerView.run_analysis()` execution                                          |
| `DimerView.run_analysis (20x cache hit)`             | 20 repeated `DimerView.run_analysis()` cache hits                                     |
| `update_line_numbers 30 kb (10x)`                    | Recomputing line-number gutter indices across a 30 kb wrapped template (10 calls)     |
| `TemplateInput._update_status_bar 30 kb (50x)`       | Updating the template status bar across 50 cursor positions on a 30 kb template       |

## How to Run It

Activate the repository virtual environment first:

```bash
source .venv/bin/activate
```

Run all benchmark suites (`core`, `io`, and `gui`) with default settings
(`repeats=3`):

```bash
python scripts/benchmark.py
```

### Command-Line Options

| Option                      | Default | Description                                                                                     |
| --------------------------- | ------- | ----------------------------------------------------------------------------------------------- |
| `--suite {all,core,io,gui}` | `all`   | Run a specific suite or all suites                                                              |
| `--repeats N`               | `3`     | Minimum timed repetitions per workload (fast micro-benchmarks run at least 5 times)             |
| `--profile`                 | `false` | Print top-12 `cProfile` hotspots sorted by `tottime` for 30 kb multiplex PCR and 5 kb 2D design |
| `--json`                    | `false` | Emit structured JSON output instead of a formatted text table                                   |

### Examples

Run only the core benchmarks with 5 repetitions:

```bash
python scripts/benchmark.py --suite core --repeats 5
```

Inspect top `cProfile` hotspots alongside core timings:

```bash
python scripts/benchmark.py --suite core --repeats 1 --profile
```

Export GUI benchmark timings as JSON:

```bash
python scripts/benchmark.py --suite gui --repeats 1 --json
```

## Notes on Interpreting Results

- Look at **Best (ms)** when comparing algorithmic changes on quiet hardware, as
  it minimises background OS scheduling noise; use **Mean (ms)** to observe
  consistency across runs.
- `PyYAML CSafeLoader available: True` at the top of the table confirms that
  libyaml C bindings are active. In browser WebAssembly (Pyodide) builds where
  libyaml is unavailable, `parse_state_yaml()` and `serialise_state()`
  automatically fall back to pure-Python `SafeLoader` and `SafeDumper`.
- Headless GUI timings reflect Python-side view, layout, and canvas generation
  work and do not include Flutter engine rasterisation or IPC latency.
- **GUI suite failure semantics**: The GUI analysis benchmarks (`PCRView` and
  `DimerView` workloads) assert that each analysis reports success. If any of
  them fails, the script raises a `RuntimeError`
  (`GUI benchmark analysis failed: {workload}`) and exits with a non-zero status
  instead of reporting timings for that suite.
