#!/usr/bin/env python
# Copyright (C) 2026 AmplifyP Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Benchmark AmplifyP computational, YAML IO, and headless GUI workloads."""

from __future__ import annotations

import argparse
import cProfile
import io
import json
import pstats
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

_SRC_DIR = Path(__file__).resolve().parent.parent / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

import flet as ft  # noqa: E402
import yaml  # noqa: E402

from amplifyp.dimer import (  # noqa: E402
    PrimerDimerGenerator,
    calculate_dimer_stats,
)
from amplifyp.dna import (  # noqa: E402
    DNA,
    GLOBAL_COMPLEMENT_TABLE,
    DNAType,
    Primer,
    clean_sequence,
    format_sequence,
)
from amplifyp.gui.settings import GUISettings  # noqa: E402
from amplifyp.gui.user_data import GUIInput  # noqa: E402
from amplifyp.gui.views.designer_2d.grid_2d_results import (  # noqa: E402
    Grid2DResultsView,
)
from amplifyp.gui.views.dimer.dimer_view import DimerView  # noqa: E402
from amplifyp.gui.views.input.input_view import InputView  # noqa: E402
from amplifyp.gui.views.input.template.formatter import (  # noqa: E402
    update_line_numbers,
)
from amplifyp.gui.views.pcr.pcr_view import PCRView  # noqa: E402
from amplifyp.gui.views.settings.settings_view import SettingsView  # noqa: E402
from amplifyp.io import (  # noqa: E402
    ProjectState,
    ReactionInput,
    parse_state_yaml,
    serialise_state,
)
from amplifyp.melting import calculate_tm  # noqa: E402
from amplifyp.pcr import PCR  # noqa: E402
from amplifyp.primer_designer_1d import PrimerDesigner1D  # noqa: E402
from amplifyp.primer_designer_2d import PrimerDesigner2D  # noqa: E402
from amplifyp.repliconf import Repliconf  # noqa: E402

_BASES = "ACGT"
_DEGEN_ALPHABET = "A" * 20 + "C" * 20 + "G" * 20 + "T" * 20 + "MRWSYKN"


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    """Measured runtime for a single benchmark workload."""

    suite: str
    name: str
    best_ms: float
    mean_ms: float
    repeats: int


def _make_dna(length: int, seed: int = 42, alphabet: str = _BASES) -> str:
    """Generate a deterministic synthetic DNA sequence via a linear LCG."""
    state = seed & 0xFFFFFFFF
    chars: list[str] = []
    append = chars.append
    alpha_len = len(alphabet)
    for _ in range(length):
        state = (state * 1664525 + 1013904223) & 0xFFFFFFFF
        append(alphabet[(state >> 16) % alpha_len])
    return "".join(chars)


def _measure(
    suite: str,
    name: str,
    fn: Callable[[], Any],
    repeats: int,
) -> BenchmarkResult:
    """Run ``fn`` ``repeats`` times and record best and mean timings in ms."""
    samples: list[float] = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        samples.append((time.perf_counter() - t0) * 1000.0)
    return BenchmarkResult(
        suite=suite,
        name=name,
        best_ms=min(samples),
        mean_ms=sum(samples) / len(samples),
        repeats=repeats,
    )


def _run_core_benchmarks(
    repeats: int, *, with_profile: bool = False
) -> tuple[list[BenchmarkResult], list[tuple[str, str]]]:
    """Benchmark PCR simulation, primer designers, dimers, Tm, and DNA."""
    results: list[BenchmarkResult] = []
    profiles: list[tuple[str, str]] = []

    seq_5kb = _make_dna(5_000, seed=1)
    seq_30kb = _make_dna(30_000, seed=2)
    tpl_5kb_lin = DNA(seq_5kb, DNAType.LINEAR, "pTest5kbLin")
    tpl_5kb_circ = DNA(seq_5kb, DNAType.CIRCULAR, "pTest5kb")
    tpl_30kb_lin = DNA(seq_30kb, DNAType.LINEAR, "Lambda30kb")

    fwd_40 = Primer(seq_5kb[500:540], name="Fwd40")
    rev_40_seq = seq_5kb[1800:1840][::-1].translate(GLOBAL_COMPLEMENT_TABLE)
    rev_40 = Primer(rev_40_seq, name="Rev40")

    primers_20: list[Primer] = []
    for idx in range(10):
        pos_f = 1000 + idx * 2500
        pos_r = pos_f + 800
        primers_20.append(Primer(seq_30kb[pos_f : pos_f + 22], name=f"F{idx}"))
        rc = seq_30kb[pos_r : pos_r + 22][::-1].translate(
            GLOBAL_COMPLEMENT_TABLE
        )
        primers_20.append(Primer(rc, name=f"R{idx}"))

    def run_pcr_5kb() -> None:
        p = PCR(tpl_5kb_circ)
        p.add_primer(Primer(fwd_40.seq[:20], name="F"))
        p.add_primer(Primer(rev_40.seq[-20:], name="R"))
        p.predict_amplicons()

    def run_pcr_30kb_multiplex() -> None:
        p = PCR(tpl_30kb_lin)
        p.add_primers(primers_20)
        p.predict_amplicons()

    results.append(
        _measure(
            "core",
            "PCR: 2 primers on 5 kb circular plasmid",
            run_pcr_5kb,
            max(repeats, 5),
        )
    )
    results.append(
        _measure(
            "core",
            "PCR: 20 multiplex primers on 30 kb linear",
            run_pcr_30kb_multiplex,
            repeats,
        )
    )

    primers_21 = [*primers_20, Primer("ACGTMRWSYKVHDBNACGT", name="Degen1")]

    def run_repliconf_21() -> None:
        for primer in primers_21:
            Repliconf(tpl_30kb_lin, primer).search()

    results.append(
        _measure(
            "core",
            "Repliconf.search: 21 primers on 30 kb linear",
            run_repliconf_21,
            repeats,
        )
    )

    if with_profile:
        pr = cProfile.Profile()
        pr.enable()
        run_pcr_30kb_multiplex()
        pr.disable()
        buf = io.StringIO()
        pstats.Stats(pr, stream=buf).sort_stats("tottime").print_stats(12)
        profiles.append(("PCR: 20 multiplex primers on 30 kb", buf.getvalue()))

    results.append(
        _measure(
            "core",
            "Designer 1D: 23 steps (no template)",
            lambda: PrimerDesigner1D(fwd_40, min_length=18),
            max(repeats, 5),
        )
    )
    results.append(
        _measure(
            "core",
            "Designer 1D: 23 steps + 5 kb circular",
            lambda: PrimerDesigner1D(
                fwd_40, min_length=18, template=tpl_5kb_circ
            ),
            repeats,
        )
    )
    dna_30kb_primer = DNA(seq_30kb[1000:1040])
    results.append(
        _measure(
            "core",
            "Designer 1D: 23 steps + 30 kb linear",
            lambda: PrimerDesigner1D(
                dna_30kb_primer, min_length=18, template=tpl_30kb_lin
            ),
            repeats,
        )
    )

    results.append(
        _measure(
            "core",
            "Designer 2D: 23x23 pairs (no template)",
            lambda: PrimerDesigner2D(fwd_40, 18, rev_40, 18),
            repeats,
        )
    )
    results.append(
        _measure(
            "core",
            "Designer 2D: 23x23 pairs + 5 kb linear",
            lambda: PrimerDesigner2D(
                fwd_40, 18, rev_40, 18, template=tpl_5kb_lin
            ),
            repeats,
        )
    )
    results.append(
        _measure(
            "core",
            "Designer 2D: 23x23 pairs + 5 kb circular",
            lambda: PrimerDesigner2D(
                fwd_40, 18, rev_40, 18, template=tpl_5kb_circ
            ),
            repeats,
        )
    )

    if with_profile:
        pr2 = cProfile.Profile()
        pr2.enable()
        PrimerDesigner2D(fwd_40, 18, rev_40, 18, template=tpl_5kb_circ)
        pr2.disable()
        buf2 = io.StringIO()
        pstats.Stats(pr2, stream=buf2).sort_stats("tottime").print_stats(12)
        profiles.append(
            ("Designer 2D: 23x23 pairs + 5 kb circular", buf2.getvalue())
        )

    pool_100 = [
        Primer(_make_dna(25, seed=100 + idx), name=f"P{idx}")
        for idx in range(100)
    ]
    degen_pool_100 = [
        Primer(
            _make_dna(25, seed=300 + idx, alphabet=_DEGEN_ALPHABET),
            name=f"D{idx}",
        )
        for idx in range(100)
    ]

    def run_dimer_stats_10k() -> None:
        seqs = [(p.seq_upper, len(p)) for p in degen_pool_100]
        for s1, n1 in seqs:
            for s2, n2 in seqs:
                calculate_dimer_stats(s1, s2, n1, n2)

    def run_dimer_pool() -> None:
        gen = PrimerDimerGenerator()
        for primer in pool_100:
            gen.add_primer(primer)
        gen.analyse_primers()

    results.append(
        _measure(
            "core",
            "calculate_dimer_stats: 10,000 pairs",
            run_dimer_stats_10k,
            repeats,
        )
    )
    results.append(
        _measure(
            "core",
            "PrimerDimerGenerator: 100 primers (5,050 pairs)",
            run_dimer_pool,
            repeats,
        )
    )
    results.append(
        _measure(
            "core",
            "calculate_tm: 1,000 calls (SantaLucia/Owczarzy)",
            lambda: [calculate_tm(p) for p in pool_100 * 10],
            max(repeats, 5),
        )
    )

    wrapped_60kb = "\n".join(
        seq_30kb[i : i + 80] for i in range(0, len(seq_30kb), 80)
    )
    wrapped_60kb = f"{wrapped_60kb}\n{wrapped_60kb}"
    results.append(
        _measure(
            "core",
            "clean_sequence: 60 kb wrapped FASTA (100x)",
            lambda: [clean_sequence(wrapped_60kb) for _ in range(100)],
            repeats,
        )
    )

    def run_dna_ops_500() -> None:
        for _ in range(500):
            mid = tpl_30kb_lin[100:600]
            _ = fwd_40 + mid + rev_40.reverse_complement()

    results.append(
        _measure(
            "core",
            "DNA ops: slice, add, reverse_complement (500x)",
            run_dna_ops_500,
            max(repeats, 5),
        )
    )

    return results, profiles


def _run_io_benchmarks(repeats: int) -> list[BenchmarkResult]:
    """Benchmark YAML serialisation/parsing and GUISettings caching."""
    seq_30kb = _make_dna(30_000, seed=2)
    primers = [
        {"name": f"P{i}", "seq": _make_dna(25, seed=200 + i), "active": True}
        for i in range(50)
    ]
    state = ProjectState(
        input=ReactionInput(
            template=seq_30kb,
            template_circular=False,
            primers=primers,
        )
    )
    state_dict = state.to_dict()
    yaml_str = state.to_yaml()
    gui_settings = GUISettings()

    return [
        _measure(
            "io",
            "yaml.safe_dump (SafeDumper, 30 kb + 50 primers)",
            lambda: yaml.safe_dump(state_dict, sort_keys=False),
            repeats,
        ),
        _measure(
            "io",
            "serialise_state() (CSafeDumper, 30 kb + 50 primers)",
            lambda: serialise_state(state_dict),
            repeats,
        ),
        _measure(
            "io",
            "ProjectState.to_yaml() (30 kb + 50 primers)",
            lambda: state.to_yaml(),
            repeats,
        ),
        _measure(
            "io",
            "yaml.safe_load (SafeLoader, 30 kb + 50 primers)",
            lambda: ProjectState.from_dict(yaml.safe_load(yaml_str)),
            repeats,
        ),
        _measure(
            "io",
            "parse_state_yaml() (30 kb + 50 primers)",
            lambda: parse_state_yaml(yaml_str),
            repeats,
        ),
        _measure(
            "io",
            "GUISettings.get_replication_settings() (100x)",
            lambda: [
                gui_settings.get_replication_settings() for _ in range(100)
            ],
            repeats,
        ),
        _measure(
            "io",
            "GUISettings.get_primer_dimer_settings() (100x)",
            lambda: [
                gui_settings.get_primer_dimer_settings() for _ in range(100)
            ],
            repeats,
        ),
    ]


def _run_gui_benchmarks(repeats: int) -> list[BenchmarkResult]:
    """Benchmark GUI views using a mocked headless Flet page."""
    results: list[BenchmarkResult] = []
    page = MagicMock(spec=ft.Page)
    page.width = 1200
    page.height = 800
    page.web = False
    settings = GUISettings()

    results.append(
        _measure(
            "gui",
            "SettingsView.__init__ (359 controls)",
            lambda: SettingsView(page, settings),
            repeats,
        )
    )

    fwd_p = Primer("ACGTACGTACGTACGTACGTACGT")
    rev_p = Primer("TGCATGCATGCATGCATGCATGCA")
    designer = PrimerDesigner2D(
        fwd_p,
        12,
        rev_p,
        12,
        template=None,
        replication_settings=settings.get_replication_settings(),
    )
    grid = Grid2DResultsView(
        settings=settings,
        page=page,
        on_select_step_callback=lambda _step: None,
    )
    grid.update_grid(designer)
    grid_keys = list(grid._cell_containers.keys())
    steps = designer.all_steps

    def click_100_cells() -> None:
        n_keys = len(grid_keys)
        for i in range(100):
            grid._on_cell_click(steps[0], grid_keys[i % n_keys])

    results.append(
        _measure(
            "gui",
            f"Grid2DResultsView._on_cell_click (100x, {len(grid_keys)} cells)",
            click_100_cells,
            repeats,
        )
    )

    seq_30kb = _make_dna(30_000, seed=42)
    p1 = seq_30kb[500:520]
    p2 = seq_30kb[2500:2520][::-1].translate(GLOBAL_COMPLEMENT_TABLE)
    input_data = GUIInput(
        template=seq_30kb,
        template_circular=False,
        primers=[
            {"name": "p1", "seq": p1, "active": True},
            {"name": "p2", "seq": p2, "active": True},
        ],
    )

    pcr_view = PCRView(page, input_data, settings)
    pcr_view.run_pcr()

    def run_pcr_cold() -> None:
        pcr_view._cached_state_key = None
        pcr_view.run_pcr()

    results.append(
        _measure("gui", "PCRView.run_pcr (cold 30 kb)", run_pcr_cold, repeats)
    )
    results.append(
        _measure(
            "gui",
            "PCRView.run_pcr (20x resize cache hit)",
            lambda: [
                pcr_view.run_pcr(keep_cards=True, refresh_cards=False)
                for _ in range(20)
            ],
            repeats,
        )
    )

    multi_primers: list[dict[str, Any]] = []
    for idx in range(5):
        f_start = 1000 + idx * 4000
        r_start = f_start + 500
        f_seq = seq_30kb[f_start : f_start + 20]
        r_seq = seq_30kb[r_start : r_start + 20][::-1].translate(
            GLOBAL_COMPLEMENT_TABLE
        )
        multi_primers.append({"name": f"MF{idx}", "seq": f_seq, "active": True})
        multi_primers.append({"name": f"MR{idx}", "seq": r_seq, "active": True})
    cards_input = GUIInput(
        template=seq_30kb,
        template_circular=False,
        primers=multi_primers,
    )
    cards_view = PCRView(page, cards_input, settings)
    cards_view.run_pcr()

    def run_open_all_cards() -> None:
        cards_view.result_list.controls = [
            c
            for c in cards_view.result_list.controls
            if not isinstance(c, ft.Card)
        ]
        cards_view.open_all_cards()

    results.append(
        _measure(
            "gui",
            "PCRView.open_all_cards (multiple amplicons)",
            run_open_all_cards,
            repeats,
        )
    )

    dimer_view = DimerView(page, input_data, settings)
    dimer_view.run_analysis()

    def run_dimer_cold() -> None:
        dimer_view._cached_state_key = None
        dimer_view.run_analysis()

    results.append(
        _measure(
            "gui", "DimerView.run_analysis (cold)", run_dimer_cold, repeats
        )
    )
    results.append(
        _measure(
            "gui",
            "DimerView.run_analysis (20x cache hit)",
            lambda: [dimer_view.run_analysis() for _ in range(20)],
            repeats,
        )
    )

    input_view = InputView(page, input_data, settings)
    tmpl_input = input_view.template_input
    tmpl_input.template_sequence.value = format_sequence(seq_30kb, 60)
    tmpl_input._cleaned_len = len(seq_30kb)
    tmpl_input._is_focused = True

    def run_line_numbers_10x() -> None:
        for _ in range(10):
            update_line_numbers(tmpl_input, update=False)

    results.append(
        _measure(
            "gui",
            "update_line_numbers 30 kb (10x)",
            run_line_numbers_10x,
            repeats,
        )
    )

    def cursor_50_moves() -> None:
        for idx in range(50):
            pos = 1000 + idx * 100
            tmpl_input._update_status_bar(
                ft.TextSelection(base_offset=pos, extent_offset=pos),
                update=False,
            )

    results.append(
        _measure(
            "gui",
            "TemplateInput._update_status_bar 30 kb (50x)",
            cursor_50_moves,
            repeats,
        )
    )

    return results


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line options for the benchmark runner."""
    parser = argparse.ArgumentParser(
        description="Benchmark AmplifyP computational, IO, and GUI workloads.",
    )
    parser.add_argument(
        "--suite",
        choices=("all", "core", "io", "gui"),
        default="all",
        help="Benchmark suite to run (default: all).",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Minimum number of timed repetitions per workload (default: 3).",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Print cProfile breakdowns for multiplex PCR and 2D design.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Emit benchmark measurements as JSON.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run selected benchmark suites and format output."""
    args = _parse_args(argv)
    repeats = max(1, args.repeats)
    all_results: list[BenchmarkResult] = []
    profiles: list[tuple[str, str]] = []

    if args.suite in ("all", "core"):
        core_results, profiles = _run_core_benchmarks(
            repeats, with_profile=args.profile
        )
        all_results.extend(core_results)
    if args.suite in ("all", "io"):
        all_results.extend(_run_io_benchmarks(repeats))
    if args.suite in ("all", "gui"):
        all_results.extend(_run_gui_benchmarks(repeats))

    if args.json_output:
        payload = {
            "csafeloader_available": hasattr(yaml, "CSafeLoader"),
            "results": [
                {
                    "suite": r.suite,
                    "name": r.name,
                    "best_ms": round(r.best_ms, 4),
                    "mean_ms": round(r.mean_ms, 4),
                    "repeats": r.repeats,
                }
                for r in all_results
            ],
        }
        print(json.dumps(payload, indent=2))
        return 0

    print(f"PyYAML CSafeLoader available: {hasattr(yaml, 'CSafeLoader')}")
    print("=" * 80)
    print(f"{'Workload':<54} {'Best (ms)':>10} {'Mean (ms)':>10}")
    print("-" * 80)
    current_suite = ""
    for res in all_results:
        if res.suite != current_suite:
            current_suite = res.suite
            print(f"[{current_suite.upper()}]")
        print(f"  {res.name:<52} {res.best_ms:10.2f} {res.mean_ms:10.2f}")

    for title, stats_text in profiles:
        print("\n" + "=" * 80)
        print(f"cProfile top hotspots: {title}")
        print(stats_text.rstrip())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
