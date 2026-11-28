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

"""Primer melting temperature calculation module.

This module provides functions to calculate the melting temperature (Tm) of
DNA primers using the Nearest-Neighbour thermodynamics model, including IUPAC
degenerate base support and min/mean/max Tm range estimation.
"""

import math
from dataclasses import dataclass
from enum import StrEnum
from itertools import pairwise, product
from typing import Final

from .dna import GLOBAL_COMPLEMENT_TABLE, Primer
from .errors import InsufficientThermodynamicDataError
from .settings import (
    GLOBAL_TM_SETTINGS,
    TMSettings,
)

# Thermodynamic parameters (SantaLucia 1998)
# Enthalpy (dH) in cal/mol, Entropy (dS) in cal/(K*mol)
# Keys are dinucleotides (5'->3')
# Source: SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and
# oligonucleotide DNA
# nearest-neighbour thermodynamics". PNAS, 95(4), 1460-1465.
NN_THERMO_DATA: Final[dict[str, tuple[float, float]]] = {
    "AA": (-7900, -22.2),
    "TT": (-7900, -22.2),
    "AT": (-7200, -20.4),
    "TA": (-7200, -21.3),
    "CA": (-8500, -22.7),
    "TG": (-8500, -22.7),
    "GT": (-8400, -22.4),
    "AC": (-8400, -22.4),
    "CT": (-7800, -21.0),
    "AG": (-7800, -21.0),
    "GA": (-8200, -22.2),
    "TC": (-8200, -22.2),
    "CG": (-10600, -27.2),
    "GC": (-9800, -24.4),
    "GG": (-8000, -19.9),
    "CC": (-8000, -19.9),
}

# Optimisation: Pre-compute tuple keys for faster lookup in loops
# pairwise(seq) produces tuples of characters, e.g. ('A', 'A')
_NN_THERMO_DATA_TUPLE: Final[dict[tuple[str, str], tuple[float, float]]] = {
    (k[0], k[1]): v for k, v in NN_THERMO_DATA.items()
}

# Standard IUPAC DNA base expansions covering Nucleotides.PRIMER
IUPAC_BASES: Final[dict[str, tuple[str, ...]]] = {
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
    "M": ("A", "C"),
    "R": ("A", "G"),
    "W": ("A", "T"),
    "S": ("C", "G"),
    "Y": ("C", "T"),
    "K": ("G", "T"),
    "V": ("A", "C", "G"),
    "H": ("A", "C", "T"),
    "D": ("A", "G", "T"),
    "B": ("C", "G", "T"),
    "N": ("A", "C", "G", "T"),
}

_STANDARD_BASES: Final[frozenset[str]] = frozenset("ACGT")

_IUPAC_GC_FRACTION: Final[dict[str, float]] = {
    base: sum(1.0 for c in expanded if c in "GC") / len(expanded)
    for base, expanded in IUPAC_BASES.items()
}

_TERMINAL_INIT_IUPAC: Final[dict[str, tuple[float, float]]] = {
    base: (
        sum(100.0 if c in "GC" else 2300.0 for c in expanded) / len(expanded),
        sum(-2.8 if c in "GC" else 4.1 for c in expanded) / len(expanded),
    )
    for base, expanded in IUPAC_BASES.items()
}

# Mean nearest-neighbour parameters for dinucleotides containing at least one
# degenerate IUPAC base. Canonical ACGT pairs stay in _NN_THERMO_DATA_TUPLE so
# tests patching _NN_THERMO_DATA_TUPLE are not shadowed.
_IUPAC_NN_THERMO_DATA_TUPLE: Final[
    dict[tuple[str, str], tuple[float, float]]
] = {
    (b1, b2): (
        sum(_NN_THERMO_DATA_TUPLE[(c1, c2)][0] for c1 in exp1 for c2 in exp2)
        / (len(exp1) * len(exp2)),
        sum(_NN_THERMO_DATA_TUPLE[(c1, c2)][1] for c1 in exp1 for c2 in exp2)
        / (len(exp1) * len(exp2)),
    )
    for b1, exp1 in IUPAC_BASES.items()
    for b2, exp2 in IUPAC_BASES.items()
    if b1 not in _STANDARD_BASES or b2 not in _STANDARD_BASES
}

DEFAULT_MAX_VARIANTS: Final[int] = 4096


@dataclass(frozen=True, slots=True)
class TmRange:
    """Melting temperature range across IUPAC primer variants.

    Attributes:
        min_tm: Minimum melting temperature in degrees Celsius.
        mean_tm: Mean melting temperature in degrees Celsius.
        max_tm: Maximum melting temperature in degrees Celsius.
        is_approximate: True if ``min_tm`` and ``max_tm`` were estimated using
            ``O(L)`` dynamic-programming bounds because ``variant_count``
            exceeded ``max_variants``.
        variant_count: Number of concrete IUPAC variants represented (``1`` when
            non-degenerate or using a non-expanding model).
    """

    min_tm: float
    mean_tm: float
    max_tm: float
    is_approximate: bool = False
    variant_count: int = 1


def _calculate_tm_sl_owczarzy_core(
    seq: str,
    settings: TMSettings,
    log_redundancy_fold: float = 0.0,
    is_degenerate: bool = False,
) -> float:
    """Compute SantaLucia 1998 / Owczarzy 2008 Tm from an uppercase IUPAC seq.

    Degenerate bases use mean nearest-neighbour and terminal initiation
    parameters across their IUPAC expansion, expected GC fraction for divalent
    salt correction, and equimolar per-variant strand concentration
    (``dna_conc / redundancy_fold`` via ``log_redundancy_fold``).
    """
    n = len(seq)
    if n < 2:
        return 0.0

    # Constants
    R: Final[float] = 1.987  # cal/(K*mol)

    # Initiation with terminal corrections (SantaLucia 1998, averaged for IUPAC)
    # Init w/ term G.C: dH=0.1 kcal, dS=-2.8 eu
    # Init w/ term A.T: dH=2.3 kcal, dS=4.1 eu
    left_dh, left_ds = _TERMINAL_INIT_IUPAC[seq[0]]
    right_dh, right_ds = _TERMINAL_INIT_IUPAC[seq[-1]]
    dh: float = left_dh + right_dh
    ds: float = left_ds + right_ds

    # Symmetry correction for self-complementary canonical sequences only
    is_self_comp = (not is_degenerate) and (
        seq == seq[::-1].translate(GLOBAL_COMPLEMENT_TABLE)
    )
    if is_self_comp:
        ds += -1.4
    conc_factor = 1.0 if is_self_comp else 4.0

    # Nearest neighbour steps
    for dinuc in pairwise(seq):
        val = _NN_THERMO_DATA_TUPLE.get(dinuc)
        if val is None:
            val = _IUPAC_NN_THERMO_DATA_TUPLE.get(dinuc)
        if val is not None:
            dh += val[0]
            ds += val[1]

    # Salt Correction (Owczarzy 2008)
    # References:
    # Owczarzy et al. (2008). "Predicting Stability of DNA Duplexes in
    # Solutions Containing Magnesium and Monovalent Cations". Biochemistry,
    # 47(19), 5336-5353.

    # 1. Calculate concentrations
    # Monovalent: Na+, K+, Tris+ (mM -> M)
    # Clamp salt inputs to non-negative values
    mono_mM = max(0.0, settings.monovalent_salt_conc)
    mono_M = mono_mM * 1e-3

    # Divalent: Mg2+ (mM -> M)
    # Deduct dNTPs which chelate Mg2+ in a 1:1 ratio to find free Mg2+
    div_mM = max(
        0.0,
        max(0.0, settings.divalent_salt_conc) - max(0.0, settings.dntp_conc),
    )
    div_M = div_mM * 1e-3

    # 2. Determine mode (Monovalent only, Mixed, or Divalent dominant)
    # Ratio R = sqrt([Mg2+]) / [Mon+]
    if math.isclose(mono_M, 0.0, abs_tol=1e-9):
        ratio = 999.0  # Large number, Divalent dominant
    else:
        ratio = math.sqrt(div_M) / mono_M

    # Calculate Tm_1M_Na (Kelvin) using equimolar per-species concentration
    total_dna_conc_M = settings.dna_conc * 1e-9
    if total_dna_conc_M <= 0:
        total_dna_conc_M = 50e-9

    log_conc_term = (
        math.log(total_dna_conc_M / conc_factor) - log_redundancy_fold
    )
    denom_1M = ds + R * log_conc_term
    if math.isclose(denom_1M, 0.0, abs_tol=1e-9):  # pragma: no cover
        raise InsufficientThermodynamicDataError(
            "Denominator is zero in Tm calculation"
        )

    if math.isclose(dh, 0.0, abs_tol=1e-9):
        raise InsufficientThermodynamicDataError(
            "Invalid sequence: lacks standard thermodynamic base pairs"
        )

    tm_1m_K = dh / denom_1M
    if math.isclose(tm_1m_K, 0.0, abs_tol=1e-9):
        raise InsufficientThermodynamicDataError("Calculated base Tm is zero")

    # Now apply corrections
    if math.isclose(div_M, 0.0, abs_tol=1e-9):
        # Monovalent only (SantaLucia 1998)
        if mono_M > 0:
            ds_corr = 0.368 * (n - 1) * math.log(mono_M)
            denom_corr = ds + ds_corr + R * log_conc_term
            if math.isclose(denom_corr, 0.0, abs_tol=1e-9):
                raise InsufficientThermodynamicDataError(
                    "Denominator with salt correction is zero"
                )
            tm_final_K = dh / denom_corr
        else:
            tm_final_K = tm_1m_K

    elif ratio < 0.22:
        # Monovalent dominant (use SantaLucia correction with [Mon+])
        ds_corr = 0.368 * (n - 1) * math.log(mono_M)
        denom_corr = ds + ds_corr + R * log_conc_term
        if math.isclose(denom_corr, 0.0, abs_tol=1e-9):
            raise InsufficientThermodynamicDataError(
                "Denominator with salt correction is zero"
            )
        tm_final_K = dh / denom_corr

    else:
        # Mixed (0.22 <= R < 6.0) or Magnesium Dominant (R >= 6.0)
        # Owczarzy 2008 Eq 16 coefficients
        a = 3.92e-5
        b = -9.11e-6
        c = 6.26e-5
        d = 1.42e-5
        e = -4.82e-4
        f = 5.25e-4
        g = 8.31e-5

        # Expected fraction of GC across IUPAC bases
        fgc = sum(_IUPAC_GC_FRACTION[base] for base in seq) / n
        log_mg = math.log(div_M)

        # Correction term
        # 1/Tm(Mg) = 1/Tm(1M) + a + b*ln(Mg) + fgc*(c + d*ln(Mg)) +
        #            (1/(2*(N-1))) * (e + f*ln(Mg) + g*(ln(Mg))^2)
        corr = (
            a
            + b * log_mg
            + fgc * (c + d * log_mg)
            + (1 / (2 * (n - 1))) * (e + f * log_mg + g * (log_mg**2))
        )

        # Note: If Mixed mode, Owczarzy 2008 sometimes specifies different
        # coeffs or eq. But Eq 16 is often used generally for Mg presence in
        # simplified implementations. Since exact mixed mode coeffs are
        # complex/variable, we stick to Eq 16 which accounts for Mg effects
        # well.
        tm_inv = (1.0 / tm_1m_K) + corr
        if math.isclose(tm_inv, 0.0, abs_tol=1e-9):
            raise InsufficientThermodynamicDataError(
                "Inverse Tm with salt correction is zero"
            )
        tm_final_K = 1.0 / tm_inv

    return tm_final_K - 273.15


def _select_dp_sequence(seq: str, t_ref_k: float, pick_max: bool) -> str:
    """Select an extremal canonical sequence over IUPAC choices via O(L) DP.

    When ``pick_max`` is False, minimises free energy ``dH - t_ref_k * dS``
    (yielding the strongest-binding / highest-Tm variant). When True,
    maximises free energy (yielding the weakest-binding / lowest-Tm variant).
    """
    n = len(seq)
    term_dg: dict[str, float] = {
        base: dh - t_ref_k * ds
        for base, (dh, ds) in _TERMINAL_INIT_IUPAC.items()
    }
    first_choices = IUPAC_BASES[seq[0]]
    scores: dict[str, float] = {b: term_dg[b] for b in first_choices}
    backpointers: list[dict[str, str]] = []

    for idx in range(1, n):
        choices = IUPAC_BASES[seq[idx]]
        next_scores: dict[str, float] = {}
        bp: dict[str, str] = {}
        is_last = idx == n - 1
        for curr in choices:
            end_dg = term_dg[curr] if is_last else 0.0
            best_prev = ""
            best_val = float("-inf") if pick_max else float("inf")
            for prev, prev_score in scores.items():
                dh, ds = _NN_THERMO_DATA_TUPLE[(prev, curr)]
                cand = prev_score + (dh - t_ref_k * ds) + end_dg
                better = cand > best_val if pick_max else cand < best_val
                if better:
                    best_val = cand
                    best_prev = prev
            next_scores[curr] = best_val
            bp[curr] = best_prev
        scores = next_scores
        backpointers.append(bp)

    end_base = (max if pick_max else min)(scores, key=lambda b: scores[b])
    chosen = [end_base]
    curr_base = end_base
    for bp in reversed(backpointers):
        curr_base = bp[curr_base]
        chosen.append(curr_base)
    chosen.reverse()
    return "".join(chosen)


def calculate_tm_santalucia_1998_owczarzy_2008(
    primer: Primer, settings: TMSettings = GLOBAL_TM_SETTINGS
) -> float:
    """Calculate the melting temperature (Tm) of a primer sequence.

    Uses the Nearest-Neighbour model with SantaLucia 1998 thermodynamic
    parameters and Owczarzy 2008 salt corrections. Supports all 15 standard
    IUPAC DNA bases (``A``, ``C``, ``G``, ``T``, ``M``, ``R``, ``W``, ``S``,
    ``Y``, ``K``, ``V``, ``H``, ``D``, ``B``, ``N``) by averaging nearest-
    neighbour and terminal initiation parameters across each degenerate
    position and dividing ``dna_conc`` by the primer redundancy fold.

    Formula used:
        Tm = (dH / (dS + R * ln(Ct/4))) - 273.15

    Where:
        dH: Total enthalpy (cal/mol)
        dS: Total entropy (cal/(K*mol))
        R: Gas constant (1.987 cal/(K*mol))
        Ct: Per-species oligonucleotide concentration (M)

    Salt correction is applied to the entropy term (dS) as per SantaLucia 1998:
        dS_corrected = dS_1M + 0.368 * (N - 1) * ln([Na+])

    Args:
        primer: The primer object containing the sequence.
        settings: Melting settings containing concentrations.

    Returns:
        The melting temperature in degrees Celsius.
    """
    seq = primer.seq_upper
    if len(seq) < 2:
        return 0.0

    if any(base not in IUPAC_BASES for base in seq):
        raise InsufficientThermodynamicDataError(
            "Sequence contains non-IUPAC bases for which "
            "nearest-neighbour parameters are not available"
        )

    is_degenerate = any(base not in _STANDARD_BASES for base in seq)
    log_redundancy_fold = (
        sum(math.log(len(IUPAC_BASES[base])) for base in seq)
        if is_degenerate
        else 0.0
    )
    return _calculate_tm_sl_owczarzy_core(
        seq,
        settings,
        log_redundancy_fold=log_redundancy_fold,
        is_degenerate=is_degenerate,
    )


def calculate_tm_lander_amplify4(
    primer: Primer,
    settings: TMSettings = GLOBAL_TM_SETTINGS,
) -> float:
    """Calculate Tm using the original Amplify4 algorithm.

    This method is a direct port of the `calcTm` method from the Swift
    codebase (Primer.swift). It uses its own set of entropy and enthalpy
    tables (typically 5x5 matrices) and specific correction factors.

    Args:
        primer: The primer object containing the sequence.
        settings: Melting settings containing concentrations and matrices.

    Returns:
        The melting temperature in degrees Celsius. Returns 0.0 if the
        sequence length is less than 1.
    """
    seq = primer.seq_upper
    seq_len = len(seq)
    if seq_len < 1:
        return 0.0

    entropy = settings.entropy
    enthalpy = settings.enthalpy

    entr: float = 108.0
    enth: float = 0.0

    # Sum neighbours
    # Note: entropy/enthalpy tables in Swift are accessed as [y][x]
    # where x is current base, y is next base.
    for x, y in pairwise(seq):
        try:
            entr += entropy[y, x]
            enth += enthalpy[y, x]
        except KeyError as exc:
            raise InsufficientThermodynamicDataError(
                "Sequence contains bases without Amplify4 thermodynamic "
                f"parameters: {(x, y)}"
            ) from exc

    # Scaling applied in Swift code
    entr = -entr * 0.1
    enth = -enth * 0.1

    # Corrections
    # DNAConc in settings is usually nM (default 50).
    # Swift formula: 1.987 * log(DNAConc/4e9)
    # If DNAConc=50, 50/4e9 = 1.25e-8 => 12.5 nM concentration assumption?
    # Or maybe it treats input as raw number?
    # We use settings.DNAConc exactly as Swift does.
    dna_conc_val = settings.dna_conc
    if dna_conc_val <= 0:
        dna_conc_val = 50.0
    log_dna = 1.987 * math.log(dna_conc_val / 4.0e9)

    # Salt: 16.6 * log(saltConc/1000) / log(10.0)
    # saltConc is typically mM (default 50).
    salt_conc_val = settings.monovalent_salt_conc
    if salt_conc_val <= 0:
        salt_conc_val = 50.0  # prevent log error if 0

    log_salt = 16.6 * math.log10(salt_conc_val / 1000.0)

    denom = entr + log_dna
    if math.isclose(denom, 0.0, abs_tol=1e-9):
        raise InsufficientThermodynamicDataError(
            "Denominator is zero in Tm calculation (Amplify4)"
        )

    tm = (enth * 1000.0) / denom - 273.15 + log_salt
    return tm


class TmMethod(StrEnum):
    """Supported melting temperature calculation methods."""

    SANTALUCIA_1998_OWCZARZY_2008 = "SantaLucia 1998 / Owczarzy 2008 (Default)"
    LANDER_AMPLIFY4 = "Lander / Amplify 4"


def calculate_tm(
    primer: Primer,
    settings: TMSettings = GLOBAL_TM_SETTINGS,
    method: TmMethod | str = TmMethod.SANTALUCIA_1998_OWCZARZY_2008,
    fallback_to_amplify4: bool = True,
) -> float:
    """Calculate the melting temperature (Tm) of a primer.

    Args:
        primer: The primer object containing the sequence.
        settings: Melting settings containing concentrations and matrices.
        method: Calculation method to use. Accepts a `TmMethod` member or its
            string value.
        fallback_to_amplify4: If True and the SantaLucia/Owczarzy model raises
            `InsufficientThermodynamicDataError` (for example, on unsupported
            non-IUPAC bases), fall back to `calculate_tm_lander_amplify4`.

    Returns:
        The melting temperature in degrees Celsius.
    """
    resolved_method = TmMethod(method)
    if resolved_method == TmMethod.LANDER_AMPLIFY4:
        return calculate_tm_lander_amplify4(primer, settings)
    try:
        return calculate_tm_santalucia_1998_owczarzy_2008(primer, settings)
    except InsufficientThermodynamicDataError:
        if fallback_to_amplify4:
            return calculate_tm_lander_amplify4(primer, settings)
        raise


def calculate_tm_range(
    primer: Primer,
    settings: TMSettings = GLOBAL_TM_SETTINGS,
    method: TmMethod | str = TmMethod.SANTALUCIA_1998_OWCZARZY_2008,
    fallback_to_amplify4: bool = True,
    max_variants: int = DEFAULT_MAX_VARIANTS,
) -> TmRange:
    """Calculate min, mean, and max Tm across concrete IUPAC variants.

    For the SantaLucia/Owczarzy model, concrete ``ACGT`` variants are
    enumerated exhaustively when the IUPAC redundancy fold is at most
    ``max_variants``, with each variant evaluated at equimolar per-species
    concentration (``dna_conc / variant_count``). Above ``max_variants``,
    ``mean_tm`` comes from the mean IUPAC approximation while
    ``min_tm``/``max_tm`` are estimated via ``O(L)`` weakest/strongest-base
    dynamic-programming bounds.

    For ``TmMethod.LANDER_AMPLIFY4`` (or fallback to Amplify 4), a single
    point estimate is returned with ``variant_count=1``.

    Args:
        primer: The primer object containing the sequence.
        settings: Melting settings containing concentrations and matrices.
        method: Calculation method to use.
        fallback_to_amplify4: If True and the SantaLucia/Owczarzy model raises
            ``InsufficientThermodynamicDataError``, fall back to Amplify 4.
        max_variants: Maximum IUPAC redundancy fold for exhaustive enumeration.

    Returns:
        A ``TmRange`` containing min, mean, max temperatures and metadata.

    Raises:
        ValueError: If ``max_variants`` is < 1, or ``method`` is unknown.
    """
    if max_variants < 1:
        raise ValueError("max_variants must be at least 1")

    resolved_method = TmMethod(method)
    if resolved_method == TmMethod.LANDER_AMPLIFY4:
        tm = calculate_tm_lander_amplify4(primer, settings)
        return TmRange(
            min_tm=tm,
            mean_tm=tm,
            max_tm=tm,
            is_approximate=False,
            variant_count=1,
        )

    seq = primer.seq_upper
    if len(seq) < 2:
        variant_count = (
            math.prod(len(IUPAC_BASES[b]) for b in seq)
            if all(b in IUPAC_BASES for b in seq)
            else 1
        )
        return TmRange(
            min_tm=0.0,
            mean_tm=0.0,
            max_tm=0.0,
            is_approximate=False,
            variant_count=variant_count,
        )

    try:
        if any(base not in IUPAC_BASES for base in seq):
            raise InsufficientThermodynamicDataError(
                "Sequence contains non-IUPAC bases for which "
                "nearest-neighbour parameters are not available"
            )

        is_degenerate = any(base not in _STANDARD_BASES for base in seq)
        if not is_degenerate:
            tm = _calculate_tm_sl_owczarzy_core(
                seq,
                settings,
                log_redundancy_fold=0.0,
                is_degenerate=False,
            )
            return TmRange(
                min_tm=tm,
                mean_tm=tm,
                max_tm=tm,
                is_approximate=False,
                variant_count=1,
            )

        variant_count = math.prod(len(IUPAC_BASES[base]) for base in seq)
        log_redundancy_fold = sum(
            math.log(len(IUPAC_BASES[base])) for base in seq
        )

        if variant_count <= max_variants:
            tms = [
                _calculate_tm_sl_owczarzy_core(
                    "".join(combo),
                    settings,
                    log_redundancy_fold=log_redundancy_fold,
                    is_degenerate=False,
                )
                for combo in product(*(IUPAC_BASES[base] for base in seq))
            ]
            return TmRange(
                min_tm=min(tms),
                mean_tm=sum(tms) / len(tms),
                max_tm=max(tms),
                is_approximate=False,
                variant_count=variant_count,
            )

        mean_tm = _calculate_tm_sl_owczarzy_core(
            seq,
            settings,
            log_redundancy_fold=log_redundancy_fold,
            is_degenerate=True,
        )
        t_ref_k = max(200.0, mean_tm + 273.15)
        weakest_seq = _select_dp_sequence(seq, t_ref_k, pick_max=True)
        strongest_seq = _select_dp_sequence(seq, t_ref_k, pick_max=False)
        weakest_tm = _calculate_tm_sl_owczarzy_core(
            weakest_seq,
            settings,
            log_redundancy_fold=log_redundancy_fold,
            is_degenerate=False,
        )
        strongest_tm = _calculate_tm_sl_owczarzy_core(
            strongest_seq,
            settings,
            log_redundancy_fold=log_redundancy_fold,
            is_degenerate=False,
        )
        return TmRange(
            min_tm=min(weakest_tm, mean_tm),
            mean_tm=mean_tm,
            max_tm=max(strongest_tm, mean_tm),
            is_approximate=True,
            variant_count=variant_count,
        )
    except InsufficientThermodynamicDataError:
        if fallback_to_amplify4:
            tm = calculate_tm_lander_amplify4(primer, settings)
            return TmRange(
                min_tm=tm,
                mean_tm=tm,
                max_tm=tm,
                is_approximate=False,
                variant_count=1,
            )
        raise
