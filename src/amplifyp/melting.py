# Copyright (C) 2026 Fufu Fang
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
DNA primers using the Nearest-Neighbor thermodynamics model.
"""

import math
from typing import Final

from .dna import Primer
from .settings import (
    DEFAULT_AMPLIFY4_TM_SETTINGS,
    Amplify4TMSettings,
    TMSettings,
)

# Thermodynamic parameters (SantaLucia 1998)
# Enthalpy (dH) in cal/mol, Entropy (dS) in cal/(K*mol)
# Keys are dinucleotides (5'->3')
# Source: SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and
# oligonucleotide DNA
# nearest-neighbor thermodynamics". PNAS, 95(4), 1460-1465.
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


def calculate_tm(primer: Primer, settings: TMSettings) -> float:
    """Calculate the melting temperature (Tm) of a primer sequence.

    Uses the Nearest-Neighbor model with SantaLucia 1998 thermodynamic
    parameters and salt corrections.

    Formula used:
        Tm = (dH / (dS + R * ln(Ct/4))) - 273.15

    Where:
        dH: Total enthalpy (cal/mol)
        dS: Total entropy (cal/(K*mol))
        R: Gas constant (1.987 cal/(K*mol))
        Ct: Total oligonucleotide concentration (M)

    Salt correction is applied to the entropy term (dS) as per SantaLucia 1998:
        dS_corrected = dS_1M + 0.368 * (N - 1) * ln([Na+])

    Args:
        primer: The primer object containing the sequence.
        settings: Melting settings containing concentrations.

    Returns:
        The melting temperature in degrees Celsius.
    """
    seq = primer.seq.upper()
    n = len(seq)
    if n < 2:
        return 0.0

    # Constants
    R: Final[float] = 1.987  # cal/(K*mol)

    # Calculate enthalpy and entropy
    dh: float = 0.0
    ds: float = 0.0

    # Initiation with terminal corrections (SantaLucia 1998)
    # Init w/ term G.C: dH=0.1 kcal, dS=-2.8 eu
    # Init w/ term A.T: dH=2.3 kcal, dS=4.1 eu
    # Values converted to cal: 100, 2300, etc.
    # Note: paper assumes unit kcal/mol. Dictionary above is in cal/mol.

    # Terminal corrections
    left_term = seq[0]
    right_term = seq[-1]

    # Left end
    if left_term in "GC":
        dh += 100
        ds += -2.8
    elif left_term in "AT":
        dh += 2300
        ds += 4.1

    # Right end
    if right_term in "GC":
        dh += 100
        ds += -2.8
    elif right_term in "AT":
        dh += 2300
        ds += 4.1

    # Nearest neighbor steps
    for i in range(n - 1):
        dinuc = seq[i : i + 2]
        if dinuc in NN_THERMO_DATA:
            val = NN_THERMO_DATA[dinuc]
            dh += val[0]
            ds += val[1]
        else:
            # Skip invalid dinucleotides (e.g. Ns)
            pass

    # Salt Correction (Owczarzy 2008)
    # References:
    # Owczarzy et al. (2008). "Predicting Stability of DNA Duplexes in
    # Solutions Containing Magnesium and Monovalent Cations". Biochemistry,
    # 47(19), 5336-5353.

    # 1. Calculate concentrations
    # Monovalent: Na+, K+, Tris+ (mM -> M)
    mono_mM = settings.monovalent_salt_conc
    mono_M = mono_mM * 1e-3

    # Divalent: Mg2+ (mM -> M)
    # We ignore dNTPs for now (which chelate Mg2+) as per standard simple
    # inputs,
    # or assume free Mg2+ is provided.
    div_mM = settings.divalent_salt_conc
    div_M = div_mM * 1e-3

    # 2. Determine mode (Monovalent only, Mixed, or Divalent dominant)
    # Ratio R = sqrt([Mg2+]) / [Mon+]
    if mono_M == 0:
        ratio = 999.0  # Large number, Divalent dominant
    else:
        ratio = math.sqrt(div_M) / mono_M

    # Base params for entropy correction (SantaLucia 1998 for Monovalent)
    # dS_corr = dS_1M + 0.368 * (N-1) * ln([Na+])

    # We need to compute the 1M Na+ background Tm first for Owczarzy methods
    # But current code adds corrections to dS directly.
    # SantaLucia 1998 works by correcting dS.
    # Owczarzy 2008 works by correcting 1/Tm.

    # Let's calculate the "base" Tm (at 1M Na+) first
    # 1M Na+ means NO entropy correction term (ln(1) = 0)

    # Calculate Tm_1M_Na (Kelvin)
    total_dna_conc_M = settings.dna_conc * 1e-9
    if total_dna_conc_M <= 0:
        total_dna_conc_M = 50e-9

    denom_1M = ds + R * math.log(total_dna_conc_M / 4.0)
    if denom_1M == 0:
        return 0.0

    tm_1m_K = dh / denom_1M

    # Now apply corrections
    if div_M == 0:
        # Monovalent only (SantaLucia 1998)
        # Apply strict SantaLucia 1998 entropy correction
        # This effectively modifies the denom.
        if mono_M > 0:
            ds_corr = 0.368 * (n - 1) * math.log(mono_M)
            denom_corr = ds + ds_corr + R * math.log(total_dna_conc_M / 4.0)
            tm_final_K = dh / denom_corr
        else:
            tm_final_K = tm_1m_K

    elif ratio < 0.22:
        # Monovalent dominant (use Monovalent logic)
        # Owczarzy recommends using the Monovalent formula with
        # effective salt? Or just strict Monovalent?
        # Paper says "Monovalent ion dominant", implies treating roughly as Na+.
        # We will use the SantaLucia correction with [Mon+]
        ds_corr = 0.368 * (n - 1) * math.log(mono_M)
        denom_corr = ds + ds_corr + R * math.log(total_dna_conc_M / 4.0)
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

        # Calculate fraction of GC
        fgc = (seq.count("G") + seq.count("C")) / n

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
        # well. Note also Von Ahsen et al (2001) suggests Na_eq = Mon +
        # 3.8*sqrt(Mg) for mixed. But we use the Owczarzy Eq 16 which is
        # specific for Mg.

        tm_inv = (1.0 / tm_1m_K) + corr
        tm_final_K = 1.0 / tm_inv

    return tm_final_K - 273.15


def calculate_tm_amplify4(
    primer: Primer,
    amplify4_settings: Amplify4TMSettings = DEFAULT_AMPLIFY4_TM_SETTINGS,
) -> float:
    """Calculate Tm using the original Amplify4 algorithm.

    This method is a direct port of the ``calcTm`` method from the Swift
    codebase (Primer.swift). It uses its own set of entropy and enthalpy
    tables (typically 5x5 matrices) and specific correction factors.

    Args:
        primer: The primer object containing the sequence.
        amplify4_settings: Amplify4-specific settings (tables, etc.).

    Returns:
        The melting temperature in degrees Celsius. Returns 0.0 if the
        sequence length is less than 1.
    """
    seq = primer.seq.upper()
    seq_len = len(seq)
    if seq_len < 1:
        return 0.0

    # Effective length limit
    n = min(amplify4_settings.effective_primer_length, seq_len)

    # Map bases to indices: G=0, A=1, T=2, C=3, Other=4
    # (Matches Amplify4: Primer.swift)
    seqn: list[int] = []
    for c in seq:
        if c == "G":
            seqn.append(0)
        elif c == "A":
            seqn.append(1)
        elif c == "T":
            seqn.append(2)
        elif c == "C":
            seqn.append(3)
        else:
            seqn.append(4)

    entropy = amplify4_settings.entropy
    enthalpy = amplify4_settings.enthalpy

    entr: float = 108.0
    enth: float = 0.0

    # Sum neighbors
    # Note: entropy/enthalpy tables in Swift are accessed as [y][x]
    # where x is current base, y is next base.
    for i in range(n - 1):
        x = seqn[i]
        y = seqn[i + 1]
        entr += entropy[y][x]
        enth += enthalpy[y][x]

    # Scaling applied in Swift code
    entr = -entr * 0.1
    enth = -enth * 0.1

    # Corrections
    # DNAConc in settings is usually nM (default 50).
    # Swift formula: 1.987 * log(DNAConc/4e9)
    # If DNAConc=50, 50/4e9 = 1.25e-8 => 12.5 nM concentration assumption?
    # Or maybe it treats input as raw number?
    # We use settings.DNAConc exactly as Swift does.
    dna_conc_val = amplify4_settings.dna_conc
    log_dna = 1.987 * math.log(dna_conc_val / 4.0e9)

    # Salt: 16.6 * log(saltConc/1000) / log(10.0)
    # saltConc is typically mM (default 50).
    salt_conc_val = amplify4_settings.monovalent_salt_conc
    if salt_conc_val <= 0:
        salt_conc_val = 50.0  # prevent log error if 0

    log_salt = 16.6 * math.log10(salt_conc_val / 1000.0)

    tm = (enth * 1000.0) / (entr + log_dna) - 273.15 + log_salt
    return tm
