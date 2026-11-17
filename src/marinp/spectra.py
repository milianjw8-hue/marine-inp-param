"""Convert droplet-freezing assay results to atmospheric INP concentrations.

Uses the singular (time-independent) description of Vali (1971):
    K(T) = -ln(1 - f(T)) / V_drop             [per volume of liquid]
    INP(T) = K(T) * dilution * V_liquid / V_air [per volume of air]

f is the cumulative frozen fraction of droplets at temperature T.
"""

from __future__ import annotations

import numpy as np


def frozen_fraction_to_inp(frozen_fraction, v_drop_uL: float, v_liquid_mL: float, v_air_L: float, dilution: float = 1.0):
    """INP concentration [L^-1 of air] from cumulative frozen fraction.

    Parameters
    ----------
    frozen_fraction : array, 0-1 cumulative frozen fraction at each temperature
    v_drop_uL       : droplet volume [uL]
    v_liquid_mL     : volume of the suspension the filter was washed into [mL]
    v_air_L         : volume of air sampled through the filter [L]
    dilution        : dilution factor of the suspension relative to the wash (>= 1)

    Fully frozen points (f = 1) are undefined and returned as NaN (lower limit only).
    """
    f = np.asarray(frozen_fraction, dtype=float)
    f = np.where((f >= 0) & (f < 1), f, np.nan)
    k_per_mL = -np.log1p(-f) / (v_drop_uL * 1e-3)
    return k_per_mL * dilution * v_liquid_mL / v_air_L


def detection_limit(n_drops: int, v_drop_uL: float, v_liquid_mL: float, v_air_L: float, dilution: float = 1.0):
    """INP concentration [L^-1] corresponding to one frozen droplet."""
    return frozen_fraction_to_inp(1.0 / n_drops, v_drop_uL, v_liquid_mL, v_air_L, dilution)
