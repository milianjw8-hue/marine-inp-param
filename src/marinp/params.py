"""Published marine INP parameterizations.

W15  Wilson, T. W., et al. (2015), Nature, 525, 234-238, doi:10.1038/nature14986
     INPs per gram TOC = exp(11.2186 - 0.4459 * T), T in degC.
     Fitted to sea-surface-microlayer data; used in models between about -27 and -6.5 degC.

M18  McCluskey, C. S., et al. (2018), J. Geophys. Res. Atmos., 123, 6196-6212,
     doi:10.1029/2017JD028033
     ice-nucleation active site density per unit SSA surface area
     n_s(T) = exp(-0.545 * (T - 273.15) + 1.0125), T in K, n_s in m^-2.

Both return NaN outside their nominal temperature range unless `extrapolate=True`.
"""

from __future__ import annotations

import warnings

import numpy as np

W15_A = 11.2186  # ln(g^-1)
W15_B = 0.4459  # degC^-1
W15_T_RANGE = (-27.0, -6.5)  # degC, as applied in global models (e.g. EC-Earth3-AerChem)

M18_SLOPE = 0.545  # K^-1
M18_INTERCEPT = 1.0125
# Nominal range is not fixed here; check McCluskey et al. (2018) before use outside the observed range.
M18_T_RANGE = (-30.0, -5.0)  # degC, permissive default

OM_OC_DEFAULT = 2.0  # organic matter : organic carbon mass ratio (sensitivity range 1.2-2.6)


def _mask(T_C, rng, extrapolate):
    T_C = np.asarray(T_C, dtype=float)
    if extrapolate:
        return T_C, np.ones_like(T_C, dtype=bool)
    ok = (T_C >= rng[0]) & (T_C <= rng[1])
    return T_C, ok


def nm_wilson2015(T_C, extrapolate: bool = False):
    """W15 cumulative INPs per gram of TOC [g^-1]."""
    T_C, ok = _mask(T_C, W15_T_RANGE, extrapolate)
    nm = np.exp(W15_A - W15_B * T_C)
    return np.where(ok, nm, np.nan)


def inp_wilson2015(T_C, toc_g_m3, extrapolate: bool = False):
    """W15 INP concentration [L^-1] from TOC mass concentration [g C m^-3]."""
    return nm_wilson2015(T_C, extrapolate) * np.asarray(toc_g_m3, dtype=float) * 1e-3


def ns_mccluskey2018(T_C, extrapolate: bool = False):
    """M18 ice-nucleation active site density [m^-2] of sea spray aerosol."""
    T_C, ok = _mask(T_C, M18_T_RANGE, extrapolate)
    T_K = T_C + 273.15
    ns = np.exp(-M18_SLOPE * (T_K - 273.15) + M18_INTERCEPT)
    return np.where(ok, ns, np.nan)


def inp_mccluskey2018(T_C, ssa_surface_m2_m3, extrapolate: bool = False):
    """M18 INP concentration [L^-1] from SSA surface area concentration [m^2 m^-3]."""
    return ns_mccluskey2018(T_C, extrapolate) * np.asarray(ssa_surface_m2_m3, dtype=float) * 1e-3


def moa_to_toc(moa_g_m3, om_oc: float = OM_OC_DEFAULT):
    """Marine organic aerosol mass -> TOC mass, TOC = MOA / (OM:OC)."""
    if not 1.0 <= om_oc <= 3.0:
        warnings.warn(f"OM:OC = {om_oc} is outside the usual 1.2-2.6 range", stacklevel=2)
    return np.asarray(moa_g_m3, dtype=float) / om_oc


def ug_to_g(x_ug_m3):
    """Convenience: ug m^-3 -> g m^-3."""
    return np.asarray(x_ug_m3, dtype=float) * 1e-6
