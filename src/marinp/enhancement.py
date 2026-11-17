"""Warm-temperature biogenic enhancement added to the W15 baseline.

Observed marine INP spectra often exceed W15 in a warm window (about -10 to -22 degC),
which is commonly attributed to heat-labile, biological INPs. This module represents that
excess as an ADDITIVE term in the same per-gram-TOC space as W15:

    n_total(T) = n_W15(T) + n_add(T)
    n_add(T)   = 10**(a + b * (T - T_ref)) * W(T)

W(T) is a smooth window that is ~1 between T_cold and T_warm and decays to 0 outside,
built from two logistic ramps of width `delta` (degC). a and b are fitted to data.
No default values of a and b are implied to be physical; they come from your observations.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import least_squares

from .params import nm_wilson2015


@dataclass
class WarmWindow:
    a: float = 7.0  # log10 g^-1 at T_ref (placeholder until fitted)
    b: float = -0.15  # log10 g^-1 per degC (negative: more INP when colder)
    T_warm: float = -10.0
    T_cold: float = -22.0
    T_ref: float = -16.0
    delta: float = 1.0
    a_se: float = float("nan")
    b_se: float = float("nan")
    n_fit: int = 0

    def window(self, T_C):
        T_C = np.asarray(T_C, dtype=float)
        lo = 1.0 / (1.0 + np.exp(-(T_C - self.T_cold) / self.delta))
        hi = 1.0 / (1.0 + np.exp((T_C - self.T_warm) / self.delta))
        return lo * hi

    def nm_add(self, T_C):
        """Additive INP per gram TOC [g^-1]."""
        T_C = np.asarray(T_C, dtype=float)
        return 10.0 ** (self.a + self.b * (T_C - self.T_ref)) * self.window(T_C)

    def to_dict(self):
        return asdict(self)


def nm_enhanced(T_C, ww: WarmWindow, extrapolate: bool = True):
    return nm_wilson2015(T_C, extrapolate=extrapolate) + ww.nm_add(T_C)


def inp_enhanced(T_C, toc_g_m3, ww: WarmWindow, extrapolate: bool = True):
    """W15 + warm-window term, INP concentration [L^-1]."""
    return nm_enhanced(T_C, ww, extrapolate) * np.asarray(toc_g_m3, dtype=float) * 1e-3


def fit_warm_enhancement(T_C, inp_obs_L, toc_g_m3, template: WarmWindow | None = None) -> WarmWindow:
    """Fit a and b of the warm-window term to observed INP spectra.

    The cost is the log10 misfit of the full model (W15 + term) against observations,
    using only points inside [T_cold - 2*delta, T_warm + 2*delta] with positive INP and TOC.
    Standard errors come from the Gauss-Newton covariance at the optimum.
    """
    ww = WarmWindow(**(template.to_dict() if template else {}))
    T = np.asarray(T_C, dtype=float)
    y = np.asarray(inp_obs_L, dtype=float)
    toc = np.broadcast_to(np.asarray(toc_g_m3, dtype=float), T.shape)
    sel = (
        np.isfinite(T) & np.isfinite(y) & np.isfinite(toc) & (y > 0) & (toc > 0)
        & (T >= ww.T_cold - 2 * ww.delta) & (T <= ww.T_warm + 2 * ww.delta)
    )
    if sel.sum() < 3:
        raise ValueError("need at least 3 positive observations inside the warm window")
    T, y, toc = T[sel], y[sel], toc[sel]
    nm_obs = y * 1e3 / toc  # g^-1
    base = nm_wilson2015(T, extrapolate=True)

    def resid(p):
        ww.a, ww.b = p
        return np.log10(base + ww.nm_add(T)) - np.log10(nm_obs)

    # start from the median excess over W15
    excess = np.clip(nm_obs - base, 1e-30, None)
    a0 = float(np.median(np.log10(excess)))
    res = least_squares(resid, x0=[a0, ww.b], method="lm")
    ww.a, ww.b = float(res.x[0]), float(res.x[1])
    dof = max(sel.sum() - 2, 1)
    s2 = (res.fun**2).sum() / dof
    try:
        cov = np.linalg.inv(res.jac.T @ res.jac) * s2
        ww.a_se, ww.b_se = (float(x) for x in np.sqrt(np.diag(cov)))
    except np.linalg.LinAlgError:
        pass
    ww.n_fit = int(sel.sum())
    return ww
