"""Compare SCAM6 M-PACE experiments (cam6/scam6_mpace_matrix.sh).

Needs `pip install "marinp[scam]"` (xarray, netCDF4).

Case names follow mpace_<scheme>_hm<on|off>. For each case the time mean (optionally over a window)
of the column and the mixed-phase layer is reported:

  lwp_gm2, iwp_gm2        liquid and ice water path (TGCLDLWP, TGCLDIWP), g m-2
  liquid_fraction         LWP / (LWP + IWP)
  ni_max_L                maximum of the time-mean in-cloud-ish ice number profile, L-1 (NUMICE * rho)
  marfrz_col, dstfrz_col  column-integrated marine and dust immersion freezing rates, m-2 s-1
  marine_share            marine / (marine + dust + BC) immersion freezing

Then each scheme is compared with the control (scheme 'none') at the same HM setting, and each HM-off
case with its HM-on twin, to separate the INP effect from the rime-splintering effect.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

RAIR = 287.04
G = 9.80616
NAME = re.compile(r"mpace_(?P<scheme>[a-z0-9]+)_hm(?P<hm>on|off)")


def open_case(path):
    import xarray as xr

    p = Path(path)
    files = sorted(p.glob("*.cam.h0.*.nc")) if p.is_dir() else [p]
    if not files:
        raise FileNotFoundError(f"no *.cam.h0.*.nc in {p}")
    return xr.open_mfdataset(files, combine="by_coords") if len(files) > 1 else xr.open_dataset(files[0])


def _column_integral(rate, ds):
    """Integrate a (time, lev, ...) rate in m-3 s-1 over height using hydrostatic layer thickness."""
    if "PDELDRY" in ds:
        dp = ds["PDELDRY"]
    elif "hyai" in ds and "PS" in ds:
        pint = ds["hyai"] * ds["P0"] + ds["hybi"] * ds["PS"]
        dp = pint.diff("ilev").rename(ilev="lev").assign_coords(lev=ds["lev"])
    else:
        return np.nan
    rho = _pressure(ds) / (RAIR * ds["T"])
    dz = dp / (rho * G)
    return (rate * dz).sum("lev")


def _pressure(ds):
    if "hyam" in ds and "PS" in ds:
        return ds["hyam"] * ds["P0"] + ds["hybm"] * ds["PS"]
    return ds["lev"] * 100.0


def summarize(ds, start=None, end=None) -> dict:
    if start or end:
        ds = ds.sel(time=slice(start, end))
    out = {}
    lwp = float(ds["TGCLDLWP"].mean()) * 1e3 if "TGCLDLWP" in ds else np.nan
    iwp = float(ds["TGCLDIWP"].mean()) * 1e3 if "TGCLDIWP" in ds else np.nan
    out["lwp_gm2"], out["iwp_gm2"] = lwp, iwp
    out["liquid_fraction"] = lwp / (lwp + iwp) if (lwp + iwp) > 0 else np.nan
    if "NUMICE" in ds and "T" in ds:
        rho = _pressure(ds) / (RAIR * ds["T"])
        ni = (ds["NUMICE"] * rho * 1e-3).mean("time")       # # kg-1 * kg m-3 -> m-3 -> L-1
        out["ni_max_L"] = float(ni.max())
    rates = {}
    for key, var in (("marfrz_col", "MARFRZIMM"), ("dstfrz_col", "DSTFREZIMM"), ("bcfrz_col", "BCFREZIMM")):
        if var in ds:
            rates[key] = float(_column_integral(ds[var], ds).mean())
            out[key] = rates[key]
    tot = sum(v for v in rates.values() if np.isfinite(v))
    out["marine_share"] = rates.get("marfrz_col", np.nan) / tot if tot > 0 else np.nan
    return out


def compare(cases: dict, start=None, end=None) -> pd.DataFrame:
    """cases: {case_name: path}. Returns one row per case plus differences to the matching control."""
    rows = []
    for name, path in cases.items():
        m = NAME.search(name)
        if not m:
            raise ValueError(f"case name {name} does not match mpace_<scheme>_hm<on|off>")
        ds = open_case(path)
        rows.append({"case": name, **m.groupdict(), **summarize(ds, start, end)})
    t = pd.DataFrame(rows)
    metrics = ["lwp_gm2", "iwp_gm2", "liquid_fraction", "ni_max_L"]
    ctrl = t[t.scheme == "none"].set_index("hm")
    for c in metrics:
        if c in t:
            t[f"d_{c}_vs_none"] = [r[c] - ctrl[c].get(r["hm"], np.nan) for _, r in t.iterrows()]
    on = t[t.hm == "on"].set_index("scheme")
    for c in metrics:
        if c in t:
            t[f"d_{c}_hm_off_minus_on"] = [r[c] - on[c].get(r["scheme"], np.nan) if r["hm"] == "off" else np.nan
                                          for _, r in t.iterrows()]
    return t.sort_values(["hm", "scheme"]).reset_index(drop=True)


def find_cases(root) -> dict:
    """Map case names to their run or archive folder under root (first match wins)."""
    out = {}
    for p in sorted(Path(root).glob("mpace_*_hm*")):
        for sub in (p / "run", p / "atm" / "hist", p):
            if sub.is_dir() and list(sub.glob("*.cam.h0.*.nc")):
                out[p.name] = sub
                break
    return out
