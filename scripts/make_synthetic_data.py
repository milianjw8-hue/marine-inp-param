"""Generate SYNTHETIC marine INP spectra for testing (not real observations).

Each sample (filter) has a TOC mass concentration, a chlorophyll-like bloom indicator and a
sample-level biological factor. INP spectra follow W15 + a known warm-window term scaled by
the biological factor, with measurement noise and a detection limit.

Output CSV columns: sample_id, T_C, inp_L, toc_g_m3, log10_chl
Usage: python scripts/make_synthetic_data.py --out data/synthetic_spectra.csv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from marinp.enhancement import WarmWindow, nm_enhanced

TRUTH = WarmWindow(a=8.5, b=-0.12)  # synthetic "true" warm term, used only to make test data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/synthetic_spectra.csv")
    ap.add_argument("--n-samples", type=int, default=80)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    temps = np.arange(-6.0, -26.5, -1.0)
    rows = []
    for s in range(args.n_samples):
        log10_toc = rng.normal(-6.7, 0.35)  # g C m^-3 (about 0.2 ug C m^-3)
        bio = rng.normal(0.0, 0.4)  # log10 factor on the warm term
        log10_chl = 0.5 + 1.5 * bio + rng.normal(0, 0.2)  # bloom indicator correlated with bio
        ww = WarmWindow(a=TRUTH.a + bio, b=TRUTH.b)
        nm = nm_enhanced(temps, ww)
        inp = nm * 10**log10_toc * 1e-3 * 10 ** rng.normal(0, 0.15, temps.size)
        for T, y in zip(temps, inp):
            if y >= 1e-4:  # detection limit, L^-1
                rows.append({"sample_id": f"S{s:03d}", "T_C": T, "inp_L": y,
                             "toc_g_m3": 10**log10_toc, "log10_chl": log10_chl})
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    out.with_suffix(".truth.json").write_text(json.dumps(TRUTH.to_dict(), indent=2))
    print(f"{len(rows)} points from {args.n_samples} synthetic samples -> {out}")


if __name__ == "__main__":
    main()
