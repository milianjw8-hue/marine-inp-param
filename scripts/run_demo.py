"""Fit the warm-window term and train/validate the ML emulator on an INP spectra CSV.

Input CSV columns: sample_id, T_C, inp_L, toc_g_m3 (+ optional predictors, e.g. log10_chl)
Usage:
  python scripts/run_demo.py --data data/synthetic_spectra.csv --out outputs/demo \
      --features T_C log10_toc log10_chl
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from marinp.emulator import INPEmulator
from marinp.enhancement import WarmWindow, fit_warm_enhancement
from marinp.plots import parity, spectrum_per_toc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="outputs/demo")
    ap.add_argument("--features", nargs="+", default=["T_C", "log10_toc"])
    ap.add_argument("--mode", choices=["residual", "direct"], default="residual")
    ap.add_argument("--t-warm", type=float, default=-10.0)
    ap.add_argument("--t-cold", type=float, default=-22.0)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.data)

    ww = fit_warm_enhancement(df["T_C"], df["inp_L"], df["toc_g_m3"],
                              template=WarmWindow(T_warm=args.t_warm, T_cold=args.t_cold))
    (out / "warm_term_fit.json").write_text(json.dumps(ww.to_dict(), indent=2))

    truth_path = Path(args.data).with_suffix(".truth.json")
    truth = WarmWindow(**json.loads(truth_path.read_text())) if truth_path.exists() else None
    spectrum_per_toc(df, ww, out / "spectrum_per_toc.png", truth=truth)

    emu = INPEmulator(features=args.features, mode=args.mode)
    cv_df, scores = emu.cross_validate(df)
    cv_df.to_csv(out / "cv_predictions.csv", index=False)
    (out / "cv_scores.json").write_text(json.dumps(scores, indent=2))
    parity(cv_df, scores, out / "parity.png")

    print("warm-window fit:", {k: round(v, 3) if isinstance(v, float) else v for k, v in ww.to_dict().items()})
    for name, s in scores.items():
        print(f"{name:9s} log10 bias {s['log10_bias']:+.2f}  RMSE {s['log10_rmse']:.2f}  "
              f"within x2 {100 * s['within_x2']:.0f}%  within x10 {100 * s['within_x10']:.0f}%")


if __name__ == "__main__":
    main()
