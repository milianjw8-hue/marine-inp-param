import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from marinp.emulator import INPEmulator, log_scores

ROOT = Path(__file__).resolve().parents[1]


def _data(tmp_path):
    out = tmp_path / "spectra.csv"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "make_synthetic_data.py"), "--out", str(out),
                    "--n-samples", "40"], check=True)
    return pd.read_csv(out)


def test_log_scores_perfect():
    s = log_scores([1, 10, 100], [1, 10, 100])
    assert s["log10_rmse"] == 0 and s["within_x2"] == 1


def test_emulator_beats_w15_in_grouped_cv(tmp_path):
    df = _data(tmp_path)
    emu = INPEmulator(features=["T_C", "log10_toc", "log10_chl"], max_iter=150)
    cv, scores = emu.cross_validate(df, n_splits=4)
    assert scores["emulator"]["log10_rmse"] < scores["W15"]["log10_rmse"]
    assert np.isfinite(cv["inp_pred"]).all()


def test_monotonic_in_temperature(tmp_path):
    df = _data(tmp_path)
    emu = INPEmulator(features=["T_C", "log10_toc"], mode="direct", max_iter=150).fit(df)
    grid = pd.DataFrame({"T_C": np.arange(-26, -6, 1.0), "toc_g_m3": 2e-7})
    pred = emu.predict(grid)
    assert np.all(np.diff(pred) <= 1e-12)  # warmer -> not more INP
