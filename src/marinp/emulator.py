"""Physics-constrained machine-learning emulator of marine INP spectra.

The model predicts log10 INP [L^-1] from temperature and aerosol/ocean predictors with a
gradient-boosted tree ensemble whose response is constrained to be
  * monotonically non-increasing in temperature (colder -> more INP), and
  * monotonically non-decreasing in log10 TOC (more organic carbon -> more INP).

Two target modes:
  "direct"   learn log10 INP directly
  "residual" learn log10(INP_obs / INP_W15); prediction = W15 * 10**residual (hybrid physics-ML)

Validation is grouped by spectrum (sample) so that points from the same filter never appear
in both training and test folds.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

from .params import inp_wilson2015


def log_scores(obs, pred) -> dict:
    """Scores in log10 space (INP is approximately log-normal)."""
    obs, pred = np.asarray(obs, float), np.asarray(pred, float)
    ok = np.isfinite(obs) & np.isfinite(pred) & (obs > 0) & (pred > 0)
    d = np.log10(pred[ok]) - np.log10(obs[ok])
    if d.size == 0:
        return {"n": 0}
    return {
        "n": int(d.size),
        "log10_bias": float(d.mean()),
        "log10_rmse": float(np.sqrt((d**2).mean())),
        "within_x2": float((np.abs(d) <= np.log10(2)).mean()),
        "within_x10": float((np.abs(d) <= 1).mean()),
        "r_log10": float(np.corrcoef(np.log10(obs[ok]), np.log10(pred[ok]))[0, 1]) if d.size > 2 else np.nan,
    }


@dataclass
class INPEmulator:
    features: list[str] = field(default_factory=lambda: ["T_C", "log10_toc"])
    mode: str = "residual"
    max_iter: int = 300
    learning_rate: float = 0.05
    max_leaf_nodes: int = 15
    min_samples_leaf: int = 10
    l2_regularization: float = 1.0
    random_state: int = 0

    def _constraints(self):
        cst = []
        for f in self.features:
            if f == "T_C":
                cst.append(-1)
            elif f == "log10_toc":
                cst.append(1 if self.mode == "direct" else 0)
            else:
                cst.append(0)
        return cst

    def _model(self):
        return HistGradientBoostingRegressor(
            max_iter=self.max_iter,
            learning_rate=self.learning_rate,
            max_leaf_nodes=self.max_leaf_nodes,
            min_samples_leaf=self.min_samples_leaf,
            l2_regularization=self.l2_regularization,
            monotonic_cst=self._constraints(),
            random_state=self.random_state,
        )

    @staticmethod
    def prepare(df: pd.DataFrame) -> pd.DataFrame:
        """Add derived columns: log10_toc and log10_inp, W15 baseline."""
        out = df.copy()
        out["log10_toc"] = np.log10(out["toc_g_m3"])
        out["inp_w15"] = inp_wilson2015(out["T_C"], out["toc_g_m3"], extrapolate=True)
        if "inp_L" in out:
            out["log10_inp"] = np.log10(out["inp_L"])
        return out

    def _target(self, df):
        if self.mode == "direct":
            return df["log10_inp"].to_numpy()
        return (df["log10_inp"] - np.log10(df["inp_w15"])).to_numpy()

    def fit(self, df: pd.DataFrame):
        df = self.prepare(df)
        df = df[np.isfinite(df["log10_inp"])]
        self.model_ = self._model().fit(df[self.features].to_numpy(), self._target(df))
        return self

    def predict(self, df: pd.DataFrame) -> np.ndarray:
        df = self.prepare(df)
        y = self.model_.predict(df[self.features].to_numpy())
        if self.mode == "direct":
            return 10.0**y
        return df["inp_w15"].to_numpy() * 10.0**y

    def cross_validate(self, df: pd.DataFrame, group_col: str = "sample_id", n_splits: int = 5) -> tuple[pd.DataFrame, dict]:
        """Grouped K-fold CV. Returns (predictions table, scores dict incl. W15 baseline)."""
        df = self.prepare(df)
        df = df[np.isfinite(df["log10_inp"])].reset_index(drop=True)
        pred = np.full(len(df), np.nan)
        for tr, te in GroupKFold(n_splits=n_splits).split(df, groups=df[group_col]):
            m = INPEmulator(**{k: getattr(self, k) for k in self.__dataclass_fields__})
            m.fit(df.iloc[tr])
            pred[te] = m.predict(df.iloc[te])
        df["inp_pred"] = pred
        scores = {"emulator": log_scores(df["inp_L"], pred), "W15": log_scores(df["inp_L"], df["inp_w15"])}
        return df, scores
