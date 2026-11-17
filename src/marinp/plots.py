"""Standard figures: INP-per-gram-TOC spectrum and emulator parity plot."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .enhancement import WarmWindow, nm_enhanced  # noqa: E402
from .params import W15_T_RANGE, nm_wilson2015  # noqa: E402


def spectrum_per_toc(df, ww: WarmWindow, path, truth: WarmWindow | None = None):
    """Observed INP per gram TOC vs temperature with W15 and W15 + warm-window fit."""
    T = np.linspace(-28, -5, 231)
    nm_obs = df["inp_L"] * 1e3 / df["toc_g_m3"]
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    ax.scatter(df["T_C"], nm_obs, s=8, c="0.6", alpha=0.6, label="observations")
    ax.plot(T, nm_wilson2015(T, extrapolate=True), "k--", lw=1.2, label="W15 (extrapolated outside range)")
    Tin = T[(T >= W15_T_RANGE[0]) & (T <= W15_T_RANGE[1])]
    ax.plot(Tin, nm_wilson2015(Tin), "k-", lw=1.6, label="W15")
    ax.plot(T, nm_enhanced(T, ww), "-", color="#d62728", lw=2,
            label=f"W15 + warm term (fit: a={ww.a:.2f}, b={ww.b:.3f})")
    if truth is not None:
        ax.plot(T, nm_enhanced(T, truth), ":", color="#1f77b4", lw=1.5, label="synthetic truth")
    ax.axvspan(ww.T_cold, ww.T_warm, color="#d62728", alpha=0.07)
    ax.set_yscale("log")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("INP per gram TOC (g$^{-1}$)")
    ax.invert_xaxis()
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def parity(cv_df, scores, path):
    """Predicted vs observed INP for W15 and the emulator (grouped CV predictions)."""
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2), sharex=True, sharey=True)
    lo = np.nanmin(cv_df["inp_L"]) / 3
    hi = np.nanmax(cv_df["inp_L"]) * 3
    for ax, col, name in zip(axes, ["inp_w15", "inp_pred"], ["W15", "emulator"]):
        sc = ax.scatter(cv_df["inp_L"], cv_df[col], s=6, alpha=0.5, c=cv_df["T_C"], cmap="viridis")
        x = np.array([lo, hi])
        ax.plot(x, x, "k-", lw=0.8)
        ax.plot(x, x * 10, "k:", lw=0.8)
        ax.plot(x, x / 10, "k:", lw=0.8)
        s = scores[name]
        ax.set_title(f"{name}: RMSE(log10)={s['log10_rmse']:.2f}, ×10 {100 * s['within_x10']:.0f}%", fontsize=9)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel("observed INP (L$^{-1}$)")
        ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel("predicted INP (L$^{-1}$)")
    fig.tight_layout()
    fig.colorbar(sc, ax=axes, label="Temperature (°C)", shrink=0.85)
    fig.savefig(path, dpi=140)
    plt.close(fig)
