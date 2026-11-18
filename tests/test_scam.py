"""SCAM6 comparison on fake history files, and a compile check of the CAM interface with stub modules."""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

xr = pytest.importorskip("xarray")
from marinp.scam import compare, find_cases  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

EFFECT = {"none": 0.0, "m18": 0.1, "w15": 0.2, "w15warm": 0.5}   # fake glaciation strength


def _fake_case(d: Path, scheme: str, hm: str):
    nt, nl = 24, 10
    lev = np.linspace(500, 1000, nl)
    hyai, hybi = np.zeros(nl + 1), np.linspace(0.45, 1.0, nl + 1)
    hyam, hybm = 0.5 * (hyai[1:] + hyai[:-1]), 0.5 * (hybi[1:] + hybi[:-1])
    e = EFFECT[scheme] + (0.0 if hm == "on" else -0.05)
    ni = np.full((nt, nl, 1), 1e3 * (1 + 4 * e))          # kg-1
    t = np.full((nt, nl, 1), 260.0)
    mar = np.full((nt, nl, 1), 0.0 if scheme == "none" else 1e-3 * (1 + e))
    ds = xr.Dataset({"T": (("time", "lev", "ncol"), t), "NUMICE": (("time", "lev", "ncol"), ni),
                     "TGCLDLWP": (("time", "ncol"), np.full((nt, 1), 0.15 * (1 - e))),
                     "TGCLDIWP": (("time", "ncol"), np.full((nt, 1), 0.02 * (1 + 2 * e))),
                     "MARFRZIMM": (("time", "lev", "ncol"), mar),
                     "DSTFREZIMM": (("time", "lev", "ncol"), np.full((nt, nl, 1), 1e-3)),
                     "BCFREZIMM": (("time", "lev", "ncol"), np.full((nt, nl, 1), 1e-4)),
                     "PS": (("time", "ncol"), np.full((nt, 1), 101000.0)), "P0": 100000.0,
                     "hyam": ("lev", hyam), "hybm": ("lev", hybm), "hyai": ("ilev", hyai), "hybi": ("ilev", hybi)},
                    coords={"time": pd.date_range("2004-10-05", periods=nt, freq="h"), "lev": lev,
                            "ilev": np.linspace(450, 1000, nl + 1)})
    (d / "run").mkdir(parents=True)
    ds.to_netcdf(d / "run" / f"{d.name}.cam.h0.2004-10-05-07171.nc")


def test_compare_separates_inp_and_hm(tmp_path):
    for s in EFFECT:
        for hm in ("on", "off"):
            _fake_case(tmp_path / f"mpace_{s}_hm{hm}", s, hm)
    cases = find_cases(tmp_path)
    assert len(cases) == 8
    t = compare(cases).set_index("case")
    assert t.loc["mpace_w15warm_hmon", "d_lwp_gm2_vs_none"] < t.loc["mpace_w15_hmon", "d_lwp_gm2_vs_none"] < 0
    assert t.loc["mpace_none_hmon", "d_lwp_gm2_vs_none"] == 0
    assert t.loc["mpace_w15_hmoff", "d_ni_max_L_hm_off_minus_on"] < 0
    assert np.isnan(t.loc["mpace_w15_hmon", "d_ni_max_L_hm_off_minus_on"])
    assert t.loc["mpace_none_hmon", "marine_share"] == 0
    assert 0 < t.loc["mpace_w15warm_hmon", "marine_share"] < 1
    assert t["liquid_fraction"].between(0, 1).all()


@pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not installed")
def test_cam_interface_compiles_with_stubs(tmp_path):
    (tmp_path / "marinp_in").write_text(
        "&marinp_nl\n marinp_scheme='w15warm', moa_sfc=0.2, moa_h=800., ww_a=6.5, ww_b=-0.12, disable_hm=.true.\n/\n")
    exe = tmp_path / "drv"
    subprocess.run(["gfortran", "-fcheck=all", "-J", str(tmp_path), str(ROOT / "tests" / "fortran_stubs" / "stubs.F90"),
                    str(ROOT / "fortran" / "marine_inp_mod.F90"),
                    str(ROOT / "cam6" / "SourceMods" / "src.cam" / "marine_inp_cam.F90"),
                    str(ROOT / "tests" / "fortran_stubs" / "drv.F90"), "-o", str(exe)], check=True)
    out = subprocess.run([str(exe)], cwd=tmp_path, capture_output=True, text=True, check=True).stdout
    vals = out.split("=")[-1].split()
    ninp = float(vals[1])
    from marinp.enhancement import WarmWindow
    from marinp.params import nm_wilson2015

    ww = WarmWindow(a=6.5, b=-0.12)
    toc = 0.2 * np.exp(-200 / 800) * 1e-6 / 2.0
    assert ninp == pytest.approx((nm_wilson2015(-15.0) + ww.nm_add(-15.0)) * toc, rel=1e-10)
    assert vals[2] == "T"


def test_patches_reference_interface():
    het = (ROOT / "cam6" / "patches" / "hetfrz_classnuc_cam.F90.patch").read_text()
    mg = (ROOT / "cam6" / "patches" / "micro_mg2_0.F90.patch").read_text()
    assert "call marinp_frzimm(" in het and "marinp_cam_readnl()" in het and "'MARFRZIMM'" in het
    assert "marinp_disable_hm" in mg
