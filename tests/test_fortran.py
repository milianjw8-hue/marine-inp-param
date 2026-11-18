"""The Fortran module (fortran/marine_inp_mod.F90) must give the same numbers as the Python package."""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from marinp.enhancement import WarmWindow, nm_enhanced
from marinp.params import nm_wilson2015, ns_mccluskey2018

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not installed")


@pytest.fixture(scope="module")
def table(tmp_path_factory):
    d = tmp_path_factory.mktemp("f90")
    exe = d / "driver"
    subprocess.run(["gfortran", "-O0", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-J", str(d),
                    str(ROOT / "fortran" / "marine_inp_mod.F90"), str(ROOT / "fortran" / "test_driver.F90"),
                    "-o", str(exe)], check=True)
    ww = WarmWindow(a=6.5, b=-0.12, T_warm=-10.0, T_cold=-22.0, T_ref=-16.0, delta=1.5)
    out = subprocess.run([str(exe)], input=f"{ww.a} {ww.b} {ww.T_warm} {ww.T_cold} {ww.T_ref} {ww.delta}\n",
                         capture_output=True, text=True, check=True).stdout
    return np.loadtxt(out.splitlines()), ww


def test_w15_and_m18_match(table):
    t, ww = table
    T = t[:, 0]
    py_w15 = np.nan_to_num(nm_wilson2015(T), nan=0.0)   # Fortran returns 0 outside the range
    py_m18 = np.nan_to_num(ns_mccluskey2018(T), nan=0.0)
    assert np.allclose(t[:, 1], py_w15, rtol=1e-12, atol=0)
    assert np.allclose(t[:, 2], py_m18, rtol=1e-12, atol=0)
    assert (t[(T < -27) | (T > -6.5), 1] == 0).all()


def test_warm_window_matches(table):
    t, ww = table
    T = t[:, 0]
    assert np.allclose(t[:, 3], ww.window(T), rtol=1e-12, atol=1e-300)
    assert np.allclose(t[:, 4], ww.nm_add(T), rtol=1e-12, atol=1e-300)
    py = np.nan_to_num(nm_wilson2015(T), nan=0.0) + ww.nm_add(T)
    assert np.allclose(t[:, 5], py, rtol=1e-12)
    # the extrapolating Python form agrees where W15 is in range
    inside = (T >= -27) & (T <= -6.5)
    assert np.allclose(t[inside, 5], nm_enhanced(T[inside], ww), rtol=1e-12)


def test_immersion_rate(table):
    t, ww = table
    T = t[:, 0]
    n_inp = (np.nan_to_num(nm_wilson2015(T), nan=0.0) + ww.nm_add(T)) * 1e-7      # m-3
    expect = np.maximum(np.minimum(n_inp, 5e7) - 1e3, 0) / 1800.0
    assert np.allclose(t[:, 6], expect, rtol=1e-12, atol=1e-300)
    assert (t[:, 6] >= 0).all()
