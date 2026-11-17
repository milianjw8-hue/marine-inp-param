import numpy as np
import pytest

from marinp.enhancement import WarmWindow, fit_warm_enhancement, inp_enhanced
from marinp.params import inp_wilson2015
from marinp.spectra import detection_limit, frozen_fraction_to_inp


def test_window_shape():
    ww = WarmWindow()
    assert ww.window(-16.0) == pytest.approx(1.0, abs=0.01)
    assert ww.window(-5.0) < 0.01
    assert ww.window(-28.0) < 0.01


def test_enhanced_never_below_w15():
    T = np.arange(-26, -7, 0.5)
    ww = WarmWindow(a=8.0, b=-0.1)
    assert np.all(inp_enhanced(T, 1e-7, ww) >= inp_wilson2015(T, 1e-7, extrapolate=True))


def test_fit_recovers_known_parameters():
    rng = np.random.default_rng(3)
    truth = WarmWindow(a=8.3, b=-0.10)
    T = np.repeat(np.arange(-24.0, -8.0, 1.0), 20)
    toc = 10 ** rng.normal(-6.7, 0.3, T.size)
    y = inp_enhanced(T, toc, truth) * 10 ** rng.normal(0, 0.1, T.size)
    fit = fit_warm_enhancement(T, y, toc)
    assert fit.a == pytest.approx(8.3, abs=0.05)
    assert fit.b == pytest.approx(-0.10, abs=0.01)
    assert fit.a_se < 0.05


def test_fit_needs_data():
    with pytest.raises(ValueError):
        fit_warm_enhancement([-5.0, -30.0], [1.0, 1.0], 1e-7)


def test_vali_conversion():
    # 50 % frozen, 1 uL drops, 10 mL wash, 1000 L air -> ln2 / 1e-3 mL * 10 mL / 1000 L
    inp = frozen_fraction_to_inp(0.5, v_drop_uL=1.0, v_liquid_mL=10.0, v_air_L=1000.0)
    assert inp == pytest.approx(np.log(2) / 1e-3 * 10 / 1000)
    assert np.isnan(frozen_fraction_to_inp(1.0, 1.0, 10.0, 1000.0))
    assert detection_limit(100, 1.0, 10.0, 1000.0) > 0
