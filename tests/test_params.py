import numpy as np
import pytest

from marinp.params import inp_mccluskey2018, inp_wilson2015, moa_to_toc, nm_wilson2015, ns_mccluskey2018


def test_w15_published_form():
    # exp(11.2186 - 0.4459*T), T in degC
    assert nm_wilson2015(-15.0) == pytest.approx(np.exp(11.2186 + 0.4459 * 15.0))


def test_w15_units_to_per_litre():
    toc = 2e-7  # g C m^-3
    expected = np.exp(11.2186 + 0.4459 * 20.0) * toc / 1e3
    assert inp_wilson2015(-20.0, toc) == pytest.approx(expected)


def test_w15_range_masking():
    assert np.isnan(nm_wilson2015(-30.0))
    assert np.isnan(nm_wilson2015(-3.0))
    assert np.isfinite(nm_wilson2015(-30.0, extrapolate=True))


def test_m18_published_form():
    # exp(-0.545*(T_K-273.15) + 1.0125)
    assert ns_mccluskey2018(-20.0) == pytest.approx(np.exp(0.545 * 20.0 + 1.0125))
    assert inp_mccluskey2018(-20.0, 1e-5) == pytest.approx(np.exp(0.545 * 20.0 + 1.0125) * 1e-5 / 1e3)


def test_monotonic_in_temperature():
    T = np.arange(-26, -7, 1.0)
    assert np.all(np.diff(nm_wilson2015(T)) < 0)


def test_moa_to_toc():
    assert moa_to_toc(2.0) == pytest.approx(1.0)
    with pytest.warns(UserWarning):
        moa_to_toc(1.0, om_oc=5.0)
