"""marinp: marine biogenic ice-nucleating particle (INP) parameterizations and ML emulators.

Units used throughout
---------------------
Temperature          T_C  in degC (negative for supercooled conditions)
TOC mass             g C m^-3 of air
SSA surface area     m^2 m^-3 of air
INP concentration    L^-1 of air (standard air conditions are NOT applied)
"""

__version__ = "0.2.0"

from .params import inp_wilson2015, inp_mccluskey2018, nm_wilson2015, ns_mccluskey2018, moa_to_toc  # noqa: F401
from .enhancement import WarmWindow, fit_warm_enhancement, inp_enhanced  # noqa: F401
from .spectra import frozen_fraction_to_inp  # noqa: F401
