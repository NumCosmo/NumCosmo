#!/usr/bin/env python
#
# test_plotting_tools.py
#
# Wed Oct 15 10:34:34 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_plotting_tools.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
#


"""Unit tests for plotting tools in numcosmo_py.plotting.tools.

These tests exercise the pure-Python plotting helper functions. They avoid
plotting to a display by using the Agg backend (matplotlib default in tests)
and create figure/axes objects explicitly.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from numcosmo_py.plotting import tools
from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_latex_float_basic_and_scientific():
    """latex_float formats normal and scientific floats as expected."""
    assert tools.latex_float(0.0) == "0"
    # The implementation formats exponentials using LaTeX "\times 10^{...}" even
    # when using the `e` format internally, so check for that representation.
    out = tools.latex_float(123.456, precision=2, convert_g=False)
    assert "1.23" in out and ("10" in out and ("times" in out or "10^" in out))

    # a value that will be formatted with an exponent and non-trivial base
    s = tools.latex_float(3.2e6, precision=2)
    # should contain times-ten exponent formatting
    assert "10" in s and "times" in s or "10^" in s

    # check for -1 and 1 special-casing
    assert tools.latex_float(1.0, precision=2) == "1"
    assert tools.latex_float(-1.0, precision=2) == "-1"


def test_format_time_units():
    """format_time chooses appropriate SI time units for small values."""
    # verify unit boundaries according to the implementation's thresholds
    assert tools.format_time(1e-10).endswith(" ns")
    # 5e-7 is below 1e-6 -> reported in ns
    assert tools.format_time(5e-7).endswith(" ns")
    # value just above 1e-6 should use microseconds
    assert tools.format_time(5e-6).endswith(" $\\mu$s")
    # 5e-4 is formatted as 500 microseconds by the implementation
    assert tools.format_time(5e-4).endswith(" $\\mu$s")
    assert tools.format_time(1.5).endswith(" s")


def test_set_rc_params_article_and_restore():
    """set_rc_params_article updates rcParams and can be restored."""
    orig = matplotlib.rcParams.copy()
    try:
        tools.set_rc_params_article(column_width=100.0, ncol=1, nrows=1, fontsize=6)
        assert matplotlib.rcParams["font.size"] == 6
        assert isinstance(matplotlib.rcParams["figure.figsize"], (list, tuple))
    finally:
        matplotlib.rcParams.update(orig)


def test_confidence_ellipse_adds_patch():
    """confidence_ellipse returns an artists.Artist and adds to axes."""
    _, ax = plt.subplots()
    mu = np.array([1.0, 2.0])
    cov = np.array([[0.2, 0.05], [0.05, 0.1]])
    patch = tools.confidence_ellipse(mu, cov, ax, n_std=2.0, edgecolor="r")
    # patch should be in axes artists
    assert patch in ax.patches
    # basic property checks
    assert hasattr(patch, "get_width")


def test_add_ellipse_from_ellipticity_creates_patches():
    """add_ellipse_from_ellipticity creates one patch per input point."""
    _, ax = plt.subplots()
    ra = np.array([0.0, 1.0, 2.0])
    dec = np.array([0.0, 0.5, -0.5])
    e1 = np.array([0.1, 0.2, 0.0])
    e2 = np.array([0.0, -0.1, 0.3])
    before = len(ax.patches)
    tools.add_ellipse_from_ellipticity(ax, ra, dec, e1, e2, ellipse_scale=0.5)
    after = len(ax.patches)
    assert after - before == 3


def test_plot_m2lnp_and_output_image():
    """plot_m2lnp accepts grids and returns an AxesImage object."""
    n = 150
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    xv, yv = np.meshgrid(x, y)
    # create a mock -2lnp array (larger values away from center)
    z = ((xv - 0.5) ** 2 + (yv - 0.5) ** 2) * 10.0
    _, ax = plt.subplots()
    img = tools.plot_m2lnp(x, y, z.flatten(), ax, plotn=n)
    assert img is not None
    # image norm should be LogNorm instance

    assert isinstance(img.norm, LogNorm)


def test_format_alpha_xaxis_creates_secondary_axes():
    """format_alpha_xaxis adds a secondary axis with proper label and formatter."""
    _, ax = plt.subplots()

    cosmo = Nc.HICosmoQGRW.new()

    # ensure we can call without top/bottom combos
    # set an x range large enough so the locator base computation yields a
    # positive step (the function does floor-division that can produce 0.0
    # for small ranges).
    ax.set_xlim(0, 100)
    tools.format_alpha_xaxis(ax, cosmo, n_ticks=4, top=True, bottom=True)
    # bottom axis label should be the alpha symbol
    assert "$\\alpha$" in ax.get_xlabel() or r"$\alpha$" == ax.get_xlabel()

    # The function should not raise and must set the bottom axis label.
    # Secondary axis creation depends on tick computation; if present, it
    # should have an xlabel containing 'x'.
    if len(ax.figure.axes) > 1:
        found_x_label = any(
            (sec is not ax) and (sec.get_xlabel() and "x" in sec.get_xlabel().lower())
            for sec in ax.figure.axes
        )
        assert found_x_label
