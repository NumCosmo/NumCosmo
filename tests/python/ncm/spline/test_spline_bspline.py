"""Tests for NcmSplineBSpline."""

import math

import numpy as np
import pytest

from numcosmo_py import Ncm

Ncm.cfg_init()


def _gauss(x):
    return np.exp(-0.5 * ((x - 1200.0) / 300.0) ** 2)


def _build(order, n, f=_gauss, a=0.0, b=2400.0):
    xs = np.linspace(a, b, n)
    ys = f(xs)
    return xs, Ncm.SplineBSpline.new_full(
        order,
        Ncm.Vector.new_array(xs.tolist()),
        Ncm.Vector.new_array(ys.tolist()),
        True,
    )


def _max_midpoint_error(xs, sp, f=_gauss):
    mid = 0.5 * (xs[:-1] + xs[1:])
    return max(abs(sp.eval(float(m)) - f(m)) for m in mid)


def test_order_property():
    """Order round-trips and defaults to 8."""
    assert Ncm.SplineBSpline.new(8).get_order() == 8
    sbs = Ncm.SplineBSpline.new(4)
    assert sbs.get_order() == 4
    sbs.set_order(6)
    assert sbs.get_order() == 6


def test_interpolates_the_samples():
    """A knot value is reproduced exactly, whatever the order."""
    for order in (4, 6, 8):
        xs, sp = _build(order, 200)
        for i in (0, 37, 100, 199):
            assert sp.eval(float(xs[i])) == pytest.approx(_gauss(xs[i]), abs=1e-13)


@pytest.mark.parametrize(
    "order,n,tol",
    [
        (4, 100, 1.0e-06),
        (6, 100, 1.0e-08),
        (8, 100, 1.0e-10),
        (4, 1000, 1.0e-10),
        (6, 1000, 1.0e-14),
        (8, 1000, 1.0e-14),
    ],
)
def test_accuracy_improves_with_order(order, n, tol):
    """Interpolation error at the measured level; see the class docs table."""
    xs, sp = _build(order, n)
    assert _max_midpoint_error(xs, sp) < tol


def test_cubic_cannot_reach_machine_precision():
    """The reason this class exists: order 4 stalls where order 8 does not.

    A cubic spline is C^2 and its error floor sits far above machine precision at any
    sample density, so data supplied as a fixed table caps the accuracy of everything
    downstream of it.
    """
    xs4, sp4 = _build(4, 2000)
    xs8, sp8 = _build(8, 2000)
    err4 = _max_midpoint_error(xs4, sp4)
    err8 = _max_midpoint_error(xs8, sp8)

    assert err4 > 1.0e-13
    assert err8 < 1.0e-14
    assert err4 / err8 > 100.0


def test_derivatives_against_closed_form():
    xs, sp = _build(8, 1000)
    x0 = 1350.0
    g = _gauss(x0)
    d1 = -(x0 - 1200.0) / 300.0**2 * g
    d2 = ((x0 - 1200.0) ** 2 / 300.0**4 - 1.0 / 300.0**2) * g

    assert sp.eval_deriv(x0) == pytest.approx(d1, rel=1.0e-11)
    assert sp.eval_deriv2(x0) == pytest.approx(d2, rel=1.0e-9)


def test_integral_against_closed_form():
    """Gaussian over a finite range: erf, not the infinite-range normalisation.

    The endpoints sit 4 sigma out, so the two truncated tails are worth 4.8e-2 -- far
    above the tolerance here, which is why the closed form must be the finite one.
    """
    xs, sp = _build(8, 1000)
    sigma, mu = 300.0, 1200.0
    exact = (
        sigma
        * math.sqrt(math.pi / 2.0)
        * (
            math.erf((2400.0 - mu) / (sigma * math.sqrt(2.0)))
            - math.erf((0.0 - mu) / (sigma * math.sqrt(2.0)))
        )
    )

    assert sp.eval_integ(0.0, 2400.0) == pytest.approx(exact, rel=1.0e-12)


def test_reprepare_after_order_change():
    """Changing the order rebuilds the workspace rather than reusing a stale one."""
    xs = np.linspace(0.0, 2400.0, 500)
    xv = Ncm.Vector.new_array(xs.tolist())
    yv = Ncm.Vector.new_array(_gauss(xs).tolist())
    sp = Ncm.SplineBSpline.new_full(4, xv, yv, True)
    err_cubic = _max_midpoint_error(xs, sp)

    sp.set_order(8)
    sp.prepare()
    err_oct = _max_midpoint_error(xs, sp)

    assert err_oct < err_cubic / 100.0


def test_serialization_round_trip():
    """The class is registered, so it survives a serialize/deserialize cycle.

    prepare() is required after deserialisation: NcmSpline carries the samples but not
    the fitted coefficients, and this is the framework-wide convention, not something
    specific to this class -- NcmSplineCubicNotaknot segfaults if you skip it.
    """
    xs, sp = _build(6, 200)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    dup = ser.dup_obj(sp)
    dup.prepare()

    assert isinstance(dup, Ncm.SplineBSpline)
    assert dup.get_order() == 6
    assert dup.eval(1234.5) == pytest.approx(sp.eval(1234.5), rel=1.0e-14)
