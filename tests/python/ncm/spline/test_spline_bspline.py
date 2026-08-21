"""Tests for NcmSplineBSpline."""

import math
import os
import subprocess
import sys

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


def _run_child(code: str) -> subprocess.CompletedProcess:
    """Run `code` in a child process, for errors that abort rather than raise.

    NCM_FFTW_* is dropped from the environment: a test that sets an invalid
    planner value during the session would otherwise abort the child at init,
    before it reaches the code under test.
    """
    env = {k: v for k, v in os.environ.items() if not k.startswith("NCM_FFTW")}
    return subprocess.run(
        [sys.executable, "-c", code],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


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


@pytest.mark.parametrize(
    "n,tol,expected_order",
    [
        (100, 1.0e-06, 4),
        (100, 1.0e-08, 6),
        (100, 1.0e-10, 8),
        (500, 1.0e-06, 4),
        (500, 1.0e-10, 6),
        (2000, 1.0e-06, 4),
        (2000, 1.0e-14, 6),
    ],
)
def test_order_derived_from_tolerance(n, tol, expected_order):
    """The order is an output: the lowest one that meets the request, not a default."""
    xs = np.linspace(0.0, 2400.0, n)
    sp = Ncm.SplineBSpline.new_tol(tol, 0.0)
    sp.set_array(xs.tolist(), _gauss(xs).tolist(), True)

    assert sp.get_order() == expected_order


@pytest.mark.parametrize("n,tol", [(100, 1.0e-06), (500, 1.0e-10), (2000, 1.0e-14)])
def test_achieved_error_is_honest(n, tol):
    """The reported estimate matches the true interpolation error, and meets the request.

    Callers use this to predict their own accuracy: an integral of the spline inherits
    roughly this error times its cancellation ratio.
    """
    xs = np.linspace(0.0, 2400.0, n)
    sp = Ncm.SplineBSpline.new_tol(tol, 0.0)
    sp.set_array(xs.tolist(), _gauss(xs).tolist(), True)

    true_err = _max_midpoint_error(xs, sp)
    assert true_err <= tol
    assert sp.get_achieved_error() == pytest.approx(true_err, rel=0.5)


def test_impossible_request_fails_loudly():
    """A request the samples cannot support must abort, never silently return worse.

    100 samples top out near 2.7e-11; asking for 1e-13 has no valid answer, and quietly
    returning a 100x worse spline would hide the one fact the caller needs.
    """
    xs = np.linspace(0.0, 2400.0, 100)
    ys = _gauss(xs)
    code = (
        "import numpy as np\n"
        "from numcosmo_py import Ncm\n"
        "Ncm.cfg_init()\n"
        f"xs = np.linspace(0.0, 2400.0, {len(xs)})\n"
        "ys = np.exp(-0.5*((xs-1200.0)/300.0)**2)\n"
        "sp = Ncm.SplineBSpline.new_tol(1.0e-13, 0.0)\n"
        "sp.set_array(xs.tolist(), ys.tolist(), True)\n"
    )
    proc = _run_child(code)

    assert proc.returncode != 0
    assert "cannot support a requested interpolation error" in proc.stderr


def test_manual_order_is_not_overridden():
    """reltol = 0 (the default) leaves the explicitly chosen order alone."""
    xs = np.linspace(0.0, 2400.0, 2000)
    sp = Ncm.SplineBSpline.new_full(
        10,
        Ncm.Vector.new_array(xs.tolist()),
        Ncm.Vector.new_array(_gauss(xs).tolist()),
        True,
    )

    assert sp.get_order() == 10
    assert sp.get_achieved_error() == 0.0


def test_copy_empty_preserves_the_order():
    """copy_empty carries the order over.

    The order is the only configuration an empty spline has, so a copy that
    silently reverted to the default would interpolate at a different accuracy
    than the spline it was copied from.
    """
    _, sp = _build(6, 200)
    empty = sp.copy_empty()

    assert isinstance(empty, Ncm.SplineBSpline)
    assert empty.get_order() == 6
    assert empty.get_len() == 0

    # A distinct object: filling the copy leaves the original untouched.
    xs = np.linspace(0.0, 2400.0, 200)
    empty.set_array(xs.tolist(), np.cos(xs / 400.0).tolist(), True)

    assert empty.eval(700.0) == pytest.approx(math.cos(700.0 / 400.0), rel=1.0e-9)
    assert sp.eval(700.0) == pytest.approx(_gauss(700.0), rel=1.0e-9)


@pytest.mark.parametrize("order,tol", [(4, 1.0e-8), (6, 1.0e-4)])
def test_deriv_nmax_matches_the_closed_form(order, tol):
    """The highest non-zero derivative is the degree: for x^(order-1) it is (order-1)!.

    The tolerance is per order because this is the worst-conditioned quantity the
    spline reports -- measured 3e-10 at order 4, 4e-6 at order 6 and 6e-2 at
    order 8, on 200 samples.
    """
    deg = order - 1
    xs = np.linspace(0.0, 2.0, 200)
    sp = Ncm.SplineBSpline.new_full(
        order,
        Ncm.Vector.new_array(xs.tolist()),
        Ncm.Vector.new_array((xs**deg).tolist()),
        True,
    )
    exact = float(math.factorial(deg))

    for x in (0.3, 0.9, 1.5):
        assert sp.eval_deriv_nmax(x) == pytest.approx(exact, rel=tol)


@pytest.mark.parametrize("order,accessor", [(2, "eval_deriv"), (3, "eval_deriv2")])
def test_deriv_nmax_agrees_with_the_dedicated_accessor(order, accessor):
    """At order 2 and 3 the top derivative is the 1st and 2nd, which have their own
    accessors -- the same number by an independent route, so it must match exactly.
    """
    _, sp = _build(order, 200)

    for x in (300.0, 900.0, 1500.0):
        assert sp.eval_deriv_nmax(x) == getattr(sp, accessor)(x)


def test_name_reports_the_order():
    """The name carries the order, and is rebuilt when the order changes."""
    # The vfunc takes the instance explicitly through introspection, hence do_name(sp).
    sp = Ncm.SplineBSpline.new(4)
    assert sp.do_name(sp) == "NcmSplineBSpline[order 4]"

    sp.set_order(8)
    assert sp.do_name(sp) == "NcmSplineBSpline[order 8]"


def test_name_identifies_the_spline_in_the_min_size_error():
    """What the name is for: naming the offending spline in ncm_spline_set's error.

    Order 8 needs 8 samples, so 5 is rejected -- and the message has to say which
    spline refused, since a program can hold several.
    """
    code = (
        "import numpy as np\n"
        "from numcosmo_py import Ncm\n"
        "Ncm.cfg_init()\n"
        "sp = Ncm.SplineBSpline.new(8)\n"
        "xs = np.linspace(0.0, 1.0, 5)\n"
        "sp.set_array(xs.tolist(), (xs**2).tolist(), True)\n"
    )
    proc = _run_child(code)

    assert proc.returncode != 0
    assert "min size for [NcmSplineBSpline[order 8]] is 8" in proc.stderr
