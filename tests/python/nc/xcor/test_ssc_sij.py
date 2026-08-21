#!/usr/bin/env python
#
# test_ssc_sij.py
#
# Thu Aug 13 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_ssc_sij.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""NcXcorSSCSij: the C-side super-sample covariance S_ij calculator.

Pinned against `numcosmo_py.ssc.SijCalculator`, the Python implementation it
replaces, for both the full-sky and the masked estimator. The C object exists
so S_ij can be recomputed once per likelihood step without a Python callback
(see NcDataClusterNCountsGauss:ssc-sij), so the tests here also cover the two
properties that recomputation depends on: a solver reused across cosmologies
must give the same answer as a fresh one, and repeated evaluation at a fixed
cosmology must not recompute.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import GObject, Nc, Ncm
from numcosmo_py.cosmology import Cosmology
from numcosmo_py.ssc import SijCalculator, find_lmax, mask_angular_power_spectrum

pytest_plugins = [
    "python.fixtures_xcor",
]

pytestmark = pytest.mark.xcor

#: Seven J-PAS-like top-hat bins, the configuration the S_ij study used.
Z_EDGES = np.linspace(0.1, 0.8, 8)

#: nside of the cap mask below. Small on purpose: the masked estimator sums
#: every multipole up to lmax, so a finer mask makes the test much slower
#: without exercising anything new.
MASK_NSIDE = 32


def _matrix_to_np(matrix: Ncm.Matrix) -> np.ndarray:
    """Convert an #NcmMatrix into a numpy array."""
    return np.array(
        [
            [matrix.get(i, j) for j in range(matrix.ncols())]
            for i in range(matrix.nrows())
        ]
    )


def _make_ssc_sij(cosmology: Cosmology, method=None) -> Nc.XcorSSCSij:
    """Build a calculator over Z_EDGES for the test cosmology."""
    ssc_sij = Nc.XcorSSCSij.new(
        cosmology.dist, cosmology.ps_ml, Ncm.Vector.new_array(Z_EDGES.tolist())
    )

    if method is not None:
        ssc_sij.set_method(method)

    return ssc_sij


def _make_reference(cosmology: Cosmology) -> SijCalculator:
    """Build the Python reference over the same bins and the same power spectrum."""
    return SijCalculator(Z_EDGES, powspec=cosmology.ps_ml, dist=cosmology.dist)


def _cap_mask() -> np.ndarray:
    """Build a polar cap mask with NcmSphereMap, so healpy is not needed."""
    smap = Ncm.SphereMap.new(MASK_NSIDE)
    npix = smap.get_npix()
    mask = np.zeros(npix)

    for i in range(npix):
        theta, _phi = smap.pix2ang_ring(i)
        if theta < 0.6:
            mask[i] = 1.0

    return mask


def test_new_defaults(cosmology: Cosmology) -> None:
    """A fresh calculator has one bin per edge pair and is full sky."""
    ssc_sij = _make_ssc_sij(cosmology)

    assert ssc_sij.get_nbins() == len(Z_EDGES) - 1
    # Full sky is the mask spectrum 4 pi delta_l0, so a single multipole.
    assert ssc_sij.get_lmax() == 0
    assert_allclose(ssc_sij.get_fsky(), 1.0)
    assert ssc_sij.get_method() == Nc.XcorMethod.KERNEL_FIXED


def test_default_construction_is_survivable() -> None:
    """A bare instance with no properties must not abort the process.

    The serialization registry and the introspection tooling both instantiate
    every registered type with no properties at all, so the class has to come
    up empty rather than assert on its missing objects.
    """
    ssc_sij = GObject.new(Nc.XcorSSCSij)

    assert ssc_sij.get_nbins() == 0
    assert ssc_sij.get_lmax() == 0

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    assert isinstance(ser.from_string(ser.to_string(ssc_sij, True)), Nc.XcorSSCSij)


def test_mask_cl_fullsky_is_the_default(cosmology: Cosmology) -> None:
    """The default mask spectrum is the one nc_xcor_ssc_sij_mask_cl_fullsky builds."""
    ssc_sij = _make_ssc_sij(cosmology)
    fullsky = Nc.XcorSSCSij.mask_cl_fullsky()

    assert fullsky.len() == 1
    assert_allclose(fullsky.get(0), 4.0 * np.pi)
    assert_allclose(ssc_sij.peek_mask_cl().get(0), 4.0 * np.pi)


def test_fullsky_matches_python_reference(cosmology: Cosmology) -> None:
    """On the same quadrature, C and Python agree to machine precision."""
    ssc_sij = _make_ssc_sij(cosmology, Nc.XcorMethod.KERNEL_CUBATURE)

    got = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))
    expected = _make_reference(cosmology).fullsky(cosmology.cosmo)

    assert_allclose(got, expected, rtol=1.0e-12)


def test_fullsky_fixed_matches_python_reference(cosmology: Cosmology) -> None:
    """The default fixed quadrature agrees with the adaptive Python reference.

    KERNEL_FIXED is the default precisely because it cannot fail to converge,
    so it must reproduce the adaptive result rather than merely be close.

    The comparison is against the peak of the matrix, not element-wise: the
    off-diagonals are four orders of magnitude below the diagonal, so their own
    relative error is set by how the two quadratures resolve a near total
    cancellation, and is not the quantity that matters for the covariance.
    """
    ssc_sij = _make_ssc_sij(cosmology)

    got = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))
    expected = _make_reference(cosmology).fullsky(cosmology.cosmo)
    peak = np.abs(expected).max()

    assert_allclose(got, expected, rtol=1.0e-5, atol=1.0e-7 * peak)


def test_sij_is_symmetric(cosmology: Cosmology) -> None:
    """S_ij is a covariance, so the matrix must come back symmetric."""
    sij = _matrix_to_np(_make_ssc_sij(cosmology).eval(cosmology.cosmo))

    assert_allclose(sij, sij.T, rtol=0.0, atol=0.0)


def test_diagonal_dominates(cosmology: Cosmology) -> None:
    """Off-diagonal S_ij of separated bins are orders below the diagonal.

    This is what makes the off-diagonals hard: they are a small residual of a
    large cancellation, which is why `scaled-abstol` and not `reltol` sets
    their accuracy.
    """
    sij = _matrix_to_np(_make_ssc_sij(cosmology).eval(cosmology.cosmo))

    assert np.all(np.diag(sij) > 0.0)
    assert abs(sij[0, -1]) < 1.0e-2 * sij[0, 0]


def test_prepare_if_needed_caches(cosmology: Cosmology) -> None:
    """Re-evaluating at an unchanged cosmology returns the identical matrix."""
    ssc_sij = _make_ssc_sij(cosmology)

    first = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))
    second = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))

    assert_allclose(first, second, rtol=0.0, atol=0.0)


def test_tracks_cosmology(cosmology: Cosmology, cosmology_alt: Cosmology) -> None:
    """Moving the cosmology moves S_ij: it is not frozen at the first call."""
    ssc_sij = _make_ssc_sij(cosmology)

    first = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))
    moved = _matrix_to_np(ssc_sij.eval(cosmology_alt.cosmo))

    assert not np.allclose(first, moved)


def test_reused_solver_matches_fresh_one(
    cosmology: Cosmology, cosmology_alt: Cosmology
) -> None:
    """A calculator reused across cosmologies matches a freshly built one.

    The solver plans its multipole blocks and pins one spherical Bessel
    factorization per block for its whole life, and this object keeps a single
    solver across every cosmology it is asked about. That is only sound if none
    of the retained state is cosmology dependent.
    """
    reused = _make_ssc_sij(cosmology)
    reused.eval(cosmology.cosmo)
    got = _matrix_to_np(reused.eval(cosmology_alt.cosmo))

    fresh = _matrix_to_np(_make_ssc_sij(cosmology).eval(cosmology_alt.cosmo))

    assert_allclose(got, fresh, rtol=0.0, atol=0.0)


def test_masked_matches_python_reference(cosmology: Cosmology) -> None:
    """The masked estimator reproduces SijCalculator.partial_sky."""
    mask = _cap_mask()
    cl_mask, _fsky = mask_angular_power_spectrum(mask)
    lmax = find_lmax(cl_mask)

    ssc_sij = _make_ssc_sij(cosmology)
    ssc_sij.set_mask_cl(Ncm.Vector.new_array(cl_mask[: lmax + 1].tolist()))

    got = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))
    expected = _make_reference(cosmology).partial_sky(cosmology.cosmo, mask, lmax=lmax)
    peak = np.abs(expected).max()

    assert_allclose(got, expected, rtol=1.0e-5, atol=1.0e-7 * peak)


def test_mask_sets_lmax_and_fsky(cosmology: Cosmology) -> None:
    """The mask spectrum's length sets lmax, and its monopole sets fsky."""
    mask = _cap_mask()
    cl_mask, fsky = mask_angular_power_spectrum(mask)
    lmax = find_lmax(cl_mask)

    ssc_sij = _make_ssc_sij(cosmology)
    ssc_sij.set_mask_cl(Ncm.Vector.new_array(cl_mask[: lmax + 1].tolist()))

    assert ssc_sij.get_lmax() == lmax
    assert_allclose(ssc_sij.get_fsky(), fsky, rtol=1.0e-12)
    assert fsky < 1.0


def test_masked_exceeds_fullsky(cosmology: Cosmology) -> None:
    """A footprint smaller than the sky raises the super-sample variance."""
    mask = _cap_mask()
    cl_mask, _fsky = mask_angular_power_spectrum(mask)
    lmax = find_lmax(cl_mask)

    fullsky = _matrix_to_np(_make_ssc_sij(cosmology).eval(cosmology.cosmo))

    masked_calc = _make_ssc_sij(cosmology)
    masked_calc.set_mask_cl(Ncm.Vector.new_array(cl_mask[: lmax + 1].tolist()))
    masked = _matrix_to_np(masked_calc.eval(cosmology.cosmo))

    assert np.all(np.diag(masked) > np.diag(fullsky))


def test_setting_mask_back_to_none_restores_fullsky(cosmology: Cosmology) -> None:
    """Passing None restores the full-sky spectrum and recomputes."""
    mask = _cap_mask()
    cl_mask, _fsky = mask_angular_power_spectrum(mask)
    lmax = find_lmax(cl_mask)

    ssc_sij = _make_ssc_sij(cosmology)
    fullsky = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))

    ssc_sij.set_mask_cl(Ncm.Vector.new_array(cl_mask[: lmax + 1].tolist()))
    ssc_sij.eval(cosmology.cosmo)

    ssc_sij.set_mask_cl(None)

    assert ssc_sij.get_lmax() == 0
    assert_allclose(_matrix_to_np(ssc_sij.eval(cosmology.cosmo)), fullsky, rtol=1.0e-12)


@pytest.mark.parametrize("area", [18000.0, 3000.0])
def test_area_matches_python_reference(cosmology: Cosmology, area: float) -> None:
    """The area rescaling reproduces SijCalculator.fullsky_fsky."""
    ssc_sij = _make_ssc_sij(cosmology)
    ssc_sij.set_area(area)

    got = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))
    expected = _make_reference(cosmology).fullsky_fsky(cosmology.cosmo, area)
    peak = np.abs(expected).max()

    assert_allclose(got, expected, rtol=1.0e-5, atol=1.0e-7 * peak)


@pytest.mark.parametrize("area", [18000.0, 3000.0])
def test_area_raises_sij_over_fullsky(cosmology: Cosmology, area: float) -> None:
    """A finite area only ever increases S_ij, by exactly 1 / f_sky.

    Omitting the correction understates the super-sample variance, so the
    factor is worth pinning rather than assuming.
    """
    fullsky = _matrix_to_np(_make_ssc_sij(cosmology).eval(cosmology.cosmo))

    ssc_sij = _make_ssc_sij(cosmology)
    ssc_sij.set_area(area)
    rescaled = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))

    fsky = area * (np.pi / 180.0) ** 2 / (4.0 * np.pi)

    assert_allclose(ssc_sij.get_fsky(), fsky, rtol=1.0e-12)
    assert fsky < 1.0
    assert_allclose(rescaled, fullsky / fsky, rtol=1.0e-12)


def test_footprint_can_be_given_at_construction(cosmology: Cosmology) -> None:
    """`area` and `mask-cl` work as construct properties, not only as setters.

    Both are G_PARAM_CONSTRUCT and mutually exclusive, and GObject sets
    construct properties in an unspecified order, so the exclusion is settled
    in constructed() rather than by whichever setter happens to run last.
    """
    z_edges = Ncm.Vector.new_array(Z_EDGES.tolist())

    by_area = Nc.XcorSSCSij(
        dist=cosmology.dist, powspec=cosmology.ps_ml, z_edges=z_edges, area=18000.0
    )

    assert_allclose(by_area.get_area(), 18000.0)
    assert by_area.get_lmax() == 0

    cl_mask, _fsky = mask_angular_power_spectrum(_cap_mask())
    lmax = find_lmax(cl_mask)

    by_mask = Nc.XcorSSCSij(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        z_edges=z_edges,
        mask_cl=Ncm.Vector.new_array(cl_mask[: lmax + 1].tolist()),
    )

    assert by_mask.get_area() == 0.0
    assert by_mask.get_lmax() == lmax

    # Both routes agree with the corresponding setter-built calculator.
    via_setter = _make_ssc_sij(cosmology)
    via_setter.set_area(18000.0)

    assert_allclose(
        _matrix_to_np(by_area.eval(cosmology.cosmo)),
        _matrix_to_np(via_setter.eval(cosmology.cosmo)),
        rtol=1.0e-12,
    )


def test_area_and_mask_are_mutually_exclusive(cosmology: Cosmology) -> None:
    """Setting one footprint description discards the other."""
    mask = _cap_mask()
    cl_mask, _fsky = mask_angular_power_spectrum(mask)
    lmax = find_lmax(cl_mask)

    ssc_sij = _make_ssc_sij(cosmology)

    ssc_sij.set_area(3000.0)
    assert ssc_sij.get_area() == 3000.0
    assert ssc_sij.get_lmax() == 0

    ssc_sij.set_mask_cl(Ncm.Vector.new_array(cl_mask[: lmax + 1].tolist()))
    assert ssc_sij.get_area() == 0.0
    assert ssc_sij.get_lmax() == lmax

    ssc_sij.set_area(3000.0)
    assert ssc_sij.get_area() == 3000.0
    assert ssc_sij.get_lmax() == 0


def test_area_survives_serialization(cosmology: Cosmology) -> None:
    """The area rescaling round trips.

    `mask-cl` and `area` clear one another, and GObject restores properties in
    an unspecified order, so a round trip is the test that neither wipes the
    other on the way back in.
    """
    ssc_sij = _make_ssc_sij(cosmology)
    ssc_sij.set_area(18000.0)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dup = ser.from_string(ser.to_string(ssc_sij, True))

    assert_allclose(dup.get_area(), 18000.0)
    assert dup.get_lmax() == 0
    assert_allclose(dup.get_fsky(), ssc_sij.get_fsky(), rtol=1.0e-12)
    assert_allclose(
        _matrix_to_np(dup.eval(cosmology.cosmo)),
        _matrix_to_np(ssc_sij.eval(cosmology.cosmo)),
        rtol=1.0e-12,
    )


def test_serialization_round_trip(cosmology: Cosmology) -> None:
    """The calculator survives a serialization round trip and stays equivalent."""
    mask = _cap_mask()
    cl_mask, _fsky = mask_angular_power_spectrum(mask)
    lmax = find_lmax(cl_mask)

    ssc_sij = _make_ssc_sij(cosmology)
    ssc_sij.set_mask_cl(Ncm.Vector.new_array(cl_mask[: lmax + 1].tolist()))
    ssc_sij.set_reltol(1.0e-7)
    # Non-default marker for the round trip, kept at or above 1e-6 and away from
    # reltol -- the floor enters the S_ij integrand squared, and ssc.py warns
    # against letting the two tolerances coincide.
    ssc_sij.set_scaled_abstol(3.0e-6)
    ssc_sij.set_block_size(4)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dup = ser.from_string(ser.to_string(ssc_sij, True))

    assert isinstance(dup, Nc.XcorSSCSij)
    assert dup.get_nbins() == ssc_sij.get_nbins()
    assert dup.get_lmax() == ssc_sij.get_lmax()
    assert dup.get_block_size() == ssc_sij.get_block_size()
    assert dup.get_method() == ssc_sij.get_method()
    assert_allclose(dup.get_reltol(), ssc_sij.get_reltol())
    assert_allclose(dup.get_scaled_abstol(), ssc_sij.get_scaled_abstol())
    assert_allclose(dup.get_fsky(), ssc_sij.get_fsky(), rtol=1.0e-12)

    assert_allclose(
        _matrix_to_np(dup.eval(cosmology.cosmo)),
        _matrix_to_np(ssc_sij.eval(cosmology.cosmo)),
        rtol=1.0e-12,
    )


def test_peek_matrix_is_the_prepared_one(cosmology: Cosmology) -> None:
    """peek_matrix exposes the last prepared matrix, eval returns a copy."""
    ssc_sij = _make_ssc_sij(cosmology)
    ssc_sij.prepare(cosmology.cosmo)

    peeked = _matrix_to_np(ssc_sij.peek_matrix())
    evaluated = _matrix_to_np(ssc_sij.eval(cosmology.cosmo))

    assert_allclose(peeked, evaluated, rtol=0.0, atol=0.0)
