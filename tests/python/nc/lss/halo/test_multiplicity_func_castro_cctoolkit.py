#!/usr/bin/env python
#
# test_multiplicity_func_castro_cctoolkit.py
#
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

"""Compare NcMultiplicityFuncCastro against the CCToolkit reference.

How the truth tables were produced
----------------------------------
The multiplicity function is a pure function of the peak height ``nu``, the
variance slope ``s = dlnsigma/dlnR``, and the background quantities
``Omega_m(z)`` (plus ``Omega_de(z_ta)`` and ``w(z_ta)`` for the 2025 model).
The stored references were produced by evaluating CCToolkit's
``cctoolkit/hmf.py`` (https://github.com/TiagoBsCastro/CCToolkit,
``multiplicity_function_castro23`` / ``multiplicity_function_castro25``) at
exactly those quantities:

* ``s`` is the analytic slope of a power-law spectrum, ``-(n+3)/2``, which is
  what ``_make_psf`` below builds. It is handed to the reference directly.
* ``Omega_m(z)``, ``Omega_de(z_ta)`` and ``w(z_ta)`` come from the same
  ``_make_cosmo`` cosmology reproduced here, with ``z_ta`` from the same
  Einstein-de Sitter relation.
* ``nu = delta_c / sigma`` with the Kitayama & Suto ``delta_c(Omega_m(z))``.

Run this module as a script to rewrite the tables.

What the tolerances measure
---------------------------
The model is a closed-form function of ``sigma``, the slope and the background,
so ``eval_full`` evaluates it with the slope supplied directly. Compared that
way the two implementations agree to 1.6e-15, i.e. to round-off, and that is
what ``test_*_matches_cctoolkit`` asserts.

The separate ``test_*_through_filter`` cases run the same comparison through
``eval``, which derives the slope from the filter instead. They are far looser
because the residual is then dominated by how well the filter reproduces the
analytic power-law slope: it is sensitive to the infrared range, and a k range
starting at 1e-6 rather than 1e-10 degrades the agreement from 5e-5 to 2e-3,
entirely through the shallow-slope case. Those cases check the wiring, not the
model.
"""

import math

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

HALO_FINDERS = [
    Nc.MultiplicityFuncCastroHaloFinder.ROCKSTAR,
    Nc.MultiplicityFuncCastroHaloFinder.AHF,
    Nc.MultiplicityFuncCastroHaloFinder.SUBFIND,
    Nc.MultiplicityFuncCastroHaloFinder.VELOCIRAPTOR,
]

_C23_TABLE = "truth_tables/halo/nc_multiplicity_func_castro_c23.bin"
_C25_TABLE = "truth_tables/halo/nc_multiplicity_func_castro_c25.bin"

# eval_full takes the slope directly, leaving only round-off.
RTOL = 1.0e-13
# Through the filter the slope is rebuilt from a tabulated power law; measured
# at 5.2e-5 worst case with the settings below. See the module docstring.
RTOL_FILTER = 2.0e-4


def _load_table(name: str) -> np.ndarray:
    """Load a truth table as an (nrows, ncols) array."""
    path = Ncm.cfg_get_data_filename(name, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(matrix.nrows(), matrix.ncols())


C23_GOLDEN = _load_table(_C23_TABLE)
C25_GOLDEN = _load_table(_C25_TABLE)


def _make_cosmo(w: float = -1.0) -> Nc.HICosmo:
    """Return the cosmology the truth tables were generated with."""
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_by_name("Omegac", 0.25)
    cosmo.param_set_by_name("Omegab", 0.05)
    cosmo.param_set_by_name("Omegak", 0.0)
    cosmo.param_set_by_name("w", w)
    return cosmo


def _make_psf(slope: float, cosmo: Nc.HICosmo) -> Ncm.PowspecFilter:
    """Return a filter whose variance slope is the requested one.

    A power law with index n gives dlnsigma/dlnR = -(n+3)/2 exactly.
    """
    n_index = -2.0 * slope - 3.0
    # The infrared cut dominates the slope accuracy for shallow spectra; see
    # the module docstring.
    ka = np.geomspace(1.0e-10, 1.0e3, 6000)
    za = np.linspace(0.0, 5.0, 20, dtype=np.float64)
    za_v, ka_v = np.meshgrid(za, ka)
    lnPk = np.log(1.0e-9 * ka_v**n_index) + 0.0 * za_v

    Pk2d = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(za.tolist()),
        y_vector=Ncm.Vector.new_array(np.log(ka).tolist()),
        z_matrix=Ncm.Matrix.new_array(lnPk.flatten().tolist(), len(za)),
    )
    ps = Ncm.PowspecSpline2d.new(Pk2d)
    ps.prepare()

    psf = Ncm.PowspecFilter.new(ps, Ncm.PowspecFilterType.TOPHAT)
    psf.set_reltol(1.0e-9)
    psf.set_best_lnr0()
    psf.prepare(cosmo)
    return psf


@pytest.mark.parametrize("row", C23_GOLDEN)
def test_c23_matches_cctoolkit(row: np.ndarray) -> None:
    """The 2023 calibration reproduces the CCToolkit reference to round-off."""
    finder_idx, slope, z, sigma, expected = row
    mc = Nc.MultiplicityFuncCastro.new_full(
        Nc.MultiplicityFuncCastroModel.C23, HALO_FINDERS[int(finder_idx)]
    )

    got = mc.eval_full(_make_cosmo(), sigma, slope, z)
    assert got == pytest.approx(expected, rel=RTOL)


@pytest.mark.parametrize("row", C25_GOLDEN)
def test_c25_matches_cctoolkit(row: np.ndarray) -> None:
    """The 2025 calibration reproduces the CCToolkit reference to round-off."""
    w, slope, z, sigma, expected = row
    mc = Nc.MultiplicityFuncCastro.new_full(
        Nc.MultiplicityFuncCastroModel.C25,
        Nc.MultiplicityFuncCastroHaloFinder.ROCKSTAR,
    )

    got = mc.eval_full(_make_cosmo(w), sigma, slope, z)
    assert got == pytest.approx(expected, rel=RTOL)


@pytest.mark.parametrize("row", C23_GOLDEN[::12])
def test_c23_through_filter(row: np.ndarray) -> None:
    """Going through the filter reproduces eval_full, checking the wiring."""
    finder_idx, slope, z, sigma, expected = row
    cosmo = _make_cosmo()
    psf = _make_psf(slope, cosmo)

    mc = Nc.MultiplicityFuncCastro.new_full(
        Nc.MultiplicityFuncCastroModel.C23, HALO_FINDERS[int(finder_idx)]
    )
    mc.set_psf(psf)

    got = mc.eval(cosmo, sigma, math.log(8.0), z)
    assert got == pytest.approx(expected, rel=RTOL_FILTER)


@pytest.mark.parametrize("row", C25_GOLDEN[::12])
def test_c25_through_filter(row: np.ndarray) -> None:
    """Going through the filter reproduces eval_full, checking the wiring."""
    w, slope, z, sigma, expected = row
    cosmo = _make_cosmo(w)
    psf = _make_psf(slope, cosmo)

    mc = Nc.MultiplicityFuncCastro.new_full(
        Nc.MultiplicityFuncCastroModel.C25,
        Nc.MultiplicityFuncCastroHaloFinder.ROCKSTAR,
    )
    mc.set_psf(psf)

    got = mc.eval(cosmo, sigma, math.log(8.0), z)
    assert got == pytest.approx(expected, rel=RTOL_FILTER)


def _regenerate() -> None:  # pragma: no cover - developer tool
    """Rewrite both truth tables, evaluating CCToolkit at our own inputs.

    Requires CCToolkit on the path. ``hmf.py`` is pure numpy/scipy, so it is
    loaded directly to avoid the package __init__ pulling in CAMB.
    """
    # pylint: disable=too-many-locals
    import importlib.util  # pylint: disable=import-outside-toplevel

    spec = importlib.util.spec_from_file_location(
        "cc_hmf", "../CCToolkit/cctoolkit/hmf.py"
    )
    hmf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hmf)

    finders = [
        hmf.best_fit_values_ROCKSTAR,
        hmf.best_fit_values_AHF,
        hmf.best_fit_values_SUBFIND,
        hmf.best_fit_values_VELOCIraptor,
    ]
    slopes = sorted(set(C23_GOLDEN[:, 1]))
    zs = sorted(set(C23_GOLDEN[:, 2]))
    sigmas = sorted(set(C23_GOLDEN[:, 3]))
    ws = sorted(set(C25_GOLDEN[:, 0]))

    def delta_c(cosmo: Nc.HICosmo, z: float) -> float:
        om = cosmo.E2Omega_m(z) / cosmo.E2(z)
        return (
            3.0
            / 20.0
            * (12.0 * math.pi) ** (2.0 / 3.0)
            * (1.0 + 0.012299 * math.log10(om))
        )

    rows23 = []
    cosmo = _make_cosmo()
    for i, bf in enumerate(finders):
        for slope in slopes:
            for z in zs:
                om = cosmo.E2Omega_m(z) / cosmo.E2(z)
                d_c = delta_c(cosmo, z)
                for sigma in sigmas:
                    f = float(
                        hmf.multiplicity_function_castro23(d_c / sigma, slope, om, bf)
                    )
                    rows23.append([float(i), slope, z, sigma, f])

    rows25 = []
    for w in ws:
        cosmo = _make_cosmo(w)
        for slope in slopes:
            for z in zs:
                om = cosmo.E2Omega_m(z) / cosmo.E2(z)
                d_c = delta_c(cosmo, z)
                z_ta = (1.0 + z) * 2.0 ** (2.0 / 3.0) - 1.0
                om_de = cosmo.E2Omega_de(z_ta) / cosmo.E2(z_ta)
                w_ta = cosmo.w_de(z_ta)
                for sigma in sigmas:
                    f = float(
                        hmf.multiplicity_function_castro25(
                            d_c / sigma,
                            slope,
                            om,
                            om_de,
                            w_ta,
                            hmf.best_fit_values_castro25,
                        )
                    )
                    rows25.append([w, slope, z, sigma, f])

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    for name, rows in ((_C23_TABLE, rows23), (_C25_TABLE, rows25)):
        arr = np.array(rows, dtype=np.float64)
        matrix = Ncm.Matrix.new_array(arr.flatten().tolist(), arr.shape[1])
        ser.to_binfile(matrix, f"data/{name}")
        print(f"wrote data/{name}: {arr.shape[0]} x {arr.shape[1]}")


if __name__ == "__main__":  # pragma: no cover
    _regenerate()
