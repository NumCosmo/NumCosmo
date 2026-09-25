#
# measure.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# measure.py
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

"""NumCosmo measured against its own certified xcor references.

:func:`radial` computes the radial integral of every cell of a certified window table
with ``NcmSBesselIntegratorLevin``, the way ``test_xcor_window_truth_table.py`` does,
and returns the measurements as a :class:`~numcosmo_py.certified.results.Campaign`.
"""

import pathlib
import typing

import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.certified import results, tables
from numcosmo_py.certified.xcor import kernels
from numcosmo_py.certified.xcor import windows as W


def _window(meta: dict) -> W.Window:
    """Return the window a table group's metadata defines."""
    params = tuple(
        tuple(meta[name]) if isinstance(meta[name], list) else meta[name]
        for name in W.FORM_PARAMS[meta["form"]]
    )

    return W.Window(meta["form"], params, "")


def _integrate(
    kernel,
    ells: np.ndarray,
    ks: np.ndarray,
    deriv: int,
    reltol: float | None,
) -> np.ndarray:
    """Return the library's I_ell(k) of one window at each (ell, k) row."""
    supports = [kernel.get_comp_support(c) for c in range(kernel.get_n_comps())]
    block = Ncm.Vector.new(1)
    measured = np.empty(len(ells))

    for ell in np.unique(ells):
        sbi = Ncm.SBesselIntegratorLevin.new(int(ell), int(ell))

        if reltol is not None:
            sbi.set_reltol(reltol)

        for row in np.flatnonzero(ells == ell):
            total = 0.0

            for comp, (low, high) in enumerate(supports):
                sbi.integrate_deriv(
                    lambda chi, _k, w=kernel, c=comp, lo=low, hi=high: w.eval_W_comp(
                        c, min(max(chi, lo), hi)
                    ),
                    low,
                    high,
                    float(ks[row]),
                    deriv,
                    block,
                )
                total += block.get(0)

            measured[row] = total

    return measured


def peaks(table: Ncm.VarDict) -> dict[tuple[str, int], float]:
    """Return the largest certified |I| of each window and multipole of a table.

    The scale a radial deviation is measured against: the peak of I_ell(k) sits near
    x_nu = 1, so only a table covering that point supplies it.
    """
    out: dict[tuple[str, int], float] = {}

    for name, (_, cells) in tables.groups(table).items():
        ells = np.asarray(cells["ell"], dtype=np.int64)
        values = np.abs(np.asarray(cells["value"]))

        for ell in np.unique(ells):
            out[(name, int(ell))] = float(values[ells == ell].max())

    return out


def radial(
    table_path: pathlib.Path,
    code_version: str,
    peak_table: Ncm.VarDict,
    peak_name: str,
    reltol: float | None = None,
) -> results.Campaign:
    """Measure NumCosmo's radial integral against a certified window table.

    Each component is integrated over its own support with chi clamped to it, as the
    library's own path does; the kernel is built from the table's definition of the
    window, not from the vocabulary entry of the same name. The deviation is scaled by
    the peak of the same window and multipole in ``peak_table``, which must cover
    x_nu = 1: a table of the tail alone has no peak, and its own largest value would
    turn roundoff into an apparent error.

    :param table_path: a window table in the :mod:`numcosmo_py.certified.tables` layout.
    :param code_version: the NumCosmo commit or version being measured.
    :param peak_table: the table the scale is taken from (see :func:`peaks`).
    :param peak_name: its name, recorded in the campaign's configuration.
    :param reltol: the integrator's relative tolerance; its default when omitted.
    :return: the campaign, not yet written.
    """
    table = tables.read(table_path)
    meta = tables.metadata(table)
    deriv = int(meta["bessel_deriv"])
    scales = peaks(peak_table)
    config: dict[str, typing.Any] = {
        "integrator": "NcmSBesselIntegratorLevin",
        "bessel_deriv": deriv,
        "scale": "peak of the window and multipole",
        "scale_table": peak_name,
    }

    if reltol is not None:
        config["reltol"] = reltol

    campaign = results.Campaign(table_path, "numcosmo", code_version, config)
    cosmo, dist, ps = kernels.make_cosmo_bits()

    for name, (gmeta, cells) in tables.groups(table).items():
        kernel = kernels.build_kernel(_window(gmeta), dist, ps)
        kernel.set_l_limber(-1)
        kernel.prepare(cosmo)
        ells = np.asarray(cells["ell"], dtype=np.int64)
        ks = np.asarray(cells["k"])

        campaign.add(
            name,
            {
                "ell": ells,
                "k": ks,
                "reference": np.asarray(cells["value"]),
                "radius": np.asarray(cells["radius"]),
                "measured": _integrate(kernel, ells, ks, deriv, reltol),
                "scale": np.array([scales[(name, int(e))] for e in ells]),
            },
        )

    return campaign
