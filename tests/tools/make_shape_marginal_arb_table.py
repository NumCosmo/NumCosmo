#!/usr/bin/env python
#
# make_shape_marginal_arb_table.py
#
# Thu Nov 12 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_shape_marginal_arb_table.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Write the certified table of the shape marginal for the data-file release.

Reads the values of

    ln P (eps_obs | g) = ln int_{|chi_I| < 1} P_2D (chi_I) N_2 (eps_obs - f_g (chi_I)) d^2 chi_I

certified by ``tests/tools/nc_galaxy_shape_marginal_arb.c``, one point per line with the
whitespace-separated columns ``conv pop p1 p2 std_noise g1 g2 e1 e2 lnP``: the shear
convention (``trace`` or ``trace_det``), the population (``gauss`` with sigma = p1, or
``beta`` with alpha = p1 and beta = p2), the noise, the shear, the observed ellipticity and
the certified value, ``nan`` where the tool did not certify the point. Every certified
point is kept and its value rounded to the nearest double.

Writes an NcmObjDictStr with one entry, ``points``: an NcmMatrix with one row per point
and the columns ``conv, pop, p1, p2, std_noise, g1, g2, e1, e2, lnP``, where ``conv`` is
the value of NcGalaxyWLObsEllipConv and ``pop`` is 0 for the truncated Gaussian and 1 for
the Beta population. The file is the asset ``nc_galaxy_shape_marginal_arb_cert.bin`` of
the data-file release; its consumer is
``tests/c/nc/lss/galaxy/test_nc_galaxy_shape_factor_spectral_cert.c``.

    python3 tests/tools/make_shape_marginal_arb_table.py points.txt
"""

import argparse
import math
from pathlib import Path

from numcosmo_py import Nc, Ncm

CONV = {
    "trace": int(Nc.GalaxyWLObsEllipConv.TRACE),
    "trace_det": int(Nc.GalaxyWLObsEllipConv.TRACE_DET),
}
POP = {"gauss": 0, "beta": 1}


def main() -> None:
    """Write the table."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("points", type=Path, help="certified points, one per line")
    parser.add_argument(
        "--output", type=Path, default=Path("nc_galaxy_shape_marginal_arb_cert.bin")
    )
    args = parser.parse_args()

    data: list[float] = []
    n_rows = 0
    n_skipped = 0
    for line in args.points.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        conv, pop, *values = fields
        ln_p = float(values[-1])
        if not math.isfinite(ln_p):
            n_skipped += 1
            continue
        data.extend([float(CONV[conv]), float(POP[pop])])
        data.extend(float(x) for x in values)
        n_rows += 1

    Ncm.cfg_init()
    table = Ncm.ObjDictStr.new()
    table.add("points", Ncm.Matrix.new_array(data, 10))
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    ser.dict_str_to_binfile(table, str(args.output))
    print(f"{n_rows} certified points written to {args.output}, {n_skipped} skipped")


if __name__ == "__main__":
    main()
