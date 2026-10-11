#!/usr/bin/env python3
#
# make_xcor_window_dense_bin.py
#
# Sun Oct 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_xcor_window_dense_bin.py
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

"""Write the dense certified windows as one NcmSerialize file for the C tests.

The certified values are the gzipped JSON tables of the dense window run, one per
Bessel-derivative order d = 0, 1, 2, with the convention of
``make_xcor_window_truth_table.py``:

    I_ell^(d)(k) = int W(chi) j_ell^(d)(k chi) dchi,  chi in Mpc, k in 1/Mpc,

with W the normalized analytic window, on 50 k per decade over the library closure's
domain of each (shape, ell), certified to a relative radius of 1e-25. The three orders
share one k grid.

The output is one NcmObjDictStr holding:

- ``kernel/<shape>/d<d>``: the analytic kernel of each shape, as the Python window tests
  build it, with #NcXcorKernelRadial:bessel-deriv set to d (construct-only, so set in its
  serialized form);
- ``table/<shape>|<ell>``: an NcmMatrix, one row per k, columns (k, I^(0), I^(1), I^(2)),
  as doubles.

It is a file of the NumCosmo data-file release, too large for the repository:

    python3 tests/tools/make_xcor_window_dense_bin.py <dir with the json.gz> \\
        xcor_window_ilk_dense.bin
"""

from __future__ import annotations

import gzip
import json
import pathlib
import sys

sys.path.insert(
    0, str(pathlib.Path(__file__).resolve().parents[1] / "python" / "nc" / "xcor")
)

from numcosmo_py import Nc, Ncm  # noqa: E402

from test_xcor_window_truth_table import _make_kernel  # noqa: E402

DERIVS = (0, 1, 2)


def _cells(path: pathlib.Path) -> tuple[dict, dict]:
    """Every certified (shape, ell) cell of one table, finished or partial, and the shapes."""
    payload = json.load(gzip.open(path, "rt"))
    cells = {}
    shapes = {}

    for shape, s in payload["shapes"].items():
        shapes[shape] = s
        for i, ell in enumerate(payload["ells"]):
            cells[(shape, ell)] = (s["kvals"][i], s["table"][i])

    for key, cell in payload.get("partial", {}).items():
        shape, ell = key.split("|")
        shapes.setdefault(shape, {"shape": cell["shape_form"], "ctor": cell["ctor"]})
        cells[(shape, int(ell))] = (cell["kvals"], cell["table"])

    return cells, shapes


def main() -> None:
    """Read the three tables and write the NcmSerialize file."""
    src = pathlib.Path(sys.argv[1])
    out = sys.argv[2]

    Ncm.cfg_init()
    tables = [_cells(src / f"xcor_window_ilk_dense_d{d}.json.gz") for d in DERIVS]
    shapes = tables[0][1]
    keys = sorted(set.intersection(*(set(cells) for cells, _ in tables)))

    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(5.0)
    dist.prepare(cosmo)
    ps = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.LCDM
    )

    root = Ncm.ObjDictStr.new()
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    n_shapes = 0
    for shape in sorted({shape for shape, _ in keys}):
        s = shapes[shape]
        kernel = _make_kernel(
            s["shape"], dist, ps, Ncm.SBesselIntegratorLevin.new(0, 8), s["ctor"]
        )
        text = ser.to_string(kernel, True)
        if text.count("'bessel-deriv': <uint32 0>") != 1:
            raise ValueError(f"{shape}: no bessel-deriv in the serialized kernel")
        n_shapes += 1
        for d in DERIVS:
            ser.reset(True)
            root.add(
                f"kernel/{shape}/d{d}",
                ser.from_string(
                    text.replace(
                        "'bessel-deriv': <uint32 0>", f"'bessel-deriv': <uint32 {d}>"
                    )
                ),
            )

    for shape, ell in keys:
        kvals = tables[0][0][(shape, ell)][0]
        for cells, _ in tables[1:]:
            if cells[(shape, ell)][0] != kvals:
                raise ValueError(f"{shape} ell = {ell}: the orders' k grids differ")

        m = Ncm.Matrix.new(len(kvals), 1 + len(DERIVS))
        for j, k in enumerate(kvals):
            m.set(j, 0, k)
            for d, (cells, _) in enumerate(tables):
                m.set(j, 1 + d, float(cells[(shape, ell)][1][j]))
        root.add(f"table/{shape}|{ell}", m)

    ser.reset(True)
    ser.dict_str_to_binfile(root, out)
    print(f"{len(keys)} (shape, ell) tables of {n_shapes} shapes written to {out}")


if __name__ == "__main__":
    main()
