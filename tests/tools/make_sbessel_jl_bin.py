#!/usr/bin/env python3
#
# make_sbessel_jl_bin.py
#
# Wed Oct 14 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_sbessel_jl_bin.py
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

"""Write an sbessel truth table (gzipped JSON) as an NcmVarDict binfile for the C tests.

The JSON holds

    table[l][i] = int_lower-bound^upper-bound f(x) j_l(kvals[i] x) dx

for the f named by the file (Gaussian or rational, with ``center`` and ``std``), with
``lvals`` and ``kvals`` its grid and ``convention`` the statement above. The binfile
has the same keys: ``center``, ``std``, ``lower-bound``, ``upper-bound`` (doubles),
``kvals`` (double array), ``lvals`` (int array), ``table`` (NcmMatrix, one row per
multipole, one column per k) and ``convention``, plus ``source`` (the JSON file name),
``source-sha256`` and ``generator`` (this script):

    python3 tests/tools/make_sbessel_jl_bin.py gauss_jl_500.json.gz gauss_jl_500.bin

The output is a file of the data-file release; the JSON sources are kept outside the
repository.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import pathlib
import sys

from numcosmo_py import Ncm


def main() -> None:
    """Convert one table."""
    src = pathlib.Path(sys.argv[1])
    out = sys.argv[2]

    Ncm.cfg_init()
    payload = json.load(gzip.open(src, "rt"))
    lvals = [int(ell) for ell in payload["lvals"]]
    kvals = [float(k) for k in payload["kvals"]]
    rows = payload["table"]

    if lvals != list(range(lvals[0], lvals[0] + len(lvals))):
        raise ValueError("lvals must be consecutive")
    if len(rows) != len(lvals) or any(len(row) != len(kvals) for row in rows):
        raise ValueError("table shape does not match lvals x kvals")

    table = Ncm.Matrix.new(len(lvals), len(kvals))
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            table.set(i, j, float(value))

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    vd = Ncm.VarDict.new()
    for key in ("center", "std", "lower-bound", "upper-bound"):
        vd.set_double(key, float(payload[key]))
    vd.set_double_array("kvals", kvals)
    vd.set_int_array("lvals", lvals)
    vd.set_object("table", ser, table)
    vd.set_string("convention", payload["convention"])
    vd.set_string("source", src.name)
    vd.set_string("source-sha256", hashlib.sha256(src.read_bytes()).hexdigest())
    vd.set_string("generator", "tests/tools/make_sbessel_jl_bin.py")

    ser.var_dict_to_variant_file(vd, out, True)
    print(f"{len(lvals)} x {len(kvals)} table from {src.name} written to {out}")


if __name__ == "__main__":
    main()
