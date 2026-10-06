#!/usr/bin/env python
#
# make_sparam_desc_fixtures.py
#
# Sat Oct 3 10:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_sparam_desc_fixtures.py
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

"""Freeze cosmologies with modified parameter descriptions, written by NumCosmo 0.27.

Run this once against NumCosmo 0.27, where Yp is still a parameter of the
cosmology classes. NcmModel:sparam-array keys each modified description by
parameter index, and removing Yp shifted every later index, so these files are
the input the index-to-name migration has to read. Re-running it against a newer
library does not reproduce them.

Writes, into data/truth_tables/bbn/:
  <case>.obj    GVariant text
  <case>.bin    GVariant binary
  <case>.yaml   YAML
  sparam_desc.json   each parameter's description and value, by name
"""

import json
import os

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "truth_tables", "bbn"
)


def modify(cosmo, name, lower, upper, scale, free):
    """Change the description of parameter @name so that it is serialized."""
    ok, idx = cosmo.param_index_from_name(name)
    assert ok
    cosmo.param_set_lower_bound(idx, lower)
    cosmo.param_set_upper_bound(idx, upper)
    cosmo.param_set_scale(idx, scale)
    cosmo.param_set_ftype(idx, Ncm.ParamType.FREE if free else Ncm.ParamType.FIXED)


def de_cpl():
    """NcHICosmoDECpl with descriptions modified after the Yp index, and on Yp."""
    cosmo = Nc.HICosmoDECpl.new()
    modify(cosmo, "ENnu", 2.0, 4.0, 0.11, False)
    modify(cosmo, "Omegab", 0.01, 0.09, 0.0012, True)
    modify(cosmo, "w0", -2.5, -0.3, 0.013, True)
    modify(cosmo, "w1", -0.9, 0.9, 0.017, False)
    modify(cosmo, "Yp", 0.2, 0.3, 0.0021, False)
    return cosmo


def lcdm():
    """NcHICosmoLCDM with only ENnu modified, the index 5 slot of Omegab after Yp left."""
    cosmo = Nc.HICosmoLCDM.new()
    modify(cosmo, "ENnu", 2.5, 3.5, 0.07, True)
    return cosmo


CASES = {"de_cpl_desc": de_cpl(), "lcdm_desc": lcdm()}


def describe(cosmo):
    """Each parameter's description and value, by name, without Yp."""
    out = {}
    for i in range(cosmo.sparam_len()):
        name = cosmo.param_name(i)
        if name == "Yp":
            continue
        out[name] = {
            "lower": cosmo.param_get_lower_bound(i),
            "upper": cosmo.param_get_upper_bound(i),
            "scale": cosmo.param_get_scale(i),
            "free": cosmo.param_get_ftype(i) == Ncm.ParamType.FREE,
            "value": cosmo.param_get(i),
        }
    return out


def main():
    ser = Ncm.Serialize.new(0)
    descs = {}

    for name, obj in CASES.items():
        ser.reset(True)
        ser.to_file(obj, os.path.join(OUT, f"{name}.obj"))
        ser.reset(True)
        ser.to_binfile(obj, os.path.join(OUT, f"{name}.bin"))
        ser.reset(True)
        ser.to_yaml_file(obj, os.path.join(OUT, f"{name}.yaml"))
        descs[name] = {"type": obj.__gtype__.name, "params": describe(obj)}

    with open(os.path.join(OUT, "sparam_desc.json"), "w", encoding="utf-8") as f:
        json.dump(descs, f, indent=2, sort_keys=True)
        f.write("\n")


if __name__ == "__main__":
    main()
