#!/usr/bin/env python
"""Freeze pre-NcBBN cosmologies as serialization fixtures.

Run this against a library that still carries ``Yp`` as a parameter of the
cosmology classes, i.e. *before* the NcBBN submodel migration. The files it
writes are the "old style" inputs the migration has to keep reading; once
NcBBN lands this script no longer reproduces them, and it should not be
re-run to "refresh" the fixtures -- that would defeat their purpose.

Writes, into data/truth_tables/bbn/:
  <case>.obj    GVariant text
  <case>.bin    GVariant binary
  <case>.yaml   YAML
  golden.json   what each case must still evaluate to after the migration
"""

import json
import os

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "truth_tables", "bbn"
)


def set_yp(cosmo, value, free):
    """Set Yp and its fit type on a cosmology."""
    ok, idx = cosmo.param_index_from_name("Yp")
    assert ok
    cosmo.param_set(idx, value)
    cosmo.param_set_ftype(
        idx, Ncm.ParamType.FREE if free else Ncm.ParamType.FIXED
    )


def de_xcdm(yp, free):
    """An NcHICosmoDEXcdm with a deliberately off-default baryon density."""
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.param_set_by_name("H0", 67.36)
    cosmo.param_set_by_name("Omegab", 0.04930)
    cosmo.param_set_by_name("Omegac", 0.26440)
    cosmo.param_set_by_name("ENnu", 3.046)
    set_yp(cosmo, yp, free)
    return cosmo


def lcdm(yp):
    """NcHICosmoLCDM never had the BBN branch: Yp is always the parameter."""
    cosmo = Nc.HICosmoLCDM.new()
    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_by_name("Omegab", 0.045)
    set_yp(cosmo, yp, False)
    return cosmo


def mset_with_submodels():
    """A realistic saved run: cosmo + reionization + primordial submodels."""
    cosmo = de_xcdm(0.24, free=False)
    cosmo.add_submodel(Nc.HIReionCamb.new())
    cosmo.add_submodel(Nc.HIPrimPowerLaw.new())
    return Ncm.MSet.new_array([cosmo])


def peek_cosmo(obj):
    """The cosmology inside obj, whether it is an mset or a cosmology."""
    if isinstance(obj, Ncm.MSet):
        return obj.peek(Nc.HICosmo.id())
    return obj


CASES = {
    "de_xcdm_yp_fixed": de_xcdm(0.24, free=False),
    "de_xcdm_yp_free": de_xcdm(0.2521, free=True),
    "lcdm_yp": lcdm(0.2478),
    "mset_de_xcdm": mset_with_submodels(),
}


def main():
    ser = Ncm.Serialize.new(0)
    golden = {}

    for name, obj in CASES.items():
        ser.reset(True)
        ser.to_file(obj, os.path.join(OUT, f"{name}.obj"))
        ser.reset(True)
        ser.to_binfile(obj, os.path.join(OUT, f"{name}.bin"))
        ser.reset(True)
        ser.to_yaml_file(obj, os.path.join(OUT, f"{name}.yaml"))

        cosmo = peek_cosmo(obj)
        ok, idx = cosmo.param_index_from_name("Yp")
        assert ok
        golden[name] = {
            "type": cosmo.__gtype__.name,
            "Yp_4He": Nc.HICosmo.Yp_4He(cosmo),
            "Yp_param": cosmo.param_get(idx),
            "Yp_free": cosmo.param_get_ftype(idx) == Ncm.ParamType.FREE,
            "Omega_b0h2": Nc.HICosmo.Omega_b0h2(cosmo),
            # Neff goes through E2Press_mnu, which NcHICosmoLCDM does not
            # implement -- and LCDM never had the BBN branch anyway.
            "Neff": (
                Nc.HICosmo.Neff(cosmo) if isinstance(cosmo, Nc.HICosmoDE) else None
            ),
            "submodels": [
                cosmo.peek_submodel(i).__gtype__.name
                for i in range(cosmo.get_submodel_len())
            ],
        }
        print(f"{name:20s} Yp_4He={golden[name]['Yp_4He']:.10f} "
              f"free={golden[name]['Yp_free']} "
              f"submodels={golden[name]['submodels']}")

    with open(os.path.join(OUT, "golden.json"), "w", encoding="utf-8") as f:
        json.dump(golden, f, indent=2, sort_keys=True)
        f.write("\n")


if __name__ == "__main__":
    main()
