#!/usr/bin/env python
"""Extract test photometric-redshift densities P(z) from the HSC PDR1 weak-lensing catalogs.

The fixture holds, for each of the five NcGalaxyWLObsCatalogId catalogs, the P(z)
splines of the NGAL galaxies with the most separate supported peaks (runs of positive
knot values separated by zeros), ties broken by galaxy index. Densities with several
peaks separated by zero density are the ones an inverse-CDF sampler steps over.

Output: data/hsc_pdr1_pz_multipeak.bin, an NcmObjDictStr mapping
"<catalog nick>:<galaxy index>" to the catalog's NcmSpline. Needs the catalogs, which
Nc.GalaxyWLObs.new_from_catalog_id() downloads.
"""

import os

import numpy as np

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..",
    "..",
    "data",
    "hsc_pdr1_pz_multipeak.bin",
)

NGAL = 2


def n_peaks(pz: Ncm.Spline) -> int:
    """Number of runs of positive knot values of pz."""
    positive = np.array(pz.peek_yv().dup_array()) > 0.0
    return int(positive[0]) + int(np.sum(positive[1:] & ~positive[:-1]))


def main() -> None:
    """Write the fixture."""
    ods = Ncm.ObjDictStr.new()

    for cid in range(len(Nc.GalaxyWLObsCatalogId.__enum_values__)):
        catalog = Nc.GalaxyWLObsCatalogId(cid)
        obs = Nc.GalaxyWLObs.new_from_catalog_id(catalog)
        peaks = np.array([n_peaks(obs.peek_pz(i)) for i in range(obs.len())])
        chosen = sorted(np.lexsort((np.arange(len(peaks)), -peaks))[:NGAL])

        for i in chosen:
            ods.add(f"{catalog.value_nick}:{int(i)}", obs.peek_pz(int(i)))

        print(
            f"{catalog.value_nick}: galaxies {[int(i) for i in chosen]}, "
            f"peaks {[int(peaks[i]) for i in chosen]}"
        )

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    ser.dict_str_to_binfile(ods, OUT)
    print(f"{ods.len()} P(z) written to {os.path.normpath(OUT)}")


if __name__ == "__main__":
    main()
