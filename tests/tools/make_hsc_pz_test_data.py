#!/usr/bin/env python
"""Extract test photometric-redshift densities P(z) from the HSC PDR1 weak-lensing catalogs.

The fixture holds, for each of the five NcGalaxyWLObsCatalogId catalogs, the P(z)
tables of the NGAL galaxies with the most separate supported peaks (runs of positive
table values separated by zeros), ties broken by galaxy index. Densities with several
peaks separated by zero density are the ones an inverse-CDF sampler steps over.

Output: data/hsc_pdr1_pz_multipeak.bin, an NcmObjDictStr mapping
"<catalog nick>:<galaxy index>" to the NcmSplineCubicNotaknot of the table. Needs the
catalogs, which Nc.GalaxyWLObs.new_from_catalog_id() downloads.

The catalogs hold each P(z) as a NcGalaxyPzBSplineSqrt, not as its table. The table
is rebuilt from it: its redshifts are equally spaced over the range of the density, as
many as the control points of its B-spline, and its values are the density there. A
table value was zero where the rebuilt one is below ZERO_THRESHOLD of the largest: the
tables are cut at 0.005 of their largest value, and a rebuilt value at a zero is below
1e-28 of it. Such values are written as zero.

The committed fixture was written from the tables themselves. With --check the script
compares the rebuilt fixture with it and writes nothing: the keys are equal, the zeros
are the same, and the redshifts and values agree within CHECK_TOLERANCE of the largest
value of each table.
"""

import argparse
import os
import sys

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

# A rebuilt value below this fraction of the largest one is a zero of the table.
ZERO_THRESHOLD = 1.0e-20

# Largest difference accepted by --check; the largest found is 1.8e-14.
CHECK_TOLERANCE = 1.0e-13


def table(pz: Nc.GalaxyPzBSplineSqrt) -> tuple[np.ndarray, np.ndarray]:
    """Rebuild the table of pz: its redshifts and its values, zeros included."""
    z_min, z_max = pz.get_range()
    z = np.linspace(z_min, z_max, pz.peek_bspline().peek_control_points().len())
    p = np.array([pz.eval_p(z_j) for z_j in z])
    p[p < ZERO_THRESHOLD * p.max()] = 0.0
    return z, p


def n_peaks(p: np.ndarray) -> int:
    """Number of runs of positive values of the table p."""
    positive = p > 0.0
    return int(positive[0]) + int(np.sum(positive[1:] & ~positive[:-1]))


def build() -> Ncm.ObjDictStr:
    """Build the fixture from the catalogs."""
    ods = Ncm.ObjDictStr.new()

    for cid in range(len(Nc.GalaxyWLObsCatalogId.__enum_values__)):
        catalog = Nc.GalaxyWLObsCatalogId(cid)
        obs = Nc.GalaxyWLObs.new_from_catalog_id(catalog)
        tables = [table(obs.peek_pz(i)) for i in range(obs.len())]
        peaks = np.array([n_peaks(p) for _, p in tables])
        chosen = sorted(np.lexsort((np.arange(len(peaks)), -peaks))[:NGAL])

        for i in chosen:
            z, p = tables[i]
            ods.add(
                f"{catalog.value_nick}:{int(i)}",
                Ncm.SplineCubicNotaknot.new_full(
                    Ncm.Vector.new_array(z.tolist()),
                    Ncm.Vector.new_array(p.tolist()),
                    True,
                ),
            )

        print(
            f"{catalog.value_nick}: galaxies {[int(i) for i in chosen]}, "
            f"peaks {[int(peaks[i]) for i in chosen]}"
        )

    return ods


def check(ods: Ncm.ObjDictStr) -> bool:
    """Compare the rebuilt fixture with the committed one."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    committed = ser.dict_str_from_binfile(OUT)
    equal = sorted(ods.keys()) == sorted(committed.keys())
    largest = 0.0

    if not equal:
        print(f"keys differ: {sorted(ods.keys())} and {sorted(committed.keys())}")
        return False

    for key in sorted(ods.keys()):
        new = ods.peek(key)
        old = committed.peek(key)
        z_new = np.array(new.peek_xv().dup_array())
        p_new = np.array(new.peek_yv().dup_array())
        z_old = np.array(old.peek_xv().dup_array())
        p_old = np.array(old.peek_yv().dup_array())

        if len(z_new) != len(z_old) or not np.array_equal(p_new == 0.0, p_old == 0.0):
            print(f"{key}: the number of points or the zeros differ")
            equal = False
            continue

        diff = max(np.abs(z_new - z_old).max(), np.abs(p_new - p_old).max())
        largest = max(largest, diff / p_old.max())

    print(f"largest difference: {largest:.3g} of the largest value of a table")
    return equal and largest <= CHECK_TOLERANCE


def main() -> int:
    """Write the fixture, or compare it with the committed one."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="do not write; exit with status 1 if the committed fixture differs",
    )
    args = parser.parse_args()
    ods = build()

    if args.check:
        return 0 if check(ods) else 1

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    ser.dict_str_to_binfile(ods, OUT)
    print(f"{ods.len()} P(z) written to {os.path.normpath(OUT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
