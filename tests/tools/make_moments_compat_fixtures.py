#!/usr/bin/env python
"""Freeze the truncated-series shape factors as references for the exact ones.

Run this against a library that still carries NcGalaxyShapeFactorTiltedSeries
and NcGalaxyShapeFactorMomentSeries, i.e. *before* their removal. They were
the only independent check on NcGalaxyShapeFactorMomentsTilt and
NcGalaxyShapeFactorMomentsGauss at small shear -- for the TRACE_DET
convention, the only check besides the closed-form derivation itself -- and
the values written here keep that check reproducible once they are gone. It
should not be re-run to "refresh" the fixtures: after the removal it cannot
be, and before it, doing so would defeat their purpose.

Both old classes are truncated series in the shear, so they agree with the
exact classes only up to their truncation error, which grows with g. The
tolerance stored for each (class, convention, g) is ten times the largest
gap measured here, so the tests assert "still the same object", not an
equality the two constructions never had.

Writes data/truth_tables/moments/golden.json.
"""

import json
import os

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "truth_tables", "moments"
)
TRUNC_ORDER = 9
CONVS = {"trace": Nc.GalaxyWLObsEllipConv.TRACE, "trace-det": Nc.GalaxyWLObsEllipConv.TRACE_DET}
SIGMA_POPS = (0.2, 0.4658)
STD_NOISES = (0.02, 0.1)
SHEARS = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
EPS_OBS = ((0.05, 0.02), (0.3, -0.1), (-0.4, 0.25), (0.6, 0.3), (0.1, -0.7), (0.85, 0.0))
TOL_FACTOR = 10.0
TOL_FLOOR = 1.0e-8


def build_mset(sigma_pop):
    """Lens models plus the population the marginal reads; mirrored by the tests."""
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    pop = Nc.GalaxyShapePopGauss.new()
    pop["sigma"] = sigma_pop
    hms.param_set_by_name("log10MDelta", 14.0)
    hp.param_set_by_name("z", 0.2)
    hp.prepare(cosmo)
    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop):
        mset.set(model)
    mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())
    mset.set(Nc.GalaxyRedshiftObsGauss.new())
    return mset, pop


def make_data(gsf, mset, std_noise):
    """A single galaxy's data, returned with the fragments it references."""
    posf = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    pos_data = Nc.GalaxyPositionFactorData.new(posf, mset)
    zf = Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0)
    z_data = Nc.GalaxyRedshiftFactorData.new(zf, mset)
    data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)
    gsf.data_set(data, 0.0, 0.0, std_noise, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)
    gsf.prepare_data_array(mset, [data], True, True)
    return data, (pos_data, z_data)


def main():
    """Evaluate every case and write golden.json."""
    pairs = {
        "tilted_series": (
            lambda conv: Nc.GalaxyShapeFactorTiltedSeries.new(conv, TRUNC_ORDER),
            Nc.GalaxyShapeFactorMomentsTilt.new,
        ),
        "moment_series": (
            lambda conv: Nc.GalaxyShapeFactorMomentSeries.new(conv, TRUNC_ORDER),
            Nc.GalaxyShapeFactorMomentsGauss.new,
        ),
    }
    cases = []
    gap = {old: {conv: {str(g): 0.0 for g in SHEARS} for conv in CONVS} for old in pairs}

    for conv_name, conv in CONVS.items():
        for sigma_pop in SIGMA_POPS:
            mset, pop = build_mset(sigma_pop)
            for std_noise in STD_NOISES:
                objs = {}
                for old, (make_old, make_new) in pairs.items():
                    gsf_old, gsf_new = make_old(conv), make_new(conv)
                    objs[old] = (gsf_old, make_data(gsf_old, mset, std_noise),
                                 gsf_new, make_data(gsf_new, mset, std_noise))
                for g in SHEARS:
                    for e1, e2 in EPS_OBS:
                        case = {"conv": conv_name, "sigma_pop": sigma_pop, "std_noise": std_noise,
                                "g": g, "e1": e1, "e2": e2}
                        for old, (gsf_old, (d_old, _), gsf_new, (d_new, _)) in objs.items():
                            v_old = gsf_old.eval_ln_marginal(pop, d_old, g, 0.0, e1, e2)
                            v_new = gsf_new.eval_ln_marginal(pop, d_new, g, 0.0, e1, e2)
                            case[old] = v_old
                            gap[old][conv_name][str(g)] = max(gap[old][conv_name][str(g)], abs(v_old - v_new))
                        cases.append(case)

    tol = {old: {conv: {g: max(TOL_FACTOR * v, TOL_FLOOR) for g, v in per_g.items()}
                 for conv, per_g in per_conv.items()}
           for old, per_conv in gap.items()}

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "golden.json"), "w", encoding="utf-8") as f:
        json.dump({"trunc_order": TRUNC_ORDER, "tol_factor": TOL_FACTOR, "tol_floor": TOL_FLOOR,
                   "measured_gap": gap, "tol": tol, "cases": cases}, f, indent=1)
    print(f"wrote {len(cases)} cases to {os.path.join(OUT, 'golden.json')}")
    for old, per_conv in gap.items():
        for conv, per_g in per_conv.items():
            print(f"  {old:14s} {conv:10s} max gap by g: " + ", ".join(f"{g}:{v:.1e}" for g, v in per_g.items()))


if __name__ == "__main__":
    main()
