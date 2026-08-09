"""Step 9: real NcXcorSolver vs CCL -- register all N kernels + all pairs,
plan_blocks(8), solve().

Same cosmology/kernel setup as the v2 prototype it replaced, for direct
comparability with the plan doc's §9.7/§9.8 tables.

Two things this reports that the original version of this script did not, both
because collapsing them lost information the plan doc then repeated:

  * timing split into cold (first solve, which builds the per-block
    factorizations) and warm (every later solve). Only the warm number
    describes an MCMC, where solve() runs once per likelihood evaluation.
  * deviation from CCL resolved per ell, and split at ell=8. A single max over
    ell=[2,25] is set entirely by ell=2 and says nothing about the rest.

`--tuned` (reltol=1e-2 / cheb-reltol=1e-2) is now refused by
_nc_xcor_check_kernel_tolerance(): it asks the outer k-integral for more
precision than that closure carries. §9.8 had already concluded those settings
should not be used; the check enforces it.
"""

import itertools
import os
import sys
import time

import numpy as np
import pyccl

from numcosmo_py import Nc, Ncm
from numcosmo_py.ccl.nc_ccl import create_nc_obj, CCLParams

Ncm.cfg_init()
CCLParams.set_default_params()

ccl_cosmo = pyccl.Cosmology(
    Omega_c=0.25,
    Omega_b=0.05,
    Neff=3.046,
    h=0.7,
    sigma8=0.9,
    n_s=0.96,
    Omega_k=0.0,
    w0=-1.0,
    wa=0.0,
    m_nu=[0.0, 0.0, 0.0],
    transfer_function="eisenstein_hu",
    matter_power_spectrum="linear",
)
cosmology = create_nc_obj(ccl_cosmo, dist_z_max=15.0)
cosmo, dist, ps_ml = cosmology.cosmo, cosmology.dist, cosmology.ps_ml


def _arg_lmax(default=25):
    for a in sys.argv[1:]:
        if a.startswith("--lmax="):
            return int(a.split("=", 1)[1])
    return default


# ell=[2,25] is 3 blocks at block size 8, so it cannot exercise more than 3
# threads (§9.9). Pass --lmax=200 or --lmax=400 for a range with enough blocks
# to say anything about threading.
lmin, lmax = 2, _arg_lmax()
ells = np.arange(lmin, lmax + 1, dtype=float)

# First ten unchanged, so N=5/N=10 stay comparable with the plan doc's tables.
mus_pool = [
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    1.1,
    1.2,
    1.3,
    1.4,
    1.5,
    1.6,
    1.7,
    1.8,
    1.9,
    2.0,
    2.1,
]


def build_kernels(n, mus, sigma=0.08, z_len=1000, tuned=True):
    kernels, dndz_list = [], []
    for mu in mus[:n]:
        z_a = np.linspace(max(mu - 8 * sigma, 0.0), mu + 8 * sigma, z_len)
        nz_a = np.exp(-((z_a - mu) ** 2) / sigma**2 / 2.0) / np.sqrt(
            2.0 * np.pi * sigma**2
        )
        dndz_list.append((z_a, nz_a))
        integrator = Ncm.SBesselIntegratorLevin.new(lmin, lmax)
        if tuned:
            integrator.set_reltol(1.0e-2)
            integrator.set_cheb_reltol(1.0e-2)
            integrator.set_max_order(64)
        k = Nc.XcorKernelWeakLensing(
            dist=dist,
            powspec=ps_ml,
            dndz=Ncm.SplineCubicNotaknot.new_full(
                Ncm.Vector.new_array(z_a), Ncm.Vector.new_array(nz_a), True
            ),
            nbar=3.0,
            intr_shear=7.0,
            integrator=integrator,
        )
        k.set_l_limber(-1)
        k.prepare(cosmo)
        kernels.append(k)
    return kernels, dndz_list


def run(n, mus, tuned, reps=4):
    nc_kernels, dndz_list = build_kernels(n, mus, tuned=tuned)

    xc = Nc.Xcor.new(dist, ps_ml, Nc.XcorMethod.KERNEL_CUBATURE)
    xc.prepare(cosmo)

    solver = Nc.XcorSolver.new()
    kernel_ids = [solver.register_kernel(k) for k in nc_kernels]
    pairs = list(itertools.combinations_with_replacement(range(n), 2))

    for i, j in pairs:
        solver.request_cl(kernel_ids[i], kernel_ids[j], lmin, lmax)

    solver.plan_blocks(8)

    nc_t = []
    for _ in range(reps):
        t0 = time.perf_counter()
        solver.solve(xc, cosmo)
        nc_t.append(time.perf_counter() - t0)

    nc_pairs = {
        (i, j): np.array(solver.get_result(r).dup_array())
        for r, (i, j) in enumerate(pairs)
    }

    ccl_tracers = [pyccl.WeakLensingTracer(ccl_cosmo, dndz=dz) for dz in dndz_list]
    psp = ccl_cosmo.get_linear_power()
    ccl_t = []

    for _ in range(reps):
        t0 = time.perf_counter()
        ccl_pairs = {
            (i, j): pyccl.angular_cl(
                ccl_cosmo,
                ccl_tracers[i],
                ccl_tracers[j],
                ells,
                l_limber=int(lmax) + 1,
                p_of_k_a=psp,
                p_of_k_a_lin=psp,
            )
            for i, j in pairs
        }
        ccl_t.append(time.perf_counter() - t0)

    n_pairs = len(pairs)
    n_threads = os.environ.get("OMP_NUM_THREADS", "?")
    nc_warm, ccl_warm = min(nc_t[1:]), min(ccl_t[1:])

    print(
        f"N={n} kernels, {n_pairs} pairs, ell=[{lmin},{lmax}], tuned={tuned}, "
        f"OMP_NUM_THREADS={n_threads}, OMP_THREAD_LIMIT={os.environ.get('OMP_THREAD_LIMIT', '?')}"
    )
    print(
        f"  NumCosmo solve(): cold {nc_t[0]*1e3:8.1f} ms | warm {nc_warm*1e3:8.1f} ms"
    )
    print(
        f"  CCL angular_cl  : cold {ccl_t[0]*1e3:8.1f} ms | warm {ccl_warm*1e3:8.1f} ms"
    )
    print(
        f"  ratio (NumCosmo/CCL): cold {nc_t[0]/ccl_t[0]:5.2f}x | warm {nc_warm/ccl_warm:5.2f}x"
    )

    dev = np.array(
        [
            [abs(nc_pairs[p][e] / ccl_pairs[p][e] - 1) for p in pairs]
            for e in range(len(ells))
        ]
    )
    print("  deviation vs CCL by ell (max over pairs):")

    for e, l in enumerate(ells):
        if l <= 10 or l % 5 == 0:
            print(f"    ell={int(l):3d}  {dev[e].max():.3e}")

    lo = ells < 8
    print(f"  ell <  8: max {dev[lo].max():.3e}  median {np.median(dev[lo]):.3e}")
    print(f"  ell >= 8: max {dev[~lo].max():.3e}  median {np.median(dev[~lo]):.3e}")

    return nc_t, ccl_t


def sanity_check_solve_vs_compute(n, mus, tuned):
    """Compare solve() output against direct nc_xcor_compute() for the
    same kernels, bypassing CCL entirely -- isolates whether a
    solve()-vs-compute() discrepancy explains the poor CCL agreement."""
    nc_kernels, _ = build_kernels(n, mus, tuned=tuned)

    xc = Nc.Xcor.new(dist, ps_ml, Nc.XcorMethod.KERNEL_CUBATURE)
    xc.prepare(cosmo)

    solver = Nc.XcorSolver.new()
    kernel_ids = [solver.register_kernel(k) for k in nc_kernels]
    pairs = list(itertools.combinations_with_replacement(range(n), 2))

    for i, j in pairs:
        solver.request_cl(kernel_ids[i], kernel_ids[j], lmin, lmax)

    solver.plan_blocks(8)
    solver.solve(xc, cosmo)

    max_dev = 0.0
    for r, (i, j) in enumerate(pairs):
        expected = Ncm.Vector.new(lmax - lmin + 1)
        xc.compute(nc_kernels[i], nc_kernels[j], cosmo, lmin, lmax, expected)
        got = np.array(solver.get_result(r).dup_array())
        exp = np.array(expected.dup_array())
        dev = np.max(np.abs(got / exp - 1))
        max_dev = max(max_dev, dev)
        print(f"  pair ({i},{j}): max dev solve() vs compute() = {dev:.4e}")

    print(f"solve() vs compute() max dev over all pairs: {max_dev:.4e}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 5
    tuned = "--tuned" in sys.argv
    if "--sanity" in sys.argv:
        sanity_check_solve_vs_compute(n, mus_pool, tuned)
    else:
        run(n, mus_pool, tuned)
