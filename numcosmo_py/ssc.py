#
# ssc.py
#
# Mon Aug 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# ssc.py
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

r"""Super-sample covariance $S_{ij}$ matrices computed with NumCosmo only.

This module replaces the `PySSC` code path used so far in `numcosmo_py`: it
needs neither `scipy` nor `healpy`, and it uses the non-Limber `NcXcorSolver`
tier (the "UltraLevin" spherical-Bessel machinery) for every radial integral.

The two quantities `PySSC` provides are:

*Full sky* --- with $C^{ij}_\ell$ the angular power spectrum of the
volume-averaged matter density contrast over redshift bins $i$ and $j$,

$$S_{ij} = \frac{C^{ij}_0}{4\pi} .$$

*Partial sky* --- with $C^{\rm mask}_\ell$ the angular power spectrum of the
survey mask and $f_{\rm sky} = \sqrt{C^{\rm mask}_0 / 4\pi}$,

$$S_{ij} = \frac{1}{(4\pi f_{\rm sky})^2}
    \sum_{\ell=0}^{\ell_{\rm max}} (2\ell+1) C^{\rm mask}_\ell C^{ij}_\ell .$$

Both reduce to the same $C^{ij}_\ell$ engine, which here is
#NcXcorKernelClusterTophat plus #NcXcorSolver in its non-Limber mode.

Long partial-sky runs report progress through the `progress` callback of
:meth:`SijCalculator.compute_cl`, which fires once per solved multipole chunk.

Precision
---------

The off-diagonal $S_{ij}$ of well-separated bins is a small residual of a large
cancellation: $|S_{06}|$ is four orders of magnitude below $S_{00}$ for J-PAS
bins. Its accuracy is set by `scaled_abstol`, the absolute floor of the adaptive
refinement building the $U_i(k)$ spline --- **not** by `reltol`, which is
inactive while `scaled_abstol` binds. #NcXcorKernel defaults to
`scaled_abstol = 1e-4`, and at that value tightening `reltol` changes nothing at
all. Measured against a converged independent quadrature, for seven J-PAS bins:

=====================  =========  =========  =========  =========  =======
`scaled_abstol`        $S_{00}$   $S_{01}$   $S_{04}$   $S_{06}$   time
=====================  =========  =========  =========  =========  =======
`1e-4` (NcXcorKernel)     -0.07%     -0.59%     -19.6%     -59.1%   0.04 s
`1e-6` (used here)        0.006%    -0.014%      0.28%      -5.2%   0.13 s
`1e-8`                    0.001%    -0.003%      0.02%     -0.06%   0.36 s
=====================  =========  =========  =========  =========  =======

Read that table for what it is: a *relative* error on elements that are
themselves four orders of magnitude below $S_{00}$. It is the wrong quantity to
tune against, because those elements contribute to the covariance in absolute
terms, not relative ones, and `-59.1%` of a negligible number is still
negligible. The quantity that matters is what reaches the parameter
uncertainties. Measured on a cap footprint, $z \in [0.1, 0.8]$, seven bins,
18000 deg$^2$, against the `1e-8` result:

=====================  ============  ===============  ==================
`scaled_abstol`        max |dS/S|    max |dS|/S_{00}  max |dsigma/sigma|
=====================  ============  ===============  ==================
`1e-4` (NcXcorKernel)      1.4e-02          3.0e-04              0.028%
`1e-5` (used here)         1.7e-03          1.2e-04              0.009%
`1e-6`                     8.7e-05          1.4e-06              0.0001%
=====================  ============  ===============  ==================

`sigma` is the Fisher uncertainty on $(\Omega_c, w, \ln 10^{10}A_s)$ for the
full forecast likelihood, measured against `1e-7`. Every setting in the table
is converged for anything the forecast reports: the effects it studies are
5--50%, so even `1e-4` has three orders of magnitude of margin. Do not tighten
this chasing the off-diagonal percentages in the first table.

What does move with `scaled_abstol` is cost, and it matters in one place. With
a frozen $S_{ij}$ the matrix is built once per experiment and the choice is
irrelevant (a second either way). With `--vary-fitting-sij` it is rebuilt at
every likelihood step, so it sits directly on the chain's critical path. One
rebuild, cap, 18000 deg$^2$, seven bins:

=====================  ==========  ====================
`scaled_abstol`        S rebuild   relative to `1e-6`
=====================  ==========  ====================
`1e-4`                    0.28 s   3.8x faster
`1e-5` (used here)        0.51 s   2.1x faster
`1e-6`                    1.06 s   ---
`1e-7`                    2.03 s   1.9x slower
=====================  ==========  ====================

Hence `1e-5`: accurate to `0.009%` on the reported uncertainties, and twice as
cheap as `1e-6` in the varying case. The accuracy that `1e-6` and below buy is
real but unusable, and in the varying case it is paid for on every step.

There is a floor under all of this, and `1e-6` is it. `scaled_abstol` is a
fraction of the peak of `W(k)`, but the quantity integrated is `k^2 W_i W_j`,
so it enters *squared*: `1e-6` here is `1e-12` on the integrand, already past
what the outer integral carries, and `1e-8` would be `1e-16`, below double
precision. Below `1e-6` there is nothing left to buy at any price --
`nc_xcor_kernel_set_scaled_abstol` warns there, see
`NC_XCOR_KERNEL_MIN_USEFUL_SCALED_ABSTOL`. The `1e-7` and `1e-8` columns in the
tables above are measurements of that, not options.

Keep `scaled_abstol` away from `reltol`
--------------------------------------

`scaled_abstol` must not equal the outer $k$-integral's `reltol`. The p-adaptive
`g_error` in `ncm_integral_nd_eval` ("p-adaptive methods report failure when
they run out of Clenshaw-Curtis levels") is not monotonic in `scaled_abstol`:
it fires at that one coincidence and nowhere near it. Generating the
configuration that first hit it (15 knots, cap, 3000 deg$^2$, $w = -0.8$, with
`reltol = 1e-6`) at each floor:

=====================  ====  ====  ====  ====  ====  ====  ====
`scaled_abstol`        1e-4  1e-5  3e-6  1e-6  3e-7  1e-7  1e-8
=====================  ====  ====  ====  ====  ====  ====  ====
p-adaptive failures       0     0     0     1     0     0     0
=====================  ====  ====  ====  ====  ====  ====  ====

The refinement stops exactly at the level the outer rule is trying to resolve,
so the non-smoothness it leaves behind sits precisely at `pcubature`'s
convergence threshold; a factor of three either way breaks the tie. This is the
mismatch `_nc_xcor_check_kernel_tolerance` describes, but that guard compares
`reltol` against the *Levin* tolerances, not against `scaled_abstol`, so it
does not catch it.

So the offset from `reltol` is deliberate, and it is taken in the cheap
direction: `DEFAULT_SCALED_ABSTOL = 1e-5` against `DEFAULT_RELTOL = 1e-6`.
Moving the other way, to `1e-7`, would clear the coincidence just as well but
cost 1.9x per rebuild for accuracy that is already unusable. Do not "tidy" the
two constants to the same value.

The C-side counterpart must move with it. `NcXcorSSCSij` carries its own
`NC_XCOR_SSC_SIJ_DEFAULT_SCALED_ABSTOL` (`nc_xcor_ssc_sij.c`) and pushes it onto
the kernels it builds, so it is what the varying path actually uses --- not this
constant. The two are kept equal on purpose: `create_ssc_sij_calculator()`
promises that a fixed and a varying run "differ only in whether $S_{ij}$ follows
the cosmology, not in how it is computed", and comment 15 of the covariance
paper is precisely the comparison between them. Changing one alone silently
breaks that.

The failure is not fatal in any case: it falls back to h-adaptive subdivision,
which is what let the affected run finish. The offset removes the retry rather
than the crash.

`adaptive_epsilon` was verified not to bind at any of these settings, and
`reltol` is inert while `scaled_abstol` binds.
"""

from typing import Callable, Sequence
import time

import numpy as np
from numpy.typing import NDArray

from numcosmo_py import Nc, Ncm

__all__ = [
    "ProgressCallback",
    "SijCalculator",
    "find_lmax",
    "mask_angular_power_spectrum",
]

#: Signature of the progress callback: (done, total, elapsed seconds, message).
ProgressCallback = Callable[[int, int, float, str], None]

#: Default multipole block size handed to `NcXcorSolver.plan_blocks()`. Eight is
#: the empirical sweet spot measured in `dev-notes/xcor_ultralevin_batching_plan.md`
#: section 1.3, not `ell_cache_max`.
DEFAULT_BLOCK_SIZE = 8

#: Relative tolerance for the `U_i(k)` spline and the outer `k` integral. Not the
#: knob that limits accuracy (see `DEFAULT_SCALED_ABSTOL`), and it cannot be
#: tightened much: at `1e-7` the p-adaptive cubature runs out of Clenshaw-Curtis
#: levels on the cross integrand for `l > 0` and aborts.
DEFAULT_RELTOL = 1.0e-6

#: Absolute floor for the adaptive refinement of the `U_i(k)` spline. One order
#: tighter than #NcXcorKernel's own `1e-4` default, and the knob that limits
#: off-diagonal accuracy -- see the module docstring. Deliberately offset from
#: `DEFAULT_RELTOL`: making the two equal is the one setting that trips the
#: p-adaptive cubature, and the offset is taken upwards because tightening costs
#: ~2x per rebuild when `--vary-fitting-sij` puts S on the likelihood's critical
#: path, for accuracy already far below anything a forecast reports.
#:
#: Must be kept equal to `NC_XCOR_SSC_SIJ_DEFAULT_SCALED_ABSTOL`
#: (`nc_xcor_ssc_sij.c`), which is what the varying path actually uses.
DEFAULT_SCALED_ABSTOL = 1.0e-5


def print_progress(done: int, total: int, elapsed: float, message: str) -> None:
    """Report progress on stdout, one line per call.

    Suitable as the `progress` argument of :meth:`SijCalculator.compute_cl`.

    :param done: Units of work completed.
    :param total: Total units of work.
    :param elapsed: Seconds since the computation started.
    :param message: Short description of the step just finished.
    """
    frac = done / total if total > 0 else 1.0
    eta = elapsed * (1.0 - frac) / frac if frac > 0.0 else float("nan")
    print(
        f"[Sij] {100.0 * frac:5.1f}%  {done}/{total}  "
        f"elapsed {elapsed:7.1f}s  eta {eta:7.1f}s  {message}",
        flush=True,
    )


def _nside_from_npix(npix: int) -> int:
    """Recover a HEALPix `nside` from a map length, validating it."""
    nside = int(round(np.sqrt(npix / 12.0)))

    if nside < 1 or 12 * nside * nside != npix:
        raise ValueError(f"{npix} is not a valid HEALPix map length.")

    return nside


def mask_angular_power_spectrum(
    mask: NDArray[np.float64],
    lmax: int | None = None,
    mask2: NDArray[np.float64] | None = None,
    iter_n: int = 3,
) -> tuple[NDArray[np.float64], float]:
    r"""Compute the angular power spectrum of a HEALPix mask and its $f_{\rm sky}$.

    Uses #NcmSphereMap, which reproduces `healpy.anafast` to machine precision
    for the same number of iterations, so no `healpy` dependency is needed.

    :param mask: HEALPix map of the mask, in RING ordering.
    :param lmax: Maximum multipole; defaults to twice the map's `nside`.
    :param mask2: Optional second mask, for a cross mask spectrum.
    :param iter_n: Number of analysis iterations, matching `healpy`'s default of 3.
    :return: A tuple `(cl_mask, fsky)` with `cl_mask` of length `lmax + 1`.
    """
    nside = _nside_from_npix(mask.size)

    # Validated before any set_map(): NcmSphereMap asserts the map length
    # against its own npix, so a mismatch there aborts the process rather than
    # raising.
    if mask2 is not None:
        nside2 = _nside_from_npix(mask2.size)

        if nside2 != nside:
            raise ValueError(
                f"mask and mask2 must have the same nside, got {nside} and "
                f"{nside2}."
            )

    smap = Ncm.SphereMap.new(nside)
    smap.set_map(mask)

    if lmax is None:
        lmax = 2 * nside

    smap.set_lmax(lmax)
    smap.set_iter(iter_n)
    smap.prepare_alm()

    if mask2 is None:
        cl_mask = np.array([smap.get_Cl(ell) for ell in range(lmax + 1)])
    else:
        smap2 = Ncm.SphereMap.new(nside)
        smap2.set_map(mask2)
        smap2.set_lmax(lmax)
        smap2.set_iter(iter_n)
        smap2.prepare_alm()
        cl_mask = np.array(smap.compute_cross_Cl(smap2).dup_array())

    fsky = np.sqrt(cl_mask[0] / (4.0 * np.pi))

    return cl_mask, fsky


def find_lmax(cl_mask: NDArray[np.float64], var_tol: float = 0.05) -> int:
    r"""Find the smallest $\ell_{\rm max}$ capturing the mask variance to `var_tol`.

    The mask variance is $\sum_\ell (2\ell+1) C^{\rm mask}_\ell / 4\pi$; this
    returns the smallest truncation reproducing it to within `var_tol`. Same
    criterion as `PySSC.find_lmax`.

    :param cl_mask: Angular power spectrum of the mask.
    :param var_tol: Relative tolerance on the recovered variance.
    :return: The chosen maximum multipole.
    """
    ell = np.arange(cl_mask.size)
    summand = (2.0 * ell + 1.0) * cl_mask / (4.0 * np.pi)
    var_target = summand.sum()

    if var_target <= 0.0:
        raise ValueError("mask has non-positive total variance.")

    cumulative = np.cumsum(summand)
    within = np.abs(cumulative - var_target) / var_target <= var_tol
    # np.argmax on a boolean array returns the first True, or 0 if none are True.
    if not within.any():
        return int(ell[-1])

    return int(np.argmax(within))


class SijCalculator:
    r"""Computes SSC $S_{ij}$ matrices for a set of top-hat redshift bins.

    Each bin becomes an #NcXcorKernelClusterTophat, whose radial window is
    normalized to unit integral in comoving distance, so $C^{ij}_\ell$ is the
    angular power spectrum of the volume-averaged density contrast --- exactly
    the quantity `PySSC` builds from tabulated kernels with `convention=0`.

    Every kernel is put in permanent non-Limber mode (`l_limber = -1`): the
    Limber approximation is meaningless at the low multipoles that dominate
    $S_{ij}$, and it makes the cross spectrum of two disjoint bins vanish.
    """

    def __init__(
        self,
        z_edges: Sequence[float] | NDArray[np.float64],
        powspec: Ncm.Powspec | None = None,
        dist: Nc.Distance | None = None,
        block_size: int = DEFAULT_BLOCK_SIZE,
        reltol: float = DEFAULT_RELTOL,
        scaled_abstol: float = DEFAULT_SCALED_ABSTOL,
    ) -> None:
        """Build the per-bin kernels.

        :param z_edges: Redshift bin edges; `len(z_edges) - 1` bins are created.
        :param powspec: Linear matter power spectrum; defaults to an
            Eisenstein--Hu transfer function over `k` in `[1e-6, 1e3]` 1/Mpc.
        :param dist: Distance object; defaults to `Nc.Distance.new(3.0)`.
        :param block_size: Multipole block size for #NcXcorSolver.
        :param reltol: Relative tolerance for the kernel spline and the outer
            `k` integral. Cannot go below the integrator's `cheb-reltol`
            (`1e-8` by default), which caps the achievable precision.
        :param scaled_abstol: Absolute floor for the adaptive refinement of the
            `U_i(k)` spline. **This, not `reltol`, is what limits the accuracy
            of the off-diagonal `S_ij`** --- see the module docstring.
        """
        z_edges = np.asarray(z_edges, dtype=np.float64)

        if z_edges.ndim != 1 or z_edges.size < 2:
            raise ValueError(
                "z_edges must be a 1-dimensional array of at least 2 edges."
            )
        if np.any(np.diff(z_edges) <= 0.0):
            raise ValueError("z_edges must be strictly increasing.")

        if dist is None:
            dist = Nc.Distance.new(float(z_edges[-1]) + 1.0)
        if powspec is None:
            powspec = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH.new())
            powspec.require_kmin(1.0e-6)
            powspec.require_kmax(1.0e3)

        self.z_edges = z_edges
        self.dist = dist
        self.powspec = powspec
        self.block_size = block_size
        self.reltol = reltol
        self.scaled_abstol = scaled_abstol
        self._mask_cache: tuple[tuple, tuple[NDArray[np.float64], float]] | None = None

        self.kernels = [
            Nc.XcorKernelClusterTophat(
                dist=dist,
                powspec=powspec,
                z_lower=float(z_lower),
                z_upper=float(z_upper),
                integrator=Ncm.SBesselIntegratorLevin.new(0, block_size),
            )
            for z_lower, z_upper in zip(z_edges[:-1], z_edges[1:])
        ]

        for kernel in self.kernels:
            kernel.set_l_limber(-1)
            kernel.set_reltol(reltol)
            kernel.set_scaled_abstol(scaled_abstol)

        self.xcor = Nc.Xcor.new(dist, powspec, Nc.XcorMethod.KERNEL_CUBATURE)
        self.xcor.set_reltol(reltol)

    @property
    def nbins(self) -> int:
        """Number of redshift bins."""
        return len(self.kernels)

    def prepare(self, cosmo: Nc.HICosmo) -> None:
        """Prepare every kernel and the #NcXcor object for `cosmo`.

        :param cosmo: The cosmological model.
        """
        self.dist.prepare_if_needed(cosmo)
        self.powspec.prepare_if_needed(cosmo)

        for kernel in self.kernels:
            kernel.prepare(cosmo)

        self.xcor.prepare(cosmo)

    def compute_cl(
        self,
        cosmo: Nc.HICosmo,
        lmin: int = 0,
        lmax: int = 0,
        progress: ProgressCallback | None = None,
        chunk_size: int | None = None,
    ) -> NDArray[np.float64]:
        r"""Compute $C^{ij}_\ell$ for every bin pair over `[lmin, lmax]`.

        The multipole range is solved in chunks so that progress can be reported
        while the computation runs. Chunking costs nothing: distinct multipole
        blocks share no reusable operator state anyway.

        :param cosmo: The cosmological model.
        :param lmin: Lowest multipole.
        :param lmax: Highest multipole, included.
        :param progress: Optional callback invoked after each chunk, receiving
            `(done, total, elapsed, message)`. :func:`print_progress` is a ready
            implementation.
        :param chunk_size: Multipoles solved per progress step; defaults to four
            solver blocks.
        :return: Array of shape `(nbins, nbins, lmax - lmin + 1)`.
        """
        if lmax < lmin:
            raise ValueError(f"lmax ({lmax}) must not be smaller than lmin ({lmin}).")

        if chunk_size is None:
            chunk_size = 4 * self.block_size

        self.prepare(cosmo)

        nbins = self.nbins
        nell = lmax - lmin + 1
        cl = np.zeros((nbins, nbins, nell))
        pairs = [(i, j) for i in range(nbins) for j in range(i, nbins)]

        chunks = [
            (lo, min(lo + chunk_size - 1, lmax))
            for lo in range(lmin, lmax + 1, chunk_size)
        ]
        start = time.monotonic()

        for n_done, (chunk_lmin, chunk_lmax) in enumerate(chunks, start=1):
            solver = Nc.XcorSolver.new()
            ids = [solver.register_kernel(kernel) for kernel in self.kernels]

            for i, j in pairs:
                solver.request_cl(ids[i], ids[j], chunk_lmin, chunk_lmax)

            solver.plan_blocks(self.block_size)
            solver.solve(self.xcor, cosmo)

            offset = chunk_lmin - lmin
            for request, (i, j) in enumerate(pairs):
                values = np.array(solver.get_result(request).dup_array())
                cl[i, j, offset : offset + values.size] = values
                cl[j, i, offset : offset + values.size] = values

            if progress is not None:
                progress(
                    n_done,
                    len(chunks),
                    time.monotonic() - start,
                    f"ell {chunk_lmin}-{chunk_lmax}",
                )

        return cl

    def fullsky(self, cosmo: Nc.HICosmo) -> NDArray[np.float64]:
        r"""Compute the full-sky $S_{ij} = C^{ij}_0 / 4\pi$.

        :param cosmo: The cosmological model.
        :return: The $S_{ij}$ matrix, shape `(nbins, nbins)`.
        """
        return self.compute_cl(cosmo, 0, 0)[:, :, 0] / (4.0 * np.pi)

    def fullsky_fsky(self, cosmo: Nc.HICosmo, area: float) -> NDArray[np.float64]:
        r"""Compute the full-sky $S_{ij}$ rescaled by $1/f_{\rm sky}$.

        This is the crude finite-area correction, not a real mask deconvolution;
        use :meth:`partial_sky` when the footprint's shape matters.

        :param cosmo: The cosmological model.
        :param area: Survey area in square degrees.
        :return: The rescaled $S_{ij}$ matrix.
        """
        sqd_fullsky = 4.0 * np.pi * (180.0 / np.pi) ** 2

        if not 0.0 < area < sqd_fullsky:
            raise ValueError(
                f"area must lie in (0, {sqd_fullsky:.1f}) square degrees, got {area}."
            )

        fsky = area * (np.pi / 180.0) ** 2 / (4.0 * np.pi)

        return self.fullsky(cosmo) / fsky

    def _mask_spectrum(
        self,
        mask: NDArray[np.float64],
        mask2: NDArray[np.float64] | None,
    ) -> tuple[NDArray[np.float64], float]:
        """Return the mask spectrum, reusing the cached one when unchanged.

        The mask spectrum does not depend on cosmology, so a chain varying
        $S_{ij}$ must not pay for a spherical harmonic transform per step.

        :param mask: HEALPix map of the mask.
        :param mask2: Optional second mask.
        :return: A tuple `(cl_mask, fsky)`.
        """
        key = (
            mask.shape,
            mask.tobytes(),
            None if mask2 is None else mask2.tobytes(),
        )

        if self._mask_cache is None or self._mask_cache[0] != key:
            self._mask_cache = (key, mask_angular_power_spectrum(mask, mask2=mask2))

        return self._mask_cache[1]

    def partial_sky(
        self,
        cosmo: Nc.HICosmo,
        mask: NDArray[np.float64],
        var_tol: float = 0.05,
        lmax: int | None = None,
        mask2: NDArray[np.float64] | None = None,
        progress: ProgressCallback | None = None,
    ) -> NDArray[np.float64]:
        r"""Compute the partial-sky $S_{ij}$ for a HEALPix footprint.

        :param cosmo: The cosmological model.
        :param mask: HEALPix map of the mask, in RING ordering.
        :param var_tol: Relative mask-variance tolerance used to pick
            $\ell_{\rm max}$ when `lmax` is not given.
        :param lmax: Explicit maximum multipole, overriding `var_tol`.
        :param mask2: Optional second mask, for a cross mask spectrum.
        :param progress: Optional progress callback, see :meth:`compute_cl`.
        :return: The $S_{ij}$ matrix, shape `(nbins, nbins)`.
        """
        cl_mask, fsky = self._mask_spectrum(mask, mask2)

        if lmax is None:
            lmax = find_lmax(cl_mask, var_tol)

        if lmax >= cl_mask.size:
            raise ValueError(
                f"lmax ({lmax}) exceeds the mask spectrum's range ({cl_mask.size - 1})."
            )

        cl_mask = cl_mask[: lmax + 1]
        ell = np.arange(lmax + 1)

        cl = self.compute_cl(cosmo, 0, lmax, progress=progress)
        weights = (2.0 * ell + 1.0) * cl_mask

        sij = np.tensordot(cl, weights, axes=([2], [0])) / (4.0 * np.pi * fsky) ** 2

        # np.tensordot is typed as returning floating[Any] whatever it is given,
        # so the dtype has to be restated rather than inferred. Both operands
        # are float64 already, making this a no-op at runtime.
        return np.asarray(sij, dtype=np.float64)
