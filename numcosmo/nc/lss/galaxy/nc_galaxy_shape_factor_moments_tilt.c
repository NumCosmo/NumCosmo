/***************************************************************************
 *            nc_galaxy_shape_factor_moments_tilt.c
 *
 *  Mon Sep 22 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_moments_tilt.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyShapeFactorMomentsTilt:
 *
 * Exponential tilt of the exact zero-shear marginal, solved exactly and
 * interpolated on a dyadic Chebyshev mesh in $\hat g=\min(g,1/g)$.
 *
 * The density is the exact zero-shear marginal tilted by a quadratic
 * statistic,
 * $$
 *   \ln P(\chi_\mathrm{obs}\mid g) = \ln P_0(|\chi_\mathrm{obs}|)
 *     + \lambda_1 x + \lambda_2 x^2 + \lambda_3 y^2 - W(\lambda),
 * $$
 * with $(x,y)$ the observed ellipticity rotated by $-\arg g$ and
 * $T=(x,x^2,y^2)$ the sufficient statistic. $\lambda$ is fixed by matching
 * the model's moments of $T$ to the exact moments $t(\hat g)$ of the
 * lensed, noisy marginal, $\nabla W(\lambda) = t$, and $W$ is the exact
 * log-normalizer, so the model integrates to one by construction.
 *
 * The target moments are closed form in both ellipticity conventions. Under
 * TRACE the map depends on $g$ only through $\delta = 2g/(1+g^2)$ and the
 * circle averages follow from the Moebius lemma; under TRACE_DET the map is
 * a disc automorphism, so $\mathrm{E}[x] = \hat g$ exactly and the second
 * moments reduce to one radial average. #NcGalaxyShapeFactorMomentsGauss
 * uses the same moments for a matched Gaussian.
 *
 * The moment equations are solved exactly, by a damped Newton iteration at
 * Chebyshev-Lobatto nodes, and the solution is interpolated. Solving them
 * order by order in $\delta$ instead gives a Taylor section of the inverse
 * map whose radius of convergence is $\hat g\simeq 0.62$, so it is
 * inaccurate at large shear.
 *
 * ## Accuracy against the exact marginal
 *
 * The tilt is an approximation to the exact lensed, noisy marginal: it
 * matches three moments, and #NcGalaxyShapeFactorMomentsTilt:accuracy-gate
 * bounds only the interpolation of the solved tilt, not its distance from
 * the exact marginal. At $g = 0$ the model is the exact zero-shear
 * marginal. Measured on the HSC-SSP PDR1 catalog HWL16a-002 (TRACE
 * convention, Gaussian population with $\sigma = 0.3$, the catalog's own
 * ellipticities and noise), against #NcGalaxyShapeFactorQuad at $g = 0.03$
 * on 188 galaxies with $|\chi_\mathrm{obs}| < 1$: $|\Delta\ln P|$ up to
 * $0.26$, median $0.016$. On the 14161 galaxies of that catalog with
 * $|\chi_\mathrm{obs}| < 1$, with the same redshift nodes, the best-fit
 * $\log_{10} M$ is $0.030$ dex below the one from
 * #NcGalaxyShapeFactorSeriesLensed at order 12 (which agrees with
 * #NcGalaxyShapeFactorQuad to $6\times10^{-5}$ per galaxy), $0.32$ of its
 * standard error. For an exact marginal use #NcGalaxyShapeFactorSeriesLensed
 * at order 12 or higher (Gauss and Beta populations) or
 * #NcGalaxyShapeFactorQuad ($|\chi_\mathrm{obs}| < 1$ only).
 *
 * ## The mesh
 *
 * $\lambda_1$ has a boundary layer at $\hat g=1$ whose width scales as
 * $\sigma_\nu$ and whose height scales as $1/(2\sigma_\nu^2)$. A single
 * Chebyshev series on $[0,1]$ equioscillates in absolute error over the
 * whole interval, so it needs an order growing like $2/\sigma_\nu$, and
 * below $\sigma_\nu\simeq 0.015$ no order reaches the accuracy gate at all:
 * past order 64 the construction is limited by the conditioning of the
 * solve rather than by resolution.
 *
 * The representation used here is a dyadic mesh of panels
 * $p_j = 1 - 2^{-j}$, $j = 0,\dots,K-1$, each carrying its own Chebyshev
 * series, with
 * $$
 *   K = \max\left(2, \mathrm{round}\left(\log_2(1/\sigma_\nu)\right)\right).
 * $$
 * The layer is resolved once $2^{-(K-1)}\lesssim\sigma_\nu$, and the error
 * does not decrease for larger $K$, so $K$ is fixed by this rule and is not
 * a parameter to tune. The worst $|\ln Z|$ over three population widths is
 * $9.8\times10^{-4}$ at $\sigma_\nu=0.008$, $5.3\times10^{-4}$ at $0.015$,
 * $3.8\times10^{-5}$ at $0.05$ and $3.0\times10^{-6}$ at $0.40$, inside
 * the #NcGalaxyShapeFactorMomentsTilt:accuracy-gate at every width. The
 * mesh needs fewer Newton solves than a single series, and its evaluation
 * degree is five to eight instead of sixteen to sixty-four.
 *
 * The panel index costs one frexp() and no search, because panel $j$ covers
 * $1-\hat g\in(2^{-(j+1)},2^{-j}]$, i.e. $j$ is the exponent field of the
 * double $1-\hat g$. This is why the dyadic mesh is used instead of an
 * adaptively placed one with fewer panels.
 *
 * ## Local degree
 *
 * Each panel chooses its own degree at build time. Chebyshev-Lobatto nodes
 * nest exactly under doubling (the degree-$d$ nodes are the even-indexed
 * degree-$2d$ nodes), so refining from $d$ to $2d$ adds $d$ solves and
 * discards none. The stopping test is the signed combination
 * $$
 *   \Xi(m) = \sum_{k>m}\left| t_1 a_k + t_2 b_k + t_3 c_k - e_k \right|,
 * $$
 * where $a,b,c,e$ are the Chebyshev coefficients of
 * $\lambda_1,\lambda_2,\lambda_3,W$ and $t$ is the target moment vector at
 * the panel midpoint. This is the linear combination the evaluation forms,
 * and it determines the model-mass error. The termwise bound
 * $\sum(|t_1a_k|+|t_2b_k|+\dots)$ is a rigorous upper bound on it but
 * overestimates by a median factor of 45, which would select degrees far
 * higher than needed. $\Xi$ is accurate to a factor of about three and is
 * not an upper bound, so the test is run against the gate divided by ten.
 *
 * ## Error reporting
 *
 * A shear above the table's upper end is outside the range the table was
 * built for. A Chebyshev panel diverges outside its own interval, so the
 * evaluation does not extrapolate: it is counted, warned about once per
 * instance, and returns zero probability, which #NcDataClusterWLFactor
 * routes through its existing low-probability path. Newton failures during
 * a build are counted separately and read back with
 * nc_galaxy_shape_factor_moments_tilt_get_solve_error_count(); the count
 * is per distinct table, not per galaxy, since one table serves every
 * galaxy sharing its key.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_moments_tilt.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_moments_private.h"
#include "nc/background/nc_hicosmo.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/halo/nc_halo_mass_summary.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "ncm/core/ncm_obj_array.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_memory_pool.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Radial nodes of the source-disc quadrature inside the Newton solve. The
 * source disc has no noise edge layer, so 32 nodes resolve it. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_SRC_NR 32

#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NPHI_MIN 16

/* Largest angular node count of the Newton solve. The verification pass
 * re-evaluates at twice the working count and re-solves if the residual is
 * not acceptable there; at the cap that pass cannot run, so under-resolution
 * is not detected there. The largest lambda_1 of the small-noise end needs
 * 2048 nodes, so the cap is 4096. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NPHI_MAX 4096

/* Newton controls.
 *
 * TOL is the residual the iteration stops at; OK is the residual it accepts.
 * The moments are quadratures of exp(lambda . T), and the residual their
 * cancellation can resolve has a floor between 1.3e-9 and 2.4e-9, flat in
 * |lambda| over the range the small-noise end reaches (1.3e-9 at
 * |lambda| = 284, 1.8e-9 at |lambda| = 7900), so the acceptance is a
 * constant above the floor and not a function of the iterate. A residual
 * of 1e-8 on moments of order one is five orders below the accuracy gate.
 *
 * An acceptance at the floor turns converged solves into reported failures,
 * and a failed node is replaced by its predecessor's lambda, which is a
 * wrong value at that node.
 */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_TOL 1.0e-12
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_OK 1.0e-8

#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_MAXIT 40
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_MAXHALF 40

/* Armijo constant of the line search on phi = W - lambda.t. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_ARMIJO_C1 1.0e-4

/* ln P_0 panel: 64 nodes on a window localized around the mass. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_LNP0_NNODES 64
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_LNP0_DELTA 40.0

/* Most bytes of a galaxy's table data_prefetch() fetches. A restricted
 * table is below this and is fetched whole; a full-mesh one is cut at this
 * size, which leaves out its last panels, the ones the fewest evaluations
 * reach. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_PREFETCH_CAP 4096

/* Smallest panel count of the mesh. The count is K = round(log2(1/sigma_nu))
 * bounded below by this value and above by
 * NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS (8); K exceeds the cap for
 * sigma_nu < 2^-8.5 ~ 0.0028, and such a sigma_nu is clamped to the cap
 * and warned about once per instance. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_MIN_PANELS 2

/* Margin of the stopping test against the accuracy gate. The signed tail
 * tracks the true error to a factor of about three and is not an upper
 * bound on it, so the test is run against gate / SAFETY. */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_GATE_SAFETY 10.0

/*
 * ---- Gauss-Legendre node tables on [-1,1], built once ----
 *
 * Both node counts are compile-time constants, so the tables depend on
 * nothing per galaxy. gsl_integration_glfixed_table_alloc() allocates and
 * must not be called per evaluation.
 */
static gdouble _moments_tilt_gl32_x[NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_SRC_NR];
static gdouble _moments_tilt_gl32_w[NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_SRC_NR];

static void
_moments_tilt_gl_init (void)
{
  static gsize init = 0;

  if (g_once_init_enter (&init))
  {
    _nc_galaxy_shape_factor_moments_gl_fill (_moments_tilt_gl32_x, _moments_tilt_gl32_w,
                                             NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_SRC_NR);
    g_once_init_leave (&init, 1);
  }
}

/*
 * ---- The analytic-noise source-disc quadrature ----
 *
 * Computes Z, E[T] and Cov[T] under q = P_0 e^{lam.T}/Z, with the noise
 * integral done in closed form (this file's class doc comment), and writes
 * E[T] to @E, Cov[T] to @C and ln Z to @lnZ. @r/@wr are the source-disc
 * radial nodes and their weights with P_pop folded in (so sum wr = 1 and
 * lam = 0 gives Z = 1 exactly); @cphi/@s2phi are the angular tables.
 * Returns FALSE outside the natural domain u,v > 0 and when Z is not finite
 * and positive.
 *
 * The exponent is evaluated twice rather than stored: for fixed r it is
 * the quadratic K0 + K1 c + K2 c^2 in c = cos(phi), three flops, while
 * storing it would need an nr x nphi scratch array that grows with the
 * angular refinement.
 *
 * @C may be NULL, in which case the fourth-order sums (E[x^3], E[x^4],
 * E[x^2 y^2], E[y^4]) are skipped; they are more than half the inner
 * loop's arithmetic. The line search needs only E[T] to test a trial step,
 * so it passes NULL.
 */
static gboolean
_moments_tilt_moments (const gdouble *r, const gdouble *wr, const guint nr,
                       const gdouble *cphi, const gdouble *s2phi, const guint nphi,
                       const gdouble sn2, const gdouble *lam,
                       gdouble *E, gdouble *C, gdouble *lnZ)
{
  const gboolean want_cov = (C != NULL);
  const gdouble l1        = lam[0];
  const gdouble l2        = lam[1];
  const gdouble l3        = lam[2];
  const gdouble u         = 1.0 - 2.0 * l2 * sn2;
  const gdouble v         = 1.0 - 2.0 * l3 * sn2;
  gdouble t2, w2, m0, l1_u, l2_u, l3_v, wphi;
  gdouble smax = -G_MAXDOUBLE;
  gdouble Z    = 0.0;
  gdouble s_x  = 0.0, s_x2 = 0.0, s_y2 = 0.0;
  gdouble s_x3 = 0.0, s_x4 = 0.0, s_xy2 = 0.0, s_x2y2 = 0.0, s_y4 = 0.0;
  guint i, j;

  if ((u <= 0.0) || (v <= 0.0))
    return FALSE;

  t2   = sn2 / u;
  w2   = sn2 / v;
  m0   = l1 * sn2 / u;
  l1_u = l1 / u;
  l2_u = l2 / u;
  l3_v = l3 / v;
  wphi = 1.0 / nphi;

  /* Pass 1: the exponent's maximum over the grid. Taken from the grid
   * itself rather than from a continuous bound, so the largest term is
   * exactly 1 at every angular resolution. */
  for (i = 0; i < nr; i++)
  {
    const gdouble ri = r[i];
    const gdouble r2 = ri * ri;
    const gdouble K0 = l3_v * r2;
    const gdouble K1 = l1_u * ri;
    const gdouble K2 = l2_u * r2 - K0;

    for (j = 0; j < nphi; j++)
    {
      const gdouble c   = cphi[j];
      const gdouble lnM = (K2 * c + K1) * c + K0;

      if (lnM > smax)
        smax = lnM;
    }
  }

  /* Pass 2: accumulate. */
  for (i = 0; i < nr; i++)
  {
    const gdouble ri = r[i];
    const gdouble r2 = ri * ri;
    const gdouble K0 = l3_v * r2;
    const gdouble K1 = l1_u * ri;
    const gdouble K2 = l2_u * r2 - K0;
    const gdouble mr = ri / u;
    const gdouble pr = r2 / (v * v);
    const gdouble Wr = wr[i] * wphi;

    for (j = 0; j < nphi; j++)
    {
      const gdouble c   = cphi[j];
      const gdouble lnM = (K2 * c + K1) * c + K0;
      const gdouble Wg  = Wr * exp (lnM - smax);
      const gdouble m   = mr * c + m0;
      const gdouble m2  = m * m;
      const gdouble p2  = pr * s2phi[j];
      const gdouble Ex  = m;
      const gdouble Ex2 = m2 + t2;
      const gdouble Ey2 = p2 + w2;

      Z    += Wg;
      s_x  += Wg * Ex;
      s_x2 += Wg * Ex2;
      s_y2 += Wg * Ey2;

      if (want_cov)
      {
        const gdouble Ex3 = m * (m2 + 3.0 * t2);
        const gdouble Ex4 = m2 * (m2 + 6.0 * t2) + 3.0 * t2 * t2;
        const gdouble Ey4 = p2 * (p2 + 6.0 * w2) + 3.0 * w2 * w2;

        s_x3   += Wg * Ex3;
        s_x4   += Wg * Ex4;
        s_xy2  += Wg * Ex * Ey2;
        s_x2y2 += Wg * Ex2 * Ey2;
        s_y4   += Wg * Ey4;
      }
    }
  }

  if (!(Z > 0.0) || !gsl_finite (Z))
    return FALSE;

  E[0] = s_x / Z;
  E[1] = s_x2 / Z;
  E[2] = s_y2 / Z;

  if (want_cov)
  {
    C[0] = s_x2 / Z - E[0] * E[0];
    C[1] = s_x3 / Z - E[0] * E[1];
    C[2] = s_xy2 / Z - E[0] * E[2];
    C[3] = C[1];
    C[4] = s_x4 / Z - E[1] * E[1];
    C[5] = s_x2y2 / Z - E[1] * E[2];
    C[6] = C[2];
    C[7] = C[5];
    C[8] = s_y4 / Z - E[2] * E[2];
  }

  /* The r-independent terms are added after the quadrature, not folded
   * into the exponent. A term folded in and also added here would be a
   * constant shift of W, hence of ln P, indistinguishable from model
   * error. */
  *lnZ = log (Z) + smax + 0.5 * l1 * l1 * sn2 / u - 0.5 * log (u * v);

  return TRUE;
}

/* Solves the 3x3 system A x = b by Gaussian elimination with partial
 * pivoting; returns FALSE on a zero pivot. Cov(T) is positive definite in
 * exact arithmetic, but at a Newton iterate far from the solution the
 * quadrature's version of it can be indefinite, where a Cholesky
 * factorization fails. */
static gboolean
_moments_tilt_solve3 (const gdouble *A_in, const gdouble *b_in, gdouble *x)
{
  gdouble A[9];
  gdouble b[3];
  guint i, j, k;

  memcpy (A, A_in, 9 * sizeof (gdouble));
  memcpy (b, b_in, 3 * sizeof (gdouble));

  for (k = 0; k < 3; k++)
  {
    guint piv     = k;
    gdouble pivot = fabs (A[k * 3 + k]);

    for (i = k + 1; i < 3; i++)
    {
      if (fabs (A[i * 3 + k]) > pivot)
      {
        pivot = fabs (A[i * 3 + k]);
        piv   = i;
      }
    }

    if (!(pivot > 0.0))
      return FALSE;

    if (piv != k)
    {
      for (j = 0; j < 3; j++)
      {
        const gdouble tmp = A[k * 3 + j];

        A[k * 3 + j]   = A[piv * 3 + j];
        A[piv * 3 + j] = tmp;
      }

      {
        const gdouble tmp = b[k];

        b[k]   = b[piv];
        b[piv] = tmp;
      }
    }

    for (i = k + 1; i < 3; i++)
    {
      const gdouble f = A[i * 3 + k] / A[k * 3 + k];

      for (j = k; j < 3; j++)
        A[i * 3 + j] -= f * A[k * 3 + j];

      b[i] -= f * b[k];
    }
  }

  for (i = 3; i-- > 0;)
  {
    gdouble s = b[i];

    for (j = i + 1; j < 3; j++)
      s -= A[i * 3 + j] * x[j];

    x[i] = s / A[i * 3 + i];
  }

  return TRUE;
}

static void
_moments_tilt_set_phi (gdouble *cphi, gdouble *s2phi, const guint nphi)
{
  guint j;

  for (j = 0; j < nphi; j++)
  {
    const gdouble phi = 2.0 * M_PI * j / nphi;
    const gdouble s   = sin (phi);

    cphi[j]  = cos (phi);
    s2phi[j] = s * s;
  }
}

/*
 * Angular node count of the Newton solve at @lam: the smallest power of two
 * in [NPHI_MIN, NPHI_MAX] not below k sqrt(A), A = |lam_1|/u. The measured
 * constant needed runs from 10.5 to 26.2 across ghat and the catalog's
 * noise corners. Convergence in n_phi is spectral, so under-resolution has
 * no symptom in the residual (8 sqrt(A) loses eight digits); the solve
 * therefore also verifies its answer at 2 n_phi.
 */
static guint
_moments_tilt_nphi_for (const gdouble *lam, const gdouble sn2)
{
  const gdouble u = MAX (1.0 - 2.0 * lam[1] * sn2, 1.0e-6);
  const gdouble A = fabs (lam[0]) / u;

  /* k = 16 halves the build and holds lambda to 4e-5 for sigma_nu >= 0.02;
   * below that it loses four orders against k = 32, so the small-noise end
   * uses the larger constant. */
  const gdouble k = (sn2 >= 4.0e-4) ? 16.0 : 32.0;
  const gdouble n = k * sqrt (MAX (A, 1.0));
  guint nphi      = NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NPHI_MIN;

  while ((nphi < (guint) ceil (n)) && (nphi < NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NPHI_MAX))
    nphi *= 2;

  return nphi;
}

/*
 * Damped Newton for lam with E_lam[T] = @target, started from @lam; writes
 * the solution to @lam and W(lam) to @W. Returns TRUE if the final residual
 * is below NEWTON_OK, FALSE if a moment evaluation or the 3x3 solve fails
 * or the residual is above NEWTON_OK at the iteration cap; on FALSE at the
 * cap @lam and @W still hold a matching pair. @cphi/@s2phi are scratch of
 * NPHI_MAX entries whose filled length is tracked by @nphi_io across calls,
 * so the continuation along ghat rebuilds the angular tables only when the
 * resolution rule changes them.
 */
static gboolean
_moments_tilt_newton (const gdouble *r, const gdouble *wr, const guint nr,
                      gdouble *cphi, gdouble *s2phi, guint *nphi_io,
                      const gdouble sn2, const gdouble *target, const guint maxit,
                      gdouble *lam, gdouble *W)
{
  gdouble E[3], C[9], F[3], step[3], trial[3], E2[3];
  gdouble lnZ  = 0.0, lnZ2 = 0.0;
  gdouble fmax      = G_MAXDOUBLE;
  gdouble fmax_prev = G_MAXDOUBLE;
  guint nphi        = _moments_tilt_nphi_for (lam, sn2);
  guint nphi_floor  = NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NPHI_MIN;
  guint it, pass;

  if (nphi != *nphi_io)
  {
    _moments_tilt_set_phi (cphi, s2phi, nphi);
    *nphi_io = nphi;
  }

  /* Two passes: the solve, then one verification at doubled angular
   * resolution that re-solves if the residual there is above NEWTON_OK. */
  for (pass = 0; pass < 2; pass++)
  {
    for (it = 0; it < maxit; it++)
    {
      gdouble t = 1.0;
      guint h, a;

      if (!_moments_tilt_moments (r, wr, nr, cphi, s2phi, nphi, sn2, lam, E, C, &lnZ))
        return FALSE;

      fmax = 0.0;

      for (a = 0; a < 3; a++)
      {
        F[a] = E[a] - target[a];
        fmax = MAX (fmax, fabs (F[a]));
      }

      if (fmax < NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_TOL)
        break;

      /* Stop once the residual is below NEWTON_OK and has not halved since
       * the previous iteration. TOL is below the floor the moment
       * quadrature reaches at large |lambda|, so without this test the last
       * nodes of a small-noise table run to the iteration cap, each
       * iteration costing a full nr x nphi moment pass plus the line
       * search's own. The test requires the residual to be below NEWTON_OK
       * first, so it never stops a solve short of acceptance. */
      if ((fmax >= 0.5 * fmax_prev) && (fmax < NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_OK))
        break;

      fmax_prev = fmax;

      if (!_moments_tilt_solve3 (C, F, step))
        return FALSE;

      /* Line search: halves the step until the Armijo condition holds on
       * the convex objective phi(lambda) = W(lambda) - lambda.t, whose
       * gradient is the residual F. Both this and a test on ||F|| cost one
       * first-moment pass, but only the Armijo condition on phi has the
       * global convergence of damped Newton on a strictly convex function:
       * a residual-norm test can stall at a point where ||F|| cannot
       * decrease along the Newton ray while phi still can.
       *
       * F.step = F^T C^{-1} F > 0 for a positive-definite covariance, so
       * -step is a descent direction. If F.step <= 0 numerically, the
       * residual-norm test is used instead. */
      {
        gdouble phi_cur = lnZ;
        gdouble gTs     = 0.0;

        for (a = 0; a < 3; a++)
        {
          phi_cur -= lam[a] * target[a];
          gTs     += F[a] * step[a];
        }

        for (h = 0; h < NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_MAXHALF; h++)
        {
          for (a = 0; a < 3; a++)
            trial[a] = lam[a] - t * step[a];

          if (_moments_tilt_moments (r, wr, nr, cphi, s2phi, nphi, sn2, trial, E2, NULL, &lnZ2))
          {
            if (gTs > 0.0)
            {
              gdouble phi_try = lnZ2;

              for (a = 0; a < 3; a++)
                phi_try -= trial[a] * target[a];

              if (phi_try <= phi_cur - NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_ARMIJO_C1 * t * gTs)
                break;
            }
            else
            {
              gdouble f2 = 0.0;

              for (a = 0; a < 3; a++)
                f2 = MAX (f2, fabs (E2[a] - target[a]));

              if (f2 < fmax)
                break;
            }
          }

          t *= 0.5;
        }
      }

      for (a = 0; a < 3; a++)
        lam[a] -= t * step[a];

      {
        /* The node count follows the resolution rule at the new iterate,
         * bounded below by nphi_floor, the count the verification pass
         * set: a smaller count would undo the refinement. */
        const guint n2 = MAX (_moments_tilt_nphi_for (lam, sn2), nphi_floor);

        if (n2 != nphi)
        {
          nphi = n2;
          _moments_tilt_set_phi (cphi, s2phi, nphi);
          *nphi_io = nphi;
        }
      }
    }

    /* At the iteration cap lam took a step after lnZ and fmax were
     * computed, so both are recomputed at the lam that is returned. */
    if (it == maxit)
    {
      guint a;

      if (!_moments_tilt_moments (r, wr, nr, cphi, s2phi, nphi, sn2, lam, E, NULL, &lnZ))
        return FALSE;  /* LCOV_EXCL_LINE */

      fmax = 0.0;

      for (a = 0; a < 3; a++)
        fmax = MAX (fmax, fabs (E[a] - target[a]));
    }

    if (pass == 0)
    {
      const guint n2 = MIN (2 * nphi, NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NPHI_MAX);

      if (n2 == nphi)
        break;

      _moments_tilt_set_phi (cphi, s2phi, n2);
      nphi       = n2;
      nphi_floor = n2;
      *nphi_io   = n2;

      /* The same lambda as the pass that just converged, on a finer
       * angular rule: u, v > 0 still hold, so it is feasible. */
      if (!_moments_tilt_moments (r, wr, nr, cphi, s2phi, nphi, sn2, lam, E2, NULL, &lnZ2))
        return FALSE;  /* LCOV_EXCL_LINE */

      {
        gdouble f2 = 0.0;
        guint a;

        for (a = 0; a < 3; a++)
          f2 = MAX (f2, fabs (E2[a] - target[a]));

        if (f2 <= NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_OK)
        {
          /* The residual holds at the finer rule: keep lam and the refined
           * lnZ. */
          lnZ = lnZ2;
          break;
        }
      }
    }
  }

  *W = lnZ;

  return fmax < NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_OK;
}

/*
 * ---- Stage 2: the dyadic Chebyshev mesh in ghat ----
 *
 * Panel j spans [1 - 2^-j, 1 - 2^-(j+1)], except the last, whose upper end
 * is the table's own @top. Each panel carries its own Chebyshev series of
 * its own degree, and the four series (lambda_1, lambda_2, lambda_3, W)
 * share one interleaved coefficient block so that the evaluation reads each
 * cache line once instead of striding the array four times.
 *
 * off[j] indexes @coef in doubles and includes the factor of four:
 * off[0] = 0 and off[j+1] = off[j] + 4 * (deg[j] + 1). Build and evaluation
 * must agree on this layout: a mismatch reads a valid but wrong part of
 * the block and does not run off its end.
 */

/* Number of panels for the noise level @sn: K = round(log2(1/sn)) clamped
 * to [MIN_PANELS, NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS]; *@clamped is
 * set when the cap binds, and the caller warns. The error does not
 * decrease above K, so K is a stopping rule and not a tuning parameter.
 * The floor of two gives two orders of margin at the large-sn end for one
 * redundant panel. */
static guint
_moments_tilt_panels_for_sn (const gdouble sn, gboolean *clamped)
{
  const gdouble k = (sn > 0.0) ? floor (-log2 (sn) + 0.5) : (gdouble) NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS;

  *clamped = (k > (gdouble) NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS);

  if (k < (gdouble) NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_MIN_PANELS)
    return NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_MIN_PANELS;

  if (*clamped)
    return NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS;

  return (guint) k;
}

/*
 * Size of one coefficient vector @v = (a_k, b_k, c_k, e_k), measured as the
 * largest error it can contribute to ln P: the signed combination the
 * evaluation forms, maximized over four statistic vectors. The four are
 * the target moments @t_mid at the panel midpoint and the corners of the
 * observable disc, |chi_obs| = 1 along and across the gauge axis, where
 * the statistic (x, x^2, y^2) reaches its extremes.
 *
 * The midpoint form alone measures the model-mass error, which is what the
 * accuracy gate is stated in, but it does not determine the resolution: the
 * exact solve satisfies dW/dghat = t . dlambda/dghat, so over a panel where
 * t varies little the coefficients of W cancel t . lambda term by term. The
 * mass error is then small while lambda is unresolved, and the density
 * shape is wrong between nodes. On a Gaussian population of width 0.4658
 * at sigma_nu = 0.015 the midpoint form alone stops a panel at degree 7
 * whose worst |ln Z| is 9.3 nats.
 *
 * The termwise bound |t1 a_k| + |t2 b_k| + |t3 c_k| + |e_k| also avoids
 * the cancellation and is a rigorous upper bound, but it overestimates by
 * a median factor of 45 and selects degrees far above what the function
 * needs.
 */

static gdouble
_moments_tilt_probe_max (const gdouble *v, const gdouble *t_mid)
{
  const gdouble a = v[0];
  const gdouble b = v[1];
  const gdouble c = v[2];
  const gdouble e = v[3];
  gdouble m       = fabs (t_mid[0] * a + t_mid[1] * b + t_mid[2] * c - e);

  m = MAX (m, fabs (a + b - e));
  m = MAX (m, fabs (-a + b - e));
  m = MAX (m, fabs (c - e));

  return m;
}

/*
 * Writes to @lam the extrapolation of lambda to @ghat from the @n_hist
 * previous converged nodes: constant from one, linear from two, quadratic
 * (Lagrange) from three. The quadratic start reduces the Newton iteration
 * count of a table by about a third against starting from the previous
 * node.
 *
 * The history is a shift register whose most recent entry is always slot 2,
 * so a partially filled one is read from the end: slots 3 - @n_hist to 2.
 */
static void
_moments_tilt_warm_start (const gdouble *g_hist, const gdouble *lam_hist, const guint n_hist, const gdouble ghat, gdouble *lam)
{
  guint a;

  g_assert_cmpuint (n_hist, >=, 1);

  if (n_hist == 1)
  {
    for (a = 0; a < 3; a++)
      lam[a] = lam_hist[6 + a];

    return;
  }

  if (n_hist == 2)
  {
    const gdouble d12 = g_hist[2] - g_hist[1];
    const gdouble w   = (d12 != 0.0) ? (ghat - g_hist[1]) / d12 : 0.0;

    for (a = 0; a < 3; a++)
      lam[a] = lam_hist[3 + a] + w * (lam_hist[6 + a] - lam_hist[3 + a]);

    return;
  }

  {
    const gdouble x0 = g_hist[0], x1 = g_hist[1], x2 = g_hist[2];
    const gdouble d01 = x0 - x1, d02 = x0 - x2, d12 = x1 - x2;

    /* The abscissae are distinct memo entries, never closer than the memo's
     * lookup tolerance, so no divisor vanishes. */
    {
      const gdouble l0 = ((ghat - x1) * (ghat - x2)) / (d01 * d02);
      const gdouble l1 = ((ghat - x0) * (ghat - x2)) / (-d01 * d12);
      const gdouble l2 = ((ghat - x0) * (ghat - x1)) / (d02 * d12);

      for (a = 0; a < 3; a++)
        lam[a] = l0 * lam_hist[a] + l1 * lam_hist[3 + a] + l2 * lam_hist[6 + a];
    }
  }
}

/*
 * Clamps lambda_2 and lambda_3 of a starting point to 0.45 / sigma_nu^2,
 * inside the natural domain of W, which needs 1 - 2 lambda_2 sigma_nu^2 > 0
 * and the same for lambda_3. An extrapolated start can overshoot that
 * bound, and the first moment evaluation at an infeasible point returns
 * FALSE, which fails the solve.
 */
static void
_moments_tilt_clamp_start (gdouble *lam, const gdouble sn2)
{
  const gdouble cap = 0.45 / sn2;

  lam[1] = MIN (lam[1], cap);
  lam[2] = MIN (lam[2], cap);
}

/*
 * ---- The table build, on the shared NcmSpectral-driven builder ----
 *
 * This class supplies the node callback: the exact target moments at the
 * folded shear, then the damped Newton solve for (lambda, W). The builder
 * places the nodes, refines, converges and transforms; see
 * _nc_galaxy_shape_factor_moments_build().
 */
typedef struct _NcGalaxyShapeFactorMomentsTiltBuild
{
  NcGalaxyWLObsEllipConv ellip_conv;
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  gdouble sn2;
  GArray *r_arr;
  GArray *p_arr;
  gdouble *cphi;
  gdouble *s2phi;
  gdouble *src_r;
  gdouble *src_w;
  guint nr;
  guint nphi_cur;
  guint newton_maxit;
} NcGalaxyShapeFactorMomentsTiltBuild;

/*
 * Writes to @lam the starting point for the node at @ghat, from the solved
 * nodes around it. Between two solved neighbors (every refinement node) it
 * interpolates linearly. Above all of them (the first level of a panel) it
 * extrapolates from the three below with _moments_tilt_warm_start().
 */
static void
_moments_tilt_start_from_memo (const NcGalaxyShapeFactorMomentsMemo *memo, const guint below, const gdouble ghat, gdouble *lam)
{
  const guint n = _nc_galaxy_shape_factor_moments_memo_len (memo);
  guint a;

  /* ghat = 0 is the first node of every build and is never solved (see
   * _moments_tilt_node_cb()), so a node being solved always has one below. */
  g_assert (below != NC_GALAXY_SHAPE_FACTOR_MOMENTS_MEMO_NONE);

  if (below + 1 < n)
  {
    const gdouble x0  = _nc_galaxy_shape_factor_moments_memo_x (memo, below);
    const gdouble x1  = _nc_galaxy_shape_factor_moments_memo_x (memo, below + 1);
    const gdouble *v0 = _nc_galaxy_shape_factor_moments_memo_vals (memo, below);
    const gdouble *v1 = _nc_galaxy_shape_factor_moments_memo_vals (memo, below + 1);
    const gdouble w   = (x1 > x0) ? (ghat - x0) / (x1 - x0) : 0.0;

    for (a = 0; a < 3; a++)
      lam[a] = v0[a] + w * (v1[a] - v0[a]);

    return;
  }

  {
    gdouble g_hist[3]   = { 0.0, 0.0, 0.0 };
    gdouble lam_hist[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    const guint n_hist  = MIN (below + 1, 3);
    guint i;

    /* Shift-register layout: the most recent entry is always slot 2. */
    for (i = 0; i < n_hist; i++)
    {
      const guint src  = below - (n_hist - 1 - i);
      const guint slot = 3 - n_hist + i;
      const gdouble *v = _nc_galaxy_shape_factor_moments_memo_vals (memo, src);

      g_hist[slot]           = _nc_galaxy_shape_factor_moments_memo_x (memo, src);
      lam_hist[3 * slot + 0] = v[0];
      lam_hist[3 * slot + 1] = v[1];
      lam_hist[3 * slot + 2] = v[2];
    }

    _moments_tilt_warm_start (g_hist, lam_hist, n_hist, ghat, lam);
  }
}

static gboolean
_moments_tilt_node_cb (gpointer user_data, const gdouble ghat, const NcGalaxyShapeFactorMomentsMemo *memo, gdouble *vals)
{
  NcGalaxyShapeFactorMomentsTiltBuild *b = (NcGalaxyShapeFactorMomentsTiltBuild *) user_data;
  const guint below                      = _nc_galaxy_shape_factor_moments_memo_below (memo, ghat);
  gdouble target[3], lam[3];
  gdouble W = 0.0;
  gboolean ok;

  /* lambda(0) = 0 and W(0) = 0 exactly: never solved. */
  if (ghat <= 0.0)
  {
    vals[0] = vals[1] = vals[2] = vals[3] = 0.0;

    return TRUE;
  }

  _nc_galaxy_shape_factor_moments_exact_t (b->ellip_conv, b->pop, b->pop_data, ghat, b->sn2, b->r_arr, &b->p_arr, target);

  _moments_tilt_start_from_memo (memo, below, ghat, lam);
  _moments_tilt_clamp_start (lam, b->sn2);

  ok = _moments_tilt_newton (b->src_r, b->src_w, b->nr, b->cphi, b->s2phi, &b->nphi_cur, b->sn2, target, b->newton_maxit, lam, &W);

  if (!ok)
  {
    /* Retry from the nearest solved node below, the start a plain
     * continuation would use. */
    const gdouble *v = _nc_galaxy_shape_factor_moments_memo_vals (memo, below);

    lam[0] = v[0];
    lam[1] = v[1];
    lam[2] = v[2];
    ok     = _moments_tilt_newton (b->src_r, b->src_w, b->nr, b->cphi, b->s2phi, &b->nphi_cur, b->sn2, target, b->newton_maxit, lam, &W);

    if (!ok)
    {
      /* Keep the neighbor's value rather than a diverged iterate, which
       * would spread through the whole panel. */
      vals[0] = v[0];
      vals[1] = v[1];
      vals[2] = v[2];
      vals[3] = v[3];

      return FALSE;
    }
  }

  vals[0] = lam[0];
  vals[1] = lam[1];
  vals[2] = lam[2];
  vals[3] = W;

  return ok;
}

static void
_moments_tilt_panel_ctx_cb (gpointer user_data, const gdouble lo, const gdouble hi, gdouble *ctx)
{
  NcGalaxyShapeFactorMomentsTiltBuild *b = (NcGalaxyShapeFactorMomentsTiltBuild *) user_data;

  _nc_galaxy_shape_factor_moments_exact_t (b->ellip_conv, b->pop, b->pop_data, 0.5 * (lo + hi), b->sn2, b->r_arr, &b->p_arr, ctx);
}

static gdouble
_moments_tilt_probe_cb (gpointer user_data, const gdouble *v, const gdouble *ctx)
{
  return _moments_tilt_probe_max (v, ctx);
}

static NcGalaxyShapeFactorMomentsTable *
_moments_tilt_build_table (NcmSpectral *spectral, const NcGalaxyWLObsEllipConv ellip_conv, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                           const gdouble sn, const guint n_panels, const gdouble top,
                           const gdouble gate, const guint max_degree, const guint newton_maxit)
{
  const guint nq    = NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES;
  const guint nr    = NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_SRC_NR;
  const gdouble tol = gate / NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_GATE_SAFETY;
  NcGalaxyShapeFactorMomentsTiltBuild b;
  NcGalaxyShapeFactorMomentsTable *table;
  guint i;

  _moments_tilt_gl_init ();

  b.ellip_conv = ellip_conv;

  b.pop          = pop;
  b.pop_data     = pop_data;
  b.sn2          = sn * sn;
  b.r_arr        = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), MAX (nq, nr));
  b.p_arr        = NULL;
  b.cphi         = g_new0 (gdouble, NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NPHI_MAX);
  b.s2phi        = g_new0 (gdouble, NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NPHI_MAX);
  b.src_r        = g_new0 (gdouble, nr);
  b.src_w        = g_new0 (gdouble, nr);
  b.nr           = nr;
  b.nphi_cur     = 0;
  b.newton_maxit = newton_maxit;

  /* Source-disc nodes and weights, independent of the shear, so built once
   * for the whole table. sum src_w = int P_pop dr = 1. */
  g_array_set_size (b.r_arr, nr);
  {
    gdouble *r_data = (gdouble *) b.r_arr->data;

    for (i = 0; i < nr; i++)
      r_data[i] = 0.5 * (_moments_tilt_gl32_x[i] + 1.0);

    nc_galaxy_shape_pop_eval_p_array (pop, pop_data, b.r_arr, &b.p_arr);

    for (i = 0; i < nr; i++)
    {
      b.src_r[i] = r_data[i];
      b.src_w[i] = 0.5 * _moments_tilt_gl32_w[i] * g_array_index (b.p_arr, gdouble, i);
    }
  }

  g_array_set_size (b.r_arr, nq);

  table = _nc_galaxy_shape_factor_moments_build (spectral, 4, n_panels, top, max_degree, tol, tol,
                                                 &_moments_tilt_node_cb, &_moments_tilt_panel_ctx_cb,
                                                 &_moments_tilt_probe_cb, &b);

  g_array_unref (b.r_arr);

  if (b.p_arr != NULL)
    g_array_unref (b.p_arr);

  g_free (b.cphi);
  g_free (b.s2phi);
  g_free (b.src_r);
  g_free (b.src_w);

  return table;
}

/*
 * Evaluates lambda_1 x + lambda_2 x^2 + lambda_3 y^2 - W on panel @j at the
 * Chebyshev argument @t by one fused Clenshaw recurrence: the four series
 * share the argument and the weights (x, x^2, y^2, -1) do not depend on the
 * summation index, so the combination is formed inside the recurrence.
 * That is 9d flops instead of 12d, over one pass of the interleaved block
 * instead of four strided ones.
 */
static gdouble
_moments_tilt_clenshaw_fused (const NcGalaxyShapeFactorMomentsTable *table, const guint j,
                              const gdouble t, const gdouble x, const gdouble x2, const gdouble y2)
{
  const gdouble *coef = &table->coef[table->off[j]];
  const guint d       = table->deg[j];
  const gdouble t2    = 2.0 * t;
  gdouble b1 = 0.0, b2 = 0.0;
  guint k;

  /* b0 = t2 b1 + (ck - b2): ck and b2 are known before b1 is, so only
   * the multiply-add is on the loop-carried dependency chain. */
  for (k = d; k >= 1; k--)
  {
    const gdouble ck = coef[4 * k + 0] * x + coef[4 * k + 1] * x2 + coef[4 * k + 2] * y2 - coef[4 * k + 3];
    const gdouble b0 = t2 * b1 + (ck - b2);

    b2 = b1;
    b1 = b0;
  }

  return t * b1 - b2 + (coef[0] * x + coef[1] * x2 + coef[2] * y2 - coef[3]);
}

/*
 * ---- ln P_0 ----
 *
 * P_0(R) = (1/2 pi sn^2) int_0^1 P_pop^NC(r) exp(-(R-r)^2/2sn^2)
 *          I0_scaled(R r/sn^2) dr, in NumCosmo's own radial-marginal
 * convention (nc_galaxy_shape_pop_eval_p(): the disc measure 2 pi r
 * already folded in).
 *
 * The window is localized around the mass, and the exponents are combined
 * in log space. Any error here is a per-galaxy additive constant in ln P
 * and cancels from every derivative in g, so a generic quadrature is used
 * rather than a closed form.
 */
static gdouble _moments_tilt_lnP0_x[NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_LNP0_NNODES];
static gdouble _moments_tilt_lnP0_w[NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_LNP0_NNODES];

static void
_moments_tilt_lnP0_nodes_init (void)
{
  static gsize init = 0;

  if (g_once_init_enter (&init))
  {
    _nc_galaxy_shape_factor_moments_gl_fill (_moments_tilt_lnP0_x, _moments_tilt_lnP0_w,
                                             NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_LNP0_NNODES);
    g_once_init_leave (&init, 1);
  }
}

static gdouble
_moments_tilt_ln_P0 (NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, gdouble R, gdouble sn)
{
  const guint n_nodes = NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_LNP0_NNODES;
  const gdouble sn2   = sn * sn;
  gdouble r_lo, r_hi, half, mid;
  GArray *r_arr;
  GArray *p_arr = NULL;
  gdouble *r_data;
  gdouble expo[NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_LNP0_NNODES];
  gdouble max_expo = -G_MAXDOUBLE;
  gdouble sum      = 0.0;
  guint i;

  _moments_tilt_lnP0_nodes_init ();

  if (R <= 1.0)
  {
    const gdouble half_width = sn * sqrt (2.0 * NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_LNP0_DELTA);

    r_lo = MAX (0.0, R - half_width);
    r_hi = MIN (1.0, R + half_width);
  }
  else
  {
    const gdouble d          = R - 1.0;
    const gdouble half_width = sqrt (d * d + 2.0 * sn2 * NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_LNP0_DELTA);

    r_lo = MAX (0.0, R - half_width);
    r_hi = 1.0;
  }

  half  = 0.5 * (r_hi - r_lo);
  mid   = 0.5 * (r_hi + r_lo);
  r_arr = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_nodes);
  g_array_set_size (r_arr, n_nodes);
  r_data = (gdouble *) r_arr->data;

  for (i = 0; i < n_nodes; i++)
    r_data[i] = mid + half * _moments_tilt_lnP0_x[i];

  nc_galaxy_shape_pop_eval_p_array (pop, pop_data, r_arr, &p_arr);

  for (i = 0; i < n_nodes; i++)
  {
    const gdouble r  = r_data[i];
    const gdouble dr = R - r;
    const gdouble z  = R * r / sn2;

    expo[i]  = -0.5 * dr * dr / sn2 + log (gsl_sf_bessel_I0_scaled (z));
    max_expo = MAX (max_expo, expo[i]);
  }

  {
    const gdouble *p_data = (const gdouble *) p_arr->data;

    for (i = 0; i < n_nodes; i++)
      sum += _moments_tilt_lnP0_w[i] * p_data[i] * exp (expo[i] - max_expo);
  }

  g_array_unref (r_arr);
  g_array_unref (p_arr);

  return max_expo + log (half * sum) - log (2.0 * M_PI * sn2);
}

struct _NcGalaxyShapeFactorMomentsTilt
{
  NcGalaxyShapeFactor parent_instance;
};

/*
 * ---- Instance-level table cache ----
 *
 * A table is a function of (population parameters, sigma_nu, e_rms, panel
 * count) only: it does not depend on the shear or on the observed
 * ellipticity, and NcGalaxyShapePopData's only per-galaxy input is e_rms.
 * The population enters as a generation counter rather than as part of
 * every key, so a population change drops the whole cache at once instead
 * of leaving one stale generation per step.
 *
 * The cache is held by the instance, not by the per-galaxy ldata: every
 * galaxy sharing a key reads the same table, and NcDataClusterWLFactor
 * rebuilds every NcGalaxyShapeFactorData on resample, which would discard
 * per-galaxy copies. Each galaxy holds a reference, not a copy, so the
 * whole-cache flush above is safe.
 *
 * Concurrency: tab_cache is protected by cache_lock and dp_box by
 * range_lock (see the private struct). NcmFitESMCMC's threaded modes
 * deep-copy the dataset per thread, so each thread owns a distinct factor
 * instance; the locks cover nc_data_cluster_wl_factor_data_prepare()
 * building tables for several galaxies of one instance at once.
 */

typedef struct _NcGalaxyShapeFactorMomentsTiltPrivate
{
  NcGalaxyWLObsEllipConv ellip_conv;
  gdouble accuracy_gate;
  guint max_degree;
  guint newton_maxit;
  guint64 pop_hash;

  /* (sn, e_rms, n_panels) -> table, valid for pop_hash == tab_cache_hash. */
  GHashTable *tab_cache;
  guint64 tab_cache_hash;

  gboolean strict_solve;
  gboolean restrict_range;

  /* Cached by prepare() for the reachable-range bound: the models the bound
   * needs, plus a private copy of the density profile held at the upper
   * bound of the mass. The bound is evaluated on the copy so that the
   * profile in the mset is not modified. */
  NcHICosmo *cosmo;
  NcWLSurfaceMassDensity *smd;
  NcHaloDensityProfile *dp_box;
  gdouble z_cl;
  gdouble delta_R;
  gboolean range_ready;

  /* Whether the per-galaxy bound holds for the current free parameters:
   * only the mass, the halo center and the shape population may be free
   * (see _moments_tilt_range_holds). @range_gen is incremented when it
   * changes, so tables cut to a bound are refreshed. */
  gboolean range_holds;
  guint range_gen;

  /* One NcmSpectral per concurrent build: it holds per-instance plans and
   * buffers. */
  NcmMemoryPool *spectral_pool;

  /* tab_cache is shared by every galaxy of the instance, and
   * nc_data_cluster_wl_factor_data_prepare() may build tables for several
   * at once. dp_box is shared too, and a density profile updates state
   * lazily on use. Each has its own lock; neither lock is held while a
   * table is built. */
  GMutex cache_lock;
  GMutex range_lock;

  /* Serialized table set. @pending_tables holds what a load brought in,
   * not yet validated: the stamp can only be checked once the population is
   * known, at prepare(), and properties are restored in file order, so
   * neither property may depend on the other having been set. */
  NcmObjArray *pending_tables;
  gchar *tables_stamp;
  NcGalaxyShapePop *pop_ref;
  gint stamp_warned;

  /* Memoized PROP_TABLES payload: the same NcmVector instances are returned
   * across gets as long as @tab_cache is unchanged (@tables_generation
   * tracks that). Registering them with an NcmSerialize (see
   * nc_galaxy_shape_factor_moments_tilt_register_shared()) lets every
   * ESMCMC worker thread's clone reuse the identical, read-only vectors
   * instead of re-serializing and rebuilding every table, which is
   * otherwise the dominant per-thread setup cost. */
  NcmObjArray *shared_tables_cache;
  guint64 shared_tables_cache_gen;
  guint64 tables_generation;

  /* Counters and one-shot warning flags, accessed only through
   * g_atomic_int_*; reset only by the public reset functions. */
  gint solve_error_count;
  gint solve_warned;
  gint table_build_count;
  gint range_error_count;
  gint range_warned;
  gint panels_warned;
} NcGalaxyShapeFactorMomentsTiltPrivate;

/*
 * Per-galaxy scratch, split in two: the table depends only on
 * (pop_hash, sn, panel count) and ln_P0 additionally on R, so an R-only
 * change never re-runs the Newton solves.
 */
typedef struct _NcGalaxyShapeFactorMomentsTiltLData
{
  NcGalaxyShapeFactorMomentsTable *table;
  gdouble sn_seen;
  guint64 pop_hash_seen_tab;
  gboolean tab_valid;

  /* Reachable range, filled by data_prepare(). With range_valid FALSE the
   * galaxy never saw data_prepare() and uses the full mesh, which is
   * correct and slower. */
  gdouble ghat_max;
  gboolean range_valid;
  guint range_gen_seen;
  guint n_panels_seen;

  gdouble ln_P0;
  gdouble R2_seen;
  guint64 pop_hash_seen_p0;
  gboolean p0_valid;

  /* Bytes of @table the evaluation reads, set with @table: what
   * data_prefetch() fetches without reading the table header. */
  gsize table_span;
} NcGalaxyShapeFactorMomentsTiltLData;

enum
{
  PROP_0,
  PROP_ACCURACY_GATE,
  PROP_MAX_DEGREE,
  PROP_NEWTON_MAXIT,
  PROP_RESTRICT_RANGE,
  PROP_STRICT_SOLVE,
  PROP_TABLES,
  PROP_TABLES_STAMP,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorMomentsTilt, nc_galaxy_shape_factor_moments_tilt, NC_TYPE_GALAXY_SHAPE_FACTOR)

/*
 * ---- Serializing the table set ----
 *
 * A table is a function of its key only, so a run that has built the tables
 * can pass them to the next run. They are stored as an NcmObjArray of
 * NcmVector, one per table, each carrying its own key and layout:
 *
 *   [ sigma_nu, e_rms, n_panels, top, n_failed, deg_0 .. deg_{n-1}, coef .. ]
 *
 * Only the degrees are stored, never the offsets: off[] is a prefix sum of
 * the degrees, and storing it too would put the "offsets count doubles"
 * convention on both sides of the format boundary, where the two copies
 * could disagree.
 *
 * Validity is checked by a stamp, not a hash. The pkey counters the
 * in-process caches use are per-process and carry no meaning after a load,
 * so the stamp is built from values: the format tag, the build settings,
 * and the serialized population the tables were solved against. A mismatch
 * drops the whole set and rebuilds; it is not an error.
 */
#define NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_TABLES_FORMAT "moments-tilt-v1"

static gchar *
_moments_tilt_stamp (NcGalaxyShapeFactorMomentsTiltPrivate * const self)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *pop_desc   = ncm_serialize_to_string (ser, G_OBJECT (self->pop_ref), TRUE);
  gchar *stamp      = g_strdup_printf ("%s|ellip-conv=%d|gate=%.17g|max-degree=%u|newton-maxit=%u|pop=%s",
                                       NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_TABLES_FORMAT,
                                       (gint) self->ellip_conv, self->accuracy_gate, self->max_degree, self->newton_maxit, pop_desc);

  g_free (pop_desc);
  ncm_serialize_free (ser);

  return stamp;
}

/* Allocator of the NcmSpectral pool. The pool keeps the returned pointer in
 * a slot of its own and hands out the slot: ncm_memory_pool_get() returns
 * an NcmSpectral **. */
static gpointer
_moments_tilt_spectral_alloc (gpointer userdata)
{
  return ncm_spectral_new_with_max_level (6);
}

static void
_moments_tilt_spectral_free (gpointer p)
{
  ncm_spectral_free (NCM_SPECTRAL (p));
}

static void
nc_galaxy_shape_factor_moments_tilt_init (NcGalaxyShapeFactorMomentsTilt *gsfmt)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);

  self->ellip_conv = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE;

  self->accuracy_gate = 1.0e-3;
  self->max_degree    = 16;
  self->newton_maxit  = 0;
  self->pop_hash      = 0;
  self->tab_cache     = g_hash_table_new_full (&_nc_galaxy_shape_factor_moments_key_hash, &_nc_galaxy_shape_factor_moments_key_equal,
                                               &g_free, (GDestroyNotify) & _nc_galaxy_shape_factor_moments_table_unref);
  self->tab_cache_hash = 0;
  self->strict_solve   = FALSE;
  self->restrict_range = TRUE;
  self->cosmo          = NULL;
  self->smd            = NULL;
  self->dp_box         = NULL;
  self->z_cl           = 0.0;
  self->delta_R        = 0.0;
  self->range_ready    = FALSE;
  self->range_holds    = FALSE;
  self->range_gen      = 0;
  self->spectral_pool  = ncm_memory_pool_new (&_moments_tilt_spectral_alloc, NULL, &_moments_tilt_spectral_free);
  g_mutex_init (&self->cache_lock);
  g_mutex_init (&self->range_lock);
  self->pending_tables          = NULL;
  self->tables_stamp            = NULL;
  self->pop_ref                 = NULL;
  self->stamp_warned            = 0;
  self->shared_tables_cache     = NULL;
  self->shared_tables_cache_gen = 0;
  self->tables_generation       = 0;
  self->solve_error_count       = 0;
  self->solve_warned            = 0;
  self->table_build_count       = 0;
  self->range_error_count       = 0;
  self->range_warned            = 0;
  self->panels_warned           = 0;
}

/* Reads the ellipticity convention, a construct-only property of the
 * parent, once at construction. The target moments branch on it, and a
 * table built under one convention is wrong under the other. */
static void
_nc_galaxy_shape_factor_moments_tilt_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_moments_tilt_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorMomentsTiltPrivate * const self =
      nc_galaxy_shape_factor_moments_tilt_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (object));

    self->ellip_conv = nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (object));
  }
}

static void
_nc_galaxy_shape_factor_moments_tilt_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorMomentsTilt *gsfmt              = NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (object);
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);

  switch (prop_id)
  {
    case PROP_ACCURACY_GATE:
      self->accuracy_gate = g_value_get_double (value);
      break;
    case PROP_MAX_DEGREE:
      self->max_degree = g_value_get_uint (value);
      break;
    case PROP_NEWTON_MAXIT:
      self->newton_maxit = g_value_get_uint (value);
      break;
    case PROP_RESTRICT_RANGE:
      self->restrict_range = g_value_get_boolean (value);
      break;
    case PROP_STRICT_SOLVE:
      self->strict_solve = g_value_get_boolean (value);
      break;
    case PROP_TABLES:
    {
      NcmObjArray *oa = g_value_get_boxed (value);

      g_clear_pointer (&self->pending_tables, ncm_obj_array_unref);

      if (oa != NULL)
        self->pending_tables = ncm_obj_array_ref (oa);

      break;
    }
    case PROP_TABLES_STAMP:
      g_clear_pointer (&self->tables_stamp, g_free);
      self->tables_stamp = g_value_dup_string (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/*
 * Returns a reference to the memoized "tables" payload, rebuilding it only
 * when @tab_cache has changed (@tables_generation) since the last build.
 * The NcmVector instances inside are the same across calls as long as the
 * generation does not change, so
 * nc_galaxy_shape_factor_moments_tilt_register_shared() can register them
 * in an NcmSerialize: a later ncm_serialize_to_variant() of the same
 * instance finds each vector already named and emits a reference instead of
 * re-serializing it.
 */
static NcmObjArray *
_nc_galaxy_shape_factor_moments_tilt_get_shared_tables (NcGalaxyShapeFactorMomentsTiltPrivate * const self)
{
  NcmObjArray *oa;

  g_mutex_lock (&self->cache_lock);

  if ((self->shared_tables_cache == NULL) || (self->shared_tables_cache_gen != self->tables_generation))
  {
    NcmObjArray *built = ncm_obj_array_new ();
    GHashTableIter iter;
    gpointer k, v;

    g_hash_table_iter_init (&iter, self->tab_cache);

    while (g_hash_table_iter_next (&iter, &k, &v))
    {
      NcmVector *vec = _nc_galaxy_shape_factor_moments_table_to_vector ((const NcGalaxyShapeFactorMomentsKey *) k,
                                                                        (const NcGalaxyShapeFactorMomentsTable *) v);

      ncm_obj_array_add (built, G_OBJECT (vec));
      ncm_vector_free (vec);
    }

    g_clear_pointer (&self->shared_tables_cache, ncm_obj_array_unref);
    self->shared_tables_cache     = built;
    self->shared_tables_cache_gen = self->tables_generation;
  }

  oa = ncm_obj_array_ref (self->shared_tables_cache);

  g_mutex_unlock (&self->cache_lock);

  return oa;
}

static void
_nc_galaxy_shape_factor_moments_tilt_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorMomentsTilt *gsfmt              = NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (object);
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);

  switch (prop_id)
  {
    case PROP_ACCURACY_GATE:
      g_value_set_double (value, self->accuracy_gate);
      break;
    case PROP_MAX_DEGREE:
      g_value_set_uint (value, self->max_degree);
      break;
    case PROP_NEWTON_MAXIT:
      g_value_set_uint (value, self->newton_maxit);
      break;
    case PROP_RESTRICT_RANGE:
      g_value_set_boolean (value, self->restrict_range);
      break;
    case PROP_STRICT_SOLVE:
      g_value_set_boolean (value, self->strict_solve);
      break;
    case PROP_TABLES:
      g_value_take_boxed (value, _nc_galaxy_shape_factor_moments_tilt_get_shared_tables (self));
      break;
    case PROP_TABLES_STAMP:
      g_value_take_string (value, (self->pop_ref != NULL) ? _moments_tilt_stamp (self) : NULL);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_moments_tilt_ldata_destroy (gpointer p)
{
  NcGalaxyShapeFactorMomentsTiltLData *ldata = (NcGalaxyShapeFactorMomentsTiltLData *) p;

  g_clear_pointer (&ldata->table, _nc_galaxy_shape_factor_moments_table_unref);
  g_free (ldata);
}

static void
_nc_galaxy_shape_factor_moments_tilt_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

/* Invalidates the per-galaxy table and ln_P0 when a catalog row is read: a
 * per-galaxy population reads its moments from the row, so a new row
 * changes them without any model pkey moving. This covers only row reads;
 * the sn and R2 value checks in peek_table() cover
 * nc_galaxy_shape_factor_data_set() and gen(). */
static void
_nc_galaxy_shape_factor_moments_tilt_ldata_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorMomentsTiltLData *ldata = (NcGalaxyShapeFactorMomentsTiltLData *) data->ldata;

  ldata->tab_valid = FALSE;
  ldata->p0_valid  = FALSE;
}

static void
_nc_galaxy_shape_factor_moments_tilt_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_moments_tilt_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorMomentsTiltLData *ldata = g_new0 (NcGalaxyShapeFactorMomentsTiltLData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_moments_tilt_ldata_destroy;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_moments_tilt_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_moments_tilt_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_moments_tilt_ldata_required_columns;
}

/*
 * Whether the reachable-range bound of data_prepare() holds for the free
 * parameters of @mset: TRUE when only log10MDelta, the halo center (ra,
 * dec) and shape-population parameters are free. The bound covers
 * log10MDelta up to its upper bound and the halo center over its box; the
 * shape population changes only the tables, which are keyed on it. Every
 * other free parameter (cosmology, z_cl, cDelta and the other profile and
 * mass-summary parameters, the surface mass density, the redshift models)
 * moves the shear or the redshift support outside what the bound covers.
 */
static gboolean
_moments_tilt_range_holds (NcmMSet *mset)
{
  const guint n_models = ncm_mset_nmodels (mset);
  guint i, pid;

  /* Read from the models, not from the free-parameter map, which the
   * caller may not have prepared. */
  for (i = 0; i < n_models; i++)
  {
    NcmModel *model = ncm_mset_peek (mset, ncm_mset_get_mid_array_pos (mset, i));
    const guint len = ncm_model_len (model);

    if (NC_IS_GALAXY_SHAPE_POP (model))
      continue;

    for (pid = 0; pid < len; pid++)
    {
      if (ncm_model_param_get_ftype (model, pid) != NCM_PARAM_TYPE_FREE)
        continue;

      if (NC_IS_HALO_POSITION (model) && ((pid == NC_HALO_POSITION_RA) || (pid == NC_HALO_POSITION_DEC)))
        continue;

      if (NC_IS_HALO_MASS_SUMMARY (model))
      {
        const gchar *name = ncm_model_param_name (model, pid);

        if (g_strcmp0 (name, "log10MDelta") == 0)
          continue;
      }

      return FALSE;
    }
  }

  return TRUE;
}

static void
_nc_galaxy_shape_factor_moments_tilt_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self =
    nc_galaxy_shape_factor_moments_tilt_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (gsf));

  /* No capability gate: only nc_galaxy_shape_pop_eval_p_array() and
   * nc_galaxy_shape_pop_moment_2k() are used, which every
   * NcGalaxyShapePop provides. */
  self->pop_hash = nc_galaxy_shape_factor_get_pop_hash (gsf);

  {
    NcGalaxyShapePop *pop = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

    if (pop != self->pop_ref)
    {
      nc_galaxy_shape_pop_clear (&self->pop_ref);
      self->pop_ref = (pop != NULL) ? nc_galaxy_shape_pop_ref (pop) : NULL;
    }
  }

  if (self->tab_cache_hash != self->pop_hash)
  {
    g_hash_table_remove_all (self->tab_cache);
    self->tab_cache_hash = self->pop_hash;
    self->tables_generation++;
  }

  /* Adopt the tables a load brought in, once, and only if their stamp
   * matches this population and these build settings. This runs after the
   * generation flush above, so an adopted set is not dropped at once. */
  if (G_UNLIKELY (self->pending_tables != NULL))
  {
    NcmObjArray *pending = self->pending_tables;
    gchar *stamp_now     = (self->pop_ref != NULL) ? _moments_tilt_stamp (self) : NULL;
    const gboolean match = (stamp_now != NULL) && (self->tables_stamp != NULL) &&
                           (g_strcmp0 (stamp_now, self->tables_stamp) == 0);
    guint adopted = 0;
    guint i;

    self->pending_tables = NULL;

    for (i = 0; match && (i < ncm_obj_array_len (pending)); i++)
    {
      GObject *obj = ncm_obj_array_peek (pending, i);
      NcGalaxyShapeFactorMomentsKey key;
      NcGalaxyShapeFactorMomentsTable *table;

      if (!NCM_IS_VECTOR (obj))
        continue;

      if (_nc_galaxy_shape_factor_moments_table_from_vector (NCM_VECTOR (obj), &key, &table))
      {
        NcGalaxyShapeFactorMomentsKey *key_copy = g_new (NcGalaxyShapeFactorMomentsKey, 1);

        *key_copy = key;
        g_hash_table_insert (self->tab_cache, key_copy, table);
        adopted++;
      }
    }

    if (adopted > 0)
      self->tables_generation++;

    if (!match && (ncm_obj_array_len (pending) > 0) &&
        g_atomic_int_compare_and_exchange (&self->stamp_warned, 0, 1))
      g_warning ("NcGalaxyShapeFactorMomentsTilt: the %u stored tables were "
                 "built for a different population or different build settings, "
                 "so they were discarded and will be rebuilt. Results are "
                 "unaffected; only the time the build costs is.",
                 ncm_obj_array_len (pending));

    self->tab_cache_hash = self->pop_hash;

    ncm_obj_array_unref (pending);
    g_free (stamp_now);
  }

  self->range_ready = FALSE;

  if (self->restrict_range)
  {
    NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
    NcHaloPosition *hp          = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
    NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
    NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));

    if ((cosmo != NULL) && (hp != NULL) && (smd != NULL) && (dp != NULL))
    {
      if (self->dp_box == NULL)
      {
        NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

        self->dp_box = NC_HALO_DENSITY_PROFILE (ncm_model_dup (NCM_MODEL (dp), ser));

        ncm_serialize_free (ser);
      }

      nc_wl_surface_mass_density_prepare_if_needed (smd, cosmo);
      nc_halo_position_prepare_if_needed (hp, cosmo);

      nc_hicosmo_clear (&self->cosmo);
      nc_wl_surface_mass_density_clear (&self->smd);
      self->cosmo = nc_hicosmo_ref (cosmo);
      self->smd   = nc_wl_surface_mass_density_ref (smd);

      self->z_cl    = nc_halo_position_get_redshift (hp);
      self->delta_R = _nc_galaxy_shape_factor_moments_centre_delta (hp, cosmo);

      self->range_ready = _nc_galaxy_shape_factor_moments_set_box_corner (self->dp_box, dp);
    }
  }

  {
    const gboolean holds = self->restrict_range && _moments_tilt_range_holds (mset);

    if (holds != self->range_holds)
    {
      self->range_holds = holds;
      self->range_gen++;
    }
  }
}

/*
 * Per-galaxy reachable range: writes ghat_max to the galaxy's ldata and
 * builds its table. The dyadic mesh exists to resolve the layer at
 * ghat = 1, and a galaxy reaches that layer only if its line of sight can
 * cross the tangential critical curve somewhere in the parameter box the
 * sampler explores. Most do not, and their tables stop below the layer.
 *
 * The bound is taken at a fixed concentration and needs one convergence
 * and one shear evaluation at the corner (R_min, z_sup, M_max). At fixed
 * (R, z_s, c) both kappa and gamma grow with M and with z_s. If kappa < 1
 * at the corner, |g| = gamma / (1 - kappa) is largest there and ghat =
 * min(|g|, 1/|g|) is bounded by min(|g|, 1). If kappa >= 1 at the corner,
 * |g| passed through 1 at a smaller mass or redshift, so ghat reaches 1 and
 * the full mesh is kept. The shear is not monotonic in the concentration
 * (outside the scale radius a higher concentration lowers it), so the bound
 * is not taken over a concentration range: a free cDelta disables the
 * restriction (see _moments_tilt_range_holds).
 *
 * The inner radial cut is not used as a floor on R: no r_min weighting is
 * applied at fit time, so a galaxy selected at the fiducial center is
 * still evaluated at whatever radius the current center gives.
 *
 * The shear the table is read at is the calibrated one,
 * g = (1 + m) g_t e^{2i phi} + c (see nc_galaxy_shape_factor.c), not the
 * bare reduced shear the two calls return, so the bound includes the
 * galaxy's own m and c: |g| <= |1 + m| |g_t| + |c| by the triangle
 * inequality. Without them the bound is too small by up to m |g_t| for
 * every galaxy with m > 0, and the evaluation then returns zero
 * probability above a panel top.
 */
static inline void _nc_galaxy_shape_factor_moments_tilt_peek_table (NcGalaxyShapeFactorMomentsTiltPrivate * const self,
                                                                    NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                                    const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                                                    const NcGalaxyShapeFactorMomentsTable **table_out,
                                                                    gdouble *ln_P0_out);

static void
_nc_galaxy_shape_factor_moments_tilt_data_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self =
    nc_galaxy_shape_factor_moments_tilt_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (gsf));
  NcGalaxyShapeFactorMomentsTiltLData *ldata = (NcGalaxyShapeFactorMomentsTiltLData *) data->ldata;
  const gdouble R                            = nc_galaxy_shape_factor_data_get_radius (data);
  const gdouble R_min                        = MAX (0.0, R - self->delta_R);
  const gdouble z_sup                        = z_max;
  gdouble ghat_max                           = 1.0;

  if (!self->range_ready)
  {
    ldata->range_valid = FALSE;

    return;
  }

  if ((R_min <= 0.0) || (z_sup <= self->z_cl))
  {
    /* No background source, or a radius the center's prior can take to the
     * halo itself: nothing to bound. */
    ghat_max = (z_sup <= self->z_cl) ? 0.0 : 1.0;
  }
  else if (!gsl_finite (z_sup))
  {
    /* Redshift support with no upper end: keep the full mesh, since the
     * shear cannot be evaluated at z = inf. */
    ghat_max = 1.0;
  }
  else
  {
    gdouble kappa, g_t;

    /* Two evaluations, under range_lock because the profile copy is shared
     * across galaxies. */
    g_mutex_lock (&self->range_lock);
    kappa = nc_wl_surface_mass_density_convergence (self->smd, self->dp_box, self->cosmo, R_min, z_sup, self->z_cl, self->z_cl);
    g_t   = fabs (nc_wl_surface_mass_density_reduced_shear (self->smd, self->dp_box, self->cosmo, R_min, z_sup, self->z_cl, self->z_cl));
    g_mutex_unlock (&self->range_lock);

    /* kappa >= 1 at the corner: |g| crossed 1 inside the box. A non-finite
     * value bounds nothing. */
    if (gsl_finite (kappa) && gsl_finite (g_t) && (kappa < 1.0))
    {
      const gdouble g_cal = fabs (1.0 + data->m) * g_t + hypot (data->c1, data->c2);

      ghat_max = MIN (g_cal, 1.0);
    }
    else
    {
      ghat_max = 1.0;
    }
  }

  ldata->ghat_max    = ghat_max;
  ldata->range_valid = TRUE;

  /* Force the table to be re-acquired: the range sets its panel count, and
   * the cache is keyed on that. */
  ldata->tab_valid = FALSE;

  /* Build the table now rather than at the first evaluation, so that a
   * prepared dataset evaluates without building anything.
   * NcDataClusterWLFactor runs this step before the redshift grid: the
   * auto-node calibration probes the shape integrand, and with the range
   * already recorded the probe finds this table instead of building a
   * full-mesh one that nothing reads afterwards. */
  {
    NcGalaxyShapePop *pop = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

    if (pop != NULL)
    {
      const NcGalaxyShapeFactorMomentsTable *table;
      gdouble ln_P0;

      _nc_galaxy_shape_factor_moments_tilt_peek_table (self, pop, data,
                                                       data->epsilon_obs_1, data->epsilon_obs_2,
                                                       &table, &ln_P0);
    }
  }
}

static void
_nc_galaxy_shape_factor_moments_tilt_report_failures (NcGalaxyShapeFactorMomentsTiltPrivate * const self,
                                                      const guint n_failed, const gdouble sn, const guint n_panels)
{
  if (G_LIKELY (n_failed == 0))
    return;

  if (G_UNLIKELY (self->strict_solve))
    g_error ("NcGalaxyShapeFactorMomentsTilt: the exact moment-matching Newton "
             "solve failed at %u nodes of a %u-panel table (sigma_nu=%g). The "
             "target moments are closed form, so this means the population "
             "parameters do not describe a distribution on the unit disc, or "
             "the population is so narrow against the noise that the moment "
             "equations are too stiff to solve (seen for sigma_pop ~ 0.01 at "
             "sigma_nu ~ 0.01 under TRACE_DET).",
             n_failed, n_panels, sn);

  if (G_UNLIKELY (g_atomic_int_compare_and_exchange (&self->solve_warned, 0, 1)))
    g_warning ("NcGalaxyShapeFactorMomentsTilt: the exact moment-matching Newton "
               "solve failed at %u nodes of a %u-panel table (sigma_nu=%g); the "
               "affected panel is built from the remaining nodes and is not "
               "accurate. Warned once per instance; read the running total "
               "with nc_galaxy_shape_factor_moments_tilt_get_solve_error_count(). "
               "The count is per distinct table, not per galaxy.",
               n_failed, n_panels, sn);

  g_atomic_int_add (&self->solve_error_count, (gint) n_failed);
}

/*
 * Returns a new reference to the table for (sn, e_rms, n_panels), building
 * it if it is not cached.
 *
 * The lookup and the insertion happen under cache_lock; the build does not.
 * Two threads that miss the same key therefore both build it, and the
 * second to finish discards its own copy for the one already inserted. The
 * reference is taken while the lock is held: a concurrent replacement of
 * the same key would otherwise free the table between the lookup and the
 * ref.
 */
static NcGalaxyShapeFactorMomentsTable *
_nc_galaxy_shape_factor_moments_tilt_acquire_table (NcGalaxyShapeFactorMomentsTiltPrivate * const self,
                                                    NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                    const gdouble sn, const guint n_panels, const gdouble top)
{
  const NcGalaxyShapeFactorMomentsKey key = { sn, data->pop_data->e_rms, n_panels };
  NcGalaxyShapeFactorMomentsTable *table;
  NcGalaxyShapeFactorMomentsTable *built;

  g_mutex_lock (&self->cache_lock);
  table = (NcGalaxyShapeFactorMomentsTable *) g_hash_table_lookup (self->tab_cache, &key);

  if (table != NULL)
    table = _nc_galaxy_shape_factor_moments_table_ref (table);

  g_mutex_unlock (&self->cache_lock);

  if (table != NULL)
    return table;

  {
    NcmSpectral **sp = ncm_memory_pool_get (self->spectral_pool);

    built = _moments_tilt_build_table (*sp, self->ellip_conv, pop, data->pop_data, sn, n_panels, top,
                                       self->accuracy_gate, self->max_degree, self->newton_maxit);
    ncm_memory_pool_return (sp);
  }

  g_atomic_int_inc (&self->table_build_count);
  _nc_galaxy_shape_factor_moments_tilt_report_failures (self, built->n_failed, sn, n_panels);

  g_mutex_lock (&self->cache_lock);
  table = (NcGalaxyShapeFactorMomentsTable *) g_hash_table_lookup (self->tab_cache, &key);

  /* LCOV_EXCL_START: only when another thread built the same table meanwhile. */
  if (table != NULL)
  {
    table = _nc_galaxy_shape_factor_moments_table_ref (table);
    _nc_galaxy_shape_factor_moments_table_unref (built);
  }
  /* LCOV_EXCL_STOP */
  else
  {
    NcGalaxyShapeFactorMomentsKey *key_copy = g_new (NcGalaxyShapeFactorMomentsKey, 1);

    *key_copy = key;
    g_hash_table_insert (self->tab_cache, key_copy, built);
    table = _nc_galaxy_shape_factor_moments_table_ref (built);
    self->tables_generation++;
  }

  g_mutex_unlock (&self->cache_lock);

  return table;
}

/*
 * The two cold paths of peek_table(), kept out of line so that its checks
 * inline into every evaluation. refresh_table() acquires the galaxy's table
 * for the current population generation, range generation and std_noise;
 * refresh_ln_P0() recomputes ln_P0, which additionally depends on the
 * observed radius.
 */
G_GNUC_NO_INLINE static void
_nc_galaxy_shape_factor_moments_tilt_refresh_table (NcGalaxyShapeFactorMomentsTiltPrivate * const self,
                                                    NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorMomentsTiltLData *ldata = (NcGalaxyShapeFactorMomentsTiltLData *) data->ldata;
  const gdouble sn                           = data->std_noise;

  gboolean clamped   = FALSE;
  const guint n_full = _moments_tilt_panels_for_sn (sn, &clamped);
  guint n_panels     = n_full;
  gdouble top        = 1.0;
  NcGalaxyShapeFactorMomentsTable *table;

  /* Quantize the reachable range upward onto the mesh: J panels cover
   * [0, 1 - 2^-J], and J = K restores the full mesh. A galaxy that never
   * saw data_prepare() has no bound recorded and gets the full mesh. */
  if (self->range_holds && ldata->range_valid && (ldata->ghat_max < 1.0))
  {
    const gdouble need = ceil (log2 (1.0 / (1.0 - ldata->ghat_max)));

    n_panels = (guint) CLAMP (need, 1.0, (gdouble) n_full);
    top      = (n_panels == n_full) ? 1.0 : 1.0 - ldexp (1.0, -(gint) n_panels);
  }

  if (G_UNLIKELY (clamped && g_atomic_int_compare_and_exchange (&self->panels_warned, 0, 1)))
    g_warning ("NcGalaxyShapeFactorMomentsTilt: sigma_nu=%g needs more "
               "than the %u panels the layout holds, so the mesh stops there "
               "and the boundary layer at ghat=1 is under-resolved. Warned "
               "once per instance.",
               sn, NC_GALAXY_SHAPE_FACTOR_MOMENTS_MAX_PANELS);

  table = _nc_galaxy_shape_factor_moments_tilt_acquire_table (self, pop, data, sn, n_panels, top);

  g_clear_pointer (&ldata->table, _nc_galaxy_shape_factor_moments_table_unref);
  ldata->table             = table; /* already a reference */
  ldata->table_span        = _nc_galaxy_shape_factor_moments_table_span (table, table->n_panels, NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_PREFETCH_CAP);
  ldata->pop_hash_seen_tab = self->pop_hash;
  ldata->sn_seen           = sn;
  ldata->range_gen_seen    = self->range_gen;
  ldata->n_panels_seen     = n_panels;
  ldata->tab_valid         = TRUE;

  /* sn moved, so a previously cached ln_P0 is stale too. */
  ldata->p0_valid = FALSE;
}

G_GNUC_NO_INLINE static void
_nc_galaxy_shape_factor_moments_tilt_refresh_ln_P0 (NcGalaxyShapeFactorMomentsTiltPrivate * const self,
                                                    NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                    const gdouble R2)
{
  NcGalaxyShapeFactorMomentsTiltLData *ldata = (NcGalaxyShapeFactorMomentsTiltLData *) data->ldata;
  const gdouble sn                           = data->std_noise;

  ldata->ln_P0 = _moments_tilt_ln_P0 (pop, data->pop_data, sqrt (R2), sn);

  ldata->pop_hash_seen_p0 = self->pop_hash;
  ldata->R2_seen          = R2;
  ldata->p0_valid         = TRUE;
}

static inline void
_nc_galaxy_shape_factor_moments_tilt_peek_table (NcGalaxyShapeFactorMomentsTiltPrivate * const self,
                                                 NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                 const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                                 const NcGalaxyShapeFactorMomentsTable **table_out,
                                                 gdouble *ln_P0_out)
{
  NcGalaxyShapeFactorMomentsTiltLData *ldata = (NcGalaxyShapeFactorMomentsTiltLData *) data->ldata;
  const gdouble sn                           = data->std_noise;

  /* R = |epsilon_obs| is read from the arguments, not from
   * data->epsilon_obs_1/2: the two coincide on the fixed-nodes pipeline
   * path but not in general. */
  const gdouble R2 = gsl_pow_2 (epsilon_obs_1) + gsl_pow_2 (epsilon_obs_2);

  if (G_UNLIKELY (!ldata->tab_valid || (ldata->pop_hash_seen_tab != self->pop_hash) ||
                  (ldata->sn_seen != sn) || (ldata->range_gen_seen != self->range_gen)))
    _nc_galaxy_shape_factor_moments_tilt_refresh_table (self, pop, data);

  if (G_UNLIKELY (!ldata->p0_valid || (ldata->pop_hash_seen_p0 != self->pop_hash) ||
                  (ldata->R2_seen != R2)))
    _nc_galaxy_shape_factor_moments_tilt_refresh_ln_P0 (self, pop, data, R2);

  *table_out = ldata->table;
  *ln_P0_out = ldata->ln_P0;
}

/*
 * Returns ln P_0 + lambda_1 x + lambda_2 x^2 + lambda_3 y^2 - W, or its
 * exponential when @want_log is FALSE. (x, y) is eps_obs rotated by
 * -arg(g), and the table is read at ghat = min(|g|, 1/|g|); the fold is
 * exact, since the marginal is invariant under g -> 1/g and the table is a
 * function of ghat alone.
 *
 * The upper-end check runs before the panel index is formed. The index
 * clamps to the last panel, so a ghat above the table's end would otherwise
 * be evaluated there, outside that panel's own interval, where a Chebyshev
 * series diverges. Such a point is outside the range the table was built
 * for: it is counted, warned about once per instance, and returns zero
 * probability (GSL_NEGINF for the log).
 *
 * There is no ceiling guard on the density. A truncated tilt needs one, to
 * catch a lambda so negative that the model becomes a spurious spike; the
 * tilt here is solved, not truncated, and the density integrates to one by
 * construction of W.
 */
static gdouble
_nc_galaxy_shape_factor_moments_tilt_eval (NcGalaxyShapeFactorMomentsTiltPrivate * const self,
                                           const NcGalaxyShapeFactorMomentsTable *table, const gdouble ln_P0,
                                           const gdouble g_1, const gdouble g_2,
                                           const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                           const gboolean want_log)
{
  const gdouble g_mag  = sqrt (g_1 * g_1 + g_2 * g_2);
  const gdouble cos_pg = (g_mag > 0.0) ? g_1 / g_mag : 1.0;
  const gdouble sin_pg = (g_mag > 0.0) ? g_2 / g_mag : 0.0;
  const gdouble x      = epsilon_obs_1 * cos_pg + epsilon_obs_2 * sin_pg;
  const gdouble y      = -epsilon_obs_1 * sin_pg + epsilon_obs_2 * cos_pg;
  const gdouble ghat   = (g_mag <= 1.0) ? g_mag : 1.0 / g_mag;
  guint j;
  gdouble lnP;

  /* lambda(0) = 0 and W(0) = 0 hold exactly, and the zero-shear node is
   * never solved. A truncated series is not exact at its nodes, so the
   * zero-shear case is returned directly rather than through the
   * interpolant: at g = 0 this class is the exact zero-shear marginal, and
   * the identity holds bitwise. */
  if (g_mag == 0.0)
    return want_log ? ln_P0 : exp (ln_P0);

  if (G_UNLIKELY (ghat > table->top + NC_GALAXY_SHAPE_FACTOR_MOMENTS_TOP_TOL))
  {
    g_atomic_int_inc (&self->range_error_count);

    if (G_UNLIKELY (g_atomic_int_compare_and_exchange (&self->range_warned, 0, 1)))
      g_warning ("NcGalaxyShapeFactorMomentsTilt: ghat=%g lies above the "
                 "%g this galaxy's table covers, so the reachable range it was "
                 "built for was violated. The point is reported as zero "
                 "probability rather than extrapolated. Warned once per "
                 "instance; read the running total with "
                 "nc_galaxy_shape_factor_moments_tilt_get_range_error_count().",
                 ghat, table->top);

    return want_log ? GSL_NEGINF : 0.0;
  }

  j   = _nc_galaxy_shape_factor_moments_panel_index (table, ghat);
  lnP = ln_P0 + _moments_tilt_clenshaw_fused (table, j, _nc_galaxy_shape_factor_moments_panel_arg (table, j, ghat), x, x * x, y * y);

  return want_log ? lnP : exp (lnP);
}

static void
_nc_galaxy_shape_factor_moments_tilt_data_prefetch (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const guint stage)
{
  const NcGalaxyShapeFactorMomentsTiltLData *ldata = (const NcGalaxyShapeFactorMomentsTiltLData *) data->ldata;

  if (stage == 1)
    ncm_prefetch_span (ldata, sizeof (NcGalaxyShapeFactorMomentsTiltLData));
  else if (stage == 2)
    ncm_prefetch_span (ldata->table, ldata->table_span);
}

static gdouble
_nc_galaxy_shape_factor_moments_tilt_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self =
    nc_galaxy_shape_factor_moments_tilt_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (gsf));
  const NcGalaxyShapeFactorMomentsTable *table;
  gdouble ln_P0;

  _nc_galaxy_shape_factor_moments_tilt_peek_table (self, pop, data, epsilon_obs_1, epsilon_obs_2, &table, &ln_P0);

  return _nc_galaxy_shape_factor_moments_tilt_eval (self, table, ln_P0, g_1, g_2, epsilon_obs_1, epsilon_obs_2, FALSE);
}

static gdouble
_nc_galaxy_shape_factor_moments_tilt_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self =
    nc_galaxy_shape_factor_moments_tilt_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (gsf));
  const NcGalaxyShapeFactorMomentsTable *table;
  gdouble ln_P0;

  _nc_galaxy_shape_factor_moments_tilt_peek_table (self, pop, data, epsilon_obs_1, epsilon_obs_2, &table, &ln_P0);

  return _nc_galaxy_shape_factor_moments_tilt_eval (self, table, ln_P0, g_1, g_2, epsilon_obs_1, epsilon_obs_2, TRUE);
}

static gchar *
_nc_galaxy_shape_factor_moments_tilt_get_desc (NcGalaxyShapeFactor *gsf)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self =
    nc_galaxy_shape_factor_moments_tilt_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (gsf));
  gchar *parent_desc = NC_GALAXY_SHAPE_FACTOR_CLASS (nc_galaxy_shape_factor_moments_tilt_parent_class)->get_desc (gsf);
  gchar *desc        = g_strdup_printf ("%s, gate=%g, max_degree=%u, newton_maxit=%u", parent_desc, self->accuracy_gate, self->max_degree, self->newton_maxit);

  g_free (parent_desc);

  return desc;
}

static void
_nc_galaxy_shape_factor_moments_tilt_dispose (GObject *object)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self =
    nc_galaxy_shape_factor_moments_tilt_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (object));

  g_clear_pointer (&self->tab_cache, g_hash_table_unref);
  nc_hicosmo_clear (&self->cosmo);
  nc_wl_surface_mass_density_clear (&self->smd);
  nc_halo_density_profile_clear (&self->dp_box);
  nc_galaxy_shape_pop_clear (&self->pop_ref);
  g_clear_pointer (&self->pending_tables, ncm_obj_array_unref);
  g_clear_pointer (&self->tables_stamp, g_free);
  g_clear_pointer (&self->shared_tables_cache, ncm_obj_array_unref);

  if (self->spectral_pool != NULL)
  {
    ncm_memory_pool_free (self->spectral_pool, TRUE);
    self->spectral_pool = NULL;
  }

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_moments_tilt_parent_class)->dispose (object);
}

static void
_nc_galaxy_shape_factor_moments_tilt_finalize (GObject *object)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self =
    nc_galaxy_shape_factor_moments_tilt_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (object));

  g_mutex_clear (&self->cache_lock);
  g_mutex_clear (&self->range_lock);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_moments_tilt_parent_class)->finalize (object);
}

static void
nc_galaxy_shape_factor_moments_tilt_class_init (NcGalaxyShapeFactorMomentsTiltClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_shape_factor_moments_tilt_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_moments_tilt_get_property;
  object_class->dispose      = &_nc_galaxy_shape_factor_moments_tilt_dispose;
  object_class->finalize     = &_nc_galaxy_shape_factor_moments_tilt_finalize;
  object_class->constructed  = &_nc_galaxy_shape_factor_moments_tilt_constructed;

  /**
   * NcGalaxyShapeFactorMomentsTilt:accuracy-gate:
   *
   * Target accuracy of the interpolant, in nats of model mass
   * $|\ln Z|$. Each panel is truncated to the smallest degree whose signed
   * coefficient tail is below one tenth of this value; the factor of ten is
   * the margin for the tail being an estimate of the error and not a bound
   * on it. The default of $10^{-3}$ is about a hundred times below the
   * irreducible error of the tilt family itself, so it does not limit the
   * accuracy of a fit.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ACCURACY_GATE,
                                   g_param_spec_double ("accuracy-gate",
                                                        "Accuracy gate",
                                                        "Target |ln Z| accuracy of the interpolant",
                                                        1.0e-12, 1.0, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorMomentsTilt:max-degree:
   *
   * Largest local degree a panel may reach. The adaptive refinement doubles
   * from four, so the degrees tried are 4, 8, 16. The degree selected is
   * four to six on almost every panel, since the dyadic grading is
   * self-similar and each panel has the same shape in its own rescaled
   * coordinate, so the default of 16 is not reached.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_DEGREE,
                                   g_param_spec_uint ("max-degree",
                                                      "Maximum local degree",
                                                      "Largest Chebyshev degree a single panel may reach",
                                                      4, 64, 16,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorMomentsTilt:newton-maxit:
   *
   * Largest number of Newton iterations per pass of the solve for the tilt
   * parameters at one table node. A solve that reaches it returns the last
   * iterate with its own log-normalization, and counts as converged only if
   * its residual is below the acceptance threshold $10^{-8}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_NEWTON_MAXIT,
                                   g_param_spec_uint ("newton-maxit",
                                                      "Maximum Newton iterations",
                                                      "Largest number of Newton iterations per solve pass at a table node",
                                                      1, 1000, NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_NEWTON_MAXIT,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorMomentsTilt:restrict-range:
   *
   * Whether to build each galaxy's table only over the range of folded
   * shear its line of sight can reach inside the parameter box, instead of
   * over all of $[0,1]$.
   *
   * The mesh exists to resolve the boundary layer at $\hat g=1$, which a
   * galaxy reaches only if its line of sight can cross the tangential
   * critical curve. On a real cluster most do not, and the bound costs two
   * shear evaluations per galaxy against the dozens of Newton solves of the
   * panels it removes. It needs nc_galaxy_shape_factor_data_prepare() to
   * have run for the galaxy; without that the full mesh is used, which is
   * slower and equally correct.
   *
   * The bound is a property of the prior box, so it does not change as the
   * sampler moves, and it removes panels only if that prior is tight: a
   * halo center free over the whole field makes every galaxy reach
   * $\hat g = 1$, and every table is then the full mesh.
   *
   * The bound covers the mass (log10MDelta), the halo center (ra, dec) and
   * the shape population. If any other parameter is free in the model set
   * at prepare time, such as the concentration, a cosmological parameter,
   * the halo redshift or a redshift-model parameter, the bound does not
   * hold and every galaxy gets the full mesh, as with this property set to
   * %FALSE.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RESTRICT_RANGE,
                                   g_param_spec_boolean ("restrict-range",
                                                         "Restrict to the reachable range",
                                                         "Whether to build each table only over the reachable folded shear",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorMomentsTilt:strict-solve:
   *
   * Whether a Newton failure during a table build is a fatal error instead
   * of a warning counted in
   * nc_galaxy_shape_factor_moments_tilt_get_solve_error_count(). %FALSE by
   * default, so that an exploratory run reports and continues.
   *
   * Installed without %G_PARAM_CONSTRUCT, which would overwrite the value
   * set in the instance initializer.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_STRICT_SOLVE,
                                   g_param_spec_boolean ("strict-solve",
                                                         "Strict solve",
                                                         "Whether a Newton failure aborts instead of warning",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorMomentsTilt:tables:
   *
   * The interpolation tables this instance has built, as an
   * #NcmObjArray of #NcmVector, one per table, each carrying its own key and
   * layout. Serializing an instance includes them, so a later run started
   * from the saved instance builds no table.
   *
   * They are adopted only if #NcGalaxyShapeFactorMomentsTilt:tables-stamp
   * matches the current population and build settings. A mismatch discards
   * them and rebuilds, which changes only the time taken.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_TABLES,
                                   g_param_spec_boxed ("tables",
                                                       "Interpolation tables",
                                                       "Built interpolation tables, one NcmVector each",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorMomentsTilt:tables-stamp:
   *
   * What #NcGalaxyShapeFactorMomentsTilt:tables was built against: the
   * format tag, the build settings, and the serialized population.
   *
   * It is built from values rather than from the parameter-key counters the
   * in-process caches use: those counters are per-process and carry no
   * meaning across a save and load.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_TABLES_STAMP,
                                   g_param_spec_string ("tables-stamp",
                                                        "Table provenance stamp",
                                                        "What the stored tables were built against",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  gsf_class->data_init = &_nc_galaxy_shape_factor_moments_tilt_data_init;

  /* The table cache and the shared profile copy are locked (see the
   * private struct); everything else a galaxy touches is its own. */
  gsf_class->prepare          = &_nc_galaxy_shape_factor_moments_tilt_prepare;
  gsf_class->data_prepare     = &_nc_galaxy_shape_factor_moments_tilt_data_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_moments_tilt_eval_marginal;
  gsf_class->data_prefetch    = &_nc_galaxy_shape_factor_moments_tilt_data_prefetch;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_moments_tilt_eval_ln_marginal;
  gsf_class->get_desc         = &_nc_galaxy_shape_factor_moments_tilt_get_desc;
}

/**
 * nc_galaxy_shape_factor_moments_tilt_new:
 * @ellip_conv: the ellipticity convention #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorMomentsTilt.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorMomentsTilt
 */
NcGalaxyShapeFactorMomentsTilt *
nc_galaxy_shape_factor_moments_tilt_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  NcGalaxyShapeFactorMomentsTilt *gsfmt = g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_MOMENTS_TILT,
                                                        "ellip-conv", ellip_conv,
                                                        NULL);

  return gsfmt;
}

/**
 * nc_galaxy_shape_factor_moments_tilt_ref:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 *
 * Increases the reference count of @gsfmt by one.
 *
 * Returns: (transfer full): @gsfmt
 */
NcGalaxyShapeFactorMomentsTilt *
nc_galaxy_shape_factor_moments_tilt_ref (NcGalaxyShapeFactorMomentsTilt *gsfmt)
{
  return g_object_ref (gsfmt);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_free:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 *
 * Decreases the reference count of @gsfmt by one.
 *
 */
void
nc_galaxy_shape_factor_moments_tilt_free (NcGalaxyShapeFactorMomentsTilt *gsfmt)
{
  g_object_unref (gsfmt);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_clear:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 *
 * If *@gsfmt is not %NULL, decreases its reference count by one and sets
 * *@gsfmt to %NULL.
 *
 */
void
nc_galaxy_shape_factor_moments_tilt_clear (NcGalaxyShapeFactorMomentsTilt **gsfmt)
{
  g_clear_object (gsfmt);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_get_solve_error_count:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 *
 * Number of Chebyshev nodes at which the Newton solve failed since the last
 * reset. One table serves every galaxy sharing its key, so this counts
 * nodes of distinct tables and not galaxies affected. Under a parallel
 * nc_data_cluster_wl_factor_data_prepare(), two threads can build the same
 * table at once and one copy is discarded; the discarded build's failures
 * are counted too, so the total can exceed the count over the tables kept.
 *
 * Returns: the running count of failed nodes
 */
guint
nc_galaxy_shape_factor_moments_tilt_get_solve_error_count (NcGalaxyShapeFactorMomentsTilt *gsfmt)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);

  return (guint) g_atomic_int_get (&self->solve_error_count);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_reset_solve_error_count:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 *
 * Zeroes the failed-node counter and re-arms the one-shot warning.
 *
 */
void
nc_galaxy_shape_factor_moments_tilt_reset_solve_error_count (NcGalaxyShapeFactorMomentsTilt *gsfmt)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);

  g_atomic_int_set (&self->solve_error_count, 0);
  g_atomic_int_set (&self->solve_warned, 0);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_get_range_error_count:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 *
 * Number of evaluations whose folded shear was above the range its table
 * covers since the last reset. A non-zero value means a table was built
 * for a smaller range than the sampler visited, and those evaluations
 * returned zero probability instead of an extrapolated value.
 *
 * Returns: the running count of out-of-range evaluations
 */
guint
nc_galaxy_shape_factor_moments_tilt_get_range_error_count (NcGalaxyShapeFactorMomentsTilt *gsfmt)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);

  return (guint) g_atomic_int_get (&self->range_error_count);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_reset_range_error_count:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 *
 * Zeroes the out-of-range counter and re-arms the one-shot warning.
 *
 */
void
nc_galaxy_shape_factor_moments_tilt_reset_range_error_count (NcGalaxyShapeFactorMomentsTilt *gsfmt)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);

  g_atomic_int_set (&self->range_error_count, 0);
  g_atomic_int_set (&self->range_warned, 0);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_get_table_build_count:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 *
 * Number of tables built since the last reset, not counting tables served
 * from the instance cache. A serial run over the same population and the
 * same per-galaxy shape dispersions builds the same number, and a run
 * reusing stored tables builds none. A parallel
 * nc_data_cluster_wl_factor_data_prepare() can build a table twice when
 * two threads need it at once, keeping one; both builds are counted, so
 * the count then exceeds the number of distinct tables and varies from run
 * to run, while the results do not.
 *
 * Returns: the running count of table builds
 */
guint
nc_galaxy_shape_factor_moments_tilt_get_table_build_count (NcGalaxyShapeFactorMomentsTilt *gsfmt)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);

  return (guint) g_atomic_int_get (&self->table_build_count);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_reset_table_build_count:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 *
 * Zeroes the table-build counter.
 *
 */
void
nc_galaxy_shape_factor_moments_tilt_reset_table_build_count (NcGalaxyShapeFactorMomentsTilt *gsfmt)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);

  g_atomic_int_set (&self->table_build_count, 0);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_register_shared:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 * @ser: a #NcmSerialize
 *
 * Registers the #NcmVector of each table of this instance's current table
 * set (its #NcGalaxyShapeFactorMomentsTilt:tables property) as a named
 * instance in @ser, under names derived from the address of @gsfmt.
 *
 * The tables are functions of (population, shape dispersion, panel count)
 * only and are read-only once built (see the class documentation). A later
 * ncm_serialize_to_variant() or ncm_serialize_dup_obj() of @gsfmt through
 * @ser finds each vector already named and emits a reference instead of
 * re-serializing and rebuilding it, for every caller sharing @ser,
 * including across ncm_serialize_reset(ser, TRUE) (autosave only), which
 * #NcmFitESMCMC and #NcmFitMC run between per-worker fit clones.
 * #NcDataClusterWLFactor's register_shared() does the same for its "obs"
 * property. Without this registration every worker thread's clone
 * deep-copies every cached table, which for a large galaxy catalog
 * dominates per-thread setup time.
 *
 * Registering does not build any table: an instance with an empty cache
 * registers nothing, and a later ncm_serialize_to_variant() serializes
 * whatever tables exist by then, so callers populate the cache first (for
 * example by one serial prepare).
 *
 * May be called repeatedly, including as the cache grows: registering a
 * table already registered under the same name changes nothing, and new
 * tables are added.
 */
void
nc_galaxy_shape_factor_moments_tilt_register_shared (NcGalaxyShapeFactorMomentsTilt *gsfmt, NcmSerialize *ser)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);
  NcmObjArray *tables                                = _nc_galaxy_shape_factor_moments_tilt_get_shared_tables (self);
  guint i;

  for (i = 0; i < ncm_obj_array_len (tables); i++)
  {
    GObject *vec = ncm_obj_array_peek (tables, i);

    /* NcmSerialize's named-instance syntax accepts only [A-Za-z0-9:] inside
     * the brackets (see its is_named_regex): no "0x" prefix, no
     * underscore. */
    gchar *name = g_strdup_printf ("ncGsfMomentsTiltTab:%" G_GSIZE_FORMAT ":%u", (gsize) gsfmt, i);

    ncm_serialize_set (ser, vec, name, TRUE);
    g_free (name);
  }

  ncm_obj_array_unref (tables);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_exact_moments:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 * @pop: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapeFactorData
 * @ghat: the folded shear $\hat g=\min(g,1/g)\in[0,1]$
 * @mu: (out): $\mathbb E[x]$
 * @Ex2: (out): $\mathbb E[x^2]$
 * @Ey2: (out): $\mathbb E[y^2]$
 *
 * Computes the closed-form target moments of the exact marginal at @ghat,
 * including the noise variance, in this instance's ellipticity convention.
 * Exposed for validation. At $\hat g = 1$ both conventions' maps send every
 * source to 1, so there the moments are exactly
 * $(1, 1 + \sigma_\nu^2, \sigma_\nu^2)$; under TRACE_DET
 * $\mathrm{E}[x] = \hat g$ at every shear.
 *
 */
void
nc_galaxy_shape_factor_moments_tilt_exact_moments (NcGalaxyShapeFactorMomentsTilt *gsfmt, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble ghat, gdouble *mu, gdouble *Ex2, gdouble *Ey2)
{
  const guint nq    = NC_GALAXY_SHAPE_FACTOR_MOMENTS_MOMENT_NNODES;
  GArray *r_arr     = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), nq);
  GArray *p_arr     = NULL;
  const gdouble sn2 = data->std_noise * data->std_noise;
  gdouble t[3];

  g_array_set_size (r_arr, nq);
  _nc_galaxy_shape_factor_moments_exact_t (nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfmt)), pop, data->pop_data, ghat, sn2, r_arr, &p_arr, t);

  *mu  = t[0];
  *Ex2 = t[1];
  *Ey2 = t[2];

  g_array_unref (r_arr);

  if (p_arr != NULL)
    g_array_unref (p_arr);
}

/**
 * nc_galaxy_shape_factor_moments_tilt_eval_tilt:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 * @pop: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: first reduced-shear component
 * @g_2: second reduced-shear component
 * @lambda_1: (out): $\lambda_1$
 * @lambda_2: (out): $\lambda_2$
 * @lambda_3: (out): $\lambda_3$
 * @W: (out): the exact log-normalizer $W$
 *
 * Computes the interpolated tilt at @g_1, @g_2, building this galaxy's
 * table if needed. Exposed for validation. Zero shear returns zeros, which
 * is the exact value; a folded shear above the table's upper end also
 * returns zeros, without counting or warning.
 *
 */
void
nc_galaxy_shape_factor_moments_tilt_eval_tilt (NcGalaxyShapeFactorMomentsTilt *gsfmt, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, gdouble *lambda_1, gdouble *lambda_2, gdouble *lambda_3, gdouble *W)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);
  const gdouble g_mag                                = hypot (g_1, g_2);
  const gdouble ghat                                 = (g_mag <= 1.0) ? g_mag : 1.0 / g_mag;
  const NcGalaxyShapeFactorMomentsTable *table;
  gdouble ln_P0, lw[4];
  guint j;

  _nc_galaxy_shape_factor_moments_tilt_peek_table (self, pop, data, 0.0, 0.0, &table, &ln_P0);

  /* Zero shear returns lambda = 0, W = 0, the exact value; a folded shear
   * above the table's upper end returns the same, since the table has no
   * value there. */
  if ((g_mag == 0.0) || (ghat > table->top + NC_GALAXY_SHAPE_FACTOR_MOMENTS_TOP_TOL))
  {
    *lambda_1 = *lambda_2 = *lambda_3 = *W = 0.0;

    return;
  }

  j = _nc_galaxy_shape_factor_moments_panel_index (table, ghat);
  _nc_galaxy_shape_factor_moments_clenshaw (table, j, _nc_galaxy_shape_factor_moments_panel_arg (table, j, ghat), lw);

  *lambda_1 = lw[0];
  *lambda_2 = lw[1];
  *lambda_3 = lw[2];
  *W        = lw[3];
}

/**
 * nc_galaxy_shape_factor_moments_tilt_peek_layout:
 * @gsfmt: a #NcGalaxyShapeFactorMomentsTilt
 * @pop: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapeFactorData
 * @n_panels: (out): number of panels
 * @top: (out): upper end of the mesh
 * @degrees: (out) (element-type gdouble) (transfer full): per-panel degree
 *
 * Returns the mesh of this galaxy's table, building the table if needed.
 * Exposed so that a validation run can report the node count of a
 * configuration and check the adaptive degree against the measured
 * profile.
 *
 */
void
nc_galaxy_shape_factor_moments_tilt_peek_layout (NcGalaxyShapeFactorMomentsTilt *gsfmt, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, guint *n_panels, gdouble *top, GArray **degrees)
{
  NcGalaxyShapeFactorMomentsTiltPrivate * const self = nc_galaxy_shape_factor_moments_tilt_get_instance_private (gsfmt);
  const NcGalaxyShapeFactorMomentsTable *table;
  gdouble ln_P0;
  guint j;

  _nc_galaxy_shape_factor_moments_tilt_peek_table (self, pop, data, 0.0, 0.0, &table, &ln_P0);

  *n_panels = table->n_panels;
  *top      = table->top;
  *degrees  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), table->n_panels);

  for (j = 0; j < table->n_panels; j++)
  {
    const gdouble d = (gdouble) table->deg[j];

    g_array_append_val (*degrees, d);
  }
}

