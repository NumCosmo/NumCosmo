/***************************************************************************
 *            nc_xcor.h
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcXCor:
 *
 * Angular auto- and cross-spectra.
 *
 * The angular power spectrum between observables $A$ and $B$ with kernels $W^A$ and $W^B$ is
 * \begin{equation}
 * C_{\ell}^{AB} = \int dz W^A(z) \int dz^\prime W^B (z^\prime) \times \int dk \frac{2}{\pi} k^2 P(k, z, z^\prime) j_{\ell}(k\chi(z)) j_{\ell} (k\chi(z^\prime)).
 * \end{equation}
 * In the Limber approximation, it reduces to
 * \begin{equation}
 * C_{\ell}^{AB} = \int_0^{z_*} dz \frac{H(z)}{c \chi^2(z)} W^A(z) W^B (z) P\left(k = \frac{\ell +1/2}{\chi(z)} , z \right),
 * \end{equation}
 * where $P\left(k = \frac{\ell +1/2}{\chi(z)} , z \right)$ is the power spectrum (a #NcmPowspec) at redshift $z$ and $chi(z)$ the comoving distance (a #NcDistance).
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/nc_xcor.h"
#include "nc/xcor/nc_xcor_priv.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_DISTANCE,
  PROP_MATTER_POWER_SPECTRUM,
  PROP_METH,
  PROP_CLOSURE_TYPE,
  PROP_RELTOL,
  PROP_ELL_BATCH_SIZE,
};

G_DEFINE_TYPE (NcXcor, nc_xcor, G_TYPE_OBJECT)

static void
nc_xcor_init (NcXcor *xc)
{
  xc->ps             = NULL;
  xc->dist           = NULL;
  xc->RH             = 0.0;
  xc->meth           = NC_XCOR_METHOD_LIMBER_Z_GSL;
  xc->closure_type   = NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV;
  xc->ell_batch_size = 8;
}

static void
_nc_xcor_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcor *xc = NC_XCOR (object);

  g_return_if_fail (NC_IS_XCOR (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      nc_distance_clear (&xc->dist);
      xc->dist = g_value_dup_object (value);
      break;
    case PROP_MATTER_POWER_SPECTRUM:
      ncm_powspec_clear (&xc->ps);
      xc->ps = g_value_dup_object (value);
      break;
    case PROP_METH:
      xc->meth = g_value_get_enum (value);
      break;
    case PROP_CLOSURE_TYPE:
      nc_xcor_set_closure_type (xc, g_value_get_enum (value));
      break;
    case PROP_RELTOL:
      nc_xcor_set_reltol (xc, g_value_get_double (value));
      break;
    case PROP_ELL_BATCH_SIZE:
      nc_xcor_set_ell_batch_size (xc, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcor *xc = NC_XCOR (object);

  g_return_if_fail (NC_IS_XCOR (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      g_value_set_object (value, xc->dist);
      break;
    case PROP_MATTER_POWER_SPECTRUM:
      g_value_set_object (value, xc->ps);
      break;
    case PROP_METH:
      g_value_set_enum (value, xc->meth);
      break;
    case PROP_CLOSURE_TYPE:
      g_value_set_enum (value, nc_xcor_get_closure_type (xc));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_xcor_get_reltol (xc));
      break;
    case PROP_ELL_BATCH_SIZE:
      g_value_set_uint (value, nc_xcor_get_ell_batch_size (xc));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_dispose (GObject *object)
{
  NcXcor *xc = NC_XCOR (object);

  nc_distance_clear (&xc->dist);
  ncm_powspec_clear (&xc->ps);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_parent_class)->dispose (object);
}

static void
_nc_xcor_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_parent_class)->finalize (object);
}

static void
nc_xcor_class_init (NcXcorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  /*GObjectClass* parent_class = G_OBJECT_CLASS (klass); */

  object_class->set_property = &_nc_xcor_set_property;
  object_class->get_property = &_nc_xcor_get_property;
  object_class->dispose      = &_nc_xcor_dispose;
  object_class->finalize     = &_nc_xcor_finalize;

  /**
   * NcXcor:distance:
   *
   * This property keeps the distance object.
   */
  g_object_class_install_property (object_class,
                                   PROP_DISTANCE,
                                   g_param_spec_object ("distance",
                                                        NULL,
                                                        "Distance.",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcor:power-spec:
   *
   * This property keeps the matter power spectrum object.
   */
  g_object_class_install_property (object_class,
                                   PROP_MATTER_POWER_SPECTRUM,
                                   g_param_spec_object ("power-spec",
                                                        NULL,
                                                        "Matter power spectrum.",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcor:meth:
   *
   * This property keeps the method enumerator to compute the integrals.
   */
  g_object_class_install_property (object_class,
                                   PROP_METH,
                                   g_param_spec_enum ("meth",
                                                      NULL,
                                                      "Method.",
                                                      NC_TYPE_XCOR_METHOD,
                                                      NC_XCOR_METHOD_LIMBER_Z_GSL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcor:closure-type:
   *
   * How each kernel represents its sampled $W_\ell(k)$. See
   * #NcXcorKernelClosure.
   *
   * The choice is made here rather than on #NcXcorKernel because it is the
   * computation, not the kernel, that should be switchable as a whole: the
   * two representations are alternative fits to the same sampled function, and
   * comparing them is only meaningful when every kernel in a run uses the
   * same one. A pair may still be mixed -- %NC_XCOR_METHOD_KERNEL_EXACT
   * integrates a spline against a panel set exactly, on the common refinement
   * of the two breakpoint sets -- which is what a pair straddling the Limber
   * threshold produces whatever this is set to.
   *
   * Both closures sample the same function over the same domain and differ
   * only in what is fitted to it, so a computation can be switched over and
   * the two compared directly. Measured against the certified Arb C_ell
   * table (43 entries, 17 pairs): the Chebyshev closure is closer in 36 of
   * 43, its median deviation is 4.5x smaller at every tolerance rung, and it
   * has no catastrophic regime -- the spline at loose tolerance returns the
   * wrong sign at 31x the pair scale on a far-separated pair. On a
   * tomographic workload through #NcXcorSolver (5 Gaussian bins, 15 pairs,
   * $\ell \le 60$, one thread) it costs 1.4x at reltol $10^{-4}$ and 1.13x at
   * $10^{-6}$, for a median deviation 4 to 5 orders smaller; on hard-edged
   * cluster top-hats it costs 1.8x and the spline does not converge at all,
   * the two disagreeing by $10^{-4}$ of the pair scale with both asked for
   * $10^{-8}$.
   *
   * That is why the default is %NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV. The spline
   * stays as the independent cross-check: a deviation both closures share at
   * every tolerance is a wrong reference, not a bad fit, and that diagnostic
   * needs two structurally different representations. It is also the cheaper
   * one where the requested tolerance, rather than the representation, is
   * what binds.
   *
   * It applies to the non-Limber closure only. Under Limber each multipole is
   * supported on its own band and zero outside it, so the block's window
   * carries a step per multipole; a Chebyshev series converges on the
   * non-Limber kernel because $W_\ell(k)$ is entire in $k$, and a step is not.
   * Multipoles taken under Limber keep the spline closure whatever this is set
   * to.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_CLOSURE_TYPE,
                                   g_param_spec_enum ("closure-type",
                                                      NULL,
                                                      "Representation used for the k-space closures.",
                                                      NC_TYPE_XCOR_KERNEL_CLOSURE,
                                                      NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcor:reltol:
   *
   * This property keeps the relative tolerance.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance.",
                                                        GSL_DBL_EPSILON, 1.0e-1, NC_XCOR_PRECISION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcor:ell-batch-size:
   *
   * This property keeps the multipole batch size for the kernel-space block
   * methods. Tuned for 8 or 16; #NC_XCOR_KERNEL_MAX_ELL_BLOCK is only a hard
   * safety cap. See nc_xcor_set_ell_batch_size().
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL_BATCH_SIZE,
                                   g_param_spec_uint ("ell-batch-size",
                                                      NULL,
                                                      "Multipole batch size for the kernel-space block methods.",
                                                      1, NC_XCOR_KERNEL_MAX_ELL_BLOCK, 8,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_xcor_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @meth: a #NcXcorMethod to compute the integrals
 *
 * Two methods are available to compute integrals: independent GSL numerical integration or vector integration using Sundials's CVode algorithm.
 *
 * Returns: (transfer full): a #NcXcor
 *
 */
NcXcor *
nc_xcor_new (NcDistance *dist, NcmPowspec *ps, NcXcorMethod meth)
{
  return g_object_new (NC_TYPE_XCOR,
                       "distance", dist,
                       "power-spec", ps,
                       "meth", meth,
                       NULL);
}

/**
 * nc_xcor_ref:
 * @xc: a #NcXcor
 *
 * Increments the reference count of @xc.
 *
 * Returns: (transfer full): @xc
 */
NcXcor *
nc_xcor_ref (NcXcor *xc)
{
  return g_object_ref (xc);
}

/**
 * nc_xcor_free:
 * @xc: a #NcXcor
 *
 * Decrements the reference count of @xc, and frees it if the count reaches 0.
 *
 */
void
nc_xcor_free (NcXcor *xc)
{
  g_object_unref (xc);
}

/**
 * nc_xcor_clear:
 * @xc: a #NcXcor
 *
 * If *@xc is not %NULL, decrements the reference count of @xc, and frees it if the
 * count reaches 0.
 *
 */
void
nc_xcor_clear (NcXcor **xc)
{
  g_clear_object (xc);
}

/**
 * nc_xcor_get_meth:
 * @xc: a #NcXcor
 *
 * Returns: the #NcXcorMethod used by @xc
 */
NcXcorMethod
nc_xcor_get_meth (NcXcor *xc)
{
  return xc->meth;
}

/**
 * nc_xcor_set_closure_type:
 * @xc: a #NcXcor
 * @closure_type: a #NcXcorKernelClosure
 *
 * Sets #NcXcor:closure-type. Read when a closure is built, so one already built
 * keeps the representation it was built with.
 *
 */
void
nc_xcor_set_closure_type (NcXcor *xc, NcXcorKernelClosure closure_type)
{
  xc->closure_type = closure_type;
}

/**
 * nc_xcor_get_closure_type:
 * @xc: a #NcXcor
 *
 * Returns: the representation used for the k-space closures. See
 * #NcXcor:closure-type.
 */
NcXcorKernelClosure
nc_xcor_get_closure_type (NcXcor *xc)
{
  return xc->closure_type;
}

/**
 * nc_xcor_set_reltol:
 * @xc: a #NcXcor
 * @reltol: a relative tolerance
 *
 * Sets the relative tolerance of @xc.
 *
 */
void
nc_xcor_set_reltol (NcXcor *xc, const gdouble reltol)
{
  xc->reltol = reltol;
}

/**
 * nc_xcor_get_reltol:
 * @xc: a #NcXcor
 *
 * Returns: the relative tolerance of @xc
 */
gdouble
nc_xcor_get_reltol (NcXcor *xc)
{
  return xc->reltol;
}

/**
 * nc_xcor_set_ell_batch_size:
 * @xc: a #NcXcor
 * @ell_batch_size: multipole batch size
 *
 * Sets the multipole batch size used by the kernel-space block methods
 * (%NC_XCOR_METHOD_KERNEL_CUBATURE, %NC_XCOR_METHOD_KERNEL_EXACT and
 * %NC_XCOR_METHOD_KERNEL_GSL_BLOCK): each batch builds one k-space closure per
 * kernel and shares it across the whole batch. It does not reach
 * %NC_XCOR_METHOD_KERNEL_GSL, which fits one closure per multipole.
 *
 * The Levin machinery is tuned for 8 (the default) or 16; wider batches are
 * counterproductive, not faster. #NC_XCOR_KERNEL_MAX_ELL_BLOCK is a hard
 * safety cap on a single nc_xcor_kernel_get_eval_vectorized() call, not a
 * setting to aim for.
 *
 */
void
nc_xcor_set_ell_batch_size (NcXcor *xc, const guint ell_batch_size)
{
  g_return_if_fail ((ell_batch_size > 0) && (ell_batch_size <= (guint) NC_XCOR_KERNEL_MAX_ELL_BLOCK));

  xc->ell_batch_size = ell_batch_size;
}

/**
 * nc_xcor_get_ell_batch_size:
 * @xc: a #NcXcor
 *
 * Returns: the multipole batch size
 */
guint
nc_xcor_get_ell_batch_size (NcXcor *xc)
{
  return xc->ell_batch_size;
}

/**
 * nc_xcor_prepare:
 * @xc: a #NcXcor
 * @cosmo: a #NcHICosmo
 *
 * Prepares @xc for computation.
 *
 */
void
nc_xcor_prepare (NcXcor *xc, NcHICosmo *cosmo)
{
  nc_distance_prepare_if_needed (xc->dist, cosmo);
  ncm_powspec_prepare_if_needed (xc->ps, NCM_MODEL (cosmo));

  xc->RH = nc_hicosmo_RH_Mpc (cosmo);
}

/*
 * QUADPACK's status is a statement about certification, not about the answer:
 * GSL_EROUND means the error estimate stopped improving before reaching the
 * requested tolerance. Where the request carries a safety margin over
 * NcXcor:reltol, the margin is a goal and reltol is the requirement, so the
 * achieved error decides -- aborting on the status alone throws away a result
 * that already meets what the caller asked for.
 */
void
_nc_xcor_check_qag_status (const gchar *where, gint ret, gdouble reltol, gdouble result, gdouble err)
{
  if (ret == GSL_SUCCESS)
    return;

  if ((ret == GSL_EROUND) && (err <= reltol * fabs (result)))
    return;

  g_error ("%s: %s. Result % 22.15g with estimated error %.8e, against the "
           "requested %.8e (NcXcor:reltol %.8e). The quadrature could not "
           "certify that accuracy on this integrand; loosen NcXcor:reltol to "
           "what the integrand carries.",
           where, gsl_strerror (ret), result, err, reltol * fabs (result), reltol);
}

gboolean
_nc_xcor_kernels_limber_disjoint (NcXcorKernel *xclk1, NcXcorKernel *xclk2, gboolean isauto, guint *l_zero)
{
  gdouble zmin, zmax, zmid, zmin_2, zmax_2, zmid_2;
  gint l_limber_1, l_limber_2;

  if (isauto)
    return FALSE;  /* A kernel always overlaps itself. */

  l_limber_1 = nc_xcor_kernel_get_l_limber (xclk1);
  l_limber_2 = nc_xcor_kernel_get_l_limber (xclk2);

  /* A negative l_limber pins the kernel to the non-Limber tier at every
   * multipole, so no threshold exists. */
  if ((l_limber_1 < 0) || (l_limber_2 < 0))
    return FALSE;

  nc_xcor_kernel_get_z_range (xclk1, &zmin, &zmax, &zmid);
  nc_xcor_kernel_get_z_range (xclk2, &zmin_2, &zmax_2, &zmid_2);

  if (GSL_MAX (zmin, zmin_2) < GSL_MIN (zmax, zmax_2))
    return FALSE;

  /* Both kernels are in the Limber tier from the larger of the two thresholds
   * up; below it at least one is still non-Limber and the pair correlates. */
  *l_zero = (guint) GSL_MAX (l_limber_1, l_limber_2);

  return TRUE;
}

/*
 * Dispatches one contiguous multipole range to the configured kernel-space
 * method. Split out so nc_xcor_compute() can drive it for a sub-range without
 * duplicating the switch, which is what a range straddling a Limber threshold
 * needs.
 */
static void
_nc_xcor_kernel_space_compute (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  /* Only some methods know their own error budget -- see
   * nc_xcor_method_has_error_estimate(). The others leave NaN rather than a
   * zero that would read as "no error". */
  if (vp_err != NULL)
    ncm_vector_set_all (vp_err, GSL_NAN);

  if (xc->meth == NC_XCOR_METHOD_KERNEL_GSL)
  {
    _nc_xcor_kernel_gsl (xc, xclk1, xclk2, cosmo, lmin, lmax, isauto, vp);
  }
  else
  {
    const NcXcorKQuad *kquad = _nc_xcor_kquad_for_method (xc->meth);

    g_assert (kquad != NULL);

    _nc_xcor_kernel_space_run (xc, kquad, xclk1, xclk2, cosmo, lmin, lmax, isauto, vp, vp_err);
  }
}

/**
 * nc_xcor_compute_full:
 * @xc: a #NcXcor
 * @xclk1: a #NcXcorKernel
 * @xclk2: (nullable): a #NcXcorKernel, or %NULL for the auto spectrum
 * @cosmo: a #NcHICosmo
 * @lmin: a #guint
 * @lmax: a #guint
 * @vp: a #NcmVector
 * @vp_err: (nullable): a #NcmVector for the error estimate, or %NULL
 *
 * Computes the spectra and, when @vp_err is non-%NULL, an absolute error
 * estimate for each $C_\ell$. The estimate has the same length and units as
 * @vp, so @vp_err[i] / @vp[i] is the relative estimate.
 *
 * Only %NC_XCOR_METHOD_KERNEL_EXACT provides an estimate. Other methods set
 * @vp_err to NaN rather than zero, which would read as "no error".
 *
 * For %NC_XCOR_METHOD_KERNEL_EXACT, the outer integral is exact on both
 * closures -- GL(5) on merged knot panels of two cubics, a bilinear form in
 * the coefficients on cells where either side is a Chebyshev panel -- so the
 * quadrature contributes no error and @vp_err reports the closure fit error,
 * how well nc_xcor_kernel_get_eval_vectorized_full() fitted $W_\ell(k)$,
 * propagated through this pair's conditioning. It propagates
 * $\delta (W^A W^B) = \vert W^A \vert \delta W^B + \vert W^B \vert \delta W^A$:
 *
 * $$ \sigma_\ell \simeq \int \mathrm{d}k\, k^2 \left( \vert W^A_\ell \vert\, \delta W^B_\ell
 *    + \vert W^B_\ell \vert\, \delta W^A_\ell \right) $$
 *
 * The per-cell $\delta W^i$ is the residual the fit achieved, recorded per
 * interval on which the closure is a single polynomial -- a knot interval, or
 * a Chebyshev panel -- while #NcXcorKernel:track-fit-residual is on, which is
 * the default; see nc_xcor_kernel_integrand_peek_residuals(). With tracking
 * off, or on an interval whose refinement was never accepted, $\delta W^i$
 * falls back to the criterion the fit was asked for,
 * $\delta W^i \le \epsilon_i \vert W^i \vert + a_i W^i_\mathrm{max}$, and the
 * same propagation gives
 *
 * $$ \sigma_\ell \simeq (\epsilon_A + \epsilon_B) \int \mathrm{d}k\, k^2 \vert W^A_\ell W^B_\ell \vert
 *    + a_A W^A_{\ell,\mathrm{max}} \int \mathrm{d}k\, k^2 \vert W^B_\ell \vert
 *    + a_B W^B_{\ell,\mathrm{max}} \int \mathrm{d}k\, k^2 \vert W^A_\ell \vert $$
 *
 * Read @vp_err as an upper bound, not as a tight estimate: it is built on
 * absolute values, so it adds error that the integral cancels, and it exceeds
 * the true relative error by 1.2 to 1467 times over the measured cases.
 *
 * Two limits qualify that bound. The fit criterion is an $L^2$ norm over the
 * whole multipole block at each $k$, so a multipole that is sub-dominant
 * within its block is held only to the block's norm and its error can be
 * understated. And the outer integral runs over the intersection of the two
 * closures' fitted $k$ ranges, so @vp_err cannot account for what the
 * intersection discarded. A small @vp_err means the quadrature and the fit are
 * under control, not that the $C_\ell$ is correct.
 *
 * See `dev-notes/xcor_exact_quadrature.md` section 10 for the measured
 * calibration.
 *
 */
void
nc_xcor_compute_full (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp, NcmVector *vp_err)
{
  const guint nell           = ncm_vector_len (vp);
  const gboolean isauto      = (xclk2 == NULL) || (xclk2 == xclk1);
  const gdouble const_factor = 1.0 / gsl_pow_3 (xc->RH);
  guint i;
  gdouble zmin, zmax, zmid;

  if (nell != lmax - lmin + 1)
    g_error ("nc_xcor_compute_full: vector size does not match multipole limits");

  if ((vp_err != NULL) && (ncm_vector_len (vp_err) != nell))
    g_error ("nc_xcor_compute_full: error vector size does not match multipole limits");

  if (isauto)
    xclk2 = xclk1;

  /* Kernel-space (non-Limber) methods do each kernel's radial integral
   * separately and couple the two only through the outer k integral, so two
   * kernels with disjoint redshift support still have a non-zero cross
   * spectrum: two disjoint radial shells are correlated through the same 3D
   * field. The z-overlap short-circuit below is specific to the Limber-z tier,
   * whose C_l is a single integral over the common support, and must not be
   * applied here.
   *
   * The exception is a kernel-space method running kernels that are themselves
   * in the Limber tier: there W(k) is supported only on its own per-ell band,
   * disjoint bins have disjoint support, and the C_l is zero. Integrating it
   * anyway multiplies the two exponential extrapolation tails, which is a
   * smoothing device rather than physics, and gives a large spurious cross
   * spectrum.
   *
   * The tier is chosen per multipole, so a range straddling the threshold is
   * split: the tail from l_zero up is zeroed, the head below it is integrated
   * normally. */
  if (nc_xcor_method_is_kernel_space (xc->meth))
  {
    guint l_zero = 0;

    if (_nc_xcor_kernels_limber_disjoint (xclk1, xclk2, isauto, &l_zero) && (l_zero <= lmax))
    {
      ncm_vector_set_zero (vp);

      /* The zeroed tail is exactly zero, not merely small, so its error is
       * zero too -- the head below overwrites its own entries. */
      if (vp_err != NULL)
        ncm_vector_set_zero (vp_err);

      if (l_zero > lmin)
      {
        NcmVector *vp_head     = ncm_vector_get_subvector (vp, 0, l_zero - lmin);
        NcmVector *vp_err_head = (vp_err != NULL) ? ncm_vector_get_subvector (vp_err, 0, l_zero - lmin) : NULL;

        _nc_xcor_kernel_space_compute (xc, xclk1, xclk2, cosmo, lmin, l_zero - 1, isauto, vp_head, vp_err_head);
        ncm_vector_free (vp_head);
        ncm_vector_clear (&vp_err_head);
      }
    }
    else
    {
      _nc_xcor_kernel_space_compute (xc, xclk1, xclk2, cosmo, lmin, lmax, isauto, vp, vp_err);
    }

    return;
  }

  if (vp_err != NULL)
    ncm_vector_set_all (vp_err, GSL_NAN);

  nc_xcor_kernel_get_z_range (xclk1, &zmin, &zmax, &zmid);

  if (!isauto)
  {
    gdouble zmin_2, zmax_2, zmid_2;

    nc_xcor_kernel_get_z_range (xclk2, &zmin_2, &zmax_2, &zmid_2);
    zmin = GSL_MAX (zmin, zmin_2);
    zmax = GSL_MIN (zmax, zmax_2);
  }

  if (zmin < zmax)
  {
    switch (xc->meth)
    {
      case NC_XCOR_METHOD_LIMBER_Z_GSL:
        _nc_xcor_limber_z_gsl (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      case NC_XCOR_METHOD_LIMBER_Z_CUBATURE:
        _nc_xcor_limber_z_cubature (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }

    if (isauto)
    {
      for (i = 0; i < nell; i++)
      {
        const guint ell                = lmin + i;
        const gdouble const_factor_ell = nc_xcor_kernel_eval_limber_z_prefactor (xclk1, cosmo, ell);

        ncm_vector_mulby (vp, i, const_factor * const_factor_ell * const_factor_ell);
      }
    }
    else
    {
      for (i = 0; i < nell; i++)
      {
        const guint ell                 = lmin + i;
        const gdouble const_factor_ell1 = nc_xcor_kernel_eval_limber_z_prefactor (xclk1, cosmo, ell);
        const gdouble const_factor_ell2 = nc_xcor_kernel_eval_limber_z_prefactor (xclk2, cosmo, ell);

        ncm_vector_mulby (vp, i, const_factor * const_factor_ell1 * const_factor_ell2);
      }
    }
  }
  else
  {
    ncm_vector_set_zero (vp);
  }
}

/**
 * nc_xcor_compute:
 * @xc: a #NcXcor
 * @xclk1: a #NcXcorKernel
 * @xclk2: (nullable): a #NcXcorKernel, or %NULL for the auto spectrum
 * @cosmo: a #NcHICosmo
 * @lmin: a #guint
 * @lmax: a #guint
 * @vp: a #NcmVector
 *
 * Performs the computation of the power spectrum $C_{\ell}^{AB}$. The kernels of
 * observables A and B are @xclk1 and @xclk2. If @xclk2 is NULL, the auto power
 * spectrum is computed. The result for multipoles lmin to lmax (included) is stored in
 * the #NcmVector @vp.
 *
 * Use nc_xcor_compute_full() to also obtain a per-multipole error estimate.
 *
 */
void
nc_xcor_compute (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp)
{
  nc_xcor_compute_full (xc, xclk1, xclk2, cosmo, lmin, lmax, vp, NULL);
}

