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
  xc->closure_type   = NC_XCOR_KERNEL_CLOSURE_SPLINE;
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
   * pair, not the kernel, that constrains it. %NC_XCOR_METHOD_KERNEL_EXACT
   * integrates a pair on the common refinement of the two closures' panels,
   * which requires both to be of the same kind; a per-kernel setting could
   * express a mixed pair the exact method cannot integrate at all. The
   * pointwise methods are indifferent to it, but they read the same property,
   * so no computation mixes the two representations.
   *
   * The spline default is what every method has been calibrated against.
   * %NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV samples the same function over the same
   * domain and differs only in what is fitted to it, so a computation can be
   * switched over and the two compared directly.
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
                                                      NC_XCOR_KERNEL_CLOSURE_SPLINE,
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
 * As nc_xcor_compute(), additionally filling @vp_err with an estimate of the
 * absolute error on each $C_\ell$ -- the same length as @vp, and in the same
 * units, so @vp_err[i] / @vp[i] is the relative estimate.
 *
 * Only %NC_XCOR_METHOD_KERNEL_EXACT provides one; every other method leaves
 * NaN, which is deliberate -- a zero there would read as "no error".
 *
 * The estimate is not a quadrature error, and does not belong to the k
 * integral at all. That method integrates its spline closures exactly, so the
 * outer quadrature contributes nothing; what @vp_err reports is a
 * **kernel-building** error -- how well nc_xcor_kernel_get_eval_vectorized_full()
 * fitted $W_\ell(k)$ -- propagated through the conditioning of this particular
 * pair. The quadrature's only contribution is the amplification factor, which
 * is also the one thing that cannot be known before the pair is formed.
 *
 * Propagating $\delta (W^A W^B) = \vert W^A \vert \delta W^B + \vert W^B
 * \vert \delta W^A$ gives
 *
 * $$ \sigma_\ell \simeq \int \mathrm{d}k\, k^2 \left( \vert W^A_\ell \vert\, \delta W^B_\ell
 *    + \vert W^B_\ell \vert\, \delta W^A_\ell \right) $$
 *
 * with $\delta W^i$ the residual the closure's fit **achieved**, recorded per
 * knot interval while #NcXcorKernel:track-fit-residual is on -- which it is by
 * default. See nc_xcor_kernel_integrand_peek_residuals().
 *
 * With tracking off, or on an interval whose refinement was never accepted,
 * $\delta W^i$ falls back to the criterion the fit was *asked* for,
 * nc_xcor_kernel_integrand_set_tolerances()'s $\delta W^i \le \epsilon_i \vert
 * W^i \vert + a_i W^i_\mathrm{max}$, and the same propagation gives the older
 * three-term form:
 *
 * $$ \sigma_\ell \simeq (\epsilon_A + \epsilon_B) \int \mathrm{d}k\, k^2 \vert W^A_\ell W^B_\ell \vert
 *    + a_A W^A_{\ell,\mathrm{max}} \int \mathrm{d}k\, k^2 \vert W^B_\ell \vert
 *    + a_B W^B_{\ell,\mathrm{max}} \int \mathrm{d}k\, k^2 \vert W^A_\ell \vert $$
 *
 * The first rides on the product, so it is the one the pair's cancellation
 * amplifies. The second and third are each closure's peak-scaled floor weighted
 * by the *other* closure's true size, and they are usually the larger of the
 * two: for cluster top-hat bins at the library defaults they dominate the
 * relative term by one to two orders. That is worth stating plainly, because a
 * floor set per closure is often assumed to reach $C_\ell$ only squared. It
 * does so only where both closures sit on their floors at once; wherever just
 * one does, it is linear in $a$ and weighted by the other's real amplitude.
 *
 * ## How conservative, measured
 *
 * Refinement beats the tolerance it was given by one to three orders, and by
 * an amount that depends on the kernel, so the fallback form above is a ceiling
 * rather than an estimate. Measured against a reference built at reltol
 * $10^{-10}$, over $\ell = 2 \dots 9$ at the library defaults, worst ratio of
 * estimate to true relative error:
 *
 * | pair | true relative error | achieved | tolerance-only |
 * |---|---|---|---|
 * | top-hat, auto | 4.3e-6 to 1.3e-4 | 1.2-50x | 12-858x |
 * | top-hat, cross adjacent | 6.7e-4 to 0.13 | 2.9-16x | 35-320x |
 * | top-hat, cross separated | 0.07 to 8993 | 3.7e-4-11x | 5.1e-3-161x |
 * | Gaussian, auto | 6.5e-8 to 9.6e-7 | 68-630x | 537-7949x |
 * | Gaussian, cross separated | 5.0e-4 to 1.2 | 237-1467x | 6440-50487x |
 *
 * Using the achieved residual is worth a uniform 13 to 34 times across all
 * five, at no measurable cost -- the record is one double per knot per
 * multipole, and building it does not slow the closure down.
 *
 * Read the result as a ceiling still: a small @vp_err is a strong statement, a
 * large one warrants checking rather than despair. Both rows where the ratio
 * drops below one are separated top-hat bins whose $C_\ell$ has no digits left
 * at all, and where the estimate says so -- over the three pairs, three
 * thresholds and both forms, there is no cell where @vp_err calls a $C_\ell$
 * good that is not.
 *
 * Those figures are a **worst case**, and deliberately so: a cluster top-hat has
 * a sharp edge in $\xi$, which gives $W_\ell(k)$ a $1/k$ tail instead of an
 * exponential one. Nothing integrates across that edge -- it is declared through
 * the component's limits, see #NcXcorKernelComponent -- but the tail keeps far
 * more of k-space above the closure's floor, so the fit costs more and is
 * worse. On the same comoving shells a smooth kernel needs 161 knots against
 * the top-hat's 541, and its cross spectrum is accurate to 7.7e-4 rather than
 * 0.13 -- a factor of 165.
 *
 * That comparison also bounds what @vp_err can do, and the achieved residual
 * only half fixes it. The table above still runs from 2.9x over-conservative on
 * a top-hat cross to 1467x on a Gaussian one, because a second mechanism is
 * left: a spline's error alternates sign from panel to panel, and an estimate
 * built on $\vert \cdot \vert$ adds what the integral cancels. That
 * cancellation is near-total for a smoothly fitted kernel and partial for a
 * ragged one, which is exactly the direction of the spread. Closing it needs a
 * signed error model, not a better residual.
 *
 * ## And read them against the tolerances an application actually sets
 *
 * The table uses #NcXcorKernel's bare defaults, which exist to be cheap. A
 * caller that cares sets its own: #NcXcorSSCSij uses reltol $10^{-6}$ with
 * scaled-abstol $10^{-5}$, deliberately offset from each other -- the rationale
 * lives at that object's defaults, and is worth reading before changing either
 * here. At those, on the same top-hat bins, the diagonal is accurate to 5.9e-6
 * rather than 1.3e-4, and the adjacent-bin cross to 2.8e-3 rather than 0.13.
 *
 * The separated-bin cross stays poor in *relative* terms, but that is the wrong
 * measure for it: it is 2.4e-4 of the diagonal, so what it contributes to
 * anything built from these is far below the diagonal's own error. Judge a
 * @vp_err against the amplitude its term carries, never on its own.
 *
 * One thing pushes the other way, against the conservatism: the criterion is an
 * $L^2$ norm over the whole multipole block at each $k$, not over one
 * multipole, so a multipole that is sub-dominant within its block is held only
 * to the block's norm. For those, $\epsilon \vert W_\ell \vert$ understates the
 * fit error.
 *
 * **What it does not cover**, and it is the same classification again -- the
 * other kernel-building error, which is a range rather than a residual. The
 * outer integral runs over the intersection of
 * the two closures' fitted k-ranges, and @vp_err measures only what is inside
 * that range -- it cannot see what the intersection discarded. Two kernels
 * whose closures are fitted on different k-supports lose the non-overlapping
 * part silently, and that loss grows with separation, the same direction in
 * which the cancellation grows. A small @vp_err therefore means the quadrature
 * and the fit are in hand, not that the $C_\ell$ is right.
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

  /* The kernel-space (non-Limber) methods perform each kernel's radial
   * integral separately and couple the two only through the outer k integral,
   * so two kernels with disjoint redshift support still have a non-zero cross
   * spectrum -- two disjoint radial shells are correlated through the same 3D
   * field. The z-overlap short-circuit below is therefore specific to the
   * Limber-z tier, whose C_l is a single integral over the common support, and
   * must not be applied here.
   *
   * The exception is a kernel-space method running kernels that are themselves
   * in the Limber tier: there W(k) is supported only on its own per-ell band,
   * disjoint bins have disjoint support, and the Cl is zero. Integrating it
   * anyway multiplies the two exponential extrapolation tails, which is a
   * numerical smoothing device rather than physics, and yields a large
   * spurious cross spectrum.
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

