/***************************************************************************
 *            nc_xcor_kernel_analytic_tophat_smooth.c
 *
 *  Thu August 21 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcXcorKernelAnalyticTophatSmooth:
 *
 * Top-hat convolved with a Gaussian: what a real tomographic bin looks like.
 *
 * \begin{equation}
 *   W(\chi) = \frac{1}{N}\left[
 *   \mathrm{erf}\left(\frac{\chi_\mathrm{u} - \chi}{\sqrt{2}\sigma}\right) -
 *   \mathrm{erf}\left(\frac{\chi_\mathrm{l} - \chi}{\sqrt{2}\sigma}\right)
 *   \right] ,
 * \end{equation}
 *
 * the exact convolution of the indicator of $[\chi_\mathrm{l}, \chi_\mathrm{u}]$
 * with a Gaussian of width $\sigma$, truncated at $n\sigma$ beyond each edge and
 * renormalized there in closed form.
 *
 * This sits between #NcXcorKernelAnalyticTophat and #NcXcorKernelAnalyticGauss
 * and is the more honest model of a photo-z bin than either: bin edges are
 * sharp in true redshift but the photo-z scatter smooths them, so the window is
 * flat in the middle with rounded shoulders rather than either a hard step or a
 * single bump. $\sigma \to 0$ recovers the top-hat and
 * $\chi_\mathrm{u} \to \chi_\mathrm{l}$ the Gaussian, so it interpolates the two
 * stress cases the other shapes probe separately.
 *
 * See #NcXcorKernelAnalytic for the unit and normalization conventions.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_tophat_smooth.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelAnalyticTophatSmooth
{
  /*< private >*/
  NcXcorKernelAnalytic parent_instance;

  gdouble chi_lower;
  gdouble chi_upper;
  gdouble chi_sigma;
  gdouble n_sigma;

  gdouble chi_min;
  gdouble chi_max;
  gdouble norm;
};

enum
{
  PROP_0,
  PROP_CHI_LOWER,
  PROP_CHI_UPPER,
  PROP_CHI_SIGMA,
  PROP_N_SIGMA,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelAnalyticTophatSmooth, nc_xcor_kernel_analytic_tophat_smooth, NC_TYPE_XCOR_KERNEL_ANALYTIC)

static void
nc_xcor_kernel_analytic_tophat_smooth_init (NcXcorKernelAnalyticTophatSmooth *xckats)
{
  xckats->chi_lower = 0.0;
  xckats->chi_upper = 0.0;
  xckats->chi_sigma = 0.0;
  xckats->n_sigma   = 0.0;
  xckats->chi_min   = 0.0;
  xckats->chi_max   = 0.0;
  xckats->norm      = 0.0;
}

static void
nc_xcor_kernel_analytic_tophat_smooth_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticTophatSmooth *xckats = NC_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH (object));

  switch (prop_id)
  {
    case PROP_CHI_LOWER:
      xckats->chi_lower = g_value_get_double (value);
      break;
    case PROP_CHI_UPPER:
      xckats->chi_upper = g_value_get_double (value);
      break;
    case PROP_CHI_SIGMA:
      xckats->chi_sigma = g_value_get_double (value);
      break;
    case PROP_N_SIGMA:
      xckats->n_sigma = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_tophat_smooth_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticTophatSmooth *xckats = NC_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH (object));

  switch (prop_id)
  {
    case PROP_CHI_LOWER:
      g_value_set_double (value, xckats->chi_lower);
      break;
    case PROP_CHI_UPPER:
      g_value_set_double (value, xckats->chi_upper);
      break;
    case PROP_CHI_SIGMA:
      g_value_set_double (value, xckats->chi_sigma);
      break;
    case PROP_N_SIGMA:
      g_value_set_double (value, xckats->n_sigma);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* Antiderivative of erf: d/dt [t erf(t) + exp(-t^2)/sqrt(pi)] = erf(t). */
static gdouble
_erf_antiderivative (gdouble t)
{
  return t * erf (t) + exp (-t * t) / sqrt (M_PI);
}

static void
nc_xcor_kernel_analytic_tophat_smooth_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_tophat_smooth_parent_class)->constructed (object);
  {
    NcXcorKernelAnalyticTophatSmooth *xckats = NC_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH (object);

    if (!(xckats->chi_upper > xckats->chi_lower))
      g_error ("nc_xcor_kernel_analytic_tophat_smooth_constructed: empty bin, "
               "chi-lower %.17g >= chi-upper %.17g.", xckats->chi_lower, xckats->chi_upper);

    xckats->chi_min = GSL_MAX (0.0, xckats->chi_lower - xckats->n_sigma * xckats->chi_sigma);
    xckats->chi_max = xckats->chi_upper + xckats->n_sigma * xckats->chi_sigma;

    {
      /* int erf((c - chi)/(sqrt(2) sigma)) dchi over [A, B] is
       * sqrt(2) sigma [F(t(A)) - F(t(B))] with t(chi) = (c - chi)/(sqrt(2) sigma)
       * and F the antiderivative above -- so the whole normalization is closed
       * form, and exact on the truncated support rather than only as n -> inf. */
      const gdouble s2   = M_SQRT2 * xckats->chi_sigma;
      const gdouble tu_a = (xckats->chi_upper - xckats->chi_min) / s2;
      const gdouble tu_b = (xckats->chi_upper - xckats->chi_max) / s2;
      const gdouble tl_a = (xckats->chi_lower - xckats->chi_min) / s2;
      const gdouble tl_b = (xckats->chi_lower - xckats->chi_max) / s2;

      xckats->norm = s2 * ((_erf_antiderivative (tu_a) - _erf_antiderivative (tu_b)) -
                           (_erf_antiderivative (tl_a) - _erf_antiderivative (tl_b)));
    }

    if (!(xckats->norm > 0.0))
      g_error ("nc_xcor_kernel_analytic_tophat_smooth_constructed: window has zero mass on "
               "[%.17g, %.17g] Mpc.", xckats->chi_min, xckats->chi_max);
  }
}

static void
nc_xcor_kernel_analytic_tophat_smooth_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_tophat_smooth_parent_class)->finalize (object);
}

static guint _nc_xcor_kernel_analytic_tophat_smooth_get_n_comps (NcXcorKernelAnalytic *xcka);
static gdouble _nc_xcor_kernel_analytic_tophat_smooth_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi);
static void _nc_xcor_kernel_analytic_tophat_smooth_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);

static void
nc_xcor_kernel_analytic_tophat_smooth_class_init (NcXcorKernelAnalyticTophatSmoothClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);
  NcXcorKernelAnalyticClass *parent_class = NC_XCOR_KERNEL_ANALYTIC_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_analytic_tophat_smooth_constructed;
  object_class->finalize    = &nc_xcor_kernel_analytic_tophat_smooth_finalize;
  model_class->set_property = &nc_xcor_kernel_analytic_tophat_smooth_set_property;
  model_class->get_property = &nc_xcor_kernel_analytic_tophat_smooth_get_property;

  ncm_model_class_set_name_nick (model_class, "Analytic smoothed top-hat radial window", "AnalyticTophatSmooth");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelAnalyticTophatSmooth:chi-lower:
   *
   * Lower edge $\chi_\mathrm{l}$ of the bin before smoothing, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_LOWER,
                                   g_param_spec_double ("chi-lower",
                                                        NULL,
                                                        "Lower bin edge in Mpc",
                                                        0.0, G_MAXDOUBLE, 1000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticTophatSmooth:chi-upper:
   *
   * Upper edge $\chi_\mathrm{u}$ of the bin before smoothing, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_UPPER,
                                   g_param_spec_double ("chi-upper",
                                                        NULL,
                                                        "Upper bin edge in Mpc",
                                                        0.0, G_MAXDOUBLE, 2000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticTophatSmooth:chi-sigma:
   *
   * Width $\sigma$ of the smoothing, in Mpc: the photo-z scatter that rounds
   * the bin edges.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_SIGMA,
                                   g_param_spec_double ("chi-sigma",
                                                        NULL,
                                                        "Smoothing width in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 150.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticTophatSmooth:n-sigma:
   *
   * How far beyond each edge the window is kept, in units of $\sigma$. Part of
   * the definition, as in #NcXcorKernelAnalyticGauss.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_SIGMA,
                                   g_param_spec_double ("n-sigma",
                                                        NULL,
                                                        "Truncation beyond each edge in units of sigma",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 6.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->get_n_comps      = &_nc_xcor_kernel_analytic_tophat_smooth_get_n_comps;
  parent_class->eval_W_comp      = &_nc_xcor_kernel_analytic_tophat_smooth_eval_W_comp;
  parent_class->get_comp_support = &_nc_xcor_kernel_analytic_tophat_smooth_get_comp_support;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_kernel_analytic_tophat_smooth_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_lower: lower bin edge, in Mpc
 * @chi_upper: upper bin edge, in Mpc
 * @chi_sigma: smoothing width $\sigma$, in Mpc
 * @n_sigma: truncation beyond each edge, in units of $\sigma$
 *
 * Creates a new #NcXcorKernelAnalyticTophatSmooth.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticTophatSmooth
 */
NcXcorKernelAnalyticTophatSmooth *
nc_xcor_kernel_analytic_tophat_smooth_new (NcDistance *dist, NcmPowspec *ps, gdouble chi_lower, gdouble chi_upper, gdouble chi_sigma, gdouble n_sigma)
{
  NcXcorKernelAnalyticTophatSmooth *xckats = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH,
                                                           "dist", dist,
                                                           "powspec", ps,
                                                           "chi-lower", chi_lower,
                                                           "chi-upper", chi_upper,
                                                           "chi-sigma", chi_sigma,
                                                           "n-sigma", n_sigma,
                                                           NULL);

  return xckats;
}

/**
 * nc_xcor_kernel_analytic_tophat_smooth_new_full:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_lower: lower bin edge, in Mpc
 * @chi_upper: upper bin edge, in Mpc
 * @chi_sigma: smoothing width $\sigma$, in Mpc
 * @n_sigma: truncation beyond each edge, in units of $\sigma$
 * @sbi: a #NcmSBesselIntegrator
 *
 * Creates a new #NcXcorKernelAnalyticTophatSmooth carrying @sbi, as
 * nc_xcor_kernel_analytic_tophat_smooth_new() does not. A #NcXcorKernel only
 * accepts the non-Limber modes of nc_xcor_kernel_set_l_limber() once it holds
 * an integrator, so this is the constructor to use for them.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticTophatSmooth
 */
NcXcorKernelAnalyticTophatSmooth *
nc_xcor_kernel_analytic_tophat_smooth_new_full (NcDistance *dist, NcmPowspec *ps, gdouble chi_lower, gdouble chi_upper, gdouble chi_sigma, gdouble n_sigma, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelAnalyticTophatSmooth *xckats = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH,
                                                           "dist", dist,
                                                           "powspec", ps,
                                                           "chi-lower", chi_lower,
                                                           "chi-upper", chi_upper,
                                                           "chi-sigma", chi_sigma,
                                                           "n-sigma", n_sigma,
                                                           "integrator", sbi,
                                                           NULL);

  return xckats;
}

/**
 * nc_xcor_kernel_analytic_tophat_smooth_get_chi_lower:
 * @xckats: a #NcXcorKernelAnalyticTophatSmooth
 *
 * Returns: the lower bin edge, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_tophat_smooth_get_chi_lower (NcXcorKernelAnalyticTophatSmooth *xckats)
{
  return xckats->chi_lower;
}

/**
 * nc_xcor_kernel_analytic_tophat_smooth_get_chi_upper:
 * @xckats: a #NcXcorKernelAnalyticTophatSmooth
 *
 * Returns: the upper bin edge, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_tophat_smooth_get_chi_upper (NcXcorKernelAnalyticTophatSmooth *xckats)
{
  return xckats->chi_upper;
}

/**
 * nc_xcor_kernel_analytic_tophat_smooth_get_chi_sigma:
 * @xckats: a #NcXcorKernelAnalyticTophatSmooth
 *
 * Returns: the smoothing width $\sigma$, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_tophat_smooth_get_chi_sigma (NcXcorKernelAnalyticTophatSmooth *xckats)
{
  return xckats->chi_sigma;
}

static guint
_nc_xcor_kernel_analytic_tophat_smooth_get_n_comps (NcXcorKernelAnalytic *xcka)
{
  return 1;
}

static gdouble
_nc_xcor_kernel_analytic_tophat_smooth_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi)
{
  NcXcorKernelAnalyticTophatSmooth *xckats = NC_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH (xcka);
  const gdouble s2                         = M_SQRT2 * xckats->chi_sigma;
  const gdouble tu                         = (xckats->chi_upper - chi) / s2;
  const gdouble tl                         = (xckats->chi_lower - chi) / s2;

  if ((chi < xckats->chi_min) || (chi > xckats->chi_max))
    return 0.0;

  /* erf (tu) - erf (tl) with tu > tl. Both erf values approach the same
   * limit in either tail, so the difference cancels to the window's own size
   * -- 1e-8 relative at 6 sigma, which no relative tolerance can see past.
   * erfc keeps the small quantity small throughout. */
  if ((tu >= 0.0) && (tl >= 0.0))
    return (erfc (tl) - erfc (tu)) / xckats->norm;
  else if ((tu <= 0.0) && (tl <= 0.0))
    return (erfc (-tu) - erfc (-tl)) / xckats->norm;
  else
    return (erf (tu) - erf (tl)) / xckats->norm;
}

static void
_nc_xcor_kernel_analytic_tophat_smooth_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  NcXcorKernelAnalyticTophatSmooth *xckats = NC_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH (xcka);

  *chi_min = xckats->chi_min;
  *chi_max = xckats->chi_max;
}

