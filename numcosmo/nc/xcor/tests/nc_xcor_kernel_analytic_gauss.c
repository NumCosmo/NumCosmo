/***************************************************************************
 *            nc_xcor_kernel_analytic_gauss.c
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
 * NcXcorKernelAnalyticGauss:
 *
 * Gaussian radial window in comoving distance, truncated at a fixed number of
 * standard deviations.
 *
 * \begin{equation}
 *   W(\chi) = \frac{1}{N}
 *   \exp\left[-\frac{(\chi - \chi_0)^2}{2\sigma^2}\right],
 *   \qquad \chi \in [\chi_\mathrm{min}, \chi_\mathrm{max}] ,
 * \end{equation}
 *
 * with $\chi_\mathrm{min} = \max(0, \chi_0 - n\sigma)$,
 * $\chi_\mathrm{max} = \chi_0 + n\sigma$, and $N$ the integral of the
 * exponential over that interval, so that $\int W \mathrm{d}\chi = 1$
 * exactly on the truncated support rather than only in the $n \to \infty$
 * limit.
 *
 * This is the baseline photo-z bin shape: $\sigma \sim 300$ Mpc for a wide
 * bin, $\sigma \sim 50$ Mpc for a thin one. See #NcXcorKernelAnalytic for the
 * unit and normalization conventions.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_gauss.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelAnalyticGauss
{
  /*< private >*/
  NcXcorKernelAnalytic parent_instance;

  gdouble chi_mean;
  gdouble chi_sigma;
  gdouble n_sigma;

  gdouble chi_min;
  gdouble chi_max;
  gdouble norm;
};

enum
{
  PROP_0,
  PROP_CHI_MEAN,
  PROP_CHI_SIGMA,
  PROP_N_SIGMA,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelAnalyticGauss, nc_xcor_kernel_analytic_gauss, NC_TYPE_XCOR_KERNEL_ANALYTIC)

static void
nc_xcor_kernel_analytic_gauss_init (NcXcorKernelAnalyticGauss *xckag)
{
  xckag->chi_mean  = 0.0;
  xckag->chi_sigma = 0.0;
  xckag->n_sigma   = 0.0;
  xckag->chi_min   = 0.0;
  xckag->chi_max   = 0.0;
  xckag->norm      = 0.0;
}

static void
nc_xcor_kernel_analytic_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticGauss *xckag = NC_XCOR_KERNEL_ANALYTIC_GAUSS (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_GAUSS (object));

  switch (prop_id)
  {
    case PROP_CHI_MEAN:
      xckag->chi_mean = g_value_get_double (value);
      break;
    case PROP_CHI_SIGMA:
      xckag->chi_sigma = g_value_get_double (value);
      break;
    case PROP_N_SIGMA:
      xckag->n_sigma = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticGauss *xckag = NC_XCOR_KERNEL_ANALYTIC_GAUSS (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_GAUSS (object));

  switch (prop_id)
  {
    case PROP_CHI_MEAN:
      g_value_set_double (value, xckag->chi_mean);
      break;
    case PROP_CHI_SIGMA:
      g_value_set_double (value, xckag->chi_sigma);
      break;
    case PROP_N_SIGMA:
      g_value_set_double (value, xckag->n_sigma);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_gauss_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_gauss_parent_class)->constructed (object);
  {
    NcXcorKernelAnalyticGauss *xckag = NC_XCOR_KERNEL_ANALYTIC_GAUSS (object);
    const gdouble s2                 = M_SQRT2 * xckag->chi_sigma;

    xckag->chi_min = GSL_MAX (0.0, xckag->chi_mean - xckag->n_sigma * xckag->chi_sigma);
    xckag->chi_max = xckag->chi_mean + xckag->n_sigma * xckag->chi_sigma;

    xckag->norm = 0.5 * sqrt (2.0 * M_PI) * xckag->chi_sigma *
                  (erf ((xckag->chi_max - xckag->chi_mean) / s2) -
                   erf ((xckag->chi_min - xckag->chi_mean) / s2));

    if (!(xckag->norm > 0.0))
      g_error ("nc_xcor_kernel_analytic_gauss_constructed: window has zero mass on "
               "[%.17g, %.17g] Mpc for chi-mean %.17g, chi-sigma %.17g.",
               xckag->chi_min, xckag->chi_max, xckag->chi_mean, xckag->chi_sigma);
  }
}

static void
nc_xcor_kernel_analytic_gauss_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_gauss_parent_class)->finalize (object);
}

static guint _nc_xcor_kernel_analytic_gauss_get_n_comps (NcXcorKernelAnalytic *xcka);
static gdouble _nc_xcor_kernel_analytic_gauss_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi);
static void _nc_xcor_kernel_analytic_gauss_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);

static void
nc_xcor_kernel_analytic_gauss_class_init (NcXcorKernelAnalyticGaussClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);
  NcXcorKernelAnalyticClass *parent_class = NC_XCOR_KERNEL_ANALYTIC_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_analytic_gauss_constructed;
  object_class->finalize    = &nc_xcor_kernel_analytic_gauss_finalize;
  model_class->set_property = &nc_xcor_kernel_analytic_gauss_set_property;
  model_class->get_property = &nc_xcor_kernel_analytic_gauss_get_property;

  ncm_model_class_set_name_nick (model_class, "Analytic Gaussian radial window", "AnalyticGauss");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelAnalyticGauss:chi-mean:
   *
   * Centre $\chi_0$ of the window, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_MEAN,
                                   g_param_spec_double ("chi-mean",
                                                        NULL,
                                                        "Window centre in Mpc",
                                                        0.0, G_MAXDOUBLE, 1000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticGauss:chi-sigma:
   *
   * Standard deviation $\sigma$ of the window, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_SIGMA,
                                   g_param_spec_double ("chi-sigma",
                                                        NULL,
                                                        "Window standard deviation in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 300.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticGauss:n-sigma:
   *
   * Truncation half-width in units of $\sigma$. Part of the window's
   * definition, not a numerical tolerance: the window is exactly zero beyond
   * it and is renormalized over what remains.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_SIGMA,
                                   g_param_spec_double ("n-sigma",
                                                        NULL,
                                                        "Truncation half-width in units of sigma",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->get_n_comps      = &_nc_xcor_kernel_analytic_gauss_get_n_comps;
  parent_class->eval_W_comp      = &_nc_xcor_kernel_analytic_gauss_eval_W_comp;
  parent_class->get_comp_support = &_nc_xcor_kernel_analytic_gauss_get_comp_support;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_kernel_analytic_gauss_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_mean: window centre $\chi_0$, in Mpc
 * @chi_sigma: window standard deviation $\sigma$, in Mpc
 * @n_sigma: truncation half-width, in units of $\sigma$
 *
 * Creates a new #NcXcorKernelAnalyticGauss.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticGauss
 */
NcXcorKernelAnalyticGauss *
nc_xcor_kernel_analytic_gauss_new (NcDistance *dist, NcmPowspec *ps, gdouble chi_mean, gdouble chi_sigma, gdouble n_sigma)
{
  NcXcorKernelAnalyticGauss *xckag = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_GAUSS,
                                                   "dist", dist,
                                                   "powspec", ps,
                                                   "chi-mean", chi_mean,
                                                   "chi-sigma", chi_sigma,
                                                   "n-sigma", n_sigma,
                                                   NULL);

  return xckag;
}

/**
 * nc_xcor_kernel_analytic_gauss_new_full:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_mean: window centre $\chi_0$, in Mpc
 * @chi_sigma: window standard deviation $\sigma$, in Mpc
 * @n_sigma: truncation half-width, in units of $\sigma$
 * @sbi: a #NcmSBesselIntegrator
 *
 * Creates a new #NcXcorKernelAnalyticGauss carrying @sbi, as
 * nc_xcor_kernel_analytic_gauss_new() does not. A #NcXcorKernel only accepts
 * the non-Limber modes of nc_xcor_kernel_set_l_limber() once it holds an
 * integrator, so this is the constructor to use for them.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticGauss
 */
NcXcorKernelAnalyticGauss *
nc_xcor_kernel_analytic_gauss_new_full (NcDistance *dist, NcmPowspec *ps, gdouble chi_mean, gdouble chi_sigma, gdouble n_sigma, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelAnalyticGauss *xckag = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_GAUSS,
                                                   "dist", dist,
                                                   "powspec", ps,
                                                   "chi-mean", chi_mean,
                                                   "chi-sigma", chi_sigma,
                                                   "n-sigma", n_sigma,
                                                   "integrator", sbi,
                                                   NULL);

  return xckag;
}

/**
 * nc_xcor_kernel_analytic_gauss_get_chi_mean:
 * @xckag: a #NcXcorKernelAnalyticGauss
 *
 * Returns: the window centre $\chi_0$, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_gauss_get_chi_mean (NcXcorKernelAnalyticGauss *xckag)
{
  return xckag->chi_mean;
}

/**
 * nc_xcor_kernel_analytic_gauss_get_chi_sigma:
 * @xckag: a #NcXcorKernelAnalyticGauss
 *
 * Returns: the window standard deviation $\sigma$, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_gauss_get_chi_sigma (NcXcorKernelAnalyticGauss *xckag)
{
  return xckag->chi_sigma;
}

/**
 * nc_xcor_kernel_analytic_gauss_get_n_sigma:
 * @xckag: a #NcXcorKernelAnalyticGauss
 *
 * Returns: the truncation half-width, in units of $\sigma$.
 */
gdouble
nc_xcor_kernel_analytic_gauss_get_n_sigma (NcXcorKernelAnalyticGauss *xckag)
{
  return xckag->n_sigma;
}

static guint
_nc_xcor_kernel_analytic_gauss_get_n_comps (NcXcorKernelAnalytic *xcka)
{
  return 1;
}

static gdouble
_nc_xcor_kernel_analytic_gauss_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi)
{
  NcXcorKernelAnalyticGauss *xckag = NC_XCOR_KERNEL_ANALYTIC_GAUSS (xcka);
  const gdouble d                  = chi - xckag->chi_mean;

  if ((chi < xckag->chi_min) || (chi > xckag->chi_max))
    return 0.0;

  return exp (-0.5 * gsl_pow_2 (d / xckag->chi_sigma)) / xckag->norm;
}

static void
_nc_xcor_kernel_analytic_gauss_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  NcXcorKernelAnalyticGauss *xckag = NC_XCOR_KERNEL_ANALYTIC_GAUSS (xcka);

  *chi_min = xckag->chi_min;
  *chi_max = xckag->chi_max;
}

