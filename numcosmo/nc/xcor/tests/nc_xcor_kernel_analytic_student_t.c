/***************************************************************************
 *            nc_xcor_kernel_analytic_student_t.c
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
 * NcXcorKernelAnalyticStudentT:
 *
 * Radial window with power-law tails: a Student-t profile in comoving
 * distance.
 *
 * \begin{equation}
 *   W(\chi) = \frac{1}{D}
 *   \left[1 + \frac{1}{\nu}\left(\frac{\chi - \chi_0}{s}\right)^2\right]
 *   ^{-\frac{\nu + 1}{2}} ,
 * \end{equation}
 *
 * with $D$ fixed in closed form so that $\int W \mathrm{d}\chi = 1$ over the
 * truncated support.
 *
 * The tail falls as $|\chi - \chi_0|^{-(\nu+1)}$, not exponentially. That is
 * the point of the shape: a Gaussian window's tail is negligible a few widths
 * out, so where it is truncated barely matters and its transform decays fast;
 * a power-law tail does neither, which is how several physical kernels
 * actually behave. $\nu$ tunes the decay continuously -- $\nu = 1$ is Cauchy,
 * the heaviest useful case, and large $\nu$ returns to the Gaussian of
 * #NcXcorKernelAnalyticGauss.
 *
 * Truncation therefore discards real mass here, and deliberately so: at
 * $\nu = 1$, cutting at $10s$ leaves out about 6% of the profile. The window
 * is renormalized over exactly what is kept, so the definition stays exact and
 * both sides of a cross-code comparison compute the same thing.
 * nc_xcor_kernel_analytic_student_t_get_tail_mass() reports the fraction
 * discarded, so a spec can state it rather than discover it.
 *
 * See #NcXcorKernelAnalytic for the unit and normalization conventions.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_student_t.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelAnalyticStudentT
{
  /*< private >*/
  NcXcorKernelAnalytic parent_instance;

  gdouble chi_mean;
  gdouble chi_scale;
  gdouble nu;
  gdouble n_scale;

  gdouble chi_min;
  gdouble chi_max;
  gdouble norm;
  gdouble kept_mass;
};

enum
{
  PROP_0,
  PROP_CHI_MEAN,
  PROP_CHI_SCALE,
  PROP_NU,
  PROP_N_SCALE,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelAnalyticStudentT, nc_xcor_kernel_analytic_student_t, NC_TYPE_XCOR_KERNEL_ANALYTIC)

static void
nc_xcor_kernel_analytic_student_t_init (NcXcorKernelAnalyticStudentT *xckas)
{
  xckas->chi_mean  = 0.0;
  xckas->chi_scale = 0.0;
  xckas->nu        = 0.0;
  xckas->n_scale   = 0.0;
  xckas->chi_min   = 0.0;
  xckas->chi_max   = 0.0;
  xckas->norm      = 0.0;
  xckas->kept_mass = 0.0;
}

static void
nc_xcor_kernel_analytic_student_t_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticStudentT *xckas = NC_XCOR_KERNEL_ANALYTIC_STUDENT_T (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_STUDENT_T (object));

  switch (prop_id)
  {
    case PROP_CHI_MEAN:
      xckas->chi_mean = g_value_get_double (value);
      break;
    case PROP_CHI_SCALE:
      xckas->chi_scale = g_value_get_double (value);
      break;
    case PROP_NU:
      xckas->nu = g_value_get_double (value);
      break;
    case PROP_N_SCALE:
      xckas->n_scale = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_student_t_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticStudentT *xckas = NC_XCOR_KERNEL_ANALYTIC_STUDENT_T (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_STUDENT_T (object));

  switch (prop_id)
  {
    case PROP_CHI_MEAN:
      g_value_set_double (value, xckas->chi_mean);
      break;
    case PROP_CHI_SCALE:
      g_value_set_double (value, xckas->chi_scale);
      break;
    case PROP_NU:
      g_value_set_double (value, xckas->nu);
      break;
    case PROP_N_SCALE:
      g_value_set_double (value, xckas->n_scale);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_student_t_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_student_t_parent_class)->constructed (object);
  {
    NcXcorKernelAnalyticStudentT *xckas = NC_XCOR_KERNEL_ANALYTIC_STUDENT_T (object);
    const gdouble nu                    = xckas->nu;

    xckas->chi_min = GSL_MAX (0.0, xckas->chi_mean - xckas->n_scale * xckas->chi_scale);
    xckas->chi_max = xckas->chi_mean + xckas->n_scale * xckas->chi_scale;

    {
      /* Unit-scale variable: the profile is the standard t density up to the
       * constant below, so the kept mass is a difference of its CDF. */
      const gdouble t_lo = (xckas->chi_min - xckas->chi_mean) / xckas->chi_scale;
      const gdouble t_hi = (xckas->chi_max - xckas->chi_mean) / xckas->chi_scale;

      /* C_nu = sqrt(nu pi) Gamma(nu/2) / Gamma((nu+1)/2) turns the CDF
       * difference into the integral of the unnormalized profile. */
      const gdouble ln_C = 0.5 * log (nu * M_PI) + gsl_sf_lngamma (0.5 * nu) - gsl_sf_lngamma (0.5 * (nu + 1.0));

      xckas->kept_mass = gsl_cdf_tdist_P (t_hi, nu) - gsl_cdf_tdist_P (t_lo, nu);
      xckas->norm      = xckas->chi_scale * exp (ln_C) * xckas->kept_mass;
    }

    if (!(xckas->norm > 0.0))
      g_error ("nc_xcor_kernel_analytic_student_t_constructed: window has zero mass on "
               "[%.17g, %.17g] Mpc for chi-mean %.17g, chi-scale %.17g, nu %.17g.",
               xckas->chi_min, xckas->chi_max, xckas->chi_mean, xckas->chi_scale, nu);
  }
}

static void
nc_xcor_kernel_analytic_student_t_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_student_t_parent_class)->finalize (object);
}

static guint _nc_xcor_kernel_analytic_student_t_get_n_comps (NcXcorKernelAnalytic *xcka);
static gdouble _nc_xcor_kernel_analytic_student_t_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi);
static void _nc_xcor_kernel_analytic_student_t_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);

static void
nc_xcor_kernel_analytic_student_t_class_init (NcXcorKernelAnalyticStudentTClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);
  NcXcorKernelAnalyticClass *parent_class = NC_XCOR_KERNEL_ANALYTIC_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_analytic_student_t_constructed;
  object_class->finalize    = &nc_xcor_kernel_analytic_student_t_finalize;
  model_class->set_property = &nc_xcor_kernel_analytic_student_t_set_property;
  model_class->get_property = &nc_xcor_kernel_analytic_student_t_get_property;

  ncm_model_class_set_name_nick (model_class, "Analytic power-law-tailed radial window", "AnalyticStudentT");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelAnalyticStudentT:chi-mean:
   *
   * Centre $\chi_0$ of the window, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_MEAN,
                                   g_param_spec_double ("chi-mean",
                                                        NULL,
                                                        "Window centre in Mpc",
                                                        0.0, G_MAXDOUBLE, 1500.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticStudentT:chi-scale:
   *
   * Scale $s$ of the window, in Mpc. For $\nu > 2$ the profile has variance
   * $s^2 \nu / (\nu - 2)$; for $\nu \le 2$ it has none, which is the regime the
   * shape exists to probe.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_SCALE,
                                   g_param_spec_double ("chi-scale",
                                                        NULL,
                                                        "Window scale in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 300.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticStudentT:nu:
   *
   * Degrees of freedom $\nu$, which set the tail exponent: $W$ falls as
   * $|\chi - \chi_0|^{-(\nu+1)}$. $\nu = 1$ is Cauchy; large $\nu$ approaches a
   * Gaussian.
   */
  g_object_class_install_property (object_class,
                                   PROP_NU,
                                   g_param_spec_double ("nu",
                                                        NULL,
                                                        "Degrees of freedom",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticStudentT:n-scale:
   *
   * Truncation half-width in units of $s$. With power-law tails this discards
   * a fraction of the profile that is not negligible, so it is part of the
   * window's definition rather than a tolerance; see
   * nc_xcor_kernel_analytic_student_t_get_tail_mass().
   */
  g_object_class_install_property (object_class,
                                   PROP_N_SCALE,
                                   g_param_spec_double ("n-scale",
                                                        NULL,
                                                        "Truncation half-width in units of the scale",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->get_n_comps      = &_nc_xcor_kernel_analytic_student_t_get_n_comps;
  parent_class->eval_W_comp      = &_nc_xcor_kernel_analytic_student_t_eval_W_comp;
  parent_class->get_comp_support = &_nc_xcor_kernel_analytic_student_t_get_comp_support;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_kernel_analytic_student_t_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_mean: window centre $\chi_0$, in Mpc
 * @chi_scale: window scale $s$, in Mpc
 * @nu: degrees of freedom $\nu$; 1 is Cauchy
 * @n_scale: truncation half-width, in units of $s$
 *
 * Creates a new #NcXcorKernelAnalyticStudentT.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticStudentT
 */
NcXcorKernelAnalyticStudentT *
nc_xcor_kernel_analytic_student_t_new (NcDistance *dist, NcmPowspec *ps, gdouble chi_mean, gdouble chi_scale, gdouble nu, gdouble n_scale)
{
  NcXcorKernelAnalyticStudentT *xckas = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_STUDENT_T,
                                                      "dist", dist,
                                                      "powspec", ps,
                                                      "chi-mean", chi_mean,
                                                      "chi-scale", chi_scale,
                                                      "nu", nu,
                                                      "n-scale", n_scale,
                                                      NULL);

  return xckas;
}

/**
 * nc_xcor_kernel_analytic_student_t_new_full:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_mean: window centre $\chi_0$, in Mpc
 * @chi_scale: window scale $s$, in Mpc
 * @nu: degrees of freedom $\nu$; 1 is Cauchy
 * @n_scale: truncation half-width, in units of $s$
 * @sbi: a #NcmSBesselIntegrator
 *
 * Creates a new #NcXcorKernelAnalyticStudentT carrying @sbi, as
 * nc_xcor_kernel_analytic_student_t_new() does not. A #NcXcorKernel only
 * accepts the non-Limber modes of nc_xcor_kernel_set_l_limber() once it holds
 * an integrator, so this is the constructor to use for them.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticStudentT
 */
NcXcorKernelAnalyticStudentT *
nc_xcor_kernel_analytic_student_t_new_full (NcDistance *dist, NcmPowspec *ps, gdouble chi_mean, gdouble chi_scale, gdouble nu, gdouble n_scale, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelAnalyticStudentT *xckas = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_STUDENT_T,
                                                      "dist", dist,
                                                      "powspec", ps,
                                                      "chi-mean", chi_mean,
                                                      "chi-scale", chi_scale,
                                                      "nu", nu,
                                                      "n-scale", n_scale,
                                                      "integrator", sbi,
                                                      NULL);

  return xckas;
}

/**
 * nc_xcor_kernel_analytic_student_t_get_chi_mean:
 * @xckas: a #NcXcorKernelAnalyticStudentT
 *
 * Returns: the window centre $\chi_0$, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_student_t_get_chi_mean (NcXcorKernelAnalyticStudentT *xckas)
{
  return xckas->chi_mean;
}

/**
 * nc_xcor_kernel_analytic_student_t_get_chi_scale:
 * @xckas: a #NcXcorKernelAnalyticStudentT
 *
 * Returns: the window scale $s$, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_student_t_get_chi_scale (NcXcorKernelAnalyticStudentT *xckas)
{
  return xckas->chi_scale;
}

/**
 * nc_xcor_kernel_analytic_student_t_get_nu:
 * @xckas: a #NcXcorKernelAnalyticStudentT
 *
 * Returns: the degrees of freedom $\nu$.
 */
gdouble
nc_xcor_kernel_analytic_student_t_get_nu (NcXcorKernelAnalyticStudentT *xckas)
{
  return xckas->nu;
}

/**
 * nc_xcor_kernel_analytic_student_t_get_n_scale:
 * @xckas: a #NcXcorKernelAnalyticStudentT
 *
 * Returns: the truncation half-width, in units of $s$.
 */
gdouble
nc_xcor_kernel_analytic_student_t_get_n_scale (NcXcorKernelAnalyticStudentT *xckas)
{
  return xckas->n_scale;
}

/**
 * nc_xcor_kernel_analytic_student_t_get_tail_mass:
 * @xckas: a #NcXcorKernelAnalyticStudentT
 *
 * Gets the fraction of the untruncated profile that falls outside the support
 * and is therefore not represented. Unlike a Gaussian's, this is not
 * negligible: at $\nu = 1$ and $n = 10$ it is about 6%. It is reported so that
 * a benchmark specification can state it.
 *
 * Returns: the discarded fraction, between 0 and 1.
 */
gdouble
nc_xcor_kernel_analytic_student_t_get_tail_mass (NcXcorKernelAnalyticStudentT *xckas)
{
  return 1.0 - xckas->kept_mass;
}

static guint
_nc_xcor_kernel_analytic_student_t_get_n_comps (NcXcorKernelAnalytic *xcka)
{
  return 1;
}

static gdouble
_nc_xcor_kernel_analytic_student_t_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi)
{
  NcXcorKernelAnalyticStudentT *xckas = NC_XCOR_KERNEL_ANALYTIC_STUDENT_T (xcka);
  const gdouble t                     = (chi - xckas->chi_mean) / xckas->chi_scale;

  if ((chi < xckas->chi_min) || (chi > xckas->chi_max))
    return 0.0;

  return pow (1.0 + t * t / xckas->nu, -0.5 * (xckas->nu + 1.0)) / xckas->norm;
}

static void
_nc_xcor_kernel_analytic_student_t_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  NcXcorKernelAnalyticStudentT *xckas = NC_XCOR_KERNEL_ANALYTIC_STUDENT_T (xcka);

  *chi_min = xckas->chi_min;
  *chi_max = xckas->chi_max;
}

