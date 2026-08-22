/***************************************************************************
 *            nc_xcor_kernel_analytic_power_exp.c
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
 * NcXcorKernelAnalyticPowerExp:
 *
 * Power-law times stretched exponential, the shape of a realistic redshift
 * distribution.
 *
 * \begin{equation}
 *   W(\chi) = \frac{1}{N}
 *   \left(\frac{\chi}{\chi_0}\right)^{\alpha}
 *   \exp\left[-\left(\frac{\chi}{\chi_0}\right)^{\beta}\right],
 *   \qquad \chi \in [\chi_\mathrm{l}, \chi_\mathrm{u}] ,
 * \end{equation}
 *
 * with $N$ in closed form, an incomplete gamma function, so that
 * $\int W \mathrm{d}\chi = 1$ over exactly the interval kept.
 *
 * This is the family the survey $\mathrm{d}n/\mathrm{d}z$ parametrizations
 * belong to, written in comoving distance rather than redshift: rising as a
 * power law, falling as a stretched exponential, skewed rather than symmetric.
 * $\alpha \sim 2$, $\beta \sim 1.5$ gives an LSST-like source distribution;
 * small $\alpha$ with $\beta \sim 1$ gives the broad, low-$\chi$-weighted shape
 * of an ISW kernel. Both are emulations with the right structure, exactly
 * evaluable, not fits to a particular survey.
 *
 * See #NcXcorKernelAnalytic for the unit and normalization conventions.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/nc_xcor_kernel_analytic_power_exp.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelAnalyticPowerExp
{
  /*< private >*/
  NcXcorKernelAnalytic parent_instance;

  gdouble chi_scale;
  gdouble alpha;
  gdouble beta;
  gdouble chi_lower;
  gdouble chi_upper;

  gdouble norm;
};

enum
{
  PROP_0,
  PROP_CHI_SCALE,
  PROP_ALPHA,
  PROP_BETA,
  PROP_CHI_LOWER,
  PROP_CHI_UPPER,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelAnalyticPowerExp, nc_xcor_kernel_analytic_power_exp, NC_TYPE_XCOR_KERNEL_ANALYTIC)

static void
nc_xcor_kernel_analytic_power_exp_init (NcXcorKernelAnalyticPowerExp *xckap)
{
  xckap->chi_scale = 0.0;
  xckap->alpha     = 0.0;
  xckap->beta      = 0.0;
  xckap->chi_lower = 0.0;
  xckap->chi_upper = 0.0;
  xckap->norm      = 0.0;
}

static void
nc_xcor_kernel_analytic_power_exp_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticPowerExp *xckap = NC_XCOR_KERNEL_ANALYTIC_POWER_EXP (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_POWER_EXP (object));

  switch (prop_id)
  {
    case PROP_CHI_SCALE:
      xckap->chi_scale = g_value_get_double (value);
      break;
    case PROP_ALPHA:
      xckap->alpha = g_value_get_double (value);
      break;
    case PROP_BETA:
      xckap->beta = g_value_get_double (value);
      break;
    case PROP_CHI_LOWER:
      xckap->chi_lower = g_value_get_double (value);
      break;
    case PROP_CHI_UPPER:
      xckap->chi_upper = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_power_exp_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticPowerExp *xckap = NC_XCOR_KERNEL_ANALYTIC_POWER_EXP (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_POWER_EXP (object));

  switch (prop_id)
  {
    case PROP_CHI_SCALE:
      g_value_set_double (value, xckap->chi_scale);
      break;
    case PROP_ALPHA:
      g_value_set_double (value, xckap->alpha);
      break;
    case PROP_BETA:
      g_value_set_double (value, xckap->beta);
      break;
    case PROP_CHI_LOWER:
      g_value_set_double (value, xckap->chi_lower);
      break;
    case PROP_CHI_UPPER:
      g_value_set_double (value, xckap->chi_upper);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_power_exp_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_power_exp_parent_class)->constructed (object);
  {
    NcXcorKernelAnalyticPowerExp *xckap = NC_XCOR_KERNEL_ANALYTIC_POWER_EXP (object);

    if (!(xckap->chi_upper > xckap->chi_lower))
      g_error ("nc_xcor_kernel_analytic_power_exp_constructed: empty support, "
               "chi-lower %.17g >= chi-upper %.17g.", xckap->chi_lower, xckap->chi_upper);

    {
      /* With u = (chi/chi_0)^beta the integral is an incomplete gamma:
       *   int chi^alpha exp(-(chi/chi_0)^beta) dchi
       *     = (chi_0/beta) Gamma(s) [P(s, u_hi) - P(s, u_lo)],  s = (alpha+1)/beta,
       * which is exact and needs no quadrature. */
      const gdouble s     = (xckap->alpha + 1.0) / xckap->beta;
      const gdouble u_lo  = pow (xckap->chi_lower / xckap->chi_scale, xckap->beta);
      const gdouble u_hi  = pow (xckap->chi_upper / xckap->chi_scale, xckap->beta);
      const gdouble gamma = exp (gsl_sf_lngamma (s));

      xckap->norm = (xckap->chi_scale / xckap->beta) * gamma *
                    (gsl_sf_gamma_inc_P (s, u_hi) - gsl_sf_gamma_inc_P (s, u_lo));
    }

    if (!(xckap->norm > 0.0))
      g_error ("nc_xcor_kernel_analytic_power_exp_constructed: window has zero mass on "
               "[%.17g, %.17g] Mpc for chi-scale %.17g, alpha %.17g, beta %.17g.",
               xckap->chi_lower, xckap->chi_upper, xckap->chi_scale, xckap->alpha, xckap->beta);
  }
}

static void
nc_xcor_kernel_analytic_power_exp_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_power_exp_parent_class)->finalize (object);
}

static guint _nc_xcor_kernel_analytic_power_exp_get_n_comps (NcXcorKernelAnalytic *xcka);
static gdouble _nc_xcor_kernel_analytic_power_exp_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi);
static void _nc_xcor_kernel_analytic_power_exp_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);

static void
nc_xcor_kernel_analytic_power_exp_class_init (NcXcorKernelAnalyticPowerExpClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);
  NcXcorKernelAnalyticClass *parent_class = NC_XCOR_KERNEL_ANALYTIC_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_analytic_power_exp_constructed;
  object_class->finalize    = &nc_xcor_kernel_analytic_power_exp_finalize;
  model_class->set_property = &nc_xcor_kernel_analytic_power_exp_set_property;
  model_class->get_property = &nc_xcor_kernel_analytic_power_exp_get_property;

  ncm_model_class_set_name_nick (model_class, "Analytic power-exponential radial window", "AnalyticPowerExp");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelAnalyticPowerExp:chi-scale:
   *
   * Scale $\chi_0$ of the exponential cut-off, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_SCALE,
                                   g_param_spec_double ("chi-scale",
                                                        NULL,
                                                        "Cut-off scale in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 1200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticPowerExp:alpha:
   *
   * Power-law index $\alpha$ of the rise.
   */
  g_object_class_install_property (object_class,
                                   PROP_ALPHA,
                                   g_param_spec_double ("alpha",
                                                        NULL,
                                                        "Power-law index of the rise",
                                                        0.0, G_MAXDOUBLE, 2.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticPowerExp:beta:
   *
   * Stretch exponent $\beta$ of the fall.
   */
  g_object_class_install_property (object_class,
                                   PROP_BETA,
                                   g_param_spec_double ("beta",
                                                        NULL,
                                                        "Stretch exponent of the fall",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 1.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticPowerExp:chi-lower:
   *
   * Lower end $\chi_\mathrm{l}$ of the support, in Mpc. Part of the window's
   * definition: the tail is renormalized over what is kept.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_LOWER,
                                   g_param_spec_double ("chi-lower",
                                                        NULL,
                                                        "Lower end of the support in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 50.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticPowerExp:chi-upper:
   *
   * Upper end $\chi_\mathrm{u}$ of the support, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_UPPER,
                                   g_param_spec_double ("chi-upper",
                                                        NULL,
                                                        "Upper end of the support in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 4000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->get_n_comps      = &_nc_xcor_kernel_analytic_power_exp_get_n_comps;
  parent_class->eval_W_comp      = &_nc_xcor_kernel_analytic_power_exp_eval_W_comp;
  parent_class->get_comp_support = &_nc_xcor_kernel_analytic_power_exp_get_comp_support;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_kernel_analytic_power_exp_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_scale: cut-off scale $\chi_0$, in Mpc
 * @alpha: power-law index $\alpha$ of the rise
 * @beta: stretch exponent $\beta$ of the fall
 * @chi_lower: lower end of the support, in Mpc
 * @chi_upper: upper end of the support, in Mpc
 *
 * Creates a new #NcXcorKernelAnalyticPowerExp.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticPowerExp
 */
NcXcorKernelAnalyticPowerExp *
nc_xcor_kernel_analytic_power_exp_new (NcDistance *dist, NcmPowspec *ps, gdouble chi_scale, gdouble alpha, gdouble beta, gdouble chi_lower, gdouble chi_upper)
{
  NcXcorKernelAnalyticPowerExp *xckap = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_POWER_EXP,
                                                      "dist", dist,
                                                      "powspec", ps,
                                                      "chi-scale", chi_scale,
                                                      "alpha", alpha,
                                                      "beta", beta,
                                                      "chi-lower", chi_lower,
                                                      "chi-upper", chi_upper,
                                                      NULL);

  return xckap;
}

/**
 * nc_xcor_kernel_analytic_power_exp_new_full:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_scale: cut-off scale $\chi_0$, in Mpc
 * @alpha: power-law index $\alpha$ of the rise
 * @beta: stretch exponent $\beta$ of the fall
 * @chi_lower: lower end of the support, in Mpc
 * @chi_upper: upper end of the support, in Mpc
 * @sbi: a #NcmSBesselIntegrator
 *
 * Creates a new #NcXcorKernelAnalyticPowerExp carrying @sbi, as
 * nc_xcor_kernel_analytic_power_exp_new() does not. A #NcXcorKernel only
 * accepts the non-Limber modes of nc_xcor_kernel_set_l_limber() once it holds
 * an integrator, so this is the constructor to use for them.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticPowerExp
 */
NcXcorKernelAnalyticPowerExp *
nc_xcor_kernel_analytic_power_exp_new_full (NcDistance *dist, NcmPowspec *ps, gdouble chi_scale, gdouble alpha, gdouble beta, gdouble chi_lower, gdouble chi_upper, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelAnalyticPowerExp *xckap = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_POWER_EXP,
                                                      "dist", dist,
                                                      "powspec", ps,
                                                      "chi-scale", chi_scale,
                                                      "alpha", alpha,
                                                      "beta", beta,
                                                      "chi-lower", chi_lower,
                                                      "chi-upper", chi_upper,
                                                      "integrator", sbi,
                                                      NULL);

  return xckap;
}

/**
 * nc_xcor_kernel_analytic_power_exp_get_chi_scale:
 * @xckap: a #NcXcorKernelAnalyticPowerExp
 *
 * Returns: the cut-off scale $\chi_0$, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_power_exp_get_chi_scale (NcXcorKernelAnalyticPowerExp *xckap)
{
  return xckap->chi_scale;
}

/**
 * nc_xcor_kernel_analytic_power_exp_get_alpha:
 * @xckap: a #NcXcorKernelAnalyticPowerExp
 *
 * Returns: the power-law index $\alpha$.
 */
gdouble
nc_xcor_kernel_analytic_power_exp_get_alpha (NcXcorKernelAnalyticPowerExp *xckap)
{
  return xckap->alpha;
}

/**
 * nc_xcor_kernel_analytic_power_exp_get_beta:
 * @xckap: a #NcXcorKernelAnalyticPowerExp
 *
 * Returns: the stretch exponent $\beta$.
 */
gdouble
nc_xcor_kernel_analytic_power_exp_get_beta (NcXcorKernelAnalyticPowerExp *xckap)
{
  return xckap->beta;
}

static guint
_nc_xcor_kernel_analytic_power_exp_get_n_comps (NcXcorKernelAnalytic *xcka)
{
  return 1;
}

static gdouble
_nc_xcor_kernel_analytic_power_exp_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi)
{
  NcXcorKernelAnalyticPowerExp *xckap = NC_XCOR_KERNEL_ANALYTIC_POWER_EXP (xcka);
  const gdouble x                     = chi / xckap->chi_scale;

  if ((chi < xckap->chi_lower) || (chi > xckap->chi_upper))
    return 0.0;

  return pow (x, xckap->alpha) * exp (-pow (x, xckap->beta)) / xckap->norm;
}

static void
_nc_xcor_kernel_analytic_power_exp_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  NcXcorKernelAnalyticPowerExp *xckap = NC_XCOR_KERNEL_ANALYTIC_POWER_EXP (xcka);

  *chi_min = xckap->chi_lower;
  *chi_max = xckap->chi_upper;
}

