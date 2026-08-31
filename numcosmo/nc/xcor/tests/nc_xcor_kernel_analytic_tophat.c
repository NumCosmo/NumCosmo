/***************************************************************************
 *            nc_xcor_kernel_analytic_tophat.c
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
 * NcXcorKernelAnalyticTophat:
 *
 * Sharp-edged top-hat radial window in comoving distance.
 *
 * \begin{equation}
 *   W(\chi) = \frac{1}{\chi_\mathrm{u} - \chi_\mathrm{l}},
 *   \qquad \chi \in [\chi_\mathrm{l}, \chi_\mathrm{u}] ,
 * \end{equation}
 *
 * and zero outside. The discontinuities at the two edges are the point of the
 * shape: it is the piecewise stress case, where a solver that assumes a smooth
 * integrand has nothing to converge to and panel edges must fall on the
 * breakpoints. See #NcXcorKernelRadial for the unit and normalization
 * conventions.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_tophat.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelAnalyticTophat
{
  /*< private >*/
  NcXcorKernelRadial parent_instance;

  gdouble chi_lower;
  gdouble chi_upper;
};

enum
{
  PROP_0,
  PROP_CHI_LOWER,
  PROP_CHI_UPPER,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelAnalyticTophat, nc_xcor_kernel_analytic_tophat, NC_TYPE_XCOR_KERNEL_RADIAL)

static void
nc_xcor_kernel_analytic_tophat_init (NcXcorKernelAnalyticTophat *xckat)
{
  xckat->chi_lower = 0.0;
  xckat->chi_upper = 0.0;
}

static void
nc_xcor_kernel_analytic_tophat_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticTophat *xckat = NC_XCOR_KERNEL_ANALYTIC_TOPHAT (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_TOPHAT (object));

  switch (prop_id)
  {
    case PROP_CHI_LOWER:
      xckat->chi_lower = g_value_get_double (value);
      break;
    case PROP_CHI_UPPER:
      xckat->chi_upper = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_tophat_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticTophat *xckat = NC_XCOR_KERNEL_ANALYTIC_TOPHAT (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_TOPHAT (object));

  switch (prop_id)
  {
    case PROP_CHI_LOWER:
      g_value_set_double (value, xckat->chi_lower);
      break;
    case PROP_CHI_UPPER:
      g_value_set_double (value, xckat->chi_upper);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_tophat_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_tophat_parent_class)->constructed (object);
  {
    NcXcorKernelAnalyticTophat *xckat = NC_XCOR_KERNEL_ANALYTIC_TOPHAT (object);

    if (!(xckat->chi_upper > xckat->chi_lower))
      g_error ("nc_xcor_kernel_analytic_tophat_constructed: empty support, "
               "chi-lower %.17g >= chi-upper %.17g.", xckat->chi_lower, xckat->chi_upper);
  }
}

static void
nc_xcor_kernel_analytic_tophat_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_tophat_parent_class)->finalize (object);
}

static guint _nc_xcor_kernel_analytic_tophat_get_n_comps (NcXcorKernelRadial *xcka);
static gdouble _nc_xcor_kernel_analytic_tophat_eval_W_comp (NcXcorKernelRadial *xcka, guint comp, gdouble chi);
static void _nc_xcor_kernel_analytic_tophat_get_comp_support (NcXcorKernelRadial *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);

static void
nc_xcor_kernel_analytic_tophat_class_init (NcXcorKernelAnalyticTophatClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class            = NCM_MODEL_CLASS (klass);
  NcXcorKernelRadialClass *parent_class = NC_XCOR_KERNEL_RADIAL_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_analytic_tophat_constructed;
  object_class->finalize    = &nc_xcor_kernel_analytic_tophat_finalize;
  model_class->set_property = &nc_xcor_kernel_analytic_tophat_set_property;
  model_class->get_property = &nc_xcor_kernel_analytic_tophat_get_property;

  ncm_model_class_set_name_nick (model_class, "Analytic top-hat radial window", "AnalyticTophat");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelAnalyticTophat:chi-lower:
   *
   * Lower edge $\chi_\mathrm{l}$ of the window, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_LOWER,
                                   g_param_spec_double ("chi-lower",
                                                        NULL,
                                                        "Lower edge in Mpc",
                                                        0.0, G_MAXDOUBLE, 1000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticTophat:chi-upper:
   *
   * Upper edge $\chi_\mathrm{u}$ of the window, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_UPPER,
                                   g_param_spec_double ("chi-upper",
                                                        NULL,
                                                        "Upper edge in Mpc",
                                                        0.0, G_MAXDOUBLE, 2000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->get_n_comps      = &_nc_xcor_kernel_analytic_tophat_get_n_comps;
  parent_class->eval_W_comp      = &_nc_xcor_kernel_analytic_tophat_eval_W_comp;
  parent_class->get_comp_support = &_nc_xcor_kernel_analytic_tophat_get_comp_support;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_kernel_analytic_tophat_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_lower: lower edge $\chi_\mathrm{l}$, in Mpc
 * @chi_upper: upper edge $\chi_\mathrm{u}$, in Mpc
 *
 * Creates a new #NcXcorKernelAnalyticTophat.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticTophat
 */
NcXcorKernelAnalyticTophat *
nc_xcor_kernel_analytic_tophat_new (NcDistance *dist, NcmPowspec *ps, gdouble chi_lower, gdouble chi_upper)
{
  NcXcorKernelAnalyticTophat *xckat = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_TOPHAT,
                                                    "dist", dist,
                                                    "powspec", ps,
                                                    "chi-lower", chi_lower,
                                                    "chi-upper", chi_upper,
                                                    NULL);

  return xckat;
}

/**
 * nc_xcor_kernel_analytic_tophat_new_full:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_lower: lower edge $\chi_\mathrm{l}$, in Mpc
 * @chi_upper: upper edge $\chi_\mathrm{u}$, in Mpc
 * @sbi: a #NcmSBesselIntegrator
 *
 * Creates a new #NcXcorKernelAnalyticTophat carrying @sbi, as
 * nc_xcor_kernel_analytic_tophat_new() does not. A #NcXcorKernel only accepts
 * the non-Limber modes of nc_xcor_kernel_set_l_limber() once it holds an
 * integrator, so this is the constructor to use for them.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticTophat
 */
NcXcorKernelAnalyticTophat *
nc_xcor_kernel_analytic_tophat_new_full (NcDistance *dist, NcmPowspec *ps, gdouble chi_lower, gdouble chi_upper, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelAnalyticTophat *xckat = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_TOPHAT,
                                                    "dist", dist,
                                                    "powspec", ps,
                                                    "chi-lower", chi_lower,
                                                    "chi-upper", chi_upper,
                                                    "integrator", sbi,
                                                    NULL);

  return xckat;
}

/**
 * nc_xcor_kernel_analytic_tophat_get_chi_lower:
 * @xckat: a #NcXcorKernelAnalyticTophat
 *
 * Returns: the lower edge $\chi_\mathrm{l}$, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_tophat_get_chi_lower (NcXcorKernelAnalyticTophat *xckat)
{
  return xckat->chi_lower;
}

/**
 * nc_xcor_kernel_analytic_tophat_get_chi_upper:
 * @xckat: a #NcXcorKernelAnalyticTophat
 *
 * Returns: the upper edge $\chi_\mathrm{u}$, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_tophat_get_chi_upper (NcXcorKernelAnalyticTophat *xckat)
{
  return xckat->chi_upper;
}

static guint
_nc_xcor_kernel_analytic_tophat_get_n_comps (NcXcorKernelRadial *xcka)
{
  return 1;
}

static gdouble
_nc_xcor_kernel_analytic_tophat_eval_W_comp (NcXcorKernelRadial *xcka, guint comp, gdouble chi)
{
  NcXcorKernelAnalyticTophat *xckat = NC_XCOR_KERNEL_ANALYTIC_TOPHAT (xcka);

  if ((chi < xckat->chi_lower) || (chi > xckat->chi_upper))
    return 0.0;

  return 1.0 / (xckat->chi_upper - xckat->chi_lower);
}

static void
_nc_xcor_kernel_analytic_tophat_get_comp_support (NcXcorKernelRadial *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  NcXcorKernelAnalyticTophat *xckat = NC_XCOR_KERNEL_ANALYTIC_TOPHAT (xcka);

  *chi_min = xckat->chi_lower;
  *chi_max = xckat->chi_upper;
}

