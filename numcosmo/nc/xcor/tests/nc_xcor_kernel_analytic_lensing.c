/***************************************************************************
 *            nc_xcor_kernel_analytic_lensing.c
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
 * NcXcorKernelAnalyticLensing:
 *
 * Lensing efficiency of a top-hat source distribution: an integral of the
 * source window, and so broader and smoother than it.
 *
 * \begin{equation}
 *   W(\chi) \propto \chi \int_{\max(\chi,\chi_a)}^{\chi_b}
 *   \frac{\mathrm{d}\chi'}{\chi_b - \chi_a}
 *   \left(1 - \frac{\chi}{\chi'}\right) ,
 * \end{equation}
 *
 * which for a top-hat source on $[\chi_a, \chi_b]$ is elementary:
 * $\chi[(\chi_b - \chi_a) - \chi\ln(\chi_b/\chi_a)]$ in front of the sources and
 * $\chi[(\chi_b - \chi) - \chi\ln(\chi_b/\chi)]$ among them, both over
 * $\chi_b - \chi_a$, and zero beyond. The normalization is elementary too, so
 * the whole window is closed form.
 *
 * The point of the shape is that it is an *integral* of a source distribution:
 * it is broad, peaks at roughly half the source distance, and is smooth even
 * where the source window has hard edges. A solver tuned on bin-like windows
 * meets something qualitatively different here, and its exact value is still
 * known.
 *
 * The geometric efficiency alone is used; the $(1+z)$ and prefactor a physical
 * convergence kernel carries are evolution, and by the convention of
 * #NcXcorKernelAnalytic all evolution belongs in the window's definition
 * rather than being applied on top of it.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_lensing.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelAnalyticLensing
{
  /*< private >*/
  NcXcorKernelAnalytic parent_instance;

  gdouble chi_lower;
  gdouble chi_source_lower;
  gdouble chi_source_upper;

  gdouble norm;
};

enum
{
  PROP_0,
  PROP_CHI_LOWER,
  PROP_CHI_SOURCE_LOWER,
  PROP_CHI_SOURCE_UPPER,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelAnalyticLensing, nc_xcor_kernel_analytic_lensing, NC_TYPE_XCOR_KERNEL_ANALYTIC)

static void
nc_xcor_kernel_analytic_lensing_init (NcXcorKernelAnalyticLensing *xckal)
{
  xckal->chi_lower        = 0.0;
  xckal->chi_source_lower = 0.0;
  xckal->chi_source_upper = 0.0;
  xckal->norm             = 0.0;
}

static void
nc_xcor_kernel_analytic_lensing_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticLensing *xckal = NC_XCOR_KERNEL_ANALYTIC_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_LENSING (object));

  switch (prop_id)
  {
    case PROP_CHI_LOWER:
      xckal->chi_lower = g_value_get_double (value);
      break;
    case PROP_CHI_SOURCE_LOWER:
      xckal->chi_source_lower = g_value_get_double (value);
      break;
    case PROP_CHI_SOURCE_UPPER:
      xckal->chi_source_upper = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_lensing_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticLensing *xckal = NC_XCOR_KERNEL_ANALYTIC_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_LENSING (object));

  switch (prop_id)
  {
    case PROP_CHI_LOWER:
      g_value_set_double (value, xckal->chi_lower);
      break;
    case PROP_CHI_SOURCE_LOWER:
      g_value_set_double (value, xckal->chi_source_lower);
      break;
    case PROP_CHI_SOURCE_UPPER:
      g_value_set_double (value, xckal->chi_source_upper);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* Unnormalized efficiency: chi times the source integral. */
static gdouble
_lensing_shape (NcXcorKernelAnalyticLensing *xckal, gdouble chi)
{
  const gdouble a = xckal->chi_source_lower;
  const gdouble b = xckal->chi_source_upper;

  if (chi >= b)
    return 0.0;

  if (chi <= a)
    return chi * ((b - a) - chi * log (b / a)) / (b - a);

  return chi * ((b - chi) - chi * log (b / chi)) / (b - a);
}

static void
nc_xcor_kernel_analytic_lensing_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_lensing_parent_class)->constructed (object);
  {
    NcXcorKernelAnalyticLensing *xckal = NC_XCOR_KERNEL_ANALYTIC_LENSING (object);
    const gdouble a                    = xckal->chi_source_lower;
    const gdouble b                    = xckal->chi_source_upper;
    const gdouble c                    = xckal->chi_lower;

    if (!(b > a))
      g_error ("nc_xcor_kernel_analytic_lensing_constructed: empty source bin, "
               "chi-source-lower %.17g >= chi-source-upper %.17g.", a, b);

    if (!((c > 0.0) && (c < a)))
      g_error ("nc_xcor_kernel_analytic_lensing_constructed: chi-lower %.17g must lie "
               "strictly between the observer and the nearest source at %.17g. It is "
               "part of the window's definition, and zero would leave the Limber band "
               "unbounded.", c, a);

    {
      /* Both regimes integrate in closed form; the second uses
       * int chi^2 ln(chi) dchi = chi^3 (3 ln chi - 1) / 9. */
      const gdouble L     = log (b / a);
      const gdouble front = ((b - a) * (a * a - c * c) / 2.0 - L * (a * a * a - c * c * c) / 3.0) / (b - a);
      const gdouble among = (b * (b * b - a * a) / 2.0
                             - (b * b * b - a * a * a) / 3.0
                             + a * a * a * L / 3.0
                             + (a * a * a - b * b * b) / 9.0) / (b - a);

      xckal->norm = front + among;
    }

    if (!(xckal->norm > 0.0))
      g_error ("nc_xcor_kernel_analytic_lensing_constructed: window has zero mass for "
               "chi-lower %.17g, sources on [%.17g, %.17g] Mpc.", c, a, b);
  }
}

static void
nc_xcor_kernel_analytic_lensing_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_lensing_parent_class)->finalize (object);
}

static guint _nc_xcor_kernel_analytic_lensing_get_n_comps (NcXcorKernelAnalytic *xcka);
static gdouble _nc_xcor_kernel_analytic_lensing_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi);
static void _nc_xcor_kernel_analytic_lensing_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);

static void
nc_xcor_kernel_analytic_lensing_class_init (NcXcorKernelAnalyticLensingClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);
  NcXcorKernelAnalyticClass *parent_class = NC_XCOR_KERNEL_ANALYTIC_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_analytic_lensing_constructed;
  object_class->finalize    = &nc_xcor_kernel_analytic_lensing_finalize;
  model_class->set_property = &nc_xcor_kernel_analytic_lensing_set_property;
  model_class->get_property = &nc_xcor_kernel_analytic_lensing_get_property;

  ncm_model_class_set_name_nick (model_class, "Analytic lensing-efficiency radial window", "AnalyticLensing");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelAnalyticLensing:chi-lower:
   *
   * Observer-side end $\chi_\mathrm{l}$ of the support, in Mpc. The efficiency
   * vanishes linearly at the observer, but the support must start strictly
   * inside: at $\chi = 0$ the Limber band $\nu/\chi$ is unbounded.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_LOWER,
                                   g_param_spec_double ("chi-lower",
                                                        NULL,
                                                        "Observer-side end of the support in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 50.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticLensing:chi-source-lower:
   *
   * Near edge $\chi_a$ of the top-hat source distribution, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_SOURCE_LOWER,
                                   g_param_spec_double ("chi-source-lower",
                                                        NULL,
                                                        "Near edge of the source bin in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 2000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticLensing:chi-source-upper:
   *
   * Far edge $\chi_b$ of the top-hat source distribution, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_SOURCE_UPPER,
                                   g_param_spec_double ("chi-source-upper",
                                                        NULL,
                                                        "Far edge of the source bin in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 3000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->get_n_comps      = &_nc_xcor_kernel_analytic_lensing_get_n_comps;
  parent_class->eval_W_comp      = &_nc_xcor_kernel_analytic_lensing_eval_W_comp;
  parent_class->get_comp_support = &_nc_xcor_kernel_analytic_lensing_get_comp_support;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_kernel_analytic_lensing_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_lower: observer-side end of the support, in Mpc
 * @chi_source_lower: near edge $\chi_a$ of the source bin, in Mpc
 * @chi_source_upper: far edge $\chi_b$ of the source bin, in Mpc
 *
 * Creates a new #NcXcorKernelAnalyticLensing.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticLensing
 */
NcXcorKernelAnalyticLensing *
nc_xcor_kernel_analytic_lensing_new (NcDistance *dist, NcmPowspec *ps, gdouble chi_lower, gdouble chi_source_lower, gdouble chi_source_upper)
{
  NcXcorKernelAnalyticLensing *xckal = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_LENSING,
                                                     "dist", dist,
                                                     "powspec", ps,
                                                     "chi-lower", chi_lower,
                                                     "chi-source-lower", chi_source_lower,
                                                     "chi-source-upper", chi_source_upper,
                                                     NULL);

  return xckal;
}

/**
 * nc_xcor_kernel_analytic_lensing_new_full:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_lower: observer-side end of the support, in Mpc
 * @chi_source_lower: near edge $\chi_a$ of the source bin, in Mpc
 * @chi_source_upper: far edge $\chi_b$ of the source bin, in Mpc
 * @sbi: a #NcmSBesselIntegrator
 *
 * Creates a new #NcXcorKernelAnalyticLensing carrying @sbi, as
 * nc_xcor_kernel_analytic_lensing_new() does not. A #NcXcorKernel only accepts
 * the non-Limber modes of nc_xcor_kernel_set_l_limber() once it holds an
 * integrator, so this is the constructor to use for them.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticLensing
 */
NcXcorKernelAnalyticLensing *
nc_xcor_kernel_analytic_lensing_new_full (NcDistance *dist, NcmPowspec *ps, gdouble chi_lower, gdouble chi_source_lower, gdouble chi_source_upper, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelAnalyticLensing *xckal = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_LENSING,
                                                     "dist", dist,
                                                     "powspec", ps,
                                                     "chi-lower", chi_lower,
                                                     "chi-source-lower", chi_source_lower,
                                                     "chi-source-upper", chi_source_upper,
                                                     "integrator", sbi,
                                                     NULL);

  return xckal;
}

/**
 * nc_xcor_kernel_analytic_lensing_get_chi_source_lower:
 * @xckal: a #NcXcorKernelAnalyticLensing
 *
 * Returns: the near edge $\chi_a$ of the source bin, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_lensing_get_chi_source_lower (NcXcorKernelAnalyticLensing *xckal)
{
  return xckal->chi_source_lower;
}

/**
 * nc_xcor_kernel_analytic_lensing_get_chi_source_upper:
 * @xckal: a #NcXcorKernelAnalyticLensing
 *
 * Returns: the far edge $\chi_b$ of the source bin, in Mpc.
 */
gdouble
nc_xcor_kernel_analytic_lensing_get_chi_source_upper (NcXcorKernelAnalyticLensing *xckal)
{
  return xckal->chi_source_upper;
}

static guint
_nc_xcor_kernel_analytic_lensing_get_n_comps (NcXcorKernelAnalytic *xcka)
{
  return 1;
}

static gdouble
_nc_xcor_kernel_analytic_lensing_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi)
{
  NcXcorKernelAnalyticLensing *xckal = NC_XCOR_KERNEL_ANALYTIC_LENSING (xcka);

  if ((chi < xckal->chi_lower) || (chi > xckal->chi_source_upper))
    return 0.0;

  return _lensing_shape (xckal, chi) / xckal->norm;
}

static void
_nc_xcor_kernel_analytic_lensing_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  NcXcorKernelAnalyticLensing *xckal = NC_XCOR_KERNEL_ANALYTIC_LENSING (xcka);

  *chi_min = xckal->chi_lower;
  *chi_max = xckal->chi_source_upper;
}

