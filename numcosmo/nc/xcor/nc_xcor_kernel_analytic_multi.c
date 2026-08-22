/***************************************************************************
 *            nc_xcor_kernel_analytic_multi.c
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
 * NcXcorKernelAnalyticMulti:
 *
 * Multimodal radial window: a weighted sum of Gaussian bumps in comoving
 * distance.
 *
 * \begin{equation}
 *   W(\chi) = \frac{1}{N} \sum_i w_i
 *   \exp\left[-\frac{(\chi - \chi_i)^2}{2\sigma_i^2}\right] ,
 * \end{equation}
 *
 * with $N$ fixed in closed form so that $\int W \mathrm{d}\chi = 1$ over the
 * support actually kept.
 *
 * Each bump is truncated at $n\sigma_i$, and the bumps are grouped into
 * components by whether those intervals meet: bumps whose intervals overlap
 * form **one** component spanning their union, and a bump whose interval
 * touches no other is a component on its own. A component is then integrated
 * over exactly its own interval, so every boundary of the window is a boundary
 * of an integration domain -- see #NcXcorKernelAnalytic.
 *
 * The grouping is by support, not by peak. Bumps that merge into one
 * continuous stretch stay together and each is evaluated across the whole
 * stretch rather than cut at its own $n\sigma_i$: cutting there would put a
 * step inside a domain, which is precisely what the split exists to avoid.
 * Only a real gap, where the window is identically zero, separates components.
 *
 * This gives the protocol a window with several scales at once -- useful for a
 * two-population $\mathrm{d}n/\mathrm{d}z$, and, with the bumps pushed apart,
 * for a genuinely disconnected support.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/nc_xcor_kernel_analytic_multi.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcXcorKernelAnalyticMultiComp
{
  gdouble chi_min;
  gdouble chi_max;
  guint first;
  guint len;
} NcXcorKernelAnalyticMultiComp;

struct _NcXcorKernelAnalyticMulti
{
  /*< private >*/
  NcXcorKernelAnalytic parent_instance;

  NcmVector *chi_mean;
  NcmVector *chi_sigma;
  NcmVector *weight;
  gdouble n_sigma;

  /* Bumps sorted by the lower end of their support, and the components they
   * group into. */
  GArray *order;
  NcXcorKernelAnalyticMultiComp comps[NC_XCOR_KERNEL_ANALYTIC_MAX_COMPS];
  guint n_comps;
  gdouble norm;
};

enum
{
  PROP_0,
  PROP_CHI_MEAN,
  PROP_CHI_SIGMA,
  PROP_WEIGHT,
  PROP_N_SIGMA,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelAnalyticMulti, nc_xcor_kernel_analytic_multi, NC_TYPE_XCOR_KERNEL_ANALYTIC)

static void
nc_xcor_kernel_analytic_multi_init (NcXcorKernelAnalyticMulti *xckam)
{
  xckam->chi_mean  = NULL;
  xckam->chi_sigma = NULL;
  xckam->weight    = NULL;
  xckam->n_sigma   = 0.0;
  xckam->order     = g_array_new (FALSE, FALSE, sizeof (guint));
  xckam->n_comps   = 0;
  xckam->norm      = 0.0;
}

static void
nc_xcor_kernel_analytic_multi_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticMulti *xckam = NC_XCOR_KERNEL_ANALYTIC_MULTI (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_MULTI (object));

  switch (prop_id)
  {
    case PROP_CHI_MEAN:
      ncm_vector_clear (&xckam->chi_mean);
      xckam->chi_mean = g_value_dup_object (value);
      break;
    case PROP_CHI_SIGMA:
      ncm_vector_clear (&xckam->chi_sigma);
      xckam->chi_sigma = g_value_dup_object (value);
      break;
    case PROP_WEIGHT:
      ncm_vector_clear (&xckam->weight);
      xckam->weight = g_value_dup_object (value);
      break;
    case PROP_N_SIGMA:
      xckam->n_sigma = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_analytic_multi_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelAnalyticMulti *xckam = NC_XCOR_KERNEL_ANALYTIC_MULTI (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_ANALYTIC_MULTI (object));

  switch (prop_id)
  {
    case PROP_CHI_MEAN:
      g_value_set_object (value, xckam->chi_mean);
      break;
    case PROP_CHI_SIGMA:
      g_value_set_object (value, xckam->chi_sigma);
      break;
    case PROP_WEIGHT:
      g_value_set_object (value, xckam->weight);
      break;
    case PROP_N_SIGMA:
      g_value_set_double (value, xckam->n_sigma);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* Support of a single bump, truncated at n_sigma and clipped at the observer. */
static void
_bump_support (NcXcorKernelAnalyticMulti *xckam, guint i, gdouble *lo, gdouble *hi)
{
  const gdouble mu    = ncm_vector_get (xckam->chi_mean, i);
  const gdouble sigma = ncm_vector_get (xckam->chi_sigma, i);

  *lo = GSL_MAX (0.0, mu - xckam->n_sigma * sigma);
  *hi = mu + xckam->n_sigma * sigma;
}

/* Integral of one bump over [lo, hi], in closed form. */
static gdouble
_bump_mass (NcXcorKernelAnalyticMulti *xckam, guint i, gdouble lo, gdouble hi)
{
  const gdouble mu    = ncm_vector_get (xckam->chi_mean, i);
  const gdouble sigma = ncm_vector_get (xckam->chi_sigma, i);
  const gdouble w     = ncm_vector_get (xckam->weight, i);
  const gdouble s2    = M_SQRT2 * sigma;

  return w * 0.5 * sqrt (2.0 * M_PI) * sigma * (erf ((hi - mu) / s2) - erf ((lo - mu) / s2));
}

/*
 * Groups the bumps into components. Sorting by the lower end of the support
 * turns "do these intervals meet?" into a single sweep, and a merged component
 * is a contiguous run of the sorted order, so a component needs only its first
 * index and length.
 */
static void
_nc_xcor_kernel_analytic_multi_build_comps (NcXcorKernelAnalyticMulti *xckam)
{
  const guint n = ncm_vector_len (xckam->chi_mean);
  GArray *lows  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n);
  GArray *perm  = g_array_sized_new (FALSE, FALSE, sizeof (gsize), n);
  guint i;

  g_array_set_size (lows, n);
  g_array_set_size (perm, n);
  g_array_set_size (xckam->order, n);

  for (i = 0; i < n; i++)
  {
    gdouble lo, hi;

    _bump_support (xckam, i, &lo, &hi);
    g_array_index (lows, gdouble, i) = lo;
  }

  gsl_sort_index ((gsize *) perm->data, (gdouble *) lows->data, 1, n);

  for (i = 0; i < n; i++)
    g_array_index (xckam->order, guint, i) = (guint) g_array_index (perm, gsize, i);

  xckam->n_comps = 0;

  for (i = 0; i < n; i++)
  {
    const guint b = g_array_index (xckam->order, guint, i);
    gdouble lo, hi;

    _bump_support (xckam, b, &lo, &hi);

    if ((xckam->n_comps > 0) && (lo <= xckam->comps[xckam->n_comps - 1].chi_max))
    {
      NcXcorKernelAnalyticMultiComp *cur = &xckam->comps[xckam->n_comps - 1];

      cur->chi_max = GSL_MAX (cur->chi_max, hi);
      cur->len++;
    }
    else
    {
      if (xckam->n_comps == NC_XCOR_KERNEL_ANALYTIC_MAX_COMPS)
        g_error ("nc_xcor_kernel_analytic_multi: the bumps fall into more than %d disjoint "
                 "groups, which is more components than a kernel carries. Widen n-sigma, "
                 "move the bumps closer, or use fewer of them.",
                 NC_XCOR_KERNEL_ANALYTIC_MAX_COMPS);

      xckam->comps[xckam->n_comps].chi_min = lo;
      xckam->comps[xckam->n_comps].chi_max = hi;
      xckam->comps[xckam->n_comps].first   = i;
      xckam->comps[xckam->n_comps].len     = 1;
      xckam->n_comps++;
    }
  }

  g_array_unref (lows);
  g_array_unref (perm);
}

static void
nc_xcor_kernel_analytic_multi_constructed (GObject *object)
{
  {
    NcXcorKernelAnalyticMulti *xckam = NC_XCOR_KERNEL_ANALYTIC_MULTI (object);
    guint n, c;

    if ((xckam->chi_mean == NULL) || (xckam->chi_sigma == NULL) || (xckam->weight == NULL))
      g_error ("nc_xcor_kernel_analytic_multi_constructed: chi-mean, chi-sigma and weight "
               "must all be given at construction time.");

    n = ncm_vector_len (xckam->chi_mean);

    if ((n == 0) || (ncm_vector_len (xckam->chi_sigma) != n) || (ncm_vector_len (xckam->weight) != n))
      g_error ("nc_xcor_kernel_analytic_multi_constructed: chi-mean, chi-sigma and weight "
               "must have the same non-zero length, got %u, %u and %u.",
               n, ncm_vector_len (xckam->chi_sigma), ncm_vector_len (xckam->weight));

    _nc_xcor_kernel_analytic_multi_build_comps (xckam);

    xckam->norm = 0.0;

    for (c = 0; c < xckam->n_comps; c++)
    {
      const NcXcorKernelAnalyticMultiComp *comp = &xckam->comps[c];
      guint j;

      for (j = 0; j < comp->len; j++)
        xckam->norm += _bump_mass (xckam, g_array_index (xckam->order, guint, comp->first + j),
                                   comp->chi_min, comp->chi_max);
    }

    if (!(xckam->norm > 0.0))
      g_error ("nc_xcor_kernel_analytic_multi_constructed: window has zero mass; check the "
               "weights and widths.");
  }

  /* Chain up : end. The base builds one component per interval, so the grouping
   * above has to be settled first. */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_multi_parent_class)->constructed (object);
}

static void
nc_xcor_kernel_analytic_multi_dispose (GObject *object)
{
  NcXcorKernelAnalyticMulti *xckam = NC_XCOR_KERNEL_ANALYTIC_MULTI (object);

  ncm_vector_clear (&xckam->chi_mean);
  ncm_vector_clear (&xckam->chi_sigma);
  ncm_vector_clear (&xckam->weight);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_multi_parent_class)->dispose (object);
}

static void
nc_xcor_kernel_analytic_multi_finalize (GObject *object)
{
  NcXcorKernelAnalyticMulti *xckam = NC_XCOR_KERNEL_ANALYTIC_MULTI (object);

  g_clear_pointer (&xckam->order, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_analytic_multi_parent_class)->finalize (object);
}

static guint _nc_xcor_kernel_analytic_multi_get_n_comps (NcXcorKernelAnalytic *xcka);
static gdouble _nc_xcor_kernel_analytic_multi_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi);
static void _nc_xcor_kernel_analytic_multi_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);

static void
nc_xcor_kernel_analytic_multi_class_init (NcXcorKernelAnalyticMultiClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);
  NcXcorKernelAnalyticClass *parent_class = NC_XCOR_KERNEL_ANALYTIC_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_analytic_multi_constructed;
  object_class->dispose     = &nc_xcor_kernel_analytic_multi_dispose;
  object_class->finalize    = &nc_xcor_kernel_analytic_multi_finalize;
  model_class->set_property = &nc_xcor_kernel_analytic_multi_set_property;
  model_class->get_property = &nc_xcor_kernel_analytic_multi_get_property;

  ncm_model_class_set_name_nick (model_class, "Analytic multimodal radial window", "AnalyticMulti");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelAnalyticMulti:chi-mean:
   *
   * Centres $\chi_i$ of the bumps, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_MEAN,
                                   g_param_spec_object ("chi-mean",
                                                        NULL,
                                                        "Bump centres in Mpc",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticMulti:chi-sigma:
   *
   * Standard deviations $\sigma_i$ of the bumps, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_SIGMA,
                                   g_param_spec_object ("chi-sigma",
                                                        NULL,
                                                        "Bump standard deviations in Mpc",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticMulti:weight:
   *
   * Relative weights $w_i$ of the bumps. Only their ratios matter: the window
   * is renormalized to unit integral regardless of the overall scale.
   */
  g_object_class_install_property (object_class,
                                   PROP_WEIGHT,
                                   g_param_spec_object ("weight",
                                                        NULL,
                                                        "Relative bump weights",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelAnalyticMulti:n-sigma:
   *
   * Truncation half-width of each bump, in units of its own $\sigma_i$. Part of
   * the window's definition, and also what decides which bumps share a
   * component: it sets the intervals whose overlaps are merged.
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

  parent_class->get_n_comps      = &_nc_xcor_kernel_analytic_multi_get_n_comps;
  parent_class->eval_W_comp      = &_nc_xcor_kernel_analytic_multi_eval_W_comp;
  parent_class->get_comp_support = &_nc_xcor_kernel_analytic_multi_get_comp_support;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_kernel_analytic_multi_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_mean: (in): bump centres $\chi_i$, in Mpc
 * @chi_sigma: (in): bump standard deviations $\sigma_i$, in Mpc
 * @weight: (in): relative bump weights $w_i$
 * @n_sigma: truncation half-width of each bump, in units of its own $\sigma_i$
 *
 * Creates a new #NcXcorKernelAnalyticMulti. The three vectors must have the
 * same length.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticMulti
 */
NcXcorKernelAnalyticMulti *
nc_xcor_kernel_analytic_multi_new (NcDistance *dist, NcmPowspec *ps, NcmVector *chi_mean, NcmVector *chi_sigma, NcmVector *weight, gdouble n_sigma)
{
  NcXcorKernelAnalyticMulti *xckam = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_MULTI,
                                                   "dist", dist,
                                                   "powspec", ps,
                                                   "chi-mean", chi_mean,
                                                   "chi-sigma", chi_sigma,
                                                   "weight", weight,
                                                   "n-sigma", n_sigma,
                                                   NULL);

  return xckam;
}

/**
 * nc_xcor_kernel_analytic_multi_new_full:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @chi_mean: (in): bump centres $\chi_i$, in Mpc
 * @chi_sigma: (in): bump standard deviations $\sigma_i$, in Mpc
 * @weight: (in): relative bump weights $w_i$
 * @n_sigma: truncation half-width of each bump, in units of its own $\sigma_i$
 * @sbi: a #NcmSBesselIntegrator
 *
 * Creates a new #NcXcorKernelAnalyticMulti carrying @sbi, as
 * nc_xcor_kernel_analytic_multi_new() does not. A #NcXcorKernel only accepts
 * the non-Limber modes of nc_xcor_kernel_set_l_limber() once it holds an
 * integrator, so this is the constructor to use for them.
 *
 * Returns: (transfer full): a new #NcXcorKernelAnalyticMulti
 */
NcXcorKernelAnalyticMulti *
nc_xcor_kernel_analytic_multi_new_full (NcDistance *dist, NcmPowspec *ps, NcmVector *chi_mean, NcmVector *chi_sigma, NcmVector *weight, gdouble n_sigma, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelAnalyticMulti *xckam = g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_MULTI,
                                                   "dist", dist,
                                                   "powspec", ps,
                                                   "chi-mean", chi_mean,
                                                   "chi-sigma", chi_sigma,
                                                   "weight", weight,
                                                   "n-sigma", n_sigma,
                                                   "integrator", sbi,
                                                   NULL);

  return xckam;
}

/**
 * nc_xcor_kernel_analytic_multi_get_n_bumps:
 * @xckam: a #NcXcorKernelAnalyticMulti
 *
 * Gets the number of Gaussian bumps. This is not the number of components:
 * bumps whose supports overlap share one.
 *
 * Returns: the number of bumps.
 */
guint
nc_xcor_kernel_analytic_multi_get_n_bumps (NcXcorKernelAnalyticMulti *xckam)
{
  return ncm_vector_len (xckam->chi_mean);
}

/**
 * nc_xcor_kernel_analytic_multi_get_n_sigma:
 * @xckam: a #NcXcorKernelAnalyticMulti
 *
 * Returns: the truncation half-width, in units of each bump's own $\sigma_i$.
 */
gdouble
nc_xcor_kernel_analytic_multi_get_n_sigma (NcXcorKernelAnalyticMulti *xckam)
{
  return xckam->n_sigma;
}

/**
 * nc_xcor_kernel_analytic_multi_peek_chi_mean:
 * @xckam: a #NcXcorKernelAnalyticMulti
 *
 * Returns: (transfer none): the bump centres, in Mpc.
 */
NcmVector *
nc_xcor_kernel_analytic_multi_peek_chi_mean (NcXcorKernelAnalyticMulti *xckam)
{
  return xckam->chi_mean;
}

/**
 * nc_xcor_kernel_analytic_multi_peek_chi_sigma:
 * @xckam: a #NcXcorKernelAnalyticMulti
 *
 * Returns: (transfer none): the bump standard deviations, in Mpc.
 */
NcmVector *
nc_xcor_kernel_analytic_multi_peek_chi_sigma (NcXcorKernelAnalyticMulti *xckam)
{
  return xckam->chi_sigma;
}

/**
 * nc_xcor_kernel_analytic_multi_peek_weight:
 * @xckam: a #NcXcorKernelAnalyticMulti
 *
 * Returns: (transfer none): the relative bump weights.
 */
NcmVector *
nc_xcor_kernel_analytic_multi_peek_weight (NcXcorKernelAnalyticMulti *xckam)
{
  return xckam->weight;
}

static guint
_nc_xcor_kernel_analytic_multi_get_n_comps (NcXcorKernelAnalytic *xcka)
{
  return NC_XCOR_KERNEL_ANALYTIC_MULTI (xcka)->n_comps;
}

static gdouble
_nc_xcor_kernel_analytic_multi_eval_W_comp (NcXcorKernelAnalytic *xcka, guint comp, gdouble chi)
{
  NcXcorKernelAnalyticMulti *xckam           = NC_XCOR_KERNEL_ANALYTIC_MULTI (xcka);
  const NcXcorKernelAnalyticMultiComp *cinfo = &xckam->comps[comp];
  gdouble sum                                = 0.0;
  guint j;

  if ((chi < cinfo->chi_min) || (chi > cinfo->chi_max))
    return 0.0;

  /* Every bump of the group is evaluated across the whole stretch: their own
   * n_sigma cuts fall inside it, and honouring them there would put steps in
   * the middle of an integration domain. */
  for (j = 0; j < cinfo->len; j++)
  {
    const guint i       = g_array_index (xckam->order, guint, cinfo->first + j);
    const gdouble mu    = ncm_vector_get (xckam->chi_mean, i);
    const gdouble sigma = ncm_vector_get (xckam->chi_sigma, i);
    const gdouble w     = ncm_vector_get (xckam->weight, i);

    sum += w * exp (-0.5 * gsl_pow_2 ((chi - mu) / sigma));
  }

  return sum / xckam->norm;
}

static void
_nc_xcor_kernel_analytic_multi_get_comp_support (NcXcorKernelAnalytic *xcka, guint comp, gdouble *chi_min, gdouble *chi_max)
{
  NcXcorKernelAnalyticMulti *xckam = NC_XCOR_KERNEL_ANALYTIC_MULTI (xcka);

  *chi_min = xckam->comps[comp].chi_min;
  *chi_max = xckam->comps[comp].chi_max;
}

