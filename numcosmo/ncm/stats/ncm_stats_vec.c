/***************************************************************************
 *            ncm_stats_vec.c
 *
 *  Fri August 02 13:41:01 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_vec.c
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmStatsVec:
 *
 * Online weighted statistics for vectors.
 *
 * Maintains the weighted mean, variance, covariance, quantiles, and
 * autocorrelation diagnostics of appended samples.
 *
 * The mean is updated as $$\bar{x}_n = \bar{x}_{n-1} + (x_n -
 * \bar{x}_{n-1})\frac{w_n}{W_n},$$ where $\bar{x}_n$ is the mean of the first
 * $n$ elements, $x_n$ the $n$-th element, $w_n$ the $n$-th weight and $W_n$ the
 * sum of the first $n$ weights.
 *
 * The variance follows from $$M_n = M_{n-1} + (x_n -
 * \bar{x}_{n-1})^2w_n\frac{W_{n-1}}{W_n},$$ with $$V_n =
 * \frac{M_n}{W^\text{bias}_{n}}, \quad W^\text{bias}_{n} \equiv
 * \frac{W_n^2 - \sum^n_iw_i^2}{W_n},$$ where $W^\text{bias}_{n}$ is the bias
 * corrected weight.
 *
 * The covariance follows from $$N(x,y)_n = N(x,y)_{n-1} + (x_n -
 * \bar{x}_n)(y_n - \bar{y}_{n-1})w_n,$$ with $$Cov(x,y)_n =
 * \frac{N(x,y)_n}{W^\text{bias}_{n}}.$$
 *
 * # Using a NcmStatsVec. #
 * |[<!-- language="C" -->
 *
 * // One dimensional NcmStatsVec computing mean and variance.
 * NcmStatsVec *svec = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
 *
 * // Set and update three values of the single random variable.
 * ncm_stats_vec_set (svec, 0, 1.0);
 * ncm_stats_vec_update (svec);
 * ncm_stats_vec_set (svec, 0, 2.0);
 * ncm_stats_vec_update (svec);
 * ncm_stats_vec_set (svec, 0, 1.5);
 * ncm_stats_vec_update (svec);
 *
 * {
 *   gdouble mean = ncm_stats_vec_get_mean (svec, 0);
 *   gdouble var = ncm_stats_vec_get_var (svec, 0);
 *   ...
 * }
 *
 * ]|
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_vec.h"
#include "ncm/stats/ncm_stats_acorr.h"
#include "ncm/algebra/ncm_lapack.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_sort.h>

#include <math.h>

#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_LEN,
  PROP_TYPE,
  PROP_SAVE_X,
  PROP_SIZE,
};

struct _NcmStatsVec
{
  /*< private >*/
  GObject parent_instance;
  NcmStatsVecType t;
  NcmStatsVecUpdateFunc update;
  guint len;
  gboolean save_x;
  gdouble weight;
  gdouble weight2;
  gdouble bias_wt;
  guint nitens;
  NcmVector *x;
  NcmVector *mean;
  NcmVector *var;
  NcmMatrix *cov;
  NcmMatrix *real_cov;
  GPtrArray *saved_x;
  GPtrArray *q_array;
};

G_DEFINE_TYPE (NcmStatsVec, ncm_stats_vec, G_TYPE_OBJECT)

static void
ncm_stats_vec_init (NcmStatsVec *svec)
{
  svec->t       = NCM_STATS_VEC_TYPES_LEN;
  svec->update  = NULL;
  svec->weight  = 0.0;
  svec->weight2 = 0.0;
  svec->bias_wt = 0.0;

/*
 * Increment calculation
 *  svec->mean_inc = 0.0;
 *  svec->var_inc  = 0.0;
 *  svec->cov_inc  = 0.0;
 */
  svec->nitens   = 0;
  svec->x        = NULL;
  svec->mean     = NULL;
  svec->var      = NULL;
  svec->cov      = NULL;
  svec->real_cov = NULL;
  svec->saved_x  = NULL;
  svec->save_x   = FALSE;

  svec->q_array = g_ptr_array_new ();
  g_ptr_array_set_free_func (svec->q_array, (GDestroyNotify) gsl_rstat_quantile_free);
}

static void
_ncm_stats_vec_dispose (GObject *object)
{
  NcmStatsVec *svec = NCM_STATS_VEC (object);


  ncm_vector_clear (&svec->x);
  ncm_vector_clear (&svec->mean);
  ncm_vector_clear (&svec->var);
  ncm_matrix_clear (&svec->cov);
  ncm_matrix_clear (&svec->real_cov);

  if (svec->saved_x != NULL)
  {
    g_ptr_array_unref (svec->saved_x);
    svec->saved_x = NULL;
    svec->save_x  = FALSE;
  }

  g_clear_pointer (&svec->q_array, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_vec_parent_class)->dispose (object);
}

static void
_ncm_stats_vec_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_vec_parent_class)->finalize (object);
}

static void _ncm_stats_vec_update_from_vec_weight_cov (NcmStatsVec *svec, const gdouble w, NcmVector *x);
static void _ncm_stats_vec_update_from_vec_weight_var (NcmStatsVec *svec, const gdouble w, NcmVector *x);
static void _ncm_stats_vec_update_from_vec_weight_mean (NcmStatsVec *svec, const gdouble w, NcmVector *x);

static void
_ncm_stats_vec_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_stats_vec_parent_class)->constructed (object);
  {
    NcmStatsVec *svec = NCM_STATS_VEC (object);

    g_return_if_fail (NCM_IS_STATS_VEC (object));
    g_assert_cmpuint (svec->len, >, 0);

    if (svec->save_x)
    {
      g_assert (svec->saved_x == NULL);
      svec->saved_x = g_ptr_array_new ();
      g_ptr_array_set_free_func (svec->saved_x, (GDestroyNotify) ncm_vector_free);
    }

    switch (svec->t)
    {
      case NCM_STATS_VEC_COV:
        g_assert_cmpuint (svec->len, >, 0);
        g_assert (svec->cov == NULL);

        svec->cov    = ncm_matrix_new (svec->len, svec->len);
        svec->update = &_ncm_stats_vec_update_from_vec_weight_cov;

        ncm_matrix_set_zero (svec->cov);

        G_GNUC_FALLTHROUGH;
      case NCM_STATS_VEC_VAR:
        g_assert (svec->var == NULL);

        svec->var = ncm_vector_new (svec->len);

        if (svec->update == NULL)
          svec->update = &_ncm_stats_vec_update_from_vec_weight_var;

        ncm_vector_set_zero (svec->var);

        G_GNUC_FALLTHROUGH;
      case NCM_STATS_VEC_MEAN:
        g_assert (svec->x == NULL);
        g_assert (svec->mean == NULL);

        svec->x    = ncm_vector_new (svec->len);
        svec->mean = ncm_vector_new (svec->len);

        if (svec->update == NULL)
          svec->update = &_ncm_stats_vec_update_from_vec_weight_mean;

        ncm_vector_set_zero (svec->x);
        ncm_vector_set_zero (svec->mean);
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }
}

static void
_ncm_stats_vec_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsVec *svec = NCM_STATS_VEC (object);

  g_return_if_fail (NCM_IS_STATS_VEC (object));

  switch (prop_id)
  {
    case PROP_LEN:
      svec->len = g_value_get_uint (value);
      break;
    case PROP_TYPE:
      svec->t = g_value_get_enum (value);
      break;
    case PROP_SAVE_X:
      svec->save_x = g_value_get_boolean (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_vec_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsVec *svec = NCM_STATS_VEC (object);

  g_return_if_fail (NCM_IS_STATS_VEC (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, svec->len);
      break;
    case PROP_TYPE:
      g_value_set_enum (value, svec->t);
      break;
    case PROP_SAVE_X:
      g_value_set_boolean (value, svec->save_x);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_stats_vec_class_init (NcmStatsVecClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->dispose      = &_ncm_stats_vec_dispose;
  object_class->finalize     = &_ncm_stats_vec_finalize;
  object_class->constructed  = &_ncm_stats_vec_constructed;
  object_class->set_property = &_ncm_stats_vec_set_property;
  object_class->get_property = &_ncm_stats_vec_get_property;

  /**
   * NcmStatsVec:length:
   *
   * Number of random variables.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("length",
                                                      NULL,
                                                      "Statistics vector length",
                                                      1, G_MAXUINT32, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmStatsVec:type:
   *
   * The statistics to be calculated.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_TYPE,
                                   g_param_spec_enum ("type",
                                                      NULL,
                                                      "Statistics vector type",
                                                      NCM_TYPE_STATS_VEC_TYPE, NCM_STATS_VEC_MEAN,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmStatsVec:save-x:
   *
   * Whether to save each input vector.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SAVE_X,
                                   g_param_spec_boolean ("save-x",
                                                         NULL,
                                                         "Whether to save all input vectors",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_stats_vec_new:
 * @len: number of random variables
 * @t: type of statistics to be calculated
 * @save_x: whenever to save each vector x
 *
 * Creates a new #NcmStatsVec.
 *
 * Returns: (transfer full): a new #NcmStatsVec.
 */
NcmStatsVec *
ncm_stats_vec_new (guint len, NcmStatsVecType t, gboolean save_x)
{
  NcmStatsVec *svec = g_object_new (NCM_TYPE_STATS_VEC,
                                    "length", len,
                                    "type", t,
                                    "save-x", save_x,
                                    NULL);

  return svec;
}

/**
 * ncm_stats_vec_ref:
 * @svec: a #NcmStatsVec
 *
 * Increase the reference of @svec by one.
 *
 * Returns: (transfer full): @svec.
 */
NcmStatsVec *
ncm_stats_vec_ref (NcmStatsVec *svec)
{
  return g_object_ref (svec);
}

/**
 * ncm_stats_vec_free:
 * @svec: a #NcmStatsVec
 *
 * Decrease the reference count of @svec by one.
 *
 */
void
ncm_stats_vec_free (NcmStatsVec *svec)
{
  g_object_unref (svec);
}

/**
 * ncm_stats_vec_clear:
 * @svec: a #NcmStatsVec
 *
 * Decrease the reference count of @svec by one, and sets the pointer *svec to
 * NULL.
 *
 */
void
ncm_stats_vec_clear (NcmStatsVec **svec)
{
  g_clear_object (svec);
}

/**
 * ncm_stats_vec_reset:
 * @svec: a #NcmStatsVec
 * @rm_saved: a boolean
 *
 * Reset all data in @svec. If @rm_saved is TRUE and @svec has
 * saved data, it will be also removed from the object.
 *
 */
void
ncm_stats_vec_reset (NcmStatsVec *svec, gboolean rm_saved)
{
  if (rm_saved && svec->save_x)
  {
    g_assert (svec->saved_x != NULL);
    g_ptr_array_unref (svec->saved_x);
    svec->saved_x = g_ptr_array_new ();
    g_ptr_array_set_free_func (svec->saved_x, (GDestroyNotify) ncm_vector_free);
    g_ptr_array_set_size (svec->saved_x, 0);
  }

  svec->weight  = 0.0;
  svec->weight2 = 0.0;
  svec->bias_wt = 0.0;
  svec->nitens  = 0;

  switch (svec->t)
  {
    case NCM_STATS_VEC_COV:
      g_assert (svec->cov != NULL);
      ncm_matrix_set_zero (svec->cov);

      G_GNUC_FALLTHROUGH;
    case NCM_STATS_VEC_VAR:
      g_assert (svec->var != NULL);
      ncm_vector_set_zero (svec->var);

      G_GNUC_FALLTHROUGH;
    case NCM_STATS_VEC_MEAN:
      g_assert (svec->x != NULL);
      ncm_vector_set_zero (svec->x);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  if (svec->q_array->len == svec->len)
  {
    guint i;
    const gdouble p = ((gsl_rstat_quantile_workspace *) g_ptr_array_index (svec->q_array, 0))->p;

    g_ptr_array_set_size (svec->q_array, 0);

    for (i = 0; i < svec->len; i++)
    {
      gsl_rstat_quantile_workspace *qws_i = gsl_rstat_quantile_alloc (p);

      g_ptr_array_add (svec->q_array, qws_i);
    }
  }
}

static void
_ncm_stats_vec_update_from_vec_weight_cov (NcmStatsVec *svec, const gdouble w, NcmVector *x)
{
  const gdouble curweight = svec->weight + w;
  const guint sveclen     = svec->len;
  guint i;

  svec->nitens++;

  if (w == 0.0)
    return;

  for (i = 0; i < sveclen; i++)
  {
    guint j;
    gdouble mean_i        = ncm_vector_fast_get (svec->mean, i);
    const gdouble x_i     = ncm_vector_fast_get (x, i);
    const gdouble delta_i = x_i - mean_i;
    const gdouble R_i     = delta_i * w / curweight;
    const gdouble var     = ncm_vector_fast_get (svec->var, i);
    const gdouble dvar    = svec->weight * delta_i * R_i;

    mean_i += R_i;
    ncm_vector_fast_set (svec->mean, i, mean_i);
    ncm_vector_fast_set (svec->var, i, var + dvar);

    for (j = i + 1; j < sveclen; j++)
    {
      const gdouble x_j    = ncm_vector_fast_get (x, j);
      const gdouble mean_j = ncm_vector_fast_get (svec->mean, j);
      const gdouble dC_ij  = w * (x_i - mean_i) * (x_j - mean_j);
      const gdouble oC_ij  = ncm_matrix_get (svec->cov, i, j);
      const gdouble C_ij   = oC_ij + dC_ij;

      ncm_matrix_set (svec->cov, i, j, C_ij);
      ncm_matrix_set (svec->cov, j, i, C_ij);
    }
  }

  svec->weight   = curweight;
  svec->weight2 += w * w;
  svec->bias_wt  = 1.0 / (svec->weight - svec->weight2 / svec->weight);

  if (svec->q_array->len == svec->len)
  {
    guint i;

    for (i = 0; i < svec->len; i++)
    {
      const gdouble x_i                   = ncm_vector_fast_get (x, i);
      gsl_rstat_quantile_workspace *qws_i = g_ptr_array_index (svec->q_array, i);

      gsl_rstat_quantile_add (x_i, qws_i);
    }
  }
}

static void
_ncm_stats_vec_update_from_vec_weight_var (NcmStatsVec *svec, const gdouble w, NcmVector *x)
{
  const gdouble curweight = svec->weight + w;
  const guint sveclen     = svec->len;
  guint i;

  svec->nitens++;

  if (w == 0.0)
    return;

  for (i = 0; i < sveclen; i++)
  {
    const gdouble mean_i  = ncm_vector_fast_get (svec->mean, i);
    const gdouble x_i     = ncm_vector_fast_get (x, i);
    const gdouble delta_i = x_i - mean_i;
    const gdouble R_i     = delta_i * w / curweight;
    const gdouble var     = ncm_vector_fast_get (svec->var, i);
    const gdouble dvar    = svec->weight * delta_i * R_i;

    ncm_vector_fast_set (svec->mean, i, mean_i + R_i);
    ncm_vector_fast_set (svec->var, i, var + dvar);
  }

  svec->weight   = curweight;
  svec->weight2 += w * w;
  svec->bias_wt  = 1.0 / (svec->weight - svec->weight2 / svec->weight);

  if (svec->q_array->len == svec->len)
  {
    guint i;

    for (i = 0; i < svec->len; i++)
    {
      const gdouble x_i                   = ncm_vector_fast_get (x, i);
      gsl_rstat_quantile_workspace *qws_i = g_ptr_array_index (svec->q_array, i);

      gsl_rstat_quantile_add (x_i, qws_i);
    }
  }
}

static void
_ncm_stats_vec_update_from_vec_weight_mean (NcmStatsVec *svec, const gdouble w, NcmVector *x)
{
  const gdouble curweight = svec->weight + w;
  const guint sveclen     = svec->len;
  guint i;

  svec->nitens++;

  if (w == 0.0)
    return;

  for (i = 0; i < sveclen; i++)
  {
    const gdouble mean_i  = ncm_vector_fast_get (svec->mean, i);
    const gdouble x_i     = ncm_vector_fast_get (x, i);
    const gdouble delta_i = x_i - mean_i;
    const gdouble R_i     = delta_i * w / curweight;

    ncm_vector_fast_set (svec->mean, i, mean_i + R_i);
  }

  svec->weight   = curweight;
  svec->weight2 += w * w;
  svec->bias_wt  = 1.0 / (svec->weight - svec->weight2 / svec->weight);

  if (svec->q_array->len == svec->len)
  {
    guint i;

    for (i = 0; i < svec->len; i++)
    {
      const gdouble x_i                   = ncm_vector_fast_get (x, i);
      gsl_rstat_quantile_workspace *qws_i = g_ptr_array_index (svec->q_array, i);

      gsl_rstat_quantile_add (x_i, qws_i);
    }
  }
}

/**
 * ncm_stats_vec_update_weight:
 * @svec: a #NcmStatsVec
 * @w: The statistical weight
 *
 * Updates the statistics using @svec->x set in @svec and @weight, then reset
 * @svec->x to zero.
 *
 */
void
ncm_stats_vec_update_weight (NcmStatsVec *svec, const gdouble w)
{
  svec->update (svec, w, svec->x);

  if (svec->save_x)
  {
    NcmVector *v = ncm_vector_dup (svec->x);

    g_ptr_array_add (svec->saved_x, v);
  }

  ncm_vector_set_zero (svec->x);
}

/**
 * ncm_stats_vec_append_weight:
 * @svec: a #NcmStatsVec
 * @x: a #NcmVector to be added
 * @w: the weight of @x
 * @dup: a boolean
 *
 * Appends and updates the statistics using weight @w for the vector @x #NcmVector of same
 * size #NcmStatsVec:length and with continuous allocation. i.e., NcmVector:stride == 1.
 *
 * If @svec was created with save_x TRUE, the paramenter @dup determines if the vector
 * @x will be duplicated or if just a reference for @x will be saved.
 *
 */
void
ncm_stats_vec_append_weight (NcmStatsVec *svec, NcmVector *x, gdouble w, gboolean dup)
{
  svec->update (svec, w, x);

  if (svec->save_x)
  {
    if (dup)
      g_ptr_array_add (svec->saved_x, ncm_vector_dup (x));
    else
      g_ptr_array_add (svec->saved_x, ncm_vector_ref (x));
  }
}

/**
 * ncm_stats_vec_prepend_weight:
 * @svec: a #NcmStatsVec
 * @x: a #NcmVector to be added
 * @w: the weight of @x
 * @dup: a boolean
 *
 * Prepends and updates the statistics using the vector @x and weight @w.
 * It assumes that #NcmVector is of same size #NcmStatsVec:length and
 * with continuous allocation. i.e., NcmVector:stride == 1.
 *
 * If @svec was created with save_x TRUE, the paramenter @dup determines if the vector
 * will be duplicated or if just a reference for @x will be saved.
 *
 */
void
ncm_stats_vec_prepend_weight (NcmStatsVec *svec, NcmVector *x, gdouble w, gboolean dup)
{
  svec->update (svec, w, x);

  if (svec->save_x)
  {
    const guint cp_len  = 1;
    const guint old_len = svec->saved_x->len;
    const guint new_len = old_len + cp_len;

    g_ptr_array_set_size (svec->saved_x, new_len);
    memmove (&svec->saved_x->pdata[cp_len], svec->saved_x->pdata, sizeof (gpointer) * old_len);

    if (dup)
      g_ptr_array_index (svec->saved_x, 0) = ncm_vector_dup (x);
    else
      g_ptr_array_index (svec->saved_x, 0) = ncm_vector_ref (x);
  }
}

/**
 * ncm_stats_vec_append:
 * @svec: a #NcmStatsVec
 * @x: a #NcmVector to be added
 * @dup: a boolean
 *
 * Appends and updates the statistics using weight 1.0 for the vector @x #NcmVector of same
 * size #NcmStatsVec:length and with continuous allocation. i.e., NcmVector:stride == 1.
 *
 * If @svec was created with save_x TRUE, the paramenter @dup determines if the vector
 * @x will be duplicated or if just a reference for @x will be saved.
 *
 */
void
ncm_stats_vec_append (NcmStatsVec *svec, NcmVector *x, gboolean dup)
{
  ncm_stats_vec_append_weight (svec, x, 1.0, dup);
}

/**
 * ncm_stats_vec_prepend:
 * @svec: a #NcmStatsVec
 * @x: a #NcmVector to be added
 * @dup: a boolean
 *
 * Prepends and updates the statistics using the vector @x and weight 1.0.
 * It assumes that #NcmVector is of same size #NcmStatsVec:length and
 * with continuous allocation. i.e., NcmVector:stride == 1.
 *
 * If @svec was created with save_x TRUE, the paramenter @dup determines if the vector
 * will be duplicated or if just a reference for @x will be saved.
 *
 */
void
ncm_stats_vec_prepend (NcmStatsVec *svec, NcmVector *x, gboolean dup)
{
  ncm_stats_vec_prepend_weight (svec, x, 1.0, dup);
}

/**
 * ncm_stats_vec_append_data:
 * @svec: a #NcmStatsVec
 * @data: (element-type NcmVector): a #GPtrArray containing #NcmVector s to be added
 * @dup: a boolean
 *
 * Appends and updates the statistics using the data contained in @data and weight == 1.0.
 * It assumes that each element of @data is a #NcmVector of same size #NcmStatsVec:length and
 * with continuous allocation. i.e., NcmVector:stride == 1.
 *
 * If @svec was created with save_x TRUE, the paramenter @dup determines if the vectors
 * from @data will be duplicated or if just a reference for the current vectors in @data
 * will be saved.
 *
 */
void
ncm_stats_vec_append_data (NcmStatsVec *svec, GPtrArray *data, gboolean dup)
{
  guint i;

  for (i = 0; i < data->len; i++)
  {
    NcmVector *x = g_ptr_array_index (data, i);

    svec->update (svec, 1.0, x);

    if (svec->save_x)
    {
      if (dup)
        g_ptr_array_add (svec->saved_x, ncm_vector_dup (x));
      else
        g_ptr_array_add (svec->saved_x, ncm_vector_ref (x));
    }
  }
}

/**
 * ncm_stats_vec_prepend_data:
 * @svec: a #NcmStatsVec
 * @data: (element-type NcmVector): a #GPtrArray containing #NcmVector s to be added
 * @dup: a boolean
 *
 * Prepends and updates the statistics using the data contained in @data and weight == 1.0.
 * It assumes that each element of @data is a #NcmVector of same size #NcmStatsVec:length and
 * with continuous allocation. i.e., NcmVector:stride == 1.
 *
 * If @svec was created with save_x TRUE, the paramenter @dup determines if the vectors
 * from @data will be duplicated or if just a reference for the current vectors in @data
 * will be saved.
 *
 */
void
ncm_stats_vec_prepend_data (NcmStatsVec *svec, GPtrArray *data, gboolean dup)
{
  guint i;

  if (svec->save_x)
  {
    const guint cp_len  = data->len;
    const guint old_len = svec->saved_x->len;
    const guint new_len = old_len + cp_len;

    g_ptr_array_set_size (svec->saved_x, new_len);
    memmove (&svec->saved_x->pdata[cp_len], svec->saved_x->pdata, sizeof (gpointer) * old_len);
  }

  for (i = 0; i < data->len; i++)
  {
    NcmVector *x = g_ptr_array_index (data, i);

    svec->update (svec, 1.0, x);

    if (svec->save_x)
    {
      if (dup)
        g_ptr_array_index (svec->saved_x, i) = ncm_vector_dup (x);
      else
        g_ptr_array_index (svec->saved_x, i) = ncm_vector_ref (x);
    }
  }
}

/**
 * ncm_stats_vec_enable_quantile:
 * @svec: a #NcmStatsVec
 * @p: double $\in (0, 1)$
 *
 * Enables quantile calculation, it will calculate the $p$
 * quantile. Warning, it does not support weighted samples, the results
 * will ignores the weights.
 *
 */
void
ncm_stats_vec_enable_quantile (NcmStatsVec *svec, gdouble p)
{
  g_assert_cmpfloat (p, >, 0.0);
  g_assert_cmpfloat (p, <, 1.0);

  {
    guint i;

    g_ptr_array_set_size (svec->q_array, 0);

    for (i = 0; i < svec->len; i++)
    {
      gsl_rstat_quantile_workspace *qws_i = gsl_rstat_quantile_alloc (p);

      g_ptr_array_add (svec->q_array, qws_i);
    }
  }

  if (svec->nitens > 0)
  {
    if (!svec->save_x)
    {
      g_warning ("ncm_stats_vec_enable_quantile: Enabling quantile calculation in a non-empty NcmStatsVec,"
                 " all previous data will be ignored in the quantile.");
    }
    else
    {
      guint i;

      for (i = 0; i < svec->saved_x->len; i++)
      {
        NcmVector *x = g_ptr_array_index (svec->saved_x, i);
        guint j;

        for (j = 0; j < svec->len; j++)
        {
          const gdouble x_j                   = ncm_vector_fast_get (x, j);
          gsl_rstat_quantile_workspace *qws_j = g_ptr_array_index (svec->q_array, j);

          gsl_rstat_quantile_add (x_j, qws_j);
        }
      }
    }
  }
}

/**
 * ncm_stats_vec_disable_quantile:
 * @svec: a #NcmStatsVec
 *
 * Disables quantile calculation.
 *
 */
void
ncm_stats_vec_disable_quantile (NcmStatsVec *svec)
{
  g_ptr_array_set_size (svec->q_array, 0);
}

/**
 * ncm_stats_vec_get_quantile:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Returns the current quantile estimate configured by
 * ncm_stats_vec_enable_quantile().
 *
 * Returns: the current estimate of the quantile.
 */
gdouble
ncm_stats_vec_get_quantile (NcmStatsVec *svec, guint i)
{
  g_assert_cmpuint (i, <, svec->q_array->len);

  return gsl_rstat_quantile_get (g_ptr_array_index (svec->q_array, i));
}

/**
 * ncm_stats_vec_get_quantile_spread:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Returns the difference between the $(p + 1)/2$ and $p/2$ quantiles
 * configured by ncm_stats_vec_enable_quantile(). For $p = 0.5$ this is the
 * inter-quartile range.
 *
 * Returns: the current estimate of the quantile spread.
 */
gdouble
ncm_stats_vec_get_quantile_spread (NcmStatsVec *svec, guint i)
{
  g_assert_cmpuint (i, <, svec->q_array->len);

  {
    gsl_rstat_quantile_workspace *qws_i = g_ptr_array_index (svec->q_array, i);

    return qws_i->q[3] - qws_i->q[1];
  }
}

/**
 * ncm_stats_vec_get_quantile_all:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Returns the minimum, $p/2$, $p$, $(p + 1)/2$, and maximum quantiles
 * configured by ncm_stats_vec_enable_quantile().
 *
 * Returns: (transfer none) (array fixed-size=5): the current estimate of the quantile.
 */
gdouble *
ncm_stats_vec_get_quantile_all (NcmStatsVec *svec, guint i)
{
  g_assert_cmpuint (i, <, svec->q_array->len);

  {
    gsl_rstat_quantile_workspace *qws_i = g_ptr_array_index (svec->q_array, i);
    gdouble *q                          = qws_i->q;

    return q;
  }
}

static guint
_ncm_stats_vec_estimate_const_break_int (NcmStatsVec *svec, guint p, guint pad)
{
  g_assert_cmpuint (pad, <, ncm_stats_vec_nitens (svec));
  {
    const guint n                     = ncm_stats_vec_nitens (svec) - pad;
    const gsl_multifit_robust_type *T = gsl_multifit_robust_default;
    gsl_multifit_robust_workspace *w  = gsl_multifit_robust_alloc (T, n, 1);
    NcmMatrix *X                      = ncm_matrix_new (n, 1);
    NcmMatrix *cov                    = ncm_matrix_new (1, 1);
    NcmVector *y                      = ncm_vector_new (n);
    NcmVector *c                      = ncm_vector_new (1);
    gsl_multifit_robust_stats stats;
    gint status;
    gdouble t0, cutoff;
    guint i;

    for (i = 0; i < n; i++)
    {
      NcmVector *row_i = ncm_stats_vec_peek_row (svec, i + pad);

      ncm_vector_set (y, i, ncm_vector_get (row_i, p));
    }

    ncm_vector_set (c, 0, ncm_stats_vec_get_mean (svec, p));
    ncm_matrix_set_all (X, 1.0);

    gsl_multifit_robust_maxiter (100000, w);
    status = gsl_multifit_robust (ncm_matrix_gsl (X), ncm_vector_gsl (y), ncm_vector_gsl (c), ncm_matrix_gsl (cov), w);

    if ((status != GSL_SUCCESS) && (status != GSL_EMAXITER))
      g_error ("_ncm_stats_vec_estimate_const_break_int: error %d computing gsl_multifit_robust\n", status);

    stats = gsl_multifit_robust_statistics (w);

    t0 = ncm_vector_get (c, 0);

    cutoff = ceil (sqrt (gsl_cdf_chisq_Qinv (1.0 / n, 1.0)));

    /*printf ("# c0[%5u] % 22.15g % 22.15g % 22.15g\n", pad, t0, sqrt (gsl_cdf_chisq_Qinv (1.0 / n, 1.0)), sqrt (gsl_cdf_chisq_Pinv (1.0 / n, 1.0)));*/
    for (i = 0; i < n; i++)
    {
      NcmVector *row_i = ncm_stats_vec_peek_row (svec, i + pad);

/*
 *     printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g | % 22.15f\n",
 *             ncm_vector_get (row_i, p), t0, ncm_vector_get (row_i, p) - t0, sqrt (ncm_matrix_get (cov, 0, 0)),
 *             stats.sigma_ols, stats.sigma_mad, stats.sigma_rob, stats.sigma,
 *             (ncm_vector_get (row_i, p) - t0) / stats.sigma_rob
 *             );
 */
      if (fabs ((ncm_vector_get (row_i, p) - t0) / stats.sigma_rob) < cutoff)
        break;
    }

    gsl_multifit_robust_free (w);
    ncm_matrix_free (X);
    ncm_matrix_free (cov);
    ncm_vector_free (y);
    ncm_vector_free (c);

    return i;
  }
}

/**
 * ncm_stats_vec_estimate_const_break:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 *
 * Estimates the mean $\mu$ and standard deviation $\sigma$ of parameter @p
 * with robust regression and returns the first index $t_0$ within
 * $\alpha\sigma$ of $\mu$, where $\alpha$ satisfies
 * $$ \int_\alpha^\infty\chi_1(X)\mathrm{d}X = 1/N,$$
 * and $N$ is the size of the sample.
 *
 * Returns: $t_0$
 */
gdouble
ncm_stats_vec_estimate_const_break (NcmStatsVec *svec, guint p)
{
  guint n  = ncm_stats_vec_nitens (svec);
  guint t0 = 0;
  guint t1 = 0;

  do {
    t1  = _ncm_stats_vec_estimate_const_break_int (svec, p, t0);
    t0 += t1;

    if (t0 >= n)
    {
      t0 = n - 1;
      break;
    }
  } while (t1 > 0);

  return t0;
}

static gdouble
_ncm_stats_vec_heidel_diag_pcramer (const gdouble q)
{
  const gdouble pi_32 = sqrt (gsl_pow_3 (M_PI));
  const guint maxiter = 100;
  const gdouble ffac  = 1.0 / (pi_32 * sqrt (q));
  gdouble p           = 0.0;
  guint i             = 0;

  g_assert_cmpfloat (q, >=, 0.0);

  for (i = 0; i < maxiter; i++)
  {
    gint sig          = 0;
    const gdouble lnf = lgamma_r (i + 0.5, &sig) - lgamma_r (i + 1.0, &sig);
    const gdouble z   = sqrt (4.0 * i + 1.0);
    const gdouble u   = gsl_pow_2 (4.0 * i + 1.0) / (16.0 * q);
    const gdouble ti  = z * exp (-u + lnf + gsl_sf_bessel_lnKnu (0.25, u));

    p += ti;

    if (fabs (ti / p) < GSL_DBL_EPSILON)
      break;
  }

  return GSL_MIN (GSL_MAX (p * ffac, 0.0), 1.0);
}

/**
 * ncm_stats_vec_heidel_diag:
 * @svec: a #NcmStatsVec
 * @ntests: number of tests
 * @pvalue: required p-value
 * @bindex: (out): index of the best p-values
 * @wp: (out): worst parameter index
 * @wp_order: (out): worst parameter AR fit order
 * @wp_pvalue: (out): worst parameter p-value
 *
 * Applies the Heidelberger--Welch convergence diagnostic with @ntests
 * sequential Schruben tests. Uses 10 tests when @ntests is zero and a
 * p-value of $0.05$ when @pvalue is zero.
 *
 * Sets @bindex to the smallest index at which every parameter's p-value is
 * below $1 -$ @pvalue, and to -1 when no index qualifies. The returned vector
 * contains the p-values at @bindex, or for the full sample when no index
 * qualifies. @wp, @wp_order, and @wp_pvalue identify the
 * worst parameter at the selected index.
 *
 * See:
 *
 * - [Heidelberger (1981)](https://doi.org/10.1145/358598.358630)
 * - [Schruben (1982)](https://doi.org/10.1287/opre.30.3.569)
 * - [Heidelberger (1983)](https://doi.org/10.1287/opre.31.6.1109)
 *
 * Returns: (transfer full): a #NcmVector containing the best p-values.
 */
NcmVector *
ncm_stats_vec_heidel_diag (NcmStatsVec *svec, const guint ntests, const gdouble pvalue, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_pvalue)
{
  NcmStatsVec *chunk   = ncm_stats_vec_new (svec->len, NCM_STATS_VEC_VAR, TRUE);
  const gint half_size = svec->nitens / 2;
  const gint block     = (ntests == 0) ? ((half_size - 1) / 10 + 1) : ((half_size - 1) / (gint) ntests + 1);
  const gdouble onepv  = (pvalue == 0.0) ? 0.95 : (1.0 - pvalue);
  NcmVector *pvals     = ncm_vector_new (svec->len);
  NcmVector *Ivals     = ncm_vector_new (svec->len);
  NcmVector *spec0     = ncm_vector_new (svec->len);
  NcmVector *cumsum    = ncm_vector_new (svec->len);
  GArray *ar_order     = g_array_new (FALSE, FALSE, sizeof (guint));
  guint c_order        = 0;
  gint i;

  g_assert_cmpuint (svec->nitens, >=, 10);
  g_assert_cmpfloat (pvalue, <, 1.0);
  g_assert (svec->save_x);

  for (i = svec->nitens - 1; i >= half_size; i--)
  {
    NcmVector *row = ncm_stats_vec_peek_row (svec, i);

    ncm_stats_vec_append (chunk, row, FALSE);
  }

  {
    NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (svec->len,
                                                     NCM_STATS_ACORR_DEFAULT_MAX_LAG,
                                                     NCM_STATS_ACORR_DEFAULT_MAX_LEVELS,
                                                     NCM_STATS_ACORR_METHOD_AR);

    for (i = svec->nitens - 1; i >= half_size; i--)
      ncm_stats_acorr_update (acorr, ncm_stats_vec_peek_row (svec, i));

    for (i = 0; i < (gint) svec->len; i++)
    {
      c_order = ncm_stats_acorr_get_ar_order (acorr, i);

      ncm_vector_set (spec0, i, ncm_stats_acorr_get_spec0 (acorr, i));
      g_array_append_val (ar_order, c_order);
    }

    ncm_stats_acorr_free (acorr);
  }

  bindex[0] = -1;
  wp[0]     = 0;

  for (i = half_size - 1; i >= 0; i--)
  {
    NcmVector *row = ncm_stats_vec_peek_row (svec, i);

    ncm_stats_vec_append (chunk, row, FALSE);

    if ((i % block) == 0)
    {
      const guint nitens = svec->nitens - i;
      gint j;
      guint p;

      ncm_vector_set_zero (cumsum);
      ncm_vector_set_zero (Ivals);

      for (j = svec->nitens - 1; j >= i; j--)
      {
        NcmVector *row  = ncm_stats_vec_peek_row (svec, j);
        const gdouble n = (svec->nitens - j);

        for (p = 0; p < svec->len; p++)
        {
          const gdouble mean_p_n = ncm_stats_vec_get_mean (chunk, p) * n;
          const gdouble cumsum_p = ncm_vector_get (cumsum, p) + ncm_vector_get (row, p);

          ncm_vector_set (cumsum, p, cumsum_p);
          ncm_vector_addto (Ivals, p, gsl_pow_2 (cumsum_p - mean_p_n));
        }
      }

      {
        gdouble max_pval = 0.0;
        gint lwp         = 0;

        for (p = 0; p < svec->len; p++)
        {
          const gdouble Ival_p = ncm_vector_get (Ivals, p) / (gsl_pow_2 (nitens) * ncm_vector_get (spec0, p));
          const gdouble pval_p = _ncm_stats_vec_heidel_diag_pcramer (Ival_p);

          ncm_vector_set (Ivals, p, pval_p);

          if (pval_p > max_pval)
          {
            max_pval = pval_p;
            lwp      = p;
          }
        }

        if (max_pval <= onepv)
        {
          bindex[0] = i;
          wp[0]     = lwp;

          ncm_vector_memcpy (pvals, Ivals);
        }
      }
    }
  }

  if (bindex[0] == -1)
  {
    ncm_vector_memcpy (pvals, Ivals);
    wp[0] = ncm_vector_get_max_index (pvals);
  }

  wp_pvalue[0] = ncm_vector_get (pvals, wp[0]);
  wp_order[0]  = g_array_index (ar_order, guint, wp[0]);

  ncm_vector_clear (&spec0);
  ncm_vector_clear (&cumsum);
  ncm_vector_clear (&Ivals);

  ncm_stats_vec_clear (&chunk);

  g_array_unref (ar_order);

  return pvals;
}

/**
 * ncm_stats_vec_visual_heidel_diag:
 * @svec: a #NcmStatsVec
 * @p: vector index
 * @fi: first index
 * @mean: (out): mean
 * @var: (out): test's variance
 *
 * Computes the empirical cumulative, mean, and variance used by
 * ncm_stats_vec_heidel_diag().
 *
 * See ncm_stats_vec_heidel_diag().
 *
 * Returns: (transfer full): a #NcmVector containing the empirical cumulative distribution.
 */
NcmVector *
ncm_stats_vec_visual_heidel_diag (NcmStatsVec *svec, const guint p, const guint fi, gdouble *mean, gdouble *var)
{
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (1,
                                                   NCM_STATS_ACORR_DEFAULT_MAX_LAG,
                                                   NCM_STATS_ACORR_DEFAULT_MAX_LEVELS,
                                                   NCM_STATS_ACORR_METHOD_AR);
  const guint nitens  = svec->nitens - fi;
  gdouble cumsum      = 0.0;
  NcmVector *cumsum_v = ncm_vector_new (nitens);
  gint i, j = 0;

  g_assert_cmpuint (svec->nitens, >=, 10);
  g_assert_cmpuint (fi, <, svec->nitens);
  g_assert (svec->save_x);

  for (i = svec->nitens - 1; i >= (gint) fi; i--)
  {
    NcmVector *row      = ncm_stats_vec_peek_row (svec, i);
    const gdouble p_val = ncm_vector_get (row, p);

    cumsum += p_val;
    ncm_vector_set (cumsum_v, j, cumsum);
    j++;

    ncm_stats_acorr_update_var (acorr, 0, p_val);
  }

  mean[0] = ncm_stats_acorr_get_mean (acorr, 0);
  var[0]  = ncm_stats_acorr_get_spec0 (acorr, 0) * nitens;

  ncm_stats_acorr_free (acorr);

  return cumsum_v;
}

/**
 * ncm_stats_vec_max_ess_time:
 * @svec: a #NcmStatsVec
 * @ntests: number of tests
 * @bindex: (out): time index of the best ESS's
 * @wp: (out): worst parameter index
 * @wp_order: (out): worst parameter AR fit order
 * @wp_ess: (out): worst parameter ESS
 *
 * Calculates the time $t_m$ that maximizes the Effective Sample Size (ESS).
 * The variable @ntests control the number of divisions where the ESS
 * will be calculated, if it is zero the default 10 tests will be used.
 *
 * Returns: (transfer full): a #NcmVector containing the best ess.
 */
NcmVector *
ncm_stats_vec_max_ess_time (NcmStatsVec *svec, const guint ntests, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_ess)
{
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (svec->len,
                                                   NCM_STATS_ACORR_DEFAULT_MAX_LAG,
                                                   NCM_STATS_ACORR_DEFAULT_MAX_LEVELS,
                                                   NCM_STATS_ACORR_METHOD_AR);
  const gint size     = svec->nitens;
  const gint block    = (ntests == 0) ? ((size - 1) / 10 + 1) : ((size - 1) / (gint) ntests + 1);
  NcmVector *esss_tmp = ncm_vector_new (svec->len);
  NcmVector *esss     = ncm_vector_new (svec->len);
  gdouble max_t_ess   = 0.0;
  gint i, j = 0;

  g_assert_cmpuint (svec->nitens, >=, 10);
  g_assert (svec->save_x);

  bindex[0] = -1;

  for (i = size - 1; i >= 0; i--)
  {
    NcmVector *row_i = ncm_stats_vec_peek_row (svec, i);

    ncm_stats_acorr_update (acorr, row_i);

    if ((i == 0) || ((i % block == 0) && (j >= 99)))
    {
      gdouble min_ess = GSL_POSINF;
      guint cur_size  = size - i;
      guint lwp_order = 0;
      guint k, lwp    = 0;

      for (k = 0; k < svec->len; k++)
      {
        const gdouble ess   = ncm_stats_acorr_get_ess (acorr, k);
        const gdouble c_ess = GSL_MIN (cur_size, ess);

        lwp_order = ncm_stats_acorr_get_ar_order (acorr, k);

        ncm_vector_set (esss_tmp, k, ess);

        if (c_ess < min_ess)
        {
          min_ess = c_ess;
          lwp     = k;
        }
      }

      if (min_ess >= max_t_ess)
      {
        max_t_ess   = min_ess;
        bindex[0]   = i;
        wp_order[0] = lwp_order;
        wp[0]       = lwp;

        ncm_vector_memcpy (esss, esss_tmp);
      }
    }

    j++;
  }

  wp_ess[0] = ncm_vector_get (esss, wp[0]);

  ncm_vector_clear (&esss_tmp);
  ncm_stats_acorr_free (acorr);

  return esss;
}

/**
 * ncm_stats_vec_dup_saved_x:
 * @svec: a #NcmStatsVec
 *
 * Creates a copy of the internal saved_x array.
 *
 * Returns: (transfer full) (element-type NcmVector): a copy of the saved x array or NULL if it was not saved.
 */
GPtrArray *
ncm_stats_vec_dup_saved_x (NcmStatsVec *svec)
{
  if (svec->save_x)
  {
    GPtrArray *dup = NULL;
    guint i;

    g_assert (svec->saved_x != NULL);

    dup = g_ptr_array_new ();
    g_ptr_array_set_free_func (dup, (GDestroyNotify) ncm_vector_free);

    for (i = 0; i < svec->saved_x->len; i++)
    {
      NcmVector *v_i = ncm_vector_ref (g_ptr_array_index (svec->saved_x, i));

      g_ptr_array_add (dup, v_i);
    }

    return dup;
  }
  else
  {
    return NULL;
  }
}

/**
 * ncm_stats_vec_compute_cov_robust_diag:
 * @svec: a #NcmStatsVec
 *
 * Compute the covariance using the saved data applying a
 * a robust scale estimator for each degree of freedom.
 *
 *
 * Returns: (transfer full): A diagonal #NcmMatrix $D$ containing the estimated variances.
 */
NcmMatrix *
ncm_stats_vec_compute_cov_robust_diag (NcmStatsVec *svec)
{
  NcmMatrix *cov   = ncm_matrix_new (svec->len, svec->len);
  GArray *data     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *work     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *work_int = g_array_new (FALSE, FALSE, sizeof (gint));
  guint i;

  if (svec->nitens < 4)
    g_error ("ncm_stats_vec_compute_cov_robust_diag: too few points to estimate the covariance [%d].",
             svec->nitens);

  if (!svec->save_x)
    g_error ("ncm_stats_vec_compute_cov_robust_diag: This algorithm requires the saved data into the object.");

  g_array_set_size (data, svec->nitens);
  g_array_set_size (work, svec->nitens * 3);
  g_array_set_size (work_int, svec->nitens * 5);

  ncm_matrix_set_zero (cov);

  for (i = 0; i < svec->len; i++)
  {
    gdouble var_ii;
    guint a;

    for (a = 0; a < svec->nitens; a++)
    {
      NcmVector *theta_a      = ncm_stats_vec_peek_row (svec, a);
      const gdouble theta_a_i = ncm_vector_get (theta_a, i);

      g_array_index (data, gdouble, a) = theta_a_i;
    }

    gsl_sort (&g_array_index (data, gdouble, 0), 1, svec->nitens);
    var_ii = gsl_stats_Qn_from_sorted_data (
      &g_array_index (data, gdouble, 0),
      1, svec->nitens,
      &g_array_index (work, gdouble, 0),
      &g_array_index (work_int, gint, 0)
    );
    var_ii = gsl_pow_2 (var_ii);
    ncm_matrix_set (cov, i, i, var_ii);
  }

  g_array_unref (data);
  g_array_unref (work);
  g_array_unref (work_int);

  return cov;
}

/**
 * ncm_stats_vec_compute_cov_robust_ogk:
 * @svec: a #NcmStatsVec
 *
 * Compute the covariance matrix employing the Orthogonalized Gnanadesikan-Kettenring (OGK)
 * method. This method utilizes saved data and incorporates a robust scale estimator
 * for each degree of freedom. The OGK method provides a robust and efficient approach
 * to compute covariance, ensuring reliable estimates even in the presence of outliers
 * or skewed distributions.
 *
 * Returns: (transfer full): A diagonal #NcmMatrix $V$ containing the estimated covariance.
 */
NcmMatrix *
ncm_stats_vec_compute_cov_robust_ogk (NcmStatsVec *svec)
{
  NcmMatrix *cov     = ncm_matrix_new (svec->len, svec->len);
  NcmMatrix *E       = ncm_matrix_new (svec->len, svec->len);
  NcmMatrix *y       = ncm_matrix_new (svec->nitens, svec->len);
  NcmMatrix *z       = ncm_matrix_new (svec->len, svec->nitens);
  NcmVector *sigma_x = ncm_vector_new (svec->len);
  NcmVector *sigma_z = ncm_vector_new (svec->len);
  GArray *data       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *work       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *work_int   = g_array_new (FALSE, FALSE, sizeof (gint));
  guint a, i;

  if (svec->nitens < 4)
    g_error ("ncm_stats_vec_compute_cov_robust_diag: too few points to estimate the covariance [%d].",
             svec->nitens);

  if (!svec->save_x)
    g_error ("ncm_stats_vec_compute_cov_robust_diag: This algorithm requires the saved data into the object.");

  g_array_set_size (data, svec->nitens);
  g_array_set_size (work, svec->nitens * 3);
  g_array_set_size (work_int, svec->nitens * 5);

  for (a = 0; a < svec->nitens; a++)
  {
    NcmVector *theta_a = ncm_stats_vec_peek_row (svec, a);

    ncm_matrix_set_row (y, a, theta_a);
  }

  for (i = 0; i < svec->len; i++)
  {
    gdouble sigma_i;

    for (a = 0; a < svec->nitens; a++)
    {
      NcmVector *theta_a      = ncm_stats_vec_peek_row (svec, a);
      const gdouble theta_a_i = ncm_vector_get (theta_a, i);

      g_array_index (data, gdouble, a) = theta_a_i;
    }

    gsl_sort (&g_array_index (data, gdouble, 0), 1, svec->nitens);
    sigma_i = gsl_stats_Qn_from_sorted_data (
      &g_array_index (data, gdouble, 0),
      1, svec->nitens,
      &g_array_index (work, gdouble, 0),
      &g_array_index (work_int, gint, 0)
    );

    ncm_vector_set (sigma_x, i, sigma_i);
    ncm_matrix_mul_col (y, i, 1.0 / sigma_i);
  }

  ncm_matrix_set_identity (cov);

  for (i = 0; i < svec->len; i++)
  {
    guint j;

    for (j = i + 1; j < svec->len; j++)
    {
      gdouble s_ipj, s_imj;

      for (a = 0; a < svec->nitens; a++)
      {
        const gdouble y_a_i = ncm_matrix_get (y, a, i);
        const gdouble y_a_j = ncm_matrix_get (y, a, j);

        g_array_index (data, gdouble, a) = y_a_i + y_a_j;
      }

      gsl_sort (&g_array_index (data, gdouble, 0), 1, svec->nitens);
      s_ipj = gsl_stats_Qn_from_sorted_data (
        &g_array_index (data, gdouble, 0),
        1, svec->nitens,
        &g_array_index (work, gdouble, 0),
        &g_array_index (work_int, gint, 0)
      );

      for (a = 0; a < svec->nitens; a++)
      {
        const gdouble y_a_i = ncm_matrix_get (y, a, i);
        const gdouble y_a_j = ncm_matrix_get (y, a, j);

        g_array_index (data, gdouble, a) = y_a_i - y_a_j;
      }

      gsl_sort (&g_array_index (data, gdouble, 0), 1, svec->nitens);
      s_imj = gsl_stats_Qn_from_sorted_data (
        &g_array_index (data, gdouble, 0),
        1, svec->nitens,
        &g_array_index (work, gdouble, 0),
        &g_array_index (work_int, gint, 0)
      );
      ncm_matrix_set (cov, i, j, 0.25 * (s_ipj * s_ipj - s_imj * s_imj));
    }
  }

  {
    NcmLapackWS *lapack_work = ncm_lapack_ws_new ();
    gint neval;
    gint info;

    info = ncm_lapack_dsyevr ('V', 'A', 'U', svec->len,
                              ncm_matrix_data (cov), ncm_matrix_tda (cov),
                              0.0, 0.0, 0.0, 0.0,
                              0.0, &neval, &g_array_index (data, gdouble, 0),
                              ncm_matrix_data (E), ncm_matrix_tda (E),
                              &g_array_index (work_int, gint, 0),
                              lapack_work);

    NCM_LAPACK_CHECK_INFO ("dsyevr", info);

    ncm_lapack_ws_free (lapack_work);
  }

  {
    ncm_matrix_dgemm (z, 'N', 'T', 1.0, E, y, 0.0);

    for (i = 0; i < svec->len; i++)
    {
      gdouble sigma_z_i;

      for (a = 0; a < svec->nitens; a++)
      {
        const gdouble z_a_i = ncm_matrix_get (z, i, a);

        g_array_index (data, gdouble, a) = z_a_i;
      }

      gsl_sort (&g_array_index (data, gdouble, 0), 1, svec->nitens);
      sigma_z_i = gsl_stats_Qn_from_sorted_data (
        &g_array_index (data, gdouble, 0),
        1, svec->nitens,
        &g_array_index (work, gdouble, 0),
        &g_array_index (work_int, gint, 0)
      );

      ncm_vector_set (sigma_z, i, sigma_z_i);
    }
  }

  for (i = 0; i < svec->len; i++)
  {
    const gdouble sigma_z_i = ncm_vector_get (sigma_z, i);
    guint j;

    for (j = 0; j < svec->len; j++)
    {
      const gdouble sigma_x_j = ncm_vector_get (sigma_x, j);
      const gdouble E_ij      = ncm_matrix_get (E, i, j);
      const gdouble V_ij      = sigma_z_i * E_ij * sigma_x_j;

      ncm_matrix_set (E, i, j, V_ij);
    }
  }

  /* dsyrk writes one triangle; the caller gets a full covariance matrix. */
  ncm_matrix_dsyrk (cov, 'U', 'T', 1.0, E, 0.0);
  ncm_matrix_copy_triangle (cov, 'U');

  g_array_unref (data);
  g_array_unref (work);
  g_array_unref (work_int);
  ncm_matrix_free (y);
  ncm_matrix_free (z);
  ncm_matrix_free (E);
  ncm_vector_free (sigma_x);
  ncm_vector_free (sigma_z);

  return cov;
}

/**
 * ncm_stats_vec_peek_x:
 * @svec: a #NcmStatsVec
 *
 * Returns the vector containing the current value of the random variables.
 *
 * Returns: (transfer none): the random variables vector.
 */
NcmVector *
ncm_stats_vec_peek_x (NcmStatsVec *svec)
{
  return svec->x;
}

/**
 * ncm_stats_vec_set:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 * @x_i: the value of the @i-th variable
 *
 * Sets the value of the current @i-th random variable to @x_i.
 *
 */
void
ncm_stats_vec_set (NcmStatsVec *svec, guint i, gdouble x_i)
{
  ncm_vector_fast_set (svec->x, i, x_i);
}

/**
 * ncm_stats_vec_get:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Returns the value of the current @i-th random variable.
 *
 * Returns: @i-th random variable.
 */
gdouble
ncm_stats_vec_get (NcmStatsVec *svec, guint i)
{
  return ncm_vector_fast_get (svec->x, i);
}

/**
 * ncm_stats_vec_update:
 * @svec: a #NcmStatsVec.
 *
 * Same as ncm_stats_vec_update_weight() assuming weight equal to one.
 *
 */
void
ncm_stats_vec_update (NcmStatsVec *svec)
{
  ncm_stats_vec_update_weight (svec, 1.0);
}

/**
 * ncm_stats_vec_len:
 * @svec: a #NcmStatsVec.
 *
 * Gets @svec length.
 *
 * Returns: number of variables in @svec.
 */
guint
ncm_stats_vec_len (NcmStatsVec *svec)
{
  return svec->len;
}

/**
 * ncm_stats_vec_get_mean:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Return the current value of the variable mean, i.e., $\bar{x}_n$.
 *
 * Returns: $\bar{x}_n$.
 */
gdouble
ncm_stats_vec_get_mean (NcmStatsVec *svec, guint i)
{
  return ncm_vector_fast_get (svec->mean, i);
}

/**
 * ncm_stats_vec_get_var:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Return the current value of the variable variance, i.e., $Var_n$.
 *
 * Returns: $Var_n$.
 */
gdouble
ncm_stats_vec_get_var (NcmStatsVec *svec, guint i)
{
  g_assert (svec->t == NCM_STATS_VEC_VAR || svec->t == NCM_STATS_VEC_COV);

  return ncm_vector_fast_get (svec->var, i) * svec->bias_wt;
}

/**
 * ncm_stats_vec_get_sd:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Return the current value of the variable standard deviation,
 * i.e., $\sigma_n \equiv sqrt (Var_n)$.
 *
 * Returns: $\sigma_n$
 */
gdouble
ncm_stats_vec_get_sd (NcmStatsVec *svec, guint i)
{
  return sqrt (ncm_stats_vec_get_var (svec, i));
}

/**
 * ncm_stats_vec_get_cov:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 * @j: a variable index
 *
 * Return the current value of the variance between the @i-th and the @j-th
 * variables, i.e., $Cov_{ij}$.
 *
 * Returns: $Cov_{ij}$.
 */
gdouble
ncm_stats_vec_get_cov (NcmStatsVec *svec, guint i, guint j)
{
  g_assert (svec->t == NCM_STATS_VEC_COV);

  if (i == j)
    return ncm_stats_vec_get_var (svec, i);
  else
    return ncm_matrix_get (svec->cov, i, j) * svec->bias_wt;
}

/**
 * ncm_stats_vec_get_cor:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 * @j: a variable index
 *
 * Return the current value of the correlation between the @i-th and the @j-th
 * variables, i.e., $$Cor_{ij} \equiv \frac{Cov_{ij}}{\sigma_i\sigma_j}.$$
 *
 * Returns: $Cor_{ij}$.
 */
gdouble
ncm_stats_vec_get_cor (NcmStatsVec *svec, guint i, guint j)
{
  if (i == j)
    return 1.0;
  else
    return ncm_stats_vec_get_cov (svec, i, j) / (ncm_stats_vec_get_sd (svec, i) * ncm_stats_vec_get_sd (svec, j));
}

/**
 * ncm_stats_vec_get_weight:
 * @svec: a #NcmStatsVec
 *
 * Return the current value of the weight, for non-weighted means this is simply
 * the number of elements.
 *
 * Returns: $W_n$.
 */
gdouble
ncm_stats_vec_get_weight (NcmStatsVec *svec)
{
  return svec->weight;
}

/**
 * ncm_stats_vec_get_mean_vector:
 * @svec: a #NcmStatsVec
 * @mean: a #NcmVector
 * @offset: first parameter index
 *
 * Copy the current value of the means to the vector @mean starting from parameter @offset.
 *
 */
void
ncm_stats_vec_get_mean_vector (NcmStatsVec *svec, NcmVector *x, guint offset)
{
  g_assert (x != NULL);
  g_assert_cmpint (offset, <, svec->len);
  ncm_vector_memcpy2 (x, svec->mean, 0, offset, svec->len - offset);
}

/**
 * ncm_stats_vec_peek_mean:
 * @svec: a #NcmStatsVec
 *
 * Gets the local mean vector.
 *
 * Returns: (transfer none): the internal mean #NcmVector.
 */
NcmVector *
ncm_stats_vec_peek_mean (NcmStatsVec *svec)
{
  return svec->mean;
}

/**
 * ncm_stats_vec_get_cov_matrix:
 * @svec: a #NcmStatsVec
 * @m: a #NcmMatrix
 * @offset: first parameter index
 *
 * Copy the current value of the correlation between the variables to the
 * matrix @m starting from paramenter @offset.
 *
 */
void
ncm_stats_vec_get_cov_matrix (NcmStatsVec *svec, NcmMatrix *m, guint offset)
{
  guint i;

  g_assert (m != NULL);
  g_assert_cmpint (offset, <, svec->len);

  if (offset > 0)
  {
    NcmMatrix *m_src = ncm_matrix_get_submatrix (svec->cov, offset, offset, svec->len - offset, svec->len - offset);

    ncm_matrix_memcpy (m, m_src);
    ncm_matrix_free (m_src);
  }
  else
  {
    ncm_matrix_memcpy (m, svec->cov);
  }

  for (i = 0; i < svec->len - offset; i++)
    ncm_matrix_set (m, i, i, ncm_vector_fast_get (svec->var, i + offset));

  ncm_matrix_scale (m, svec->bias_wt);
}

/**
 * ncm_stats_vec_peek_cov_matrix:
 * @svec: a #NcmStatsVec
 * @offset: first parameter index
 *
 * Gets the internal covariance matrix starting from paramenter @offset.
 * This is the internal matrix of @svec and can change with further
 * additions to @svec. It is not guaranteed to be valid after new additions.
 *
 * Returns: (transfer none): the covariance matrix.
 */
NcmMatrix *
ncm_stats_vec_peek_cov_matrix (NcmStatsVec *svec, guint offset)
{
  gint effsize = svec->len - offset;

  g_assert_cmpint (effsize, >, 0);

  if (svec->real_cov != NULL)
  {
    if ((gint) ncm_matrix_nrows (svec->real_cov) != effsize)
    {
      ncm_matrix_free (svec->real_cov);
      svec->real_cov = ncm_matrix_new (effsize, effsize);
    }
  }
  else
  {
    svec->real_cov = ncm_matrix_new (effsize, effsize);
  }

  ncm_stats_vec_get_cov_matrix (svec, svec->real_cov, offset);

  return svec->real_cov;
}

/**
 * ncm_stats_vec_nrows:
 * @svec: a #NcmStatsVec
 *
 * Gets the number of saved rows, this function fails if the object
 * was not created with save_x == TRUE;
 *
 * Returns: the number of saved rows.
 */
guint
ncm_stats_vec_nrows (NcmStatsVec *svec)
{
  g_assert (svec->save_x);

  return svec->saved_x->len;
}

/**
 * ncm_stats_vec_nitens:
 * @svec: a #NcmStatsVec
 *
 * Gets the number of items added to the object;
 *
 * Returns: the number of items added.
 */
guint
ncm_stats_vec_nitens (NcmStatsVec *svec)
{
  return svec->nitens;
}

/**
 * ncm_stats_vec_peek_row:
 * @svec: a #NcmStatsVec
 * @i: the row's index
 *
 * The i-th data row used in the statistics, this function fails if the object
 * was not created with save_x == TRUE;
 *
 * Returns: (transfer none): the i-th data row.
 */
NcmVector *
ncm_stats_vec_peek_row (NcmStatsVec *svec, guint i)
{
  g_assert (svec->save_x);
  g_assert (i < svec->saved_x->len);

  return g_ptr_array_index (svec->saved_x, i);
}

/**
 * ncm_stats_vec_get_param_at:
 * @svec: a #NcmStatsVec
 * @i: the row's index
 * @p: the parameter's index
 *
 * Gets the p-th parameter in the i-th data row used in the statistics, this
 * function fails if the object was not created with save_x == TRUE;
 *
 * Returns: the parameter value.
 */
gdouble
ncm_stats_vec_get_param_at (NcmStatsVec *svec, guint i, guint p)
{
  g_assert (svec->save_x);
  g_assert (i < svec->nitens);

  return ncm_vector_get (g_ptr_array_index (svec->saved_x, i), p);
}

