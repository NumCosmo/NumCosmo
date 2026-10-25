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
#include "ncm/stats/ncm_stats_rows.h"
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
  GArray *saved_w;
  NcmVector *x_contiguous;
  GPtrArray *q_array;
  NcmStatsRows *rows;
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

  svec->nitens   = 0;
  svec->x        = NULL;
  svec->mean     = NULL;
  svec->var      = NULL;
  svec->cov      = NULL;
  svec->real_cov = NULL;
  svec->saved_x  = NULL;
  svec->saved_w  = NULL;
  svec->save_x   = FALSE;
  svec->rows     = NULL;

  svec->x_contiguous = NULL;

  svec->q_array = g_ptr_array_new ();
  g_ptr_array_set_free_func (svec->q_array, (GDestroyNotify) gsl_rstat_quantile_free);
}

static void
_ncm_stats_vec_dispose (GObject *object)
{
  NcmStatsVec *svec = NCM_STATS_VEC (object);

  ncm_stats_rows_clear (&svec->rows);

  ncm_vector_clear (&svec->x);
  ncm_vector_clear (&svec->x_contiguous);
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

  g_clear_pointer (&svec->saved_w, g_array_unref);

  g_clear_pointer (&svec->q_array, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_vec_parent_class)->dispose (object);
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
      svec->saved_w = g_array_new (FALSE, FALSE, sizeof (gdouble));
      g_ptr_array_set_free_func (svec->saved_x, (GDestroyNotify) ncm_vector_free);
    }

    svec->x_contiguous = ncm_vector_new (svec->len);

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
 * @t: statistics to compute
 * @save_x: whether to save each added vector
 *
 * Creates a new #NcmStatsVec, see #NcmStatsVec:length, #NcmStatsVec:type and
 * #NcmStatsVec:save-x.
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
 * Increases the reference count of @svec by one.
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
 * Decreases the reference count of @svec by one.
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
 * Decreases the reference count of *@svec by one and sets *@svec to %NULL.
 */
void
ncm_stats_vec_clear (NcmStatsVec **svec)
{
  g_clear_object (svec);
}

/**
 * ncm_stats_vec_reset:
 * @svec: a #NcmStatsVec
 * @rm_saved: whether to remove the saved rows
 *
 * Restarts the statistics, the quantiles and ncm_stats_vec_nitens(). With #NcmStatsVec:save-x
 * and @rm_saved %TRUE the saved rows are removed too. With @rm_saved %FALSE they are kept and
 * later rows are saved after them, at their absolute positions, so the rows
 * $[0, \mathrm{nitens})$, which ncm_stats_vec_peek_row(), ncm_stats_vec_get_param_at() and the
 * row-based diagnostics read, are then not the rows in the statistics.
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
    g_array_set_size (svec->saved_w, 0);
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
      ncm_vector_set_zero (svec->mean);
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

/* The input as a contiguous vector: @x itself, or its copy in x_contiguous when it is strided */
static NcmVector *
_ncm_stats_vec_contiguous (NcmStatsVec *svec, NcmVector *x)
{
  if (ncm_vector_len (x) != svec->len)
    g_error ("ncm_stats_vec: vector of length %u added to a NcmStatsVec of length %u.", ncm_vector_len (x), svec->len);

  if (ncm_vector_stride (x) == 1)
    return x;

  ncm_vector_memcpy (svec->x_contiguous, x);

  return svec->x_contiguous;
}

static void
_ncm_stats_vec_save (NcmStatsVec *svec, NcmVector *x, const gdouble w, gboolean dup, const guint pos)
{
  NcmVector *v = dup ? ncm_vector_dup (x) : ncm_vector_ref (x);

  g_ptr_array_insert (svec->saved_x, (pos == G_MAXUINT) ? -1 : (gint) pos, v);

  if (pos == G_MAXUINT)
    g_array_append_val (svec->saved_w, w);
  else
    g_array_insert_val (svec->saved_w, pos, w);
}

/**
 * ncm_stats_vec_update_weight:
 * @svec: a #NcmStatsVec
 * @w: weight, non-negative
 *
 * Adds the current vector, set by ncm_stats_vec_set(), with weight @w, then sets the
 * current vector to zero. A zero weight counts in ncm_stats_vec_nitens() but changes no
 * statistic.
 */
void
ncm_stats_vec_update_weight (NcmStatsVec *svec, const gdouble w)
{
  svec->update (svec, w, svec->x);

  if (svec->save_x)
    _ncm_stats_vec_save (svec, svec->x, w, TRUE, G_MAXUINT);

  ncm_vector_set_zero (svec->x);
}

/**
 * ncm_stats_vec_append_weight:
 * @svec: a #NcmStatsVec
 * @x: vector of length #NcmStatsVec:length, of any stride
 * @w: weight, non-negative
 * @dup: whether to save a copy of @x instead of a reference
 *
 * Adds @x with weight @w and, with #NcmStatsVec:save-x, saves it as the last row. Aborts if
 * the length of @x differs from #NcmStatsVec:length.
 */
void
ncm_stats_vec_append_weight (NcmStatsVec *svec, NcmVector *x, gdouble w, gboolean dup)
{
  svec->update (svec, w, _ncm_stats_vec_contiguous (svec, x));

  if (svec->save_x)
    _ncm_stats_vec_save (svec, x, w, dup, G_MAXUINT);
}

/**
 * ncm_stats_vec_prepend_weight:
 * @svec: a #NcmStatsVec
 * @x: vector of length #NcmStatsVec:length, of any stride
 * @w: weight, non-negative
 * @dup: whether to save a copy of @x instead of a reference
 *
 * Adds @x with weight @w and, with #NcmStatsVec:save-x, saves it as the first row. Aborts
 * if the length of @x differs from #NcmStatsVec:length.
 */
void
ncm_stats_vec_prepend_weight (NcmStatsVec *svec, NcmVector *x, gdouble w, gboolean dup)
{
  svec->update (svec, w, _ncm_stats_vec_contiguous (svec, x));

  if (svec->save_x)
    _ncm_stats_vec_save (svec, x, w, dup, 0);
}

/**
 * ncm_stats_vec_append:
 * @svec: a #NcmStatsVec
 * @x: vector of length #NcmStatsVec:length, of any stride
 * @dup: whether to save a copy of @x instead of a reference
 *
 * ncm_stats_vec_append_weight() with weight 1.
 */
void
ncm_stats_vec_append (NcmStatsVec *svec, NcmVector *x, gboolean dup)
{
  ncm_stats_vec_append_weight (svec, x, 1.0, dup);
}

/**
 * ncm_stats_vec_prepend:
 * @svec: a #NcmStatsVec
 * @x: vector of length #NcmStatsVec:length, of any stride
 * @dup: whether to save a copy of @x instead of a reference
 *
 * ncm_stats_vec_prepend_weight() with weight 1.
 */
void
ncm_stats_vec_prepend (NcmStatsVec *svec, NcmVector *x, gboolean dup)
{
  ncm_stats_vec_prepend_weight (svec, x, 1.0, dup);
}

/**
 * ncm_stats_vec_append_data:
 * @svec: a #NcmStatsVec
 * @data: (element-type NcmVector): vectors of length #NcmStatsVec:length
 * @dup: whether to save copies of the vectors instead of references
 *
 * Calls ncm_stats_vec_append() on each element of @data, in order.
 */
void
ncm_stats_vec_append_data (NcmStatsVec *svec, GPtrArray *data, gboolean dup)
{
  guint i;

  for (i = 0; i < data->len; i++)
  {
    NcmVector *x = g_ptr_array_index (data, i);

    svec->update (svec, 1.0, _ncm_stats_vec_contiguous (svec, x));

    if (svec->save_x)
      _ncm_stats_vec_save (svec, x, 1.0, dup, G_MAXUINT);
  }
}

/**
 * ncm_stats_vec_prepend_data:
 * @svec: a #NcmStatsVec
 * @data: (element-type NcmVector): vectors of length #NcmStatsVec:length
 * @dup: whether to save copies of the vectors instead of references
 *
 * Adds the elements of @data with weight 1 and, with #NcmStatsVec:save-x, saves them before
 * the existing rows in the order of @data.
 */
void
ncm_stats_vec_prepend_data (NcmStatsVec *svec, GPtrArray *data, gboolean dup)
{
  guint i;

  for (i = 0; i < data->len; i++)
  {
    NcmVector *x = g_ptr_array_index (data, i);

    svec->update (svec, 1.0, _ncm_stats_vec_contiguous (svec, x));

    if (svec->save_x)
      _ncm_stats_vec_save (svec, x, 1.0, dup, i);
  }
}

/**
 * ncm_stats_vec_enable_quantile:
 * @svec: a #NcmStatsVec
 * @p: probability, in $(0, 1)$
 *
 * Enables the running estimate of the $p$ quantile of each variable, together with the
 * $p/2$ and $(1 + p)/2$ quantiles, by the P-squared algorithm of GSL. The quantiles ignore
 * the weights, except that rows of zero weight are left out. On a non-empty @svec the
 * saved rows are replayed; without #NcmStatsVec:save-x the earlier rows are left out, with
 * a warning.
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

        /* As in the updates, a row of zero weight does not enter the quantiles */
        if (g_array_index (svec->saved_w, gdouble, i) == 0.0)
          continue;

        for (j = 0; j < svec->len; j++)
        {
          const gdouble x_j                   = ncm_vector_get (x, j);
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
 * Disables the quantile estimates.
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

/* The NcmStatsRows of the deprecated row analyses, created on first use and given the saved
 * rows and their number; each analysis resets it after use. */
static NcmStatsRows *
_ncm_stats_vec_peek_rows (NcmStatsVec *svec)
{
  g_assert (svec->save_x);

  if (svec->rows == NULL)
    svec->rows = ncm_stats_rows_new (svec->len);

  ncm_stats_rows_set_rows (svec->rows, svec->saved_x);
  ncm_stats_rows_set_nrows (svec->rows, svec->nitens);

  return svec->rows;
}

/**
 * ncm_stats_vec_estimate_const_break:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 *
 * Calls ncm_stats_rows_estimate_const_break() on the saved rows. Requires
 * #NcmStatsVec:save-x. Deprecated with #NcmStatsVec:save-x, dropped in 1.0: read the rows
 * from where they live and call ncm_stats_rows_estimate_const_break().
 *
 * Returns: $t_0$
 */
gdouble
ncm_stats_vec_estimate_const_break (NcmStatsVec *svec, guint p)
{
  NcmStatsRows *rows = _ncm_stats_vec_peek_rows (svec);
  gdouble res        = ncm_stats_rows_estimate_const_break (rows, p);

  ncm_stats_rows_reset (rows);

  return res;
}

/**
 * ncm_stats_vec_heidel_diag:
 * @svec: a #NcmStatsVec
 * @ntests: number of tests
 * @pvalue: required p-value
 * @bindex: (out): index of the best p-values
 * @wp: (out): worst parameter index
 * @wp_order: (out): worst parameter AR fit order
 * @wp_pvalue: (out): value of the worst parameter
 *
 * Calls ncm_stats_rows_heidel_diag() on the saved rows. Requires #NcmStatsVec:save-x and at
 * least 10 rows. Deprecated with #NcmStatsVec:save-x, dropped in 1.0: read the rows from
 * where they live and call ncm_stats_rows_heidel_diag().
 *
 * Returns: (transfer full): the Cramér-von Mises cumulative distribution values.
 */
NcmVector *
ncm_stats_vec_heidel_diag (NcmStatsVec *svec, const guint ntests, const gdouble pvalue, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_pvalue)
{
  NcmStatsRows *rows = _ncm_stats_vec_peek_rows (svec);
  NcmVector *res     = ncm_stats_rows_heidel_diag (rows, ntests, pvalue, bindex, wp, wp_order, wp_pvalue);

  ncm_stats_rows_reset (rows);

  return res;
}

/**
 * ncm_stats_vec_visual_heidel_diag:
 * @svec: a #NcmStatsVec
 * @p: vector index
 * @fi: first index
 * @mean: (out): mean
 * @var: (out): test's variance
 *
 * Calls ncm_stats_rows_visual_heidel_diag() on the saved rows. Requires #NcmStatsVec:save-x
 * and at least 10 rows. Deprecated with #NcmStatsVec:save-x, dropped in 1.0: read the rows
 * from where they live and call ncm_stats_rows_visual_heidel_diag().
 *
 * Returns: (transfer full): the cumulative sums.
 */
NcmVector *
ncm_stats_vec_visual_heidel_diag (NcmStatsVec *svec, const guint p, const guint fi, gdouble *mean, gdouble *var)
{
  NcmStatsRows *rows = _ncm_stats_vec_peek_rows (svec);
  NcmVector *res     = ncm_stats_rows_visual_heidel_diag (rows, p, fi, mean, var);

  ncm_stats_rows_reset (rows);

  return res;
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
 * Calls ncm_stats_rows_max_ess_time() on the saved rows. Requires #NcmStatsVec:save-x and
 * at least 10 rows. Deprecated with #NcmStatsVec:save-x, dropped in 1.0: read the rows from
 * where they live and call ncm_stats_rows_max_ess_time().
 *
 * Returns: (transfer full): a #NcmVector containing the best ess.
 */
NcmVector *
ncm_stats_vec_max_ess_time (NcmStatsVec *svec, const guint ntests, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_ess)
{
  NcmStatsRows *rows = _ncm_stats_vec_peek_rows (svec);
  NcmVector *res     = ncm_stats_rows_max_ess_time (rows, ntests, bindex, wp, wp_order, wp_ess);

  ncm_stats_rows_reset (rows);

  return res;
}

/**
 * ncm_stats_vec_dup_saved_x:
 * @svec: a #NcmStatsVec
 *
 * Creates a new array with references to the saved rows.
 *
 * Returns: (transfer full) (element-type NcmVector) (nullable): the saved rows, or %NULL without #NcmStatsVec:save-x.
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
 * Calls ncm_stats_rows_cov_robust_diag() on the saved rows. Requires #NcmStatsVec:save-x and
 * at least 4 rows. Deprecated with #NcmStatsVec:save-x, dropped in 1.0: read the rows from
 * where they live and call ncm_stats_rows_cov_robust_diag().
 *
 * Returns: (transfer full): a diagonal matrix of the variances.
 */
NcmMatrix *
ncm_stats_vec_compute_cov_robust_diag (NcmStatsVec *svec)
{
  if (!svec->save_x)
    g_error ("ncm_stats_vec_compute_cov_robust_diag: This algorithm requires the saved data into the object.");

  NcmStatsRows *rows = _ncm_stats_vec_peek_rows (svec);
  NcmMatrix *res     = ncm_stats_rows_cov_robust_diag (rows);

  ncm_stats_rows_reset (rows);

  return res;
}

/**
 * ncm_stats_vec_compute_cov_robust_ogk:
 * @svec: a #NcmStatsVec
 *
 * Calls ncm_stats_rows_cov_robust_ogk() on the saved rows. Requires #NcmStatsVec:save-x and
 * at least 4 rows. Deprecated with #NcmStatsVec:save-x, dropped in 1.0: read the rows from
 * where they live and call ncm_stats_rows_cov_robust_ogk().
 *
 * Returns: (transfer full): the covariance matrix.
 */
NcmMatrix *
ncm_stats_vec_compute_cov_robust_ogk (NcmStatsVec *svec)
{
  if (!svec->save_x)
    g_error ("ncm_stats_vec_compute_cov_robust_ogk: This algorithm requires the saved data into the object.");

  NcmStatsRows *rows = _ncm_stats_vec_peek_rows (svec);
  NcmMatrix *res     = ncm_stats_rows_cov_robust_ogk (rows);

  ncm_stats_rows_reset (rows);

  return res;
}

/**
 * ncm_stats_vec_peek_x:
 * @svec: a #NcmStatsVec
 *
 * Gets the current vector, the one ncm_stats_vec_update() adds.
 *
 * Returns: (transfer none): the current vector.
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
 * Sets the @i-th element of the current vector to @x_i.
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
 * Returns: the @i-th element of the current vector.
 */
gdouble
ncm_stats_vec_get (NcmStatsVec *svec, guint i)
{
  return ncm_vector_fast_get (svec->x, i);
}

/**
 * ncm_stats_vec_update:
 * @svec: a #NcmStatsVec
 *
 * ncm_stats_vec_update_weight() with weight 1.
 */
void
ncm_stats_vec_update (NcmStatsVec *svec)
{
  ncm_stats_vec_update_weight (svec, 1.0);
}

/**
 * ncm_stats_vec_len:
 * @svec: a #NcmStatsVec
 *
 * Returns: the number of variables.
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
 * Returns: the weighted mean $\bar{x}_n$ of the @i-th variable.
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
 * Gets the bias-corrected weighted variance $V_n$ of the @i-th variable; NaN with a single
 * row. Requires #NCM_STATS_VEC_VAR or #NCM_STATS_VEC_COV.
 *
 * Returns: $V_n$.
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
 * Returns: the standard deviation $\sqrt{V_n}$ of the @i-th variable.
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
 * Gets the bias-corrected weighted covariance of the @i-th and @j-th variables. Requires
 * #NCM_STATS_VEC_COV.
 *
 * Returns: $\mathrm{Cov}_{ij}$.
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
 * Returns: the correlation $\mathrm{Cov}_{ij}/(\sigma_i\sigma_j)$ of the @i-th and @j-th variables.
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
 * Returns: the total weight $W_n$; the number of rows when every weight is 1.
 */
gdouble
ncm_stats_vec_get_weight (NcmStatsVec *svec)
{
  return svec->weight;
}

/**
 * ncm_stats_vec_get_mean_vector:
 * @svec: a #NcmStatsVec
 * @mean: vector of length at least #NcmStatsVec:length minus @offset
 * @offset: first variable index
 *
 * Copies the means of the variables from @offset on to @mean.
 */
void
ncm_stats_vec_get_mean_vector (NcmStatsVec *svec, NcmVector *mean, guint offset)
{
  g_assert (mean != NULL);
  g_assert_cmpint (offset, <, svec->len);
  ncm_vector_memcpy2 (mean, svec->mean, 0, offset, svec->len - offset);
}

/**
 * ncm_stats_vec_peek_mean:
 * @svec: a #NcmStatsVec
 *
 * Returns: (transfer none): the vector of the means, updated in place by further rows.
 */
NcmVector *
ncm_stats_vec_peek_mean (NcmStatsVec *svec)
{
  return svec->mean;
}

/**
 * ncm_stats_vec_get_cov_matrix:
 * @svec: a #NcmStatsVec
 * @m: square matrix of order #NcmStatsVec:length minus @offset
 * @offset: first variable index
 *
 * Copies the covariance of the variables from @offset on to @m. Aborts unless @svec
 * was created with #NCM_STATS_VEC_COV.
 */
void
ncm_stats_vec_get_cov_matrix (NcmStatsVec *svec, NcmMatrix *m, guint offset)
{
  guint i;

  g_assert (m != NULL);
  g_assert_cmpint (offset, <, svec->len);

  if (svec->t != NCM_STATS_VEC_COV)
    g_error ("ncm_stats_vec_get_cov_matrix: the NcmStatsVec does not compute the covariance, create it with NCM_STATS_VEC_COV.");

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
 * Fills an internal matrix with ncm_stats_vec_get_cov_matrix() and returns it; the matrix
 * is overwritten by the next call and is not updated by further rows.
 *
 * Returns: (transfer none): the covariance matrix of the variables from @offset on.
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
 * Gets the number of saved rows, including rows kept by ncm_stats_vec_reset() with
 * @rm_saved %FALSE. Requires #NcmStatsVec:save-x.
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
 * Returns: the number of rows added since the last reset, zero-weight rows included.
 */
guint
ncm_stats_vec_nitens (NcmStatsVec *svec)
{
  return svec->nitens;
}

/**
 * ncm_stats_vec_peek_row:
 * @svec: a #NcmStatsVec
 * @i: position of the saved row
 *
 * Gets the saved row at position @i, counting every saved row; see ncm_stats_vec_reset() for
 * rows kept across a reset. Requires #NcmStatsVec:save-x.
 *
 * Returns: (transfer none): the saved row at position @i.
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
 * @i: position of the saved row
 * @p: element index
 *
 * Gets element @p of the saved row at position @i, see ncm_stats_vec_peek_row(); @i must be
 * below ncm_stats_vec_nitens(). Requires #NcmStatsVec:save-x.
 *
 * Returns: the element.
 */
gdouble
ncm_stats_vec_get_param_at (NcmStatsVec *svec, guint i, guint p)
{
  g_assert (svec->save_x);
  g_assert (i < svec->nitens);

  return ncm_vector_get (g_ptr_array_index (svec->saved_x, i), p);
}

