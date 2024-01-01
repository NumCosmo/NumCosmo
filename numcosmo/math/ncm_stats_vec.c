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
 * SECTION:ncm_stats_vec
 * @title: NcmStatsVec
 * @short_description: An online statistics vector.
 *
 * This object calculates some basic statistics (mean, variance and covariance)
 * of a set of random variables.
 *
 * The mean can be calculated online using the following formula:
 * $$\bar{x}_n = \bar{x}_{n-1} + (x_n - \bar{x}_{n-1})\frac{w_n}{W_n},$$
 * where $\bar{x}_n$ is the mean calculated using the first $n$ elements,
 * $x_n$ is the $n$-th element, $w_n$ the $n$-th weight and finally $W_n$
 * is the sum of the first $n$ weights.
 *
 * Using the expressions above we obtain the variance from as following:
 * $$M_n = M_{n-1} + (x_n - \bar{x}_{n-1})^2w_n\frac{W_{n-1}}{W_n},$$
 * where the variance of the first $n$ elements is
 * $$V_n = \frac{M_n}{W^\text{bias}_{n}}, \quad W^\text{bias}_{n} \equiv \frac{W_n^2 - \sum^n_iw_i^2}{W_n}.$$
 * In the formula above we defined the bias corrected weight $W^\text{bias}_{n}$.
 *
 * Finally, the covariance is computed through the following expression:
 * $$N(x,y)_n = N(x,y)_{n-1} + (x_n - \bar{x}_n)(y_n - \bar{y}_{n-1})w_n,$$
 * where the covariance of two variables $x$, $y$ is given by
 * $$Cov(x,y)_n = \frac{N(x,y)_n}{W^\text{bias}_{n}}.$$
 *
 * # Using a NcmStatsVec. #
 * |[<!-- language="C" -->
 *
 * // Creates a new one dimensional NcmStatsVec to calculates mean and variance.
 * NcmStatsVec *svec = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
 *
 * // Set and update three different values of the only random variable.
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
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_vec.h"
#include "math/ncm_lapack.h"
#include "math/ncm_cfg.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#ifdef HAVE_FFTW3
#include <fftw3.h>
#endif /* HAVE_FFTW3 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_sort.h>

#include <math.h>

#include "toeplitz/solvers/toeplitz.h"
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
  NcmStatsVec *tmp;
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

#ifdef HAVE_FFTW3
  guint fft_size;
  guint fft_plan_size;
  gdouble *param_data;
  fftw_complex *param_fft;
  fftw_plan param_r2c;
  fftw_plan param_c2r;
#endif /* HAVE_FFTW3 */
};

G_DEFINE_TYPE (NcmStatsVec, ncm_stats_vec, G_TYPE_OBJECT)

static void
ncm_stats_vec_init (NcmStatsVec *svec)
{
  svec->t       = NCM_STATS_VEC_TYPES_LEN;
  svec->update  = NULL;
  svec->tmp     = NULL;
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

#ifdef HAVE_FFTW3
  svec->fft_size      = 0;
  svec->fft_plan_size = 0;
  svec->param_data    = NULL;
  svec->param_fft     = NULL;
  svec->param_r2c     = NULL;
  svec->param_c2r     = NULL;
#endif /* HAVE_FFTW3 */
}

static void
_ncm_stats_vec_dispose (GObject *object)
{
  NcmStatsVec *svec = NCM_STATS_VEC (object);

  ncm_stats_vec_clear (&svec->tmp);

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
  NcmStatsVec *svec = NCM_STATS_VEC (object);

#ifdef HAVE_FFTW3
  g_clear_pointer (&svec->param_fft,  fftw_free);
  g_clear_pointer (&svec->param_data, fftw_free);
  g_clear_pointer (&svec->param_c2r,  fftw_destroy_plan);
  g_clear_pointer (&svec->param_r2c,  fftw_destroy_plan);
#endif /* HAVE_FFTW3 */

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
        g_assert_cmpuint (svec->len, >, 1);
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
   * Whenever to save each vector x through each interation.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SAVE_X,
                                   g_param_spec_boolean ("save-x",
                                                         NULL,
                                                         "Whenever to save all x vectors",
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
 * Returns the current estimate of the quantile initialized
 * through ncm_stats_vec_enable_quantile().
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
 * Returns the current estimate of the quantile spread, from the
 * probability $p$ initialized through ncm_stats_vec_enable_quantile(),
 * i.e., it returns the difference between $(p + 1)/2$ quantile
 * and the $p/2$. For example, if $p = 0.5$ then it returns the
 * interquartile range.
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

static void
_ncm_stats_vec_get_autocorr_alloc (NcmStatsVec *svec, guint size)
{
#ifdef HAVE_FFTW3
  const guint effsize = ncm_util_fact_size (2 * size);

  if (svec->tmp == NULL)
    svec->tmp = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);

  if (!svec->save_x)
    g_error ("_ncm_stats_vec_get_autocorr_alloc: NcmStatsVec must have saved data to calculate autocorrelation.");

  if (svec->fft_size < effsize)
  {
    g_clear_pointer (&svec->param_fft,  fftw_free);
    g_clear_pointer (&svec->param_data, fftw_free);
    {
      svec->param_fft  = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * (effsize / 2 + 1));
      svec->param_data = (gdouble *) fftw_malloc (sizeof (gdouble) * effsize);
    }
    svec->fft_size = effsize;
  }

  if (svec->fft_plan_size != effsize)
  {
    g_clear_pointer (&svec->param_c2r, fftw_destroy_plan);
    g_clear_pointer (&svec->param_r2c, fftw_destroy_plan);

    /*g_debug ("# _ncm_stats_vec_get_autocorr_alloc: calculating wisdown %u\n", effsize);*/
    ncm_cfg_load_fftw_wisdom ("ncm_stats_vec_autocorr_%u", effsize);

    ncm_cfg_lock_plan_fftw ();
    svec->param_r2c = fftw_plan_dft_r2c_1d (effsize, svec->param_data, svec->param_fft, /*fftw_default_flags*/ FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    svec->param_c2r = fftw_plan_dft_c2r_1d (effsize, svec->param_fft, svec->param_data, /*fftw_default_flags*/ FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    ncm_cfg_unlock_plan_fftw ();

    ncm_cfg_save_fftw_wisdom ("ncm_stats_vec_autocorr_%u", effsize);
    svec->fft_plan_size = effsize;
    /*g_debug ("# _ncm_stats_vec_get_autocorr_alloc: calculated  wisdown %u\n", effsize);*/
  }

#endif /* HAVE_FFTW3 */
}

static void
_ncm_stats_vec_get_autocov (NcmStatsVec *svec, guint p, guint subsample, guint pad)
{
#ifdef HAVE_FFTW3
  guint eff_nitens = svec->nitens / subsample - pad;

  g_assert_cmpuint (svec->nitens / subsample, >, pad);

  if (eff_nitens == 0)
    g_error ("_ncm_stats_vec_get_autocov: too few itens to calculate.");

  _ncm_stats_vec_get_autocorr_alloc (svec, eff_nitens);
  {
    guint i;
    const gdouble mean = ncm_stats_vec_get_mean (svec, p);

    memset (&svec->param_data[eff_nitens], 0, sizeof (gdouble) * (svec->fft_plan_size - eff_nitens));

    if (subsample > 1)
    {
      ncm_stats_vec_reset (svec->tmp, TRUE);

      for (i = 0; i < eff_nitens; i++)
      {
        guint j;
        gdouble e_mean = 0.0;

        for (j = 0; j < subsample; j++)
        {
          e_mean += ncm_vector_get (g_ptr_array_index (svec->saved_x, (i + pad) * subsample + j), p);
        }

        e_mean = e_mean / (1.0 * subsample);

        svec->param_data[i] = (e_mean - mean);

        ncm_stats_vec_set (svec->tmp, 0, svec->param_data[i]);
        ncm_stats_vec_update (svec->tmp);
      }
    }
    else
    {
      for (i = 0; i < eff_nitens; i++)
      {
        svec->param_data[i] = (ncm_vector_get (g_ptr_array_index (svec->saved_x, i + pad), p) - mean);
      }
    }

    fftw_execute (svec->param_r2c);

    for (i = 0; i < svec->fft_plan_size / 2 + 1; i++)
    {
      svec->param_fft[i] = svec->param_fft[i] * conj (svec->param_fft[i]);
    }

    fftw_execute (svec->param_c2r);
  }
#endif /* HAVE_FFTW3 */
}

/**
 * ncm_stats_vec_get_autocorr:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 *
 * Calculates the autocorrelation vector, the j-th element represent
 * the selfcorrelation with lag-j.
 *
 * The returning vector use the internal memory allocation and will
 * change with subsequent calls to ncm_stats_vec_get_autocorr().
 *
 * Returns: (transfer full): the autocorrelation vector.
 */
NcmVector *
ncm_stats_vec_get_autocorr (NcmStatsVec *svec, guint p)
{
#ifdef HAVE_FFTW3
  _ncm_stats_vec_get_autocov (svec, p, 1, 0);
  {
    NcmVector *autocor = ncm_vector_new_data_dup (svec->param_data, svec->nitens, 1);

    ncm_vector_scale (autocor, 1.0 / svec->param_data[0]);

    return autocor;
  }
#else
  g_error ("ncm_stats_vec_get_autocorr: recompile NumCosmo with fftw support.");

  return NULL;

#endif /* HAVE_FFTW3 */
}

/**
 * ncm_stats_vec_get_subsample_autocorr:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 * @subsample: size of the subsample ($>0$)
 *
 * Calculates the autocorrelation vector, the j-th element represent
 * the selfcorrelation with lag-j using the @subsample parameter.
 *
 * The returning vector use the internal memory allocation and will
 * change with subsequent calls to ncm_stats_vec_get_autocorr().
 *
 * Returns: (transfer full): the autocorrelation vector.
 */
NcmVector *
ncm_stats_vec_get_subsample_autocorr (NcmStatsVec *svec, guint p, guint subsample)
{
#ifdef HAVE_FFTW3
  _ncm_stats_vec_get_autocov (svec, p, subsample, 0);
  g_assert_cmpuint (svec->nitens, >=, subsample);
  {
    NcmVector *autocor = ncm_vector_new_data_dup (svec->param_data, svec->nitens / subsample, 1);

    ncm_vector_scale (autocor, 1.0 / svec->param_data[0]);

    return autocor;
  }
#else
  g_error ("ncm_stats_vec_get_subsample_autocorr: recompile NumCosmo with fftw support.");

  return NULL;

#endif /* HAVE_FFTW3 */
}

/**
 * ncm_stats_vec_fit_ar_model:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 * @order: max order
 * @ar_crit: a #NcmStatsVecARType
 * @rho: (inout) (nullable): the vector containing the ar(@p) model parameters
 * @pacf: (inout) (nullable):  the vector containing the partial autocorrelations
 * @ivar: (out): innovations variance
 * @c_order: (out): the actual order calculated
 *
 * If order is zero the value of floor $\left[10 log_{10}(s) \right]$, where $s$
 * is the number of points.
 *
 * Returns: TRUE if @c_order is equal to @order.
 */
gboolean
ncm_stats_vec_fit_ar_model (NcmStatsVec *svec, guint p, const guint order, NcmStatsVecARType ar_crit, NcmVector **rho, NcmVector **pacf, gdouble *ivar, guint *c_order)
{
#ifdef HAVE_FFTW3
  _ncm_stats_vec_get_autocov (svec, p, 1, 0);
  {
    const gint aorder          = (order == 0) ? GSL_MIN (GSL_MAX (svec->nitens - 2, 1), floor (10 * log10 (svec->nitens))) : order;
    NcmVector *M               = ncm_vector_new (2 * aorder + 1);
    const gdouble dlev_tol     = 1.0e-3;
    gboolean allocated_here[2] = {FALSE, FALSE};
    gint i;

    g_assert_cmpuint (svec->nitens, >, order + 1);

    if (*rho != NULL)
    {
      g_assert_cmpuint (ncm_vector_len (*rho), >=, aorder);
    }
    else
    {
      *rho              = ncm_vector_new (aorder);
      allocated_here[0] = TRUE;
    }

    if (*pacf != NULL)
    {
      g_assert_cmpuint (ncm_vector_len (*pacf), >=, aorder);
    }
    else
    {
      *pacf             = ncm_vector_new (aorder);
      allocated_here[1] = TRUE;
    }

    ncm_vector_fast_set (M, aorder, svec->param_data[0]);

    for (i = 0; i < aorder; i++)
    {
      const gdouble a_i = svec->param_data[i + 1];

      ncm_vector_fast_set (M, aorder + i + 1, a_i);
      ncm_vector_fast_set (M, aorder - i - 1, a_i);
    }

    d_lev_inner (ncm_vector_data (*rho),
                 ncm_vector_ptr (M, aorder),
                 aorder, svec->param_data + 1, dlev_tol, dlev_tol, 6, 0,
                 ncm_vector_data (*pacf));

    {
      const gdouble n = svec->nitens;
      gdouble var     = ivar[0] = ncm_stats_vec_get_var (svec, p);

      switch (ar_crit)
      {
        case NCM_STATS_VEC_AR_NONE:
          c_order[0] = aorder;

          for (i = 0; i < aorder; i++)
            var *= 1.0 - gsl_pow_2 (ncm_vector_get (*pacf, i));

          ivar[0] = var;

          break;
        case NCM_STATS_VEC_AR_FPE:
        {
          gdouble crit     = var;
          gdouble min_crit = crit;

          c_order[0] = 0;

          for (i = 0; i < aorder; i++)
          {
            const gdouble p = 1.0 + i;

            var *= 1.0 - gsl_pow_2 (ncm_vector_get (*pacf, i));
            crit = var * (n + p) / (n - p);

            if (crit < min_crit)
            {
              c_order[0] = i + 1;
              min_crit   = crit;
              ivar[0]    = var;
            }
          }

          break;
        }
        case NCM_STATS_VEC_AR_AIC:
        {
          gdouble crit     = n * log (var) + 2.0;
          gdouble min_crit = crit;

          c_order[0] = 0;

          for (i = 0; i < aorder; i++)
          {
            const gdouble p = 1.0 + i;

            var *= 1.0 - gsl_pow_2 (ncm_vector_get (*pacf, i));

            crit = n * log (var) + 2.0 * (p + 1.0);

            if (crit < min_crit)
            {
              c_order[0] = i + 1;
              min_crit   = crit;
              ivar[0]    = var;
            }
          }

          break;
        }
        case NCM_STATS_VEC_AR_AICC:
        {
          gdouble crit     = n * log (var) + 2.0 * n / (n - 2.0);
          gdouble min_crit = crit;

          c_order[0] = 0;

          for (i = 0; i < aorder; i++)
          {
            const gdouble p = 1.0 + i;

            var *= 1.0 - gsl_pow_2 (ncm_vector_get (*pacf, i));
            crit = n * log (var) + 2.0 * n * (p + 1.0) / (n - p - 2.0);

            if (crit < min_crit)
            {
              c_order[0] = i + 1;
              min_crit   = crit;
              ivar[0]    = var;
            }
          }

          break;
        }
        default:
          g_assert_not_reached ();
          break;
      }

      ivar[0] *= (n - 1.0) / (n - (c_order[0] + 1.0));
    }

    if (c_order[0] == 0)
    {
      if (allocated_here[0])
        ncm_vector_clear (rho);

      if (allocated_here[1])
        ncm_vector_clear (pacf);
    }
    else if (c_order[0] != (guint) aorder)
    {
      NcmVector *c_rho  = ncm_vector_get_subvector (*rho,  0, c_order[0]);
      NcmVector *c_pacf = ncm_vector_get_subvector (*pacf, 0, c_order[0]);

      ncm_vector_clear (rho);
      ncm_vector_clear (pacf);

      *rho  = c_rho;
      *pacf = c_pacf;

      d_lev_inner (ncm_vector_data (c_rho),
                   ncm_vector_ptr (M, aorder),
                   c_order[0], svec->param_data + 1, dlev_tol, dlev_tol, 6, 0,
                   ncm_vector_data (c_pacf));
    }

    ncm_vector_free (M);

    return ((guint) aorder == c_order[0]);
  }
#else
  g_error ("ncm_stats_vec_get_autocorr: recompile NumCosmo with fftw support.");

  return FALSE;

#endif /* HAVE_FFTW3 */
}

/**
 * ncm_stats_vec_ar_ess:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 * @ar_crit: a #NcmStatsVecARType
 * @spec0: (out): spectral density at zero
 * @c_order: (out): @ar_crit determined order
 *
 * Calculates the effective sample size for the parameter @p.
 *
 * Returns: the effective sample size.
 */
gdouble
ncm_stats_vec_ar_ess (NcmStatsVec *svec, guint p, NcmStatsVecARType ar_crit, gdouble *spec0, guint *c_order)
{
  NcmVector *rho = NULL, *pacf = NULL;
  gdouble ivar = 0.0;
  guint order  = 0;

  g_assert_cmpuint (p, <, svec->len);

  if (svec->nitens <= 1)
    return svec->nitens;

  while (ncm_stats_vec_fit_ar_model (svec, p, order, ar_crit, &rho, &pacf, &ivar, c_order) && (2 * c_order[0] + 1 < svec->nitens))
  {
    ncm_vector_clear (&rho);
    ncm_vector_clear (&pacf);

    order = 2 * c_order[0];
  }

  spec0[0] = ivar;

  if (c_order[0] > 0)
    spec0[0] *= 1.0 / gsl_pow_2 (1.0 - ncm_vector_sum_cpts (rho));

  ncm_vector_clear (&rho);
  ncm_vector_clear (&pacf);

  return svec->nitens * ncm_stats_vec_get_var (svec, p) / spec0[0];
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
 * Estimate mean $\mu$ and standard deviation $\sigma$ fitting the paramater @p
 * using robust regression. Computes the time $t_0$ where the parameter @p falls
 * within the $\alpha\sigma$ from $\mu$, where $\alpha$ is implicitly defined by
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
  const gint maxiter  = 100;
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
 * Applies the Heidelberger and Welch’s convergence diagnostic
 * applying @ntests Schruben tests sequentially, if @ntests == 0
 * it will use the default 10 tests. The variable @bindex will
 * contains the smallest index where all p-values are smaller than
 * @pvalue, if @pvalue is zero it used the default value of $0.05$.
 *
 * If the test is not satisfied by any index @bindex will contain
 * -1 and the return vector the p-values considering the whole system.
 *
 * See:
 * - [Heidelberger (1981)][XHeidelberger1981]
 * - [Schruben (1982)][XSchruben1982]
 * - [Heidelberger (1983)][XHeidelberger1983]
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

  for (i = 0; i < (gint) svec->len; i++)
  {
    ncm_stats_vec_ar_ess (chunk, i, NCM_STATS_VEC_AR_AICC, ncm_vector_ptr (spec0, i), &c_order);
    g_array_append_val (ar_order, c_order);
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
 * Computes the empirical cumulative and the mean used to build
 * the Heidelberger and Welch’s convergence diagnostic.
 *
 * See ncm_stats_vec_heidel_diag().
 *
 * Returns: (transfer full): a #NcmVector containing the empirical cumulative distribution.
 */
NcmVector *
ncm_stats_vec_visual_heidel_diag (NcmStatsVec *svec, const guint p, const guint fi, gdouble *mean, gdouble *var)
{
  NcmStatsVec *chunk  = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, TRUE);
  const guint nitens  = svec->nitens - fi;
  gdouble spec0       = 0.0;
  guint c_order       = 0;
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

    ncm_stats_vec_set (chunk, 0, p_val);
    ncm_stats_vec_update (chunk);
  }

  ncm_stats_vec_ar_ess (chunk, 0, NCM_STATS_VEC_AR_AICC, &spec0, &c_order);

  mean[0] = ncm_stats_vec_get_mean (chunk, 0);
  var[0]  = spec0 * nitens;

  ncm_stats_vec_clear (&chunk);

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
  NcmStatsVec *chunk  = ncm_stats_vec_new (svec->len, NCM_STATS_VEC_VAR, TRUE);
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

    ncm_stats_vec_append (chunk, row_i, FALSE);

    if ((i == 0) || ((i % block == 0) && (j >= 99)))
    {
      gdouble min_ess = GSL_POSINF;
      guint cur_size  = size - i;
      guint lwp_order = 0;
      guint k, lwp    = 0;

      for (k = 0; k < svec->len; k++)
      {
        gdouble spec0       = 0.0;
        const gdouble ess   = ncm_stats_vec_ar_ess (chunk, k, NCM_STATS_VEC_AR_AICC, &spec0, &lwp_order);
        const gdouble c_ess = GSL_MIN (cur_size, ess);

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
  ncm_stats_vec_clear (&chunk);

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
    gint ret;

    ret = gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, ncm_matrix_gsl (E), ncm_matrix_gsl (y), 0.0, ncm_matrix_gsl (z));
    NCM_TEST_GSL_RESULT ("ncm_stats_vec_compute_cov_robust_ogk", ret);

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

  {
    gint ret;

    ret = gsl_blas_dsyrk (CblasUpper, CblasTrans, 1.0, ncm_matrix_gsl (E), 0.0, ncm_matrix_gsl (cov));
    NCM_TEST_GSL_RESULT ("ncm_stats_vec_compute_cov_robust_ogk", ret);
  }

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
 * ncm_stats_vec_get_autocorr_tau:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 * @max_lag: max lag in the computation
 *
 * Calculates the integrated autocorrelation time for the parameter @p
 * using all rows of data.
 *
 * If @max_lag is 0 or larger than the current number of itens than it use
 * the current number of itens as @max_lag.
 *
 * Returns: the integrated autocorrelation time of the whole data.
 */
gdouble
ncm_stats_vec_get_autocorr_tau (NcmStatsVec *svec, const guint p, const guint max_lag)
{
#ifdef HAVE_FFTW3
  guint i;
  gdouble tau          = 0.0;
  const guint Imax_lag = (max_lag == 0) ? svec->nitens / 10 : max_lag;
  const guint Fmax_lag = (Imax_lag > 1000) ? 1000 : Imax_lag;

  _ncm_stats_vec_get_autocov (svec, p, 1, 0);

  g_assert_cmpuint (Fmax_lag, >, 0);
  g_assert_cmpuint (Fmax_lag, <, svec->nitens);

  {
    for (i = 1; i < Fmax_lag + 1; i++)
    {
      const gdouble rho_i = svec->param_data[i] / svec->param_data[0];

      tau += rho_i;
    }
  }

  tau = 1.0 + 2.0 * tau;

  return tau;

#else
  g_error ("ncm_stats_vec_get_autocorr: recompile NumCosmo with fftw support.");

  return 0.0;

#endif /* HAVE_FFTW3 */
}

/**
 * ncm_stats_vec_get_subsample_autocorr_tau:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 * @subsample: size of the subsample ($>0$)
 * @max_lag: max lag in the computation
 *
 * Calculates the integrated autocorrelation time for the parameter @p
 * using the @subsample parameter.
 *
 * Returns: the integrated autocorrelation time of data with @subsample.
 */
gdouble
ncm_stats_vec_get_subsample_autocorr_tau (NcmStatsVec *svec, const guint p, const guint subsample, const guint max_lag)
{
#ifdef HAVE_FFTW3
  guint i;
  gdouble tau          = 0.0;
  guint eff_nitens     = svec->nitens / subsample;
  const guint Imax_lag = (max_lag == 0) ? eff_nitens / 10 : max_lag;
  const guint Fmax_lag = (Imax_lag > 1000) ? 1000 : Imax_lag;

  _ncm_stats_vec_get_autocov (svec, p, subsample, 0);

  g_assert_cmpuint (Fmax_lag, >, 0);
  g_assert_cmpuint (Fmax_lag, <, eff_nitens);

  for (i = 1; i < Fmax_lag + 1; i++)
  {
    const gdouble rho_i = svec->param_data[i] / svec->param_data[0];

    tau += rho_i;
  }

  tau = 1.0 + 2.0 * tau;

  return tau;

#else
  g_error ("ncm_stats_vec_get_autocorr: recompile NumCosmo with fftw support.");

  return 0.0;

#endif /* HAVE_FFTW3 */
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
 * Same as ncm_stats_vec_update_weight() assuming weigth equal to one.
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
 * Gets the number of itens added to the object;
 *
 * Returns: the number of itens added.
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

