/***************************************************************************
 *            ncm_stats_vec.c
 *
 *  Fri August 02 13:41:01 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_vec.c
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
#include "math/ncm_cfg.h"
#include "ncm_enum_types.h"

#include <math.h>
#include <gsl/gsl_math.h>

#include "math/gsl_rstat.h"
#include "math/rquantile.c"

enum
{
  PROP_0,
  PROP_LEN,
  PROP_TYPE,
  PROP_SAVE_X,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmStatsVec, ncm_stats_vec, G_TYPE_OBJECT);

static void
ncm_stats_vec_init (NcmStatsVec *svec)
{
  svec->t        = NCM_STATS_VEC_TYPES_LEN;
  svec->weight   = 0.0;
  svec->weight2  = 0.0;
  svec->bias_wt  = 0.0;
/*
 * Increment calculation
  svec->mean_inc = 0.0;
  svec->var_inc  = 0.0;
  svec->cov_inc  = 0.0;
 */
  svec->nitens   = 0;
  svec->x        = NULL;
  svec->mean     = NULL;
  svec->var      = NULL;
  svec->cov      = NULL;
  svec->real_cov = NULL;
  svec->saved_x  = NULL;
  svec->save_x   = FALSE;

  svec->q_array  = g_ptr_array_new ();
  g_ptr_array_set_free_func (svec->q_array, (GDestroyNotify) gsl_rstat_quantile_free);
  
#ifdef NUMCOSMO_HAVE_FFTW3
  svec->fft_size       = 0;
  svec->fft_plan_size  = 0;
  svec->param_data     = NULL;
  svec->param_fft      = NULL;
  svec->param_r2c      = NULL;
  svec->param_c2r      = NULL;
#endif /* NUMCOSMO_HAVE_FFTW3 */
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
    svec->save_x = FALSE;
  }

  g_clear_pointer (&svec->q_array, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_vec_parent_class)->dispose (object);
}

static void
_ncm_stats_vec_finalize (GObject *object)
{
  NcmStatsVec *svec = NCM_STATS_VEC (object);

#ifdef NUMCOSMO_HAVE_FFTW3
  g_clear_pointer (&svec->param_fft,  fftw_free);
  g_clear_pointer (&svec->param_data, fftw_free);
  g_clear_pointer (&svec->param_c2r,  fftw_destroy_plan);
  g_clear_pointer (&svec->param_r2c,  fftw_destroy_plan);
#endif /* NUMCOSMO_HAVE_FFTW3 */

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_vec_parent_class)->finalize (object);
}

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
      g_ptr_array_set_free_func (svec->saved_x, (GDestroyNotify)ncm_vector_free);
    }

    switch (svec->t)
    {
      case NCM_STATS_VEC_COV:
        g_assert_cmpuint (svec->len, >, 1);
        g_assert (svec->cov == NULL);
        svec->cov = ncm_matrix_new (svec->len, svec->len);
        ncm_matrix_set_zero (svec->cov);
      case NCM_STATS_VEC_VAR:
        g_assert (svec->var == NULL);
        svec->var = ncm_vector_new (svec->len);
        ncm_vector_set_zero (svec->var);
      case NCM_STATS_VEC_MEAN:
        g_assert (svec->x == NULL);
        g_assert (svec->mean == NULL);
        svec->x = ncm_vector_new (svec->len);
        svec->mean = ncm_vector_new (svec->len);
        ncm_vector_set_zero (svec->x);
        ncm_vector_set_zero (svec->mean);
        break;
      case NCM_STATS_VEC_TYPES_LEN:
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_stats_vec_class_init (NcmStatsVecClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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
 * Returns: (transfer full): FIXME
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
    g_ptr_array_set_free_func (svec->saved_x, (GDestroyNotify)ncm_vector_free);
    g_ptr_array_set_size (svec->saved_x, 0);
  }

  svec->weight   = 0.0;
  svec->weight2  = 0.0;
  svec->bias_wt  = 0.0;
#ifdef NCM_STATS_VEC_INC
  svec->mean_inc = 0.0;
  svec->var_inc  = 0.0;
  svec->cov_inc  = 0.0;
#endif /* NCM_STATS_VEC_INC */
  svec->nitens   = 0;

  switch (svec->t)
  {
    case NCM_STATS_VEC_COV:
      g_assert (svec->cov != NULL);
      ncm_matrix_set_zero (svec->cov);
    case NCM_STATS_VEC_VAR:
      g_assert (svec->var != NULL);
      ncm_vector_set_zero (svec->var);
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
    const gdouble p = ((gsl_rstat_quantile_workspace *)g_ptr_array_index (svec->q_array, 0))->p;

    g_ptr_array_set_size (svec->q_array, 0);
    
    for (i = 0; i < svec->len; i++)
    {
      gsl_rstat_quantile_workspace *qws_i = gsl_rstat_quantile_alloc (p);
      g_ptr_array_add (svec->q_array, qws_i);
    }
  }
}

static void
_ncm_stats_vec_update_from_vec_weight (NcmStatsVec *svec, const gdouble w, NcmVector *x)
{
  const gdouble curweight = svec->weight + w;
  guint i;
#ifdef NCM_STATS_VEC_INC
  gdouble inc = 0.0;
  svec->mean_inc = 0.0;
  svec->var_inc  = 0.0;
  svec->cov_inc  = 0.0;
#endif /* NCM_STATS_VEC_INC */

  switch (svec->t)
  {
    case NCM_STATS_VEC_MEAN:
    {
      for (i = 0; i < svec->len; i++)
      {
        const gdouble mean_i = ncm_vector_fast_get (svec->mean, i);
        const gdouble x_i = ncm_vector_fast_get (x, i);
        const gdouble delta_i = x_i - mean_i;
        const gdouble R_i = delta_i * w / curweight;
        ncm_vector_fast_set (svec->mean, i, mean_i + R_i);

#ifdef NCM_STATS_VEC_INC
        inc = sqrt (curweight) / (fabs (mean_i / R_i) + 1.0);
        svec->mean_inc = GSL_MAX (svec->mean_inc, inc);
#endif /* NCM_STATS_VEC_INC */
      }
      break;
    }
    case NCM_STATS_VEC_VAR:
    {
      for (i = 0; i < svec->len; i++)
      {
        const gdouble mean_i  = ncm_vector_fast_get (svec->mean, i);
        const gdouble x_i     = ncm_vector_fast_get (x, i);
        const gdouble delta_i = x_i - mean_i;
        const gdouble R_i     = delta_i * w / curweight;
        const gdouble var     = ncm_vector_fast_get (svec->var, i);
        const gdouble dvar    = svec->weight * delta_i * R_i;

        ncm_vector_fast_set (svec->mean, i, mean_i + R_i);
        ncm_vector_fast_set (svec->var, i, var + dvar);

#ifdef NCM_STATS_VEC_INC
        inc = sqrt (curweight) / (fabs (mean_i / R_i) + 1.0);
        svec->mean_inc = GSL_MAX (svec->mean_inc, inc);

        inc = 1.0 / (fabs (var / dvar) + 1.0);
        svec->var_inc = GSL_MAX (svec->var_inc, inc);
#endif /* NCM_STATS_VEC_INC */
      }
      break;
    }
    case NCM_STATS_VEC_COV:
    {
      for (i = 0; i < svec->len; i++)
      {
        guint j;
        gdouble mean_i = ncm_vector_fast_get (svec->mean, i);
        const gdouble x_i     = ncm_vector_fast_get (x, i);
        const gdouble delta_i = x_i - mean_i;
        const gdouble R_i     = delta_i * w / curweight;
        const gdouble var     = ncm_vector_fast_get (svec->var, i);
        const gdouble dvar    = svec->weight * delta_i * R_i;

#ifdef NCM_STATS_VEC_INC
        inc = sqrt (curweight) / (fabs (mean_i / R_i) + 1.0);
        svec->mean_inc = GSL_MAX (svec->mean_inc, inc);
        inc = 1.0 / (fabs (var / dvar) + 1.0);
        svec->var_inc = GSL_MAX (svec->var_inc, inc);
#endif /* NCM_STATS_VEC_INC */

        mean_i += R_i;
        ncm_vector_fast_set (svec->mean, i, mean_i);
        ncm_vector_fast_set (svec->var, i, var + dvar);
        for (j = i + 1; j < svec->len; j++)
        {
          const gdouble x_j = ncm_vector_fast_get (x, j);
          const gdouble mean_j = ncm_vector_fast_get (svec->mean, j);
          const gdouble dC_ij = w * (x_i - mean_i) * (x_j - mean_j);
          const gdouble oC_ij = ncm_matrix_get (svec->cov, i, j);
          const gdouble C_ij = oC_ij + dC_ij;
          ncm_matrix_set (svec->cov, i, j, C_ij);
          ncm_matrix_set (svec->cov, j, i, C_ij);

#ifdef NCM_STATS_VEC_INC
          inc = 1.0 / (fabs (oC_ij / dC_ij) + 1.0);
          svec->cov_inc = GSL_MAX (svec->cov_inc, inc);
#endif /* NCM_STATS_VEC_INC */
        }
      }
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }

  svec->nitens++;
  svec->weight = curweight;
  svec->weight2 += w * w;
  svec->bias_wt = 1.0 / (svec->weight - svec->weight2 / svec->weight);

  if (svec->q_array->len == svec->len)
  {
    guint i;

    for (i = 0; i < svec->len; i++)
    {
      const gdouble x_i = ncm_vector_fast_get (x, i);
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
  _ncm_stats_vec_update_from_vec_weight (svec, w, svec->x);
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
  _ncm_stats_vec_update_from_vec_weight (svec, w, x);
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
  _ncm_stats_vec_update_from_vec_weight (svec, w, x);

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
    _ncm_stats_vec_update_from_vec_weight (svec, 1.0, x);
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
    const guint cp_len = data->len;
    const guint old_len = svec->saved_x->len;
    const guint new_len = old_len + cp_len;

    g_ptr_array_set_size (svec->saved_x, new_len);
    memmove (&svec->saved_x->pdata[cp_len], svec->saved_x->pdata, sizeof (gpointer) * old_len);
  }

  for (i = 0; i < data->len; i++)
  {
    NcmVector *x = g_ptr_array_index (data, i);
    _ncm_stats_vec_update_from_vec_weight (svec, 1.0, x);
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
 * will disconsider the weights.
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
          const gdouble x_j = ncm_vector_fast_get (x, j);
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
#ifdef NUMCOSMO_HAVE_FFTW3
  guint effsize = 2 * size;
  if (!svec->save_x)
    g_error ("_ncm_stats_vec_get_autocorr_alloc: NcmStatsVec must have saved data to calculate autocorrelation.");

  effsize = exp2 (ceil (log2 (effsize)));

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
    svec->param_r2c = fftw_plan_dft_r2c_1d (effsize, svec->param_data, svec->param_fft, fftw_default_flags | FFTW_DESTROY_INPUT);
    svec->param_c2r = fftw_plan_dft_c2r_1d (effsize, svec->param_fft, svec->param_data, fftw_default_flags | FFTW_DESTROY_INPUT);
    ncm_cfg_save_fftw_wisdom ("ncm_stats_vec_autocorr_%u", effsize);
    svec->fft_plan_size = effsize;
    /*g_debug ("# _ncm_stats_vec_get_autocorr_alloc: calculated  wisdown %u\n", effsize);*/
  }

#endif /* NUMCOSMO_HAVE_FFTW3 */
}

static void
_ncm_stats_vec_get_autocov (NcmStatsVec *svec, guint p, guint subsample, guint pad)
{
#ifdef NUMCOSMO_HAVE_FFTW3
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
#endif /* NUMCOSMO_HAVE_FFTW3 */
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
 * Returns: (transfer full): a read only autocorrelation vector.
 */
NcmVector *
ncm_stats_vec_get_autocorr (NcmStatsVec *svec, guint p)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  _ncm_stats_vec_get_autocov (svec, p, 1, 0);
  {
    NcmVector *autocor = ncm_vector_new_data_static (svec->param_data, svec->nitens, 1);
    ncm_vector_scale (autocor, 1.0 / svec->param_data[0]);
    return autocor;
  }
#else
  g_error ("ncm_stats_vec_get_autocorr: recompile NumCosmo with fftw support.");
  return NULL;
#endif /* NUMCOSMO_HAVE_FFTW3 */
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
 * Returns: (transfer full): a read only autocorrelation vector.
 */
NcmVector *
ncm_stats_vec_get_subsample_autocorr (NcmStatsVec *svec, guint p, guint subsample)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  _ncm_stats_vec_get_autocov (svec, p, subsample, 0);
  g_assert_cmpuint (svec->nitens, >=, subsample);
  {
    NcmVector *autocor = ncm_vector_new_data_static (svec->param_data, svec->nitens / subsample, 1);
    ncm_vector_scale (autocor, 1.0 / svec->param_data[0]);
    return autocor;
  }
#else
  g_error ("ncm_stats_vec_get_subsample_autocorr: recompile NumCosmo with fftw support.");
  return NULL;
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

/**
 * ncm_stats_vec_get_autocorr_tau:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 * @max_lag: max lag in the computation
 * @min_rho: minimum autocorrelation to be considered in the sum
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
ncm_stats_vec_get_autocorr_tau (NcmStatsVec *svec, guint p, guint max_lag, const gdouble min_rho)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  guint i;
  gdouble tau = 0.0;

  _ncm_stats_vec_get_autocov (svec, p, 1, 0);
  if (max_lag == 0 || max_lag > svec->nitens)
    max_lag = svec->nitens;

  for (i = 1; i < max_lag; i++)
  {
    gdouble rho_i = svec->param_data[i] / svec->param_data[0];
    if (rho_i < min_rho)
      break;
    /*printf ("# self cor %u % 20.15g\n", i, rho_i);*/
    tau += rho_i;
  }

  tau = 1.0 + 2.0 * tau;

  return tau;
#else
  g_error ("ncm_stats_vec_get_autocorr: recompile NumCosmo with fftw support.");
  return 0.0;
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

/**
 * ncm_stats_vec_get_subsample_autocorr_tau:
 * @svec: a #NcmStatsVec
 * @p: parameter id
 * @subsample: size of the subsample ($>0$)
 * @max_lag: max lag in the computation
 * @min_rho: minimum autocorrelation to be considered in the sum
 *
 * Calculates the integrated autocorrelation time for the parameter @p
 * using the @subsample parameter.
 *
 * Returns: the integrated autocorrelation time of data with @subsample.
 */
gdouble
ncm_stats_vec_get_subsample_autocorr_tau (NcmStatsVec *svec, guint p, guint subsample, guint max_lag, const gdouble min_rho)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  guint i;
  gdouble tau = 0.0;
  guint eff_nitens = svec->nitens / subsample;

  _ncm_stats_vec_get_autocov (svec, p, subsample, 0);
  if (max_lag == 0 || max_lag > eff_nitens)
    max_lag = eff_nitens;

  for (i = 1; i < max_lag; i++)
  {
    gdouble rho_i = svec->param_data[i] / svec->param_data[0];
    if (rho_i < min_rho)
      break;
    tau += rho_i;
  }

  tau = 1.0 + 2.0 * tau;

  return tau;
#else
  g_error ("ncm_stats_vec_get_autocorr: recompile NumCosmo with fftw support.");
  return 0.0;
#endif /* NUMCOSMO_HAVE_FFTW3 */
}


/**
 * ncm_stats_vec_peek_x:
 * @svec: a #NcmStatsVec
 *
 * Returns the vector containing the current value of the random variables.
 *
 * Returns: (transfer none): the random variables vector.
 */
/**
 * ncm_stats_vec_set:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 * @x_i: the value of the @i-th variable
 *
 * Sets the value of the current @i-th random variable to @x_i.
 *
 */
/**
 * ncm_stats_vec_get:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Returns the value of the current @i-th random variable.
 *
 * Returns: @i-th random variable.
 */
/**
 * ncm_stats_vec_update:
 * @svec: a #NcmStatsVec.
 *
 * Same as ncm_stats_vec_update_weight() assuming weigth equal to one.
 *
 */
/**
 * ncm_stats_vec_get_mean:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Return the current value of the variable mean, i.e., $\bar{x}_n$.
 *
 * Returns: $\bar{x}_n$.
 */
/**
 * ncm_stats_vec_get_var:
 * @svec: a #NcmStatsVec
 * @i: a variable index
 *
 * Return the current value of the variable variance, i.e., $Var_n$.
 *
 * Returns: $Var_n$.
 */
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
/**
 * ncm_stats_vec_get_weight:
 * @svec: a #NcmStatsVec
 *
 * Return the current value of the weight, for non-weighted means this is simply
 * the number of elements.
 *
 * Returns: $W_n$.
 */
/**
 * ncm_stats_vec_get_mean_vector:
 * @svec: a #NcmStatsVec
 * @mean: a #NcmVector
 * @offset: first parameter index
 *
 * Copy the current value of the means to the vector @mean starting from parameter @offset.
 *
 */
/**
 * ncm_stats_vec_peek_mean:
 * @svec: a #NcmStatsVec
 * 
 * Gets the local mean vector.
 * 
 * Returns: (transfer none): the internal mean #NcmVector.
 */
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
/**
 * ncm_stats_vec_nrows:
 * @svec: a #NcmStatsVec
 *
 * Gets the number of saved rows, this function fails if the object
 * was not created with save_x == TRUE;
 *
 * Returns: the number of saved rows.
 */
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
