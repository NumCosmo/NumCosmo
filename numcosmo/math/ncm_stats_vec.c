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
 * @title: Statistics vector object
 * @short_description: An online statistics vector object.
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
 * <example>
 *   <title>Using a NcmStatsVec.</title>
 *   <programlisting>
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
 *   </programlisting>
 * </example>
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
  svec->saved_x  = NULL;
  svec->save_x   = FALSE;
}

static void
_ncm_stats_vec_dispose (GObject *object)
{
  NcmStatsVec *svec = NCM_STATS_VEC (object);
  ncm_vector_clear (&svec->x);
  ncm_vector_clear (&svec->mean);
  ncm_vector_clear (&svec->var);
  ncm_matrix_clear (&svec->cov);
  if (svec->saved_x != NULL)
  {
    g_ptr_array_unref (svec->saved_x);
    svec->saved_x = NULL;
    svec->save_x = FALSE;
  }
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_vec_parent_class)->dispose (object);
}

static void
_ncm_stats_vec_finalize (GObject *object)
{

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
 * @len: number of random variables.
 * @t: type of statistics to be calculated.
 * @save_x: whenever to save each vector x.
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
 * @svec: a #NcmStatsVec.
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
 * @svec: a #NcmStatsVec.
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
 * @svec: a #NcmStatsVec.
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
 * @svec: a #NcmStatsVec.
 * 
 * Reset all data in @svec.
 * 
 */
void 
ncm_stats_vec_reset (NcmStatsVec *svec)
{
  if (svec->save_x)
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
}

/**
 * ncm_stats_vec_update_weight:
 * @svec: a #NcmStatsVec.
 * @w: The statistical weight.
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
    g_ptr_array_add (svec->saved_x, ncm_vector_dup (svec->x));
  ncm_vector_set_zero (svec->x);
}

/**
 * ncm_stats_vec_append_data:
 * @svec: a #NcmStatsVec.
 * @data: (element-type NcmVector): a #GPtrArray containing #NcmVector s to be added.
 * @dup: a #gboolean.
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
  g_assert (svec->save_x);

  for (i = 0; i < data->len; i++)
  {
    NcmVector *x = g_ptr_array_index (data, i);
    _ncm_stats_vec_update_from_vec_weight (svec, 1.0, x);
    if (dup)
      g_ptr_array_add (svec->saved_x, ncm_vector_dup (x));
    else
      g_ptr_array_add (svec->saved_x, ncm_vector_ref (x));
  }
}

/**
 * ncm_stats_vec_prepend_data:
 * @svec: a #NcmStatsVec.
 * @data: (element-type NcmVector): a #GPtrArray containing #NcmVector s to be added.
 * @dup: a boolean.
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
  const guint cp_len = data->len;
  const guint old_len = svec->saved_x->len;
  const guint new_len = old_len + cp_len;
  g_assert (svec->save_x);

  g_ptr_array_set_size (svec->saved_x, new_len);
  memmove (&svec->saved_x->pdata[cp_len], svec->saved_x->pdata, sizeof (gpointer) * old_len);

  for (i = 0; i < data->len; i++)
  {
    NcmVector *x = g_ptr_array_index (data, i);
    _ncm_stats_vec_update_from_vec_weight (svec, 1.0, x);
    if (dup)
      g_ptr_array_index (svec->saved_x, i) = ncm_vector_dup (x);
    else
      g_ptr_array_index (svec->saved_x, i) = ncm_vector_ref (x);
  }
}

/**
 * ncm_stats_vec_peek_x:
 * @svec: a #NcmStatsVec.
 * 
 * Returns the vector containing the current value of the random variables.
 * 
 * Returns: (transfer none): the random variables vector.
 */
/**
 * ncm_stats_vec_set:
 * @svec: a #NcmStatsVec.
 * @i: a variable index.
 * @x_i: the value of the @i-th variable.
 * 
 * Sets the value of the current @i-th random variable to @x_i.
 * 
 */
/**
 * ncm_stats_vec_get:
 * @svec: a #NcmStatsVec.
 * @i: a variable index.
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
 * @svec: a #NcmStatsVec.
 * @i: a variable index.
 * 
 * Return the current value of the variable mean, i.e., $\bar{x}_n$.
 * 
 * Returns: $\bar{x}_n$.
 */
/**
 * ncm_stats_vec_get_var:
 * @svec: a #NcmStatsVec.
 * @i: a variable index.
 * 
 * Return the current value of the variable variance, i.e., $Var_n$.
 * 
 * Returns: $Var_n$.
 */
/**
 * ncm_stats_vec_get_sd:
 * @svec: a #NcmStatsVec.
 * @i: a variable index.
 * 
 * Return the current value of the variable standard deviation, 
 * i.e., $\sigma_n \equiv sqrt (Var_n)$.
 * 
 * Returns: $\sigma_n$
 */
/**
 * ncm_stats_vec_get_cov:
 * @svec: a #NcmStatsVec.
 * @i: a variable index.
 * @j: a variable index.
 * 
 * Return the current value of the variance between the @i-th and the @j-th 
 * variables, i.e., $Cov_{ij}$.
 * 
 * Returns: $Cov_{ij}$.
 */
/**
 * ncm_stats_vec_get_cor:
 * @svec: a #NcmStatsVec.
 * @i: a variable index.
 * @j: a variable index.
 * 
 * Return the current value of the correlation between the @i-th and the @j-th 
 * variables, i.e., $$Cor_{ij} \equiv \frac{Cov_{ij}}{\sigma_i\sigma_j}.$$
 * 
 * Returns: $Cor_{ij}$.
 */
/**
 * ncm_stats_vec_get_weight:
 * @svec: a #NcmStatsVec.
 * 
 * Return the current value of the weight, for non-weighted means this is simply 
 * the number of elements.
 * 
 * Returns: $W_n$.
 */
/**
 * ncm_stats_vec_get_mean_vector:
 * @svec: a #NcmStatsVec.
 * @mean: a #NcmVector.
 * @offset: first parameter index.
 * 
 * Copy the current value of the means to the vector @mean starting from parameter @offset. 
 * 
 */
/**
 * ncm_stats_vec_get_cov_matrix:
 * @svec: a #NcmStatsVec.
 * @m: a #NcmMatrix.
 * @offset: first parameter index.
 * 
 * Copy the current value of the correlation between the variables to the 
 * matrix @m starting from paramenter @offset. 
 * 
 */
/**
 * ncm_stats_vec_peek_row:
 * @svec: a #NcmStatsVec.
 * @i: the row's index.
 * 
 * The i-th data row used in the statistics, this function fails if the object 
 * was not created with save_x == TRUE;  
 * 
 * Returns: (transfer none): the i-th data row.
 */
