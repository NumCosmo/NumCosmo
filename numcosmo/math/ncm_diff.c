/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_diff.c
 *
 *  Fri July 21 12:59:36 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_diff.c
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_diff
 * @title: NcmDiff
 * @short_description: Numerical differentiation object
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_diff.h"

struct _NcmDiffPrivate
{
  guint maxorder;
  gdouble rs;
  gdouble roff_pad;
	gdouble ini_h;
  GPtrArray *central_tables;
  GPtrArray *forward_tables;
  GPtrArray *backward_tables;
};

typedef struct _NcmDiffTable
{
  NcmVector *h;
  NcmVector *lambda;
} NcmDiffTable;

NcmDiffTable *_ncm_diff_table_new (const guint n);
NcmDiffTable *_ncm_diff_table_dup (NcmDiffTable *dtable);
void _ncm_diff_table_free (gpointer dtable_ptr);

enum
{
  PROP_0,
  PROP_MAXORDER,
  PROP_RS,
  PROP_ROFF_PAD,
	PROP_INI_H,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmDiff, ncm_diff, G_TYPE_OBJECT);

static void
ncm_diff_init (NcmDiff *diff)
{
  diff->priv           = ncm_diff_get_instance_private (diff);
  diff->priv->maxorder = 0;
  diff->priv->rs       = 0.0;
  diff->priv->roff_pad = 0.0;
	diff->priv->ini_h    = 0.0;

  diff->priv->central_tables  = g_ptr_array_new ();
  diff->priv->forward_tables  = g_ptr_array_new ();
  diff->priv->backward_tables = g_ptr_array_new ();

  g_ptr_array_set_free_func (diff->priv->central_tables,  &_ncm_diff_table_free);
  g_ptr_array_set_free_func (diff->priv->forward_tables,  &_ncm_diff_table_free);
  g_ptr_array_set_free_func (diff->priv->backward_tables, &_ncm_diff_table_free);
}

static void
_ncm_diff_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDiff *diff = NCM_DIFF (object);
  g_return_if_fail (NCM_IS_DIFF (object));

  switch (prop_id)
  {
    case PROP_MAXORDER:
      ncm_diff_set_max_order (diff, g_value_get_uint (value));
      break;
    case PROP_RS:
      ncm_diff_set_richardson_step (diff, g_value_get_double (value));
      break;
    case PROP_ROFF_PAD:
      ncm_diff_set_round_off_pad (diff, g_value_get_double (value));
      break;
    case PROP_INI_H:
      ncm_diff_set_ini_h (diff, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_diff_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDiff *diff = NCM_DIFF (object);
  g_return_if_fail (NCM_IS_DIFF (object));

  switch (prop_id)
  {
    case PROP_MAXORDER:
      g_value_set_uint (value, ncm_diff_get_max_order (diff));
      break;
    case PROP_RS:
      g_value_set_double (value, ncm_diff_get_richardson_step (diff));
      break;
    case PROP_ROFF_PAD:
      g_value_set_double (value, ncm_diff_get_round_off_pad (diff));
      break;
    case PROP_INI_H:
      g_value_set_double (value, ncm_diff_get_ini_h (diff));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _ncm_diff_build_diff_tables (NcmDiff *diff);

static void
_ncm_diff_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_diff_parent_class)->constructed (object);
  {
    NcmDiff *diff = NCM_DIFF (object);

    _ncm_diff_build_diff_tables (diff);
  }
}

static void
_ncm_diff_dispose (GObject *object)
{
  NcmDiff *diff = NCM_DIFF (object);

  g_clear_pointer (&diff->priv->central_tables,  g_ptr_array_unref);
  g_clear_pointer (&diff->priv->forward_tables,  g_ptr_array_unref);
  g_clear_pointer (&diff->priv->backward_tables, g_ptr_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_diff_parent_class)->dispose (object);
}

static void
_ncm_diff_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_diff_parent_class)->finalize (object);
}

static void
ncm_diff_class_init (NcmDiffClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_diff_set_property;
  object_class->get_property = &_ncm_diff_get_property;
  object_class->constructed  = &_ncm_diff_constructed;
  object_class->dispose      = &_ncm_diff_dispose;
  object_class->finalize     = &_ncm_diff_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MAXORDER,
                                   g_param_spec_uint ("max-order",
                                                      NULL,
                                                      "Maximum order",
                                                      1, G_MAXUINT, 30,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RS,
                                   g_param_spec_double ("richardson-step",
                                                      NULL,
                                                      "Richardson extrapolation step",
                                                      1.1, G_MAXDOUBLE, 2.0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ROFF_PAD,
                                   g_param_spec_double ("round-off-pad",
                                                      NULL,
                                                      "Round off padding",
                                                      1.1, G_MAXDOUBLE, 1.0e2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_INI_H,
	                                 g_param_spec_double ("ini-h",
	                                                      NULL,
	                                                      "Initial h",
	                                                      GSL_DBL_EPSILON, G_MAXDOUBLE, pow (GSL_DBL_EPSILON, 1.0 / 8.0),
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

NcmDiffTable *
_ncm_diff_table_new (const guint n)
{
  NcmDiffTable *dtable = g_new (NcmDiffTable, 1);
  dtable->h            = ncm_vector_new (n);
  dtable->lambda       = ncm_vector_new (n);

  return dtable;
}

NcmDiffTable *
_ncm_diff_table_dup (NcmDiffTable *dtable)
{
  const guint n            = ncm_vector_len (dtable->h);
  NcmDiffTable *dtable_dup = _ncm_diff_table_new (n);

  ncm_vector_memcpy (dtable_dup->h,      dtable->h);
  ncm_vector_memcpy (dtable_dup->lambda, dtable->lambda);

  return dtable_dup;
}

void 
_ncm_diff_table_free (gpointer dtable_ptr)
{
  NcmDiffTable *dtable = (NcmDiffTable *) dtable_ptr; 

  ncm_vector_free (dtable->h);
  ncm_vector_free (dtable->lambda);
  g_free (dtable);
}

static void
_ncm_diff_build_diff_table (NcmDiff *diff, GPtrArray *tables, const guint maxorder, gdouble (*g) (const guint, gpointer), gpointer user_data)
{
  if (tables->len == maxorder)
    return;
  else if (tables->len > maxorder)
  {
    g_ptr_array_set_size (tables, maxorder);
  }
  else
  {
    guint k;

    for (k = tables->len; k < maxorder; k++)
    {
      const guint n = k + 2;
      if (k == 0)
      {
        const gdouble g0 = g (0, user_data);
        const gdouble g1 = g (1, user_data);
        NcmDiffTable *dtable = _ncm_diff_table_new (n);

        ncm_vector_set (dtable->h, 0, 1.0 / g0);
        ncm_vector_set (dtable->h, 1, 1.0 / g1);

        ncm_vector_set (dtable->lambda, 0, g0 / (g0 - g1));
        ncm_vector_set (dtable->lambda, 1, g1 / (g1 - g0));

        g_ptr_array_add (tables, dtable);
      }
      else
      {
        const guint l         = n - 1;
        const gdouble gl      = g (l, user_data);
        NcmDiffTable *ldtable = g_ptr_array_index (tables, k - 1);
        NcmDiffTable *dtable  = _ncm_diff_table_new (n);
        gdouble lambda_l      = 1.0;
        guint i;
        
        for (i = 0; i < l; i++)
        {
          const gdouble gi       = g (i, user_data);
          const gdouble hi       = ncm_vector_get (ldtable->h, i);
          const gdouble lambda_i = ncm_vector_get (ldtable->lambda, i) / (1.0 - gl / gi);

          ncm_vector_set (dtable->h,      i, hi);
          ncm_vector_set (dtable->lambda, i, lambda_i);

          lambda_l *= 1.0 / (1.0 - gi / gl);
        }

        ncm_vector_set (dtable->h,      l, 1.0 / gl);
        ncm_vector_set (dtable->lambda, l, lambda_l);

        g_ptr_array_add (tables, dtable);
      }      
    }
  }
}

static gdouble
_ncm_diff_central_g (const guint k, gpointer user_data)
{
  NcmDiff *diff    = NCM_DIFF (user_data);
  const gdouble g  = pow (diff->priv->rs, 2 * k);

  return +g;
}

static gdouble
_ncm_diff_forward_g (const guint k, gpointer user_data)
{
  NcmDiff *diff    = NCM_DIFF (user_data);
  const gdouble g  = pow (diff->priv->rs, k);

  return +g;
}

static gdouble
_ncm_diff_backward_g (const guint k, gpointer user_data)
{
  NcmDiff *diff   = NCM_DIFF (user_data);
  const gdouble g = pow (diff->priv->rs, k);

  return -g;
}

static void 
_ncm_diff_build_diff_tables (NcmDiff *diff)
{
  _ncm_diff_build_diff_table (diff, diff->priv->central_tables,  diff->priv->maxorder, _ncm_diff_central_g,  diff);
  _ncm_diff_build_diff_table (diff, diff->priv->forward_tables,  diff->priv->maxorder, _ncm_diff_forward_g,  diff);
  _ncm_diff_build_diff_table (diff, diff->priv->backward_tables, diff->priv->maxorder, _ncm_diff_backward_g, diff);
}

/**
 * ncm_diff_new:
 * 
 * Creates a new #NcmDiff object.
 * 
 * Returns: a new #NcmDiff.
 */
NcmDiff *
ncm_diff_new (void)
{
  NcmDiff *diff = g_object_new (NCM_TYPE_DIFF,
                                NULL);
  return diff;
}

/**
 * ncm_diff_ref:
 * @diff: a #NcmDiff
 *
 * Increase the reference of @diff by one.
 *
 * Returns: (transfer full): @diff.
 */
NcmDiff *
ncm_diff_ref (NcmDiff *diff)
{
  return g_object_ref (diff);
}

/**
 * ncm_diff_free:
 * @diff: a #NcmDiff
 *
 * Decrease the reference count of @diff by one.
 *
 */
void
ncm_diff_free (NcmDiff *diff)
{
  g_object_unref (diff);
}

/**
 * ncm_diff_clear:
 * @diff: a #NcmDiff
 *
 * Decrease the reference count of @diff by one, and sets the pointer *@diff to
 * NULL.
 *
 */
void
ncm_diff_clear (NcmDiff **diff)
{
  g_clear_object (diff);
}

/**
 * ncm_diff_get_max_order:
 * @diff: a #NcmDiff
 *
 * Gets the maximum order used when calculating the derivatives.
 * 
 * Returns: the maximum order.
 */
guint 
ncm_diff_get_max_order (NcmDiff *diff)
{
  return diff->priv->maxorder;
}

/**
 * ncm_diff_get_richardson_step:
 * @diff: a #NcmDiff
 *
 * Gets the current Richardson step used in the tables.
 * 
 * Returns: the maximum order.
 */
gdouble 
ncm_diff_get_richardson_step (NcmDiff *diff)
{
  return diff->priv->rs;
}

/**
 * ncm_diff_get_round_off_pad:
 * @diff: a #NcmDiff
 *
 * Gets the current round-off padding used in calculations.
 * 
 * Returns: the round-off padding.
 */
gdouble 
ncm_diff_get_round_off_pad (NcmDiff *diff)
{
  return diff->priv->roff_pad;
}

/**
 * ncm_diff_get_ini_h:
 * @diff: a #NcmDiff
 *
 * Gets the current initial step used in calculations.
 * 
 * Returns: the initial step.
 */
gdouble 
ncm_diff_get_ini_h (NcmDiff *diff)
{
  return diff->priv->ini_h;
}

/**
 * ncm_diff_set_max_order:
 * @diff: a #NcmDiff
 * @maxorder: the new maximum order
 *
 * Sets the maximum order used when calculating the derivatives to @maxorder.
 * 
 */
void 
ncm_diff_set_max_order (NcmDiff *diff, const guint maxorder)
{
  g_assert_cmpuint (maxorder, >, 0);

  if (maxorder != diff->priv->maxorder)
  {
    diff->priv->maxorder = maxorder;
    if (diff->priv->central_tables->len > 0)
    {
      _ncm_diff_build_diff_tables (diff);
    }
  }  
}

/**
 * ncm_diff_set_richardson_step:
 * @diff: a #NcmDiff
 * @rs: the new Richardson step
 *
 * Sets the Richardson step used in the tables.
 * 
 */
void 
ncm_diff_set_richardson_step (NcmDiff *diff, const gdouble rs)
{
  g_assert_cmpfloat (rs, >, 1.1);

  if (rs != diff->priv->rs)
  {
    diff->priv->rs = rs;
    if (diff->priv->central_tables->len > 0)
    {
      _ncm_diff_build_diff_tables (diff);
    }
  }  
}

/**
 * ncm_diff_set_round_off_pad:
 * @diff: a #NcmDiff
 * @roff_pad: the new round-off padding
 *
 * Sets the round-off padding used in the calculations.
 * 
 */
void 
ncm_diff_set_round_off_pad (NcmDiff *diff, const gdouble roff_pad)
{
  g_assert_cmpfloat (roff_pad, >, 1.1);
  diff->priv->roff_pad = roff_pad;
}

/**
 * ncm_diff_set_ini_h:
 * @diff: a #NcmDiff
 * @ini_h: the new initial step
 *
 * Sets the initial step used in the calculations.
 * 
 */
void 
ncm_diff_set_ini_h (NcmDiff *diff, const gdouble ini_h)
{
  g_assert_cmpfloat (ini_h, >, GSL_DBL_EPSILON);
  diff->priv->ini_h = ini_h;
}

/**
 * ncm_diff_log_central_tables:
 * @diff: a #NcmDiff
 * 
 * Logs all central tables.
 * 
 */
void
ncm_diff_log_central_tables (NcmDiff *diff)
{
  guint i;

  for (i = 0; i < diff->priv->central_tables->len; i++)
  {
    NcmDiffTable *dtable = g_ptr_array_index (diff->priv->central_tables, i);

    ncm_message ("# NcmDiff[central]  order: %u\n", i + 1);
    ncm_vector_log_vals (dtable->h,      "# NcmDiff[central]  h:      ", "% 22.15g", TRUE);
    ncm_vector_log_vals (dtable->lambda, "# NcmDiff[central]  lambda: ", "% 22.15g", TRUE);
  }
}

/**
 * ncm_diff_log_forward_tables:
 * @diff: a #NcmDiff
 * 
 * Logs all central tables.
 * 
 */
void
ncm_diff_log_forward_tables (NcmDiff *diff)
{
  guint i;

  for (i = 0; i < diff->priv->forward_tables->len; i++)
  {
    NcmDiffTable *dtable = g_ptr_array_index (diff->priv->forward_tables, i);

    ncm_message ("# NcmDiff[forward]  order: %u\n", i + 1);
    ncm_vector_log_vals (dtable->h,      "# NcmDiff[forward]  h:      ", "% 22.15g", TRUE);
    ncm_vector_log_vals (dtable->lambda, "# NcmDiff[forward]  lambda: ", "% 22.15g", TRUE);
  }
}

/**
 * ncm_diff_log_backward_tables:
 * @diff: a #NcmDiff
 * 
 * Logs all central tables.
 * 
 */
void
ncm_diff_log_backward_tables (NcmDiff *diff)
{
  guint i;

  for (i = 0; i < diff->priv->backward_tables->len; i++)
  {
    NcmDiffTable *dtable = g_ptr_array_index (diff->priv->backward_tables, i);

    ncm_message ("# NcmDiff[backward] order: %u\n", i + 1);
    ncm_vector_log_vals (dtable->h,      "# NcmDiff[backward] h:      ", "% 22.15g", TRUE);
    ncm_vector_log_vals (dtable->lambda, "# NcmDiff[backward] lambda: ", "% 22.15g", TRUE);
  }
}

typedef void (*NcmDiffStepAlgo) (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff);
typedef void (*NcmDiffHessianStepAlgo) (NcmDiff *diff, NcmDiffFuncNto1 f, gpointer user_data, const guint a, const gdouble x, const gdouble hx, const guint b, const gdouble y, const gdouble hy, NcmVector *x_v, const gdouble fval, gdouble *df, gdouble *roff);

static void
_ncm_diff_rf_d1_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  /*ncm_vector_log_vals (x_v,  "x_v  ", "% 22.15g", TRUE);*/
  ncm_vector_addto (x_v, a, h);

  /*ncm_vector_log_vals (x_v,  "xh_v ", "% 22.15g", TRUE);*/
  /*ncm_vector_log_vals (f_v,  "f_v  ", "% 22.15g", TRUE);*/
  f (x_v, df, user_data);

  /*ncm_vector_log_vals (yh_v, "yh_v ", "% 22.15g", TRUE);*/

  ncm_vector_memcpy (roff, df);
  ncm_vector_sub_round_off (roff, f_v);

  ncm_vector_sub   (df, f_v);
  ncm_vector_scale (df, 1.0 / h);

  ncm_vector_mul (roff, df);

  NCM_UNUSED (yh1_v);
  NCM_UNUSED (yh2_v);

  /*
   printf ("# t = %u\n", t);
   ncm_vector_log_vals (df,   "df   ", "% 22.15g", TRUE);
   ncm_vector_log_vals (roff, "roff ", "% 22.15g", TRUE);
   */
}

static void
_ncm_diff_rc_d1_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  ncm_vector_addto (x_v, a, h);

  f (x_v, df, user_data);

  ncm_vector_set (x_v, a, x - h);

  f (x_v, yh1_v, user_data);
  
  ncm_vector_memcpy (roff, df);
  ncm_vector_sub_round_off (roff, yh1_v);

  ncm_vector_sub   (df, yh1_v);
  ncm_vector_scale (df, 0.5 / h);

  ncm_vector_mul (roff, df);

  NCM_UNUSED (yh2_v);
}

static void
_ncm_diff_rc_d2_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  ncm_vector_addto (x_v, a, h);

  f (x_v, df, user_data);

  ncm_vector_set (x_v, a, x - h);

  f (x_v, yh1_v, user_data);

  ncm_vector_add (df, yh1_v);
  ncm_vector_scale (df, 0.5);
  
  ncm_vector_memcpy (roff, df);
  ncm_vector_sub_round_off (roff, f_v);

  ncm_vector_sub   (df, f_v);
  ncm_vector_scale (df, 2.0 / (h * h));

  ncm_vector_mul (roff, df);

  NCM_UNUSED (yh2_v);
}

void 
_ncm_diff_rf_Hessian_step (NcmDiff *diff, NcmDiffFuncNto1 f, gpointer user_data, const guint a, const gdouble x, const gdouble hx, const guint b, const gdouble y, const gdouble hy, NcmVector *x_v, const gdouble fval, gdouble *df, gdouble *roff)
{
  gdouble f_hx, f_hy, f_hxhy;
    
  ncm_vector_addto (x_v, a, hx);
  f_hx = f (x_v, user_data);

  ncm_vector_set (x_v, a, x);
  ncm_vector_addto (x_v, b, hy);
  f_hy = f (x_v, user_data);

  ncm_vector_addto (x_v, a, hx);
  f_hxhy = f (x_v, user_data);

  df[0] = ((fval + f_hxhy) - (f_hx + f_hy)) / (hx * hy);

  {
    const gdouble s1 = fval + f_hxhy;
    const gdouble s2 = f_hx + f_hy;

    const gdouble d1 = s1 - s2;

    if (G_UNLIKELY (d1 == 0.0))
      roff[0] = 1.0;
    else
    {
      const gdouble abs_s1  = fabs (s1);
      const gdouble abs_s2  = fabs (s2);
      const gdouble max_s12 = GSL_MAX (abs_s1, abs_s2);
      roff[0] = fabs (max_s12 * GSL_DBL_EPSILON / d1);
    }
  }
}

static GArray *
ncm_diff_by_step_algo (NcmDiff *diff, NcmDiffStepAlgo step_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  GPtrArray *tables  = NULL;
  GPtrArray *dfs     = g_ptr_array_new ();
  GPtrArray *roffs   = g_ptr_array_new ();
  GArray *f_a        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *yh_a       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x_v     = NULL;
  NcmVector *f_v     = NULL; 
  NcmVector *yh1_v   = NULL;
  NcmVector *yh2_v   = NULL;
  NcmVector *dfb     = NULL;
  NcmVector *dfr     = NULL;
  NcmVector *roffb   = NULL;
  NcmVector *roffr   = NULL;
  NcmVector *err     = NULL;
  NcmVector *err_err = NULL;
  NcmVector *ferr    = NULL;
  NcmVector *df_best = NULL;
  GArray *df         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *not_conv   = g_array_new (FALSE, FALSE, sizeof (guchar));
  const guint nvar   = x_a->len;
  NcmMatrix *Eerr_m  = NULL;
  NcmMatrix *df_m;
  guint a;

  g_array_set_size (df, dim * nvar);
  df_m = ncm_matrix_new_array (df, dim);

  if (po == 0)
    tables = diff->priv->forward_tables;
  else
    tables = diff->priv->central_tables;

  if (Eerr != NULL)
  {
    *Eerr = g_array_new (FALSE, FALSE, sizeof (gdouble));
    g_array_set_size (*Eerr, dim * nvar);
    Eerr_m = ncm_matrix_new_array (*Eerr, dim);
  }

  g_ptr_array_set_free_func (dfs,   (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (roffs, (GDestroyNotify) ncm_vector_free);
  
  g_array_set_size (f_a,      dim);
  g_array_set_size (yh_a,     dim);
  g_array_set_size (not_conv, dim);

  x_v   = ncm_vector_new_array (x_a);
  f_v   = ncm_vector_new_array (f_a);
  yh1_v = ncm_vector_new_array (yh_a);
  yh2_v = ncm_vector_new_array (yh_a);

  dfb     = ncm_vector_new (dim);
  dfr     = ncm_vector_new (dim);

  roffb   = ncm_vector_new (dim);
  roffr   = ncm_vector_new (dim);

  err     = ncm_vector_new (dim);
  err_err = ncm_vector_new (dim);
  ferr    = ncm_vector_new (dim);
  df_best = ncm_vector_new (dim);
  
  g_array_unref (f_a);
  g_array_unref (yh_a);

  f (x_v, f_v, user_data);
  
  for (a = 0; a < x_a->len; a++)
  {
    const gdouble x       = g_array_index (x_a, gdouble, a);
    const gdouble scale   = (x == 0.0) ? 1.0 : fabs (x);
    const gdouble h0      = diff->priv->ini_h * scale;
    const guint ntry_conv = 3;
    NcmDiffTable *ldtable = NULL;
    guint order_index;
    guint t = 0;

    ncm_vector_set_all (ferr, GSL_POSINF);
    g_ptr_array_set_size (dfs, 0);
    g_ptr_array_set_size (roffs, 0);

    memset (not_conv->data, ntry_conv, not_conv->len);
        
    for (order_index = 0; order_index < diff->priv->maxorder; order_index++)
    {
      const guint nt       = order_index + 2;
      NcmDiffTable *dtable = g_ptr_array_index (tables, order_index);
      guint i;

      for (; t < nt; t++)
      {
        const gdouble ho      = h0 * ((po == 0) ? ncm_vector_get (dtable->h, t) : sqrt (ncm_vector_get (dtable->h, t)));
        volatile gdouble temp = x + ho;
        const gdouble h       = temp - x;

        NcmVector *df_t   = ncm_vector_new (dim);
        NcmVector *roff_t = ncm_vector_new (dim);

        step_algo (diff, f, user_data, a, x, h, x_v, f_v, yh1_v, yh2_v, df_t, roff_t);

        g_ptr_array_add (dfs,   df_t);
        g_ptr_array_add (roffs, roff_t);
        
        ncm_vector_set (x_v, a, x);
      }

      if (ldtable == NULL)
      {
        ncm_vector_memcpy (dfb, g_ptr_array_index (dfs, 0));
        ncm_vector_memcpy (dfr, g_ptr_array_index (dfs, 0));

        ncm_vector_memcpy (roffb, g_ptr_array_index (roffs, 0));
        ncm_vector_memcpy (roffr, g_ptr_array_index (roffs, 0));
        
        ncm_vector_memcpy (df_best, dfb);
        
        ncm_vector_scale (dfr, ncm_vector_get (dtable->lambda, 0));
        ncm_vector_axpy  (dfr, ncm_vector_get (dtable->lambda, 1), g_ptr_array_index (dfs, 1));

        ncm_vector_scale (roffr, ncm_vector_get (dtable->lambda, 0));
        ncm_vector_axpy  (roffr, ncm_vector_get (dtable->lambda, 1), g_ptr_array_index (roffs, 1));
      }
      else
      {
        ncm_vector_memcpy (dfr, g_ptr_array_index (dfs, 0));
        ncm_vector_scale (dfr, ncm_vector_get (dtable->lambda, 0));

        ncm_vector_memcpy (roffr, g_ptr_array_index (roffs, 0));
        ncm_vector_scale (roffr, ncm_vector_get (dtable->lambda, 0));

        for (i = 1; i < nt; i++)
        {
          ncm_vector_axpy (dfr,   ncm_vector_get (dtable->lambda,  i), g_ptr_array_index (dfs, i));
          ncm_vector_axpy (roffr, ncm_vector_get (dtable->lambda,  i), g_ptr_array_index (roffs, i));
        }
      }
      
      ncm_vector_memcpy (err, dfr);
      ncm_vector_memcpy (err_err, dfr);

      ncm_vector_sub (err, dfb);
      ncm_vector_cmp (err_err, dfb);
      {
        gboolean improve = FALSE;
        for (i = 0; i < dim; i++)
        {
          const gdouble err_i     = fabs (ncm_vector_get (err, i));
          const gdouble ferr_i    = ncm_vector_get (ferr, i);
          const gdouble roffb_i   = fabs (ncm_vector_get (roffb, i)) * diff->priv->roff_pad;
          const gdouble roffr_i   = fabs (ncm_vector_get (roffr, i)) * diff->priv->roff_pad;
          
          const gdouble terr_i    = GSL_MAX (err_i, GSL_MAX (roffb_i, roffr_i));
          const gdouble err_err_i = ncm_vector_get (err_err, i);

          gdouble df_best_i = ncm_vector_get (df_best, i);
          gdouble cerr_i    = ferr_i;

#define NOT_CONV (g_array_index (not_conv, guchar, i))

          if (NOT_CONV && (err_err_i < 1.0e-3))
            NOT_CONV--;
          else
          {
            if (err_err_i > 1.0e-3)
            {
              NOT_CONV = ntry_conv;
              if (err_err_i > 1.0)
                ncm_vector_set (ferr, i, GSL_POSINF);
            }
          }
          
          if ((terr_i < ferr_i) && !NOT_CONV)
          {
            df_best_i = ncm_vector_get (dfr, i);
            ncm_vector_set (df_best, i, df_best_i);
            ncm_vector_set (ferr, i, terr_i);
            cerr_i = terr_i;

            improve = TRUE;
          }

          if (NOT_CONV || (roffr_i < cerr_i))
            improve = TRUE;
/*          
          printf ("[%3u, %3u, %3u] !conv %u % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g improve: %s\n", 
                  nt, i, a, NOT_CONV, ferr_i, err_i, roffb_i, roffr_i,
                  err_err_i, improve ? "T" : "F");
*/
        }

        if (!improve)
          break;
      }
/*
      ncm_vector_log_vals (dfb,     "dfb  ", "% 22.15g", TRUE);
      ncm_vector_log_vals (dfr,     "dfr  ", "% 22.15g", TRUE);
      ncm_vector_log_vals (err,     "err  ", "% 22.15e", TRUE);

      ncm_vector_log_vals (df_best, "df   ", "% 22.15g", TRUE);
      ncm_vector_log_vals (ferr,    "ferr ", "% 22.15e", TRUE);
*/

      ncm_vector_memcpy (dfb, dfr);
      ncm_vector_memcpy (roffb, roffr);
      ldtable = dtable;
    }

    {
      NcmVector *df_a = ncm_matrix_get_row (df_m, a);
      ncm_vector_memcpy (df_a, df_best);
      ncm_vector_free (df_a);

      if (Eerr_m != NULL)
      {
        NcmVector *Eerr_a = ncm_matrix_get_row (Eerr_m, a);
        ncm_vector_memcpy (Eerr_a, ferr);
        ncm_vector_free (Eerr_a);
      }
    }
  }

  if (Eerr_m != NULL)
  {
    ncm_matrix_scale (Eerr_m, NCM_DIFF_ERR_PAD);
  }

  {
    g_array_unref (not_conv);
    
    g_ptr_array_unref (dfs);
    g_ptr_array_unref (roffs);

    ncm_vector_clear (&x_v);
    ncm_vector_clear (&f_v);
    ncm_vector_clear (&yh1_v);
    ncm_vector_clear (&yh2_v);

    ncm_vector_clear (&dfb);
    ncm_vector_clear (&dfr);

    ncm_vector_clear (&roffb);
    ncm_vector_clear (&roffr);

    ncm_vector_clear (&err);
    ncm_vector_clear (&err_err);
    
    ncm_vector_clear (&ferr);
    ncm_vector_clear (&df_best);

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

static GArray *
ncm_diff_Hessian_by_step_algo (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, guint po, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  GPtrArray *tables = NULL;
  GArray *dfs       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *roffs     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x_v    = NULL;
  GArray *df        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint nvar  = x_a->len;
  gdouble fval      = 0.0;
  NcmMatrix *Eerr_m = NULL;
  NcmMatrix *df_m;
  guint a;

  g_array_set_size (df, nvar * nvar);
  df_m = ncm_matrix_new_array (df, nvar); 

  if (po == 0)
    tables = diff->priv->forward_tables;
  else
    tables = diff->priv->central_tables;

  if (Eerr != NULL)
  {
    *Eerr = g_array_new (FALSE, FALSE, sizeof (gdouble));
    g_array_set_size (*Eerr, nvar * nvar);
    Eerr_m = ncm_matrix_new_array (*Eerr, nvar);
  }

  x_v  = ncm_vector_new_array (x_a);

  fval = f (x_v, user_data);

  for (a = 0; a < x_a->len; a++)
  {
    guint b;
    for (b = a + 1; b < x_a->len; b++)
    {
      const guint ntry_conv = 3;
      const gdouble x       = g_array_index (x_a, gdouble, a);
      const gdouble y       = g_array_index (x_a, gdouble, b);
      const gdouble scale_x = (x == 0.0) ? 1.0 : fabs (x);
      const gdouble scale_y = (y == 0.0) ? 1.0 : fabs (y);
      const gdouble hx0     = diff->priv->ini_h * scale_x;
      const gdouble hy0     = diff->priv->ini_h * scale_y;
      NcmDiffTable *ldtable = NULL;
      gdouble ferr          = GSL_POSINF;
      gdouble err           = 0.0;
      gdouble err_err       = 0.0;
      gdouble df_best       = 0.0;
      gdouble dfb           = 0.0;
      gdouble dfr           = 0.0;
      gdouble roffb         = 0.0;
      gdouble roffr         = 0.0;        
      guint t               = 0;
      guint not_converging  = ntry_conv;
      guint order_index;

      g_array_set_size (dfs, 0);
      g_array_set_size (roffs, 0);
      
      for (order_index = 0; order_index < diff->priv->maxorder; order_index++)
      {
        const guint nt        = order_index + 2;
        NcmDiffTable *dtable  = g_ptr_array_index (tables, order_index);
        guint i;

        for (; t < nt; t++)
        {
          const gdouble hxo    = hx0 * ((po == 0) ? ncm_vector_get (dtable->h, t) : sqrt (ncm_vector_get (dtable->h, t)));
          const gdouble hyo    = hy0 * ((po == 0) ? ncm_vector_get (dtable->h, t) : sqrt (ncm_vector_get (dtable->h, t)));
          volatile gdouble t_x = x + hxo;
          const gdouble hx     = t_x - x;
          volatile gdouble t_y = y + hyo;
          const gdouble hy     = t_y - y;

          gdouble df_t   = 0.0;
          gdouble roff_t = 0.0;

          Hstep_algo (diff, f, user_data, a, x, hx, b, y, hy, x_v, fval, &df_t, &roff_t);

          g_array_append_val (dfs,   df_t);
          g_array_append_val (roffs, roff_t);

          ncm_vector_set (x_v, a, x);
          ncm_vector_set (x_v, b, y);
        }

        if (ldtable == NULL)
        {
          const gdouble lambda0 = ncm_vector_get (dtable->lambda, 0);
          const gdouble lambda1 = ncm_vector_get (dtable->lambda, 1);
          
          dfb   = g_array_index (dfs, gdouble, 0);
          dfr   = dfb * lambda0 + g_array_index (dfs, gdouble, 1) * lambda1;

          roffb = g_array_index (roffs, gdouble, 0);
          roffr = roffb * lambda0 + g_array_index (roffs, gdouble, 1) * lambda1;

          df_best = dfb;
        }
        else
        {
          dfr   = 0.0;
          roffr = 0.0;

          for (i = 0; i < nt; i++)
          {
            const gdouble lambda_i = ncm_vector_get (dtable->lambda,  i);
            const gdouble df_i     = g_array_index (dfs, gdouble, i);
            const gdouble roff_i   = g_array_index (roffs, gdouble, i);
            
            dfr   += lambda_i * df_i;
            roffr += lambda_i * roff_i;
          }
        }

        err     = fabs (dfr - dfb);
        err_err = (dfr == 0.0) ? ((dfb == 0.0) ? 0.0 : fabs (dfb)) : ((dfb == 0.0) ? fabs (dfr) : fabs ((dfr - dfb) / GSL_MIN (fabs (dfr), fabs (dfb))));

        {
          gboolean improve = FALSE;

          const gdouble Eroffb = fabs (roffb) * diff->priv->roff_pad;
          const gdouble Eroffr = fabs (roffr) * diff->priv->roff_pad;
          const gdouble terr   = GSL_MAX (err, GSL_MAX (Eroffb, Eroffr));
          gdouble cerr         = ferr;

          if (not_converging && (err_err < 1.0e-3))
            not_converging--;
          else 
          {
            if (err_err > 1.0e-3)
            {
              not_converging = ntry_conv;
              if (err_err > 1.0)
                ferr = GSL_POSINF;
            }
          }
          
/*          
          printf ("[%3u, %3u, %3u] !conv %u % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g", 
                  nt, a, b, not_converging, ferr, err, Eroffb, Eroffr, terr, err_err);
*/
          if ((terr < ferr) && !not_converging)
          {
            df_best = dfr;
            ferr    = terr;
            cerr    = terr;
            improve = TRUE;

            /*printf (" -UP-");*/
          }
/*
          else
            printf ("     ");
*/
          if (not_converging || (Eroffr < cerr))
            improve = TRUE;

          /*printf (" improve: %s\n", improve ? "T" : "F");*/

          if (!improve)
            break;
        }
        
        dfb     = dfr;
        roffb   = roffr;
        ldtable = dtable;
      }

      ncm_matrix_set (df_m, a, b, df_best);
      ncm_matrix_set (df_m, b, a, df_best);
      
      if (Eerr_m != NULL)
      {
        ncm_matrix_set (Eerr_m, a, b, ferr);
        ncm_matrix_set (Eerr_m, b, a, ferr);
      }
    }
  }

  if (Eerr_m != NULL)
  {
    ncm_matrix_scale (Eerr_m, NCM_DIFF_ERR_PAD);
  }
  
  {
    g_array_unref (dfs);
    g_array_unref (roffs);

    ncm_vector_clear (&x_v);

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

/**
 * ncm_diff_rf_d1_N_to_M:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the first derivative of @f: $\partial_i f$ using the forward method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N\to \mathbb{R}^M$,
 * where $N = $ length of @x_a and $M = $ @dim.
 * 
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rf_d1_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  return ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, 0, x_a, dim, f, user_data, Eerr);
}

/**
 * ncm_diff_rc_d1_N_to_M:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the first derivative of @f: $\partial_i f$ using the central method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N\to \mathbb{R}^M$,
 * where $N = $ length of @x_a and $M = $ @dim.
 * 
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d1_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, 1, x_a, dim, f, user_data, Eerr);
}

/**
 * ncm_diff_rc_d2_N_to_M:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the second derivative of @f: $\partial_i^2 f$ using the central method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N\to \mathbb{R}^M$,
 * where $N = $ length of @x_a and $M = $ @dim.
 * 
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d2_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, dim, f, user_data, Eerr);
}

typedef struct _NcmDiffFuncParams
{
  NcmDiffFunc1toM f_1_to_M;
  NcmDiffFuncNto1 f_N_to_1;
  NcmDiffFunc1to1 f_1_to_1;
  gpointer user_data;
} NcmDiffFuncParams;

static void 
_ncm_diff_trans_1_to_M (NcmVector *x, NcmVector *y, gpointer user_data)
{
  NcmDiffFuncParams *fp = (NcmDiffFuncParams *) user_data;
  fp->f_1_to_M (ncm_vector_get (x, 0), y, fp->user_data);
}

static void 
_ncm_diff_trans_N_to_1 (NcmVector *x, NcmVector *y, gpointer user_data)
{
  NcmDiffFuncParams *fp = (NcmDiffFuncParams *) user_data;
  ncm_vector_set (y, 0, fp->f_N_to_1 (x, fp->user_data));
}

static void 
_ncm_diff_trans_1_to_1 (NcmVector *x, NcmVector *y, gpointer user_data)
{
  NcmDiffFuncParams *fp = (NcmDiffFuncParams *) user_data;
  ncm_vector_set (y, 0, fp->f_1_to_1 (ncm_vector_get (x, 0), fp->user_data));
}

/**
 * ncm_diff_rf_d1_1_to_M:
 * @diff: a #NcmDiff
 * @x: function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the first derivative of @f: $\partial_i f$ using the forward method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 * 
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rf_d1_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {f, NULL, NULL, user_data};
  GArray *x_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *df_a;
  
  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;
  
  df_a =  ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, 0, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);

  g_array_unref (x_a);

  return df_a;
}

/**
 * ncm_diff_rc_d1_1_to_M:
 * @diff: a #NcmDiff
 * @x: function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the first derivative of @f: $\partial_i f$ using the central method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 * 
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d1_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {f, NULL, NULL, user_data};
  GArray *x_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *df_a;
  
  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, 1, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);
  g_array_unref (x_a);

  return df_a;
}

/**
 * ncm_diff_rc_d2_1_to_M:
 * @diff: a #NcmDiff
 * @x: function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the second derivative of @f: $\partial_i^2 f$ using the central method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 * 
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d2_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {f, NULL, NULL, user_data};
  GArray *x_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *df_a;
  
  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;
  
  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);
  g_array_unref (x_a);

  return df_a;
}

/**
 * ncm_diff_rf_d1_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the first derivative of @f: $\partial_i f$ using the forward method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 * 
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rf_d1_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};
  
  return ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, 0, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
}

/**
 * ncm_diff_rc_d1_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the first derivative of @f: $\partial_i f$ using the central method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 * 
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d1_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};

  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
}

/**
 * ncm_diff_rc_d2_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the second derivative of @f: $\partial_i^2 f$ using the central method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 * 
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d2_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};
  
  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
}

/**
 * ncm_diff_rf_Hessian_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 * 
 * Calculates the Hessian of @f $\partial_i\partial_j f$ using the forward method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 * 
 * Returns: (transfer full) (array) (element-type double): The Hessian of @f at @x_a.
 */
GArray *
ncm_diff_rf_Hessian_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};
  GArray *dEerr = NULL;

  GArray *diag = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, &dEerr);
  GArray *res  = ncm_diff_Hessian_by_step_algo (diff, _ncm_diff_rf_Hessian_step, 0, x_a, f, user_data, Eerr);

  guint i;

  g_assert_cmpuint (diag->len * diag->len, ==, res->len);
  
  for (i = 0; i < diag->len; i++)
  {
    g_array_index (res, gdouble, i * diag->len + i) = g_array_index (diag, gdouble, i);
    if (Eerr != NULL)
    {
      g_array_index (*Eerr, gdouble, i * diag->len + i) = g_array_index (dEerr, gdouble, i);
    }
  }

  g_array_unref (dEerr);
  g_array_unref (diag);

  return res;
}

/**
 * ncm_diff_rf_d1_1_to_1:
 * @diff: a #NcmDiff
 * @x: function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @err: (out) (nullable): estimated error
 * 
 * Calculates the first derivative of @f: $\partial_i f$ using the forward method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R} \to \mathbb{R}$.
 * 
 * Returns: The derivative of @f at @x.
 */
gdouble
ncm_diff_rf_d1_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err)
{
  NcmDiffFuncParams fp = {NULL, NULL, f, user_data};
  GArray *x_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *Eerr = NULL;
  GArray *df_a;
  gdouble df;
  
  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;
  
  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, 0, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

  df = g_array_index (df_a, gdouble, 0);

  g_array_unref (x_a);
  g_array_unref (df_a);

  if (err != NULL)
    *err = g_array_index (Eerr, gdouble, 0);
  g_array_unref (Eerr);

  return df;
}

/**
 * ncm_diff_rc_d1_1_to_1:
 * @diff: a #NcmDiff
 * @x: function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @err: (out) (nullable): estimated error
 * 
 * Calculates the first derivative of @f: $\partial_i f$ using the central method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R} \to \mathbb{R}$.
 * 
 * Returns: The derivative of @f at @x.
 */
gdouble
ncm_diff_rc_d1_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err)
{
  NcmDiffFuncParams fp = {NULL, NULL, f, user_data};
  GArray *x_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *Eerr = NULL;
  GArray *df_a;
  gdouble df;
  
  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;
  
  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, 1, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

  df = g_array_index (df_a, gdouble, 0);

  g_array_unref (x_a);
  g_array_unref (df_a);

  if (err != NULL)
    *err = g_array_index (Eerr, gdouble, 0);
  g_array_unref (Eerr);

  return df;
}

/**
 * ncm_diff_rc_d2_1_to_1:
 * @diff: a #NcmDiff
 * @x: function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @err: (out) (nullable): estimated error
 * 
 * Calculates the second derivative of @f: $\partial_i^2 f$ using the central method plus 
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R} \to \mathbb{R}$.
 * 
 * Returns: The derivative of @f at @x.
 */
gdouble
ncm_diff_rc_d2_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err)
{
  NcmDiffFuncParams fp = {NULL, NULL, f, user_data};
  GArray *x_a  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *Eerr = NULL;
  GArray *df_a;
  gdouble df;
  
  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;
  
  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

  df = g_array_index (df_a, gdouble, 0);

  g_array_unref (x_a);
  g_array_unref (df_a);

  if (err != NULL)
    *err = g_array_index (Eerr, gdouble, 0);
  g_array_unref (Eerr);

  return df;
}
