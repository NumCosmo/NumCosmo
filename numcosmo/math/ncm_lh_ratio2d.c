/***************************************************************************
 *            ncm_lh_ratio2d.c
 *
 *  Mon Jun 11 13:28:00 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_lh_ratio2d
 * @title: Likelihood Ratio 2D
 * @short_description: Likelihood ratio object for bidimensional analysis.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm_lh_ratio2d.h"

#include "math/ncm_c.h"
#include "math/ncm_cfg.h"
#include "math/ncm_matrix.h"
#include "math/util.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_eigen.h>

enum
{
  PROP_0,
  PROP_FIT,
  PROP_PI1,
  PROP_PI2,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmLHRatio2d, ncm_lh_ratio2d, G_TYPE_OBJECT);
G_DEFINE_BOXED_TYPE (NcmLHRatio2dRegion, ncm_lh_ratio2d_region, ncm_lh_ratio2d_region_dup, ncm_lh_ratio2d_region_free);

static void
ncm_lh_ratio2d_init (NcmLHRatio2d *lhr2d)
{
  lhr2d->fit         = NULL;
  lhr2d->constrained = NULL;
  lhr2d->rng         = ncm_rng_new (NULL);
  lhr2d->pi[0].mid   = -1;
  lhr2d->pi[0].pid   = 0;
  lhr2d->pi[1].mid   = -1;
  lhr2d->pi[1].pid   = 0;
  lhr2d->chisquare   = 0.0;
  lhr2d->mtype       = NCM_FIT_RUN_MSGS_NONE;
  lhr2d->rtype       = NCM_LH_RATIO2D_ROOT_BRACKET;
  lhr2d->e_vec       = ncm_matrix_new (2, 2);
  lhr2d->e_val       = ncm_vector_new (2);
  lhr2d->r        = 0.0;
  lhr2d->theta    = 0.0;
  lhr2d->shift[0] = 0.0;
  lhr2d->shift[1] = 0.0;
  lhr2d->angular  = FALSE;
}

static void _ncm_lh_ratio2d_prepare_coords (NcmLHRatio2d *lhr2d);

static void
ncm_lh_ratio2d_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_lh_ratio2d_parent_class)->constructed (object);
  {
    NcmLHRatio2d *lhr2d = NCM_LH_RATIO2D (object);
    NcmMSet *mset = ncm_mset_dup (lhr2d->fit->mset);
    gint i;
    g_assert (lhr2d->fit->fstate->is_best_fit);

    for (i = 0; i < 2; i++)
    {
      g_assert_cmpint (lhr2d->pi[i].mid, >=, 0);
      g_assert_cmpint (lhr2d->pi[i].mid,  <, NCM_MODEL_MAX_ID);

      if (ncm_mset_peek (lhr2d->fit->mset, lhr2d->pi[i].mid) == NULL)
        g_error ("ncm_lh_ratio2d_constructed: cannot use parameter[%d:%u], model not set\n", 
                 lhr2d->pi[i].mid, lhr2d->pi[i].pid);

      if (ncm_mset_param_get_ftype (lhr2d->fit->mset, lhr2d->pi[0].mid, lhr2d->pi[0].pid) != NCM_PARAM_TYPE_FREE)
        g_error ("ncm_lh_ratio2d_constructed: cannot find for a non fitted parameter[%d:%u].", 
                 lhr2d->pi[i].mid, lhr2d->pi[i].pid);

      ncm_mset_param_set_ftype (mset,
                                lhr2d->pi[i].mid, lhr2d->pi[i].pid, 
                                NCM_PARAM_TYPE_FIXED);

      lhr2d->lb[i] = ncm_mset_param_get_lower_bound (lhr2d->fit->mset, lhr2d->pi[i].mid, lhr2d->pi[i].pid);
      lhr2d->ub[i] = ncm_mset_param_get_upper_bound (lhr2d->fit->mset, lhr2d->pi[i].mid, lhr2d->pi[i].pid);
      lhr2d->bf[i] = ncm_mset_param_get (lhr2d->fit->mset, lhr2d->pi[i].mid, lhr2d->pi[i].pid);
    }
    
    lhr2d->constrained = ncm_fit_copy_new (lhr2d->fit, lhr2d->fit->lh, mset,
                                           lhr2d->fit->grad.gtype);
    ncm_mset_free (mset);
    _ncm_lh_ratio2d_prepare_coords (lhr2d);
  }
}

static void
ncm_lh_ratio2d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmLHRatio2d *lhr2d = NCM_LH_RATIO2D (object);
  g_return_if_fail (NCM_IS_LH_RATIO2D (object));

  switch (prop_id)
  {
    case PROP_FIT:
      lhr2d->fit = g_value_dup_object (value);
      break;
    case PROP_PI1:
    {
      NcmMSetPIndex *pi = g_value_get_boxed (value);
      lhr2d->pi[0] = *pi;
      break;
    }
    case PROP_PI2:
    {
      NcmMSetPIndex *pi = g_value_get_boxed (value);
      lhr2d->pi[1] = *pi;
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_lh_ratio2d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmLHRatio2d *lhr2d = NCM_LH_RATIO2D (object);
  g_return_if_fail (NCM_IS_LH_RATIO2D (object));

  switch (prop_id)
  {
    case PROP_PI1:
    {
      NcmMSetPIndex *pi = ncm_mset_pindex_new (lhr2d->pi[0].mid, lhr2d->pi[0].pid);
      g_value_take_boxed (value, pi);
      break;
    }
    case PROP_PI2:
    {
      NcmMSetPIndex *pi = ncm_mset_pindex_new (lhr2d->pi[1].mid, lhr2d->pi[1].pid);
      g_value_take_boxed (value, pi);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_lh_ratio2d_dispose (GObject *object)
{
  NcmLHRatio2d *lhr2d = NCM_LH_RATIO2D (object);

  ncm_fit_clear (&lhr2d->fit);
  ncm_fit_clear (&lhr2d->constrained);
  ncm_matrix_clear (&lhr2d->e_vec);
  ncm_vector_clear (&lhr2d->e_val);
  ncm_rng_clear (&lhr2d->rng);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_lh_ratio2d_parent_class)->dispose (object);
}

static void
ncm_lh_ratio2d_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_lh_ratio2d_parent_class)->finalize (object);
}

static void
ncm_lh_ratio2d_class_init (NcmLHRatio2dClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &ncm_lh_ratio2d_constructed;
  object_class->set_property = &ncm_lh_ratio2d_set_property;
  object_class->get_property = &ncm_lh_ratio2d_get_property;
  object_class->dispose      = &ncm_lh_ratio2d_dispose;
  object_class->finalize     = &ncm_lh_ratio2d_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_object ("fit",
                                                        NULL,
                                                        "NcmFit object",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_WRITABLE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_PI1,
                                   g_param_spec_boxed ("pi1",
                                                       NULL,
                                                       "First param index",
                                                       NCM_TYPE_MSET_PINDEX,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_PI2,
                                   g_param_spec_boxed ("pi2",
                                                       NULL,
                                                       "Second param index",
                                                       NCM_TYPE_MSET_PINDEX,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


/**
 * ncm_lh_ratio2d_new:
 * @fit: FIXME
 * @pi1: FIXME
 * @pi2: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmLHRatio2d *
ncm_lh_ratio2d_new (NcmFit *fit, NcmMSetPIndex *pi1, NcmMSetPIndex *pi2)
{
  return g_object_new (NCM_TYPE_LH_RATIO2D, 
                       "fit", fit,
                       "pi1", pi1,
                       "pi2", pi2,
                       NULL);
  
}

/**
 * ncm_lh_ratio2d_free:
 * @lhr2d: FIXME
 *
 * FIXME
 *
 */
void 
ncm_lh_ratio2d_free (NcmLHRatio2d *lhr2d)
{
  g_object_unref (lhr2d);
}

/**
 * ncm_lh_ratio2d_clear:
 * @lhr2d: FIXME
 *
 * FIXME
 *
 */
void 
ncm_lh_ratio2d_clear (NcmLHRatio2d **lhr2d)
{
  g_clear_object (lhr2d);
}

static void
_ncm_lh_ratio2d_prepare_coords (NcmLHRatio2d *lhr2d)
{
  const guint fpi1 = ncm_mset_fparam_get_fpi (lhr2d->fit->mset, lhr2d->pi[0].mid, lhr2d->pi[0].pid);
  const guint fpi2 = ncm_mset_fparam_get_fpi (lhr2d->fit->mset, lhr2d->pi[1].mid, lhr2d->pi[1].pid);
  NcmMatrix *cov = ncm_matrix_new (2, 2);
  gsl_eigen_symmv_workspace *vw = gsl_eigen_symmv_alloc (2);

  g_assert_cmpint (fpi1, >=, 0);
  g_assert_cmpint (fpi2, >=, 0);
  g_assert_cmpint (fpi1, !=, fpi2);

  ncm_matrix_set (cov, 0, 0, ncm_fit_covar_fparam_cov (lhr2d->fit, fpi1, fpi1));
  ncm_matrix_set (cov, 0, 1, ncm_fit_covar_fparam_cov (lhr2d->fit, fpi1, fpi2));
  ncm_matrix_set (cov, 1, 0, ncm_fit_covar_fparam_cov (lhr2d->fit, fpi2, fpi1));
  ncm_matrix_set (cov, 1, 1, ncm_fit_covar_fparam_cov (lhr2d->fit, fpi2, fpi2));

  gsl_eigen_symmv (NCM_MATRIX_GSL (cov), 
                   ncm_vector_gsl (lhr2d->e_val), 
                   NCM_MATRIX_GSL (lhr2d->e_vec), 
                   vw);
  gsl_eigen_symmv_free (vw);

  {
    gdouble sigma_e[2];
    gdouble e_vec1[2];
    gdouble e_vec2[2];
    gint ei1 = 0, ei2 = 1;
    if (ncm_vector_get (lhr2d->e_val, 0) < ncm_vector_get (lhr2d->e_val, 1))
    { 
      ei1 = 1; 
      ei2 = 0;
      gsl_vector_swap_elements (ncm_vector_gsl (lhr2d->e_val), 0, 1);
    }

    sigma_e[0] = sqrt (ncm_vector_get (lhr2d->e_val, 0));
    sigma_e[1] = sqrt (ncm_vector_get (lhr2d->e_val, 1));

    e_vec1[0] = ncm_matrix_get (lhr2d->e_vec, 0, ei1);
    e_vec1[1] = ncm_matrix_get (lhr2d->e_vec, 1, ei1);

    e_vec2[0] = ncm_matrix_get (lhr2d->e_vec, 0, ei2);
    e_vec2[1] = ncm_matrix_get (lhr2d->e_vec, 1, ei2);      

    ncm_matrix_set (lhr2d->e_vec, 0, 0, e_vec1[0] * sigma_e[0]);
    ncm_matrix_set (lhr2d->e_vec, 1, 0, e_vec1[1] * sigma_e[0]);

    ncm_matrix_set (lhr2d->e_vec, 0, 1, e_vec2[0] * sigma_e[1]);
    ncm_matrix_set (lhr2d->e_vec, 1, 1, e_vec2[1] * sigma_e[1]);    
  }

  ncm_matrix_free (cov);
}

/**
 * ncm_lh_ratio2d_set_pindex:
 * @lhr2d: a #NcmLHRatio2d.
 * @pi1: FIXME
 * @pi2: FIXME
 *
 * FIXME
 *
 */
void 
ncm_lh_ratio2d_set_pindex (NcmLHRatio2d *lhr2d, NcmMSetPIndex *pi1, NcmMSetPIndex *pi2)
{
  gint i;
  NcmMSetPIndex pi[2] = {*pi1, *pi2};
  for (i = 0; i < 2; i++)
  {
    if (ncm_mset_param_get_ftype (lhr2d->fit->mset, pi[i].mid, pi[i].pid) != NCM_PARAM_TYPE_FREE)
      g_error ("ncm_lh_ratio2d_set_pindex: cannot find bounds for a non fitted parameter[%d:%u].", 
               pi[i].mid, pi[i].pid);

    ncm_mset_param_set_ftype (lhr2d->constrained->mset, 
                              lhr2d->pi[i].mid, lhr2d->pi[i].pid, 
                              NCM_PARAM_TYPE_FREE);
    lhr2d->pi[i] = pi[i];
    ncm_mset_param_set_ftype (lhr2d->constrained->mset, 
                              lhr2d->pi[i].mid, lhr2d->pi[i].pid, 
                              NCM_PARAM_TYPE_FIXED);

    lhr2d->lb[i] = ncm_mset_param_get_lower_bound (lhr2d->fit->mset, pi[i].mid, pi[i].pid);
    lhr2d->ub[i] = ncm_mset_param_get_upper_bound (lhr2d->fit->mset, pi[i].mid, pi[i].pid);
    lhr2d->bf[i] = ncm_mset_param_get (lhr2d->fit->mset, pi[i].mid, pi[i].pid);
  }
  _ncm_lh_ratio2d_prepare_coords (lhr2d);
}

static gboolean _ncm_lh_ratio2d_log_dot = FALSE;

static void
ncm_lh_ratio2d_log_start (NcmLHRatio2d *lhr2d, gdouble clevel)
{
  if (lhr2d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    ncm_cfg_msg_sepa ();
    g_message ("# Likelihood ratio confidence region at %2.3f%%, bestfit [% 12.8g % 12.8g]:\n", 
               100.0 * clevel, lhr2d->bf[0], lhr2d->bf[1]);
    if (lhr2d->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("#");
      _ncm_lh_ratio2d_log_dot = TRUE;
    }
  }
}

static void
ncm_lh_ratio2d_log_param_val (NcmLHRatio2d *lhr2d, gdouble r, gdouble val)
{
  if (lhr2d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (lhr2d->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
    {
      if (!_ncm_lh_ratio2d_log_dot)
      {
        g_message ("#");
        _ncm_lh_ratio2d_log_dot = TRUE;
      }
      g_message (".");
    }
    else
      g_message ("#  parameter % 12.8g likelihood ratio % 12.8g.\n", r, val);
  }
}

static void
ncm_lh_ratio2d_log_root_start (NcmLHRatio2d *lhr2d, gdouble pl, gdouble pu)
{
  if (lhr2d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (_ncm_lh_ratio2d_log_dot)
    {
      g_message ("\n");
      _ncm_lh_ratio2d_log_dot = FALSE;
    }
    g_message ("#  looking root in interval [% 12.8g % 12.8g]:\n", pl, pu);
    if (lhr2d->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("#");
      _ncm_lh_ratio2d_log_dot = TRUE;
    }
  }
}

static void
ncm_lh_ratio2d_log_root_step (NcmLHRatio2d *lhr2d, gdouble pl, gdouble pu)
{
  if (lhr2d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (lhr2d->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
      g_message (".");
    else
      g_message ("#  parameter root bounds [% 12.8g % 12.8g].\n", pl, pu);
  }
}

static void
ncm_lh_ratio2d_log_root_finish (NcmLHRatio2d *lhr2d, gdouble pr, gdouble prec)
{
  if (lhr2d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (_ncm_lh_ratio2d_log_dot)
    {
      g_message ("\n");
      _ncm_lh_ratio2d_log_dot = FALSE;
    }
    g_message ("#  root found at % 12.8g with precision %1.8e.\n", pr, prec);
  }  
}

static void
ncm_lh_ratio2d_log_border_found (NcmLHRatio2d *lhr2d, gdouble r)
{
  if (lhr2d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (_ncm_lh_ratio2d_log_dot)
    {
      g_message ("\n");
      _ncm_lh_ratio2d_log_dot = FALSE;
    }
    g_message ("#  border found at % 12.8g.\n", r);
  }  
}

static gboolean 
_ncm_lh_ratio2d_inside_interval (gdouble *p, const gdouble lb, const gdouble ub, const gdouble prec)
{
  if (p[0] < lb)
  {
    if (ncm_cmp (p[0], lb, 1e-4) == 0)
      p[0] = lb;
    else
      return FALSE;
  }
  else if (p[0] > ub)
  {
    if (ncm_cmp (p[0], ub, 1e-4) == 0)
      p[0] = ub;
    else
      return FALSE;
  }
  return TRUE;
}

static void
ncm_lh_ratio2d_tofparam (NcmLHRatio2d *lhr2d, gdouble *p1, gdouble *p2)
{
  gdouble alpha, beta;

  alpha = lhr2d->shift[0] + lhr2d->r * cos (lhr2d->theta);
  beta  = lhr2d->shift[1] + lhr2d->r * sin (lhr2d->theta);

  *p1 = lhr2d->bf[0] + alpha * ncm_matrix_get (lhr2d->e_vec, 0, 0) + beta * ncm_matrix_get (lhr2d->e_vec, 0, 1);
  *p2 = lhr2d->bf[1] + alpha * ncm_matrix_get (lhr2d->e_vec, 1, 0) + beta * ncm_matrix_get (lhr2d->e_vec, 1, 1);
}

static gdouble
ncm_lh_ratio2d_f (gdouble x, gpointer ptr)
{
  NcmLHRatio2d *lhr2d = NCM_LH_RATIO2D (ptr);
  gdouble p[2];
  gboolean skip = FALSE;
  if (lhr2d->angular)
    lhr2d->theta = x;
  else
    lhr2d->r = x;

  ncm_lh_ratio2d_tofparam (lhr2d, &p[0], &p[1]);
  
  skip = skip || !_ncm_lh_ratio2d_inside_interval (&p[0], lhr2d->lb[0], lhr2d->ub[0], 1e-4);
  skip = skip || !_ncm_lh_ratio2d_inside_interval (&p[1], lhr2d->lb[1], lhr2d->ub[1], 1e-4);

  if (skip)
    return HUGE_VAL;
  
  ncm_mset_param_set_pi (lhr2d->constrained->mset, lhr2d->pi, p, 2);
  ncm_fit_run (lhr2d->constrained, NCM_FIT_RUN_MSGS_NONE);

  lhr2d->niter     += lhr2d->constrained->fstate->niter;
  lhr2d->func_eval += lhr2d->constrained->fstate->func_eval;
  lhr2d->grad_eval += lhr2d->constrained->fstate->grad_eval;

  {
    const gdouble m2lnL_const = ncm_fit_state_get_m2lnL_curval (lhr2d->constrained->fstate);
    const gdouble m2lnL = ncm_fit_state_get_m2lnL_curval (lhr2d->fit->fstate);
    return m2lnL_const - (m2lnL + lhr2d->chisquare);
  }
}

static gdouble
ncm_lh_ratio2d_root_brent (NcmLHRatio2d *lhr2d, gdouble x0, gdouble x)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gdouble prec = 1e-5, x1 = x;

  F.function = &ncm_lh_ratio2d_f;
  F.params   = lhr2d;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x0, x1);

  ncm_lh_ratio2d_log_root_start (lhr2d, x0, x);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    if (status)
    {
      g_warning ("%s", gsl_strerror (status));
      gsl_root_fsolver_free (s);
      return GSL_NAN;
    }

    x = gsl_root_fsolver_root (s);
    x0 = gsl_root_fsolver_x_lower (s);
    x1 = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x0, x1, 0, prec);

    ncm_lh_ratio2d_log_root_step (lhr2d, x0, x1);

    if (!gsl_finite (ncm_lh_ratio2d_f (x, lhr2d)))
    {
      g_debug ("Ops");
      x = GSL_NAN;
      break;
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);
  ncm_lh_ratio2d_log_root_finish (lhr2d, x, prec);

  return x;
}

static gdouble
ncm_lh_ratio2d_numdiff_df (gdouble x, gpointer p)
{
  gsl_function F;
  gdouble res, err;
  F.function = &ncm_lh_ratio2d_f;
  F.params = p;
  res = ncm_numdiff_1 (&F, x, x * 1e-5, &err);

  return res;
}

static void
ncm_lh_ratio2d_numdiff_fdf (gdouble x, gpointer p, gdouble *y, gdouble *dy)
{
  *dy = ncm_lh_ratio2d_numdiff_df (x, p);
  *y = ncm_lh_ratio2d_f (x, p);
  return;
}

static gdouble
ncm_lh_ratio2d_root_steffenson (NcmLHRatio2d *lhr2d, gdouble x0, gdouble x1)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  gsl_function_fdf F;
  gdouble prec = 1e-5;
  gdouble x = (x0 + x1) * 0.5;

  F.f = &ncm_lh_ratio2d_f;
  F.df = &ncm_lh_ratio2d_numdiff_df;
  F.fdf = &ncm_lh_ratio2d_numdiff_fdf;
  F.params = lhr2d;

  T = gsl_root_fdfsolver_steffenson;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &F, x);

  ncm_lh_ratio2d_log_root_start (lhr2d, x0, x);
  
  do
  {
    iter++;
    status = gsl_root_fdfsolver_iterate (s);

    if (status)
    {
      g_warning ("%s", gsl_strerror (status));
      gsl_root_fdfsolver_free (s);
      return GSL_NAN;
    }

    x0 = x;
    x = gsl_root_fdfsolver_root (s);
    status = gsl_root_test_delta (x, x0, 0, prec);

    ncm_lh_ratio2d_log_root_step (lhr2d, x, x0);

    if (!gsl_finite (ncm_lh_ratio2d_f (x, lhr2d)))
    {
      g_debug ("Ops");
      x = GSL_NAN;
      break;
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  ncm_lh_ratio2d_log_root_finish (lhr2d, x, prec);
    
  gsl_root_fdfsolver_free (s);
  return x;
}

/**
 * ncm_lh_ratio2d_points_add:
 * @points: FIXME
 * @x: FIXME
 * @y: FIXME
 * @theta: FIXME
 * @p1: FIXME
 * @p2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
ncm_lh_ratio2d_points_add (NcmLHRatio2d *lhr2d, GList **points)
{
  NcmLHRatio2dPoint *p_new = g_slice_new (NcmLHRatio2dPoint);
  p_new->x = lhr2d->shift[0] + lhr2d->r * cos (lhr2d->theta);
  p_new->y = lhr2d->shift[1] + lhr2d->r * sin (lhr2d->theta);
  p_new->theta = lhr2d->theta;

  ncm_lh_ratio2d_tofparam (lhr2d, &p_new->p1, &p_new->p2);
  
  *points = g_list_append (*points, p_new);
}

/**
 * ncm_lh_ratio2d_points_exists: (skip)
 * @points: FIXME
 * @x: FIXME
 * @y: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_lh_ratio2d_points_exists (GList **points, gdouble x, gdouble y, guint n)
{
  GList *pos = g_list_last (*points);
  while (n != 0)
  {
    NcmLHRatio2dPoint *crp1, *crp2;
    gdouble t1, t2;
    if (pos == NULL) break;
    crp1 = (NcmLHRatio2dPoint *) pos->data;
    pos = g_list_previous (pos);
    if (pos == NULL) break;
    crp2 = (NcmLHRatio2dPoint *) pos->data;
    t1 = (x - crp1->x) / (crp2->x - crp1->x);
    t2 = (y - crp1->y) / (crp2->y - crp1->y);
    if (fabs (t1 - t2) < 1e-1 && t1 >= 0.0 && t1 <= 1.0)
      return TRUE;
    n--;
  }
  return FALSE;
}

/**
 * ncm_lh_ratio2d_points_print: (skip)
 * @points: FIXME
 * @out: FIXME
 *
 * FIXME
 *
 */
void
ncm_lh_ratio2d_points_print (GList *points, FILE *out)
{
  GList *spoints = points;
  points = g_list_first (spoints);
  while (points)
  {
    NcmLHRatio2dPoint *crp = (NcmLHRatio2dPoint *) points->data;
    fprintf (out, "% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", crp->p1, crp->p2, crp->x, crp->y, crp->theta);fflush(out);
    points = g_list_next (points);
  }
  points = g_list_first (spoints);
  {
    NcmLHRatio2dPoint *crp = (NcmLHRatio2dPoint *) points->data;
    fprintf (out, "% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", crp->p1, crp->p2, crp->x, crp->y, crp->theta);fflush(out);
  }
}

/**
 * ncm_lh_ratio2d_points_free: (skip)
 * @points: FIXME
 *
 * FIXME
 *
 */
void
ncm_lh_ratio2d_points_free (GList *points)
{
  points = g_list_first (points);
  while (points)
  {
    NcmLHRatio2dPoint *crp = (NcmLHRatio2dPoint *) points->data;
    g_slice_free (NcmLHRatio2dPoint, crp);
    points = g_list_next (points);
  }
  g_list_free (points);
}

/**
 * ncm_lh_ratio2d_points_to_region: (skip)
 * @points: FIXME
 * @clevel: FIXME
 *
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmLHRatio2dRegion *
ncm_lh_ratio2d_points_to_region (GList *points, gdouble clevel)
{
  GList *spoints = points;
  NcmLHRatio2dRegion *rg = g_slice_new0 (NcmLHRatio2dRegion);
  guint i = 0;
  
  points = g_list_first (spoints);
  rg->np = g_list_length (points) + 1;
  rg->p1 = ncm_vector_new (rg->np);
  rg->p2 = ncm_vector_new (rg->np);
  rg->clevel = clevel;
  
  while (points)
  {
    NcmLHRatio2dPoint *crp = (NcmLHRatio2dPoint *) points->data;
    ncm_vector_set (rg->p1, i, crp->p1);
    ncm_vector_set (rg->p2, i, crp->p2);
    i++;
    points = g_list_next (points);
  }
  points = g_list_first (spoints);
  {
    NcmLHRatio2dPoint *crp = (NcmLHRatio2dPoint *) points->data;
    ncm_vector_set (rg->p1, i, crp->p1);
    ncm_vector_set (rg->p2, i, crp->p2);
    i++;
  }
  g_assert_cmpint (i, ==, rg->np);

  return rg;
}

#define TIMEOUT 90.0
#define RESCALE (0.5)
#define NMAXTRIES 40

static void
_ncm_lh_ratio2d_set_angular_interval (NcmLHRatio2d *lhr2d, const gdouble da, gdouble *theta0, gdouble *theta1)
{
  gdouble val0 = ncm_lh_ratio2d_f (*theta0, lhr2d);
  gdouble val1 = ncm_lh_ratio2d_f (*theta1, lhr2d);
  g_assert (*theta1 > *theta0);
  g_assert (da > 0.0);

  if (val0 < 0.0)
  {
    do
    {
      ncm_lh_ratio2d_log_param_val (lhr2d, ncm_c_radian_to_degree (ncm_c_radian_0_2pi (*theta0)), val0);
      *theta1 = *theta0;
      *theta0 -= da;
    } while ((val0 = ncm_lh_ratio2d_f (*theta0, lhr2d)) < 0);
  }
  else if (val1 > 0.0)
  {
    do
    {
      ncm_lh_ratio2d_log_param_val (lhr2d, ncm_c_radian_to_degree (ncm_c_radian_0_2pi (*theta1)), val1);
      *theta0 = *theta1;
      *theta1 += da;
    } while ((val1 = ncm_lh_ratio2d_f (*theta1, lhr2d)) > 0);
  }
}

/**
 * ncm_lh_ratio2d_conf_region:
 * @lhr2d: a @NcmFit
 * @clevel: FIXME
 * @expected_np: Expected number of points, if lesser than 1 it uses the default value of 100.
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmLHRatio2dRegion *
ncm_lh_ratio2d_conf_region (NcmLHRatio2d *lhr2d, gdouble clevel, gdouble expected_np, NcmFitRunMsgs mtype)
{
  GTimer *iter_timer = g_timer_new ();
  gdouble total_time = 0.0;
  gdouble init_x, init_y;
  gint i, counter = -1;
  GList *points = NULL, *final_points = NULL;
  gboolean completed = FALSE;
  gdouble second_try = FALSE;
  static gdouble (*root) (NcmLHRatio2d *, gdouble, gdouble);
  if (expected_np <= 1.0)
    expected_np = 100.0;
  
  lhr2d->mtype = mtype;
  
  ncm_lh_ratio2d_log_start (lhr2d, clevel);

  lhr2d->chisquare = gsl_cdf_chisq_Qinv (1.0 - clevel, 2);
  lhr2d->shift[0] = 0.0;
  lhr2d->shift[1] = 0.0;
  lhr2d->theta = gsl_rng_uniform (lhr2d->rng->r) * 2.0 * M_PI;
  
  switch (lhr2d->rtype)
  {
    case NCM_LH_RATIO2D_ROOT_BRACKET:
      root = ncm_lh_ratio2d_root_brent;
      break;
    case NCM_LH_RATIO2D_ROOT_NUMDIFF:
      root = ncm_lh_ratio2d_root_steffenson;
      break;
    default:
      g_assert_not_reached ();
      break;
  }
  
  while (!completed)
  {
    gdouble val, r0, r, scale;
    gdouble theta0, theta1;

    lhr2d->angular = FALSE;

    scale = sqrt (lhr2d->chisquare);
    r0 = 0.0;
    r = scale;

    while ((val = ncm_lh_ratio2d_f (r, lhr2d)) < 0)
    {
      ncm_lh_ratio2d_log_param_val (lhr2d, r, val);
      r0 = r;
      r += scale;
    }
    r = root (lhr2d, r0, r);
    ncm_lh_ratio2d_log_border_found (lhr2d, r);

    completed = TRUE;
    if (second_try)
      completed = TRUE;

    ncm_lh_ratio2d_points_add (lhr2d, &points);
    
    lhr2d->shift[0] = lhr2d->r * cos (lhr2d->theta);
    lhr2d->shift[1] = lhr2d->r * sin (lhr2d->theta);

    init_x = lhr2d->shift[0];
    init_y = lhr2d->shift[1];

    lhr2d->angular = TRUE;
    lhr2d->r = (2.0 * M_PI * sqrt (lhr2d->chisquare) / expected_np);

    theta0 = lhr2d->theta - M_PI * 0.25;
    theta1 = lhr2d->theta + M_PI * 0.25;

    i = 0;
    total_time += g_timer_elapsed (iter_timer, NULL);
    g_timer_start (iter_timer);

    while (TRUE)
    {
      _ncm_lh_ratio2d_set_angular_interval (lhr2d, M_PI * 0.5, &theta0, &theta1);
      lhr2d->theta = root (lhr2d, theta0, theta1);

      theta0 = lhr2d->theta - M_PI * 0.25;
      theta1 = lhr2d->theta + M_PI * 0.25;

      lhr2d->theta = ncm_c_radian_0_2pi (lhr2d->theta);

      ncm_lh_ratio2d_points_add (lhr2d, &points);

      lhr2d->shift[0] += lhr2d->r * cos (lhr2d->theta);
      lhr2d->shift[1] += lhr2d->r * sin (lhr2d->theta);

      counter--;
      if (!counter)
        break;

      {
        gboolean near_x = fabs (init_x - lhr2d->shift[0]) < lhr2d->r;
        gboolean near_y = fabs (init_y - lhr2d->shift[1]) < lhr2d->r;
        if (i > 10 && near_x && near_y && counter < 0)
        {
          g_message ("#  Start found at [%d], ending...\n",i);
          completed = TRUE;
          break;
          counter = 5;
        }
      }
      total_time += g_timer_elapsed (iter_timer, NULL);
      g_timer_start (iter_timer);
      i++;
    }
    
    ncm_lh_ratio2d_points_add (lhr2d, &points);
    if (final_points == NULL)
    {
      final_points = points;
      points = NULL;
    }
    else
    {
      final_points = g_list_concat (final_points, g_list_reverse (points));
      points = NULL;
    }
    if (!completed)
    {
      second_try = TRUE;
      g_message ("#  Trying in another direction\n");
    }
  }

  total_time += g_timer_elapsed (iter_timer, NULL);
  g_timer_destroy (iter_timer);
  
  {
    NcmLHRatio2dRegion *rg = ncm_lh_ratio2d_points_to_region (final_points, clevel);
    ncm_lh_ratio2d_points_free (final_points);
    return rg;
  }
}

/**
 * ncm_lh_ratio2d_fisher_border:
 * @lhr2d: a #NcmFit.
 * @clevel: FIXME
 * @expected_np:  Expected number of points, if lesser than 1 it uses the default value of 600.
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmLHRatio2dRegion *
ncm_lh_ratio2d_fisher_border (NcmLHRatio2d *lhr2d, gdouble clevel, gdouble expected_np, NcmFitRunMsgs mtype)
{
  GList *points = NULL;
  gdouble theta;
  gdouble step;
  if (expected_np <= 1.0)
    expected_np = 600.0;

  step = 2.0 * M_PI / expected_np;

  lhr2d->chisquare = gsl_cdf_chisq_Qinv (1.0 - clevel, 2);
  lhr2d->r = sqrt (lhr2d->chisquare);

  lhr2d->shift[0] = 0.0;
  lhr2d->shift[1] = 0.0;

  for (theta = 0.0; theta <= 2.0 * M_PI; theta += step)
  {
    lhr2d->theta = theta;
    ncm_lh_ratio2d_points_add (lhr2d, &points);
  }

  {
    NcmLHRatio2dRegion *rg = ncm_lh_ratio2d_points_to_region (points, clevel);
    ncm_lh_ratio2d_points_free (points);
    return rg;
  }
}

/**
 * ncm_lh_ratio2d_region_dup:
 * @rg: a #NcmLHRatio2dRegion.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmLHRatio2dRegion *
ncm_lh_ratio2d_region_dup (NcmLHRatio2dRegion *rg)
{
  NcmLHRatio2dRegion *rg_dup = g_slice_new0 (NcmLHRatio2dRegion);
  rg_dup->np = rg->np;
  rg_dup->clevel = rg->clevel;
  rg_dup->p1 = ncm_vector_dup (rg->p1);
  rg_dup->p2 = ncm_vector_dup (rg->p2);
  
  return rg_dup;
}

/**
 * ncm_lh_ratio2d_region_free:
 * @rg: a #NcmLHRatio2dRegion.
 *
 * FIXME
 *
 */
void 
ncm_lh_ratio2d_region_free (NcmLHRatio2dRegion *rg)
{
  ncm_vector_clear (&rg->p1);
  ncm_vector_clear (&rg->p2);
  g_slice_free (NcmLHRatio2dRegion, rg);
}

/**
 * ncm_lh_ratio2d_region_clear:
 * @rg: a #NcmLHRatio2dRegion.
 *
 * FIXME
 *
 */
void 
ncm_lh_ratio2d_region_clear (NcmLHRatio2dRegion **rg)
{
  g_clear_pointer (rg, &ncm_lh_ratio2d_region_free);
}

/**
 * ncm_lh_ratio2d_region_print:
 * @rg: FIXME
 * @out: FIXME
 *
 * FIXME
 *
 */
void
ncm_lh_ratio2d_region_print (NcmLHRatio2dRegion *rg, FILE *out)
{
  guint i;
  for (i = 0; i < rg->np; i++)
  {
    const gdouble p1 = ncm_vector_get (rg->p1, i);
    const gdouble p2 = ncm_vector_get (rg->p2, i);
    fprintf (out, "% -19.15g % -19.15g\n", p1, p2);fflush(out);
  }
}

