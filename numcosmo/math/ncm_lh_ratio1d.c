/***************************************************************************
 *            ncm_lh_ratio1d.c
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
 * SECTION:ncm_lh_ratio1d
 * @title: NcmLHRatio1d
 * @short_description: Likelihood ratio for one dimensional parameter analysis.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_lh_ratio1d.h"

#include "math/ncm_c.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>

enum
{
  PROP_0,
  PROP_FIT,
  PROP_PI,
  PROP_CONSTRAINT,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmLHRatio1d, ncm_lh_ratio1d, G_TYPE_OBJECT);

static void
ncm_lh_ratio1d_init (NcmLHRatio1d *lhr1d)
{
  lhr1d->fit         = NULL;
  lhr1d->constrained = NULL;
  lhr1d->pi.mid      = -1;
  lhr1d->pi.pid      = 0;
  lhr1d->constraint  = NULL;
  lhr1d->chisquare   = 0.0;
  lhr1d->mtype       = NCM_FIT_RUN_MSGS_NONE;
  lhr1d->rtype       = NCM_LH_RATIO1D_ROOT_BRACKET;
}

static void
ncm_lh_ratio1d_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_lh_ratio1d_parent_class)->constructed (object);
  {
    NcmLHRatio1d *lhr1d = NCM_LH_RATIO1D (object);
    NcmSerialize *ser = ncm_serialize_global ();
    NcmMSet *mset = ncm_mset_dup (lhr1d->fit->mset, ser);

    ncm_serialize_free (ser);
    g_assert_cmpint (lhr1d->pi.mid, >=, 0);
    
    g_assert (lhr1d->fit->fstate->is_best_fit);

    if (ncm_mset_peek (lhr1d->fit->mset, lhr1d->pi.mid) == NULL)
      g_error ("ncm_lh_ratio1d_constructed: cannot use parameter[%d:%u], model not set.", 
               lhr1d->pi.mid, lhr1d->pi.pid);
    
    if (ncm_mset_param_get_ftype (lhr1d->fit->mset, lhr1d->pi.mid, lhr1d->pi.pid) != NCM_PARAM_TYPE_FREE)
      g_error ("ncm_lh_ratio1d_constructed: cannot find for a non fitted parameter[%d:%u].", 
               lhr1d->pi.mid, lhr1d->pi.pid);

    ncm_mset_param_set_ftype (mset,
                              lhr1d->pi.mid, lhr1d->pi.pid, 
                              NCM_PARAM_TYPE_FIXED);

    lhr1d->constrained = ncm_fit_copy_new (lhr1d->fit, lhr1d->fit->lh, mset,
                                           lhr1d->fit->grad.gtype);
    ncm_mset_free (mset);

    lhr1d->lb = ncm_mset_param_get_lower_bound (lhr1d->fit->mset, lhr1d->pi.mid, lhr1d->pi.pid);
    lhr1d->ub = ncm_mset_param_get_upper_bound (lhr1d->fit->mset, lhr1d->pi.mid, lhr1d->pi.pid);
    lhr1d->bf = ncm_mset_param_get (lhr1d->fit->mset, lhr1d->pi.mid, lhr1d->pi.pid);
  }
}

static void
ncm_lh_ratio1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmLHRatio1d *lhr1d = NCM_LH_RATIO1D (object);
  g_return_if_fail (NCM_IS_LH_RATIO1D (object));

  switch (prop_id)
  {
    case PROP_FIT:
      lhr1d->fit = g_value_dup_object (value);
      break;
    case PROP_PI:
    {
      NcmMSetPIndex *pi = g_value_get_boxed (value);
      lhr1d->pi = *pi;
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_lh_ratio1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmLHRatio1d *lhr1d = NCM_LH_RATIO1D (object);
  g_return_if_fail (NCM_IS_LH_RATIO1D (object));

  switch (prop_id)
  {
    case PROP_PI:
    {
      NcmMSetPIndex *pi = ncm_mset_pindex_new (lhr1d->pi.mid, lhr1d->pi.pid);
      g_value_take_boxed (value, pi);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_lh_ratio1d_dispose (GObject *object)
{
  NcmLHRatio1d *lhr1d = NCM_LH_RATIO1D (object);

  ncm_fit_clear (&lhr1d->fit);
  ncm_fit_clear (&lhr1d->constrained);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_lh_ratio1d_parent_class)->dispose (object);
}

static void
ncm_lh_ratio1d_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_lh_ratio1d_parent_class)->finalize (object);
}

static void
ncm_lh_ratio1d_class_init (NcmLHRatio1dClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &ncm_lh_ratio1d_constructed;
  object_class->set_property = &ncm_lh_ratio1d_set_property;
  object_class->get_property = &ncm_lh_ratio1d_get_property;
  object_class->dispose      = &ncm_lh_ratio1d_dispose;
  object_class->finalize     = &ncm_lh_ratio1d_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_object ("fit",
                                                        NULL,
                                                        "NcmFit object",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_WRITABLE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_PI,
                                   g_param_spec_boxed ("pi",
                                                       NULL,
                                                       "Param index",
                                                       NCM_TYPE_MSET_PINDEX,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_PI,
                                   g_param_spec_object ("constraint",
                                                        NULL,
                                                        "Constraint",
                                                        NCM_TYPE_MSET_FUNC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


/**
 * ncm_lh_ratio1d_new:
 * @fit: FIXME
 * @pi: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmLHRatio1d *
ncm_lh_ratio1d_new (NcmFit *fit, NcmMSetPIndex *pi)
{
  return g_object_new (NCM_TYPE_LH_RATIO1D, 
                       "fit", fit,
                       "pi", pi,
                       NULL);
  
}

/**
 * ncm_lh_ratio1d_free:
 * @lhr1d: FIXME
 *
 * FIXME
 *
 */
void 
ncm_lh_ratio1d_free (NcmLHRatio1d *lhr1d)
{
  g_object_unref (lhr1d);
}

/**
 * ncm_lh_ratio1d_clear:
 * @lhr1d: FIXME
 *
 * FIXME
 *
 */
void 
ncm_lh_ratio1d_clear (NcmLHRatio1d **lhr1d)
{
  g_clear_object (lhr1d);
}

/**
 * ncm_lh_ratio1d_set_pindex:
 * @lhr1d: a #NcmLHRatio1d.
 * @pi: FIXME
 *
 * FIXME
 *
 */
void 
ncm_lh_ratio1d_set_pindex (NcmLHRatio1d *lhr1d, NcmMSetPIndex *pi)
{
  if (ncm_mset_param_get_ftype (lhr1d->fit->mset, pi->mid, pi->pid) != NCM_PARAM_TYPE_FREE)
    g_error ("ncm_lh_ratio1d_set_pindex: cannot find bounds for a non fitted parameter[%d:%u].", 
             pi->mid, pi->pid);

  ncm_mset_param_set_ftype (lhr1d->constrained->mset, 
                            lhr1d->pi.mid, lhr1d->pi.pid, 
                            NCM_PARAM_TYPE_FREE);
  lhr1d->pi = *pi;
  ncm_mset_param_set_ftype (lhr1d->constrained->mset, 
                            lhr1d->pi.mid, lhr1d->pi.pid, 
                            NCM_PARAM_TYPE_FIXED);

  lhr1d->lb = ncm_mset_param_get_lower_bound (lhr1d->fit->mset, pi->mid, pi->pid);
  lhr1d->ub = ncm_mset_param_get_upper_bound (lhr1d->fit->mset, pi->mid, pi->pid);
  lhr1d->bf = ncm_mset_param_get (lhr1d->fit->mset, pi->mid, pi->pid);
}

static gboolean _ncm_lh_ratio1d_log_dot = FALSE;

static void
ncm_lh_ratio1d_log_start (NcmLHRatio1d *lhr1d, gdouble clevel)
{
  if (lhr1d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    ncm_cfg_msg_sepa ();
    g_message ("# Likelihood ratio bounds at %2.3f%%, bestfit % 12.8g:\n", 
               100.0 * clevel, lhr1d->bf);
    if (lhr1d->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("#");
      _ncm_lh_ratio1d_log_dot = TRUE;
    }
  }
}

static void
ncm_lh_ratio1d_log_param_val (NcmLHRatio1d *lhr1d, gdouble p, gdouble val)
{
  if (lhr1d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (lhr1d->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
    {
      if (!_ncm_lh_ratio1d_log_dot)
      {
        g_message ("#");
        _ncm_lh_ratio1d_log_dot = TRUE;
      }
      g_message (".");
    }
    else
      g_message ("#  parameter % 12.8g likelihood ratio % 12.8g.\n", 
                 lhr1d->bf + p, val);
  }
}

static void
ncm_lh_ratio1d_log_root_start (NcmLHRatio1d *lhr1d, gdouble pl, gdouble pu)
{
  if (lhr1d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (_ncm_lh_ratio1d_log_dot)
    {
      g_message ("\n");
      _ncm_lh_ratio1d_log_dot = FALSE;
    }
    g_message ("#  looking root in interval [% 12.8g % 12.8g]:\n", 
               lhr1d->bf + pl, lhr1d->bf + pu);
    if (lhr1d->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("#");
      _ncm_lh_ratio1d_log_dot = TRUE;
    }
  }
}

static void
ncm_lh_ratio1d_log_root_step (NcmLHRatio1d *lhr1d, gdouble pl, gdouble pu)
{
  if (lhr1d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (lhr1d->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
      g_message (".");
    else
      g_message ("#  parameter root bounds [% 12.8g % 12.8g].\n", 
                 lhr1d->bf + pl, lhr1d->bf + pu);
  }
}

static void
ncm_lh_ratio1d_log_root_finish (NcmLHRatio1d *lhr1d, gdouble pr, gdouble prec)
{
  if (lhr1d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (_ncm_lh_ratio1d_log_dot)
    {
      g_message ("\n");
      _ncm_lh_ratio1d_log_dot = FALSE;
    }
    g_message ("#  root found at % 12.8g with precision %1.8e.\n", 
               lhr1d->bf + pr, prec);
  }  
}

static void
ncm_lh_ratio1d_log_finish (NcmLHRatio1d *lhr1d, gdouble pl, gdouble pu)
{
  if (lhr1d->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (_ncm_lh_ratio1d_log_dot)
    {
      g_message ("\n");
      _ncm_lh_ratio1d_log_dot = FALSE;
    }
    g_message ("#  lower and upper bounds found [% 12.8g % 12.8g].\n", 
               lhr1d->bf + pl, lhr1d->bf + pu);
    g_message ("#  iteration            [%06d]\n", lhr1d->niter);
    g_message ("#  function evaluations [%06d]\n", lhr1d->func_eval);
    g_message ("#  gradient evaluations [%06d]\n", lhr1d->grad_eval);
    
  }  
}

static gdouble
ncm_lh_ratio1d_f (gdouble x, gpointer ptr)
{
  NcmLHRatio1d *lhr1d = NCM_LH_RATIO1D (ptr);
  gdouble p = lhr1d->bf + x;

  p = GSL_MAX (p, lhr1d->lb);
  p = GSL_MIN (p, lhr1d->ub);

  ncm_mset_param_set (lhr1d->constrained->mset, lhr1d->pi.mid, lhr1d->pi.pid, p);

  ncm_fit_run (lhr1d->constrained, NCM_FIT_RUN_MSGS_NONE);

  lhr1d->niter     += lhr1d->constrained->fstate->niter;
  lhr1d->func_eval += lhr1d->constrained->fstate->func_eval;
  lhr1d->grad_eval += lhr1d->constrained->fstate->grad_eval;

  if (p == lhr1d->lb)
  {
    g_warning ("reached lower bound stoping...");
    return 0.0;
  }
  if (p == lhr1d->ub)
  {
    g_warning ("reached upper bound stoping...");
    return 0.0;
  }

  {
    const gdouble m2lnL_const = ncm_fit_state_get_m2lnL_curval (lhr1d->constrained->fstate);
    const gdouble m2lnL = ncm_fit_state_get_m2lnL_curval (lhr1d->fit->fstate);
    return m2lnL_const - (m2lnL + lhr1d->chisquare);
  }
}

static gdouble
ncm_lh_ratio1d_root_brent (NcmLHRatio1d *lhr1d, gdouble x0, gdouble x)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gdouble prec = 1e-5, x1 = x;

  F.function = &ncm_lh_ratio1d_f;
  F.params   = lhr1d;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x0, x1);

  ncm_lh_ratio1d_log_root_start (lhr1d, x0, x);

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

    ncm_lh_ratio1d_log_root_step (lhr1d, x0, x1);

    if (!gsl_finite (ncm_lh_ratio1d_f (x, lhr1d)))
    {
      g_debug ("Ops");
      x = GSL_NAN;
      break;
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);
  ncm_lh_ratio1d_log_root_finish (lhr1d, x, prec);

  return x;
}

static gdouble
ncm_lh_ratio1d_numdiff_df (gdouble x, gpointer p)
{
  gsl_function F;
  gdouble res, err;
  F.function = &ncm_lh_ratio1d_f;
  F.params = p;
  res = ncm_numdiff_1 (&F, x, x * 1e-5, &err);

  return res;
}

static void
ncm_lh_ratio1d_numdiff_fdf (gdouble x, gpointer p, gdouble *y, gdouble *dy)
{
  *dy = ncm_lh_ratio1d_numdiff_df (x, p);
  *y = ncm_lh_ratio1d_f (x, p);
  return;
}


static gdouble
ncm_lh_ratio1d_root_steffenson (NcmLHRatio1d *lhr1d, gdouble x0, gdouble x1)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  gsl_function_fdf F;
  gdouble prec = 1e-5;
  gdouble x = (x0 + x1) * 0.5;

  F.f = &ncm_lh_ratio1d_f;
  F.df = &ncm_lh_ratio1d_numdiff_df;
  F.fdf = &ncm_lh_ratio1d_numdiff_fdf;
  F.params = lhr1d;

  T = gsl_root_fdfsolver_steffenson;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &F, x);

  ncm_lh_ratio1d_log_root_start (lhr1d, x0, x);
  
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

    ncm_lh_ratio1d_log_root_step (lhr1d, x, x0);

    if (!gsl_finite (ncm_lh_ratio1d_f (x, lhr1d)))
    {
      g_debug ("Ops");
      x = GSL_NAN;
      break;
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  ncm_lh_ratio1d_log_root_finish (lhr1d, x, prec);
    
  gsl_root_fdfsolver_free (s);
  return x;
}

#define NCM_LH_RATIO1D_SCALE_INCR (1.1)

/**
 * ncm_lh_ratio1d_find_bounds:
 * @lhr1d: a #NcmLHRatio1d.
 * @clevel: The confidence level (0,1).
 * @mtype: FIXME
 * @lb: (out): Lower bound
 * @ub: (out): Upper bound 
 * 
 * FIXME
 * 
 */
void 
ncm_lh_ratio1d_find_bounds (NcmLHRatio1d *lhr1d, gdouble clevel, NcmFitRunMsgs mtype, gdouble *lb, gdouble *ub)
{
  gdouble scale, r, r_min, r_max, val;
  static gdouble (*root) (NcmLHRatio1d *lhr1d, gdouble x0, gdouble x);

  g_assert (clevel > 0.0 && clevel < 1.0);

  lhr1d->chisquare = gsl_cdf_chisq_Qinv (1.0 - clevel, 1.0);
  scale = sqrt (lhr1d->chisquare) * 
    ncm_fit_covar_sd (lhr1d->fit, lhr1d->pi.mid, lhr1d->pi.pid);
  r = 0.0;
  r_min = -scale;
  r_max =  scale;
    
  switch (lhr1d->rtype)
  {
    case NCM_LH_RATIO1D_ROOT_BRACKET:
      root = ncm_lh_ratio1d_root_brent;
      break;
    case NCM_LH_RATIO1D_ROOT_NUMDIFF:
      root = ncm_lh_ratio1d_root_steffenson;
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  lhr1d->mtype = mtype;

  ncm_lh_ratio1d_log_start (lhr1d, clevel);
  
  while ((val = ncm_lh_ratio1d_f (r_min, lhr1d)) < 0.0)
  {
    ncm_lh_ratio1d_log_param_val (lhr1d, r_min, val);
    r = r_min;
    r_min *= NCM_LH_RATIO1D_SCALE_INCR;
  }
  r_min = root (lhr1d, r_min, r);

  r = 0.0;
  while ((val = ncm_lh_ratio1d_f (r_max, lhr1d)) < 0.0)
  {
    ncm_lh_ratio1d_log_param_val (lhr1d, r_max, val);
    r = r_max;
    r_max *= NCM_LH_RATIO1D_SCALE_INCR;
  }
  r_max = root (lhr1d, r, r_max);

  *lb = r_min;
  *ub = r_max;

  ncm_lh_ratio1d_log_finish (lhr1d, r_min, r_max);
}

