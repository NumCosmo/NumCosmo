/***************************************************************************
 *            ncm_fit_gsl_ls.c
 *
 *  Mon Jun 11 12:04:33 2007
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
 * SECTION:ncm_fit_gsl_ls
 * @title: Least Squares Algorithims -- GSL
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_gsl_ls.h"

#include <gsl/gsl_blas.h>

G_DEFINE_TYPE (NcmFitGSLLS, ncm_fit_gsl_ls, NCM_TYPE_FIT);

static void
ncm_fit_gsl_ls_init (NcmFitGSLLS *fit_gsl_ls)
{
  fit_gsl_ls->ls = NULL; 
}

static gint ncm_fit_gsl_ls_f (const gsl_vector *x, gpointer p, gsl_vector *f);
static gint ncm_fit_gsl_ls_df (const gsl_vector *x, gpointer p, gsl_matrix *J);
static gint ncm_fit_gsl_ls_fdf (const gsl_vector *x, gpointer p, gsl_vector *f, gsl_matrix *J);

static void
_ncm_fit_gsl_ls_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_gsl_ls_parent_class)->constructed (object);
  {
    NcmFitGSLLS *fit_gsl_ls = NCM_FIT_GSL_LS (object);
    NcmFit *fit = NCM_FIT (fit_gsl_ls);
    const gsl_multifit_fdfsolver_type *T;

    fit_gsl_ls->f.f      = &ncm_fit_gsl_ls_f;
    fit_gsl_ls->f.df     = &ncm_fit_gsl_ls_df;
    fit_gsl_ls->f.fdf    = &ncm_fit_gsl_ls_fdf;
    fit_gsl_ls->f.p      = fit->fparam_len;
    fit_gsl_ls->f.n      = fit->data_len;
    fit_gsl_ls->f.params = fit;

    T = gsl_multifit_fdfsolver_lmsder;
    fit_gsl_ls->ls = gsl_multifit_fdfsolver_alloc (T, fit_gsl_ls->f.n, fit_gsl_ls->f.p);
  }
}

static void
ncm_fit_gsl_ls_finalize (GObject *object)
{
  NcmFitGSLLS *fit_gsl_ls = NCM_FIT_GSL_LS (object);

  gsl_multifit_fdfsolver_free (fit_gsl_ls->ls);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_gsl_ls_parent_class)->finalize (object);
}

static NcmFit *_ncm_fit_gsl_ls_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
static gboolean _ncm_fit_gsl_ls_run (NcmFit *fit, NcmFitRunMsgs mtype);
static const gchar *_ncm_fit_gsl_ls_get_desc (NcmFit *fit);

static void
ncm_fit_gsl_ls_class_init (NcmFitGSLLSClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmFitClass* fit_class     = NCM_FIT_CLASS (klass);

  object_class->constructed  = &_ncm_fit_gsl_ls_constructed;
  object_class->finalize     = &ncm_fit_gsl_ls_finalize;

  fit_class->copy_new = &_ncm_fit_gsl_ls_copy_new;
  fit_class->run      = &_ncm_fit_gsl_ls_run;
  fit_class->get_desc = &_ncm_fit_gsl_ls_get_desc;

}

static NcmFit *
_ncm_fit_gsl_ls_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  return ncm_fit_gsl_ls_new (lh, mset, gtype);
}

static gdouble ncm_fit_gsl_ls_test_grad (NcmFit *fit);

#define _NCM_FIT_GSL_LS_MIN_PREC_RETRY (1e-3)

gboolean
_ncm_fit_gsl_ls_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitGSLLS *fit_gsl_ls = NCM_FIT_GSL_LS (fit);
  gint status;
  gdouble prec = 1e-11;

  ncm_mset_fparams_get_vector (fit->mset, fit->x);
  gsl_multifit_fdfsolver_set (fit_gsl_ls->ls, &fit_gsl_ls->f, ncm_vector_gsl (fit->x));
  
  do
  {
    fit->niter++;
    status = gsl_multifit_fdfsolver_iterate (fit_gsl_ls->ls);

    if (fit->niter == 1 && !gsl_finite (gsl_blas_dnrm2(fit_gsl_ls->ls->f)))
    {
      ncm_mset_fparams_set_vector (fit->mset, fit->x);
      return FALSE;
    }
    prec = ncm_fit_gsl_ls_test_grad (fit);

    if (status)
    {
      if (mtype > NCM_FIT_RUN_MSGS_NONE)
        ncm_fit_log_step_error (fit, gsl_strerror (status));
      if (status == GSL_CONTINUE)
        break;
    }

    //status = gsl_multifit_test_gradient (fit->grad, prec);
    status = gsl_multifit_test_delta (fit_gsl_ls->ls->dx, fit_gsl_ls->ls->x, 1e-13, 1e-13);

    {
      gdouble sqrt_m2lnL = gsl_blas_dnrm2 (fit_gsl_ls->ls->f);
      ncm_fit_log_step (fit, sqrt_m2lnL * sqrt_m2lnL);
    }
  }
  while (status == GSL_CONTINUE && fit->niter < fit->maxeval);

  fit->sqrt_m2lnL = gsl_blas_dnrm2 (fit_gsl_ls->ls->f);
  fit->m2lnL = fit->sqrt_m2lnL * fit->sqrt_m2lnL;
  fit->m2lnL_dof = fit->m2lnL / fit->dof;
  fit->m2lnL_prec = prec;

  {
    NcmVector *_x = ncm_vector_new_gsl_static (fit_gsl_ls->ls->x);
    NcmVector *_f = ncm_vector_new_gsl_static (fit_gsl_ls->ls->f);
    NcmMatrix *_J = ncm_matrix_new_gsl_static (fit_gsl_ls->ls->J);
    ncm_fit_save (fit, _x, _f, _J);
    ncm_vector_free (_x);
    ncm_vector_free (_f);
    ncm_matrix_free (_J);
  }

  return TRUE;
}

static gint
ncm_fit_gsl_ls_f (const gsl_vector *x, gpointer p, gsl_vector *f)
{
  NcmFit *fit = NCM_FIT (p);
  NcmVector *fv = ncm_vector_new_gsl_static (f);
  
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  ncm_fit_ls_f (fit, fv);
  
  ncm_vector_free (fv);
  return GSL_SUCCESS;
}

static gint
ncm_fit_gsl_ls_df (const gsl_vector *x, gpointer p, gsl_matrix *J)
{
  NcmFit *fit = NCM_FIT (p);
  NcmMatrix *Jm = ncm_matrix_new_gsl_static (J);
  
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  ncm_fit_ls_J (fit, Jm);

  ncm_matrix_free (Jm);
  return GSL_SUCCESS;
}

static gint
ncm_fit_gsl_ls_fdf (const gsl_vector *x, gpointer p, gsl_vector *f, gsl_matrix *J)
{
  NcmFit *fit = NCM_FIT (p);
  NcmVector *fv = ncm_vector_new_gsl_static (f);
  NcmMatrix *Jm = ncm_matrix_new_gsl_static (J);
  
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  ncm_fit_ls_f_J (fit, fv, Jm);
  
  ncm_vector_free (fv);
  ncm_matrix_free (Jm);
  return GSL_SUCCESS;
}

static const gchar *
_ncm_fit_gsl_ls_get_desc (NcmFit *fit)
{
  static gchar *desc = NULL;
  if (desc == NULL)
  {
    NcmFitGSLLS *fit_gsl_ls = NCM_FIT_GSL_LS (fit);
    desc = g_strdup_printf ("GSL Least Squares:%s", gsl_multifit_fdfsolver_name (fit_gsl_ls->ls));
  }
  return desc;
}

static gdouble
ncm_fit_gsl_ls_test_grad (NcmFit *fit)
{
  NcmFitGSLLS *fit_gsl_ls = NCM_FIT_GSL_LS (fit);
  guint i;
  gdouble final_error = 0.0;
  gsl_matrix *J = fit_gsl_ls->ls->J;
  gsl_vector *f = fit_gsl_ls->ls->f;
  
  for (i = 0; i < J->size2; i++)
  {
    guint j;
    gdouble sum = 0.0, max_abs, min, max;
    gsl_vector_view J_col = gsl_matrix_column (J, i);
    gsl_vector_memcpy (ncm_vector_gsl (fit->f), f);
    gsl_vector_mul (ncm_vector_gsl (fit->f), &J_col.vector);
    for (j = 0; j < ncm_vector_len (fit->f); j++)
      sum += ncm_vector_get (fit->f, j);
    gsl_vector_minmax (ncm_vector_gsl (fit->f), &min, &max);
    max_abs = GSL_MAX (fabs(min), fabs(max));
    final_error += fabs (sum / max_abs);
  }

  return final_error;
}

/**
 * ncm_fit_gsl_ls_new:
 * @lh: FIXME
 * @mset: FIXME
 * @gtype: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmFit *
ncm_fit_gsl_ls_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  return g_object_new (NCM_TYPE_FIT_GSL_LS, 
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       NULL
                       );
}
