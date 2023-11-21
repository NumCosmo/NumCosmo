/***************************************************************************
 *            ncm_fit_gsl_ls.c
 *
 *  Mon Jun 11 12:04:33 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * @title: NcmFitGSLLS
 * @short_description: Best-fit finder -- GSL least squares algorithms.
 *
 * This object implements a best-fit finder using the GSL least squares
 * algorithms. It is a subclass of #NcmFit.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_gsl_ls.h"
#include "math/ncm_fit_state.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmFitGSLLS
{
  /*< private >*/
  NcmFit parent_instance;
  gsl_multifit_fdfsolver *ls;
  gsl_multifit_function_fdf f;
  const gsl_multifit_fdfsolver_type *T;
};


G_DEFINE_TYPE (NcmFitGSLLS, ncm_fit_gsl_ls, NCM_TYPE_FIT)

static void
ncm_fit_gsl_ls_init (NcmFitGSLLS *fit_gsl_ls)
{
  fit_gsl_ls->ls = NULL;
  fit_gsl_ls->T  = gsl_multifit_fdfsolver_lmsder;
}

static gint _ncm_fit_gsl_ls_f (const gsl_vector *x, gpointer p, gsl_vector *f);
static gint _ncm_fit_gsl_ls_df (const gsl_vector *x, gpointer p, gsl_matrix *J);
static gint _ncm_fit_gsl_ls_fdf (const gsl_vector *x, gpointer p, gsl_vector *f, gsl_matrix *J);

static void
_ncm_fit_gsl_ls_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_gsl_ls_parent_class)->constructed (object);
  {
    NcmFitGSLLS *fit_gsl_ls = NCM_FIT_GSL_LS (object);
    NcmFit *fit             = NCM_FIT (fit_gsl_ls);
    NcmFitState *fstate     = ncm_fit_peek_state (fit);

    if (ncm_fit_state_get_fparam_len (fstate) > 0)
    {
      fit_gsl_ls->f.f      = &_ncm_fit_gsl_ls_f;
      fit_gsl_ls->f.df     = &_ncm_fit_gsl_ls_df;
      fit_gsl_ls->f.fdf    = &_ncm_fit_gsl_ls_fdf;
      fit_gsl_ls->f.p      = ncm_fit_state_get_fparam_len (fstate);
      fit_gsl_ls->f.n      = ncm_fit_state_get_data_len (fstate);
      fit_gsl_ls->f.params = fit;

      fit_gsl_ls->ls = gsl_multifit_fdfsolver_alloc (fit_gsl_ls->T,
                                                     fit_gsl_ls->f.n,
                                                     fit_gsl_ls->f.p);
    }
  }
}

static void
ncm_fit_gsl_ls_finalize (GObject *object)
{
  NcmFitGSLLS *fit_gsl_ls = NCM_FIT_GSL_LS (object);

  g_clear_pointer (&fit_gsl_ls->ls, gsl_multifit_fdfsolver_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_gsl_ls_parent_class)->finalize (object);
}

static NcmFit *_ncm_fit_gsl_ls_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
static void _ncm_fit_gsl_ls_reset (NcmFit *fit);
static gboolean _ncm_fit_gsl_ls_run (NcmFit *fit, NcmFitRunMsgs mtype);
static const gchar *_ncm_fit_gsl_ls_get_desc (NcmFit *fit);

static void
ncm_fit_gsl_ls_class_init (NcmFitGSLLSClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmFitClass *fit_class     = NCM_FIT_CLASS (klass);

  object_class->constructed = &_ncm_fit_gsl_ls_constructed;
  object_class->finalize    = &ncm_fit_gsl_ls_finalize;

  fit_class->copy_new = &_ncm_fit_gsl_ls_copy_new;
  fit_class->reset    = &_ncm_fit_gsl_ls_reset;
  fit_class->run      = &_ncm_fit_gsl_ls_run;
  fit_class->get_desc = &_ncm_fit_gsl_ls_get_desc;

  fit_class->is_least_squares = TRUE;
}

static NcmFit *
_ncm_fit_gsl_ls_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  NCM_UNUSED (fit);

  return ncm_fit_gsl_ls_new (lh, mset, gtype);
}

static void
_ncm_fit_gsl_ls_reset (NcmFit *fit)
{
  /* Chain up : start */
  NCM_FIT_CLASS (ncm_fit_gsl_ls_parent_class)->reset (fit);
  {
    NcmFitGSLLS *fit_gsl_ls = NCM_FIT_GSL_LS (fit);
    NcmFitState *fstate     = ncm_fit_peek_state (fit);

    if ((fit_gsl_ls->f.p != ncm_fit_state_get_fparam_len (fstate)) || (fit_gsl_ls->f.n != ncm_fit_state_get_data_len (fstate)))
    {
      g_clear_pointer (&fit_gsl_ls->ls, gsl_multifit_fdfsolver_free);

      if (ncm_fit_state_get_fparam_len (fstate) > 0)
      {
        fit_gsl_ls->f.p = ncm_fit_state_get_fparam_len (fstate);
        fit_gsl_ls->f.n = ncm_fit_state_get_data_len (fstate);
        fit_gsl_ls->ls  = gsl_multifit_fdfsolver_alloc (fit_gsl_ls->T,
                                                        fit_gsl_ls->f.n,
                                                        fit_gsl_ls->f.p);
      }
    }
  }
}

#define _NCM_FIT_GSL_LS_MIN_PREC_RETRY (1e-3)

gboolean
_ncm_fit_gsl_ls_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitGSLLS *fit_gsl_ls = NCM_FIT_GSL_LS (fit);
  NcmFitState *fstate     = ncm_fit_peek_state (fit);
  NcmMSet *mset           = ncm_fit_peek_mset (fit);
  gint status, info = 0;

  if (ncm_fit_equality_constraints_len (fit) || ncm_fit_inequality_constraints_len (fit))
    g_error ("_ncm_fit_gsl_ls_run: GSL algorithms do not support constraints.");

  g_assert (ncm_fit_state_get_fparam_len (fstate) != 0);

  ncm_mset_fparams_get_vector (mset, ncm_fit_state_peek_fparams (fstate));
  gsl_multifit_fdfsolver_set (fit_gsl_ls->ls, &fit_gsl_ls->f, ncm_vector_gsl (ncm_fit_state_peek_fparams (fstate)));

  status = gsl_multifit_fdfsolver_driver (fit_gsl_ls->ls,
                                          ncm_fit_get_maxiter (fit),
                                          ncm_fit_get_params_reltol (fit),
                                          ncm_fit_get_m2lnL_reltol (fit),
                                          ncm_fit_get_m2lnL_reltol (fit),
                                          &info
                                         );

  {
    NcmVector *_x = ncm_vector_new_gsl_static (fit_gsl_ls->ls->x);
    NcmVector *_f = ncm_vector_new_gsl_static (fit_gsl_ls->ls->f);
    NcmMatrix *_J = ncm_matrix_new (fit_gsl_ls->f.n, fit_gsl_ls->f.p);

    gsl_multifit_fdfsolver_jac (fit_gsl_ls->ls, ncm_matrix_gsl (_J));

    ncm_fit_params_set_vector (fit, _x);
    ncm_fit_state_set_params_prec (fstate, ncm_fit_get_params_reltol (fit));
    ncm_fit_state_set_ls (fstate, _f, _J);
    ncm_fit_state_set_niter (fstate, gsl_multifit_fdfsolver_niter (fit_gsl_ls->ls));

    ncm_vector_free (_x);
    ncm_vector_free (_f);
    ncm_matrix_free (_J);
  }

  if (status == GSL_SUCCESS)
    return TRUE;
  else
    return FALSE;
}

static gint
_ncm_fit_gsl_ls_f (const gsl_vector *x, gpointer p, gsl_vector *f)
{
  NcmFit *fit   = NCM_FIT (p);
  NcmMSet *mset = ncm_fit_peek_mset (fit);
  NcmVector *fv = ncm_vector_new_gsl_static (f);

  ncm_fit_params_set_gsl_vector (fit, x);

  if (!ncm_mset_params_valid (mset))
    return GSL_EDOM;

  ncm_fit_log_step (fit);
  ncm_fit_ls_f (fit, fv);

  ncm_vector_free (fv);

  return GSL_SUCCESS;
}

static gint
_ncm_fit_gsl_ls_df (const gsl_vector *x, gpointer p, gsl_matrix *J)
{
  NcmFit *fit   = NCM_FIT (p);
  NcmMSet *mset = ncm_fit_peek_mset (fit);
  NcmMatrix *Jm = ncm_matrix_new_gsl_static (J);

  ncm_fit_params_set_gsl_vector (fit, x);

  if (!ncm_mset_params_valid (mset))
    return GSL_EDOM;

  ncm_fit_log_step (fit);
  ncm_fit_ls_J (fit, Jm);

  ncm_matrix_free (Jm);

  return GSL_SUCCESS;
}

static gint
_ncm_fit_gsl_ls_fdf (const gsl_vector *x, gpointer p, gsl_vector *f, gsl_matrix *J)
{
  NcmFit *fit   = NCM_FIT (p);
  NcmMSet *mset = ncm_fit_peek_mset (fit);
  NcmVector *fv = ncm_vector_new_gsl_static (f);
  NcmMatrix *Jm = ncm_matrix_new_gsl_static (J);

  ncm_fit_params_set_gsl_vector (fit, x);

  if (!ncm_mset_params_valid (mset))
    return GSL_EDOM;

  ncm_fit_log_step (fit);
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

    desc = g_strdup_printf ("GSL Least Squares:%s",
                            fit_gsl_ls->ls != NULL ? gsl_multifit_fdfsolver_name (fit_gsl_ls->ls) : "not-set");
  }

  return desc;
}

/**
 * ncm_fit_gsl_ls_new:
 * @lh: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @gtype: a #NcmFitGradType
 *
 * Creates a new #NcmFitGSLLS object with the given likelihood, model set and
 * gradient type.
 *
 * Returns: (transfer full): a new #NcmFitGSLLS object.
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

