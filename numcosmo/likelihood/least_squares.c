/***************************************************************************
 *            least_squares.c
 *
 *  Mon Jun 11 12:04:33 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:least_squares
 * @title: Least Squares Algorithims -- GSL
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "likelihood/least_squares.h"

#include <gsl/gsl_blas.h>

typedef struct _NcHICosmoQPieceWise
{
  GList *models;
  GList *params;
  gint npieces;
  gdouble z_f;
  gdouble piece;
} NcHICosmoQPieceWise;

static int nc_residual_multifit_nlin_f (const gsl_vector *x, gpointer p, gsl_vector *f);
static int nc_residual_multifit_nlin_df (const gsl_vector *x, gpointer p, gsl_matrix *J);
static int nc_residual_multifit_nlin_fdf (const gsl_vector *x, gpointer p, gsl_vector *f, gsl_matrix *J);
static gdouble nc_least_squares_test_grad (NcmFit *fit);

/***************************************************************************
 *
 *
 ****************************************************************************/

#define _NC_MIN_PREC_RETRY (1e-3)

/**
 * ncm_fit_run_ls:
 * @fit: FIXME
 * @niters: FIXME
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_run_ls (NcmFit *fit, gint niters, NcmFitRunMsgs mtype)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_function_fdf f;
  gint status;
  gdouble prec = 1e-16;
  guint retries = 50;
  gboolean noretry = TRUE;

  if (fit->fparam_len == 0)
  {
    ncm_fit_m2lnL_val (fit, &fit->m2lnL);
    fit->sqrt_m2lnL = sqrt (fit->m2lnL);
    fit->m2lnL_dof = fit->m2lnL / fit->dof;
    return TRUE;
  }

  ncm_mset_fparams_get_vector (fit->mset, fit->x);

  f.f      = &nc_residual_multifit_nlin_f;
  f.df     = &nc_residual_multifit_nlin_df;
  f.fdf    = &nc_residual_multifit_nlin_fdf;
  f.p      = fit->fparam_len;
  f.n      = fit->data_len;
  f.params = fit;

  if (fit->least_squares == NULL)
  {
    T = gsl_multifit_fdfsolver_lmsder;
//    T = gsl_multifit_fdfsolver_lmder;
    fit->least_squares = gsl_multifit_fdfsolver_alloc (T, f.n, f.p);

  }
  fit->solver_name = gsl_multifit_fdfsolver_name (fit->least_squares);
  gsl_multifit_fdfsolver_set (fit->least_squares, &f, ncm_vector_gsl (fit->x));
  fit->mtype = mtype;
  ncm_fit_log_start (fit);

  do
  {
    fit->niter++;
    status = gsl_multifit_fdfsolver_iterate (fit->least_squares);

    if (fit->niter == 1 && !gsl_finite(gsl_blas_dnrm2(fit->least_squares->f)))
    {
      ncm_mset_fparams_set_vector (fit->mset, fit->x);
      //nc_params_print_all (fit->cp, stdout);
      //g_debug ("Init params infinity");
      return FALSE;
    }
    //gsl_multifit_gradient (fit->least_squares->J, fit->least_squares->f, fit->grad);
    prec = nc_least_squares_test_grad (fit);

    if (status)
    {
      if (mtype > NCM_FIT_RUN_MSGS_NONE)
        ncm_fit_log_step_error (fit, gsl_strerror (status));
      if (status == GSL_CONTINUE)
        break;
      if(!noretry && (fabs(gsl_blas_dasum(ncm_vector_gsl(fit->df))) < _NC_MIN_PREC_RETRY))
        noretry = TRUE;
/*
      if (noretry || !(retries--))
        prec = gsl_blas_dasum(fit->grad) * 1.001;
      else
      {
        gsl_vector_memcpy (, fit->least_squares->x);
        gsl_multifit_fdfsolver_set (fit->least_squares, &f, );
      }
*/
    }

    //status = gsl_multifit_test_gradient (fit->grad, prec);
    status = gsl_multifit_test_delta (fit->least_squares->dx, fit->least_squares->x, 1e-13, 1e-13);

    if (status == GSL_ETOL && FALSE)
    {
      if(!noretry && (fabs(gsl_blas_dasum(ncm_vector_gsl(fit->df))) < _NC_MIN_PREC_RETRY))
        noretry = TRUE;
      if (noretry || !(retries--))
        prec = gsl_blas_dasum(ncm_vector_gsl(fit->df)) * 1.001;
      else
      {
        gsl_vector_memcpy (ncm_vector_gsl (fit->x), fit->least_squares->x);
        gsl_multifit_fdfsolver_set (fit->least_squares, &f, ncm_vector_gsl(fit->x));
      }
      status = GSL_CONTINUE;
    }

    {
      gdouble sqrt_m2lnL = gsl_blas_dnrm2 (fit->least_squares->f);
      ncm_fit_log_step (fit, sqrt_m2lnL * sqrt_m2lnL);
    }
  }
  while (status == GSL_CONTINUE && fit->niter < niters);

  fit->sqrt_m2lnL = gsl_blas_dnrm2 (fit->least_squares->f);
  fit->m2lnL = fit->sqrt_m2lnL * fit->sqrt_m2lnL;
  fit->m2lnL_dof = fit->m2lnL / fit->dof;
  fit->m2lnL_prec = prec;

  ncm_fit_log_end (fit);

  ncm_fit_save (fit,
                ncm_vector_new_gsl (fit->least_squares->x),
                ncm_vector_new_gsl (fit->least_squares->f),
                ncm_matrix_new_gsl (fit->least_squares->J)
                );

  return TRUE;
}

/***************************************************************************
 *
 *
 ****************************************************************************/

static int
nc_residual_multifit_nlin_f (const gsl_vector *x, gpointer p, gsl_vector *f)
{
  NcmFit *fit = NCM_FIT (p);
  NcmVector *fv = ncm_vector_new_gsl (f);
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  ncm_fit_leastsquares_f (fit, fv);
  ncm_vector_free (fv);
  return GSL_SUCCESS;
}

/***************************************************************************
 *
 *
 ****************************************************************************/

static int
nc_residual_multifit_nlin_df (const gsl_vector *x, gpointer p, gsl_matrix *J)
{
  NcmFit *fit = NCM_FIT (p);
  NcmMatrix *Jm = ncm_matrix_new_gsl (J);
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  fit->ls_J (fit, Jm);
  ncm_matrix_free (Jm);
  return GSL_SUCCESS;
}

/***************************************************************************
 *
 *
 ****************************************************************************/

static int
nc_residual_multifit_nlin_fdf (const gsl_vector *x, gpointer p, gsl_vector *f, gsl_matrix *J)
{
  NcmFit *fit = NCM_FIT (p);
  NcmVector *fv = ncm_vector_new_gsl (f);
  NcmMatrix *Jm = ncm_matrix_new_gsl (J);
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  fit->ls_f_J (fit, fv, Jm);
  ncm_vector_free (fv);
  ncm_matrix_free (Jm);
  return GSL_SUCCESS;
}

/***************************************************************************
 *
 *
 ****************************************************************************/

static gdouble
nc_least_squares_test_grad (NcmFit *fit)
{
  guint i;
  gdouble final_error = 0.0;
  gsl_matrix *J = fit->least_squares->J;
  gsl_vector *f = fit->least_squares->f;
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
    final_error += fabs(sum / max_abs);
  }
  return final_error;
}
