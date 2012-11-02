/***************************************************************************
 *            levmar.c
 *
 *  Wed Feb 24 21:20:09 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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
 * SECTION:levmar
 * @title: Least Squares Algorithims -- Levmar
 * @short_description: Interface for Levenberg-Marquardt nonlinear least squares algorithm library
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "likelihood/levmar.h"

#include <gsl/gsl_blas.h>
#ifdef NUMCOSMO_HAVE_LEVMAR
#ifdef NC_LEVMAR_NEED_PREFIX
#include <levmar/levmar.h>
#else
#include <levmar.h>
#endif /* NC_LEVMAR_NEED_PREFIX */

static void nc_residual_levmar_f (gdouble *p, gdouble *hx, gint m, gint n, gpointer adata);
static void nc_residual_levmar_J (gdouble *p, gdouble *j, gint m, gint n, gpointer adata);

/***************************************************************************
 *
 *
 ****************************************************************************/

#undef _ALGO_NAME
#define _ALGO_NAME "liblevmar:dlevmar_der"

/**
 * ncm_fit_run_levmar_der:
 * @fit: FIXME
 * @niters: FIXME
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_run_levmar_der (NcmFit *fit, gint niters, NcmFitRunMsgs mtype)
{
  gdouble info[LM_INFO_SZ];
  //  gdouble opts[5];
  gint ret;

  if (fit->fparam_len == 0)
  {
	ncm_fit_m2lnL_val (fit, &fit->m2lnL);
	fit->sqrt_m2lnL = sqrt (fit->m2lnL);
	fit->m2lnL_dof = fit->m2lnL / fit->dof;
	return TRUE;
  }

  if (fit->levmar == NULL)
	fit->levmar = g_slice_alloc (LM_DIF_WORKSZ (fit->fparam_len, fit->data_len) * sizeof(gdouble));

  ncm_mset_fparams_get_vector (fit->mset, fit->x);

  fit->solver_name = _ALGO_NAME;
  fit->mtype = mtype;
  ncm_fit_log_start (fit);

  g_assert (ncm_vector_gsl (fit->f)->stride == 1 &&
            ncm_vector_gsl (fit->x)->stride == 1 &&
            NCM_MATRIX_GSL (fit->covar)->tda == NCM_MATRIX_NCOLS (fit->covar));

  ret = dlevmar_der (
                     &nc_residual_levmar_f, &nc_residual_levmar_J,
                     ncm_vector_gsl (fit->x)->data, ncm_vector_gsl (fit->f)->data, fit->fparam_len, fit->data_len,
                     niters, NULL, info, fit->levmar, NCM_MATRIX_GSL (fit->covar)->data, fit
                     );
  if (ret < 0)
	ncm_fit_log_step_error (fit, "(%d)", ret);


  fit->sqrt_m2lnL = sqrt(info[1]);
  fit->m2lnL = info[1];
  fit->m2lnL_dof = info[1] / fit->dof;
  fit->m2lnL_prec = info[2] / info[1];
  fit->niter = info[5];

  ncm_fit_log_end (fit);

  ncm_mset_fparams_set_vector (fit->mset, fit->x);

  return TRUE;
}

#undef _ALGO_NAME
#define _ALGO_NAME "liblevmar:dlevmar_dif"

/**
 * ncm_fit_run_levmar_dif:
 * @fit: FIXME
 * @niters: FIXME
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_run_levmar_dif (NcmFit *fit, gint niters, NcmFitRunMsgs mtype)
{
  gdouble info[LM_INFO_SZ];
  //  gdouble opts[5];
  gint ret;

  if (fit->fparam_len == 0)
  {
	ncm_fit_m2lnL_val (fit, &fit->m2lnL);
	fit->sqrt_m2lnL = sqrt (fit->m2lnL);
	fit->m2lnL_dof = fit->m2lnL / fit->dof;
	return TRUE;
  }

  if (fit->levmar == NULL)
	fit->levmar = g_slice_alloc (LM_DIF_WORKSZ (fit->fparam_len, fit->data_len) * sizeof(gdouble));

  ncm_mset_fparams_get_vector (fit->mset, fit->x);

  fit->solver_name = _ALGO_NAME;
  fit->mtype = mtype;
  ncm_fit_log_start (fit);

  g_assert (ncm_vector_gsl (fit->f)->stride == 1 &&
            ncm_vector_gsl (fit->x)->stride == 1 &&
            NCM_MATRIX_GSL (fit->covar)->tda == NCM_MATRIX_NCOLS (fit->covar));

  ret = dlevmar_dif (
                     &nc_residual_levmar_f,
                     ncm_vector_gsl (fit->x)->data, ncm_vector_gsl (fit->f)->data, fit->fparam_len, fit->data_len,
                     niters, NULL, info, fit->levmar, NCM_MATRIX_GSL (fit->covar)->data, fit
                     );

  if (ret < 0)
	ncm_fit_log_step_error (fit, "(%d)", ret);

  fit->sqrt_m2lnL = sqrt(info[1]);
  fit->m2lnL = info[1];
  fit->m2lnL_dof = info[1] / fit->dof;
  fit->m2lnL_prec = info[2] / info[1];
  fit->niter = info[5];

  ncm_fit_log_end (fit);

  ncm_mset_fparams_set_vector (fit->mset, fit->x);

  return TRUE;
}

/***************************************************************************
 *
 *
 ****************************************************************************/

static void
nc_residual_levmar_f (gdouble *p, gdouble *hx, gint m, gint n, gpointer adata)
{
  NcmFit *fit = NCM_FIT (adata);
  NcmVector *f = ncm_vector_new_data_static (hx, n, 1);

  fit->niter++;
  ncm_mset_fparams_set_array (fit->mset, p);
  ncm_fit_leastsquares_f (fit, f);

  ncm_fit_log_step (fit, gsl_pow_2(gsl_blas_dnrm2 (ncm_vector_gsl (f))));
  ncm_vector_free (f);
}

static void
nc_residual_levmar_J (gdouble *p, gdouble *j, gint m, gint n, gpointer adata)
{
  NcmFit *fit = NCM_FIT (adata);
  NcmMatrix *J = ncm_matrix_new_data_static (j, n, m);

  ncm_mset_fparams_set_array (fit->mset, p);

  fit->ls_J (fit, J);
  ncm_matrix_free (J);
}

#endif /* NUMCOSMO_HAVE_LEVMAR */
