/***************************************************************************
 *            multimin.c
 *
 *  Mon Jun 11 12:08:20 2007
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
 * SECTION:multimin
 * @title: Non-linear Minimization -- GSL
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "likelihood/multimin.h"

#include <gsl/gsl_blas.h>

/***************************************************************************
 * Calculate the best fit
 *
 ****************************************************************************/

static gdouble nc_residual_multimin_f (const gsl_vector *x, gpointer p);
static void nc_residual_multimin_df (const gsl_vector *x, gpointer p, gsl_vector *df);
static void nc_residual_multimin_fdf (const gsl_vector *x, gpointer p, gdouble *f, gsl_vector *df);

/**
 * ncm_fit_run_mm:
 * @fit: FIXME
 * @niters: FIXME
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_run_mm (NcmFit *fit, gint niters, NcmFitRunMsgs mtype)
{
  const gsl_multimin_fdfminimizer_type *T[] = {
	gsl_multimin_fdfminimizer_vector_bfgs2,
	gsl_multimin_fdfminimizer_conjugate_fr,
	gsl_multimin_fdfminimizer_conjugate_pr
  };
  gsl_multimin_fdfminimizer *s;
  gsl_multimin_function_fdf f;
  gdouble err_a = 1.0e300;
  gdouble err_b = 1.0e-1;
  gint status;
  gdouble prec = NCM_FIT_DEFAULT_M2LNL_RELTOL;
  gint i;

  if (fit->fparam_len == 0)
  {
	ncm_fit_m2lnL_val (fit, &fit->m2lnL);
	fit->sqrt_m2lnL = sqrt (fit->m2lnL);
	fit->m2lnL_dof = fit->m2lnL / fit->dof;
	return TRUE;
  }

  ncm_mset_fparams_get_vector (fit->mset, fit->x);

  for (i = 0; i < fit->fparam_len; i++)
  {
	gdouble pscale = ncm_mset_fparam_get_scale (fit->mset, i);
	err_a = GSL_MIN (err_a, pscale);
  }

  f.f   = &nc_residual_multimin_f;
  f.df  = &nc_residual_multimin_df;
  f.fdf = &nc_residual_multimin_fdf;
  f.n   = fit->fparam_len;
  f.params = fit;

  s = gsl_multimin_fdfminimizer_alloc (T[0], fit->fparam_len);
  gsl_multimin_fdfminimizer_set (s, &f, ncm_vector_gsl (fit->x), err_a, err_b);

  fit->solver_name = gsl_multimin_fdfminimizer_name (s);
  fit->mtype = mtype;
  ncm_fit_log_start (fit);

  do
  {
	gdouble pscale;
	fit->niter++;
	status = gsl_multimin_fdfminimizer_iterate (s);
	pscale = prec * (s->f != 0.0 ? s->f : 1.0);

	if (fit->niter == 1 && !gsl_finite(s->f))
	{
	  ncm_mset_fparams_set_vector (fit->mset, fit->x);
	  return FALSE;
	}

	if (status == GSL_ENOPROG)
	{
	  if (mtype > NCM_FIT_RUN_MSGS_NONE)
		ncm_fit_log_step_error (fit, gsl_strerror (status));
	  status = GSL_SUCCESS;
	}
	else
	  status = gsl_multimin_test_gradient (s->gradient, pscale);
	ncm_fit_log_step (fit, s->f);
  }
  while ( (status == GSL_CONTINUE) && (fit->niter < niters) );

  fit->m2lnL = s->f;
  fit->sqrt_m2lnL = sqrt(s->f);
  fit->m2lnL_dof = fit->m2lnL / fit->dof;
  fit->m2lnL_prec = gsl_blas_dnrm2 (s->gradient) / s->f;

  ncm_fit_log_end (fit);

  ncm_mset_fparams_set_gsl_vector (fit->mset, s->x);
  gsl_multimin_fdfminimizer_free (s);

  return TRUE;
}


/***************************************************************************
 *
 *
 ****************************************************************************/

static gdouble
nc_residual_multimin_f (const gsl_vector *x, gpointer p)
{
  NcmFit *fit = NCM_FIT (p);
  gdouble result;
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  ncm_fit_m2lnL_val (fit, &result);
  return result;
}

/***************************************************************************
 *
 *
 ****************************************************************************/

static void
nc_residual_multimin_df (const gsl_vector *x, gpointer p, gsl_vector *df)
{
  NcmFit *fit = NCM_FIT (p);
  NcmVector *dfv = ncm_vector_new_gsl (df);
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  fit->m2lnL_grad (fit, dfv);
  ncm_vector_free (dfv);
  return;
}

/***************************************************************************
 *
 *
 ****************************************************************************/

static void
nc_residual_multimin_fdf (const gsl_vector *x, gpointer p, gdouble *f, gsl_vector *df)
{
  NcmFit *fit = NCM_FIT (p);
  NcmVector *dfv = ncm_vector_new_gsl (df);
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  fit->m2lnL_val_grad (fit, f, dfv);
  ncm_vector_free (dfv);
  return;
}
