/***************************************************************************
 *            multimin_simplex.c
 *
 *  Tue Jun 19 10:55:06 2007
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
 * SECTION:multimin_simplex
 * @title: Non-linear Simplex Minimization -- GSL
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <string.h>
#include <glib.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>


/***************************************************************************
 * Calculate the best fit
 *
 ****************************************************************************/

static gdouble nc_residual_multimin_f (const gsl_vector *x, gpointer p);

/**
 * FIXME
 */
gboolean
ncm_fit_run_sp (NcmFit *fit, gint niters, NcmFitRunMsgs mtype)
{
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s;
  gsl_multimin_function f;
  gsl_vector *ss;
  gint status;
  gdouble prec = NCM_FIT_DEFAULT_M2LNL_RELTOL;
  gdouble last_size = 1e300;
  gulong still_count = 0;

  if (fit->fparam_len == 0)
  {
	ncm_fit_m2lnL_val (fit, &fit->m2lnL);
	fit->sqrt_m2lnL = sqrt (fit->m2lnL);
	fit->m2lnL_dof = fit->m2lnL / fit->dof;
	return TRUE;
  }

  f.f      = &nc_residual_multimin_f;
  f.n      = fit->fparam_len;
  f.params = fit;

  ss = gsl_vector_alloc (fit->fparam_len);

  ncm_mset_fparams_get_vector (fit->mset, fit->x);

  {
	gint i;
	for (i = 0; i < ncm_mset_fparams_len (fit->mset); i++)
	{
	  gdouble pscale = ncm_mset_fparam_get_scale (fit->mset, i);
	  gsl_vector_set (ss, i, pscale * 1e-3);
	}
  }

  s = gsl_multimin_fminimizer_alloc (T, fit->fparam_len);
  gsl_multimin_fminimizer_set (s, &f, ncm_vector_gsl (fit->x), ss);

  fit->solver_name = gsl_multimin_fminimizer_name (s);
  fit->mtype = mtype;
  ncm_fit_log_start (fit);
  do
  {
	gdouble size;
	fit->niter++;
	status = gsl_multimin_fminimizer_iterate (s);

	if (fit->niter == 1 && !gsl_finite(s->fval))
	{
	  ncm_mset_fparams_set_vector (fit->mset, fit->x);
	  return FALSE;
	}

	if (status)
	{
	  if (mtype > NCM_FIT_RUN_MSGS_NONE)
		ncm_fit_log_step_error (fit, gsl_strerror (status));
	}

	size = gsl_multimin_fminimizer_size (s);

	if (size == last_size)
	{
	  if (++still_count == 3)
	  {
		ncm_fit_log_step_error (fit, "size do not improve [prec: %8.5e size: %8.5e]", fit->solver_name, prec, size);
		status = GSL_SUCCESS;
	  }
	  else
		status = GSL_CONTINUE;
	}
	else
	{
	  still_count = 0;
	  last_size = size;
	  status = gsl_multimin_test_size (size, prec);
	  if (status == GSL_ETOL)
		status = GSL_CONTINUE;
	}
	ncm_fit_log_step (fit, s->fval);
  }
  while ( (status == GSL_CONTINUE) && (fit->niter < niters) );

  ncm_mset_fparams_set_gsl_vector (fit->mset, s->x);

  fit->m2lnL = s->fval;
  fit->sqrt_m2lnL = sqrt(s->fval);
  fit->m2lnL_dof = fit->m2lnL / fit->dof;
  fit->m2lnL_prec = last_size;

  ncm_fit_log_end (fit);

  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (ss);

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
