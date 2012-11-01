/***************************************************************************
 *            nc_nlopt.c
 *
 *  Sat Apr  3 16:07:02 2010
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
 * SECTION:nc_nlopt
 * @title: NLopt interface
 * @short_description: Interface for NLopt optmization library
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#ifdef NUMCOSMO_HAVE_NLOPT

#include "likelihood/nc_nlopt.h"

typedef gdouble (*_NcmFitNLOptOldFunc) (gint n, const gdouble *x, gdouble *grad, gpointer userdata);
static gdouble _ncm_fit_nlopt_func (guint n, const gdouble *x, gdouble *grad, gpointer userdata);

/***************************************************************************
 *
 *
 ****************************************************************************/

#undef _ALGO_NAME
#define _ALGO_NAME "NLopt"

/**
 * ncm_fit_run_nlopt: (skip)
 * @fit: a #NcmFit
 * @algo: FIXME
 * @niters: FIXME
 * @mtype: a #NcmFitRunMsgs
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_run_nlopt (NcmFit *fit, nlopt_algorithm algo, gint niters, NcmFitRunMsgs mtype)
{
  gint i;
  gdouble minf;

  switch (algo)
  {
	case NLOPT_LN_COBYLA:
	  fit->solver_name = _ALGO_NAME":COBYLA";
	  break;
	case NLOPT_LN_BOBYQA:
	  fit->solver_name = _ALGO_NAME":BOBYQA";
	  break;
	case NLOPT_LN_NEWUOA:
	  fit->solver_name = _ALGO_NAME":NEWUOA";
	  break;
	case NLOPT_LN_NEWUOA_BOUND:
	  fit->solver_name = _ALGO_NAME":NEWUOA Bound";
	  break;
	case NLOPT_LN_PRAXIS:
	  fit->solver_name = _ALGO_NAME":PRAXIS";
	  break;
	case NLOPT_LN_NELDERMEAD:
	  fit->solver_name = _ALGO_NAME":NelderMead";
	  break;
	case NLOPT_LN_SBPLX:
	  fit->solver_name = _ALGO_NAME":Sbplx";
	  break;
	case NLOPT_LD_MMA:
	  fit->solver_name = _ALGO_NAME":Method of Moving Asymptotes";
	  break;
	case NLOPT_LD_LBFGS:
	  fit->solver_name = _ALGO_NAME":Low-storage BFGS";
	  break;
	case NLOPT_LD_TNEWTON_PRECOND_RESTART:
	  fit->solver_name = _ALGO_NAME":TNewtonPrecondRestart";
	  break;
	case NLOPT_LD_TNEWTON_PRECOND:
	  fit->solver_name = _ALGO_NAME":TNewtonPrecond";
	  break;
	case NLOPT_LD_TNEWTON_RESTART:
	  fit->solver_name = _ALGO_NAME":TNewtonRestart";
	  break;
	case NLOPT_LD_TNEWTON:
	  fit->solver_name = _ALGO_NAME":TNewton";
	  break;
	case NLOPT_LD_VAR1:
	  fit->solver_name = _ALGO_NAME":Shifted limited-memory variable-metric rank 1";
	  break;
	case NLOPT_LD_VAR2:
	  fit->solver_name = _ALGO_NAME":Shifted limited-memory variable-metric rank 2";
	  break;
	default:
	  fit->solver_name = _ALGO_NAME;
  }

  fit->mtype = mtype;

  if (fit->fparam_len == 0)
  {
	ncm_fit_m2lnL_val (fit, &fit->m2lnL);
	fit->sqrt_m2lnL = sqrt (fit->m2lnL);
	fit->m2lnL_dof = fit->m2lnL / fit->dof;
	return TRUE;
  }

  ncm_mset_fparams_get_vector (fit->mset, fit->x);
  fit->mtype = mtype;
  ncm_fit_log_start (fit);

  {
	guint free_params_len = ncm_mset_fparams_len (fit->mset);
	gdouble *lb = g_slice_alloc (sizeof(gdouble) * fit->fparam_len);
	gdouble *ub = g_slice_alloc (sizeof(gdouble) * fit->fparam_len);
	gdouble *pabs = g_slice_alloc (sizeof(gdouble) * fit->fparam_len);
	gdouble *pscale = g_slice_alloc (sizeof(gdouble) * fit->fparam_len);

	for (i = 0; i < free_params_len; i++)
	{
	  lb[i] = ncm_mset_fparam_get_lower_bound (fit->mset, i);
	  ub[i] = ncm_mset_fparam_get_upper_bound (fit->mset, i);
	  pabs[i] = ncm_mset_fparam_get_abstol (fit->mset, i);
	  pscale[i] = ncm_mset_fparam_get_scale (fit->mset, i);
	}

#ifdef HAVE_NLOPT_2_2
	{
	  gboolean nlopt_alloc = FALSE;
	  nlopt_result ret;

	  if (fit->nlopt == NULL)
		nlopt_alloc = TRUE;
	  else if (fit->nlopt_algo != algo)
	  {
		nlopt_destroy (fit->nlopt);
		nlopt_alloc = TRUE;
	  }
	  if (nlopt_alloc)
	  {
		fit->nlopt = nlopt_create(algo, fit->fparam_len);
		fit->nlopt_algo = algo;
		nlopt_set_min_objective (fit->nlopt, &_ncm_fit_nlopt_func, fit);
		nlopt_set_lower_bounds (fit->nlopt, lb);
		nlopt_set_upper_bounds (fit->nlopt, ub);
		nlopt_set_ftol_rel (fit->nlopt, fit->m2lnL_prec_target);
		if (NCM_FIT_DEFAULT_M2LNL_ABSTOL != 0)
		  nlopt_set_ftol_abs (fit->nlopt, NCM_FIT_DEFAULT_M2LNL_ABSTOL);
		nlopt_set_xtol_rel (fit->nlopt, fit->params_prec_target);
		nlopt_set_xtol_abs (fit->nlopt, pabs);
		nlopt_set_maxeval (fit->nlopt, niters);
		nlopt_set_initial_step (fit->nlopt, pscale);
	  }

	  ret = nlopt_optimize (fit->nlopt, ncm_vector_gsl (fit->x)->data, &minf);
	  fit->m2lnL_prec = nlopt_get_ftol_rel (fit->nlopt);
	  fit->params_prec = nlopt_get_xtol_rel (fit->nlopt);

	  if (ret < 0)
		ncm_fit_log_step_error (fit, "(%d)", ret);
	}
#else
	{
	  gint ret;
	  ret = nlopt_minimize (algo,
	                        fit->fparam_len,
	                        (_NcmFitNLOptOldFunc)&_ncm_fit_nlopt_func, fit,
	                        lb, ub,
	                        fit->x->data,
	                        &minf,
	                        -HUGE_VAL,
	                        fit->m2lnL_prec, NCM_FIT_DEFAULT_M2LNL_ABSTOL,
	                        fit->param_prec, pabs,
	                        niters, 0);
	  if (ret < 0)
		ncm_fit_log_step_error (fit, "(%d)", ret);
	}
#endif

	g_slice_free1 (sizeof(gdouble) * fit->fparam_len, lb);
	g_slice_free1 (sizeof(gdouble) * fit->fparam_len, ub);
	g_slice_free1 (sizeof(gdouble) * fit->fparam_len, pabs);
	g_slice_free1 (sizeof(gdouble) * fit->fparam_len, pscale);
  }

  fit->sqrt_m2lnL = sqrt(minf);
  fit->m2lnL = minf;
  fit->m2lnL_dof = minf / fit->dof;

  ncm_fit_log_end (fit);

  ncm_mset_fparams_set_vector (fit->mset, fit->x);

  return TRUE;
}

/***************************************************************************
 *
 *
 ****************************************************************************/

static gdouble
_ncm_fit_nlopt_func (guint n, const gdouble *x, gdouble *grad, gpointer userdata)
{
  NcmFit *fit = NCM_FIT (userdata);
  gdouble m2lnL;

  fit->niter++;
  ncm_mset_fparams_set_array (fit->mset, x);

  if (grad)
  {
	NcmVector *gradv = ncm_vector_new_data_static (grad, n, 1);
	fit->m2lnL_val_grad (fit, &m2lnL, gradv);
	ncm_vector_free (gradv);
  }
  else
	ncm_fit_m2lnL_val (fit, &m2lnL);

  ncm_fit_log_step (fit, m2lnL);
  if (!gsl_finite (m2lnL))
	return GSL_POSINF;
  return m2lnL;
}

#endif /* NUMCOSMO_HAVE_NLOPT */
