/***************************************************************************
 *            ncm_fit.c
 *
 *  Sat Aug 16 16:22:13 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:ncm_fit
 * @title: Cosmological Fit
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <string.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_deriv.h>

G_DEFINE_TYPE (NcmFit, ncm_fit, G_TYPE_OBJECT);

/**
 * ncm_fit_new:
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @type: a #NcmFitType.
 * @gtype: a #NcmFitGradType.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFit *
ncm_fit_new (NcLikelihood *lh, NcmMSet *mset, NcmFitType type, NcmFitGradType gtype)
{
  NcmFit *fit = g_object_new (NCM_TYPE_FIT, NULL);
  gint n = nc_dataset_get_n (lh->ds);

  g_assert (nc_dataset_all_init (lh->ds));

  fit->lh = nc_likelihood_copy (lh);
  fit->mset = ncm_mset_copy_ref (mset);
  fit->type = type;
  fit->gtype = gtype;

  ncm_mset_prepare_fparam_map (fit->mset);

  fit->data_len = n + nc_likelihood_priors_length (lh);
  fit->fparam_len = ncm_mset_fparam_len (fit->mset);
  fit->dof = fit->data_len - fit->fparam_len;

  g_assert (fit->data_len > 0);

  if (fit->fparam_len > 0)
  {
	fit->covar     = ncm_matrix_new (fit->fparam_len, fit->fparam_len);
	fit->hessian   = ncm_matrix_new (fit->fparam_len, fit->fparam_len);
	fit->x         = ncm_vector_new (fit->fparam_len);
	fit->df        = ncm_vector_new (fit->fparam_len);
	fit->J         = ncm_matrix_new (fit->data_len, fit->fparam_len);
  }

  fit->jackpv = ncm_vector_new (n);
  fit->f      = ncm_vector_new (fit->data_len);
  fit->bs     = gsl_vector_uint_alloc (fit->data_len);

  if ((gtype == NCM_FIT_GRAD_ANALYTICAL) && !nc_likelihood_has_m2lnL_grad (lh))
	g_error ("Likelihood do not support analytical gradient, try to use a numerical algorithm.");

  switch (gtype)
  {
	case NCM_FIT_GRAD_ANALYTICAL:
	  fit->ls_J = &ncm_fit_leastsquares_J;
	  fit->ls_f_J = &ncm_fit_leastsquares_f_J;
	  fit->m2lnL_grad = &ncm_fit_m2lnL_grad;
	  fit->m2lnL_val_grad = &ncm_fit_m2lnL_val_grad;
	  fit->diff_name = "Analytical gradient";
	  break;
	case NCM_FIT_GRAD_NUMDIFF_FORWARD:
	  fit->ls_J = &ncm_fit_numdiff_forward_leastsquares_J;
	  fit->ls_f_J = &ncm_fit_numdiff_forward_leastsquares_f_J;
	  fit->m2lnL_grad = &ncm_fit_numdiff_forward_m2lnL_grad;
	  fit->m2lnL_val_grad = &ncm_fit_numdiff_forward_m2lnL_val_grad;
	  fit->diff_name = "Numerical differentiantion (forward)";
	  break;
	case NCM_FIT_GRAD_NUMDIFF_CENTRAL:
	  fit->ls_J = &ncm_fit_numdiff_central_leastsquares_J;
	  fit->ls_f_J = &ncm_fit_numdiff_central_leastsquares_f_J;
	  fit->m2lnL_grad = &ncm_fit_numdiff_central_m2lnL_grad;
	  fit->m2lnL_val_grad = &ncm_fit_numdiff_central_m2lnL_val_grad;
	  fit->diff_name = "Numerical differentiantion (central)";
	  break;
	case NCM_FIT_GRAD_NUMDIFF_ACCURATE:
	  fit->ls_J = &ncm_fit_numdiff_central_leastsquares_J;
	  fit->ls_f_J = &ncm_fit_numdiff_central_leastsquares_f_J;
	  fit->m2lnL_grad = &ncm_fit_numdiff_accurate_m2lnL_grad;
	  fit->m2lnL_val_grad = &ncm_fit_numdiff_accurate_m2lnL_val_grad;
	  fit->diff_name = "Numerical differentiantion (Richardson extrapolation)";
	  break;
	default:
	  g_error ("Invalid gtype %d", gtype);
	  break;
  }

  return fit;
}

/**
 * ncm_fit_ref:
 * @fit: a #NcmFit.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmFit *
ncm_fit_ref (NcmFit *fit)
{
  return g_object_ref (fit);
}

/**
 * ncm_fit_copy:
 * @fit: a #NcmFit.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmFit *
ncm_fit_copy (NcmFit *fit)
{
  return ncm_fit_new (fit->lh, fit->mset, fit->type, fit->gtype);
}

/**
 * ncm_fit_free:
 * @fit: a #NcmFit.
 *
 * FIXME
 *
 */
void
ncm_fit_free (NcmFit *fit)
{
  g_object_unref (fit);
}

/**
 * ncm_fit_save:
 * @fit: a #NcmFit.
 * @x: a #NcmVector.
 * @f: a #NcmVector.
 * @J: a #NcmMatrix.
 *
 * FIXME
 */
void
ncm_fit_save (NcmFit *fit, NcmVector *x, NcmVector *f, NcmMatrix *J)
{
  ncm_mset_fparams_set_vector (fit->mset, x);
  fit->sqrt_m2lnL = gsl_blas_dnrm2(ncm_vector_gsl(f));
  fit->m2lnL = fit->sqrt_m2lnL*fit->sqrt_m2lnL;
  fit->m2lnL_dof = fit->m2lnL / fit->dof;
  ncm_vector_memcpy (fit->f, f);

  gsl_blas_dgemv (CblasTrans, 2.0, NCM_MATRIX_GSL(fit->J), ncm_vector_gsl(fit->f), 0.0, ncm_vector_gsl(fit->df));
  ncm_matrix_memcpy (fit->J, J);
}

/**
 * ncm_fit_covar_leastsquares_calc:
 * @fit: a #NcmFit.
 *
 * FIXME
 */
void
ncm_fit_covar_leastsquares_calc (NcmFit *fit)
{
  gsl_multifit_covar (NCM_MATRIX_GSL (fit->J), 0.0, NCM_MATRIX_GSL (fit->covar));
}

/**
 * ncm_fit_covar_fparam_var:
 * @fit: a #NcmFit.
 * @fpi: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_fparam_var (NcmFit *fit, guint fpi)
{
  return ncm_matrix_get (fit->covar, fpi, fpi);
}

/**
 * ncm_fit_covar_fparam_sd:
 * @fit: a #NcmFit.
 * @fpi: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_fparam_sd (NcmFit *fit, guint fpi)
{
  return sqrt (ncm_fit_covar_fparam_var (fit, fpi));
}

/**
 * ncm_fit_covar_fparam_cov:
 * @fit: a #NcmFit.
 * @fpi1: FIXME
 * @fpi2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_fparam_cov (NcmFit *fit, guint fpi1, guint fpi2)
{
  return ncm_matrix_get (fit->covar, fpi1, fpi2);
}

/**
 * ncm_fit_covar_fparam_cor:
 * @fit: a #NcmFit.
 * @fpi1: FIXME
 * @fpi2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_fparam_cor (NcmFit *fit, guint fpi1, guint fpi2)
{
  return ncm_fit_covar_fparam_cov (fit, fpi1, fpi2) / (ncm_fit_covar_fparam_sd (fit, fpi1) * ncm_fit_covar_fparam_sd (fit, fpi2));
}

/**
 * ncm_fit_covar_var:
 * @fit: a #NcmFit.
 * @gmid: a #NcmModelID.
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_var (NcmFit *fit, NcmModelID gmid, guint pid)
{
  gint fpi = ncm_mset_fparam_get_fpi (fit->mset, gmid, pid);
  if (fpi < 0)
	g_error ("Parameter (%d:%u) was not fit.", gmid, pid);
  return ncm_fit_covar_fparam_var (fit, fpi);
}

/**
 * ncm_fit_covar_sd:
 * @fit: a #NcmFit.
 * @gmid: a #NcmModelID.
 * @pid: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_sd (NcmFit *fit, NcmModelID gmid, guint pid)
{
  return sqrt (ncm_fit_covar_var (fit, gmid, pid));
}

/**
 * ncm_fit_covar_cov:
 * @fit: a #NcmFit
 * @gmid1: a #NcmModelID.
 * @pid1: FIXME
 * @gmid2: a #NcmModelID.
 * @pid2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_cov (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2)
{
  gint fpi1 = ncm_mset_fparam_get_fpi (fit->mset, gmid1, pid1);
  gint fpi2 = ncm_mset_fparam_get_fpi (fit->mset, gmid2, pid2);

  if (fpi1 < 0 || fpi1 < 0)
	g_error ("Parameters (%d:%u, %d:%u) were not fit.", gmid1, pid1, gmid2, pid2);

  return ncm_fit_covar_fparam_cov (fit, fpi1, fpi2);
}

/**
 * ncm_fit_covar_cor:
 * @fit: a #NcmFit
 * @gmid1: a #NcmModelID.
 * @pid1: FIXME
 * @gmid2: a #NcmModelID.
 * @pid2: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_covar_cor (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2)
{
  gint fpi1 = ncm_mset_fparam_get_fpi (fit->mset, gmid1, pid1);
  gint fpi2 = ncm_mset_fparam_get_fpi (fit->mset, gmid2, pid2);

  if (fpi1 < 0 || fpi1 < 0)
	g_error ("Parameters (%d:%u, %d:%u) were not fit.", gmid1, pid1, gmid2, pid2);

  return ncm_fit_covar_fparam_cov (fit, fpi1, fpi2) / (ncm_fit_covar_fparam_sd (fit, fpi1) * ncm_fit_covar_fparam_sd (fit, fpi2));
}

/**
 * ncm_fit_run:
 * @fit: a #NcmFit
 * @niters: FIXME
 * @mtype: a #NcmFitRunMsgs
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_fit_run (NcmFit *fit, gint niters, NcmFitRunMsgs mtype)
{
  gboolean ret;

  if (fit->dof < -1 && fit->type == NCM_FIT_TYPE_LEAST_SQUARES)
	g_error ("Cannot use algorithm NCM_FIT_TYPE_LEAST_SQUARES with dof = %d\n", fit->dof);

  fit->niter = 0;
  fit->n_func_eval = 0;
  fit->n_grad_eval = 0;

  switch (fit->type)
  {
	case NCM_FIT_TYPE_LEAST_SQUARES:
	  ret = ncm_fit_run_ls (fit, niters, mtype);
	  break;
	case NCM_FIT_TYPE_MULTIMIN:
	  ret = ncm_fit_run_mm (fit, niters, mtype);
	  break;
	case NCM_FIT_TYPE_SIMPLEX:
	  ret = ncm_fit_run_sp (fit, niters, mtype);
	  break;
#ifdef NUMCOSMO_HAVE_LEVMAR
	case NCM_FIT_TYPE_LEVMAR_DER:
	  ret = ncm_fit_run_levmar_der (fit, niters, mtype);
	  break;
	case NCM_FIT_TYPE_LEVMAR_DIF:
	  ret = ncm_fit_run_levmar_dif (fit, niters, mtype);
	  break;
#else
	case NCM_FIT_TYPE_LEVMAR_DER:
	case NCM_FIT_TYPE_LEVMAR_DIF:
	  g_error ("Levmar library not supported, recompile numcosmo with Levmar support.");
	  break;
#endif /* NUMCOSMO_HAVE_LEVMAR */
#ifdef NUMCOSMO_HAVE_NLOPT
	case NCM_FIT_TYPE_NLOPT_LN_COBYLA:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LN_COBYLA, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LN_BOBYQA:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LN_BOBYQA, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LN_NEWUOA:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LN_NEWUOA, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LN_NEWUOA_BOUND:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LN_NEWUOA_BOUND, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LN_PRAXIS:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LN_PRAXIS, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LN_NELDERMEAD:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LN_NELDERMEAD, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LN_SBPLX:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LN_SBPLX, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LD_MMA:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LD_MMA, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LD_LBFGS:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LD_LBFGS, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LD_TNEWTON_PRECOND_RESTART:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LD_TNEWTON_PRECOND_RESTART, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LD_TNEWTON_PRECOND:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LD_TNEWTON_PRECOND, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LD_TNEWTON_RESTART:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LD_TNEWTON_RESTART, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LD_TNEWTON:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LD_TNEWTON, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LD_VAR1:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LD_VAR1, niters, mtype);
	  break;
	case NCM_FIT_TYPE_NLOPT_LD_VAR2:
	  ret = ncm_fit_run_nlopt (fit, NLOPT_LD_VAR2, niters, mtype);
	  break;
#else
	case NCM_FIT_TYPE_NLOPT_LN_COBYLA:
	case NCM_FIT_TYPE_NLOPT_LN_BOBYQA:
	case NCM_FIT_TYPE_NLOPT_LN_NEWUOA:
	case NCM_FIT_TYPE_NLOPT_LN_NEWUOA_BOUND:
	case NCM_FIT_TYPE_NLOPT_LN_PRAXIS:
	case NCM_FIT_TYPE_NLOPT_LN_NELDERMEAD:
	case NCM_FIT_TYPE_NLOPT_LN_SBPLX:
	case NCM_FIT_TYPE_NLOPT_LD_MMA:
	case NCM_FIT_TYPE_NLOPT_LD_LBFGS:
	case NCM_FIT_TYPE_NLOPT_LD_TNEWTON_PRECOND_RESTART:
	case NCM_FIT_TYPE_NLOPT_LD_TNEWTON_PRECOND:
	case NCM_FIT_TYPE_NLOPT_LD_TNEWTON_RESTART:
	case NCM_FIT_TYPE_NLOPT_LD_TNEWTON:
	case NCM_FIT_TYPE_NLOPT_LD_VAR1:
	case NCM_FIT_TYPE_NLOPT_LD_VAR2:
	  g_error ("NLopt library not supported, recompile numcosmo with NLopt support.");
	  break;
#endif /* NUMCOSMO_HAVE_NLOPT */
	default:
	  g_assert_not_reached();
  }
  return ret;
}

/**
 * ncm_fit_gen_bootstrap:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_gen_bootstrap (NcmFit *fit)
{
  gsl_rng *rand = ncm_get_rng();
  guint i;
  fit->bootstrap = TRUE;
  for (i = 0; i < ncm_vector_len (fit->f); i++)
  {
	gint boot_el = gsl_rng_uniform_int (rand, ncm_vector_len (fit->f));
	gsl_vector_uint_set (fit->bs, i, boot_el);
  }
  return;
}

/**
 * ncm_fit_log_start:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_log_start (NcmFit *fit)
{
  g_timer_start (fit->timer);
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
	g_message ("# Interating using:\n");
	g_message ("#  - solver:            %s\n", fit->solver_name);
	g_message ("#  - differentiation:   %s\n", fit->diff_name);

	if (fit->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
	  g_message ("#");
  }

  fit->niter = 0;
  fit->n_func_eval = 0;
  fit->n_grad_eval = 0;
  return;
}

/**
 * ncm_fit_log_step_error:
 * @fit: a #NcmFit
 * @strerror: FIXME
 * @...: FIXME
 *
 * FIXME
 */
void
ncm_fit_log_step_error (NcmFit *fit, const gchar *strerror, ...)
{
  va_list ap;

  va_start (ap, strerror);
  if (fit->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
  {
	gchar *errmsg = g_strdup_vprintf (strerror, ap);
	g_message ("\n#  [%s] error = %s\n#", fit->solver_name, errmsg);
	g_free (errmsg);
  }
  else if (fit->mtype == NCM_FIT_RUN_MSGS_FULL)
  {
	gchar *errmsg = g_strdup_vprintf (strerror, ap);
	g_message ("#  [%s] error = %s\n", fit->solver_name, errmsg);
	g_free (errmsg);
  }
  va_end (ap);
  return;
}

/**
 * ncm_fit_log_end:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_log_end (NcmFit *fit)
{
  g_timer_stop (fit->timer);
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
	if (fit->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
	  g_message ("\n");
	g_message ("#  Minimum found with precision: % 8.5e (|df|/f or |dx|)\n", fit->m2lnL_prec);
  }
  ncm_fit_log_state (fit, fit->m2lnL);
  return;
}

/**
 * ncm_fit_log_state:
 * @fit: a #NcmFit
 * @m2lnL: FIXME
 *
 * FIXME
 */
void
ncm_fit_log_state (NcmFit *fit, gdouble m2lnL)
{
  gint i;
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
	gdouble elap_sec = g_timer_elapsed (fit->timer, NULL);
	gulong elap_min = elap_sec / 60;
	gulong elap_hour = elap_min / 60;
	gulong elap_day = elap_hour / 24;

	elap_sec = fmod (elap_sec, 60);
	elap_min = elap_min % 60;
	elap_hour =  elap_hour % 24;
	g_message ("#  Elapsed time: %02lu days, %02lu:%02lu:%010.7f\n", elap_day, elap_hour, elap_min, elap_sec);
	g_message ("#  iteration            [%06d]\n", fit->niter);
	g_message ("#  function evaluations [%06d]\n", fit->n_func_eval);
	g_message ("#  gradient evaluations [%06d]\n", fit->n_grad_eval);
	g_message ("#  degrees of freedom   [%06d]\n", fit->dof);
	g_message ("#  m2lnL     = %20.15g\n", m2lnL);
	g_message ("#  Fit parameters:\n#    ");
	for (i = 0; i < ncm_mset_fparam_len (fit->mset); i++)
	  g_message ("[% -20.15g] ", ncm_mset_fparam_get (fit->mset, i));
	g_message ("\n");
  }
  return;
}

/**
 * ncm_fit_log_step:
 * @fit: a #NcmFit
 * @m2lnL: FIXME
 *
 * FIXME
 */
void
ncm_fit_log_step (NcmFit *fit, gdouble m2lnL)
{
  if (fit->mtype == NCM_FIT_RUN_MSGS_FULL)
	ncm_fit_log_state (fit, m2lnL);
  else if (fit->mtype == NCM_FIT_RUN_MSGS_SIMPLE)
	g_message (".");
  return;
}

/**
 * ncm_fit_log_finish:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_log_finish (NcmFit *fit)
{
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
	g_message ("#  m2lnL/dof = %20.15g\n", fit->m2lnL_dof);
	g_message ("#  |m2lnL-dof|/sqrt(2*dof) = %20.15g,\n", fabs(fit->m2lnL-fit->dof)/sqrt(2.0 * fit->dof));
	g_message ("#  GoF_tt = %4.2f%% = (%4.2f + %4.2f)%%; GoF = %4.2f%%\n",
	           gsl_cdf_chisq_Q (fit->dof + fabs(fit->dof - fit->m2lnL), (gint)fit->dof)*100.0 +
	           gsl_cdf_chisq_P (fit->dof - fabs(fit->dof - fit->m2lnL), (gint)fit->dof)*100.0,
	           gsl_cdf_chisq_P (fit->dof - fabs(fit->dof - fit->m2lnL), (gint)fit->dof)*100.0,
	           gsl_cdf_chisq_Q (fit->dof + fabs(fit->dof - fit->m2lnL), (gint)fit->dof)*100.0,
	           gsl_cdf_chisq_Q (fit->m2lnL, (gint)fit->dof)*100.0);
  }
  return;
}

/**
 * ncm_fit_log_info:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_log_info (NcmFit *fit)
{
  nc_dataset_log_info (fit->lh->ds);
  ncm_mset_pretty_log (fit->mset);
  return;
}

/**
 * ncm_fit_log_covar:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_log_covar (NcmFit *fit)
{
  gint i, j;
  guint name_size = ncm_mset_max_fparam_name (fit->mset);
  guint free_params_len = ncm_mset_fparam_len (fit->mset);
  gchar *box = "---------------";

  g_message ("# Fitted parameters covariance matrix\n");
  g_message ("#                                        ");
  for (i = 0; i < name_size; i++) g_message (" ");

  for (i = 0; i < free_params_len; i++)
	i ? g_message ("%s",box) : g_message ("-%s",box);
  if (i)
	g_message ("\n");

  for (i = 0; i < free_params_len; i++)
  {
	NcmMSetPIndex *pi = ncm_mset_fparam_get_pi (fit->mset, i);
	const gchar *pname = ncm_mset_fparam_name (fit->mset, i);
	g_message ("# %*s[%02d%02d] = % -12.4g +/- % -12.4g |",
	           name_size, pname, pi->gmid, pi->pid,
	           ncm_mset_fparam_get (fit->mset, i),
	           ncm_fit_covar_fparam_sd (fit, i));
	for (j = 0; j < free_params_len; j++)
	{
	  g_message (" % -12.4g |", ncm_fit_covar_fparam_cor (fit, i, j));
	}
	g_message ("\n");
  }
  g_message ("#                                        ");
  for (i = 0; i < name_size; i++) g_message (" ");
  for (i = 0; i < free_params_len; i++)
	i ? g_message ("%s",box) : g_message ("-%s",box);
  if (i)
	g_message ("\n");

  return;
}

/**
 * ncm_fit_fisher_matrix_print:
 * @fit: a #NcmFit
 * @out: name of the file
 * @header: pointer to the command line
 *
 * This function print the command line (first line, commented), the cosmological
 * parameters' names which were fitted (second line, commented) and the Fisher Matrix.
 *
 */
void
ncm_fit_fishermatrix_print (NcmFit *fit, FILE *out, gchar *header)
{
  gint i, j;
  guint name_size = ncm_mset_max_param_name (fit->mset);
  guint free_params_len = ncm_mset_fparam_len (fit->mset);

  if (header != NULL)
	fprintf (out, "# %s\n# ", header);
  else
	fprintf (out, "# ");

  for (i = 0; i < free_params_len; i++)
  {
	const gchar *pname = ncm_mset_fparam_name (fit->mset, i);
	fprintf (out, "%*s[%02d] ", name_size, pname, i);
  }
  fprintf (out, "\n");

  for (i = 0; i < free_params_len; i++)
  {
	for (j = 0; j < free_params_len; j++)
	{
	  fprintf (out, " % -20.15g", ncm_fit_covar_fparam_cov (fit, i, j));
	}
	fprintf (out, "\n");
  }

  return;
}

/**
 * ncm_fit_m2lnL_val:
 * @fit: a #NcmFit
 * @m2lnL: (out): FIXME
   *
 * FIXME
 */
void
ncm_fit_m2lnL_val (NcmFit *fit, gdouble *m2lnL)
{
  nc_likelihood_m2lnL_val (fit->lh, fit->mset, m2lnL);
  fit->n_func_eval++;
  return;
}

/**
 * ncm_fit_data_m2lnL_val:
 * @fit: a #NcmFit
 * @data_m2lnL: (out): FIXME
   *
 * FIXME
 */
void
ncm_fit_data_m2lnL_val (NcmFit *fit, gdouble *data_m2lnL)
{
  nc_likelihood_data_m2lnL_val (fit->lh, fit->mset, data_m2lnL);
  return;
}

/**
 * ncm_fit_priors_m2lnL_val:
 * @fit: a #NcmFit
 * @priors_m2lnL: (out): FIXME
   *
 * FIXME
 */
void
ncm_fit_priors_m2lnL_val (NcmFit *fit, gdouble *priors_m2lnL)
{
  nc_likelihood_priors_m2lnL_val (fit->lh, fit->mset, priors_m2lnL);
  return;
}

/**
 * ncm_fit_m2lnL_grad:
 * @fit: a #NcmFit
 * @df: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_m2lnL_grad (NcmFit *fit, NcmVector *df)
{
  nc_likelihood_m2lnL_grad (fit->lh, fit->mset, df);
}

/**
 * ncm_fit_numdiff_forward_m2lnL_grad:
 * @fit: a #NcmFit
 * @grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_numdiff_forward_m2lnL_grad (NcmFit *fit, NcmVector *grad)
{
  gdouble m2lnL;
  ncm_fit_numdiff_forward_m2lnL_val_grad (fit, &m2lnL, grad);
}

/**
 * ncm_fit_numdiff_forward_m2lnL_val_grad:
 * @fit: a #NcmFit
 * @m2lnL: (out): FIXME
   * @grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_numdiff_forward_m2lnL_val_grad (NcmFit *fit, gdouble *m2lnL, NcmVector *grad)
{
  gint i;
  guint free_params_len = ncm_mset_max_fparam_name (fit->mset);

  ncm_fit_m2lnL_val (fit, m2lnL);

  for (i = 0; i < free_params_len; i++)
  {
	const gdouble p = ncm_mset_fparam_get (fit->mset, i);
	const gdouble p_scale = ncm_mset_fparam_get_scale (fit->mset, i);
	const gdouble h = (p_scale * NCM_FIT_NUMDIFF_SCALE);
	const gdouble pph = p + h;
	const gdouble one_h = 1.0 / h;
	gdouble m2lnL_pph;
	ncm_mset_fparam_set (fit->mset, i, pph);
	ncm_fit_m2lnL_val (fit, &m2lnL_pph);
	ncm_vector_set (grad, i, (m2lnL_pph - *m2lnL) * one_h);
	ncm_mset_fparam_set (fit->mset, i, p);
  }
}

/**
 * ncm_fit_numdiff_central_m2lnL_grad:
 * @fit: a #NcmFit
 * @grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_numdiff_central_m2lnL_grad (NcmFit *fit, NcmVector *grad)
{
  gint i;
  guint free_params_len = ncm_mset_max_fparam_name (fit->mset);

  for (i = 0; i < free_params_len; i++)
  {
	const gdouble p = ncm_mset_fparam_get (fit->mset, i);
	const gdouble p_scale = ncm_mset_fparam_get_scale (fit->mset, i);
	const gdouble h = (p_scale * NCM_FIT_NUMDIFF_SCALE);
	const gdouble pph = p + h;
	const gdouble pmh = p - h;
	const gdouble one_2h = 1.0 / (2.0 * h);
	gdouble m2lnL_pph, m2lnL_pmh;

	ncm_mset_fparam_set (fit->mset, i, pph);
	ncm_fit_m2lnL_val (fit, &m2lnL_pph);

	ncm_mset_fparam_set (fit->mset, i, pmh);
	ncm_fit_m2lnL_val (fit, &m2lnL_pmh);

	ncm_vector_set (grad, i, (m2lnL_pph - m2lnL_pmh) * one_2h);
	ncm_mset_fparam_set (fit->mset, i, p);
  }
}

/**
 * ncm_fit_numdiff_central_m2lnL_val_grad:
 * @fit: a #NcmFit
 * @m2lnL: FIXME
 * @grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_numdiff_central_m2lnL_val_grad (NcmFit *fit, gdouble *m2lnL, NcmVector *grad)
{
  ncm_fit_m2lnL_val (fit, m2lnL);
  ncm_fit_numdiff_central_m2lnL_grad (fit, grad);
}

typedef struct __ncm_fit_numdiff_1
{
  NcmFit *fit;
  guint n;
} __ncm_fit_numdiff_1;

static gdouble
_ncm_fit_numdiff_1_m2lnL (gdouble x, gpointer userdata)
{
  gdouble res;
  __ncm_fit_numdiff_1 *nd = (__ncm_fit_numdiff_1 *)userdata;
  ncm_mset_fparam_set (nd->fit->mset, nd->n, x);
  ncm_fit_m2lnL_val (nd->fit, &res);
  return res;
}

/**
 * ncm_fit_numdiff_accurate_m2lnL_grad:
 * @fit: a #NcmFit
 * @grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_numdiff_accurate_m2lnL_grad (NcmFit *fit, NcmVector *grad)
{
  gsl_function F;
  __ncm_fit_numdiff_1 nd;
  guint free_params_len = ncm_mset_max_fparam_name (fit->mset);
  gint i;

  nd.fit = fit;
  F.params = &nd;
  F.function = &_ncm_fit_numdiff_1_m2lnL;

  for (i = 0; i < free_params_len; i++)
  {
	const gdouble p = ncm_mset_fparam_get (fit->mset, i);
	const gdouble p_scale = ncm_mset_fparam_get_scale (fit->mset, i);
	gdouble err, diff;
	nd.n = i;
	diff = ncm_numdiff_1 (&F, p, p_scale, &err);
	ncm_vector_set (grad, i, diff);
	ncm_mset_fparam_set (fit->mset, i, p);
  }
}

/**
 * ncm_fit_numdiff_accurate_m2lnL_val_grad:
 * @fit: a #NcmFit
 * @m2lnL: (out): FIXME
   * @grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_numdiff_accurate_m2lnL_val_grad (NcmFit *fit, gdouble *m2lnL, NcmVector *grad)
{
  ncm_fit_m2lnL_val (fit, m2lnL);
  ncm_fit_numdiff_accurate_m2lnL_grad (fit, grad);
}

typedef struct __ncm_fit_numdiff_2
{
  NcmFit *fit;
  guint n1;
  guint n2;
  gdouble v;
  gdouble p1_scale;
  gdouble p2_scale;
} _ncm_fit_numdiff_2;

static gdouble
_ncm_fit_numdiff_2_m2lnL (gdouble u, gpointer userdata)
{
  gdouble res;
  _ncm_fit_numdiff_2 *nd = (_ncm_fit_numdiff_2 *)userdata;
  if (nd->n1 == nd->n2)
	ncm_mset_fparam_set (nd->fit->mset, nd->n1, u);
  else
  {
	ncm_mset_fparam_set (nd->fit->mset, nd->n1, (u + nd->v) / nd->p1_scale);
	ncm_mset_fparam_set (nd->fit->mset, nd->n2, (u - nd->v) / nd->p2_scale);
  }
  ncm_fit_m2lnL_val (nd->fit, &res);
  return res;
}

/**
 * ncm_fit_numdiff_m2lnL_hessian:
 * @fit: a #NcmFit
 * @H: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_numdiff_m2lnL_hessian (NcmFit *fit, NcmMatrix *H)
{
  gsl_function F;
  _ncm_fit_numdiff_2 nd;
  gint i, j;
  gdouble fx;
  const gdouble target_err = 1e-5;
  guint free_params_len = ncm_mset_fparams_len (fit->mset);

  ncm_fit_m2lnL_val (fit, &fx);

  nd.fit = fit;
  F.params = &nd;
  F.function = &_ncm_fit_numdiff_2_m2lnL;

  /* Diagonal */
  for (i = 0; i < free_params_len; i++)
  {
	const gdouble p = ncm_mset_fparam_get (fit->mset, i);
	const gdouble p_scale = ncm_mset_fparam_get_scale (fit->mset, i);
	gdouble err, diff;
	nd.n1 = i;
	nd.n2 = i;
	diff = ncm_numdiff_2_err (&F, &fx, p, p_scale, target_err, &err);
	if (fabs(err / diff) > target_err)
	  g_warning ("ncm_fit_numdiff_m2lnL_hessian: effective error (% 20.15e) bigger then (% 20.15e)", fabs(err / diff), target_err);
	ncm_matrix_set (H, i, i, diff);
	ncm_mset_fparam_set (fit->mset, i, p);
  }

  for (i = 0; i < free_params_len; i++)
  {
	for (j = i + 1; j < free_params_len; j++)
	{
	  const gdouble p1_scale = 1.0 / ncm_mset_fparam_get_scale (fit->mset, i);
	  const gdouble p2_scale = 1.0 / ncm_mset_fparam_get_scale (fit->mset, j);
	  const gdouble p1 = ncm_mset_fparam_get (fit->mset, i);
	  const gdouble p2 = ncm_mset_fparam_get (fit->mset, j);
	  const gdouble u = (p1_scale * p1 + p2_scale * p2) / 2.0;
	  const gdouble v = (p1_scale * p1 - p2_scale * p2) / 2.0;
	  const gdouble u_scale = 1.0;
	  gdouble err, diff;
	  nd.n1 = i;
	  nd.n2 = j;
	  nd.v = v;
	  nd.p1_scale = p1_scale;
	  nd.p2_scale = p2_scale;
	  diff = ncm_numdiff_2_err (&F, &fx, u, u_scale, target_err, &err);
	  if (fabs(err / diff) > target_err)
		g_warning ("ncm_fit_numdiff_m2lnL_hessian: effective error (% 20.15e) bigger then (% 20.15e)", fabs(err / diff), target_err);
	  ncm_matrix_set (H, i, j,
	                  0.5 * ( p1_scale * p2_scale * diff -
	                         (p2_scale / p1_scale) * ncm_matrix_get (H, i, i) -
	                         (p1_scale / p2_scale) * ncm_matrix_get (H, j, j)
	                         )
	                  );
	  ncm_matrix_set (H, j, i, ncm_matrix_get (H, i, j));
	  //printf ("d2[%d %d] = % 20.15g\n", i, j, (diff - gsl_matrix_get (H, k, k) - gsl_matrix_get (H, l, l)) / 2.0);
	  ncm_mset_fparam_set (fit->mset, i, p1);
	  ncm_mset_fparam_set (fit->mset, j, p2);
	}
  }
}

/**
 * ncm_fit_numdiff_forward_leastsquares_J:
 * @fit: a #NcmFit
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_numdiff_forward_leastsquares_J (NcmFit *fit, NcmMatrix *J)
{
  gint i;
  guint free_params_len = ncm_mset_max_fparam_name (fit->mset);

  ncm_fit_leastsquares_f (fit, fit->f);

  for (i = 0; i < free_params_len; i++)
  {
	NcmVector *J_col_i = ncm_matrix_get_col (J, i);
	const gdouble p = ncm_mset_fparam_get (fit->mset, i);
	const gdouble p_scale = ncm_mset_fparam_get_scale (fit->mset, i);
	const gdouble h = (p_scale * NCM_FIT_NUMDIFF_SCALE);
	const gdouble pph = p + h;
	const gdouble one_h = 1.0 / h;
	gint k;
	ncm_mset_fparam_set (fit->mset, i, pph);
	ncm_fit_leastsquares_f (fit, J_col_i);
	for (k = 0; k < ncm_vector_len (fit->f); k++)
	  ncm_vector_set (J_col_i, k,
	                  (ncm_vector_get (J_col_i, k) - ncm_vector_get (fit->f, k)) * one_h
	                  );
	ncm_mset_fparam_set (fit->mset, i, p);
	ncm_vector_free (J_col_i);
  }
}

/**
 * ncm_fit_numdiff_forward_leastsquares_f_J:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_numdiff_forward_leastsquares_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  ncm_fit_numdiff_forward_leastsquares_J (fit, J);
  ncm_vector_memcpy (f, fit->f);
}

/**
 * ncm_fit_numdiff_central_leastsquares_J:
 * @fit: a #NcmFit
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_numdiff_central_leastsquares_J (NcmFit *fit, NcmMatrix *J)
{
  gint i;
  guint free_params_len = ncm_mset_max_fparam_name (fit->mset);

  for (i = 0; i < free_params_len; i++)
  {
	  NcmVector *J_col_i = ncm_matrix_get_col (J, i);
	  const gdouble p = ncm_mset_fparam_get (fit->mset, i);
	  const gdouble p_scale = ncm_mset_fparam_get_scale (fit->mset, i);
	  const gdouble h = (p_scale * NCM_FIT_NUMDIFF_SCALE);
	  const gdouble pph = p + h;
	  const gdouble pmh = p - h;
	  const gdouble one_2h = 1.0 / (2.0 * h);
	  gint k;

	  ncm_mset_fparam_set (fit->mset, i, pph);
	  ncm_fit_leastsquares_f (fit, J_col_i);

	  ncm_mset_fparam_set (fit->mset, i, pmh);
	  ncm_fit_leastsquares_f (fit, fit->f);

	  for (k = 0; k < ncm_vector_len (fit->f); k++)
		ncm_vector_set (J_col_i, k,
		                (ncm_vector_get (J_col_i, k) - ncm_vector_get (fit->f, k)) * one_2h
		                );
	  ncm_mset_fparam_set (fit->mset, i, p);
	  ncm_vector_free (J_col_i);
  }
}

/**
 * ncm_fit_numdiff_central_leastsquares_f_J:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_numdiff_central_leastsquares_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  ncm_fit_leastsquares_f (fit, f);
  ncm_fit_numdiff_central_leastsquares_J (fit, J);
}

/**
 * ncm_fit_numdiff_m2lnL_covar:
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_fit_numdiff_m2lnL_covar (NcmFit *fit)
{
  gint ret;
  ncm_fit_numdiff_m2lnL_hessian (fit, fit->hessian);
  ncm_matrix_memcpy (fit->covar, fit->hessian);
  ncm_matrix_scale (fit->covar, 0.5);
  ret = gsl_linalg_cholesky_decomp (NCM_MATRIX_GSL (fit->covar));
  NC_TEST_GSL_RESULT("ncm_fit_numdiff_m2lnL_covar[gsl_linalg_cholesky_decomp]", ret);
  ret = gsl_linalg_cholesky_invert (NCM_MATRIX_GSL (fit->covar));
  NC_TEST_GSL_RESULT("ncm_fit_numdiff_m2lnL_covar[gsl_linalg_cholesky_invert]", ret);
}

/**
 * ncm_fit_m2lnL_val_grad:
 * @fit: a #NcmFit
 * @result: (out): FIXME
   * @df: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_m2lnL_val_grad (NcmFit *fit, gdouble *result, NcmVector *df)
{
  fit->n_func_eval++;
  fit->n_grad_eval++;

  nc_likelihood_m2lnL_val_grad (fit->lh, fit->mset, result, df);

  return;
}

typedef struct _FitDProb
{
  NcmMSetPIndex pi;
  NcmFit *fit_val;
  NcmFit *fit;
} FitDProb;

static gdouble
fit_dprob(gdouble val, gpointer p)
{
  FitDProb *dprob_arg = (FitDProb *)p;
  ncm_mset_param_set (dprob_arg->fit_val->mset, dprob_arg->pi.gmid, dprob_arg->pi.pid, val);
  ncm_fit_run (dprob_arg->fit_val, NC_BF_MAX_ITER, NCM_FIT_RUN_MSGS_NONE);
  if (dprob_arg->fit->mtype > NCM_FIT_RUN_MSGS_NONE)
	g_message (".");
  return exp (-(dprob_arg->fit_val->m2lnL - dprob_arg->fit->m2lnL) / 2.0);
}

/**
 * ncm_fit_prob:
 * @fit: a #NcmFit.
 * @gmid: a #NcmModelID.
 * @pid: FIXME
 * @a: FIXME
 * @b: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_prob (NcmFit *fit, NcmModelID gmid, guint pid, gdouble a, gdouble b)
{
  NcmFit *fit_val;
  NcmMSet *mset_val = ncm_mset_copy_all (fit->mset);
  gsl_integration_workspace **w;
  FitDProb dprob_arg;
  gdouble result, error;
  gint error_code;
  gsl_function F;

  ncm_mset_param_set_ftype (mset_val, gmid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_new (fit->lh, mset_val, fit->type, fit->gtype);

  dprob_arg.pi.gmid = gmid;
  dprob_arg.pi.pid  = pid;
  dprob_arg.fit     = fit;
  dprob_arg.fit_val = fit_val;

  F.function = &fit_dprob;
  F.params = &dprob_arg;

  w = nc_integral_get_workspace();
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
	g_message ("#");
  error_code = gsl_integration_qags (&F, a, b, 1e-10, NC_INT_ERROR, NC_INT_PARTITION, *w, &result, &error);
  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
	g_message ("\n");
  NC_TEST_GSL_RESULT ("ncm_fit_prob", error_code);
  ncm_memory_pool_return (w);

  ncm_fit_free (fit_val);
  ncm_mset_free (mset_val);
  return result;
}

/**
 * ncm_fit_dprob:
 * @fit: a #NcmFit.
 * @gmid: a #NcmModelID.
 * @pid: FIXME
 * @a: FIXME
 * @b: FIXME
 * @step: FIXME
 * @norm: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
ncm_fit_dprob (NcmFit *fit, NcmModelID gmid, guint pid, gdouble a, gdouble b, gdouble step, gdouble norm)
{
  NcmFit *fit_val;
  NcmMSet *mset_val = ncm_mset_copy_all (fit->mset);
  FitDProb dprob_arg;
  gdouble point;

  ncm_mset_param_set_ftype (mset_val, gmid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_new (fit->lh, mset_val, fit->type, fit->gtype);

  dprob_arg.pi.gmid = gmid;
  dprob_arg.pi.pid  = pid;
  dprob_arg.fit     = fit;
  dprob_arg.fit_val = fit_val;

  for (point = a; point <= b; point += step)
  {
	g_message ("%g %g\n", point, fit_dprob (point, &dprob_arg) / norm);
  }

  ncm_fit_free (fit_val);
  ncm_mset_free (mset_val);

  return;
}

/**
 * ncm_fit_lr_test_range:
 * @fit: a #NcmFit
 * @gmid: a #NcmModelID.
 * @pid: FIXME
 * @start: FIXME
 * @stop: FIXME
 * @step: FIXME
 *
 * FIXME
 */
void
ncm_fit_lr_test_range (NcmFit *fit, NcmModelID gmid, guint pid, gdouble start, gdouble stop, gdouble step)
{
  NcmFit *fit_val;
  NcmMSet *mset_val = ncm_mset_copy_all (fit->mset);
  gdouble walk;

  ncm_mset_param_set_ftype (mset_val, gmid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_new (fit->lh, mset_val, fit->type, fit->gtype);

  for (walk = start; walk <= stop; walk += step)
  {
	ncm_mset_param_set (fit_val->mset, gmid, pid, walk);
	ncm_fit_run (fit_val, NC_BF_MAX_ITER, NCM_FIT_RUN_MSGS_NONE);
	//g_message ("%g %g %g\n", walk, (fit_val->m2lnL - fit->m2lnL), gsl_cdf_chisq_Q (fit_val->m2lnL - fit->m2lnL, 1));
	g_message ("%g %g %g %g %g\n", walk,
	           (fit_val->m2lnL - fit->m2lnL),
	           gsl_ran_chisq_pdf (fit_val->m2lnL - fit->m2lnL, 1),
	           gsl_cdf_chisq_Q (fit_val->m2lnL - fit->m2lnL, 1),
	           gsl_cdf_ugaussian_Q (sqrt(fit_val->m2lnL - fit->m2lnL))
	           );
  }

  ncm_fit_free (fit_val);
  ncm_mset_free (mset_val);
  return;
}

/**
 * ncm_fit_lr_test:
 * @fit: a #NcmFit.
 * @gmid: a #NcmModelID.
 * @pid: FIXME
 * @val: FIXME
 * @dof: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_lr_test (NcmFit *fit, NcmModelID gmid, guint pid, gdouble val, gint dof)
{
  NcmFit *fit_val;
  NcmMSet *mset_val = ncm_mset_copy_all (fit->mset);
  gdouble result;

  ncm_mset_param_set_ftype (mset_val, gmid, pid, NCM_PARAM_TYPE_FIXED);
  fit_val = ncm_fit_new (fit->lh, mset_val, fit->type, fit->gtype);

  ncm_mset_param_set (fit_val->mset, gmid, pid, val);
  ncm_fit_run (fit_val, NC_BF_MAX_ITER, NCM_FIT_RUN_MSGS_NONE);
  result = gsl_cdf_chisq_Q (fit_val->m2lnL - fit->m2lnL, dof);

  ncm_fit_free (fit_val);
  ncm_mset_free (mset_val);
  return result;
}

/*
 * ncm_fit_chisq_test:
 * @fit: a #NcmFit
 * @bins: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 *
 gdouble
 ncm_fit_chisq_test (NcmFit *fit, size_t bins)
 {
   NcLikelihood *lh = fit->lh;
   NcDataSet *ds = lh->ds;
   NcData *data = nc_dataset_get_data (ds, 0);
   NcFunction distance = nc_function_get (data->res_type);
   gsl_histogram *obs = gsl_histogram_alloc (bins);
   gsl_histogram *exp = gsl_histogram_alloc (bins);
   gsl_vector_view y_exp;
   gdouble y_min, y_max;
   gint i;
   if (data == NULL)
   return GSL_NAN;

   ncm_distance_thread_pool (distance.f, fit->cp, data->x, fit->f, NULL, NULL);
   y_exp = gsl_vector_subvector (fit->f, 0, data->x->size);
   y_min = GSL_MIN (gsl_vector_min (data->y), gsl_vector_min (&y_exp.vector))*0.999;
   y_max = GSL_MAX (gsl_vector_max (data->y), gsl_vector_max (&y_exp.vector))*1.001;

   gsl_histogram_set_ranges_uniform (obs, y_min, y_max);
   gsl_histogram_set_ranges_uniform (exp, y_min, y_max);

   for (i = 0; i < data->x->size; i++)
   {
	 gsl_histogram_increment (obs, gsl_vector_get (data->y, i));
	 gsl_histogram_increment (exp, gsl_vector_get (&y_exp.vector, i));
	 }

	 gsl_histogram_sub (obs, exp);
	 //  if (gsl_histogram_min_val (exp) == 0.0)
	 //    g_error ("Cannot chisq test");
	 gsl_histogram_mul (obs, obs);
	 gsl_histogram_div (obs, exp);
	 //gsl_histogram_fprintf (fit->log, obs, "%f", "%f");
	 g_message ("\n");
	 //gsl_histogram_fprintf (fit->log, exp, "%f", "%f");
	 g_message ("\n");
	 g_message ("BLA [%g, %g] %zu %g [%g]\n", y_min, y_max, bins, gsl_histogram_sum (obs), gsl_histogram_sum (exp));
	 return gsl_histogram_sum (obs);
	 }
	 */

/**
 * ncm_fit_leastsquares_f:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 *
 * FIXME
 */
void
ncm_fit_leastsquares_f (NcmFit *fit, NcmVector *f)
{
  if (fit->bootstrap)
  {
	guint i;
	nc_likelihood_leastsquares_f (fit->lh, fit->mset, fit->f);
	for (i = 0; i < ncm_vector_len (f); i++)
	  ncm_vector_set (f, i, ncm_vector_get (fit->f, gsl_vector_uint_get (fit->bs, i)));
  }
  else
	nc_likelihood_leastsquares_f (fit->lh, fit->mset, f);
  fit->n_func_eval++;
  return;
}

/**
 * ncm_fit_leastsquares_J:
 * @fit: a #NcmFit
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_leastsquares_J (NcmFit *fit, NcmMatrix *J)
{
  nc_likelihood_leastsquares_J (fit->lh, fit->mset, J);

  if (fit->bootstrap)
  {
	guint i;
	nc_likelihood_leastsquares_J (fit->lh, fit->mset, fit->J);
	for (i = 0; i < NCM_MATRIX_NROWS (J); i++)
	{
	  NcmVector *fit_J_row = ncm_matrix_get_row (fit->J, gsl_vector_uint_get (fit->bs, i));
	  NcmVector *J_row = ncm_matrix_get_row (J, i);
	  ncm_vector_memcpy (J_row, fit_J_row);
	  ncm_vector_free (fit_J_row);
	  ncm_vector_free (J_row);
	}
  }
  fit->n_grad_eval++;

  return;
}

/**
 * ncm_fit_leastsquares_f_J:
 * @fit: a #NcmFit
 * @f: a #NcmVector
 * @J: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_leastsquares_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  nc_likelihood_leastsquares_f_J (fit->lh, fit->mset, f, J);

  if (fit->bootstrap)
  {
	guint i;
	nc_likelihood_leastsquares_J (fit->lh, fit->mset, fit->J);
	for (i = 0; i < NCM_MATRIX_NROWS (J); i++)
	{
	  NcmVector *fit_J_row = ncm_matrix_get_row (fit->J, gsl_vector_uint_get (fit->bs, i));
	  NcmVector *J_row = ncm_matrix_get_row (J, i);
	  ncm_vector_memcpy (J_row, fit_J_row);
	  ncm_vector_free (fit_J_row);
	  ncm_vector_free (J_row);
	}
  }
  fit->n_func_eval++;
  fit->n_grad_eval++;

  return;
}

/**
 * ncm_fit_function_error:
 * @fit: a #NcmFit
 * @func: a #NcmMSetFunc
 * @z: FIXME
 * @pretty_print: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_function_error (NcmFit *fit, NcmMSetFunc *func, gdouble z, gboolean pretty_print)
{
  g_assert_not_reached ();
  /*
  gdouble result;
  gint ret;
  guint free_params_len = ncm_mset_max_fparam_name (fit->mset);
  NcmVector *v = ncm_vector_new (free_params_len);

  ncm_mset_func_numdiff_fparams (func, );
  NCM_FUNC_DF (func, fit->mset, fit->pt, z, fit->df);

  ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (fit->covar), ncm_vector_gsl (fit->df), 0.0, ncm_vector_gsl (tmp1));
  NC_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
  ret = gsl_blas_ddot (ncm_vector_gsl (fit->df), ncm_vector_gsl (tmp1), &result);
  NC_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);

  if (pretty_print)
	g_message ("# % -12.4g +/- % -12.4g\n", NCM_FUNC_F(func, fit->mset, z), sqrt(result));

  ncm_fit_params_return_tmp_vector (fit->pt, tmp1);

  return sqrt(result);
  */
}

/**
 * ncm_fit_function_cov:
 * @fit: a #NcmFit
 * @func1: a #NcmMSetFunc
 * @z1: FIXME
 * @func2: a #NcmMSetFunc
 * @z2: FIXME
 * @pretty_print: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_function_cov (NcmFit *fit, NcmMSetFunc *func1, gdouble z1, NcmMSetFunc *func2, gdouble z2, gboolean pretty_print)
{
  g_assert_not_reached ();
  /*
  gdouble result, cor, s1, s2;
  gint ret;
  NcmVector *tmp1 = ncm_fit_params_get_tmp_vector (fit->pt, fit->mset);
  NcmVector *tmp2 = ncm_fit_params_get_tmp_vector (fit->pt, fit->mset);

  NCM_FUNC_DF (func1, fit->mset, fit->pt, z1, tmp1);
  NCM_FUNC_DF (func2, fit->mset, fit->pt, z2, tmp2);

  ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (fit->covar), ncm_vector_gsl (tmp1), 0.0, ncm_vector_gsl (fit->df));
  NC_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
  ret = gsl_blas_ddot (ncm_vector_gsl (fit->df), ncm_vector_gsl (tmp1), &result);
  NC_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);
  s1 = sqrt(result);

  ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (fit->covar), ncm_vector_gsl (tmp2), 0.0, ncm_vector_gsl (fit->df));
  NC_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
  ret = gsl_blas_ddot (ncm_vector_gsl (fit->df), ncm_vector_gsl (tmp2), &result);
  NC_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);
  s2 = sqrt(result);

  ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (fit->covar), ncm_vector_gsl (tmp1), 0.0, ncm_vector_gsl (fit->df));
  NC_TEST_GSL_RESULT("ncm_fit_function_error[covar.v]", ret);
  ret = gsl_blas_ddot (ncm_vector_gsl (fit->df), ncm_vector_gsl (tmp2), &result);
  NC_TEST_GSL_RESULT("ncm_fit_function_error[v.covar.v]", ret);
  cor = result / (s1*s2);

  if (pretty_print)
  {
	g_message ("# % -12.4g\t| % -12.4g\n", NCM_FUNC_F(func1, fit->mset, z1), NCM_FUNC_F(func2, fit->mset, z2));
	g_message ("#---------------------------------------------\n");
	g_message ("# % -12.4g\t | % -12.4g\t| % -12.4g\n", s1, 1.0, cor);
	g_message ("# % -12.4g\t | % -12.4g\t| % -12.4g\n", s2, cor, 1.0);
	g_message ("#---------------------------------------------\n");
  }

  ncm_fit_params_return_tmp_vector (fit->pt, tmp1);
  ncm_fit_params_return_tmp_vector (fit->pt, tmp2);

  return sqrt(result);
  */
}

/**
 * ncm_fit_montecarlo_matrix:
 * @fit: a #NcmFit.
 * @mset: a #NcmMSet.
 * @maxiter: FIXME
 * @ni: FIXME
 * @nf: FIXME
 * @mtype: a #NcmFitRunMsgs.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMatrix *
ncm_fit_montecarlo_matrix (NcmFit *fit, NcmMSet *mset, guint maxiter, guint ni, guint nf, NcmFitRunMsgs mtype)
{
  const guint n = nf - ni;
  guint free_params_len = ncm_mset_fparams_len (fit->mset);
  NcmMatrix *param_matrix = ncm_matrix_new (n, free_params_len);
  NcmVector *bf = ncm_vector_new (free_params_len);
  gint i;
  gdouble mean_time = 0.0, total_sec;
  gulong total_min, total_hour, total_day;
  g_assert (nf > ni);

  ncm_mset_prepare_fparam_map (mset);
  ncm_mset_fparams_get_vector (mset, bf);
  ncm_mset_fparams_set_vector (fit->mset, bf);

  if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
	g_message ( "#  Calculating [%06d] montecarlo fit\n", n);
  for (i = 0; i < ni; i++)
	nc_dataset_resample (fit->lh->ds, fit->mset, FALSE);
  for (; i < nf; i++)
  {
	NcmVector *bf_i = ncm_matrix_get_row (param_matrix, i - ni);

	ncm_mset_fparams_set_vector (fit->mset, bf);
	nc_dataset_resample (fit->lh->ds, fit->mset, FALSE);
	ncm_fit_run (fit, maxiter, mtype);
	ncm_mset_fparams_get_vector (fit->mset, bf_i);

	if (fit->mtype > NCM_FIT_RUN_MSGS_NONE)
	{
	  gdouble elap_sec = g_timer_elapsed (fit->timer, NULL);
	  gulong elap_min = elap_sec / 60;
	  gulong elap_hour = elap_min / 60;
	  gulong elap_day = elap_hour / 24;

	  mean_time = (i * mean_time + elap_sec) / (i + 1.0);
	  elap_sec = fmod (elap_sec, 60);
	  elap_min = elap_min % 60;
	  elap_hour = elap_hour % 24;

	  total_sec = mean_time * (n - 1 - i);
	  total_min = total_sec / 60;
	  total_hour = total_min / 60;
	  total_day = total_hour / 24;
	  total_sec = fmod (total_sec, 60);
	  total_min = total_min % 60;
	  total_hour = total_hour % 24;

	  g_message ( "# finish[%d] took: %02lu days, %02lu:%02lu:%010.7f | %02lu days, %02lu:%02lu:%010.7f left\n", i, elap_day, elap_hour, elap_min, elap_sec,
	             total_day, total_hour, total_min, total_sec);
	}
	ncm_vector_free (bf_i);
  }
  ncm_vector_free (bf);
  return param_matrix;
}

/**
 * ncm_fit_montecarlo_matrix_print:
 * @fit: a #NcmFit
 * @param_matrix: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_montecarlo_matrix_print (NcmFit *fit, NcmMatrix *param_matrix)
{
  gint i;

  for (i = 0; i < NCM_MATRIX_NROWS (param_matrix); i++)
  {
	gint j;
	for (j = 0; j < NCM_MATRIX_NCOLS (param_matrix); j++)
	  g_message ("% -20.15g ", ncm_matrix_get (param_matrix, i, j));
	g_message ("\n");
  }
}

/**
 * ncm_fit_montecarlo_matrix_mean_covar:
 * @fit: a #NcmFit
 * @param_matrix: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_fit_montecarlo_matrix_mean_covar (NcmFit *fit, NcmMatrix *param_matrix)
{
  gint i;
  for (i = 0; i < fit->fparam_len; i++)
  {
	gint j;
	NcmVector *p_i = ncm_matrix_get_col (param_matrix, i);
	ncm_vector_set (fit->x, i, gsl_stats_mean (ncm_vector_gsl (p_i)->data, ncm_vector_gsl (p_i)->stride, ncm_vector_gsl (p_i)->size));
	ncm_matrix_set (fit->covar, i, i, gsl_stats_variance (ncm_vector_gsl (p_i)->data, ncm_vector_gsl (p_i)->stride, ncm_vector_gsl (p_i)->size));
	for (j = i + 1; j < fit->fparam_len; j++)
	{
	  NcmVector *p_j = ncm_matrix_get_col (param_matrix, j);
	  ncm_matrix_set (fit->covar, i, j,
	                  gsl_stats_covariance (ncm_vector_gsl (p_i)->data, ncm_vector_gsl (p_i)->stride, ncm_vector_gsl (p_j)->data, ncm_vector_gsl (p_j)->stride, ncm_vector_gsl (p_i)->size)
	                  );
	  ncm_matrix_set (fit->covar, j, i, ncm_matrix_get (fit->covar, i, j));
	  ncm_vector_free (p_j);
	}
	ncm_vector_free (p_i);
  }
  ncm_mset_fparams_set_vector (fit->mset, fit->x);
}


static void
ncm_fit_init (NcmFit *fit)
{
  fit->covar   = NULL;
  fit->hessian = NULL;
  fit->x       = NULL;
  fit->jackpv  = NULL;
  fit->df      = NULL;
  fit->J       = NULL;
  fit->bootstrap = FALSE;
  fit->sqrt_m2lnL = 0.0;
  fit->m2lnL = 0.0;
  fit->m2lnL_dof = 0.0;

  fit->least_squares = NULL;
  fit->multimin = NULL;
  fit->simplex = NULL;
  fit->levmar = NULL;
#ifdef HAVE_NLOPT_2_2
  fit->nlopt = NULL;
  fit->nlopt_algo = 0;
#endif /* HAVE_NLOPT_2_2 */

  fit->niter = 0;
  fit->n_func_eval = 0;
  fit->n_grad_eval = 0;
  fit->m2lnL_prec = NCM_FIT_DEFAULT_M2LNL_RELTOL * 0.0 + 1e-8;
  fit->m2lnL_prec_target = NCM_FIT_DEFAULT_M2LNL_RELTOL * 0.0 + 1e-8;
  fit->params_prec = NC_HICOSMO_DEFAULT_PARAMS_RELTOL * 0.0 + 1e-8;
  fit->params_prec_target = NC_HICOSMO_DEFAULT_PARAMS_RELTOL * 0.0 + 1e-8;
  fit->timer = g_timer_new ();
  fit->solver_name = "[[No solver name]]";
  fit->mtype = NCM_FIT_RUN_MSGS_NONE;
}

static void
ncm_fit_dispose (GObject *object)
{
  NcmFit *fit = NCM_FIT (object);

  nc_likelihood_free (fit->lh);
  ncm_mset_free (fit->mset);

  if (fit->fparam_len > 0)
  {
	ncm_matrix_free (fit->covar);
	ncm_matrix_free (fit->hessian);
	ncm_vector_free (fit->x);
	ncm_vector_free (fit->df);
	ncm_matrix_free (fit->J);
  }
  ncm_vector_free (fit->f);
  ncm_vector_free (fit->jackpv);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_parent_class)->dispose (object);
}

static void
ncm_fit_finalize (GObject *object)
{
  NcmFit *fit = NCM_FIT (object);
  if (fit->least_squares != NULL)
	gsl_multifit_fdfsolver_free (fit->least_squares);

  if (fit->multimin != NULL)
	gsl_multimin_fdfminimizer_free (fit->multimin);

  if (fit->simplex != NULL)
	gsl_multimin_fminimizer_free (fit->simplex);

#ifdef NUMCOSMO_HAVE_NLOPT
#ifdef HAVE_NLOPT_2_2
  if (fit->nlopt)
  {
	nlopt_destroy (fit->nlopt);
  }
#endif
#endif

  gsl_vector_uint_free (fit->bs);
  g_timer_destroy (fit->timer);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_parent_class)->finalize (object);
}

static void
ncm_fit_class_init (NcmFitClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose  = ncm_fit_dispose;
  object_class->finalize = ncm_fit_finalize;
}

