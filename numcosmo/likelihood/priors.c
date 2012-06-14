/***************************************************************************
 *            priors.c
 *
 *  Wed Mar 19 12:46:46 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:priors
 * @title: Statistical Priors
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_sf_exp.h>

/***************************************************************************
 * Gaussian Prior NcFunction
 *
 *
 ****************************************************************************/

static NcGaussianPrior *
_nc_prior_gauss_new (void)
{
  NcGaussianPrior *gauss = g_slice_new (NcGaussianPrior);
  gauss->func = NULL;
  gauss->pi.gmid = 0;
  gauss->pi.pid = 0;
  gauss->z = GSL_NAN;
  gauss->mean = GSL_NAN;
  gauss->sigma = GSL_NAN;
  return gauss;
}

static void
_nc_prior_gauss_free (gpointer p)
{
  NcGaussianPrior *gauss = (NcGaussianPrior *)p;
  if (gauss->func)
	ncm_mset_func_free (gauss->func);
  g_slice_free (NcGaussianPrior, gauss);
}

static void
gaussian_prior_func_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcGaussianPrior *gp = (NcGaussianPrior *)obj;
  f[0] = (ncm_mset_func_eval1 (gp->func, mset, gp->z) - gp->mean) / gp->sigma;
}

/**
 * FIXME
 */
gboolean
nc_prior_add_gaussian_func (NcLikelihood *lh, NcmMSetFunc *func, gdouble z, gdouble mean, gdouble sigma)
{
  NcGaussianPrior *gp = _nc_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&gaussian_prior_func_f, 0, 1, gp, _nc_prior_gauss_free);

  gp->func = func;
  gp->z = z;
  gp->mean = mean;
  gp->sigma = sigma;

  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

static void
gaussian_prior_func0_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcGaussianPrior *gp = (NcGaussianPrior *)obj;
  f[0] = (ncm_mset_func_eval0 (gp->func, mset) - gp->mean) / gp->sigma;
}

/**
 * FIXME
 */
gboolean
nc_prior_add_gaussian_const_func (NcLikelihood *lh, NcmMSetFunc *func, gdouble mean, gdouble sigma)
{
  NcGaussianPrior *gp = _nc_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&gaussian_prior_func0_f, 0, 1, gp, _nc_prior_gauss_free);
  gp->func = func;
  gp->mean = mean;
  gp->sigma = sigma;

  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

/**
 * FIXME
 */
gboolean
nc_prior_add_gaussian (NcLikelihood *lh, NcGaussianPrior *gp)
{
  NcmMSetFunc *prior = ncm_mset_func_new (&gaussian_prior_func_f, 0, 1, gp, _nc_prior_gauss_free);
  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

static void
gaussian_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcGaussianPrior *gp = (NcGaussianPrior *)obj;
  gdouble mean = gp->mean;
  gdouble sigma = gp->sigma;
  f[0] = (ncm_mset_param_get (mset, gp->pi.gmid, gp->pi.pid) - mean) / sigma;
}

/**
 * FIXME
 */
gboolean
nc_prior_add_gaussian_data (NcLikelihood *lh, NcmModelID gmid, guint pid, gdouble mean, gdouble sigma)
{
  NcGaussianPrior *gp = _nc_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&gaussian_prior_f, 0, 1, gp, _nc_prior_gauss_free);
  gp->pi.gmid = gmid;
  gp->pi.pid = pid;
  gp->mean = mean;
  gp->sigma = sigma;

  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

static void
positive_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcGaussianPrior *gp = (NcGaussianPrior *)obj;
  f[0] = ncm_mset_param_get (mset, gp->pi.gmid, gp->pi.pid) > 0.0 ? 0.0 : GSL_POSINF;
}

/**
 * FIXME
 */
gboolean
nc_prior_add_positive (NcLikelihood *lh, NcmModelID gmid, guint pid)
{
  NcGaussianPrior *gp = _nc_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&positive_prior_f, 0, 1, gp, _nc_prior_gauss_free);
  gp->pi.gmid = gmid;
  gp->pi.pid = pid;
  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

#define HUGE_EXPONENT_NUMBER (20)

static void
oneside_a_inf_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcGaussianPrior *gp = (NcGaussianPrior *)obj;
  gdouble a = gp->mean;
  gdouble s = gp->sigma;
  gdouble p = ncm_mset_param_get (mset, gp->pi.gmid, gp->pi.pid);
  f[0] = exp (2.0 * HUGE_EXPONENT_NUMBER / s * ((a - p) + s / 2.0));
}

/**
 * FIXME
 */
gboolean
nc_prior_add_oneside_a_inf_param (NcLikelihood *lh, NcmModelID gmid, guint pid, gdouble a, gdouble s)
{
  NcGaussianPrior *gp = _nc_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&oneside_a_inf_prior_f, 0, 1, gp, _nc_prior_gauss_free);
  gp->pi.gmid = gmid;
  gp->pi.pid = pid;
  gp->mean = a;
  gp->sigma = s;

  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

static void
oneside_a_inf_prior_func_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcGaussianPrior *gp = (NcGaussianPrior *)obj;
  gdouble a = gp->mean;
  gdouble s = gp->sigma;
  gdouble p = ncm_mset_func_eval1 (gp->func, mset, gp->z);
  f[0] = exp (2.0 * HUGE_EXPONENT_NUMBER / s * ((a - p) + s / 2.0));
}

/**
 * FIXME
 */
gboolean
nc_prior_add_oneside_a_inf_func (NcLikelihood *lh, NcmMSetFunc *func, gdouble z, gdouble mean, gdouble s)
{
  NcGaussianPrior *gp = _nc_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&oneside_a_inf_prior_func_f, 0, 1, gp, _nc_prior_gauss_free);
  gp->func = func;
  gp->z = z;
  gp->mean = mean;
  gp->sigma = s;

  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

static void
oneside_a_inf_prior_const_func_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcGaussianPrior *gp = (NcGaussianPrior *)obj;
  gdouble a = gp->mean;
  gdouble s = gp->sigma;
  gdouble p = ncm_mset_func_eval0 (gp->func, mset);
  f[0] = exp (2.0 * HUGE_EXPONENT_NUMBER / s * ((a - p) + s / 2.0));
}

/**
 * FIXME
 */
gboolean
nc_prior_add_oneside_a_inf_const_func (NcLikelihood *lh, NcmMSetFunc *func, gdouble mean, gdouble s)
{
  NcGaussianPrior *gp = _nc_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&oneside_a_inf_prior_const_func_f, 0, 1, gp, _nc_prior_gauss_free);
  gp->func = func;
  gp->mean = mean;
  gp->sigma = s;

  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

/**
 * FIXME
 */
gboolean
nc_prior_add_oneside_a_inf (NcLikelihood *lh, NcGaussianPrior *gp)
{
  NcmMSetFunc *prior = ncm_mset_func_new (&oneside_a_inf_prior_func_f, 0, 1, gp, _nc_prior_gauss_free);
  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

static void
twoside_a_b_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcGaussianPrior *gp = (NcGaussianPrior *)obj;
  gdouble a = gp->mean;
  gdouble b = gp->z;
  gdouble s = gp->sigma;
  gdouble p = ncm_mset_param_get (mset, gp->pi.gmid, gp->pi.pid);
  f[0] = exp( 2.0 * HUGE_EXPONENT_NUMBER / s * ((a - p) + s / 2.0) )+
	exp( 2.0 * HUGE_EXPONENT_NUMBER / s * ((p - b) + s / 2.0) );
}

/**
 * FIXME
 */
gboolean
nc_prior_add_twoside_a_b (NcLikelihood *lh, NcmModelID gmid, guint pid, gdouble a, gdouble b, gdouble s)
{
  NcGaussianPrior *gp = _nc_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&twoside_a_b_prior_f, 0, 1, gp, _nc_prior_gauss_free);
  gp->pi.gmid = gmid;
  gp->pi.pid = pid;
  gp->mean = a;
  gp->z = b;
  gp->sigma = s;

  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}

static void
topological_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcTopologicalPrior *tp = (NcTopologicalPrior *)obj;
  NcHICosmo *model = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  gdouble Omega_k = nc_hicosmo_Omega_k (model);
  gdouble sqrt_Omega_k = sqrt (fabs (Omega_k));
  gint k = fabs (Omega_k) < NC_ZERO_LIMIT ? 0 : (Omega_k > 0.0 ? -1 : 1);
  gdouble z = tp->z;
  gdouble mean = tp->mean;
  gdouble sigma = tp->sigma;
  gdouble cd = nc_distance_comoving (tp->dist, model, z);
  if (!gsl_finite (cd))
	f[0] = GSL_POSINF;
  else if (k <= 0)
	f[0] = GSL_POSINF;
  else
	f[0] = (sqrt_Omega_k * cd - mean) / sigma;
}

/**
 * nc_prior_topological_new: (skip)
 * @z: FIXME
 * @alpha: FIXME
 * @sigma_alpha: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcTopologicalPrior *
nc_prior_topological_new (gdouble z, gdouble alpha, gdouble sigma_alpha, gint n)
{
  NcTopologicalPrior *tp = g_slice_new (NcTopologicalPrior);
  nc_prior_topological_set (tp, z, alpha, sigma_alpha, n);
  return tp;
}

/**
 * nc_prior_topological_free:
 * @tp: a #NcTopologicalPrior
 *
 * FIXME
 */
void
nc_prior_topological_free (NcTopologicalPrior *tp)
{
  nc_distance_free (tp->dist);
  g_slice_free (NcTopologicalPrior, tp);
  return;
}

/**
 * NcTopologicalPrior:
 * @tp: a #NcTopologicalPrior
 * @z: FIXME
 * @alpha: FIXME
 * @sigma_alpha: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_prior_topological_set (NcTopologicalPrior *tp, gdouble z, gdouble alpha, gdouble sigma_alpha, gint n)
{
  tp->z = z;
  tp->mean = atan(tan(M_PI/n)/cos(alpha));
  tp->sigma = sin(alpha)*tan(M_PI/n) / (cos(alpha)*cos(alpha) + tan(M_PI/n)*tan(M_PI/n));
  tp->sigma = sqrt(tp->sigma*tp->sigma*sigma_alpha*sigma_alpha);
  if (n == 2)
	tp->sigma = sigma_alpha / alpha;
  return TRUE;
}

/**
 * nc_prior_add_topological:
 * @lh: a #NcLikelihood
 * @z: FIXME
 * @alpha: FIXME
 * @sigma_alpha: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_prior_add_topological (NcLikelihood *lh, gdouble z, gdouble alpha, gdouble sigma_alpha, gint n)
{
  NcTopologicalPrior *tp = g_slice_new (NcTopologicalPrior);
  NcmMSetFunc *prior = ncm_mset_func_new (&topological_prior_f, 0, 1, tp, (GDestroyNotify)nc_prior_topological_free);
  nc_prior_topological_set (tp, z, alpha, sigma_alpha, n);
  lh->priors = g_list_append (lh->priors, prior);
  return TRUE;
}
