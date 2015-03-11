/***************************************************************************
 *            ncm_priors.c
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
 * SECTION:ncm_priors
 * @title: NcmPrior
 * @short_description: General statistical priors.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_priors.h"

static NcmPriorGauss *
_ncm_prior_gauss_new (void)
{
  NcmPriorGauss *gauss = g_slice_new (NcmPriorGauss);
  gauss->func = NULL;
  gauss->pi.mid = 0;
  gauss->pi.pid = 0;
  gauss->z = GSL_NAN;
  gauss->mean = GSL_NAN;
  gauss->sigma = GSL_NAN;
  return gauss;
}

static void
_ncm_prior_gauss_free (gpointer p)
{
  NcmPriorGauss *gauss = (NcmPriorGauss *)p;
  if (gauss->func)
    ncm_mset_func_free (gauss->func);
  g_slice_free (NcmPriorGauss, gauss);
}

static void
gaussian_prior_func_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmPriorGauss *gp = (NcmPriorGauss *)obj;
  NCM_UNUSED (x);
  f[0] = (ncm_mset_func_eval1 (gp->func, mset, gp->z) - gp->mean) / gp->sigma;
}

/**
 * ncm_prior_add_gaussian_func:
 * @lh: FIXME
 * @func: FIXME
 * @z: FIXME
 * @mean: FIXME
 * @sigma: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_gaussian_func (NcmLikelihood *lh, NcmMSetFunc *func, gdouble z, gdouble mean, gdouble sigma)
{
  NcmPriorGauss *gp = _ncm_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&gaussian_prior_func_f, 0, 1, gp, _ncm_prior_gauss_free);

  gp->func = ncm_mset_func_ref (func);
  gp->z = z;
  gp->mean = mean;
  gp->sigma = sigma;

  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

static void
gaussian_prior_func0_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmPriorGauss *gp = (NcmPriorGauss *)obj;
  NCM_UNUSED (x);
  f[0] = (ncm_mset_func_eval0 (gp->func, mset) - gp->mean) / gp->sigma;
}

/**
 * ncm_prior_add_gaussian_const_func:
 * @lh: FIXME
 * @func: FIXME
 * @mean: FIXME
 * @sigma: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_gaussian_const_func (NcmLikelihood *lh, NcmMSetFunc *func, gdouble mean, gdouble sigma)
{
  NcmPriorGauss *gp = _ncm_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&gaussian_prior_func0_f, 0, 1, gp, _ncm_prior_gauss_free);
  gp->func = ncm_mset_func_ref (func);
  gp->mean = mean;
  gp->sigma = sigma;

  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

/**
 * ncm_prior_add_gaussian:
 * @lh: FIXME
 * @gp: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_gaussian (NcmLikelihood *lh, NcmPriorGauss *gp)
{
  NcmMSetFunc *prior = ncm_mset_func_new (&gaussian_prior_func_f, 0, 1, gp, _ncm_prior_gauss_free);
  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

static void
gaussian_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmPriorGauss *gp = (NcmPriorGauss *)obj;
  gdouble mean = gp->mean;
  gdouble sigma = gp->sigma;
  NCM_UNUSED (x);
  f[0] = (ncm_mset_param_get (mset, gp->pi.mid, gp->pi.pid) - mean) / sigma;
}

/**
 * ncm_prior_add_gaussian_data:
 * @lh: FIXME
 * @mid: FIXME
 * @pid: FIXME
 * @mean: FIXME
 * @sigma: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_gaussian_data (NcmLikelihood *lh, NcmModelID mid, guint pid, gdouble mean, gdouble sigma)
{
  NcmPriorGauss *gp = _ncm_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&gaussian_prior_f, 0, 1, gp, _ncm_prior_gauss_free);
  gp->pi.mid = mid;
  gp->pi.pid = pid;
  gp->mean = mean;
  gp->sigma = sigma;

  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

static void
positive_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmPriorGauss *gp = (NcmPriorGauss *)obj;
  NCM_UNUSED (x);
  f[0] = ncm_mset_orig_param_get (mset, gp->pi.mid, gp->pi.pid) > 0.0 ? 0.0 : GSL_POSINF;
}

/**
 * ncm_prior_add_positive:
 * @lh: FIXME
 * @mid: FIXME
 * @pid: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_positive (NcmLikelihood *lh, NcmModelID mid, guint pid)
{
  NcmPriorGauss *gp = _ncm_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&positive_prior_f, 0, 1, gp, _ncm_prior_gauss_free);
  gp->pi.mid = mid;
  gp->pi.pid = pid;
  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

#define HUGE_EXPONENT_NUMBER (20)

static void
oneside_a_inf_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmPriorGauss *gp = (NcmPriorGauss *)obj;
  gdouble a = gp->mean;
  gdouble s = gp->sigma;
  gdouble p = ncm_mset_param_get (mset, gp->pi.mid, gp->pi.pid);
  NCM_UNUSED (x);
  f[0] = exp (2.0 * HUGE_EXPONENT_NUMBER / s * ((a - p) + s / 2.0));
}

/**
 * ncm_prior_add_oneside_a_inf_param:
 * @lh: FIXME
 * @mid: FIXME
 * @pid: FIXME
 * @a: FIXME
 * @s: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_oneside_a_inf_param (NcmLikelihood *lh, NcmModelID mid, guint pid, gdouble a, gdouble s)
{
  NcmPriorGauss *gp = _ncm_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&oneside_a_inf_prior_f, 0, 1, gp, _ncm_prior_gauss_free);
  gp->pi.mid = mid;
  gp->pi.pid = pid;
  gp->mean = a;
  gp->sigma = s;

  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

static void
oneside_a_inf_prior_func_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmPriorGauss *gp = (NcmPriorGauss *)obj;
  gdouble a = gp->mean;
  gdouble s = gp->sigma;
  gdouble p = ncm_mset_func_eval1 (gp->func, mset, gp->z);
  NCM_UNUSED (x);
  f[0] = exp (2.0 * HUGE_EXPONENT_NUMBER / s * ((a - p) + s / 2.0));
}

/**
 * ncm_prior_add_oneside_a_inf_func:
 * @lh: FIXME
 * @func: FIXME
 * @z: FIXME
 * @mean: FIXME
 * @s: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_oneside_a_inf_func (NcmLikelihood *lh, NcmMSetFunc *func, gdouble z, gdouble mean, gdouble s)
{
  NcmPriorGauss *gp = _ncm_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&oneside_a_inf_prior_func_f, 0, 1, gp, _ncm_prior_gauss_free);
  gp->func = ncm_mset_func_ref (func);
  gp->z = z;
  gp->mean = mean;
  gp->sigma = s;

  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

static void
oneside_a_inf_prior_const_func_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmPriorGauss *gp = (NcmPriorGauss *)obj;
  gdouble a = gp->mean;
  gdouble s = gp->sigma;
  gdouble p = ncm_mset_func_eval0 (gp->func, mset);
  NCM_UNUSED (x);
  f[0] = exp (2.0 * HUGE_EXPONENT_NUMBER / s * ((a - p) + s / 2.0));
}

/**
 * ncm_prior_add_oneside_a_inf_const_func:
 * @lh: FIXME
 * @func: FIXME
 * @mean: FIXME
 * @s: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_oneside_a_inf_const_func (NcmLikelihood *lh, NcmMSetFunc *func, gdouble mean, gdouble s)
{
  NcmPriorGauss *gp = _ncm_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&oneside_a_inf_prior_const_func_f, 0, 1, gp, _ncm_prior_gauss_free);
  gp->func = ncm_mset_func_ref (func);
  gp->mean = mean;
  gp->sigma = s;

  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

/**
 * ncm_prior_add_oneside_a_inf:
 * @lh: FIXME
 * @gp: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_oneside_a_inf (NcmLikelihood *lh, NcmPriorGauss *gp)
{
  NcmMSetFunc *prior = ncm_mset_func_new (&oneside_a_inf_prior_func_f, 0, 1, gp, _ncm_prior_gauss_free);
  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

static void
twoside_a_b_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmPriorGauss *gp = (NcmPriorGauss *)obj;
  gdouble a = gp->mean;
  gdouble b = gp->z;
  gdouble s = gp->sigma;
  gdouble p = ncm_mset_param_get (mset, gp->pi.mid, gp->pi.pid);
  NCM_UNUSED (x);
  f[0] = exp( 2.0 * HUGE_EXPONENT_NUMBER / s * ((a - p) + s / 2.0) )+
	exp( 2.0 * HUGE_EXPONENT_NUMBER / s * ((p - b) + s / 2.0) );
}

/**
 * ncm_prior_add_twoside_a_b:
 * @lh: FIXME
 * @mid: FIXME
 * @pid: FIXME
 * @a: FIXME
 * @b: FIXME
 * @s: FIXME
 *
 * FIXME
 * 
 */
void
ncm_prior_add_twoside_a_b (NcmLikelihood *lh, NcmModelID mid, guint pid, gdouble a, gdouble b, gdouble s)
{
  NcmPriorGauss *gp = _ncm_prior_gauss_new ();
  NcmMSetFunc *prior = ncm_mset_func_new (&twoside_a_b_prior_f, 0, 1, gp, _ncm_prior_gauss_free);
  gp->pi.mid = mid;
  gp->pi.pid = pid;
  gp->mean = a;
  gp->z = b;
  gp->sigma = s;

  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}
