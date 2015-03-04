/***************************************************************************
 *            nc_hicosmo_priors.c
 *
 *  Thu November 22 17:22:03 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * SECTION:nc_hicosmo_priors
 * @title: NcHICosmoPrior
 * @short_description: Collection of priors for NcHICosmo models.
 * 
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hicosmo_priors.h"

static void
_nc_hicosmo_prior_top_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcHICosmoPriorTop *tp = (NcHICosmoPriorTop *)obj;
  NcHICosmo *model = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  gdouble Omega_k = nc_hicosmo_Omega_k (model);
  gdouble sqrt_Omega_k = sqrt (fabs (Omega_k));
  gint k = fabs (Omega_k) < NCM_ZERO_LIMIT ? 0 : (Omega_k > 0.0 ? -1 : 1);
  gdouble z = tp->z;
  gdouble mean = tp->mean;
  gdouble sigma = tp->sigma;
  gdouble cd = nc_distance_comoving (tp->dist, model, z);

  NCM_UNUSED (x);
  
  if (!gsl_finite (cd))
    f[0] = GSL_POSINF;
  else if (k <= 0)
    f[0] = GSL_POSINF;
  else
    f[0] = (sqrt_Omega_k * cd - mean) / sigma;
}

/**
 * nc_hicosmo_prior_top_new: (skip)
 * @z: FIXME
 * @alpha: FIXME
 * @sigma_alpha: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoPriorTop *
nc_hicosmo_prior_top_new (gdouble z, gdouble alpha, gdouble sigma_alpha, gint n)
{
  NcHICosmoPriorTop *tp = g_slice_new (NcHICosmoPriorTop);
  nc_hicosmo_prior_top_set (tp, z, alpha, sigma_alpha, n);
  return tp;
}

/**
 * nc_hicosmo_prior_top_free:
 * @tp: FIXME
 *
 * FIXME
 */
void
nc_hicosmo_prior_top_free (NcHICosmoPriorTop *tp)
{
  nc_distance_free (tp->dist);
  g_slice_free (NcHICosmoPriorTop, tp);
}

/**
 * nc_hicosmo_prior_top_set:
 * @tp: FIXME
 * @z: FIXME
 * @alpha: FIXME
 * @sigma_alpha: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_prior_top_set (NcHICosmoPriorTop *tp, gdouble z, gdouble alpha, gdouble sigma_alpha, gint n)
{
  tp->z = z;
  tp->mean = atan(tan(M_PI/n)/cos(alpha));
  tp->sigma = sin(alpha)*tan(M_PI/n) / (cos(alpha)*cos(alpha) + tan(M_PI/n)*tan(M_PI/n));
  tp->sigma = sqrt(tp->sigma*tp->sigma*sigma_alpha*sigma_alpha);
  if (n == 2)
    tp->sigma = sigma_alpha / alpha;
}

/**
 * nc_hicosmo_prior_top_add:
 * @lh: a #NcmLikelihood
 * @z: FIXME
 * @alpha: FIXME
 * @sigma_alpha: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_prior_top_add (NcmLikelihood *lh, gdouble z, gdouble alpha, gdouble sigma_alpha, gint n)
{
  NcHICosmoPriorTop *tp = g_slice_new (NcHICosmoPriorTop);
  NcmMSetFunc *prior = ncm_mset_func_new (&_nc_hicosmo_prior_top_f, 0, 1, tp, (GDestroyNotify)nc_hicosmo_prior_top_free);
  nc_hicosmo_prior_top_set (tp, z, alpha, sigma_alpha, n);
  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}
