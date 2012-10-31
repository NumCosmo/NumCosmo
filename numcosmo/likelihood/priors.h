/***************************************************************************
 *            priors.h
 *
 *  Wed Mar 19 12:46:57 2008
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

#ifndef _NC_PRIORS_H
#define _NC_PRIORS_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/likelihood/likelihood.h>

G_BEGIN_DECLS

typedef struct _NcGaussianPrior NcGaussianPrior;

/**
 * NcGaussianPrior:
 *
 * FIXME
 */
struct _NcGaussianPrior
{
  /*< private >*/
  NcmMSetFunc *func;
  NcmMSetPIndex pi;
  gdouble z;
  gdouble mean;
  gdouble sigma;
};

typedef struct _NcTopologicalPrior NcTopologicalPrior;

/**
 * NcTopologicalPrior:
 *
 * FIXME
 */
struct _NcTopologicalPrior
{
  /*< private >*/
  NcDistance *dist;
  gdouble z;
  gdouble mean;
  gdouble sigma;
};

NcTopologicalPrior *nc_prior_topological_new (gdouble z, gdouble alpha, gdouble sigma_alpha, gint n);
void nc_prior_topological_free (NcTopologicalPrior *tp);
gboolean nc_prior_topological_set (NcTopologicalPrior *tp, gdouble z, gdouble alpha, gdouble sigma_alpha, gint n);
gboolean nc_prior_add_topological (NcLikelihood *lh, gdouble z, gdouble alpha, gdouble sigma_alpha, gint n);

gboolean nc_prior_add_oneside_a_inf_param (NcLikelihood *lh, NcmModelID gmid, guint pid, gdouble a, gdouble s);
gboolean nc_prior_add_oneside_a_inf_const_func (NcLikelihood *lh, NcmMSetFunc *func, gdouble mean, gdouble s);
gboolean nc_prior_add_oneside_a_inf_func (NcLikelihood *lh, NcmMSetFunc *func, gdouble z, gdouble mean, gdouble s);
gboolean nc_prior_add_oneside_a_inf (NcLikelihood *lh, NcGaussianPrior *gp);
gboolean nc_prior_add_twoside_a_b (NcLikelihood *lh, NcmModelID gmid, guint pid, gdouble a, gdouble b, gdouble s);
gboolean nc_prior_add_positive (NcLikelihood *lh, NcmModelID gmid, guint pid);
gboolean nc_prior_add_gaussian (NcLikelihood *lh, NcGaussianPrior *gp);
gboolean nc_prior_add_gaussian_data (NcLikelihood *lh, NcmModelID gmid, guint pid, gdouble mean, gdouble sigma);
gboolean nc_prior_add_gaussian_func (NcLikelihood *lh, NcmMSetFunc *func, gdouble z, gdouble mean, gdouble sigma);
gboolean nc_prior_add_gaussian_const_func (NcLikelihood *lh, NcmMSetFunc *func, gdouble mean, gdouble sigma);

G_END_DECLS

#endif /* _PRIORS_H */
