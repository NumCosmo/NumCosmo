/***************************************************************************
 *            ncm_priors.h
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

#ifndef _NCM_PRIORS_H_
#define _NCM_PRIORS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/math/ncm_likelihood.h>

G_BEGIN_DECLS

typedef struct _NcmPriorGauss NcmPriorGauss;

/**
 * NcmPriorGauss:
 *
 * FIXME
 */
struct _NcmPriorGauss
{
  /*< private >*/
  NcmMSetFunc *func;
  NcmMSetPIndex pi;
  gdouble z;
  gdouble mean;
  gdouble sigma;
};

void ncm_prior_add_oneside_a_inf_param (NcmLikelihood *lh, NcmModelID gmid, guint pid, gdouble a, gdouble s);
void ncm_prior_add_oneside_a_inf_const_func (NcmLikelihood *lh, NcmMSetFunc *func, gdouble mean, gdouble s);
void ncm_prior_add_oneside_a_inf_func (NcmLikelihood *lh, NcmMSetFunc *func, gdouble z, gdouble mean, gdouble s);
void ncm_prior_add_oneside_a_inf (NcmLikelihood *lh, NcmPriorGauss *gp);
void ncm_prior_add_twoside_a_b (NcmLikelihood *lh, NcmModelID gmid, guint pid, gdouble a, gdouble b, gdouble s);
void ncm_prior_add_positive (NcmLikelihood *lh, NcmModelID gmid, guint pid);
void ncm_prior_add_gaussian (NcmLikelihood *lh, NcmPriorGauss *gp);
void ncm_prior_add_gaussian_data (NcmLikelihood *lh, NcmModelID gmid, guint pid, gdouble mean, gdouble sigma);
void ncm_prior_add_gaussian_func (NcmLikelihood *lh, NcmMSetFunc *func, gdouble z, gdouble mean, gdouble sigma);
void ncm_prior_add_gaussian_const_func (NcmLikelihood *lh, NcmMSetFunc *func, gdouble mean, gdouble sigma);

G_END_DECLS

#endif /* _NCM_PRIORS_H_ */
