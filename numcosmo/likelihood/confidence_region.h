/***************************************************************************
 *            confidence-region.h
 *
 *  Fri Aug 15 15:22:57 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NC_CONFIDENCE_REGION_H
#define _NC_CONFIDENCE_REGION_H

#include <glib.h>

G_BEGIN_DECLS

/**
 * NcConfidenceRegionSearchType:
 * @NC_CONFIDENCE_REGION_SEARCH_1D: FIXME
 * @NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_RADIUS: FIXME
 * @NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_ANGLE: FIXME
 * @NC_CONFIDENCE_REGION_SEARCH_CONSTRAIN_1D: FIXME
 *
 * FIXME
 */
typedef enum _NcConfidenceRegionSearchType
{
  NC_CONFIDENCE_REGION_SEARCH_1D = 0,
  NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_RADIUS,
  NC_CONFIDENCE_REGION_SEARCH_ELLIPTIC_ANGLE,
  NC_CONFIDENCE_REGION_SEARCH_CONSTRAIN_1D
} NcConfidenceRegionSearchType;

typedef struct _NcConfidenceRegion NcConfidenceRegion;

#define NCM_FIT_CR_MAX_DIM 2

/**
 * NcConfidenceRegion:
 *
 * FIXME
 */
struct _NcConfidenceRegion
{
  /*< private >*/
  NcConfidenceRegionSearchType search_type;
  NcmFit *bestfit;
  NcmFit *constrained;
  struct _NcGaussianPrior *prior;
  NcmMSetPIndex pi[NCM_FIT_CR_MAX_DIM];
  gdouble shift[NCM_FIT_CR_MAX_DIM];
  NcmVector *covar_ev;
  NcmMatrix *covar_orto;
  gulong total_func_eval;
  guint n;
  gdouble theta;
  gdouble r;
  gdouble chi2;
  gdouble chi2_min;
  gboolean minimize;
  gboolean inv_order;
  GList *points;
};

typedef struct _NcConfidenceRegion2dPoint NcConfidenceRegion2dPoint;

/**
 * NcConfidenceRegion2dPoint:
 *
 * FIXME
 */
struct _NcConfidenceRegion2dPoint
{
  /*< private >*/
  gdouble x;
  gdouble y;
  gdouble theta;
  gdouble p1;
  gdouble p2;
};

#define NC_CR_X(cr,r,theta) ((cr)->shift[0] + (r) * cos(theta))
#define NC_CR_Y(cr,r,theta) ((cr)->shift[1] + (r) * sin(theta))
#define NC_CR_PVAL(cr,n) (ncm_mset_param_get((cr)->fit->mset,NC_CR_PNUM((cr),(n))))
#define NC_CR_COVAR_EV(cr,n) (*ncm_vector_ptr((cr)->covar_ev,(n)))

NcConfidenceRegion *nc_confidence_region_new_1d (NcmFit *fit, NcmModelID gmid, guint pid);
NcConfidenceRegion *nc_confidence_region_new (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2);
gboolean nc_confidence_region_free (NcConfidenceRegion *cr);

NcmMSetPIndex *ncm_fit_cr_get_pi (NcConfidenceRegion *cr, guint n);
NcmMSetPIndex *ncm_fit_cr_get_pi_array (NcConfidenceRegion *cr);

gdouble ncm_fit_cr_get_param (NcConfidenceRegion *cr, guint n);
gdouble ncm_fit_cr_get_bf_param (NcConfidenceRegion *cr, guint n);

gdouble ncm_fit_cr_get_param_shift (NcConfidenceRegion *cr, guint n);
void ncm_fit_cr_set_param_shift (NcConfidenceRegion *cr, guint n, gdouble s);

gdouble ncm_fit_cr_get_var (NcConfidenceRegion *cr, guint n);
gdouble ncm_fit_cr_get_sd (NcConfidenceRegion *cr, guint n);

G_END_DECLS

#endif /* CONFIDENCE_REGION_H */
