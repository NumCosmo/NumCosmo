/***************************************************************************
 *            ncm_prior_flat.h
 *
 *  Wed August 03 16:58:26 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_flat.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_PRIOR_FLAT_H_
#define _NCM_PRIOR_FLAT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_prior.h>

G_BEGIN_DECLS

#define NCM_TYPE_PRIOR_FLAT (ncm_prior_flat_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmPriorFlat, ncm_prior_flat, NCM, PRIOR_FLAT, NcmPrior)

typedef gdouble (*NcmPriorFlatMean) (NcmPriorFlat *pf, NcmMSet *mset);

struct _NcmPriorFlatClass
{
  /*< private >*/
  NcmPriorClass parent_class;
  NcmPriorFlatMean mean;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[17];
};

NcmPriorFlat *ncm_prior_flat_ref (NcmPriorFlat *pf);
void ncm_prior_flat_free (NcmPriorFlat *pf);
void ncm_prior_flat_clear (NcmPriorFlat **pf);

void ncm_prior_flat_set_x_low (NcmPriorFlat *pf, const gdouble x_low);
void ncm_prior_flat_set_x_upp (NcmPriorFlat *pf, const gdouble x_upp);
void ncm_prior_flat_set_scale (NcmPriorFlat *pf, const gdouble scale);
void ncm_prior_flat_set_var (NcmPriorFlat *pf, const gdouble var);
void ncm_prior_flat_set_h0 (NcmPriorFlat *pf, const gdouble h0);

gdouble ncm_prior_flat_get_x_low (NcmPriorFlat *pf);
gdouble ncm_prior_flat_get_x_upp (NcmPriorFlat *pf);
gdouble ncm_prior_flat_get_scale (NcmPriorFlat *pf);
gdouble ncm_prior_flat_get_var (NcmPriorFlat *pf);
gdouble ncm_prior_flat_get_h0 (NcmPriorFlat *pf);

G_END_DECLS

#endif /* _NCM_PRIOR_FLAT_H_ */

