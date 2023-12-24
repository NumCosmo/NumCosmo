/***************************************************************************
 *            ncm_prior_flat_param.h
 *
 *  Wed August 03 10:55:33 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_flat_param.h
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

#ifndef _NCM_PRIOR_FLAT_PARAM_H_
#define _NCM_PRIOR_FLAT_PARAM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_prior_flat.h>

G_BEGIN_DECLS

#define NCM_TYPE_PRIOR_FLAT_PARAM (ncm_prior_flat_param_get_type ())

G_DECLARE_FINAL_TYPE (NcmPriorFlatParam, ncm_prior_flat_param, NCM, PRIOR_FLAT_PARAM, NcmPriorFlat)

NcmPriorFlatParam *ncm_prior_flat_param_new (NcmModelID mid, guint pid, gdouble x_low, gdouble x_upp, gdouble scale);
NcmPriorFlatParam *ncm_prior_flat_param_new_pindex (const NcmMSetPIndex *pi, gdouble x_low, gdouble x_upp, gdouble scale);
NcmPriorFlatParam *ncm_prior_flat_param_new_name (NcmMSet *mset, const gchar *name, gdouble x_low, gdouble x_upp, gdouble scale);

NcmPriorFlatParam *ncm_prior_flat_param_ref (NcmPriorFlatParam *pfp);

void ncm_prior_flat_param_free (NcmPriorFlatParam *pfp);
void ncm_prior_flat_param_clear (NcmPriorFlatParam **pfp);

G_END_DECLS

#endif /* _NCM_PRIOR_FLAT_PARAM_H_ */
