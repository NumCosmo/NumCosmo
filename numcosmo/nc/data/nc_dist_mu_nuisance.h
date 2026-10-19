/***************************************************************************
 *            nc_dist_mu_nuisance.h
 *
 *  Sun Oct 18 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_dist_mu_nuisance.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_DIST_MU_NUISANCE_H_
#define _NC_DIST_MU_NUISANCE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_DIST_MU_NUISANCE (nc_dist_mu_nuisance_get_type ())

G_DECLARE_FINAL_TYPE (NcDistMuNuisance, nc_dist_mu_nuisance, NC, DIST_MU_NUISANCE, NcmModel)

/**
 * NcDistMuNuisanceSParams:
 * @NC_DIST_MU_NUISANCE_SIGMA_V: peculiar velocity dispersion $\sigma_v$ [km/s]
 *
 * #NcDistMuNuisance scalar parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_DIST_MU_NUISANCE_SPARAMS,prefix=NC_DIST_MU_NUISANCE >*/
{
  NC_DIST_MU_NUISANCE_SIGMA_V,
  /* < private > */
  NC_DIST_MU_NUISANCE_SPARAM_LEN, /*< skip >*/
} NcDistMuNuisanceSParams;

/**
 * NcDistMuNuisanceVParams:
 * @NC_DIST_MU_NUISANCE_DMU: distance modulus offset of each sample $\delta\mu_s$ [mag]
 * @NC_DIST_MU_NUISANCE_SIGMA_V_S: peculiar velocity dispersion of each sample $\sigma_{v,s}$ [km/s]
 *
 * #NcDistMuNuisance vector parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_DIST_MU_NUISANCE_VPARAMS,prefix=NC_DIST_MU_NUISANCE >*/
{
  NC_DIST_MU_NUISANCE_DMU,
  NC_DIST_MU_NUISANCE_SIGMA_V_S,
  /* < private > */
  NC_DIST_MU_NUISANCE_VPARAM_LEN, /*< skip >*/
} NcDistMuNuisanceVParams;

NCM_MSET_MODEL_DECLARE_ID (nc_dist_mu_nuisance);

NcDistMuNuisance *nc_dist_mu_nuisance_new (guint nsamples);
NcDistMuNuisance *nc_dist_mu_nuisance_new_full (guint nsamples, gboolean sigma_v_per_sample);
NcDistMuNuisance *nc_dist_mu_nuisance_ref (NcDistMuNuisance *nuisance);

void nc_dist_mu_nuisance_free (NcDistMuNuisance *nuisance);
void nc_dist_mu_nuisance_clear (NcDistMuNuisance **nuisance);

guint nc_dist_mu_nuisance_get_nsamples (NcDistMuNuisance *nuisance);
guint nc_dist_mu_nuisance_get_nsigma_v (NcDistMuNuisance *nuisance);
gdouble nc_dist_mu_nuisance_get_sigma_v (NcDistMuNuisance *nuisance);
gdouble nc_dist_mu_nuisance_get_sigma_v_sample (NcDistMuNuisance *nuisance, const gint sample);
gdouble nc_dist_mu_nuisance_get_dmu (NcDistMuNuisance *nuisance, const guint sample);
gdouble nc_dist_mu_nuisance_sigma_mu_v (NcDistMuNuisance *nuisance, const gdouble z);
gdouble nc_dist_mu_nuisance_sigma_mu_v_sample (NcDistMuNuisance *nuisance, const gdouble z, const gint sample);

#define NC_DIST_MU_NUISANCE_DEFAULT_SIGMA_V (0.0)
#define NC_DIST_MU_NUISANCE_DEFAULT_DMU (0.0)
#define NC_DIST_MU_NUISANCE_DEFAULT_SIGMA_V_S (0.0)
#define NC_DIST_MU_NUISANCE_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_DIST_MU_NUISANCE_H_ */

