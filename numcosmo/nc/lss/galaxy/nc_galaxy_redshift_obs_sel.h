/***************************************************************************
 *            nc_galaxy_redshift_obs_sel.h
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_obs_sel.h
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_GALAXY_REDSHIFT_OBS_SEL_H_
#define _NC_GALAXY_REDSHIFT_OBS_SEL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_OBS_SEL (nc_galaxy_redshift_obs_sel_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyRedshiftObsSel, nc_galaxy_redshift_obs_sel, NC, GALAXY_REDSHIFT_OBS_SEL, NcmModel)

/**
 * NcGalaxyRedshiftObsSelClass:
 *
 * The population-level distribution of a photo-z observable across the galaxy
 * population at a given true redshift. Distinct from the per-galaxy conditional
 * #NcGalaxyRedshiftObs: its scatter is a population-level model
 * parameter (not per-galaxy data), and it is conceptually a mixture over the
 * per-galaxy scatter, free to diverge from the single-galaxy kernel. Used by the
 * binning calculator to build P(z | I in W).
 */
struct _NcGalaxyRedshiftObsSelClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*eval) (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs);
  gdouble (*window_mass) (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs_lo, const gdouble obs_hi);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[16];
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_redshift_obs_sel);

NcGalaxyRedshiftObsSel *nc_galaxy_redshift_obs_sel_ref (NcGalaxyRedshiftObsSel *gsdrop);
void nc_galaxy_redshift_obs_sel_free (NcGalaxyRedshiftObsSel *gsdrop);
void nc_galaxy_redshift_obs_sel_clear (NcGalaxyRedshiftObsSel **gsdrop);

gdouble nc_galaxy_redshift_obs_sel_eval (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs);
gdouble nc_galaxy_redshift_obs_sel_window_mass (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs_lo, const gdouble obs_hi);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_OBS_SEL_H_ */
