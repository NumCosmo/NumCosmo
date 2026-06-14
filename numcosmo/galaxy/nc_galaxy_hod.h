/***************************************************************************
 *            nc_galaxy_hod.h
 *
 *  Sun Jun 14 12:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_hod.h
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

#ifndef _NC_GALAXY_HOD_H
#define _NC_GALAXY_HOD_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_HOD (nc_galaxy_hod_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyHOD, nc_galaxy_hod, NC, GALAXY_HOD, NcmModel)

struct _NcGalaxyHODClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*mean_n_central) (NcGalaxyHOD *hod, const gdouble lnM);
  gdouble (*mean_n_satellite) (NcGalaxyHOD *hod, const gdouble lnM);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[16];
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_hod);

NcGalaxyHOD *nc_galaxy_hod_ref (NcGalaxyHOD *hod);

void nc_galaxy_hod_free (NcGalaxyHOD *hod);
void nc_galaxy_hod_clear (NcGalaxyHOD **hod);

void nc_galaxy_hod_set_stochastic_central (NcGalaxyHOD *hod, gboolean stochastic_central);
gboolean nc_galaxy_hod_get_stochastic_central (NcGalaxyHOD *hod);

gdouble nc_galaxy_hod_mean_n_central (NcGalaxyHOD *hod, const gdouble lnM);
gdouble nc_galaxy_hod_mean_n_satellite (NcGalaxyHOD *hod, const gdouble lnM);

void nc_galaxy_hod_gen (NcGalaxyHOD *hod, const gdouble lnM, NcmRNG *rng, gint *n_central, gint *n_satellite);

G_END_DECLS

#endif /* _NC_GALAXY_HOD_H */
