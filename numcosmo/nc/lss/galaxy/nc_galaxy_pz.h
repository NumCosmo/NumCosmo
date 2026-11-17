/***************************************************************************
 *            nc_galaxy_pz.h
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_pz.h
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

#ifndef _NC_GALAXY_PZ_H_
#define _NC_GALAXY_PZ_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_rng.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_PZ (nc_galaxy_pz_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyPz, nc_galaxy_pz, NC, GALAXY_PZ, GObject)

struct _NcGalaxyPzClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*get_range) (NcGalaxyPz *gpz, gdouble *z_min, gdouble *z_max);
  gdouble (*eval_p) (NcGalaxyPz *gpz, const gdouble z);
  gdouble (*eval_ln_p) (NcGalaxyPz *gpz, const gdouble z);
  gdouble (*get_norm) (NcGalaxyPz *gpz);
  gdouble (*gen) (NcGalaxyPz *gpz, NcmRNG *rng);

  /* Padding to allow adding up to 13 more virtual functions without breaking ABI. */
  gpointer padding[13];
};

NcGalaxyPz *nc_galaxy_pz_ref (NcGalaxyPz *gpz);
void nc_galaxy_pz_free (NcGalaxyPz *gpz);
void nc_galaxy_pz_clear (NcGalaxyPz **gpz);

void nc_galaxy_pz_get_range (NcGalaxyPz *gpz, gdouble *z_min, gdouble *z_max);
gdouble nc_galaxy_pz_eval_p (NcGalaxyPz *gpz, const gdouble z);
gdouble nc_galaxy_pz_eval_ln_p (NcGalaxyPz *gpz, const gdouble z);
gdouble nc_galaxy_pz_get_norm (NcGalaxyPz *gpz);
gdouble nc_galaxy_pz_gen (NcGalaxyPz *gpz, NcmRNG *rng);

G_END_DECLS

#endif /* _NC_GALAXY_PZ_H_ */

