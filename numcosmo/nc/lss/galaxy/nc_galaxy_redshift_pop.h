/***************************************************************************
 *            nc_galaxy_redshift_pop.h
 *
 *  Wed Jul 31 21:09:32 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_pop.h
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#ifndef _NC_GALAXY_REDSHIFT_POP_H_
#define _NC_GALAXY_REDSHIFT_POP_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/core/ncm_rng.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_POP (nc_galaxy_redshift_pop_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyRedshiftPop, nc_galaxy_redshift_pop, NC, GALAXY_REDSHIFT_POP, NcmModel)

struct _NcGalaxyRedshiftPopClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*gen) (NcGalaxyRedshiftPop *gsdrp, NcmRNG *rng);
  gdouble (*eval) (NcGalaxyRedshiftPop *gsdrp, gdouble z);
  gdouble (*ln_eval) (NcGalaxyRedshiftPop *gsdrp, gdouble z);
  void (*set_lim) (NcGalaxyRedshiftPop *gsdrp, const gdouble z_min, const gdouble z_max);
  void (*get_lim) (NcGalaxyRedshiftPop *gsdrp, gdouble *z_min, gdouble *z_max);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[14];
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_redshift_pop);

NcGalaxyRedshiftPop *nc_galaxy_redshift_pop_ref (NcGalaxyRedshiftPop *gsdrp);

void nc_galaxy_redshift_pop_free (NcGalaxyRedshiftPop *gsdrp);
void nc_galaxy_redshift_pop_clear (NcGalaxyRedshiftPop **gsdrp);

void nc_galaxy_redshift_pop_set_lim (NcGalaxyRedshiftPop *gsdrp, const gdouble z_min, const gdouble z_max);
void nc_galaxy_redshift_pop_get_lim (NcGalaxyRedshiftPop *gsdrp, gdouble *z_min, gdouble *z_max);

gdouble nc_galaxy_redshift_pop_gen (NcGalaxyRedshiftPop *gsdrp, NcmRNG *rng);
gdouble nc_galaxy_redshift_pop_eval (NcGalaxyRedshiftPop *gsdrp, gdouble z);
gdouble nc_galaxy_redshift_pop_ln_eval (NcGalaxyRedshiftPop *gsdrp, gdouble z);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_POP_H_ */

