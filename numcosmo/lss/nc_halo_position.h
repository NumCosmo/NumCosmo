/***************************************************************************
 *            nc_halo_position.h
 *
 *  The Aug 08 08:19:02 2024
 *  Copyright  2024
 *  <code.caio@limadeoliveira.me>
 ****************************************************************************/
/*
 * nc_halo_position.h
 * Copyright (C) 2024 Caio Lima de Oliveira <code.caio@limadeoliveira.me>
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

#ifndef _NC_HALO_POSITION_H_
#define _NC_HALO_POSITION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_POSITION (nc_halo_position_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloPosition, nc_halo_position, NC, HALO_POSITION, NcmModel)

struct _NcHaloPositionClass
{
  /*< private >*/
  NcmModelClass parent_class;

  /* Padding to allow 18 virtual functions without breaking ABI */
  gpointer padding[18];
};

/**
 * NcHaloPositionSParams:
 * @NC_HALO_POSITION_RA: Right ascension
 * @NC_HALO_POSITION_DEC: Declination
 * @NC_HALO_POSITION_Z: Redshift
 *
 * Halo center model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_HALO_POSITION_SPARAMS >*/
{
  NC_HALO_POSITION_RA,
  NC_HALO_POSITION_DEC,
  NC_HALO_POSITION_Z,
  /* < private > */
  NC_HALO_POSITION_SPARAM_LEN, /*< skip >*/
} NcHaloPositionSParams;

NCM_MSET_MODEL_DECLARE_ID (nc_halo_position);

NcHaloPosition *nc_halo_position_new (NcDistance *dist);
NcHaloPosition *nc_halo_position_ref (NcHaloPosition *hp);

void nc_halo_position_free (NcHaloPosition *hp);
void nc_halo_position_clear (NcHaloPosition **hp);

#define NC_HALO_POSITION_DEFAULT_RA (0.0)
#define NC_HALO_POSITION_DEFAULT_DEC (0.0)
#define NC_HALO_POSITION_DEFAULT_Z (0.2)

#define NC_HALO_POSITION_DEFAULT_PARAMS_ABSTOL (0.0)

void nc_halo_position_polar_angles (NcHaloPosition *hp, gdouble ra, gdouble dec, gdouble *theta, gdouble *phi);
gdouble nc_halo_position_projected_radius (NcHaloPosition *hp, NcHICosmo *cosmo, gdouble theta);

G_END_DECLS

#endif /* _NC_HALO_POSITION_H_ */

