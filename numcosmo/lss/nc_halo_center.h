/***************************************************************************
 *            nc_halo_center.h
 *
 *  The Aug 08 08:19:02 2024
 *  Copyright  2024
 *  <code.caio@limadeoliveira.me>
 ****************************************************************************/
/*
 * nc_halo_center.h
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

#ifndef _NC_HALO_CENTER_H_
#define _NC_HALO_CENTER_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_CENTER (nc_halo_center_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloCenter, nc_halo_center, NC, HALO_CENTER, NcmModel)

struct _NcHaloCenterClass
{
  /*< private >*/
  NcmModelClass parent_class;

  /* Padding to allow 18 virtual functions without breaking ABI */
  gpointer padding[18];
};

/**
 * NcHaloCenterSParams:
 * @NC_HALO_CENTER_RA: Right ascension
 * @NC_HALO_CENTER_DEC: Declination
 * @NC_HALO_CENTER_Z: Redshift
 *
 * Halo center model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_HALO_CENTER_SPARAMS >*/
{
  NC_HALO_CENTER_RA,
  NC_HALO_CENTER_DEC,
  NC_HALO_CENTER_Z,
  /* < private > */
  NC_HALO_CENTER_SPARAM_LEN, /*< skip >*/
} NcHaloCenterSParams;

NCM_MSET_MODEL_DECLARE_ID (nc_halo_center);

NcHaloCenter *nc_halo_center_new (void);
NcHaloCenter *nc_halo_center_ref (NcHaloCenter *hc);

void nc_halo_center_free (NcHaloCenter *hc);
void nc_halo_center_clear (NcHaloCenter *hc);

#define NC_HALO_CENTER_DEFAULT_RA (0.0)
#define NC_HALO_CENTER_DEFAULT_DEC (0.0)
#define NC_HALO_CENTER_DEFAULT_Z (0.2)

#define NC_HALO_CENTER_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HALO_CENTER_H_ */

