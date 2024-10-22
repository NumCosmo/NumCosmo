/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_obs_redshift.h
 *
 *  Thu Aug 1 00:48:12 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift.h
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
#ifndef _NC_GALAXY_SD_OBS_REDSHIFT_H_
#define _NC_GALAXY_SD_OBS_REDSHIFT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_get_type ())
#define NC_TYPE_GALAXY_SD_OBS_REDSHIFT_DATA (nc_galaxy_sd_obs_redshift_data_get_type ())
#define NC_TYPE_GALAXY_SD_OBS_REDSHIFT_INTEGRAND (nc_galaxy_sd_obs_redshift_integrand_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDObsRedshift, nc_galaxy_sd_obs_redshift, NC, GALAXY_SD_OBS_REDSHIFT, NcmModel)
typedef struct _NcGalaxySDObsRedshiftData NcGalaxySDObsRedshiftData;

NCM_UTIL_DECLARE_CALLBACK (NcGalaxySDObsRedshiftIntegrand,
                           NC_GALAXY_SD_OBS_REDSHIFT_INTEGRAND,
                           nc_galaxy_sd_obs_redshift_integrand,
                           gdouble,
                           NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxySDObsRedshiftData * data))

#define NC_GALAXY_SD_OBS_REDSHIFT_DATA(obj) ((NcGalaxySDObsRedshiftData *) (obj))

struct _NcGalaxySDObsRedshiftClass
{
  /*< private >*/
  NcmModelClass parent_class;

  void (*gen) (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng);
  NcGalaxySDObsRedshiftIntegrand *(*integ) (NcGalaxySDObsRedshift *gsdor);
  void (*data_init) (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[15];
};

struct _NcGalaxySDObsRedshiftData
{
  gdouble z;
  gpointer ldata;
  GDestroyNotify ldata_destroy;
  void (*ldata_read_row) (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_write_row) (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_required_columns) (NcGalaxySDObsRedshiftData *data, GList *columns);
  gatomicrefcount ref_count;
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_sd_obs_redshift);

GType nc_galaxy_sd_obs_redshift_data_get_type (void) G_GNUC_CONST;

NcGalaxySDObsRedshiftData *nc_galaxy_sd_obs_redshift_data_ref (NcGalaxySDObsRedshiftData *data);
void nc_galaxy_sd_obs_redshift_data_unref (NcGalaxySDObsRedshiftData *data);
void nc_galaxy_sd_obs_redshift_data_read_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i);
void nc_galaxy_sd_obs_redshift_data_write_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i);
GList *nc_galaxy_sd_obs_redshift_data_required_columns (NcGalaxySDObsRedshiftData *data);

NcGalaxySDObsRedshift *nc_galaxy_sd_obs_redshift_ref (NcGalaxySDObsRedshift *gsdor);

void nc_galaxy_sd_obs_redshift_free (NcGalaxySDObsRedshift *gsdor);
void nc_galaxy_sd_obs_redshift_clear (NcGalaxySDObsRedshift **gsdor);

NcGalaxySDObsRedshiftIntegrand *nc_galaxy_sd_obs_redshift_integ (NcGalaxySDObsRedshift *gsdor);
NcGalaxySDObsRedshiftData *nc_galaxy_sd_obs_redshift_data_new (NcGalaxySDObsRedshift *gsdor);

#define NC_GALAXY_SD_OBS_REDSHIFT_COL_Z "z"

G_END_DECLS

#endif /* _NC_GALAXY_SD_OBS_REDSHIFT_H_ */

