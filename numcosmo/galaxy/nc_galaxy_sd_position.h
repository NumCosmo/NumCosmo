/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_position.h
 *
 *  Sat May 20 18:41:47 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position.h
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
#ifndef _NC_GALAXY_SD_POSITION_H_
#define _NC_GALAXY_SD_POSITION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/galaxy/nc_galaxy_sd_obs_redshift.h>
#include <numcosmo/galaxy/nc_galaxy_wl_obs.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_POSITION (nc_galaxy_sd_position_get_type ())
#define NC_TYPE_GALAXY_SD_POSITION_DATA (nc_galaxy_sd_position_data_get_type ())
#define NC_TYPE_GALAXY_SD_POSITION_INTEGRAND (nc_galaxy_sd_position_integrand_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDPosition, nc_galaxy_sd_position, NC, GALAXY_SD_POSITION, NcmModel)

typedef struct _NcGalaxySDPositionData NcGalaxySDPositionData;

NCM_UTIL_DECLARE_CALLBACK (NcGalaxySDPositionIntegrand,
                           NC_GALAXY_SD_POSITION_INTEGRAND,
                           nc_galaxy_sd_position_integrand,
                           gdouble,
                           NCM_UTIL_CALLBACK_ARGS (NcGalaxySDPositionData * data))


#define NC_GALAXY_SD_POSITION_DATA(obj) ((NcGalaxySDPositionData *) (obj))

struct _NcGalaxySDPositionClass
{
  /*< private >*/
  NcmModelClass parent_class;

  void (*gen) (NcGalaxySDPosition *gsdp, NcGalaxySDPositionData *data, NcmRNG *rng);
  NcGalaxySDPositionIntegrand *(*integ) (NcGalaxySDPosition *gsdp);
  gboolean (*set_ra_lim) (NcGalaxySDPosition *gsdp, const gdouble ra_min, const gdouble ra_max);
  gboolean (*set_dec_lim) (NcGalaxySDPosition *gsdp, const gdouble dec_min, const gdouble dec_max);
  gboolean (*get_ra_lim) (NcGalaxySDPosition *gsdp, gdouble *ra_min, gdouble *ra_max);
  gboolean (*get_dec_lim) (NcGalaxySDPosition *gsdp, gdouble *dec_min, gdouble *dec_max);
  NcGalaxySDPositionData *(*data_new) (NcGalaxySDPosition *gsdp, NcGalaxySDObsRedshiftData *sdz_data);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[10];
};

struct _NcGalaxySDPositionData
{
  NcGalaxySDObsRedshiftData *sdz_data;
  gdouble ra;
  gdouble dec;
  gpointer ldata;
  GDestroyNotify ldata_destroy;
  gpointer (*ldata_copy) (gpointer ldata);
  void (*ldata_read_row) (NcGalaxySDPositionData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_write_row) (NcGalaxySDPositionData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_required_columns) (NcGalaxySDPositionData *data, GList *columns);
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_sd_position);

GType nc_galaxy_sd_position_data_get_type (void) G_GNUC_CONST;

NcGalaxySDPositionData *nc_galaxy_sd_position_data_copy (NcGalaxySDPositionData *data);
void nc_galaxy_sd_position_data_free (NcGalaxySDPositionData *data);
void nc_galaxy_sd_position_data_read_row (NcGalaxySDPositionData *data, NcGalaxyWLObs *obs, const guint i);
void nc_galaxy_sd_position_data_write_row (NcGalaxySDPositionData *data, NcGalaxyWLObs *obs, const guint i);
GList *nc_galaxy_sd_position_data_required_columns (NcGalaxySDPositionData *data);

NcGalaxySDPosition *nc_galaxy_sd_position_ref (NcGalaxySDPosition *gsdp);

void nc_galaxy_sd_position_free (NcGalaxySDPosition *gsdp);
void nc_galaxy_sd_position_clear (NcGalaxySDPosition **gsdp);

gboolean nc_galaxy_sd_position_set_ra_lim (NcGalaxySDPosition *gsdp, const gdouble ra_min, const gdouble ra_max);
gboolean nc_galaxy_sd_position_get_ra_lim (NcGalaxySDPosition *gsdp, gdouble *ra_min, gdouble *ra_max);
gboolean nc_galaxy_sd_position_set_dec_lim (NcGalaxySDPosition *gsdp, const gdouble dec_min, const gdouble dec_max);
gboolean nc_galaxy_sd_position_get_dec_lim (NcGalaxySDPosition *gsdp, gdouble *dec_min, gdouble *dec_max);
NcGalaxySDPositionIntegrand *nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp);

NcGalaxySDPositionData *nc_galaxy_sd_position_data_new (NcGalaxySDPosition *gsdp, NcGalaxySDObsRedshiftData *sdz_data);

#define NC_GALAXY_SD_POSITION_COL_RA "ra"
#define NC_GALAXY_SD_POSITION_COL_DEC "dec"

G_END_DECLS

#endif /* _NC_GALAXY_SD_POSITION_H_ */

