/***************************************************************************
 *            nc_galaxy_redshift_obs.h
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_obs.h
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

#ifndef _NC_GALAXY_REDSHIFT_OBS_H_
#define _NC_GALAXY_REDSHIFT_OBS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/core/ncm_rng.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_OBS (nc_galaxy_redshift_obs_get_type ())
#define NC_TYPE_GALAXY_REDSHIFT_OBS_DATA (nc_galaxy_redshift_obs_data_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyRedshiftObs, nc_galaxy_redshift_obs, NC, GALAXY_REDSHIFT_OBS, NcmModel)
typedef struct _NcGalaxyRedshiftObsData NcGalaxyRedshiftObsData;

#define NC_GALAXY_REDSHIFT_OBS_DATA(obj) ((NcGalaxyRedshiftObsData *) (obj))

struct _NcGalaxyRedshiftObsClass
{
  /*< private >*/
  NcmModelClass parent_class;

  void (*data_init) (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data);
  gdouble (*eval) (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z);
  gdouble (*gen) (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, NcmRNG *rng);
  gdouble (*window_mass) (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, const gdouble obs_lo, const gdouble obs_hi);
  void (*get_true_z_lim) (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, gdouble *z_min, gdouble *z_max);

  /* Padding to allow adding up to 13 more virtual functions without breaking ABI. */
  gpointer padding[13];
};

/**
 * NcGalaxyRedshiftObsData:
 *
 * Per-galaxy data for the photometric-redshift observable model. Following the
 * #NcGalaxyShapePopData idiom, the whole per-galaxy photometric-redshift
 * observation (e.g. the point estimate and its scatter) lives in the opaque
 * @ldata owned by the concrete subclass: the observation's structure is defined
 * by the observable model, so the redshift calculator never accesses it directly.
 */
struct _NcGalaxyRedshiftObsData
{
  gpointer ldata;
  GDestroyNotify ldata_destroy;
  void (*ldata_read_row) (NcGalaxyRedshiftObsData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_write_row) (NcGalaxyRedshiftObsData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_required_columns) (NcGalaxyRedshiftObsData *data, GList **columns);
  gatomicrefcount ref_count;
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_redshift_obs);

GType nc_galaxy_redshift_obs_data_get_type (void) G_GNUC_CONST;

NcGalaxyRedshiftObs *nc_galaxy_redshift_obs_ref (NcGalaxyRedshiftObs *gsdre);
void nc_galaxy_redshift_obs_free (NcGalaxyRedshiftObs *gsdre);
void nc_galaxy_redshift_obs_clear (NcGalaxyRedshiftObs **gsdre);

NcGalaxyRedshiftObsData *nc_galaxy_redshift_obs_data_new (NcGalaxyRedshiftObs *gsdre);
NcGalaxyRedshiftObsData *nc_galaxy_redshift_obs_data_ref (NcGalaxyRedshiftObsData *data);
void nc_galaxy_redshift_obs_data_unref (NcGalaxyRedshiftObsData *data);

void nc_galaxy_redshift_obs_data_read_row (NcGalaxyRedshiftObsData *data, NcGalaxyWLObs *obs, const guint i);
void nc_galaxy_redshift_obs_data_write_row (NcGalaxyRedshiftObsData *data, NcGalaxyWLObs *obs, const guint i);
GList *nc_galaxy_redshift_obs_data_required_columns (NcGalaxyRedshiftObsData *data);

gdouble nc_galaxy_redshift_obs_eval (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z);
gdouble nc_galaxy_redshift_obs_gen (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, NcmRNG *rng);
gdouble nc_galaxy_redshift_obs_window_mass (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, const gdouble obs_lo, const gdouble obs_hi);
void nc_galaxy_redshift_obs_get_true_z_lim (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, gdouble *z_min, gdouble *z_max);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_OBS_H_ */

