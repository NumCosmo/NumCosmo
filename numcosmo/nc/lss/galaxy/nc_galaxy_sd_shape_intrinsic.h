/***************************************************************************
 *            nc_galaxy_sd_shape_intrinsic.h
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_intrinsic.h
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

#ifndef _NC_GALAXY_SD_SHAPE_INTRINSIC_H_
#define _NC_GALAXY_SD_SHAPE_INTRINSIC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/core/ncm_rng.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE_INTRINSIC (nc_galaxy_sd_shape_intrinsic_get_type ())
#define NC_TYPE_GALAXY_SD_SHAPE_INTRINSIC_DATA (nc_galaxy_sd_shape_intrinsic_data_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDShapeIntrinsic, nc_galaxy_sd_shape_intrinsic, NC, GALAXY_SD_SHAPE_INTRINSIC, NcmModel)
typedef struct _NcGalaxySDShapeIntrinsicData NcGalaxySDShapeIntrinsicData;

#define NC_GALAXY_SD_SHAPE_INTRINSIC_DATA(obj) ((NcGalaxySDShapeIntrinsicData *) (obj))

struct _NcGalaxySDShapeIntrinsicClass
{
  /*< private >*/
  NcmModelClass parent_class;

  void (*data_init) (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);
  void (*prepare) (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);
  gdouble (*eval_residual) (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x);
  void (*gen) (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
  gdouble (*e_rms) (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);

  /* Padding to allow adding up to 13 more virtual functions without breaking ABI. */
  gpointer padding[13];
};

/**
 * NcGalaxySDShapeIntrinsicData:
 * @e_rms: per-galaxy intrinsic shape scatter input (used by per-galaxy models;
 *   ignored by global models).
 * @jacobi_a: resolved Gauss-Jacobi weight exponent of $(1-x)$, $x = |\chi_I|^2$.
 * @jacobi_b: resolved Gauss-Jacobi weight exponent of $x$.
 *
 * Per-galaxy data for the intrinsic ellipticity distribution. Mirrors the
 * #NcGalaxySDShapeData / #NcGalaxySDPositionData idiom: typed public fields plus
 * an opaque @ldata holding the subclass-specific, varying/updatable resolved
 * parameters (filled by nc_galaxy_sd_shape_intrinsic_prepare()). The model is the
 * factory (nc_galaxy_sd_shape_intrinsic_data_new()).
 */
struct _NcGalaxySDShapeIntrinsicData
{
  gdouble e_rms;
  gdouble jacobi_a;
  gdouble jacobi_b;
  gpointer ldata;
  GDestroyNotify ldata_destroy;
  void (*ldata_read_row) (NcGalaxySDShapeIntrinsicData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_write_row) (NcGalaxySDShapeIntrinsicData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_required_columns) (NcGalaxySDShapeIntrinsicData *data, GList *columns);
  gatomicrefcount ref_count;
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_sd_shape_intrinsic);

GType nc_galaxy_sd_shape_intrinsic_data_get_type (void) G_GNUC_CONST;

NcGalaxySDShapeIntrinsic *nc_galaxy_sd_shape_intrinsic_ref (NcGalaxySDShapeIntrinsic *gsi);
void nc_galaxy_sd_shape_intrinsic_free (NcGalaxySDShapeIntrinsic *gsi);
void nc_galaxy_sd_shape_intrinsic_clear (NcGalaxySDShapeIntrinsic **gsi);

NcGalaxySDShapeIntrinsicData *nc_galaxy_sd_shape_intrinsic_data_new (NcGalaxySDShapeIntrinsic *gsi);
NcGalaxySDShapeIntrinsicData *nc_galaxy_sd_shape_intrinsic_data_ref (NcGalaxySDShapeIntrinsicData *data);
void nc_galaxy_sd_shape_intrinsic_data_unref (NcGalaxySDShapeIntrinsicData *data);

void nc_galaxy_sd_shape_intrinsic_data_read_row (NcGalaxySDShapeIntrinsicData *data, NcGalaxyWLObs *obs, const guint i);
void nc_galaxy_sd_shape_intrinsic_data_write_row (NcGalaxySDShapeIntrinsicData *data, NcGalaxyWLObs *obs, const guint i);
GList *nc_galaxy_sd_shape_intrinsic_data_required_columns (NcGalaxySDShapeIntrinsicData *data);

void nc_galaxy_sd_shape_intrinsic_prepare (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);
gdouble nc_galaxy_sd_shape_intrinsic_eval_residual (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x);
gdouble nc_galaxy_sd_shape_intrinsic_eval_p (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x);
void nc_galaxy_sd_shape_intrinsic_gen (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
gdouble nc_galaxy_sd_shape_intrinsic_e_rms (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_INTRINSIC_H_ */
