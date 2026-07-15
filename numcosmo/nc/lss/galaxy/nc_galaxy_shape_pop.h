/***************************************************************************
 *            nc_galaxy_shape_pop.h
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop.h
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

#ifndef _NC_GALAXY_SHAPE_POP_H_
#define _NC_GALAXY_SHAPE_POP_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/core/ncm_rng.h>
#include <numcosmo/ncm/algebra/ncm_laurent_series.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_get_type ())
#define NC_TYPE_GALAXY_SHAPE_POP_DATA (nc_galaxy_shape_pop_data_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyShapePop, nc_galaxy_shape_pop, NC, GALAXY_SHAPE_POP, NcmModel)
typedef struct _NcGalaxyShapePopData NcGalaxyShapePopData;

#define NC_GALAXY_SHAPE_POP_DATA(obj) ((NcGalaxyShapePopData *) (obj))

struct _NcGalaxyShapePopClass
{
  /*< private >*/
  NcmModelClass parent_class;

  void (*data_init) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
  void (*prepare) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
  gdouble (*eval_p) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble x);
  gdouble (*eval_p_rho2) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble rho2);
  void (*gen) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
  gdouble (*e_rms) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);

#ifndef NUMCOSMO_GIR_SCAN

  /*
   * Taylor-in-g analog of eval_p_rho2(): given rho2(g) = |chi_I(chi_L,g)|^2's
   * own g-Taylor coefficients (population-independent shear-map output,
   * @rho2_series), returns this population's normalized-density composition
   * P(rho2(g))'s g-Taylor coefficients. Unlike eval_p_rho2, there is no
   * sensible generic default (the composition depends entirely on the
   * population's own functional form), so the base class default just
   * errors clearly, matching NcGalaxySDShape's direct_estimate idiom for
   * the same "not every subclass needs this" situation. Consumed by
   * #NcGalaxyShapeFactorSeriesLensed. Native NcmLaurentSeries types only
   * (not introspectable), like ncm_laurent_series.h's own native interface
   * -- hence guarded out of the GIR scan.
   */
  void (*eval_p_rho2_g_series) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data,
                                NcmLaurentSeriesArena *arena, NcmLaurentSeries * const *rho2_series,
                                guint order, NcmLaurentSeries **out);
#endif /* NUMCOSMO_GIR_SCAN */

  /* Padding to allow adding up to 11 more virtual functions without breaking ABI. */
  gpointer padding[11];
};

/**
 * NcGalaxyShapePopData:
 * @e_rms: per-galaxy intrinsic shape scatter input (used by per-galaxy models;
 *   ignored by global models).
 *
 * Per-galaxy data for the intrinsic ellipticity distribution. Mirrors the
 * #NcGalaxySDShapeData / #NcGalaxySDPositionData idiom: typed public fields plus
 * an opaque @ldata holding the subclass-specific, varying/updatable resolved
 * parameters (filled by nc_galaxy_shape_pop_prepare()). The model is the
 * factory (nc_galaxy_shape_pop_data_new()).
 */
struct _NcGalaxyShapePopData
{
  gdouble e_rms;
  gpointer ldata;
  GDestroyNotify ldata_destroy;
  void (*ldata_read_row) (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_write_row) (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_required_columns) (NcGalaxyShapePopData *data, GList **columns);

  /*
   * Optional capability: NULL unless the concrete model parameterizes its
   * density through an (untruncated) Gaussian width sigma (e.g. the Gauss
   * family, Global or per-galaxy); NULL for models where that concept does
   * not apply (e.g. Beta). Consumers that need it (approximations linearizing
   * around a Gaussian, like the variance-add marginalization) go through
   * nc_galaxy_shape_pop_get_sigma(), which errors clearly if unsupported.
   */
  gdouble (*ldata_get_sigma) (NcGalaxyShapePopData *data);

  /*
   * Optional capability: NULL means the model is assumed radially symmetric
   * about chi_I=0 (its own peak, e.g. the Gauss family), i.e. mode_x=0. A
   * concrete model peaked away from the origin (e.g. NcGalaxyShapePopBeta
   * with mu far from 0) overrides this with the mode of x=|chi_I|^2. Unlike
   * ldata_get_sigma, 0 is always a meaningful default here, so consumers
   * (nc_galaxy_shape_pop_get_mode_x()) never need to error on NULL.
   */
  gdouble (*ldata_get_mode_x) (NcGalaxyShapePopData *data);

  gatomicrefcount ref_count;
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_shape_pop);

GType nc_galaxy_shape_pop_data_get_type (void) G_GNUC_CONST;

NcGalaxyShapePop *nc_galaxy_shape_pop_ref (NcGalaxyShapePop *gsp);
void nc_galaxy_shape_pop_free (NcGalaxyShapePop *gsp);
void nc_galaxy_shape_pop_clear (NcGalaxyShapePop **gsp);

NcGalaxyShapePopData *nc_galaxy_shape_pop_data_new (NcGalaxyShapePop *gsp);
NcGalaxyShapePopData *nc_galaxy_shape_pop_data_ref (NcGalaxyShapePopData *data);
void nc_galaxy_shape_pop_data_unref (NcGalaxyShapePopData *data);

void nc_galaxy_shape_pop_data_read_row (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i);
void nc_galaxy_shape_pop_data_write_row (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i);
GList *nc_galaxy_shape_pop_data_required_columns (NcGalaxyShapePopData *data);

void nc_galaxy_shape_pop_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
gdouble nc_galaxy_shape_pop_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble x);
gdouble nc_galaxy_shape_pop_eval_p_rho2 (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble rho2);
void nc_galaxy_shape_pop_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
gdouble nc_galaxy_shape_pop_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
gdouble nc_galaxy_shape_pop_get_sigma (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
gdouble nc_galaxy_shape_pop_get_mode_x (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);

#ifndef NUMCOSMO_GIR_SCAN
void nc_galaxy_shape_pop_eval_p_rho2_g_series (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data,
                                               NcmLaurentSeriesArena *arena, NcmLaurentSeries * const *rho2_series,
                                               guint order, NcmLaurentSeries **out);

#endif /* NUMCOSMO_GIR_SCAN */

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_POP_H_ */

