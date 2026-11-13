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
  gdouble (*eval_p) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r);
  void (*gen) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
  gdouble (*e_rms) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
  gdouble (*exponent_at_origin) (NcGalaxyShapePop *gsp);

  /*
   * Batched form of eval_p(): evaluates P_pop(r) at every element of @r in
   * one call, writing into @p (allocating it first if it points to NULL --
   * see nc_galaxy_shape_pop_eval_p_array()'s own doc comment). The generic
   * default just loops calling eval_p(), so every subclass is safely
   * callable through this vfunc; a subclass whose functional form batches
   * well (e.g. #NcGalaxyShapePopBeta, whose @data->ldata is invariant across
   * the whole call) can override it to amortize per-call overhead and let
   * the compiler vectorize the straight-line loop.
   */
  void (*eval_p_array) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const GArray *r, GArray **p);

  /*
   * Radial moment M_2k = <|chi_I|^2k> = int_0^1 r^2k P_pop(r) dr, where
   * P_pop(r) is the density eval_p() evaluates, normalized to
   * int_0^1 P_pop(r) dr = 1 with the radial Jacobian included (see
   * nc_galaxy_shape_pop_eval_p()). The base class default is a fixed
   * Gauss-Legendre quadrature over eval_p_array(), valid for any
   * population; a subclass whose moments have a closed form (the truncated
   * Gaussian's regularized incomplete gamma, the Beta population's
   * Beta-function ratio) overrides it. Reads the resolved parameters from
   * @data as eval_p() does, so nc_galaxy_shape_pop_prepare() must have run
   * on @data since the population's parameters last changed;
   * nc_galaxy_shape_factor_update_data_pop() and the array-prepare paths
   * that call it do so on every pass.
   */
  gdouble (*moment_2k) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const guint k);

  /*
   * The exponent alpha_1 of P_pop (r) ~ (1 - r)^alpha_1 as r -> 1, the counterpart of
   * exponent_at_origin. The base class default is an error.
   */
  gdouble (*exponent_at_one) (NcGalaxyShapePop *gsp);

  /*
   * P_pop (r) / (r^alpha_o (1 - r)^alpha_1), with alpha_o and alpha_1 the exponents of
   * exponent_at_origin and exponent_at_one, computed without the division, so that it is
   * finite at r = 0 and r = 1. Reads the resolved parameters from @data as eval_p() does.
   * Consumed by a quadrature that has the two power laws in its weight. The base class
   * default is an error.
   */
  gdouble (*eval_p_reduced) (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r);

  /* Padding to allow adding up to 8 more virtual functions without breaking ABI. */
  gpointer padding[8];
};

/**
 * NcGalaxyShapePopData:
 * @e_rms: per-galaxy intrinsic shape scatter input (used by per-galaxy models;
 *   ignored by global models).
 *
 * Per-galaxy data for the intrinsic ellipticity distribution. Mirrors the
 * #NcGalaxyShapeFactorData / #NcGalaxyPositionFactorData idiom: typed public
 * fields plus an opaque @ldata holding the subclass-specific,
 * varying/updatable resolved parameters (filled by
 * nc_galaxy_shape_pop_prepare()). The model is the factory
 * (nc_galaxy_shape_pop_data_new()).
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
gdouble nc_galaxy_shape_pop_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r);
void nc_galaxy_shape_pop_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
gdouble nc_galaxy_shape_pop_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
gdouble nc_galaxy_shape_pop_exponent_at_origin (NcGalaxyShapePop *gsp);

gdouble nc_galaxy_shape_pop_get_sigma (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);

void nc_galaxy_shape_pop_eval_p_array (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const GArray *r, GArray **p);
gdouble nc_galaxy_shape_pop_moment_2k (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const guint k);
gdouble nc_galaxy_shape_pop_exponent_at_one (NcGalaxyShapePop *gsp);
gdouble nc_galaxy_shape_pop_eval_p_reduced (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_POP_H_ */

