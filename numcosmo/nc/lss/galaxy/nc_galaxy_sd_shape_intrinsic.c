/***************************************************************************
 *            nc_galaxy_sd_shape_intrinsic.c
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_intrinsic.c
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

/**
 * NcGalaxySDShapeIntrinsic:
 *
 * Abstract model for the intrinsic galaxy ellipticity distribution.
 *
 * This small #NcmModel describes the probability density of the intrinsic
 * ellipticity modulus squared $x = |\chi_I|^2 \in [0,1]$. It is a pure
 * parametric law for which the shape calculator (#NcGalaxyWLShapeCalc) performs
 * the exact marginalization over the intrinsic ellipticity.
 *
 * The density factorizes as $P(x) = (1-x)^a\, x^b\, R(x)$, where $(1-x)^a x^b$
 * is a Gauss-Jacobi weight (so the radial integral maps naturally to
 * Gauss-Jacobi quadrature) and $R(x)$ is a smooth residual. Following the
 * #NcGalaxySDShape / #NcGalaxySDPosition idiom, the resolved per-galaxy state
 * lives in a #NcGalaxySDShapeIntrinsicData: typed public fields (@jacobi_a,
 * @jacobi_b, @e_rms) plus an opaque @ldata holding the subclass-specific,
 * varying/updatable residual parameters. nc_galaxy_sd_shape_intrinsic_prepare()
 * resolves the model parameters (and any per-galaxy @e_rms) into that data.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_sd_shape_intrinsic.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_ABSTRACT_TYPE (NcGalaxySDShapeIntrinsic, nc_galaxy_sd_shape_intrinsic, NCM_TYPE_MODEL);
G_DEFINE_BOXED_TYPE (NcGalaxySDShapeIntrinsicData, nc_galaxy_sd_shape_intrinsic_data, nc_galaxy_sd_shape_intrinsic_data_ref, nc_galaxy_sd_shape_intrinsic_data_unref); /* LCOV_EXCL_LINE */

static void
nc_galaxy_sd_shape_intrinsic_init (NcGalaxySDShapeIntrinsic *gsi)
{
}

static void
_nc_galaxy_sd_shape_intrinsic_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_shape_intrinsic_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_shape_intrinsic, NC_TYPE_GALAXY_SD_SHAPE_INTRINSIC);

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_shape_intrinsic_data_init (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  g_error ("_nc_galaxy_sd_shape_intrinsic_data_init: method not implemented.");
}

static void
_nc_galaxy_sd_shape_intrinsic_prepare (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  g_error ("_nc_galaxy_sd_shape_intrinsic_prepare: method not implemented.");
}

static gdouble
_nc_galaxy_sd_shape_intrinsic_eval_residual (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x)
{
  g_error ("_nc_galaxy_sd_shape_intrinsic_eval_residual: method not implemented.");

  return 0.0;
}

static void
_nc_galaxy_sd_shape_intrinsic_gen (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  g_error ("_nc_galaxy_sd_shape_intrinsic_gen: method not implemented.");
}

static gdouble
_nc_galaxy_sd_shape_intrinsic_e_rms (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  g_error ("_nc_galaxy_sd_shape_intrinsic_e_rms: method not implemented.");

  return 0.0;
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_sd_shape_intrinsic_class_init (NcGalaxySDShapeIntrinsicClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_sd_shape_intrinsic_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy intrinsic ellipticity distribution", "GalaxySDShapeIntrinsic");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_mset_model_register_id (model_class, "NcGalaxySDShapeIntrinsic", "Galaxy intrinsic ellipticity distribution", NULL, FALSE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->data_init     = &_nc_galaxy_sd_shape_intrinsic_data_init;
  klass->prepare       = &_nc_galaxy_sd_shape_intrinsic_prepare;
  klass->eval_residual = &_nc_galaxy_sd_shape_intrinsic_eval_residual;
  klass->gen           = &_nc_galaxy_sd_shape_intrinsic_gen;
  klass->e_rms         = &_nc_galaxy_sd_shape_intrinsic_e_rms;
}

/**
 * nc_galaxy_sd_shape_intrinsic_ref:
 * @gsi: a #NcGalaxySDShapeIntrinsic
 *
 * Increases the reference count of @gsi by one.
 *
 * Returns: (transfer full): @gsi.
 */
NcGalaxySDShapeIntrinsic *
nc_galaxy_sd_shape_intrinsic_ref (NcGalaxySDShapeIntrinsic *gsi)
{
  return g_object_ref (gsi);
}

/**
 * nc_galaxy_sd_shape_intrinsic_free:
 * @gsi: a #NcGalaxySDShapeIntrinsic
 *
 * Decreases the reference count of @gsi by one.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_free (NcGalaxySDShapeIntrinsic *gsi)
{
  g_object_unref (gsi);
}

/**
 * nc_galaxy_sd_shape_intrinsic_clear:
 * @gsi: a #NcGalaxySDShapeIntrinsic
 *
 * Decreases the reference count of *@gsi by one, and sets the pointer *@gsi to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_clear (NcGalaxySDShapeIntrinsic **gsi)
{
  g_clear_object (gsi);
}

/**
 * nc_galaxy_sd_shape_intrinsic_data_new:
 * @gsi: a #NcGalaxySDShapeIntrinsic
 *
 * Creates a new per-galaxy #NcGalaxySDShapeIntrinsicData for @gsi.
 *
 * Returns: (transfer full): a new #NcGalaxySDShapeIntrinsicData.
 */
NcGalaxySDShapeIntrinsicData *
nc_galaxy_sd_shape_intrinsic_data_new (NcGalaxySDShapeIntrinsic *gsi)
{
  NcGalaxySDShapeIntrinsicData *data = g_new0 (NcGalaxySDShapeIntrinsicData, 1);

  data->e_rms                  = 0.0;
  data->jacobi_a               = 0.0;
  data->jacobi_b               = 0.0;
  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_SD_SHAPE_INTRINSIC_GET_CLASS (gsi)->data_init (gsi, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}

/**
 * nc_galaxy_sd_shape_intrinsic_data_ref:
 * @data: a #NcGalaxySDShapeIntrinsicData
 *
 * Increases the reference count of @data by one.
 *
 * Returns: (transfer full): @data.
 */
NcGalaxySDShapeIntrinsicData *
nc_galaxy_sd_shape_intrinsic_data_ref (NcGalaxySDShapeIntrinsicData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_sd_shape_intrinsic_data_unref:
 * @data: a #NcGalaxySDShapeIntrinsicData
 *
 * Decreases the reference count of @data by one. If it reaches 0, @data is freed.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_data_unref (NcGalaxySDShapeIntrinsicData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    g_free (data);
  }
}

/**
 * nc_galaxy_sd_shape_intrinsic_data_read_row:
 * @data: a #NcGalaxySDShapeIntrinsicData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the intrinsic-distribution inputs of row @i.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_data_read_row (NcGalaxySDShapeIntrinsicData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->ldata_read_row (data, obs, i);
}

/**
 * nc_galaxy_sd_shape_intrinsic_data_write_row:
 * @data: a #NcGalaxySDShapeIntrinsicData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the intrinsic-distribution inputs of row @i.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_data_write_row (NcGalaxySDShapeIntrinsicData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->ldata_write_row (data, obs, i);
}

/**
 * nc_galaxy_sd_shape_intrinsic_data_required_columns:
 * @data: a #NcGalaxySDShapeIntrinsicData
 *
 * Returns: (element-type utf8) (transfer full): the required columns.
 */
GList *
nc_galaxy_sd_shape_intrinsic_data_required_columns (NcGalaxySDShapeIntrinsicData *data)
{
  GList *columns = NULL;

  data->ldata_required_columns (data, columns);

  return columns;
}

/**
 * nc_galaxy_sd_shape_intrinsic_prepare: (virtual prepare)
 * @gsi: a #NcGalaxySDShapeIntrinsic
 * @data: a #NcGalaxySDShapeIntrinsicData
 *
 * Resolves the model parameters (and the per-galaxy @data->e_rms, for per-galaxy
 * models) into @data: the Gauss-Jacobi exponents @data->jacobi_a, @data->jacobi_b
 * and the subclass-specific residual parameters stored in @data->ldata.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_prepare (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  NC_GALAXY_SD_SHAPE_INTRINSIC_GET_CLASS (gsi)->prepare (gsi, data);
}

/**
 * nc_galaxy_sd_shape_intrinsic_eval_residual: (virtual eval_residual)
 * @gsi: a #NcGalaxySDShapeIntrinsic
 * @data: a resolved #NcGalaxySDShapeIntrinsicData
 * @x: the squared ellipticity modulus $x = |\chi_I|^2 \in [0,1]$
 *
 * Evaluates the smooth residual $R(x)$ with $P(x) = (1-x)^a x^b R(x)$, reading
 * the resolved parameters from @data (no live parameter access).
 *
 * Returns: the residual $R(x)$.
 */
gdouble
nc_galaxy_sd_shape_intrinsic_eval_residual (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x)
{
  return NC_GALAXY_SD_SHAPE_INTRINSIC_GET_CLASS (gsi)->eval_residual (gsi, data, x);
}

/**
 * nc_galaxy_sd_shape_intrinsic_eval_p:
 * @gsi: a #NcGalaxySDShapeIntrinsic
 * @data: a resolved #NcGalaxySDShapeIntrinsicData
 * @x: the squared ellipticity modulus $x = |\chi_I|^2 \in [0,1]$
 *
 * Evaluates the full density $P(x) = (1-x)^a x^b R(x)$ of the squared intrinsic
 * ellipticity modulus (diagnostics/validation).
 *
 * Returns: the density $P(x)$.
 */
gdouble
nc_galaxy_sd_shape_intrinsic_eval_p (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x)
{
  const gdouble w = pow (1.0 - x, data->jacobi_a) * pow (x, data->jacobi_b);

  return w * nc_galaxy_sd_shape_intrinsic_eval_residual (gsi, data, x);
}

/**
 * nc_galaxy_sd_shape_intrinsic_gen: (virtual gen)
 * @gsi: a #NcGalaxySDShapeIntrinsic
 * @data: a resolved #NcGalaxySDShapeIntrinsicData
 * @rng: a #NcmRNG
 * @e_int_1: (out): first component of the sampled intrinsic ellipticity
 * @e_int_2: (out): second component of the sampled intrinsic ellipticity
 *
 * Samples an intrinsic ellipticity $\chi_I = e_1 + i e_2$ from @data.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_gen (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  NC_GALAXY_SD_SHAPE_INTRINSIC_GET_CLASS (gsi)->gen (gsi, data, rng, e_int_1, e_int_2);
}

/**
 * nc_galaxy_sd_shape_intrinsic_e_rms: (virtual e_rms)
 * @gsi: a #NcGalaxySDShapeIntrinsic
 * @data: a resolved #NcGalaxySDShapeIntrinsicData
 *
 * Computes the per-component intrinsic RMS ellipticity
 * $e_\mathrm{rms} = \sqrt{\langle |\chi_I|^2 \rangle / 2}$ implied by @data.
 *
 * Returns: the intrinsic RMS ellipticity.
 */
gdouble
nc_galaxy_sd_shape_intrinsic_e_rms (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  return NC_GALAXY_SD_SHAPE_INTRINSIC_GET_CLASS (gsi)->e_rms (gsi, data);
}
