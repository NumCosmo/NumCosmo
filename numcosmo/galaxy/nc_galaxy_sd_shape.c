/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape.c
 *
 *  Sat May 21 20:43:32 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape.c
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

/**
 * NcGalaxySDShape:
 *
 * Class describing galaxy sample shape distribution.
 *
 * This class describes a galaxy sample shape distribution. It is composed by a
 * distribution $P(s)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_wl_obs.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include "nc_enum_types.h"
#include "nc_hicosmo.h"
#include "lss/nc_halo_position.h"
#include "lss/nc_halo_density_profile.h"
#include "lss/nc_wl_surface_mass_density.h"
#include "math/ncm_rng.h"
#include "math/ncm_vector.h"

typedef struct _NcGalaxySDShapePrivate
{
  void (*apply_shear) (NcGalaxySDShape *, const NcmComplex *, const NcmComplex *, NcmComplex *);
  void (*apply_shear_inv) (NcGalaxySDShape *, const NcmComplex *, const NcmComplex *, NcmComplex *);
  gdouble (*int_to_obs_lndet_jac) (NcGalaxySDShape *, const NcmComplex *, const NcmComplex *);

  NcGalaxyWLObsEllipConv ellip_conv;
} NcGalaxySDShapePrivate;

enum
{
  PROP_0,
  PROP_ELLIP_CONV,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShape, nc_galaxy_sd_shape, NCM_TYPE_MODEL);
G_DEFINE_BOXED_TYPE (NcGalaxySDShapeData, nc_galaxy_sd_shape_data, nc_galaxy_sd_shape_data_ref, nc_galaxy_sd_shape_data_unref); /* LCOV_EXCL_LINE */
NCM_UTIL_DEFINE_CALLBACK (NcGalaxySDShapeIntegrand,
                          NC_GALAXY_SD_SHAPE_INTEGRAND,
                          nc_galaxy_sd_shape_integrand,
                          gdouble,
                          NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxySDShapeData * data),
                          NCM_UTIL_CALLBACK_ARGS (z, data))

static void
nc_galaxy_sd_shape_init (NcGalaxySDShape *gsds)
{
  NcGalaxySDShapePrivate * const self = nc_galaxy_sd_shape_get_instance_private (gsds);

  self->ellip_conv           = 0;
  self->apply_shear          = NULL;
  self->apply_shear_inv      = NULL;
  self->int_to_obs_lndet_jac = NULL;
}

static void _nc_galaxy_sd_shape_set_ellip_conv (NcGalaxySDShape *gsds, NcGalaxyWLObsEllipConv ellip_conv);

static void
_nc_galaxy_sd_shape_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShape *gsds = NC_GALAXY_SD_SHAPE (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE (gsds));

  switch (property_id)
  {
    case PROP_ELLIP_CONV:
      _nc_galaxy_sd_shape_set_ellip_conv (gsds, g_value_get_enum (value));
      break;
    default:                                                          /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                          /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_shape_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShape *gsds = NC_GALAXY_SD_SHAPE (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE (gsds));

  switch (property_id)
  {
    case PROP_ELLIP_CONV:
      g_value_set_enum (value, nc_galaxy_sd_shape_get_ellip_conv (gsds));
      break;
    default:                                                          /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                          /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_shape_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_shape_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_shape, NC_TYPE_GALAXY_SD_SHAPE);

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_shape_gen: method not implemented.");
}

static NcGalaxySDShapeIntegrand *
_nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds)
{
  g_error ("_nc_galaxy_sd_shape_integ: method not implemented.");

  return NULL;
}

static gboolean
_nc_galaxy_sd_shape_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array)
{
  g_error ("_nc_galaxy_sd_shape_prepare: method not implemented.");

  return FALSE;
}

static void
_nc_galaxy_sd_shape_data_init (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data, NcGalaxySDShapeData *data)
{
  g_error ("_nc_galaxy_sd_shape_data_new: method not implemented.");
}

static void
_nc_galaxy_sd_shape_direct_estimate (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gdouble *g1, gdouble *g2, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho)
{
  g_error ("_nc_galaxy_sd_shape_direct_estimate: method not implemented.");
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_sd_shape_class_init (NcGalaxySDShapeClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_shape_set_property;
  model_class->get_property = &_nc_galaxy_sd_shape_get_property;
  object_class->finalize    = &_nc_galaxy_sd_shape_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample shape distribution", "GalaxySDShape");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_mset_model_register_id (model_class, "NcGalaxySDShape", "Galaxy sample shape distribution", NULL, FALSE, NCM_MSET_MODEL_MAIN);

  /**
   * NcGalaxySDShape:ellip_conv:
   *
   * Weak lensing observables ellipse type #NcGalaxyWLObsEllipConv.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ELLIP_CONV,
                                   g_param_spec_enum ("ellip-conv",
                                                      "Ellipticity convention",
                                                      "Weak lensing observables ellipticity convention",
                                                      NC_TYPE_GALAXY_WL_OBS_ELLIP_CONV,
                                                      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));
  ncm_model_class_check_params_info (model_class);

  klass->gen                = &_nc_galaxy_sd_shape_gen;
  klass->integ              = &_nc_galaxy_sd_shape_integ;
  klass->prepare_data_array = &_nc_galaxy_sd_shape_prepare_data_array;
  klass->data_init          = &_nc_galaxy_sd_shape_data_init;
  klass->direct_estimate    = &_nc_galaxy_sd_shape_direct_estimate;
}

/* LCOV_EXCL_START */

/**
 * nc_galaxy_sd_shape_data_ref:
 * @data: a #NcGalaxySDShapeData
 *
 * Increases the reference count of @data by one.
 *
 */
NcGalaxySDShapeData *
nc_galaxy_sd_shape_data_ref (NcGalaxySDShapeData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/* LCOV_EXCL_STOP */

/**
 * nc_galaxy_sd_shape_data_unref:
 * @data: a #NcGalaxySDShapeData
 *
 * Decreases the reference count of @data by one. If the reference count reaches 0, the
 * data is freed.
 *
 */
void
nc_galaxy_sd_shape_data_unref (NcGalaxySDShapeData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    nc_galaxy_sd_position_data_unref (data->sdpos_data);
    g_free (data);
  }
}

static void
_nc_galaxy_sd_shape_apply_shear_trace (NcGalaxySDShape *gsds, const NcmComplex *g, const NcmComplex *chi, NcmComplex *chi_obs)
{
  complex double gn   = ncm_complex_c (g);
  complex double chin = ncm_complex_c (chi);
  complex double chin_obs;

  chin_obs = (chin + gn * (gn * conj (chin) + 2.0)) / (1.0 + gn * conj (gn) + 2.0 * creal (gn * conj (chin)));

  ncm_complex_set_c (chi_obs, chin_obs);
}

static void
_nc_galaxy_sd_shape_apply_shear_inv_trace (NcGalaxySDShape *gsds, const NcmComplex *g, const NcmComplex *chi_obs, NcmComplex *chi)
{
  complex double gn       = ncm_complex_c (g);
  complex double chin_obs = ncm_complex_c (chi_obs);
  complex double chin;

  chin = (chin_obs + gn * (gn * conj (chin_obs) - 2.0)) / (1.0 + gn * conj (gn) - 2.0 * creal (gn * conj (chin_obs)));

  ncm_complex_set_c (chi, chin);
}

static gdouble
_nc_galaxy_sd_shape_apply_shear_trace_lndet_jac (NcGalaxySDShape *gsds, const NcmComplex *g, const NcmComplex *chi_obs)
{
  complex double gn           = ncm_complex_c (g);
  complex double chin_obs     = ncm_complex_c (chi_obs);
  complex double gn_conj      = conj (gn);
  const gdouble abs_g2        = gn * gn_conj;
  const gdouble lndet_jac_den = 3.0 * log (fabs (1 - 2.0 * creal (gn_conj * chin_obs) + abs_g2));

  if (abs_g2 <= 1.0)
  {
    const gdouble lndet_jac_num = 3.0 * log1p (-abs_g2);

    return lndet_jac_num - lndet_jac_den;
  }
  else
  {
    const gdouble lndet_jac_num = 3.0 * log (abs_g2 - 1.0);

    return lndet_jac_num - lndet_jac_den;
  }
}

static void
_nc_galaxy_sd_shape_apply_shear_trace_det (NcGalaxySDShape *gsds, const NcmComplex *g, const NcmComplex *e, NcmComplex *e_obs)
{
  complex double gn = ncm_complex_c (g);
  complex double en = ncm_complex_c (e);
  complex double en_obs;

  if (cabs (gn) <= 1.0)
    en_obs = (en + gn) / (1.0 + conj (gn) * en);
  else
    en_obs = (1.0 + gn * conj (en)) / (conj (en) + conj (gn));

  ncm_complex_set_c (e_obs, en_obs);
}

static void
_nc_galaxy_sd_shape_apply_shear_inv_trace_det (NcGalaxySDShape *gsds, const NcmComplex *g, const NcmComplex *e_obs, NcmComplex *e)
{
  complex double gn     = ncm_complex_c (g);
  complex double en_obs = ncm_complex_c (e_obs);
  complex double en;

  if (cabs (gn) <= 1.0)
    en = (en_obs - gn) / (1.0 - conj (gn) * en_obs);
  else
    en = (1.0 - gn * conj (en_obs)) / (conj (en_obs) - conj (gn));

  ncm_complex_set_c (e, en);
}

static gdouble
_nc_galaxy_sd_shape_apply_shear_trace_det_lndet_jac (NcGalaxySDShape *gsds, const NcmComplex *g, const NcmComplex *e_obs)
{
  complex double gn          = ncm_complex_c (g);
  complex double en_obs      = ncm_complex_c (e_obs);
  complex double gn_conj     = conj (gn);
  complex double en_obs_conj = conj (en_obs);
  const gdouble abs_g2       = gn * gn_conj;

  if (abs_g2 <= 1.0)
  {
    const gdouble abs_en_obs2 = en_obs * en_obs_conj;
    const gdouble ln_jac_num  = 2.0 * log1p (-abs_g2);
    const gdouble ln_jac_den  = 2.0 * log1p (-2.0 * creal (gn_conj * en_obs) + abs_g2 * abs_en_obs2);

    return ln_jac_num - ln_jac_den;
  }
  else
  {
    const complex double en_obs_m_gn = en_obs - gn;
    const gdouble abs_en_obs_m_gn2   = en_obs_m_gn * conj (en_obs_m_gn);
    const gdouble ln_jac_num         = 2.0 * log (abs_g2 - 1.0);
    const gdouble ln_jac_den         = 2.0 * log (abs_en_obs_m_gn2);

    return ln_jac_num - ln_jac_den;
  }
}

static void
_nc_galaxy_sd_shape_set_ellip_conv (NcGalaxySDShape *gsds, NcGalaxyWLObsEllipConv ellip_conv)
{
  NcGalaxySDShapePrivate * const self = nc_galaxy_sd_shape_get_instance_private (gsds);

  switch (ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      self->apply_shear          = _nc_galaxy_sd_shape_apply_shear_trace;
      self->apply_shear_inv      = _nc_galaxy_sd_shape_apply_shear_inv_trace;
      self->int_to_obs_lndet_jac = _nc_galaxy_sd_shape_apply_shear_trace_lndet_jac;
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      self->apply_shear          = _nc_galaxy_sd_shape_apply_shear_trace_det;
      self->apply_shear_inv      = _nc_galaxy_sd_shape_apply_shear_inv_trace_det;
      self->int_to_obs_lndet_jac = _nc_galaxy_sd_shape_apply_shear_trace_det_lndet_jac;
      break;
    default:                                                                                             /* LCOV_EXCL_LINE */
      g_error ("nc_galaxy_sd_shape_get_ellip_conv: ellipse type %d not implemented.", self->ellip_conv); /* LCOV_EXCL_LINE */
      break;                                                                                             /* LCOV_EXCL_LINE */
  }

  self->ellip_conv = ellip_conv;
}

/**
 * nc_galaxy_sd_shape_get_ellip_conv:
 * @gsds: a #NcGalaxySDShape
 *
 * Gets the ellipse type of @gsds.
 *
 * Returns: a #NcGalaxyWLObsEllipConv
 */
NcGalaxyWLObsEllipConv
nc_galaxy_sd_shape_get_ellip_conv (NcGalaxySDShape *gsds)
{
  NcGalaxySDShapePrivate * const self = nc_galaxy_sd_shape_get_instance_private (gsds);

  return self->ellip_conv;
}

/**
 * nc_galaxy_sd_shape_apply_shear:
 * @gsds: a #NcGalaxySDShape instance
 * @g: input reduced shear as a #NcmComplex
 * @E: input intrinsic ellipticity as a #NcmComplex
 * @E_obs: output observed ellipticity as a #NcmComplex
 *
 * Applies the reduced shear @g to the intrinsic ellipticity @E, storing the resulting
 * observed ellipticity in @E_obs. The transformation depends on the
 * #NcGalaxySDShapeEllipsisType configured in @gsds.
 */
void
nc_galaxy_sd_shape_apply_shear (NcGalaxySDShape *gsds, const NcmComplex *g, const NcmComplex *E, NcmComplex *E_obs)
{
  NcGalaxySDShapePrivate * const self = nc_galaxy_sd_shape_get_instance_private (gsds);

  self->apply_shear (gsds, g, E, E_obs);
}

/**
 * nc_galaxy_sd_shape_apply_shear_inv:
 * @gsds: a #NcGalaxySDShape instance
 * @g: input reduced shear as a #NcmComplex
 * @E_obs: input observed ellipticity as a #NcmComplex
 * @E: output intrinsic ellipticity as a #NcmComplex
 *
 * Applies the inverse shear transformation using @g to recover the intrinsic
 * ellipticity @E from the observed ellipticity @E_obs. The transformation depends
 * on the #NcGalaxySDShapeEllipsisType configured in @gsds.
 */
void
nc_galaxy_sd_shape_apply_shear_inv (NcGalaxySDShape *gsds, const NcmComplex *g, const NcmComplex *E_obs, NcmComplex *E)
{
  NcGalaxySDShapePrivate * const self = nc_galaxy_sd_shape_get_instance_private (gsds);

  self->apply_shear_inv (gsds, g, E_obs, E);
}

/**
 * nc_galaxy_sd_shape_lndet_jac:
 * @gsds: a #NcGalaxySDShape instance
 * @g: input reduced shear as a #NcmComplex
 * @E_obs: input observed ellipticity as a #NcmComplex
 *
 * Computes the natural logarithm of the absolute value of the Jacobian determinant
 * of the transformation from intrinsic to observed ellipticity.
 *
 * Returns: the log-determinant of the shear Jacobian.
 */
gdouble
nc_galaxy_sd_shape_lndet_jac (NcGalaxySDShape *gsds, const NcmComplex *g, const NcmComplex *E_obs)
{
  NcGalaxySDShapePrivate * const self = nc_galaxy_sd_shape_get_instance_private (gsds);

  return self->int_to_obs_lndet_jac (gsds, g, E_obs);
}

/**
 * nc_galaxy_sd_shape_data_read_row:
 * @data: a #NcGalaxySDShapeData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the row @i from the galaxy shape data.
 *
 */
void
nc_galaxy_sd_shape_data_read_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_sd_position_data_read_row (data->sdpos_data, obs, i);
  {
    data->coord         = nc_galaxy_wl_obs_get_coord (obs);
    data->epsilon_int_1 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
    data->epsilon_int_2 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);

    data->ldata_read_row (data, obs, i);
  }
}

/**
 * nc_galaxy_sd_shape_data_write_row:
 * @data: a #NcGalaxySDShapeData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the row @i to the galaxy shape data.
 *
 */
void
nc_galaxy_sd_shape_data_write_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_sd_position_data_write_row (data->sdpos_data, obs, i);
  {
    nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i, data->epsilon_int_1);
    nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i, data->epsilon_int_2);

    data->ldata_write_row (data, obs, i);
  }
}

/**
 * nc_galaxy_sd_shape_data_required_columns:
 * @data: a #NcGalaxySDShapeData
 *
 * Returns: (element-type utf8) (transfer full): the required columns for the galaxy shape data.
 */
GList *
nc_galaxy_sd_shape_data_required_columns (NcGalaxySDShapeData *data)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2));
  data->ldata_required_columns (data, columns);

  {
    GList *sdpos_columns = nc_galaxy_sd_position_data_required_columns (data->sdpos_data);

    columns = g_list_concat (columns, sdpos_columns);
  }

  return columns;
}

/**
 * nc_galaxy_sd_shape_data_get_radius:
 * @data: a #NcGalaxySDShapeData
 *
 * Returns: the radius of the galaxy shape data.
 */
gdouble
nc_galaxy_sd_shape_data_get_radius (NcGalaxySDShapeData *data)
{
  return data->ldata_get_radius (data);
}

/**
 * nc_galaxy_sd_shape_integrand_new:
 * @func: (scope async) (closure callback_data): a #NcGalaxySDShapeIntegrandFunc
 * @callback_data_free: (scope async) (closure callback_data): a #NcGalaxySDShapeIntegrandFreeData
 * @callback_data_copy: (scope async) (closure callback_data): a #NcGalaxySDShapeIntegrandCopyData
 * @callback_data_prepare: (scope async) (closure callback_data): a #NcGalaxySDShapeIntegrandPrepareData
 * @callback_data: a gpointer
 *
 * Creates a new galaxy shape integrand.
 *
 * Returns: (transfer full): a new NcGalaxySDShapeIntegrand object.
 */
/**
 * nc_galaxy_sd_shape_integrand_copy:
 * @callback_obj: a NcGalaxySDShapeIntegrand
 *
 * Copies the integrand for the galaxy shape data.
 *
 * Returns: (transfer full): a copy of @callback_obj
 */
/**
 * nc_galaxy_sd_shape_integrand_free:
 * @callback_obj: a NcGalaxySDShapeIntegrand
 *
 * Frees the integrand for the galaxy shape data.
 *
 */
/**
 * nc_galaxy_sd_shape_integrand_prepare:
 * @callback_obj: a NcGalaxySDShapeIntegrand
 * @mset: a #NcmMSet
 *
 * Prepares the integrand for the galaxy shape data.
 *
 */

/**
 * nc_galaxy_sd_shape_ref:
 * @gsds: a #NcGalaxySDShape
 *
 * Increases the reference count of @gsds by one.
 *
 * Returns: (transfer full): @gsds.
 */
NcGalaxySDShape *
nc_galaxy_sd_shape_ref (NcGalaxySDShape *gsds)
{
  return g_object_ref (gsds);
}

/**
 * nc_galaxy_sd_shape_free:
 * @gsds: a #NcGalaxySDShape
 *
 * Decreases the reference count of @gsds by one.
 *
 */
void
nc_galaxy_sd_shape_free (NcGalaxySDShape *gsds)
{
  g_object_unref (gsds);
}

/**
 * nc_galaxy_sd_shape_clear:
 * @gsds: a #NcGalaxySDShape
 *
 * Decreases the reference count of @gsds by one, and sets the pointer *@gsds to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_clear (NcGalaxySDShape **gsds)
{
  g_clear_object (gsds);
}

/**
 * nc_galaxy_sd_shape_gen: (virtual gen)
 * @gsds: a #NcGalaxySDShape
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDShapeData
 * @rng: a #NcmRNG
 *
 * Generates a new galaxy shape. The #NcGalaxySDShapeData object @data must be
 * initialized before calling this method.
 *
 */
void
nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng)
{
  NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->gen (gsds, mset, data, rng);
}

/**
 * nc_galaxy_sd_shape_integ: (virtual integ)
 * @gsds: a #NcGalaxySDShape
 *
 * Creates a new galaxy shape integrand.
 *
 * Returns: (transfer full): a new NcGalaxySDShapeIntegrand object.
 */
NcGalaxySDShapeIntegrand *
nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->integ (gsds);
}

/**
 * nc_galaxy_sd_shape_prepare_data_array: (virtual prepare_data_array)
 * @gsds: a #NcGalaxySDShape
 * @mset: a #NcmMSet
 * @data_array: (element-type NcGalaxySDShapeData): a #GPtrArray of #NcGalaxySDShapeData
 *
 * Prepares the matrix to compute the probability density of the observaple shape.
 *
 * Returns: TRUE if the matrix was prepared, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_shape_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->prepare_data_array (gsds, mset, data_array);
}

/**
 * nc_galaxy_sd_shape_direct_estimate: (virtual direct_estimate)
 * @gsds: a #NcGalaxySDShape
 * @mset: a #NcmMSet
 * @data_array: (element-type NcGalaxySDShapeData): a #GPtrArray of #NcGalaxySDShapeData
 * @gt: (out): the reduced tangential shear
 * @gx: (out): the cross shear
 * @sigma_t: (out): the tangential scatter
 * @sigma_x: (out): the cross scatter
 * @rho: (out): the correlation
 *
 * Computes the estimate of the galaxy shape using the direct method.
 *
 */
void
nc_galaxy_sd_shape_direct_estimate (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho)
{
  NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->direct_estimate (gsds, mset, data_array, gt, gx, sigma_t, sigma_x, rho);
}

/**
 * nc_galaxy_sd_shape_data_new:
 * @gsds: a #NcGalaxySDShape
 * @sdpos_data: a #NcGalaxySDPositionData
 *
 * Creates a new galaxy shape data.
 *
 * Returns: (transfer full): a new #NcGalaxySDShapeData object.
 */
NcGalaxySDShapeData *
nc_galaxy_sd_shape_data_new (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data)
{
  NcGalaxySDShapeData *data = g_new0 (NcGalaxySDShapeData, 1);

  data->sdpos_data             = nc_galaxy_sd_position_data_ref (sdpos_data);
  data->coord                  = NC_GALAXY_WL_OBS_COORD_CELESTIAL;
  data->epsilon_int_1          = 0.0;
  data->epsilon_int_2          = 0.0;
  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->data_init (gsds, sdpos_data, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}

