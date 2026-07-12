/***************************************************************************
 *            nc_galaxy_shape_factor.c
 *
 *  Thu Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor.c
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
 * NcGalaxyShapeFactor:
 *
 * Shape likelihood-factor calculator for the weak-lensing pipeline.
 *
 * The calculator owns the whole HSM measurement engine: the intrinsic
 * ellipticity $\chi_I$ is drawn from the #NcGalaxyShapePop resolved from the
 * #NcmMSet, mapped deterministically by the reduced-shear transformation with
 * calibration bias $\tilde g = (1+m)\,g + c$, and observed with additive
 * per-galaxy Gaussian pixel noise,
 * $\epsilon_\mathrm{obs} = f_{\tilde g}(\chi_I) + n$, $n \sim N_2(0, \sigma_\mathrm{noise}^2)$.
 *
 * The likelihood factor is the marginal over the latent intrinsic shape,
 * $$P(\epsilon_\mathrm{obs} \mid z, \ldots) = \int_{|\chi_I| < 1} d^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\!\big(\epsilon_\mathrm{obs} - f_{\tilde g}(\chi_I);
 *   \sigma_\mathrm{noise}^2\big),$$
 * evaluated per galaxy at each source redshift $z$ (the outer $z$-integral
 * belongs to the orchestrator). How this two-dimensional marginal is computed
 * is the only axis subclasses vary: each subclass is one evaluation strategy
 * implementing the eval_marginal / eval_ln_marginal hooks. Everything else -
 * generation, per-galaxy geometry caches (projected radius, optzs, critical
 * surface-density nodes), frame bookkeeping and data IO - lives here, written
 * once.
 *
 * This is a calculator, not a model: it holds no fitted parameters and does
 * not live in the #NcmMSet; the models it needs (cosmology, halo position and
 * profile, surface mass density, shape population) are resolved from the mset
 * passed to each method.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop.h"
#include "nc/background/nc_hicosmo.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "ncm/stats/ncm_stats_vec.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcGalaxyShapeFactorPrivate
{
  NcGalaxyWLObsEllipConv ellip_conv;
  NcmStatsVec *obs_stats;
} NcGalaxyShapeFactorPrivate;

/*
 * Engine-owned per-galaxy geometry caches. Opaque to both introspection and
 * the integration-method subclasses: the marginal hooks receive the reduced
 * shear and the observed ellipticity already expressed in the
 * tangential/cross frame, so nothing below ever crosses the subclass seam.
 */
typedef struct _NcGalaxyShapeFactorCData
{
  gdouble epsilon_obs_t;
  gdouble epsilon_obs_x;
  gdouble c1_rot;
  gdouble c2_rot;
  gdouble radius;
  gdouble phi;
  NcWLSurfaceMassDensityOptzs optzs;
  NcWLSurfaceMassDensityCritCache *crit_cache_arr;
  NcWLSurfaceMassDensitySigmaCache sigma_cache;
  guint crit_cache_len;
} NcGalaxyShapeFactorCData;

enum
{
  PROP_0,
  PROP_ELLIP_CONV,
  PROP_LEN,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcGalaxyShapeFactor, nc_galaxy_shape_factor, G_TYPE_OBJECT)
G_DEFINE_BOXED_TYPE (NcGalaxyShapeFactorData, nc_galaxy_shape_factor_data, nc_galaxy_shape_factor_data_ref, nc_galaxy_shape_factor_data_unref); /* LCOV_EXCL_LINE */
NCM_UTIL_DEFINE_CALLBACK (NcGalaxyShapeFactorIntegrand,
                          NC_GALAXY_SHAPE_FACTOR_INTEGRAND,
                          nc_galaxy_shape_factor_integrand,
                          gdouble,
                          NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxyShapeFactorData * data),
                          NCM_UTIL_CALLBACK_ARGS (z, data))

static void
nc_galaxy_shape_factor_init (NcGalaxyShapeFactor *gsf)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);

  self->ellip_conv = 0;
  self->obs_stats  = ncm_stats_vec_new (6, NCM_STATS_VEC_COV, FALSE);
}

static void _nc_galaxy_shape_factor_set_ellip_conv (NcGalaxyShapeFactor *gsf, NcGalaxyWLObsEllipConv ellip_conv);

static void
_nc_galaxy_shape_factor_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_FACTOR (gsf));

  switch (property_id)
  {
    case PROP_ELLIP_CONV:
      _nc_galaxy_shape_factor_set_ellip_conv (gsf, g_value_get_enum (value));
      break;
    default:                                                          /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                          /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_FACTOR (gsf));

  switch (property_id)
  {
    case PROP_ELLIP_CONV:
      g_value_set_enum (value, nc_galaxy_shape_factor_get_ellip_conv (gsf));
      break;
    default:                                                          /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                          /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_dispose (GObject *object)
{
  NcGalaxyShapeFactor *gsf                = NC_GALAXY_SHAPE_FACTOR (object);
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);

  ncm_stats_vec_clear (&self->obs_stats);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_parent_class)->dispose (object);
}

static void
_nc_galaxy_shape_factor_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_parent_class)->finalize (object);
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_shape_factor_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  g_error ("_nc_galaxy_shape_factor_data_init: method not implemented.");
}

static gdouble
_nc_galaxy_shape_factor_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  g_error ("_nc_galaxy_shape_factor_eval_marginal: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_shape_factor_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  g_error ("_nc_galaxy_shape_factor_eval_ln_marginal: method not implemented.");

  return 0.0;
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_shape_factor_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  /* Default: the evaluation strategy needs no per-galaxy setup. */
}

static void
nc_galaxy_shape_factor_class_init (NcGalaxyShapeFactorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_shape_factor_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_get_property;
  object_class->dispose      = &_nc_galaxy_shape_factor_dispose;
  object_class->finalize     = &_nc_galaxy_shape_factor_finalize;

  /**
   * NcGalaxyShapeFactor:ellip-conv:
   *
   * Weak lensing observables ellipticity convention #NcGalaxyWLObsEllipConv.
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

  klass->data_init        = &_nc_galaxy_shape_factor_data_init;
  klass->prepare          = &_nc_galaxy_shape_factor_prepare;
  klass->eval_marginal    = &_nc_galaxy_shape_factor_eval_marginal;
  klass->eval_ln_marginal = &_nc_galaxy_shape_factor_eval_ln_marginal;
}

static void
_nc_galaxy_shape_factor_set_ellip_conv (NcGalaxyShapeFactor *gsf, NcGalaxyWLObsEllipConv ellip_conv)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);

  switch (ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      break;
    default:                                                                                           /* LCOV_EXCL_LINE */
      g_error ("nc_galaxy_shape_factor_set_ellip_conv: ellipse type %d not implemented.", ellip_conv); /* LCOV_EXCL_LINE */
      break;                                                                                           /* LCOV_EXCL_LINE */
  }

  self->ellip_conv = ellip_conv;
}

/**
 * nc_galaxy_shape_factor_get_ellip_conv:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Gets the ellipticity convention of @gsf.
 *
 * Returns: a #NcGalaxyWLObsEllipConv
 */
NcGalaxyWLObsEllipConv
nc_galaxy_shape_factor_get_ellip_conv (NcGalaxyShapeFactor *gsf)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);

  return self->ellip_conv;
}

/**
 * nc_galaxy_shape_factor_apply_shear:
 * @gsf: a #NcGalaxyShapeFactor
 * @g: input reduced shear as a #NcmComplex
 * @E: input intrinsic ellipticity as a #NcmComplex
 * @E_obs: output observed ellipticity as a #NcmComplex
 *
 * Applies the reduced shear @g to the intrinsic ellipticity @E, storing the
 * resulting observed ellipticity in @E_obs. The transformation depends on the
 * #NcGalaxyWLObsEllipConv configured in @gsf.
 */
void
nc_galaxy_shape_factor_apply_shear (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E, NcmComplex *E_obs)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);
  const complex double gn                 = ncm_complex_c (g);
  const complex double En                 = ncm_complex_c (E);

  switch (self->ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      ncm_complex_set_c (E_obs, nc_wl_ellipticity_apply_shear_trace_c (gn, En));
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      ncm_complex_set_c (E_obs, nc_wl_ellipticity_apply_shear_trace_det_c (gn, En));
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_galaxy_shape_factor_apply_shear_inv:
 * @gsf: a #NcGalaxyShapeFactor
 * @g: input reduced shear as a #NcmComplex
 * @E_obs: input observed ellipticity as a #NcmComplex
 * @E: output intrinsic ellipticity as a #NcmComplex
 *
 * Applies the inverse shear transformation using @g to recover the intrinsic
 * ellipticity @E from the observed ellipticity @E_obs. The transformation
 * depends on the #NcGalaxyWLObsEllipConv configured in @gsf.
 */
void
nc_galaxy_shape_factor_apply_shear_inv (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E_obs, NcmComplex *E)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);
  const complex double gn                 = ncm_complex_c (g);
  const complex double En_obs             = ncm_complex_c (E_obs);

  switch (self->ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      ncm_complex_set_c (E, nc_wl_ellipticity_apply_shear_inv_trace_c (gn, En_obs));
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      ncm_complex_set_c (E, nc_wl_ellipticity_apply_shear_inv_trace_det_c (gn, En_obs));
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_galaxy_shape_factor_lndet_jac:
 * @gsf: a #NcGalaxyShapeFactor
 * @g: input reduced shear as a #NcmComplex
 * @E_obs: input observed ellipticity as a #NcmComplex
 *
 * Computes the natural logarithm of the absolute value of the Jacobian
 * determinant of the transformation from intrinsic to observed ellipticity.
 *
 * Returns: the log-determinant of the shear Jacobian.
 */
gdouble
nc_galaxy_shape_factor_lndet_jac (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E_obs)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);
  const complex double gn                 = ncm_complex_c (g);
  const complex double En_obs             = ncm_complex_c (E_obs);

  switch (self->ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      return nc_wl_ellipticity_lndet_jac_trace_c (gn, En_obs);

    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      return nc_wl_ellipticity_lndet_jac_trace_det_c (gn, En_obs);

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

/* LCOV_EXCL_START */

/**
 * nc_galaxy_shape_factor_data_ref:
 * @data: a #NcGalaxyShapeFactorData
 *
 * Increases the reference count of @data by one.
 *
 * Returns: (transfer full): @data.
 */
NcGalaxyShapeFactorData *
nc_galaxy_shape_factor_data_ref (NcGalaxyShapeFactorData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/* LCOV_EXCL_STOP */

/**
 * nc_galaxy_shape_factor_data_unref:
 * @data: a #NcGalaxyShapeFactorData
 *
 * Decreases the reference count of @data by one. If the reference count
 * reaches 0, the data is freed.
 *
 */
void
nc_galaxy_shape_factor_data_unref (NcGalaxyShapeFactorData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    NcGalaxyShapeFactorCData *cdata = (NcGalaxyShapeFactorCData *) data->cdata;

    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    g_clear_pointer (&cdata->crit_cache_arr, g_free);
    g_free (cdata);
    nc_galaxy_shape_pop_data_unref (data->pop_data);
    nc_galaxy_position_factor_data_unref (data->pos_data);
    nc_galaxy_redshift_factor_data_unref (data->z_data);
    g_free (data);
  }
}

/**
 * nc_galaxy_shape_factor_data_read_row:
 * @data: a #NcGalaxyShapeFactorData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the row @i from @obs into @data, cascading to the upstream position
 * and redshift fragments, the population fragment and the subclass fragment.
 *
 */
void
nc_galaxy_shape_factor_data_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_position_factor_data_read_row (data->pos_data, obs, i);
  nc_galaxy_redshift_factor_data_read_row (data->z_data, obs, i);

  data->coord         = nc_galaxy_wl_obs_get_coord (obs);
  data->epsilon_int_1 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_1, i, NULL);
  data->epsilon_int_2 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_2, i, NULL);
  data->epsilon_obs_1 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1, i, NULL);
  data->epsilon_obs_2 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2, i, NULL);
  data->std_noise     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE, i, NULL);
  data->c1            = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_C1, i, NULL);
  data->c2            = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_C2, i, NULL);
  data->m             = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_M, i, NULL);

  nc_galaxy_shape_pop_data_read_row (data->pop_data, obs, i);
  data->ldata_read_row (data, obs, i);
}

/**
 * nc_galaxy_shape_factor_data_write_row:
 * @data: a #NcGalaxyShapeFactorData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the row @i of @data into @obs, cascading to the upstream position
 * and redshift fragments, the population fragment and the subclass fragment.
 *
 */
void
nc_galaxy_shape_factor_data_write_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_position_factor_data_write_row (data->pos_data, obs, i);
  nc_galaxy_redshift_factor_data_write_row (data->z_data, obs, i);

  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_1, i, data->epsilon_int_1, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_2, i, data->epsilon_int_2, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1, i, data->epsilon_obs_1, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2, i, data->epsilon_obs_2, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE, i, data->std_noise, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_C1, i, data->c1, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_C2, i, data->c2, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_M, i, data->m, NULL);

  nc_galaxy_shape_pop_data_write_row (data->pop_data, obs, i);
  data->ldata_write_row (data, obs, i);
}

/**
 * nc_galaxy_shape_factor_data_required_columns:
 * @data: a #NcGalaxyShapeFactorData
 *
 * Returns: (element-type utf8) (transfer full): the required columns for the galaxy shape data.
 */
GList *
nc_galaxy_shape_factor_data_required_columns (NcGalaxyShapeFactorData *data)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_C1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_C2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_M));

  columns = g_list_concat (columns, nc_galaxy_shape_pop_data_required_columns (data->pop_data));
  data->ldata_required_columns (data, &columns);

  columns = g_list_concat (columns, nc_galaxy_position_factor_data_required_columns (data->pos_data));
  columns = g_list_concat (columns, nc_galaxy_redshift_factor_data_required_columns (data->z_data));

  return columns;
}

/**
 * nc_galaxy_shape_factor_data_get_radius:
 * @data: a #NcGalaxyShapeFactorData
 *
 * Returns: the projected radius cached for this galaxy.
 */
gdouble
nc_galaxy_shape_factor_data_get_radius (NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorCData *cdata = (NcGalaxyShapeFactorCData *) data->cdata;

  return cdata->radius;
}

/**
 * nc_galaxy_shape_factor_ref:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Increases the reference count of @gsf by one.
 *
 * Returns: (transfer full): @gsf.
 */
NcGalaxyShapeFactor *
nc_galaxy_shape_factor_ref (NcGalaxyShapeFactor *gsf)
{
  return g_object_ref (gsf);
}

/**
 * nc_galaxy_shape_factor_free:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Decreases the reference count of @gsf by one.
 *
 */
void
nc_galaxy_shape_factor_free (NcGalaxyShapeFactor *gsf)
{
  g_object_unref (gsf);
}

/**
 * nc_galaxy_shape_factor_clear:
 * @gsf: a #NcGalaxyShapeFactor
 *
 * Decreases the reference count of @gsf by one, and sets the pointer *@gsf to
 * NULL.
 *
 */
void
nc_galaxy_shape_factor_clear (NcGalaxyShapeFactor **gsf)
{
  g_clear_object (gsf);
}

/**
 * nc_galaxy_shape_factor_data_new:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @pos_data: the upstream #NcGalaxyPositionFactorData
 * @z_data: the upstream #NcGalaxyRedshiftFactorData
 *
 * Creates a new per-galaxy shape data. The intrinsic-ellipticity fragment is
 * allocated by the #NcGalaxyShapePop resolved from @mset; the upstream
 * fragments are referenced, not copied.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorData object.
 */
NcGalaxyShapeFactorData *
nc_galaxy_shape_factor_data_new (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyPositionFactorData *pos_data, NcGalaxyRedshiftFactorData *z_data)
{
  NcGalaxyShapeFactorData *data = g_new0 (NcGalaxyShapeFactorData, 1);
  NcGalaxyShapePop *pop         = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

  g_assert_nonnull (pop);

  data->pos_data = nc_galaxy_position_factor_data_ref (pos_data);
  data->z_data   = nc_galaxy_redshift_factor_data_ref (z_data);
  data->pop_data = nc_galaxy_shape_pop_data_new (pop);
  data->coord    = NC_WL_ELLIPTICITY_FRAME_CELESTIAL;
  data->cdata    = g_new0 (NcGalaxyShapeFactorCData, 1);

  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf)->data_init (gsf, mset, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}

/**
 * nc_galaxy_shape_factor_data_set:
 * @gsf: a #NcGalaxyShapeFactor
 * @data: a #NcGalaxyShapeFactorData
 * @epsilon_obs_1: the observed ellipticity component 1
 * @epsilon_obs_2: the observed ellipticity component 2
 * @std_noise: the observational shape dispersion
 * @c1: the first additive bias parameter (in the @coord frame)
 * @c2: the second additive bias parameter (in the @coord frame)
 * @m: the multiplicative bias parameter
 * @coord: the ellipticity handedness frame #NcWLEllipticityFrame
 *
 * Sets the observed ellipticity and the per-galaxy measurement quantities.
 *
 */
void
nc_galaxy_shape_factor_data_set (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m, NcWLEllipticityFrame coord)
{
  data->epsilon_obs_1 = epsilon_obs_1;
  data->epsilon_obs_2 = epsilon_obs_2;
  data->std_noise     = std_noise;
  data->c1            = c1;
  data->c2            = c2;
  data->m             = m;
  data->coord         = coord;
}

/**
 * nc_galaxy_shape_factor_data_get:
 * @gsf: a #NcGalaxyShapeFactor
 * @data: a #NcGalaxyShapeFactorData
 * @epsilon_obs_1: (out): the observed ellipticity component 1
 * @epsilon_obs_2: (out): the observed ellipticity component 2
 * @std_noise: (out): the observational shape dispersion
 * @c1: (out): the first additive bias parameter
 * @c2: (out): the second additive bias parameter
 * @m: (out): the multiplicative bias parameter
 *
 * Gets the observed ellipticity and the per-galaxy measurement quantities.
 *
 */
void
nc_galaxy_shape_factor_data_get (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *std_noise, gdouble *c1, gdouble *c2, gdouble *m)
{
  *epsilon_obs_1 = data->epsilon_obs_1;
  *epsilon_obs_2 = data->epsilon_obs_2;
  *std_noise     = data->std_noise;
  *c1            = data->c1;
  *c2            = data->c2;
  *m             = data->m;
}

/**
 * nc_galaxy_shape_factor_gen:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @data: a #NcGalaxyShapeFactorData
 * @rng: a #NcmRNG
 *
 * Generates a galaxy shape: the intrinsic ellipticity is drawn from the
 * #NcGalaxyShapePop in @mset, sheared with the reduced shear at the galaxy's
 * (ra, dec, z) with the calibration bias applied, and observed with additive
 * Gaussian noise. The per-galaxy measurement quantities (std_noise, c1, c2, m)
 * and the frame @data->coord must be set beforehand (see
 * nc_galaxy_shape_factor_data_set()). The shape is generated in the celestial
 * frame and mapped to @data->coord (see #NcWLEllipticityFrame); since the
 * intrinsic shape and noise are isotropic the generation is frame-symmetric,
 * while the additive bias (c1, c2) is taken to be already expressed in
 * @data->coord.
 *
 */
void
nc_galaxy_shape_factor_gen (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, NcmRNG *rng)
{
  NcGalaxyShapeFactorCData *cdata              = (NcGalaxyShapeFactorCData *) data->cdata;
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  NcGalaxyShapePop *pop                        = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  const gdouble z_cl                           = nc_halo_position_get_redshift (halo_position);
  const gdouble ra                             = data->pos_data->ra;
  const gdouble dec                            = data->pos_data->dec;
  const gdouble z                              = data->z_data->z;
  const gdouble noise1                         = ncm_rng_gaussian_gen (rng, 0.0, data->std_noise);
  const gdouble noise2                         = ncm_rng_gaussian_gen (rng, 0.0, data->std_noise);
  const gdouble c1                             = data->c1;
  const gdouble c2                             = data->c2;
  const gdouble m                              = data->m;
  complex double noise                         = noise1 + I * noise2;
  complex double c                             = c1 + I * c2;
  gdouble theta                                = 0.0;
  gdouble phi                                  = 0.0;
  gdouble e_int_1, e_int_2;
  complex double e_s, e_o;
  gdouble radius;

  g_assert_nonnull (pop);

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (surface_mass_density, cosmo);

  nc_galaxy_shape_pop_prepare (pop, data->pop_data);
  nc_galaxy_shape_pop_gen (pop, data->pop_data, rng, &e_int_1, &e_int_2);
  e_s = e_int_1 + I * e_int_2;

  nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

  /* polar_angles returns the celestial position angle phi_C. The shear, the
   * intrinsic ellipticity e_s and the noise are all built here in the celestial
   * frame and then mapped into the data's frame (data->coord) so that the stored
   * observed ellipticity lands in that frame: the position angle is re-expressed
   * (phi_C -> phi_E) and the spin-2 quantities are conjugated. Both helpers are
   * the identity for CELESTIAL. e_s and noise are drawn isotropically, so the
   * parity flip leaves their distribution unchanged (the generation is
   * frame-symmetric); it is applied only to keep the deterministic pieces (phi
   * and the resulting e_o) consistent with data->coord. The additive bias c is
   * taken to be already expressed in data->coord, so it is added below without
   * any parity flip. See #NcWLEllipticityFrame. */
  phi   = nc_wl_ellipticity_celestial_to_frame_angle (data->coord, phi);
  e_s   = nc_wl_ellipticity_celestial_to_frame_c (data->coord, e_s);
  noise = nc_wl_ellipticity_celestial_to_frame_c (data->coord, noise);

  radius = nc_halo_position_projected_radius (halo_position, cosmo, theta);

  if (z > z_cl)
  {
    const gdouble gt = nc_wl_surface_mass_density_reduced_shear (surface_mass_density,
                                                                 density_profile,
                                                                 cosmo,
                                                                 radius, z, z_cl, z_cl);

    complex double g    = gt * cexp (2.0 * I * phi);
    NcmComplex cplx_E_s = NCM_COMPLEX_INIT (e_s);
    NcmComplex cplx_E_o, cplx_g;

    /* Adding bias */
    g = (1.0 + m) * g + c;

    ncm_complex_set_c (&cplx_g, g);
    nc_galaxy_shape_factor_apply_shear (gsf, &cplx_g, &cplx_E_s, &cplx_E_o);

    e_o = ncm_complex_c (&cplx_E_o);
  }
  else
  {
    e_o = e_s + c;
  }

  e_o += noise;

  data->epsilon_int_1 = creal (e_s);
  data->epsilon_int_2 = cimag (e_s);
  data->epsilon_obs_1 = creal (e_o);
  data->epsilon_obs_2 = cimag (e_o);
  cdata->radius       = radius;
  cdata->phi          = phi;

  /* Keep the tangential-frame caches (read by the integrand and fixed-node
   * likelihoods) consistent with the freshly generated observed ellipticity.
   * prepare_data_array() only refreshes these when the geometry changes, so a
   * resample - which rewrites epsilon_obs_1/2 in place without a model change -
   * would otherwise leave them stale. */
  {
    const complex double e_o_rotated  = e_o * cexp (-2.0 * I * phi);
    const complex double bias_rotated = c * cexp (-2.0 * I * phi);

    cdata->epsilon_obs_t = creal (e_o_rotated);
    cdata->epsilon_obs_x = cimag (e_o_rotated);
    cdata->c1_rot        = creal (bias_rotated);
    cdata->c2_rot        = cimag (bias_rotated);
  }
}

struct _IntegData
{
  NcGalaxyShapeFactor *gsf;
  NcGalaxyShapePop *pop;
  NcHICosmo *cosmo;
  NcHaloPosition *halo_position;
  NcWLSurfaceMassDensity *surface_mass_density;
  NcHaloDensityProfile *density_profile;
  gdouble (*marginal) (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
};

/* LCOV_EXCL_START */
static gpointer
_integ_data_copy (gpointer user_data)
{
  struct _IntegData *new_int_data = g_new0 (struct _IntegData, 1);
  struct _IntegData *int_data     = (struct _IntegData *) user_data;

  new_int_data->gsf                  = int_data->gsf;
  new_int_data->pop                  = nc_galaxy_shape_pop_ref (int_data->pop);
  new_int_data->cosmo                = nc_hicosmo_ref (int_data->cosmo);
  new_int_data->halo_position        = nc_halo_position_ref (int_data->halo_position);
  new_int_data->surface_mass_density = nc_wl_surface_mass_density_ref (int_data->surface_mass_density);
  new_int_data->density_profile      = nc_halo_density_profile_ref (int_data->density_profile);
  new_int_data->marginal             = int_data->marginal;

  return new_int_data;
}

/* LCOV_EXCL_STOP */

static void
_integ_data_free (gpointer user_data)
{
  struct _IntegData *int_data = (struct _IntegData *) user_data;

  nc_galaxy_shape_pop_free (int_data->pop);
  nc_hicosmo_free (int_data->cosmo);
  nc_halo_position_free (int_data->halo_position);
  nc_wl_surface_mass_density_free (int_data->surface_mass_density);
  nc_halo_density_profile_free (int_data->density_profile);

  g_free (int_data);
}

static void
_integ_data_prepare (gpointer user_data, NcmMSet *mset)
{
  struct _IntegData *int_data                  = (struct _IntegData *) user_data;
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  NcGalaxyShapePop *pop                        = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

  g_assert_nonnull (cosmo);
  g_assert_nonnull (halo_position);
  g_assert_nonnull (surface_mass_density);
  g_assert_nonnull (density_profile);
  g_assert_nonnull (pop);

  nc_galaxy_shape_pop_clear (&int_data->pop);
  nc_hicosmo_clear (&int_data->cosmo);
  nc_halo_position_clear (&int_data->halo_position);
  nc_wl_surface_mass_density_clear (&int_data->surface_mass_density);
  nc_halo_density_profile_clear (&int_data->density_profile);

  int_data->pop                  = nc_galaxy_shape_pop_ref (pop);
  int_data->cosmo                = nc_hicosmo_ref (cosmo);
  int_data->halo_position        = nc_halo_position_ref (halo_position);
  int_data->surface_mass_density = nc_wl_surface_mass_density_ref (surface_mass_density);
  int_data->density_profile      = nc_halo_density_profile_ref (density_profile);

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (surface_mass_density, cosmo);
}

static gdouble
_nc_galaxy_shape_factor_integ_f (gpointer callback_data, const gdouble z, NcGalaxyShapeFactorData *data)
{
  struct _IntegData *int_data     = (struct _IntegData *) callback_data;
  NcGalaxyShapeFactorCData *cdata = (NcGalaxyShapeFactorCData *) data->cdata;
  const gdouble z_cl              = nc_halo_position_get_redshift (int_data->halo_position);
  gdouble gt;

  /* Work in the tangential/cross frame (real reduced shear gt, observed
   * ellipticity epsilon_obs_t/x, calibration bias pre-rotated to c1_rot/c2_rot),
   * matching the fixed-node eval_at_nodes path. The likelihood is rotation
   * invariant, so this agrees with a sky-frame computation provided the bias is
   * rotated consistently with the observed ellipticity. */
  if (z > z_cl)
  {
    gt = nc_wl_surface_mass_density_reduced_shear_optzs (int_data->surface_mass_density,
                                                         int_data->density_profile,
                                                         int_data->cosmo,
                                                         z, z_cl, &cdata->optzs);
  }
  else
  {
    const gdouble gt0 = nc_wl_surface_mass_density_reduced_shear_optzs (int_data->surface_mass_density,
                                                                        int_data->density_profile,
                                                                        int_data->cosmo,
                                                                        z_cl * (1.0 + GSL_DBL_EPSILON), z_cl, &cdata->optzs);
    const gdouble step = exp ((z - z_cl) / 0.001);

    gt = step * gt0;
  }

  {
    /* Adding bias */
    const complex double g = (1.0 + data->m) * gt + (cdata->c1_rot + I * cdata->c2_rot);

    return int_data->marginal (int_data->gsf, int_data->pop, data,
                               creal (g), cimag (g),
                               cdata->epsilon_obs_t, cdata->epsilon_obs_x);
  }
}

/**
 * nc_galaxy_shape_factor_integ:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @use_lnp: if TRUE the integrand returns the natural logarithm of the probability density
 *
 * Creates a new shape integrand P(epsilon_obs | z, data). The marginalization
 * hook (eval_marginal or eval_ln_marginal) is resolved once here, out of the
 * per-evaluation path.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorIntegrand object.
 */
NcGalaxyShapeFactorIntegrand *
nc_galaxy_shape_factor_integ (NcGalaxyShapeFactor *gsf, NcmMSet *mset, gboolean use_lnp)
{
  NcGalaxyShapeFactorClass *klass = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
  struct _IntegData *int_data     = g_new0 (struct _IntegData, 1);
  NcGalaxyShapeFactorIntegrand *integ;

  integ = nc_galaxy_shape_factor_integrand_new (_nc_galaxy_shape_factor_integ_f,
                                                _integ_data_free,
                                                _integ_data_copy,
                                                _integ_data_prepare,
                                                int_data);

  int_data->gsf      = gsf;
  int_data->marginal = use_lnp ? klass->eval_ln_marginal : klass->eval_marginal;

  _integ_data_prepare (int_data, mset);

  return integ;
}

/**
 * nc_galaxy_shape_factor_prepare_data_array:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @data_array: (element-type NcGalaxyShapeFactorData): a #GPtrArray of #NcGalaxyShapeFactorData
 * @update_radius: whether the radius must be updated
 * @update_optzs: whether optzs must be updated
 *
 * Prepares the per-galaxy caches used by the integrand: lens geometry
 * (projected radius, tangential-frame rotation of the observed ellipticity and
 * calibration bias) and the reduced-shear optimization cache. Also prepares
 * the population fragments and runs the subclass per-galaxy prepare hook.
 *
 * Returns: TRUE if the caches were prepared.
 */
gboolean
nc_galaxy_shape_factor_prepare_data_array (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array, gboolean update_radius, gboolean update_optzs)
{
  NcGalaxyShapeFactorClass *klass              = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  NcGalaxyShapePop *pop                        = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  const gdouble z_cl                           = nc_halo_position_get_redshift (halo_position);
  NcWLSurfaceMassDensityLensCtx lens_ctx;
  gboolean pr_prefactor_ready = FALSE;
  gdouble r_s, rho_s, pr_prefactor = 0.0;
  guint i;

  g_assert_nonnull (pop);

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (surface_mass_density, cosmo);

  if (update_radius)
  {
    pr_prefactor       = nc_halo_position_projected_radius_prefactor (halo_position, cosmo);
    pr_prefactor_ready = TRUE;
  }

  if (update_optzs)
  {
    nc_wl_surface_mass_density_lens_ctx_prep (&lens_ctx, surface_mass_density, cosmo, z_cl);
    nc_halo_density_profile_r_s_rho_s (density_profile, cosmo, z_cl, &r_s, &rho_s);
  }

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxyShapeFactorData *data_i   = g_ptr_array_index (data_array, i);
    NcGalaxyShapeFactorCData *cdata_i = (NcGalaxyShapeFactorCData *) data_i->cdata;

    if ((update_radius) || (cdata_i->radius == 0.0))
    {
      const gdouble e1   = data_i->epsilon_obs_1;
      const gdouble e2   = data_i->epsilon_obs_2;
      const gdouble ra   = data_i->pos_data->ra;
      const gdouble dec  = data_i->pos_data->dec;
      complex double e_o = e1 + I * e2;
      complex double e_o_rotated;
      complex double bias_rotated;
      gdouble theta, phi;

      nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

      /* Re-express the celestial phi_C in the data frame (phi_E) so the spin-2
       * rotation e^{-2 i phi} carries the data-frame observed ellipticity into
       * the tangential/cross frame. See #NcWLEllipticityFrame. */
      phi = nc_wl_ellipticity_celestial_to_frame_angle (data_i->coord, phi);

      e_o_rotated  = e_o * cexp (-2.0 * I * phi);
      bias_rotated = (data_i->c1 + I * data_i->c2) * cexp (-2.0 * I * phi);

      if (!pr_prefactor_ready)
      {
        pr_prefactor       = nc_halo_position_projected_radius_prefactor (halo_position, cosmo);
        pr_prefactor_ready = TRUE;
      }

      cdata_i->radius        = nc_halo_position_projected_radius_from_prefactor (theta, pr_prefactor);
      cdata_i->phi           = phi;
      cdata_i->epsilon_obs_t = creal (e_o_rotated);
      cdata_i->epsilon_obs_x = cimag (e_o_rotated);
      cdata_i->c1_rot        = creal (bias_rotated);
      cdata_i->c2_rot        = cimag (bias_rotated);
    }

    if (update_optzs)
      nc_wl_surface_mass_density_reduced_shear_optzs_prep_with_lens_ctx (
        density_profile,
        &lens_ctx,
        cdata_i->radius,
        r_s,
        rho_s,
        &cdata_i->optzs
      );

    nc_galaxy_shape_pop_prepare (pop, data_i->pop_data);
    klass->prepare (gsf, mset, data_i);
  }

  return TRUE;
}

/**
 * nc_galaxy_shape_factor_prepare_data_array_at_nodes:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @data_array: (element-type NcGalaxyShapeFactorData): per-galaxy shape data
 * @z_nodes_per_galaxy: (element-type NcmVector): per-galaxy z node vectors
 * @update_radius: whether the radius must be updated
 * @update_crit: whether the critical surface density cache data must be updated
 * @update_sigma: whether the surface mass density cache must be updated
 *
 * Pre-computes cached quantities (Sigma_crit per node, Sigma(R)) for the fixed
 * quadrature path. This unconditionally (re)computes every per-galaxy node
 * cache: deciding *when* a recompute is needed is the caller's responsibility.
 *
 * Returns: TRUE if successful.
 */
gboolean
nc_galaxy_shape_factor_prepare_data_array_at_nodes (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array, const GPtrArray *z_nodes_per_galaxy, gboolean update_radius, gboolean update_crit, gboolean update_sigma)
{
  NcGalaxyShapeFactorClass *klass              = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  NcGalaxyShapePop *pop                        = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  const gdouble z_cl                           = nc_halo_position_get_redshift (halo_position);
  gboolean lens_ctx_ready                      = FALSE;
  gboolean pr_prefactor_ready                  = FALSE;
  NcWLSurfaceMassDensityLensCtx lens_ctx;
  gdouble pr_prefactor = 0.0;
  guint i;

  g_assert_nonnull (pop);

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (surface_mass_density, cosmo);

  if (update_radius)
  {
    pr_prefactor       = nc_halo_position_projected_radius_prefactor (halo_position, cosmo);
    pr_prefactor_ready = TRUE;
  }

  if (update_crit)
  {
    nc_wl_surface_mass_density_lens_ctx_prep (&lens_ctx, surface_mass_density, cosmo, z_cl);
    lens_ctx_ready = TRUE;
  }

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxyShapeFactorData *data_i   = g_ptr_array_index (data_array, i);
    NcGalaxyShapeFactorCData *cdata_i = (NcGalaxyShapeFactorCData *) data_i->cdata;
    NcmVector *z_nodes_i              = (NcmVector *) g_ptr_array_index ((GPtrArray *) z_nodes_per_galaxy, i);
    const guint n_nodes               = ncm_vector_len (z_nodes_i);

    if ((update_radius) || (cdata_i->radius == 0.0))
    {
      const gdouble e1   = data_i->epsilon_obs_1;
      const gdouble e2   = data_i->epsilon_obs_2;
      const gdouble ra   = data_i->pos_data->ra;
      const gdouble dec  = data_i->pos_data->dec;
      complex double e_o = e1 + I * e2;
      complex double e_o_rotated;
      complex double bias_rotated;
      gdouble theta, phi;

      nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

      /* Re-express the celestial phi_C in the data frame (phi_E) so the spin-2
       * rotation e^{-2 i phi} carries the data-frame observed ellipticity into
       * the tangential/cross frame. See #NcWLEllipticityFrame. */
      phi = nc_wl_ellipticity_celestial_to_frame_angle (data_i->coord, phi);

      e_o_rotated  = e_o * cexp (-2.0 * I * phi);
      bias_rotated = (data_i->c1 + I * data_i->c2) * cexp (-2.0 * I * phi);

      if (!pr_prefactor_ready)
      {
        pr_prefactor       = nc_halo_position_projected_radius_prefactor (halo_position, cosmo);
        pr_prefactor_ready = TRUE;
      }

      cdata_i->radius        = nc_halo_position_projected_radius_from_prefactor (theta, pr_prefactor);
      cdata_i->phi           = phi;
      cdata_i->epsilon_obs_t = creal (e_o_rotated);
      cdata_i->epsilon_obs_x = cimag (e_o_rotated);
      cdata_i->c1_rot        = creal (bias_rotated);
      cdata_i->c2_rot        = cimag (bias_rotated);
    }

    if (update_crit || (cdata_i->crit_cache_arr == NULL))
    {
      guint j;

      if (!lens_ctx_ready)
      {
        nc_wl_surface_mass_density_lens_ctx_prep (&lens_ctx, surface_mass_density, cosmo, z_cl);
        lens_ctx_ready = TRUE;
      }

      if (cdata_i->crit_cache_len != n_nodes)
      {
        g_free (cdata_i->crit_cache_arr);
        cdata_i->crit_cache_arr = g_new0 (NcWLSurfaceMassDensityCritCache, n_nodes);
        cdata_i->crit_cache_len = n_nodes;
      }

      for (j = 0; j < n_nodes; j++)
      {
        const gdouble z_j = ncm_vector_get (z_nodes_i, j);

        nc_wl_surface_mass_density_reduced_shear_crit_cache_prep_with_lens_ctx (
          surface_mass_density,
          cosmo,
          &lens_ctx,
          z_j,
          &cdata_i->crit_cache_arr[j]
        );
      }
    }

    if (update_sigma)
      nc_wl_surface_mass_density_reduced_shear_sigma_cache_prep (
        density_profile,
        cosmo,
        cdata_i->radius,
        z_cl,
        z_cl,
        &cdata_i->sigma_cache
      );

    nc_galaxy_shape_pop_prepare (pop, data_i->pop_data);
    klass->prepare (gsf, mset, data_i);
  }

  return TRUE;
}

/**
 * nc_galaxy_shape_factor_eval_at_nodes:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @data: a #NcGalaxyShapeFactorData
 * @z_nodes: the z values at each quadrature node
 * @out: a #NcmVector receiving the shape likelihood at each node
 *
 * Evaluates the shape likelihood P(epsilon_obs | z_j, data) at every node,
 * writing the result into @out. Requires a prior call to
 * nc_galaxy_shape_factor_prepare_data_array_at_nodes().
 *
 */
void
nc_galaxy_shape_factor_eval_at_nodes (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const NcmVector *z_nodes, NcmVector *out)
{
  NcGalaxyShapeFactorClass *klass = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
  NcHaloPosition *halo_position   = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcGalaxyShapePop *pop           = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  NcGalaxyShapeFactorCData *cdata = (NcGalaxyShapeFactorCData *) data->cdata;
  const gdouble z_cl              = nc_halo_position_get_redshift (halo_position);
  const gdouble et                = cdata->epsilon_obs_t;
  const gdouble ex                = cdata->epsilon_obs_x;
  const gdouble c1_rot            = cdata->c1_rot;
  const gdouble c2_rot            = cdata->c2_rot;
  const gdouble m                 = data->m;
  const guint n_nodes             = ncm_vector_len (z_nodes);
  guint j;

  g_assert_nonnull (pop);
  g_assert_nonnull (cdata->crit_cache_arr);

  for (j = 0; j < n_nodes; j++)
  {
    const gdouble z_j = ncm_vector_get (z_nodes, j);
    gdouble gt;

    if (z_j > z_cl)
      gt = nc_wl_surface_mass_density_reduced_shear_cache (&cdata->crit_cache_arr[j], &cdata->sigma_cache);
    else
      gt = 0.0;

    {
      const complex double g = (1.0 + m) * gt + (c1_rot + I * c2_rot);

      ncm_vector_set (out, j, klass->eval_marginal (gsf, pop, data, creal (g), cimag (g), et, ex));
    }
  }
}

/**
 * nc_galaxy_shape_factor_direct_estimate:
 * @gsf: a #NcGalaxyShapeFactor
 * @mset: a #NcmMSet
 * @data_array: (element-type NcGalaxyShapeFactorData): a #GPtrArray of #NcGalaxyShapeFactorData
 * @gt: (out): the reduced tangential shear
 * @gx: (out): the cross shear
 * @sigma_t: (out): the tangential scatter
 * @sigma_x: (out): the cross scatter
 * @rho: (out): the correlation
 *
 * Computes the inverse-variance-weighted direct estimate of the reduced shear
 * from the observed ellipticities. The per-galaxy intrinsic scatter enters
 * through the population model's rms (see nc_galaxy_shape_pop_e_rms()).
 *
 */
void
nc_galaxy_shape_factor_direct_estimate (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho)
{
  NcGalaxyShapeFactorPrivate * const self = nc_galaxy_shape_factor_get_instance_private (gsf);
  NcHICosmo *cosmo                        = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position           = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcGalaxyShapePop *pop                   = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  NcGalaxyWLObsEllipConv ellip_conv       = nc_galaxy_shape_factor_get_ellip_conv (gsf);
  guint i;

  g_assert_nonnull (pop);

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  ncm_stats_vec_reset (self->obs_stats, TRUE);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxyShapeFactorData *data_i = g_ptr_array_index (data_array, i);
    const gdouble ra                = data_i->pos_data->ra;
    const gdouble dec               = data_i->pos_data->dec;
    const gdouble e1                = data_i->epsilon_obs_1;
    const gdouble e2                = data_i->epsilon_obs_2;
    const gdouble std_noise         = data_i->std_noise;
    const gdouble m                 = data_i->m;
    const gdouble c1                = data_i->c1;
    const gdouble c2                = data_i->c2;
    complex double e_o              = e1 + I * e2;
    gdouble std_shape;
    complex double hat_g;
    gdouble var_tot, weight;
    gdouble theta, phi;

    nc_galaxy_shape_pop_prepare (pop, data_i->pop_data);
    std_shape = nc_galaxy_shape_pop_e_rms (pop, data_i->pop_data);
    var_tot   = std_shape * std_shape + std_noise * std_noise;
    weight    = 1.0 / var_tot;

    nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

    /* Re-express the celestial phi_C in the data frame (phi_E) and rotate the
     * data-frame observed ellipticity into the tangential/cross frame. The bias
     * (c1, c2) stays in the data frame and is averaged below without rotation.
     * See #NcWLEllipticityFrame. */
    phi = nc_wl_ellipticity_celestial_to_frame_angle (data_i->coord, phi);

    hat_g = e_o * cexp (-2.0 * I * phi);

    ncm_stats_vec_set (self->obs_stats, 0, creal (hat_g));
    ncm_stats_vec_set (self->obs_stats, 1, cimag (hat_g));
    ncm_stats_vec_set (self->obs_stats, 2, std_shape * std_shape);
    ncm_stats_vec_set (self->obs_stats, 3, m);
    ncm_stats_vec_set (self->obs_stats, 4, c1);
    ncm_stats_vec_set (self->obs_stats, 5, c2);

    ncm_stats_vec_update_weight (self->obs_stats, weight);
  }

  {
    const gdouble mean_gt          = ncm_stats_vec_get_mean (self->obs_stats, 0);
    const gdouble mean_gx          = ncm_stats_vec_get_mean (self->obs_stats, 1);
    const gdouble mean_sigma_true2 = ncm_stats_vec_get_mean (self->obs_stats, 2);
    const gdouble mean_m           = ncm_stats_vec_get_mean (self->obs_stats, 3);
    const gdouble mean_c1          = ncm_stats_vec_get_mean (self->obs_stats, 4);
    const gdouble mean_c2          = ncm_stats_vec_get_mean (self->obs_stats, 5);
    const gdouble R                = 1.0 - mean_sigma_true2;

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        /* epsilon convention: <epsilon> = g, no responsivity factor R */
        *gt      = (mean_gt - mean_c1) / (1.0 + mean_m);
        *gx      = (mean_gx - mean_c2) / (1.0 + mean_m);
        *sigma_t = ncm_stats_vec_get_sd (self->obs_stats, 0) / (1.0 + mean_m);
        *sigma_x = ncm_stats_vec_get_sd (self->obs_stats, 1) / (1.0 + mean_m);
        *rho     = ncm_stats_vec_get_cor (self->obs_stats, 0, 1);
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        /* distortion convention: <chi> = 2 R g, responsivity applies to both components */
        *gt      = (0.5 * mean_gt / R - mean_c1) / (1.0 + mean_m);
        *gx      = (0.5 * mean_gx / R - mean_c2) / (1.0 + mean_m);
        *sigma_t = 0.5 * ncm_stats_vec_get_sd (self->obs_stats, 0) / R / (1.0 + mean_m);
        *sigma_x = 0.5 * ncm_stats_vec_get_sd (self->obs_stats, 1) / R / (1.0 + mean_m);
        *rho     = ncm_stats_vec_get_cor (self->obs_stats, 0, 1);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }
  }
}

/**
 * nc_galaxy_shape_factor_eval_marginal: (virtual eval_marginal)
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop intrinsic-ellipticity model
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: real part of the reduced shear (bias applied)
 * @g_2: imaginary part of the reduced shear (bias applied)
 * @epsilon_obs_1: observed ellipticity component 1 (same frame as @g_1, @g_2)
 * @epsilon_obs_2: observed ellipticity component 2 (same frame as @g_1, @g_2)
 *
 * Evaluates the intrinsic-ellipticity marginal
 * P(epsilon_obs | g) = int d^2chi_I P_pop(chi_I) N_2(epsilon_obs - f_g(chi_I); std_noise^2)
 * using the subclass evaluation strategy. The shear and the observed
 * ellipticity must be expressed in the same frame; @data supplies the
 * per-galaxy noise (std_noise) and the population fragment (pop_data).
 *
 * Returns: the marginal probability density.
 */
gdouble
nc_galaxy_shape_factor_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf)->eval_marginal (gsf, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

/**
 * nc_galaxy_shape_factor_eval_ln_marginal: (virtual eval_ln_marginal)
 * @gsf: a #NcGalaxyShapeFactor
 * @pop: the #NcGalaxyShapePop intrinsic-ellipticity model
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: real part of the reduced shear (bias applied)
 * @g_2: imaginary part of the reduced shear (bias applied)
 * @epsilon_obs_1: observed ellipticity component 1 (same frame as @g_1, @g_2)
 * @epsilon_obs_2: observed ellipticity component 2 (same frame as @g_1, @g_2)
 *
 * Logarithmic counterpart of nc_galaxy_shape_factor_eval_marginal().
 *
 * Returns: the natural logarithm of the marginal probability density.
 */
gdouble
nc_galaxy_shape_factor_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf)->eval_ln_marginal (gsf, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

