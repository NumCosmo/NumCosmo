/***************************************************************************
 *            nc_cluster_mass_projection.c
 *
 *  Thu Jan 26 18:25:11 2017
 *  Copyright  2017  Cinthia Nunes de Lima and Henrique Lettieri Projection
 *  <cinthia.nlima@gmail.com>, <henrique.lettieri@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Cinthia Nunes de Lima and Henrique Lettieri Projection 2017 <cinthia.nlima@gmail.com>
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

/**
 * NcClusterMassProjection:
 *
 * Cluster mass distribution model with projection effects.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/cluster/nc_cluster_mass_projection.h"
#include "ncm/integration/ncm_integrate.h"
#include "ncm/spline/ncm_spline2d_bicubic.h"
#include "ncm/spline/ncm_spline2d.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */


#define _NC_CLUSTER_MASS_PROJECTION_DEFAULT_INT_KEY 6

typedef struct _NcClusterMassProjectionPrivate
{
  gdouble M0;
  gdouble z0;
  gdouble lnM0;
  gdouble ln1pz0;
  gdouble lnR_max;
  gdouble lnR_min;
  gboolean enable_rejection;
} NcClusterMassProjectionPrivate;


struct _NcClusterMassProjection
{
  NcClusterMass parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterMassProjection, nc_cluster_mass_projection, NC_TYPE_CLUSTER_MASS)

#define VECTOR   (NCM_MODEL (projection))
#define MU_P0    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_MU_P0))
#define MU_P1    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_MU_P1))
#define MU_P2    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_MU_P2))
#define SIGMA_P0 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_SIGMA_P0))
#define SIGMA_P1 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_SIGMA_P1))
#define SIGMA_P2 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_SIGMA_P2))
#define LNM_C0   (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_LNM_C0))
#define LNM_CZ   (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_LNM_CZ))
#define A_C0     (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_A_C0))
#define A_CZ     (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_A_CZ))
#define CUT      (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_CUT))


enum
{
  PROP_0,
  PROP_M0,
  PROP_Z0,
  PROP_LNRICHNESS_MIN,
  PROP_LNRICHNESS_MAX,
  PROP_ENABLE_REJECTION,
  PROP_SIZE,
};

static void
nc_cluster_mass_projection_init (NcClusterMassProjection *projection)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  self->M0               = 0.0;
  self->z0               = 0.0;
  self->lnM0             = 0.0;
  self->ln1pz0           = 0.0;
  self->lnR_min          = GSL_NEGINF;
  self->lnR_max          = GSL_POSINF;
  self->enable_rejection = TRUE;
}

static void
_nc_cluster_mass_projection_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (object);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  g_return_if_fail (NC_IS_CLUSTER_MASS_PROJECTION (object));

  switch (prop_id)
  {
    case PROP_M0:
      self->M0   = g_value_get_double (value);
      self->lnM0 = log (self->M0);
      break;
    case PROP_Z0:
      self->z0     = g_value_get_double (value);
      self->ln1pz0 = log1p (self->z0);
      break;
    case PROP_LNRICHNESS_MIN:
      self->lnR_min = g_value_get_double (value);
      g_assert (self->lnR_min < self->lnR_max);
      break;
    case PROP_LNRICHNESS_MAX:
      self->lnR_max = g_value_get_double (value);
      g_assert (self->lnR_min < self->lnR_max);
      break;
    case PROP_ENABLE_REJECTION:
      nc_cluster_mass_projection_set_enable_rejection (projection, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_cluster_mass_projection_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (object);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  g_return_if_fail (NC_IS_CLUSTER_MASS_PROJECTION (object));

  switch (prop_id)
  {
    case PROP_M0:
      g_value_set_double (value, self->M0);
      break;
    case PROP_Z0:
      g_value_set_double (value, self->z0);
      break;
    case PROP_LNRICHNESS_MIN:
      g_value_set_double (value, self->lnR_min);
      break;
    case PROP_LNRICHNESS_MAX:
      g_value_set_double (value, self->lnR_max);
      break;
    case PROP_ENABLE_REJECTION:
      g_value_set_boolean (value, self->enable_rejection);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_cluster_mass_projection_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_projection_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_projection_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
static gdouble _nc_cluster_mass_projection_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_projection_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
static gboolean _nc_cluster_mass_projection_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_projection_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_projection_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_projection_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);
static gdouble _nc_cluster_mass_projection_volume (NcClusterMass *clusterm);
static void _nc_cluster_mass_projection_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, NcmVector *res);

static void
nc_cluster_mass_projection_class_init (NcClusterMassProjectionClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcClusterMassClass *parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_cluster_mass_projection_set_property;
  model_class->get_property = &_nc_cluster_mass_projection_get_property;
  object_class->finalize    = &_nc_cluster_mass_projection_finalize;

  ncm_model_class_set_name_nick (model_class, "Projection Ln-normal richness distribution with projection function", "Projection");
  ncm_model_class_add_params (model_class, NC_CLUSTER_MASS_PROJECTION_SPARAM_LEN, 0, PROP_SIZE);

  /**
   * NcClusterMassProjection:M0:
   *
   * Pivot (reference) mass used to make observed and model masses
   * dimensionless. The property default and allowed range are declared in the
   * property registration below (see g_param_spec_double).
   */
  g_object_class_install_property (object_class,
                                   PROP_M0,
                                   g_param_spec_double ("M0",
                                                        NULL,
                                                        "Pivot mass",
                                                        11.0 * M_LN10, G_MAXDOUBLE, 3.0e14 / 0.71,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

/*
 * NcClusterMassProjection:Z0:
 *
 * Pivot redshift used to center redshift-dependent scaling relations. The
 * concrete default and bounds are set in the property registration below.
 */
  g_object_class_install_property (object_class,
                                   PROP_Z0,
                                   g_param_spec_double ("z0",
                                                        NULL,
                                                        "Pivot redshift",
                                                        0.0, G_MAXDOUBLE, 0.6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));



  /**
   * NcClusterMassProjection:lnRichness_min:
   *
   * Minimum observed richness (lower bound for the projection function). The
   * default and allowed range are provided in the property declaration below.
   */
  g_object_class_install_property (object_class,
                                   PROP_LNRICHNESS_MIN,
                                   g_param_spec_double ("lnRichness-min",
                                                        NULL,
                                                        "Minimum LnRichness",
                                                        0.0, G_MAXDOUBLE, M_LN10 * 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassProjection:lnRichness_max:
   *
   * Maximum observed richness (upper bound for the projection function). The
   * default and allowed range are provided in the property declaration below.
   */
  g_object_class_install_property (object_class,
                                   PROP_LNRICHNESS_MAX,
                                   g_param_spec_double ("lnRichness-max",
                                                        NULL,
                                                        "Maximum LnRichness",
                                                        0.0, G_MAXDOUBLE,  M_LN10 * 2.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassProjection:enable_rejection:
   *
   * When TRUE, generated observed richness values below the richness cut are
   * rejected; when FALSE sampling uses the truncated tail generator so values
   * below the cut are produced according to the tail distribution.
   */
  g_object_class_install_property (object_class,
                                   PROP_ENABLE_REJECTION,
                                   g_param_spec_boolean ("enable-rejection",
                                                         NULL,
                                                         "Whether rejects the sampled objects below the CUT",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  /**
   * NcClusterMassProjection:MU_P0:
   *
   * Intercept term of the mean richness model (bias in the mean). The allowed
   * range and default are set in the parameter registration below.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_MU_P0, "mu_p0", "mup0",
                              0.0,  6.0, 1.0e-1,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:MU_P1:
   *
   * Mass dependence (slope) of the mean richness model. Parameter bounds and
   * default are declared in the registration call below.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_MU_P1, "mu_p1", "mup1",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:MU_P2:
   *
   * Redshift dependence (slope) of the mean richness model. Parameter bounds
   * and defaults are provided when the parameter is registered.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_MU_P2, "mu_p2", "mup2",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:sigma_P0:
   *
   * Intercept term of the richness scatter model (standard deviation). The
   * registration call below constrains the allowed interval and provides a
   * default value.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_SIGMA_P0, "\\sigma_p0", "sigmap0",
                              1.0e-4, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:sigma_P1:
   *
   * Mass dependence (slope) of the richness scatter model. Bounds and default
   * are set in the parameter registration below.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_SIGMA_P1, "\\sigma_p1", "sigmap1",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P1,
                              NCM_PARAM_TYPE_FIXED);


/**
 * NcClusterMassProjection:sigma_P2:
 *
 * Redshift dependence (slope) of the richness scatter model. Bounds and
 * default are set in the registration below.
 */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_SIGMA_P2, "\\sigma_p2", "sigmap2",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P2,
                              NCM_PARAM_TYPE_FIXED);




/**
 * NcClusterMassProjection:CUT:
 *
 * Cut in richness.
 *
 */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_CUT, "CUT", "cut",
                              0.0,  1.0e16, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_CUT,
                              NCM_PARAM_TYPE_FIXED);


  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->P               = &_nc_cluster_mass_projection_p;
  parent_class->intP            = &_nc_cluster_mass_projection_intp;
  parent_class->intP_bin        = &_nc_cluster_mass_projection_intp_bin;
  parent_class->resample        = &_nc_cluster_mass_projection_resample;
  parent_class->P_limits        = &_nc_cluster_mass_projection_p_limits;
  parent_class->P_bin_limits    = &_nc_cluster_mass_projection_p_bin_limits;
  parent_class->N_limits        = &_nc_cluster_mass_projection_n_limits;
  parent_class->volume          = &_nc_cluster_mass_projection_volume;
  parent_class->P_vec_z_lnMobs  = &_nc_cluster_mass_projection_p_vec_z_lnMobs;
  parent_class->_obs_len        = 1;
  parent_class->_obs_params_len = 0;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);
}

static void
_nc_cluster_mass_projection_lnR_sigma (NcClusterMass *clusterm, const gdouble lnM, const gdouble z, gdouble *lnR, gdouble *sigma)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);
  const gdouble DlnM                         = lnM - self->lnM0;
  const gdouble Dln1pz                       = log1p (z) - self->ln1pz0;

  lnR[0]   = MU_P0    + MU_P1    * DlnM + MU_P2    * Dln1pz;
  sigma[0] = SIGMA_P0 + SIGMA_P1 * DlnM + SIGMA_P2 * Dln1pz;
}

void
nc_cluster_mass_projection_set_completeness (NcClusterMassProjection *projection, NcmSpline2dBicubic *completeness)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  ncm_spline2d_clear (&self->completeness);
  self->completeness = NCM_SPLINE2D (completeness);
  ncm_spline2d_prepare (self->completeness);
}

/**
 * nc_cluster_mass_projection_peek_completeness:
 * @projection: a #NcClusterMassProjection
 *
 * Get the spline for the completeness as function of
 * $\ln(M)$ and $z$.
 *
 * Returns: (transfer none): the spline for the cluster completeness.
 */
NcmSpline2d *
nc_cluster_mass_projection_peek_completeness (NcClusterMassProjection *projection)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  return self->completeness;
}

gdouble
nc_cluster_mass_projection_completeness (NcClusterMassProjection *projection, gdouble lnM, gdouble z)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  if (self->completeness == NULL)
    return 1.0;
  else
    return ncm_spline2d_eval (self->completeness, lnM, z);
}

void
nc_cluster_mass_projection_set_ipurity (NcClusterMassProjection *projection, NcmSpline2dBicubic *ipurity)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  self->ipurity = NCM_SPLINE2D (ipurity);
  ncm_spline2d_prepare (self->ipurity);
}

/**
 * nc_cluster_mass_projection_peek_ipurity:
 * @projection: a #NcClusterMassProjection
 *
 * Get the spline for the inverse of purity as function of
 * $\ln(M_obs)$ and $z$.
 *
 * Returns: (transfer none): the spline for the cluster inverse purity.
 */

NcmSpline2d *
nc_cluster_mass_projection_peek_ipurity (NcClusterMassProjection *projection)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  return self->ipurity;
}

gdouble
nc_cluster_mass_projection_ipurity (NcClusterMassProjection *projection, gdouble lnM_obs, gdouble z)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  if (self->ipurity == NULL)
    return 1.0;
  else
    return ncm_spline2d_eval (self->ipurity, lnM_obs, z);
}

static gdouble
_nc_cluster_mass_projection_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{
  NcClusterMassProjection *projection = NC_CLUSTER_MASS_PROJECTION (clusterm);
  gdouble lnR_true, sigma;

  _nc_cluster_mass_projection_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  {
    const gdouble x            = (lnM_obs[0] - lnR_true) / sigma;
    const gdouble completeness = nc_cluster_mass_projection_completeness (projection, lnM, z);
    const gdouble ipurity      = nc_cluster_mass_projection_ipurity (projection, lnM_obs[0], z);

    if (lnM_obs[0] < CUT)
      return 0.0;
    else
      return fabs (2.0 / (ncm_c_sqrt_2pi () * sigma) * exp (-0.5 * x * x) / erfc ((CUT - lnR_true) / (M_SQRT2 * sigma)) * completeness * ipurity);
  }
}

typedef struct _NcClusterMassProjectionInt
{
  NcHICosmo *cosmo;
  NcClusterMass *clusterm;
  gdouble lnM_obs;
  gdouble lnM;
  gdouble z;
  const gdouble *lnM_obs_params;
} NcClusterMassProjectionInt;

static gdouble
_nc_cluster_mass_projection_integrand (gdouble lnM_obs, gpointer userdata)
{
  NcClusterMassProjectionInt *obs_data = (NcClusterMassProjectionInt *) userdata;
  NcClusterMassProjection *projection   = NC_CLUSTER_MASS_PROJECTION (obs_data->clusterm);
  gdouble lnR_true, sigma;

  _nc_cluster_mass_projection_lnR_sigma (obs_data->clusterm, obs_data->lnM, obs_data->z, &lnR_true, &sigma);
  {
    const gdouble x       = (lnM_obs - lnR_true) / sigma;
    const gdouble ipurity = nc_cluster_mass_projection_ipurity (projection, lnM_obs, obs_data->z);

    if (lnM_obs < CUT)
      return 0.0;
    else
      return 2.0 / (ncm_c_sqrt_2pi () * sigma) * exp (-0.5 * x * x) / erfc ((CUT - lnR_true) / (M_SQRT2 * sigma)) * ipurity;
  }
}

static gdouble
_nc_cluster_mass_projection_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);
  gsl_integration_workspace **w              = ncm_integral_get_workspace ();
  NcClusterMassProjectionInt obs_data;
  gdouble intp, err, completeness;
  gsl_function F;

  obs_data.clusterm = clusterm;
  obs_data.cosmo    = cosmo;
  obs_data.lnM      = lnM;
  obs_data.z        = z;

  F.function = &_nc_cluster_mass_projection_integrand;
  F.params   = &obs_data;

  gsl_integration_qag (&F, CUT, self->lnR_max, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_MASS_PROJECTION_DEFAULT_INT_KEY, *w, &intp, &err);
  ncm_memory_pool_return (w);
  completeness = nc_cluster_mass_projection_completeness (projection, lnM, z);

  return fabs (intp * completeness);
}

static gdouble
_nc_cluster_mass_projection_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params)
{
  NcClusterMassProjection *projection = NC_CLUSTER_MASS_PROJECTION (clusterm);

  if ((lnM_obs_lower[0] < CUT) && (lnM_obs_upper[0] < CUT))
  {
    return 0.0;
  }
  else
  {
    NcClusterMassProjectionInt obs_data;
    gdouble intp_bin, err, completeness;
    gsl_function F;
    gsl_integration_workspace **w = ncm_integral_get_workspace ();

    obs_data.clusterm       = clusterm;
    obs_data.cosmo          = cosmo;
    obs_data.lnM            = lnM;
    obs_data.lnM_obs_params = lnM_obs_params;
    obs_data.z              = z;
    completeness            = nc_cluster_mass_projection_completeness (projection, lnM, z);

    F.function = &_nc_cluster_mass_projection_integrand;
    F.params   = &obs_data;

    if ((lnM_obs_lower[0] < CUT) && (lnM_obs_upper[0] >= CUT))
      gsl_integration_qag (&F, CUT, lnM_obs_upper[0], 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_MASS_PROJECTION_DEFAULT_INT_KEY, *w, &intp_bin, &err);
    else
      gsl_integration_qag (&F, lnM_obs_lower[0], lnM_obs_upper[0], 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_MASS_PROJECTION_DEFAULT_INT_KEY, *w, &intp_bin, &err);

    ncm_memory_pool_return (w);

    return fabs (intp_bin * completeness);
  }

  g_assert_not_reached ();

  return 0.0;
}

static gboolean
_nc_cluster_mass_projection_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);
  gdouble lnR_true, sigma;

  _nc_cluster_mass_projection_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  ncm_rng_lock (rng);

  if (self->enable_rejection)
  {
    lnM_obs[0] = ncm_rng_gaussian_gen (rng, lnR_true, sigma);
  }
  else
  {
    lnM_obs[0]  = ncm_rng_gaussian_tail_gen (rng, CUT - lnR_true, sigma);
    lnM_obs[0] += lnR_true;
  }

  ncm_rng_unlock (rng);

  return (lnM_obs[0] <= self->lnR_max) && (lnM_obs[0] >= CUT);
}

/**
 * ncluster_mass_projection_set_lnM_limits:
 * @projection: a #NcClusterMassProjection
 * @lnM_limits: a #NcmVector of length 2
 *
 * Set the limits for the cluster mass function.
 *
 */
void
nc_cluster_mass_projection_set_lnM_limits (NcClusterMassProjection *projection, NcmVector *lnM_limits)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  g_assert_cmpuint (ncm_vector_len (lnM_limits), ==, 2);

  ncm_vector_clear (&self->lnM_limits);
  self->lnM_limits =  ncm_vector_ref (lnM_limits);
}

static void
_nc_cluster_mass_projection_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);
  gdouble lnMl, lnMu;

  if (self->lnM_limits == NULL)
  {
    lnMl =  M_LN10 * 12.0;
    lnMu =  M_LN10 * 16.0;
  }
  else
  {
    lnMl =  ncm_vector_get (self->lnM_limits, 0);
    lnMu =  ncm_vector_get (self->lnM_limits, 1);
  }

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static void
_nc_cluster_mass_projection_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);
  gdouble lnMl, lnMu;

  if (self->lnM_limits == NULL)
  {
    lnMl =  M_LN10 * 12.0;
    lnMu =  M_LN10 * 16.0;
  }
  else
  {
    lnMl =  ncm_vector_get (self->lnM_limits, 0);
    lnMu =  ncm_vector_get (self->lnM_limits, 1);
  }

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static void
_nc_cluster_mass_projection_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);
  gdouble lnMl, lnMu;

  if (self->lnM_limits == NULL)
  {
    lnMl =  M_LN10 * 12.0;
    lnMu =  M_LN10 * 16.0;
  }
  else
  {
    lnMl =  ncm_vector_get (self->lnM_limits, 0);
    lnMu =  ncm_vector_get (self->lnM_limits, 1);
  }

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static gdouble
_nc_cluster_mass_projection_volume (NcClusterMass *clusterm)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  return (self->lnR_max - CUT);
}

static void
_nc_cluster_mass_projection_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, NcmVector *res)
{
  NcClusterMassProjection *projection          = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  const gdouble *lnM_obs_ptr = ncm_matrix_const_data (lnM_obs);
  const gdouble *z_ptr       = ncm_vector_const_data (z);
  const guint tda            = ncm_matrix_tda (lnM_obs);
  const guint sz             = ncm_vector_stride (z);
  const guint len            = ncm_vector_len (z);
  const gdouble DlnM         = lnM - self->lnM0;
  const gdouble sqrt_2pi     = ncm_c_sqrt_2pi ();
  const gdouble lnR_pre      = MU_P0    + MU_P1    * DlnM;
  const gdouble sigma_pre    = SIGMA_P0 + SIGMA_P1 * DlnM;
  const gdouble mu_p2        = MU_P2;
  const gdouble sigma_p2     = SIGMA_P2;
  const gdouble cut          = CUT;
  gdouble *res_ptr           = ncm_vector_ptr (res, 0);
  guint i;
  gdouble completeness, ipurity;

  g_assert_cmpuint (ncm_vector_stride (res), ==, 1);

  if ((tda == 1) && (sz == 1))
  {
    for (i = 0; i < len; i++)
    {
      const gdouble Dln1pz = log1p (z_ptr[i]) - self->ln1pz0;
      const gdouble lnR    = lnR_pre + mu_p2 * Dln1pz;
      const gdouble sigma  = sigma_pre + sigma_p2 * Dln1pz;
      const gdouble x      = (lnM_obs_ptr[i] - lnR) / sigma;

      completeness = nc_cluster_mass_projection_completeness (projection, lnM, z_ptr[i]);
      ipurity      = nc_cluster_mass_projection_ipurity (projection, lnM_obs_ptr[i], z_ptr[i]);

      if (lnM_obs_ptr[i] < cut)
        res_ptr[i] = 0.0;
      else
        res_ptr[i] = 2.0 / (sqrt_2pi * sigma) * exp (-0.5 * x * x) / erfc ((CUT - lnR) / (M_SQRT2 * sigma)) * completeness * ipurity;
    }
  }
  else
  {
    for (i = 0; i < len; i++)
    {
      const gdouble Dln1pz = log1p (z_ptr[i * sz]) - self->ln1pz0;
      const gdouble lnR    = lnR_pre + mu_p2 * Dln1pz;
      const gdouble sigma  = sigma_pre + sigma_p2 * Dln1pz;
      const gdouble x      = (lnM_obs_ptr[i * tda] - lnR) / sigma;
      gdouble completeness = nc_cluster_mass_projection_completeness (projection, lnM, z_ptr[i * sz]);
      gdouble ipurity      = nc_cluster_mass_projection_ipurity (projection, lnM_obs_ptr[i * tda], z_ptr[i *  sz]);

      if (lnM_obs_ptr[i * tda] < cut)
        res_ptr[i] = 0.0;
      else
        res_ptr[i] = 2.0 / (sqrt_2pi * sigma) * exp (-0.5 * x * x) / erfc ((CUT - lnR) / (M_SQRT2 * sigma)) * completeness * ipurity;
    }
  }
}

/**
 * nc_cluster_mass_projection_get_mean_richness:
 * @projection: a #NcClusterMassProjection
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the mean of the richness distribution.
 *
 */
gdouble
nc_cluster_mass_projection_get_mean_richness (NcClusterMassProjection *projection, gdouble lnM, gdouble z)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);
  const gdouble DlnM                         = lnM - self->lnM0;
  const gdouble Dln1pz                       = log1p (z) - self->ln1pz0;

  return MU_P0 + MU_P1 * DlnM + MU_P2 * Dln1pz;
}

/**
 * nc_cluster_mass_projection_get_std_richness:
 * @projection: a #NcClusterMassProjection
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the standard deviation of the richness distribution.
 *
 */
gdouble
nc_cluster_mass_projection_get_std_richness (NcClusterMassProjection *projection, gdouble lnM, gdouble z)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);
  const gdouble DlnM                         = lnM - self->lnM0;
  const gdouble Dln1pz                       = log1p (z) - self->ln1pz0;

  return SIGMA_P0 + SIGMA_P1 * DlnM + SIGMA_P2 * Dln1pz;
}

/**
 * nc_cluster_mass_projection_get_cut:
 * @projection: a #NcClusterMassProjection
 *
 * Computes the cut in richness.
 *
 * Returns: the cut in richness.
 */
gdouble
nc_cluster_mass_projection_get_cut (NcClusterMassProjection *projection)
{
  return CUT;
}

/**
 * nc_cluster_mass_projection_get_mean:
 * @projection: a #NcClusterMassProjection
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the mean of the richness distribution with the cut correction.
 *
 */
gdouble
nc_cluster_mass_projection_get_mean (NcClusterMassProjection *projection, gdouble lnM, gdouble z)
{
  const gdouble lnR_mean        = nc_cluster_mass_projection_get_mean_richness (projection, lnM, z);
  const gdouble lnR_sigma       = nc_cluster_mass_projection_get_std_richness  (projection, lnM, z);
  const gdouble A               = (CUT - lnR_mean) / lnR_sigma;
  const gdouble B               = (1.0 / (ncm_c_sqrt_2pi ())) * exp (-0.5 * (A  * A));
  const gdouble C               = 1.0 - 0.5 * (1.0 + erf (A / M_SQRT2));
  const gdouble mean_correction = (lnR_sigma * B / C);

  return lnR_mean + mean_correction;
}

/**
 * nc_cluster_mass_projection_get_std:
 * @projection: a #NcClusterMassProjection
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the standard deviation of the richness distribution with the cut correction.
 *
 */
gdouble
nc_cluster_mass_projection_get_std (NcClusterMassProjection *projection, gdouble lnM, gdouble z)
{
  const gdouble lnR_mean       = nc_cluster_mass_projection_get_mean_richness (projection, lnM, z);
  const gdouble lnR_sigma      = nc_cluster_mass_projection_get_std_richness  (projection, lnM, z);
  const gdouble A              = (CUT - lnR_mean) / lnR_sigma;
  const gdouble B              = (1.0 / (ncm_c_sqrt_2pi ())) * exp (-0.5 * (A  * A));
  const gdouble C              = 1.0 - 0.5 * (1.0 + erf (A / M_SQRT2));
  const gdouble std_correction = pow (1.0 + (A * B / C) - (B / C) * (B / C), 0.5);

  return lnR_sigma * std_correction;
}

/**
 * nc_cluster_mass_projection_set_enable_rejection:
 * @projection: a #NcClusterMassProjection
 * @on: a #NcClusterMassProjection
 *
 * Set the enable_rejection property.
 *
 */

void
nc_cluster_mass_projection_set_enable_rejection (NcClusterMassProjection *projection, gboolean on)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  self->enable_rejection = on;
}

/**
 * nc_cluster_mass_projection_get_enable_rejection:
 * @projection: a #NcClusterMassProjection
 *
 * Get if the enable_rejection property is on.
 *
 */
gboolean
nc_cluster_mass_projection_get_enable_rejection (NcClusterMassProjection *projection)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (projection);

  return self->enable_rejection;
}

