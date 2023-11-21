/***************************************************************************
 *            nc_cluster_mass_lnnormal.c
 *
 *  Fri June 22 17:48:55 2012
 *  Copyright  2012  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_mass_lnnormal
 * @title: NcClusterMassLnnormal
 * @short_description: Cluster mass ln-normal distribution.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_lnnormal.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcClusterMassLnnormal, nc_cluster_mass_lnnormal, NC_TYPE_CLUSTER_MASS)

#define VECTOR (NCM_MODEL (mlnn)->params)
#define BIAS   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_LNNORMAL_BIAS))
#define SIGMA  (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_LNNORMAL_SIGMA))

enum
{
  PROP_0,
  PROP_LNMOBS_MIN,
  PROP_LNMOBS_MAX,
  PROP_SIZE,
};

static void
nc_cluster_mass_lnnormal_init (NcClusterMassLnnormal *mlnm)
{
  mlnm->lnMobs_min = GSL_NEGINF;
  mlnm->lnMobs_max = GSL_POSINF;
}

static void
_nc_cluster_mass_lnnormal_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassLnnormal *mlnm = NC_CLUSTER_MASS_LNNORMAL (object);

  g_return_if_fail (NC_IS_CLUSTER_MASS_LNNORMAL (object));

  switch (prop_id)
  {
    case PROP_LNMOBS_MIN:
      mlnm->lnMobs_min = g_value_get_double (value);
      g_assert (mlnm->lnMobs_min < mlnm->lnMobs_max);
      break;
    case PROP_LNMOBS_MAX:
      mlnm->lnMobs_max = g_value_get_double (value);
      g_assert (mlnm->lnMobs_min < mlnm->lnMobs_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_lnnormal_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassLnnormal *mlnm = NC_CLUSTER_MASS_LNNORMAL (object);

  g_return_if_fail (NC_IS_CLUSTER_MASS_LNNORMAL (object));

  switch (prop_id)
  {
    case PROP_LNMOBS_MIN:
      g_value_set_double (value, mlnm->lnMobs_min);
      break;
    case PROP_LNMOBS_MAX:
      g_value_set_double (value, mlnm->lnMobs_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_lnnormal_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_lnnormal_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_lnnormal_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
static gdouble _nc_cluster_mass_lnnormal_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_lnnormal_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
static gboolean _nc_cluster_mass_lnnormal_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_lnnormal_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_lnnormal_p_limits_bin (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_lnnormal_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);

static void
nc_cluster_mass_lnnormal_class_init (NcClusterMassLnnormalClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcClusterMassClass *parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_cluster_mass_lnnormal_finalize;

  model_class->set_property = &_nc_cluster_mass_lnnormal_set_property;
  model_class->get_property = &_nc_cluster_mass_lnnormal_get_property;

  ncm_model_class_set_name_nick (model_class, "Ln-normal distribution", "Ln_Normal");
  ncm_model_class_add_params (model_class, 2, 0, PROP_SIZE);

  /**
   * NcClusterMassLnnormal:lnMobs_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_LNMOBS_MIN,
                                   g_param_spec_double ("lnMobs-min",
                                                        NULL,
                                                        "Minimum LnMobs",
                                                        11.0 * M_LN10, G_MAXDOUBLE, log (5.0) + 13.0 * M_LN10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassLnnormal:lnMobs_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_LNMOBS_MAX,
                                   g_param_spec_double ("lnMobs-max",
                                                        NULL,
                                                        "Maximum LnMobs",
                                                        11.0 * M_LN10, G_MAXDOUBLE, 16.0 * M_LN10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassLnnormal:bias:
   *
   * Distribution's bias.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNNORMAL_BIAS, "bias", "bias",
                              0.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNNORMAL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNNORMAL_DEFAULT_BIAS,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnnormal:sigma:
   *
   * Distribution's standard deviation, $\sigma \in [10^{-4}, 10]$.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNNORMAL_SIGMA, "sigma", "sigma",
                              1.0e-4,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNNORMAL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNNORMAL_DEFAULT_SIGMA,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->P              = &_nc_cluster_mass_lnnormal_p;
  parent_class->intP           = &_nc_cluster_mass_lnnormal_intp;
  parent_class->intP_bin       = &_nc_cluster_mass_lnnormal_intp_bin;
  parent_class->resample       = &_nc_cluster_mass_lnnormal_resample;
  parent_class->P_limits       = &_nc_cluster_mass_lnnormal_p_limits;
  parent_class->P_bin_limits   = &_nc_cluster_mass_lnnormal_p_limits_bin;
  parent_class->N_limits       = &_nc_cluster_mass_lnnormal_n_limits;
  parent_class->obs_len        = 1;
  parent_class->obs_params_len = 0;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);
}

static gdouble
_nc_cluster_mass_lnnormal_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{
  NcClusterMassLnnormal *mlnn = NC_CLUSTER_MASS_LNNORMAL (clusterm);
  const gdouble lnMobs        = lnM_obs[0];
  const gdouble sqrt2_sigma   = M_SQRT2 * SIGMA;
  const gdouble x             = (lnMobs - lnM - BIAS) / sqrt2_sigma;

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  return M_2_SQRTPI / (2.0 * M_SQRT2) * exp (-x * x) / (SIGMA);
}

static gdouble
_nc_cluster_mass_lnnormal_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcClusterMassLnnormal *mlnn = NC_CLUSTER_MASS_LNNORMAL (clusterm);
  const gdouble sqrt2_sigma   = M_SQRT2 * SIGMA;
  const gdouble x_min         = (lnM - mlnn->lnMobs_min) / sqrt2_sigma;
  const gdouble x_max         = (lnM - mlnn->lnMobs_max) / sqrt2_sigma;

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  if (x_max > 4.0)
    return -(erfc (x_min) - erfc (x_max)) / 2.0;
  else
    return (erf (x_min) - erf (x_max)) / 2.0;
}

static gdouble
_nc_cluster_mass_lnnormal_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params)
{
  NcClusterMassLnnormal *mlnn = NC_CLUSTER_MASS_LNNORMAL (clusterm);
  const gdouble sqrt2_sigma   = M_SQRT2 * SIGMA;
  const gdouble x_min         = (lnM - lnM_obs_lower[0]) / sqrt2_sigma;
  const gdouble x_max         = (lnM - lnM_obs_upper[0]) / sqrt2_sigma;

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  if (x_max > 4.0)
    return -(erfc (x_min) - erfc (x_max)) / 2.0;
  else
    return (erf (x_min) - erf (x_max)) / 2.0;
}

static gboolean
_nc_cluster_mass_lnnormal_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassLnnormal *mlnn = NC_CLUSTER_MASS_LNNORMAL (clusterm);

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  ncm_rng_lock (rng);
  lnM_obs[0] = lnM + BIAS + gsl_ran_gaussian (rng->r, SIGMA);
  ncm_rng_unlock (rng);

  return (lnM_obs[0] <= mlnn->lnMobs_max) && (lnM_obs[0] >= mlnn->lnMobs_min);
}

static void
_nc_cluster_mass_lnnormal_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassLnnormal *mlnn = NC_CLUSTER_MASS_LNNORMAL (clusterm);
  const gdouble mean          = lnM_obs[0] - BIAS;
  const gdouble lnMl          = mean - 7.0 * SIGMA;
  const gdouble lnMu          = mean + 7.0 * SIGMA;

  NCM_UNUSED (cosmo);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static void
_nc_cluster_mass_lnnormal_p_limits_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassLnnormal *mlnn = NC_CLUSTER_MASS_LNNORMAL (clusterm);
  const gdouble lnMl          = lnM_obs_lower[0] - 7.0 * SIGMA;
  const gdouble lnMu          = lnM_obs_upper[0] + 7.0 * SIGMA;

  NCM_UNUSED (cosmo);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;
}

static void
_nc_cluster_mass_lnnormal_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassLnnormal *mlnn = NC_CLUSTER_MASS_LNNORMAL (clusterm);
  const gdouble lnMl          = mlnn->lnMobs_min - 7.0 * SIGMA;
  const gdouble lnMu          = mlnn->lnMobs_max + 7.0 * SIGMA;

  NCM_UNUSED (cosmo);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

