/***************************************************************************
 *            nc_cluster_mass_vanderlinde.c
 *
 *  Tue July 3 15:21:05 2012
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
 * SECTION:nc_cluster_mass_vanderlinde
 * @title: NcClusterMassVanderlinde
 * @short_description: Sunyaev-Zel'dovich cluster mass distribution.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_vanderlinde.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcClusterMassVanderlinde, nc_cluster_mass_vanderlinde, NC_TYPE_CLUSTER_MASS);

#define VECTOR (NCM_MODEL (msz)->params)
#define A_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_VANDERLINDE_A_SZ))
#define B_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_VANDERLINDE_B_SZ))
#define C_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_VANDERLINDE_C_SZ))
#define D_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_VANDERLINDE_D_SZ))

enum
{
  PROP_0,
  PROP_SIGNIFICANCE_OBS_MIN,
  PROP_SIGNIFICANCE_OBS_MAX,
  PROP_Z0,
  PROP_M0,
  PROP_SIZE,
};

typedef struct _integrand_data
{
  NcClusterMassVanderlinde *msz;
  const gdouble *xi_params;
  gdouble z;
  gdouble lnM;
  const gdouble *xi;
  /* new: optimizing */
  gdouble lnA;
  gdouble lnM0;
  gdouble ln1pz_1pz0; /* ln ((1 + z)/(1 + z0))  */
  gdouble mu;
  gdouble erf_mu;
  gdouble D2_2; /* 2.0 * D * D */
} integrand_data;

static void
nc_cluster_mass_vanderlinde_init (NcClusterMassVanderlinde *msz)
{
  msz->signif_obs_min = 0.0;
  msz->signif_obs_max = 0.0;
  msz->z0             = 0.0;
  msz->M0             = 0.0;
}

static void
_nc_cluster_mass_vanderlinde_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (object);

  g_return_if_fail (NC_IS_CLUSTER_MASS_VANDERLINDE (object));

  switch (prop_id)
  {
    case PROP_SIGNIFICANCE_OBS_MIN:
      msz->signif_obs_min = g_value_get_double (value);
      break;
    case PROP_SIGNIFICANCE_OBS_MAX:
      msz->signif_obs_max = g_value_get_double (value);
      break;
    case PROP_Z0:
      msz->z0 = g_value_get_double (value);
      break;
    case PROP_M0:
      msz->M0 = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_vanderlinde_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (object);

  g_return_if_fail (NC_IS_CLUSTER_MASS_VANDERLINDE (object));

  switch (prop_id)
  {
    case PROP_SIGNIFICANCE_OBS_MIN:
      g_value_set_double (value, msz->signif_obs_min);
      break;
    case PROP_SIGNIFICANCE_OBS_MAX:
      g_value_set_double (value, msz->signif_obs_max);
      break;
    case PROP_Z0:
      g_value_set_double (value, msz->z0);
      break;
    case PROP_M0:
      g_value_set_double (value, msz->M0);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_vanderlinde_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_vanderlinde_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_vanderlinde_significance_m_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *xi, const gdouble *xi_params);
static gdouble _nc_cluster_mass_vanderlinde_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z);
static void _nc_cluster_mass_vanderlinde_p_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *xi, const gdouble *xi_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_vanderlinde_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);
static gboolean _nc_cluster_mass_vanderlinde_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *xi, const gdouble *xi_params, NcmRNG *rng);

static void
nc_cluster_mass_vanderlinde_class_init (NcClusterMassVanderlindeClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcClusterMassClass *parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  object_class->finalize = _nc_cluster_mass_vanderlinde_finalize;

  model_class->set_property = &_nc_cluster_mass_vanderlinde_set_property;
  model_class->get_property = &_nc_cluster_mass_vanderlinde_get_property;

  ncm_model_class_set_name_nick (model_class, "Vanderlinde et al. 2010 - SZ", "Vanderlinde_SZ");
  ncm_model_class_add_params (model_class, 4, 0, PROP_SIZE);

  /**
   * NcClusterMassVanderlinde:signif_obs_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGNIFICANCE_OBS_MIN,
                                   g_param_spec_double ("signif-obs-min",
                                                        NULL,
                                                        "Minimum observational significance",
                                                        2.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassVanderlinde:signif_obs_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGNIFICANCE_OBS_MAX,
                                   g_param_spec_double ("signif-obs-max",
                                                        NULL,
                                                        "Maximum observational significance",
                                                        2.0, G_MAXDOUBLE, 40.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassVanderlinde:z0:
   *
   * Reference redshift in the SZ signal-mass scaling relation.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_Z0,
                                   g_param_spec_double ("z0",
                                                        NULL,
                                                        "Reference redshift",
                                                        0.0, G_MAXDOUBLE, 0.6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassVanderlinde:M0:
   *
   * Reference mass (in h^(-1) * M_sun unit) in the SZ signal-mass scaling relation.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_M0,
                                   g_param_spec_double ("M0",
                                                        NULL,
                                                        "Reference mass",
                                                        1.0e13, G_MAXDOUBLE, 5.0e14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /*
   * SZ signal-mass scaling parameter: Asz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_VANDERLINDE_A_SZ, "A_{SZ}", "Asz",
                              1e-8,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_A_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /*
   * SZ signal-mass scaling parameter: Bsz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_VANDERLINDE_B_SZ, "B_{SZ}", "Bsz",
                              1e-8,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_B_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /*
   * SZ signal-mass scaling parameter: Csz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_VANDERLINDE_C_SZ, "C_{SZ}", "Csz",
                              1e-8,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_C_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /*
   * SZ signal-mass scaling parameter: Dsz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_VANDERLINDE_D_SZ, "D_{SZ}", "Dsz",
                              1e-8,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_D_SZ,
                              NCM_PARAM_TYPE_FIXED);

  parent_class->P              = &_nc_cluster_mass_vanderlinde_significance_m_p;
  parent_class->intP           = &_nc_cluster_mass_vanderlinde_intp;
  parent_class->P_limits       = &_nc_cluster_mass_vanderlinde_p_limits;
  parent_class->N_limits       = &_nc_cluster_mass_vanderlinde_n_limits;
  parent_class->resample       = &_nc_cluster_mass_vanderlinde_resample;
  parent_class->obs_len        = 1;
  parent_class->obs_params_len = 0;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);
}

static gdouble
_nc_cluster_mass_vanderlinde_significance_m_p_integrand (gdouble zeta, gpointer userdata)
{
  integrand_data *data          = (integrand_data *) userdata;
  NcClusterMassVanderlinde *msz = data->msz;
  const gdouble lnzeta          = log (zeta);
  const gdouble xi_mean         = sqrt (zeta * zeta + 3.0);
  const gdouble y               = data->xi[0] - xi_mean;
  const gdouble x               = lnzeta - data->mu;

  const gdouble result = exp (-y * y / 2.0 - x * x / data->D2_2) /
                         ((1.0 + erf (xi_mean / M_SQRT2)) * M_PI * D_SZ * zeta);

  return result;
}

static gdouble
_nc_cluster_mass_vanderlinde_significance_m_intp_integrand (gdouble zeta, gpointer userdata)
{
  integrand_data *data          = (integrand_data *) userdata;
  NcClusterMassVanderlinde *msz = data->msz;
  const gdouble lnzeta          = log (zeta);
  const gdouble xi_mean         = sqrt (zeta * zeta + 3.0);
  const gdouble a               = (xi_mean - msz->signif_obs_min);
  const gdouble x               = lnzeta - data->mu;

  const gdouble plnzeta = exp (-(x * x / (2.0 * D_SZ * D_SZ))) / (zeta * D_SZ * M_SQRT2 * M_SQRTPI);

  if (a >= 0.0)
    return plnzeta * (1.0 + erf (a / M_SQRT2)) / (1.0 + erf (xi_mean / M_SQRT2));
  else
    return plnzeta * erfc (-a / M_SQRT2) / (1.0 + erf (xi_mean / M_SQRT2));
}

static gdouble
_nc_cluster_mass_vanderlinde_significance_m_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *xi, const gdouble *xi_params)
{
  integrand_data data;
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  gdouble P, err;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  NCM_UNUSED (cosmo);

  data.msz       = msz;
  data.lnM       = lnM;
  data.z         = z;
  data.xi        = xi;
  data.xi_params = xi_params;

  /* new: optimizing */
  data.lnA        = log (A_SZ);
  data.lnM0       = log (msz->M0);
  data.ln1pz_1pz0 = log ((1.0 + z) / (1.0 + msz->z0));
  data.mu         = B_SZ * (lnM - data.lnM0) + C_SZ * data.ln1pz_1pz0 + data.lnA;
  data.D2_2       = 2.0 * D_SZ * D_SZ;

  F.function = &_nc_cluster_mass_vanderlinde_significance_m_p_integrand;
  F.params   = &data;

  {
    gdouble Pi, a, b;

    a = 0.0;
    b = xi[0];
    gsl_integration_qag (&F, a, b, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &Pi, &err);
    P = Pi;

    do {
      a  = b;
      b += xi[0];
      gsl_integration_qag (&F, a, b, P * NCM_DEFAULT_PRECISION, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &Pi, &err);
      P += Pi;
    } while (fabs (Pi / P) > NCM_DEFAULT_PRECISION);
  }

  ncm_memory_pool_return (w);

  return P;
}

static gdouble
_nc_cluster_mass_vanderlinde_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  integrand_data data;
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  gdouble P, err;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  NCM_UNUSED (cosmo);

  data.msz = msz;
  data.lnM = lnM;
  data.z   = z;

  /* new: optimizing */
  data.lnA        = log (A_SZ);
  data.lnM0       = log (msz->M0);
  data.ln1pz_1pz0 = log ((1.0 + z) / (1.0 + msz->z0));
  data.mu         = B_SZ * (lnM - data.lnM0) + C_SZ * data.ln1pz_1pz0 + data.lnA;
  data.D2_2       = 2.0 * D_SZ * D_SZ;

  F.function = &_nc_cluster_mass_vanderlinde_significance_m_intp_integrand;
  F.params   = &data;

  {
    gdouble Pi, a, b;

    a = 0.0;
    b = msz->signif_obs_max;
    gsl_integration_qag (&F, a, b, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &Pi, &err);
    P = Pi;

    do {
      a  = b;
      b += msz->signif_obs_max;
      gsl_integration_qag (&F, a, b, P * NCM_DEFAULT_PRECISION, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &Pi, &err);
      P += Pi;
    } while (fabs (Pi / P) > NCM_DEFAULT_PRECISION);
  }

  ncm_memory_pool_return (w);

  return P;
}

static gboolean
_nc_cluster_mass_vanderlinde_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *xi, const gdouble *xi_params, NcmRNG *rng)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  gdouble lnzeta, lnzeta_obs, zeta_obs, xi_mean;

  NCM_UNUSED (cosmo);
  NCM_UNUSED (xi_params);

  lnzeta = B_SZ * (lnM - log (msz->M0)) + C_SZ * log ((1.0 + z) / (1.0 + msz->z0)) + log (A_SZ);

  ncm_rng_lock (rng);
  lnzeta_obs = lnzeta + gsl_ran_gaussian (rng->r, D_SZ);

  zeta_obs = exp (lnzeta_obs);

  xi_mean = sqrt (zeta_obs * zeta_obs + 3.0);

  xi[0] = xi_mean + gsl_ran_gaussian (rng->r, 1.0);
  ncm_rng_unlock (rng);

  /*printf ("M = %e z = %.5g zeta = %.5g xi = %.5g xiobs = %.5g\n", exp(lnM), z, zeta_obs, xi_mean, xi[0]); */

  return (xi[0] >= msz->signif_obs_min);
}

static gdouble
_significance_to_mass (NcClusterMass *clusterm, gdouble z, gdouble xi)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  const gdouble zeta            = sqrt (xi * xi - 3.0);
  const gdouble lnzeta          = log (zeta);
  const gdouble lnM             = log (msz->M0) + (lnzeta - log (A_SZ) - C_SZ * (log (1.0 + z) - log (1.0 + msz->z0))) / B_SZ;

  return lnM;
}

static void
_nc_cluster_mass_vanderlinde_p_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *xi, const gdouble *xi_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  const gdouble xil             = GSL_MAX (xi[0] - 7.0, msz->signif_obs_min);
  const gdouble lnMl            = GSL_MAX (_significance_to_mass (clusterm, 2.0, xil) - 7.0 * D_SZ, log (2.0e14));
  const gdouble lnMu            = _significance_to_mass (clusterm, 0.0, xi[0] + 7.0) + 7.0 * D_SZ;

  NCM_UNUSED (cosmo);
  NCM_UNUSED (xi_params);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static void
_nc_cluster_mass_vanderlinde_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  const gdouble lnMl            = GSL_MAX (_significance_to_mass (clusterm, 2.0, msz->signif_obs_min) - 7.0 * D_SZ, log (2.0e14));
  const gdouble lnMu            = _significance_to_mass (clusterm, 0.0, msz->signif_obs_max) + 7.0 * D_SZ;

  NCM_UNUSED (cosmo);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

