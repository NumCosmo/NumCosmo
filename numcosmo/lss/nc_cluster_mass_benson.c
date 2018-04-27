/***************************************************************************
 *            nc_cluster_mass_benson.c
 *
 *  Tue July 9 14:18:11 2012
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
 * SECTION:nc_cluster_mass_benson
 * @title: NcClusterMassBenson
 * @short_description: Sunyaev-Zel'dovich cluster mass distribution.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_benson.h"
#include "math/integral.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcClusterMassBenson, nc_cluster_mass_benson, NC_TYPE_CLUSTER_MASS);

#define VECTOR (NCM_MODEL (msz)->params)
#define A_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_A_SZ))
#define B_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_B_SZ))
#define C_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_C_SZ))
#define D_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_D_SZ))

enum
{
  PROP_0,
  PROP_SIGNIFICANCE_OBS_MIN,
  PROP_SIGNIFICANCE_OBS_MAX,
  PROP_Z0,
  PROP_M0,
  PROP_SIZE,
};

static void
nc_cluster_mass_benson_init (NcClusterMassBenson *msz)
{
  msz->signif_obs_min = 0.0;
  msz->signif_obs_max = 0.0;
  msz->z0 = 0.0;
  msz->M0 = 0.0;
}

static void
_nc_cluster_mass_benson_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_BENSON (object));

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
_nc_cluster_mass_benson_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_BENSON (object));

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
_nc_cluster_mass_benson_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_benson_parent_class)->finalize (object);
}

guint _nc_cluster_mass_benson_obs_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 1; }
guint _nc_cluster_mass_benson_obs_params_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 0; }
static gdouble _nc_cluster_mass_benson_significance_m_p (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, const gdouble *xi, const gdouble *xi_params);
static gdouble _nc_cluster_mass_benson_intp (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z);
static void _nc_cluster_mass_benson_p_limits (NcClusterMass *clusterm, NcHICosmo *model, const gdouble *xi, const gdouble *xi_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_benson_n_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_lower, gdouble *lnM_upper);
static gboolean _nc_cluster_mass_benson_resample (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *xi, const gdouble *xi_params, NcmRNG *rng);


static void
nc_cluster_mass_benson_class_init (NcClusterMassBensonClass *klass)
{
  GObjectClass* object_class       = G_OBJECT_CLASS (klass);
  NcClusterMassClass* parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  parent_class->P = &_nc_cluster_mass_benson_significance_m_p;
  parent_class->intP = &_nc_cluster_mass_benson_intp;
  parent_class->P_limits = &_nc_cluster_mass_benson_p_limits;
  parent_class->N_limits = &_nc_cluster_mass_benson_n_limits;
  parent_class->resample = &_nc_cluster_mass_benson_resample;
  parent_class->obs_len = &_nc_cluster_mass_benson_obs_len;
  parent_class->obs_params_len = &_nc_cluster_mass_benson_obs_params_len;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);

  model_class->set_property = &_nc_cluster_mass_benson_set_property;
  model_class->get_property = &_nc_cluster_mass_benson_get_property;
  object_class->finalize = _nc_cluster_mass_benson_finalize;

  ncm_model_class_set_name_nick (model_class, "Benson - SZ", "Benson_SZ");
  ncm_model_class_add_params (model_class, 4, 0, PROP_SIZE);

  /**
   * NcClusterMassBenson:signif_obs_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGNIFICANCE_OBS_MIN,
                                   g_param_spec_double ("signif-obs-min",
                                                        NULL,
                                                        "Minimum obsevational significance",
                                                        2.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassBenson:signif_obs_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGNIFICANCE_OBS_MAX,
                                   g_param_spec_double ("signif-obs-max",
                                                        NULL,
                                                        "Maximum obsevational significance",
                                                        2.0, G_MAXDOUBLE, 40.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassBenson:z0:
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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); 
  /**
   * NcClusterMassBenson:M0:
   *
   * Reference mass (in h^(-1) * M_sun unit) in the SZ signal-mass scaling relation.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_M0,
                                   g_param_spec_double ("M0",
                                                        NULL,
                                                        "Reference mass",
                                                        1.0e13, G_MAXDOUBLE, 3.0e14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassBenson:Asz:
   * 
   * Slope of the SZ signal-mass scaling relation.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_A_SZ, "A_{SZ}", "Asz",
                              1e-8,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_BENSON_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_DEFAULT_A_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassBenson:Bsz:
   * 
   * SZ signal-mass scaling parameter.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_B_SZ, "B_{SZ}", "Bsz",
                              1e-8,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_BENSON_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_DEFAULT_B_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassBenson:Csz:
   * 
   * SZ signal-mass scaling parameter.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_C_SZ, "C_{SZ}", "Csz",
                              1e-8,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_BENSON_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_DEFAULT_C_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassBenson:Dsz:
   * 
   * Standard deviation of the SZ signal-mass scaling relation.
   * $D_sz \in [0.01, 2.0]$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_D_SZ, "D_{SZ}", "Dsz",
                              1e-2,  2.0, 1.0e-2,
                              NC_CLUSTER_MASS_BENSON_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_DEFAULT_D_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
}


typedef struct _integrand_data
{
  NcClusterMassBenson *msz;
  NcHICosmo *model;
  const gdouble *xi_params;
  gdouble z;
  gdouble lnM;
  const gdouble *xi;
  gdouble lnA;
  gdouble lnM0;
  gdouble lnE_E0;  /* ln (E/E0)  */
  gdouble mu;
  gdouble D2_2;
} integrand_data;

static gdouble
_nc_cluster_mass_benson_xi_mean (gdouble zeta)
{
  const gdouble xi_mean = (zeta <= 1.000001) ? zeta : sqrt (zeta * zeta + 3.0);

  return xi_mean; 
}

static gdouble
_nc_cluster_mass_benson_int_xi_cut_inf (gdouble xi_mean, gdouble xi_cut)
{
  const gdouble a = (xi_cut - xi_mean) / M_SQRT2;

  if (a < 0.0)
    return (1.0 - erf (a)) * 0.5;
  else
    return erfc (a) * 0.5;
}

static gdouble
_nc_cluster_mass_benson_significance_m_p_integrand (gdouble zeta, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassBenson *msz = data->msz;
  const gdouble lnzeta = log (zeta);
  const gdouble xi_mean = _nc_cluster_mass_benson_xi_mean (zeta);
  const gdouble y = data->xi[0] - xi_mean;
  const gdouble x = lnzeta - data->mu;
  const gdouble exp_arg = - y * y / 2.0 - x * x / data->D2_2;

  if (exp_arg < GSL_LOG_DBL_MIN)
    return 0.0;
  else
  {
    const gdouble result = exp (exp_arg) / (2.0 * M_PI * D_SZ * zeta);
    return result;
  }
}

static gdouble
_nc_cluster_mass_benson_significance_m_intp_integrand (gdouble zeta, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassBenson *msz = data->msz;
  const gdouble lnzeta = log (zeta);
  const gdouble x = lnzeta - data->mu;
  const gdouble exp_arg = -(x * x / (2.0 * D_SZ * D_SZ));
  const gdouble xi_mean = _nc_cluster_mass_benson_xi_mean (zeta);
  const gdouble plnzeta = exp (exp_arg) / (zeta * M_SQRT2 * M_SQRTPI * D_SZ);
  const gdouble xi_zeta_int = _nc_cluster_mass_benson_int_xi_cut_inf (xi_mean, msz->signif_obs_min); 

  //printf ("nhaca % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", zeta, exp_arg, xi_mean, plnzeta, xi_zeta_int);

  return plnzeta * xi_zeta_int;
}

static gdouble
_nc_cluster_mass_benson_significance_m_p (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, const gdouble *xi, const gdouble *xi_params)
{
  integrand_data data;
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  gdouble P, err;
  const gdouble E0 = nc_hicosmo_E (model, msz->z0); 
  const gdouble E = nc_hicosmo_E (model, z);
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  data.msz = msz;
  data.model = model;
  data.lnM = lnM;
  data.z = z;
  data.xi = xi;
  data.xi_params = xi_params;

  data.lnA = log (A_SZ);
  data.lnM0 = log (msz->M0);
  data.lnE_E0 = log (E / E0);
  data.mu = B_SZ * (lnM - data.lnM0) + C_SZ * data.lnE_E0 + data.lnA;
  data.D2_2 = 2.0 * D_SZ * D_SZ;

  F.function = &_nc_cluster_mass_benson_significance_m_p_integrand;
  F.params = &data;

  {
    gdouble Pi, a, b;
//    a = 0.25;
//    b = 1.0;
    a = 0.0;
    b = 1.0;
    gsl_integration_qag (&F, a, b, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &Pi, &err);
    P = Pi;
//    b = 2.0;
    b = 1.0;
    do {
      a = b;
      b += xi[0];
      gsl_integration_qag (&F, a, b, P * NCM_DEFAULT_PRECISION, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &Pi, &err); 
      P += Pi;
    } while (fabs(Pi/P) > NCM_DEFAULT_PRECISION);
  }

  ncm_memory_pool_return (w);
  
  return P;
}

static gdouble
_nc_cluster_mass_benson_intp (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z)
{
  integrand_data data;
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  gdouble P, err;
  const gdouble E0 = nc_hicosmo_E (model, msz->z0); 
  const gdouble E = nc_hicosmo_E (model, z);
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  data.msz = msz;
  data.model = model;
  data.lnM = lnM;
  data.z = z;

  data.lnA = log (A_SZ);
  data.lnM0 = log (msz->M0);
  data.lnE_E0 = log (E / E0);
  data.mu = B_SZ * (lnM - data.lnM0) + C_SZ * data.lnE_E0 + data.lnA;
  data.D2_2 = 2.0 * D_SZ * D_SZ;

  F.function = &_nc_cluster_mass_benson_significance_m_intp_integrand;
  F.params = &data;
  
  {
    gdouble Pi, a, b;
//    a = 0.25;
//    b = 1.0;
    a = 0.0;
    b = 1.0;
    gsl_integration_qag (&F, a, b, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &Pi, &err);
    P = Pi;
//    b = 2.0;
    b = 1.0;
//printf ("int_p[0,1] % 8.5g % 8.5g : % 8.5g % 8.5g\n", exp (lnM), z, P, data.mu);
    do {
      a = b;
      b += msz->signif_obs_max;
      gsl_integration_qag (&F, a, b, P * NCM_DEFAULT_PRECISION, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &Pi, &err); 
      P += Pi;
    } while (fabs(Pi/P) > NCM_DEFAULT_PRECISION);
  }

  ncm_memory_pool_return (w);

//printf ("int_p[2,-] % 8.5g % 8.5g : % 8.5g % 8.5g\n", exp (lnM), z, P, data.mu);

  
  return P;
}

static gboolean
_nc_cluster_mass_benson_resample (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *xi, const gdouble *xi_params, NcmRNG *rng)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  gdouble lnzeta, lnzeta_obs, zeta_obs, xi_mean;
  const gdouble E0 = nc_hicosmo_E (model, msz->z0);
  const gdouble E = nc_hicosmo_E (model, z);

  NCM_UNUSED (xi_params);
  
  lnzeta = B_SZ * (lnM - log (msz->M0)) + C_SZ * log (E / E0) + log (A_SZ);
 
  {
    gboolean ret;
    ncm_rng_lock (rng);
    lnzeta_obs = lnzeta + gsl_ran_gaussian (rng->r, D_SZ);
    zeta_obs = exp (lnzeta_obs);
    if (zeta_obs > 1.0 && zeta_obs < 2.0)
      ret = FALSE;
    else
    {
      xi_mean = _nc_cluster_mass_benson_xi_mean (zeta_obs);
      xi[0] = xi_mean + gsl_ran_gaussian (rng->r, 1.0);

      //printf("M = %e z = %.5g zeta = %.5g xi = %.5g xiobs = %.5g | xiobs_min = %.5g\n", exp(lnM), z, zeta_obs, xi_mean, xi[0], msz->signif_obs_min);

      ret = (xi[0] >= msz->signif_obs_min);
    }
    ncm_rng_unlock (rng);
    return ret;
  }
}

static gdouble
_significance_to_zeta (NcClusterMass *clusterm, NcHICosmo *model, gdouble z, gdouble xi)
{
  NCM_UNUSED (clusterm);
  NCM_UNUSED (model);
  NCM_UNUSED (z);
  
  return sqrt (xi * xi - 3.0);
}

static gdouble
_zeta_to_mass (NcClusterMass *clusterm, NcHICosmo *model, gdouble z, gdouble zeta)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  const gdouble E0 = nc_hicosmo_E (model, msz->z0); 
  const gdouble E = nc_hicosmo_E (model, z);
  const gdouble lnzeta = log (zeta);
  const gdouble lnM = log (msz->M0) + (lnzeta - log (A_SZ) - C_SZ * log (E / E0)) / B_SZ;

  //printf("z= %.10g xi = %.10g lnM = %.10g\n", z, xi, lnM);
  return lnM;
}

static void
_nc_cluster_mass_benson_p_limits (NcClusterMass *clusterm, NcHICosmo *model, const gdouble *xi, const gdouble *xi_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  const gdouble xil = GSL_MAX (xi[0] - 7.0, msz->signif_obs_min);
  const gdouble zetal = _significance_to_zeta (clusterm, model, 2.0, xil) - 7.0 * D_SZ;
  const gdouble lnMl = GSL_MAX (_zeta_to_mass (clusterm, model, 2.0, zetal), log (NC_CLUSTER_MASS_BENSON_M_LOWER_BOUND));

  const gdouble xiu = xi[0] + 7.0;
  const gdouble zetau = _significance_to_zeta (clusterm, model, 0.0, xiu) + 7.0 * D_SZ;
  const gdouble lnMu = _zeta_to_mass (clusterm, model, 0.0, zetau);

  NCM_UNUSED (xi_params);
  
  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static void
_nc_cluster_mass_benson_n_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);

  const gdouble xil = msz->signif_obs_min;
  const gdouble zetal = _significance_to_zeta (clusterm, model, 2.0, xil) - 7.0 * D_SZ;
  const gdouble lnMl = GSL_MAX (_zeta_to_mass (clusterm, model, 2.0, zetal), log (NC_CLUSTER_MASS_BENSON_M_LOWER_BOUND));

  const gdouble xiu = msz->signif_obs_max;
  const gdouble zetau = _significance_to_zeta (clusterm, model, 0.0, xiu) + 7.0 * D_SZ;
  const gdouble lnMu = _zeta_to_mass (clusterm, model, 0.0, zetau);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;  
}

