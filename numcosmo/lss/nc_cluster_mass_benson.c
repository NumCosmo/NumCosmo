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
 * @title: SZ Cluster Abundance Mass Distribution
 * @short_description: Sunyaev-Zel'dovich FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>
#include <gsl/gsl_integration.h>

G_DEFINE_TYPE (NcClusterMassBenson, nc_cluster_mass_benson, NC_TYPE_CLUSTER_MASS);

enum
{
  PROP_0,
  PROP_SIGNIFICANCE_OBS_MIN,
  PROP_SIGNIFICANCE_OBS_MAX,
  PROP_Z0,
  PROP_M0,
  PROP_ASZ,
  PROP_BSZ,
  PROP_CSZ,
  PROP_DSZ,
};

typedef struct _integrand_data
{
  NcClusterMassBenson *msz;
  NcHICosmo *model;
  gdouble *xi_params;
  gdouble z;
  gdouble lnM;
  gdouble *xi;
} integrand_data;

static gdouble
_p_lnzeta_lnm (NcClusterMassBenson *msz, NcHICosmo *model, gdouble lnM, gdouble z, gdouble zeta)
{
  gdouble ln_zeta = log (zeta);
  const gdouble E0 = nc_hicosmo_E (model, msz->z0); 
  const gdouble E = nc_hicosmo_E (model, z); 
  const gdouble lnM0 = log (msz->M0);

  const gdouble x = ln_zeta - msz->Bsz * (lnM - lnM0) - msz->Csz * log(E / E0) - log (msz->Asz);
   
  return (M_SQRT2 / M_SQRTPI) * exp (-(x * x / (2.0 * msz->Dsz * msz->Dsz))) / (zeta * msz->Dsz * (1.0 + erf((ln_zeta - x) / (M_SQRT2 * msz->Dsz))));
}

static gdouble
_p_significance_zeta (NcClusterMassBenson *msz, gdouble *xi, gdouble zeta)
{
  gdouble xi_mean = sqrt (zeta * zeta + 3.0);
  gdouble y = xi[0] - xi_mean;

  gdouble p_xi_zeta = (M_SQRT2 / M_SQRTPI) * exp (- y * y / 2.0) / 
	(1.0 + erf (xi_mean / M_SQRT2));

  return p_xi_zeta;
}

static gdouble
_intp_significance_zeta (NcClusterMassBenson *msz, gdouble zeta)
{
  const gdouble xi_mean = sqrt (zeta * zeta + 3.0);
  const gdouble a = (xi_mean - msz->signif_obs_min);

  if (a >= 0.0)
	return (1.0 + erf (a / M_SQRT2)) / (1.0 + erf (xi_mean / M_SQRT2));
  else
	return erfc (-a / M_SQRT2) / (1.0 + erf (xi_mean / M_SQRT2));
}

static gdouble
_nc_cluster_mass_benson_significance_m_p_integrand (gdouble zeta, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassBenson *msz = data->msz;
  const gdouble p_lnzeta_lnm = _p_lnzeta_lnm (msz, data->model, data->lnM, data->z, zeta);
  const gdouble p_xi_zeta = _p_significance_zeta (msz, data->xi, zeta);
  
  return p_lnzeta_lnm * p_xi_zeta;
}

static gdouble
_nc_cluster_mass_benson_significance_m_intp_integrand (gdouble zeta, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassBenson *msz = data->msz;
  const gdouble p_lnzeta_lnm = _p_lnzeta_lnm (msz, data->model, data->lnM, data->z, zeta);
  const gdouble intp_zeta = _intp_significance_zeta (msz, zeta);
  
  return p_lnzeta_lnm * intp_zeta;
}

static gdouble
_nc_cluster_mass_benson_significance_m_p (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *xi, gdouble *xi_params)
{
  integrand_data data;
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  gdouble P, err;
  gsl_function F;
  gsl_integration_workspace **w = nc_integral_get_workspace ();

  data.msz = msz;
  data.model = model;
  data.lnM = lnM;
  data.z = z;
  data.xi = xi;
  data.xi_params = xi_params;

  F.function = &_nc_cluster_mass_benson_significance_m_p_integrand;
  F.params = &data;

  gsl_integration_qagiu (&F, 0.0, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, *w, &P, &err);

  ncm_memory_pool_return (w);
  
  return P;
}

static gdouble
_nc_cluster_mass_benson_intp (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z)
{
  integrand_data data;
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  gdouble P, err;
  gsl_function F;
  gsl_integration_workspace **w = nc_integral_get_workspace ();

  data.msz = msz;
  data.model = model;
  data.lnM = lnM;
  data.z = z;

  F.function = &_nc_cluster_mass_benson_significance_m_intp_integrand;
  F.params = &data;

  gsl_integration_qagiu (&F, 0.0, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, *w, &P, &err);

  ncm_memory_pool_return (w);
  
  return P;
}

static gboolean
_nc_cluster_mass_benson_resample (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *xi, gdouble *xi_params)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  gsl_rng *rng = ncm_get_rng ();
  gdouble lnzeta, lnzeta_obs, zeta_obs, xi_mean;
  const gdouble E0 = nc_hicosmo_E (model, msz->z0);
  const gdouble E = nc_hicosmo_E (model, z);

  lnzeta = msz->Bsz * (lnM - log (msz->M0)) + msz->Csz * log (E / E0) + log (msz->Asz);
  
  lnzeta_obs = lnzeta + gsl_ran_gaussian (rng, msz->Dsz);
  
  zeta_obs = exp (lnzeta_obs);

  xi_mean = sqrt (zeta_obs * zeta_obs + 3.0);

  xi[0] = xi_mean + gsl_ran_gaussian (rng, 1.0);

  printf("M = %e z = %.5g zeta = %.5g xi = %.5g xiobs = %.5g\n", exp(lnM), z, zeta_obs, xi_mean, xi[0]);

  return (xi[0] >= msz->signif_obs_min);
}

static gdouble
_significance_to_mass (NcClusterMass *clusterm, NcHICosmo *model, gdouble z, gdouble xi)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  const gdouble E0 = nc_hicosmo_E (model, msz->z0); 
  const gdouble E = nc_hicosmo_E (model, z);
  const gdouble zeta = sqrt (xi * xi - 3.0);
  const gdouble lnzeta = log(zeta);
  const gdouble lnM = log (msz->M0) + (lnzeta - log(msz->Asz) - msz->Csz * log(E / E0)) / msz->Bsz;
  
  return lnM;
}

static void
_nc_cluster_mass_benson_p_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *xi, gdouble *xi_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  const gdouble xil = GSL_MAX (xi[0] - 7.0, msz->signif_obs_min);
  const gdouble lnMl = GSL_MAX (_significance_to_mass (clusterm, model, 2.0, xil) - 7.0 * msz->Dsz, log (2.0e14));
  const gdouble lnMu = _significance_to_mass (clusterm, model, 0.0, xi[0] + 7.0) + 7.0 * msz->Dsz;

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static void
_nc_cluster_mass_benson_n_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  const gdouble lnMl = GSL_MAX (_significance_to_mass (clusterm, model, 2.0, msz->signif_obs_min) - 7.0 * msz->Dsz, log (2.0e14));
  const gdouble lnMu = _significance_to_mass (clusterm, model, 0.0, msz->signif_obs_max) + 7.0 * msz->Dsz;

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;  
}

guint _nc_cluster_mass_benson_obs_len (NcClusterMass *clusterm) { return 1; }
guint _nc_cluster_mass_benson_obs_params_len (NcClusterMass *clusterm) { return 0; }

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
	case PROP_ASZ:
	  msz->Asz = g_value_get_double (value);
	  break;
	case PROP_BSZ:
	  msz->Bsz = g_value_get_double (value);
	  break;
	case PROP_CSZ:
	  msz->Csz = g_value_get_double (value);
	  break;
	case PROP_DSZ:
	  msz->Dsz = g_value_get_double (value);
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
	case PROP_ASZ:
      g_value_set_double (value, msz->Asz);
      break;
	case PROP_BSZ:
      g_value_set_double (value, msz->Bsz);
      break;
	case PROP_CSZ:
      g_value_set_double (value, msz->Csz);
      break;
	case PROP_DSZ:
      g_value_set_double (value, msz->Dsz);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_mass_benson_init (NcClusterMassBenson *msz)
{
  msz->signif_obs_min = 0.0;
  msz->signif_obs_max = 0.0;
  msz->z0 = 0.0;
  msz->M0 = 0.0;
  msz->Asz = 0.0;
  msz->Bsz = 0.0;
  msz->Csz = 0.0;
  msz->Dsz = 0.0;
}

static void
_nc_cluster_mass_benson_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_benson_parent_class)->finalize (object);
}

static void
nc_cluster_mass_benson_class_init (NcClusterMassBensonClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterMassClass* parent_class = NC_CLUSTER_MASS_CLASS (klass);

  parent_class->P = &_nc_cluster_mass_benson_significance_m_p;
  parent_class->intP = &_nc_cluster_mass_benson_intp;
  parent_class->P_limits = &_nc_cluster_mass_benson_p_limits;
  parent_class->N_limits = &_nc_cluster_mass_benson_n_limits;
  parent_class->resample = &_nc_cluster_mass_benson_resample;
  parent_class->obs_len = &_nc_cluster_mass_benson_obs_len;
  parent_class->obs_params_len = &_nc_cluster_mass_benson_obs_params_len;

  parent_class->impl = NC_CLUSTER_MASS_IMPL_ALL;
  
  object_class->finalize = _nc_cluster_mass_benson_finalize;
  object_class->set_property = &_nc_cluster_mass_benson_set_property;
  object_class->get_property = &_nc_cluster_mass_benson_get_property;

  /**
   * NcClusterMassVanderlinde:signif_obs_min:
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
   * NcClusterMassVanderlinde:signif_obs_max:
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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); 
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
                                                        1.0e13, G_MAXDOUBLE, 3.0e14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassVanderlinde:Asz:
   *
   * SZ signal-mass scaling parameter: Asz.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_ASZ,
                                   g_param_spec_double ("Asz",
                                                        NULL,
                                                        "Asz",
                                                        0.0, G_MAXDOUBLE, 5.58,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassVanderlinde:Bsz:
   *
   * SZ signal-mass scaling parameter: Bsz.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_BSZ,
                                   g_param_spec_double ("Bsz",
                                                        NULL,
                                                        "Bsz",
                                                        0.0, G_MAXDOUBLE, 1.32,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassVanderlinde:Csz:
   *
   * SZ signal-mass scaling parameter: Csz.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_CSZ,
                                   g_param_spec_double ("Csz",
                                                        NULL,
                                                        "Csz",
                                                        0.0, G_MAXDOUBLE, 0.87,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassVanderlinde:Dsz:
   *
   * SZ signal-mass scaling standard deviation: Dsz.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DSZ,
                                   g_param_spec_double ("Dsz",
                                                        NULL,
                                                        "Dsz",
                                                        0.0, G_MAXDOUBLE, 0.24,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}

