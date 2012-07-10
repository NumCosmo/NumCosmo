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

G_DEFINE_TYPE (NcClusterMassVanderlinde, nc_cluster_mass_vanderlinde, NC_TYPE_CLUSTER_MASS);

enum
{
  PROP_0,
  PROP_SIGNIFICANCE_OBS_MIN,
  PROP_SIGNIFICANCE_OBS_MAX,
  PROP_Z0,
  PROP_M0,
  PROP_A,
  PROP_B,
  PROP_C,
  PROP_D,
};

typedef struct _integrand_data
{
  NcClusterMassVanderlinde *msz;
  gdouble *xi_params;
  gdouble z;
  gdouble lnM;
  gdouble *xi;
} integrand_data;

static gdouble
_p_lnzeta_lnm (NcClusterMassVanderlinde *msz, gdouble lnM, gdouble z, gdouble zeta)
{
  gdouble ln_zeta = log (zeta);
  const gdouble lnM0 = log (msz->M0);

  gdouble x = ln_zeta - msz->B * (lnM - lnM0) - msz->C * (log (1.0 + z) - log(1.0 + msz->z0)) - log (msz->A);
   
  return (M_SQRT2 / M_SQRTPI) * exp (-(x * x / (2.0 * msz->D * msz->D))) / (zeta * msz->D * (1.0 + erf((ln_zeta - x) / (M_SQRT2 * msz->D))));
}

static gdouble
_p_significance_zeta (NcClusterMassVanderlinde *msz, gdouble *xi, gdouble zeta)
{
  gdouble xi_mean = sqrt (zeta * zeta + 3.0);
  gdouble y = xi[0] - xi_mean;

  gdouble p_xi_zeta = (M_SQRT2 / M_SQRTPI) * exp (- y * y / 2.0) / 
	(1.0 + erf (xi_mean / M_SQRT2));

  return p_xi_zeta;
}

static gdouble
_intp_significance_zeta (NcClusterMassVanderlinde *msz, gdouble zeta)
{
  const gdouble xi_mean = sqrt (zeta * zeta + 3.0);
  const gdouble a = (xi_mean - msz->signif_obs_min);

  if (a >= 0.0)
	return (1.0 + erf (a / M_SQRT2)) / (1.0 + erf (xi_mean / M_SQRT2));
  else
	return erfc (-a / M_SQRT2) / (1.0 + erf (xi_mean / M_SQRT2));
}

static gdouble
_nc_cluster_mass_vanderlinde_significance_m_p_integrand (gdouble zeta, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassVanderlinde *msz = data->msz;
  const gdouble p_lnzeta_lnm = _p_lnzeta_lnm (msz, data->lnM, data->z, zeta);
  const gdouble p_xi_zeta = _p_significance_zeta (msz, data->xi, zeta);
  
  return p_lnzeta_lnm * p_xi_zeta;
}

static gdouble
_nc_cluster_mass_vanderlinde_significance_m_intp_integrand (gdouble zeta, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassVanderlinde *msz = data->msz;
  const gdouble p_lnzeta_lnm = _p_lnzeta_lnm (msz, data->lnM, data->z, zeta);
  const gdouble intp_zeta = _intp_significance_zeta (msz, zeta);
  
  return p_lnzeta_lnm * intp_zeta;
}

static gdouble
_nc_cluster_mass_vanderlinde_significance_m_p (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *xi, gdouble *xi_params)
{
  integrand_data data;
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  gdouble P, err;
  gsl_function F;
  gsl_integration_workspace **w = nc_integral_get_workspace ();

  data.msz = msz;
  data.lnM = lnM;
  data.z = z;
  data.xi = xi;
  data.xi_params = xi_params;

  F.function = &_nc_cluster_mass_vanderlinde_significance_m_p_integrand;
  F.params = &data;

  gsl_integration_qagiu (&F, 0.0, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, *w, &P, &err);
  
  return P;
}

static gdouble
_nc_cluster_mass_vanderlinde_intp (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z)
{
  integrand_data data;
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  gdouble P, err;
  gsl_function F;
  gsl_integration_workspace **w = nc_integral_get_workspace ();

  data.msz = msz;
  data.lnM = lnM;
  data.z = z;

  F.function = &_nc_cluster_mass_vanderlinde_significance_m_intp_integrand;
  F.params = &data;

  gsl_integration_qagiu (&F, 0.0, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, *w, &P, &err);
  
  return P;
}

static gboolean
_nc_cluster_mass_vanderlinde_resample (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *xi, gdouble *xi_params)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  gsl_rng *rng = ncm_get_rng ();
  gdouble lnzeta, lnzeta_obs, zeta_obs, xi_mean;

  lnzeta = msz->B * (lnM - log (msz->M0)) + msz->C * log ((1.0 + z) / (1.0 + msz->z0)) + log (msz->A);
  
  lnzeta_obs = lnzeta + gsl_ran_gaussian (rng, msz->D);
  
  zeta_obs = exp (lnzeta_obs);

  xi_mean = sqrt (zeta_obs * zeta_obs + 3.0);

  xi[0] = xi_mean + gsl_ran_gaussian (rng, 1.0);

  printf("M = %e z = %.5g zeta = %.5g xi = %.5g xiobs = %.5g\n", exp(lnM), z, zeta_obs, xi_mean, xi[0]);

  return (xi[0] >= msz->signif_obs_min);
}

static gdouble
_significance_to_mass (NcClusterMass *clusterm, gdouble z, gdouble xi)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  
  const gdouble zeta = sqrt (xi * xi - 3.0);
  const gdouble lnzeta = log(zeta);
  const gdouble lnM = log (msz->M0) + (lnzeta - log(msz->A) - msz->C * (log(1.0 + z) - log(1.0 + msz->z0))) / msz->B;
  
  return lnM;
}

static void
_nc_cluster_mass_vanderlinde_p_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *xi, gdouble *xi_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  const gdouble xil = GSL_MAX (xi[0] - 7.0, msz->signif_obs_min);
  const gdouble lnMl = GSL_MAX (_significance_to_mass (clusterm, 2.0, xil) - 7.0 * msz->D, log (2.0e14));
  const gdouble lnMu = _significance_to_mass (clusterm, 0.0, xi[0] + 7.0) + 7.0 * msz->D;

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static void
_nc_cluster_mass_vanderlinde_n_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassVanderlinde *msz = NC_CLUSTER_MASS_VANDERLINDE (clusterm);
  const gdouble lnMl = GSL_MAX (_significance_to_mass (clusterm, 2.0, msz->signif_obs_min) - 7.0 * msz->D, log (2.0e14));
  const gdouble lnMu = _significance_to_mass (clusterm, 0.0, msz->signif_obs_max) + 7.0 * msz->D;

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;  
}

guint _nc_cluster_mass_vanderlinde_obs_len (NcClusterMass *clusterm) { return 1; }
guint _nc_cluster_mass_vanderlinde_obs_params_len (NcClusterMass *clusterm) { return 0; }

static void
_nc_cluster_mass_vanderlinde_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
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
	case PROP_A:
	  msz->A = g_value_get_double (value);
	  break;
	case PROP_B:
	  msz->B = g_value_get_double (value);
	  break;
	case PROP_C:
	  msz->C = g_value_get_double (value);
	  break;
	case PROP_D:
	  msz->D = g_value_get_double (value);
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
	case PROP_A:
      g_value_set_double (value, msz->A);
      break;
	case PROP_B:
      g_value_set_double (value, msz->B);
      break;
	case PROP_C:
      g_value_set_double (value, msz->C);
      break;
	case PROP_D:
      g_value_set_double (value, msz->D);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_mass_vanderlinde_init (NcClusterMassVanderlinde *msz)
{
  msz->signif_obs_min = 0.0;
  msz->signif_obs_max = 0.0;
  msz->z0 = 0.0;
  msz->M0 = 0.0;
  msz->A = 0.0;
  msz->B = 0.0;
  msz->C = 0.0;
  msz->D = 0.0;
}

static void
_nc_cluster_mass_vanderlinde_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_vanderlinde_parent_class)->finalize (object);
}

static void
nc_cluster_mass_vanderlinde_class_init (NcClusterMassVanderlindeClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterMassClass* parent_class = NC_CLUSTER_MASS_CLASS (klass);

  parent_class->P = &_nc_cluster_mass_vanderlinde_significance_m_p;
  parent_class->intP = &_nc_cluster_mass_vanderlinde_intp;
  parent_class->P_limits = &_nc_cluster_mass_vanderlinde_p_limits;
  parent_class->N_limits = &_nc_cluster_mass_vanderlinde_n_limits;
  parent_class->resample = &_nc_cluster_mass_vanderlinde_resample;
  parent_class->obs_len = &_nc_cluster_mass_vanderlinde_obs_len;
  parent_class->obs_params_len = &_nc_cluster_mass_vanderlinde_obs_params_len;

  parent_class->impl = NC_CLUSTER_MASS_IMPL_ALL;
  
  object_class->finalize = _nc_cluster_mass_vanderlinde_finalize;
  object_class->set_property = &_nc_cluster_mass_vanderlinde_set_property;
  object_class->get_property = &_nc_cluster_mass_vanderlinde_get_property;

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
                                                        1.0e13, G_MAXDOUBLE, 5.0e14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassVanderlinde:A:
   *
   * SZ signal-mass scaling parameter: A.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("A",
                                                        NULL,
                                                        "A",
                                                        0.0, G_MAXDOUBLE, 6.01,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassVanderlinde:B:
   *
   * SZ signal-mass scaling parameter: B.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B,
                                   g_param_spec_double ("B",
                                                        NULL,
                                                        "B",
                                                        0.0, G_MAXDOUBLE, 1.31,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassVanderlinde:C:
   *
   * SZ signal-mass scaling parameter: C.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_C,
                                   g_param_spec_double ("C",
                                                        NULL,
                                                        "C",
                                                        0.0, G_MAXDOUBLE, 1.6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassVanderlinde:D:
   *
   * SZ signal-mass scaling standard deviation: D.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_D,
                                   g_param_spec_double ("D",
                                                        NULL,
                                                        "D",
                                                        0.0, G_MAXDOUBLE, 0.21,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

