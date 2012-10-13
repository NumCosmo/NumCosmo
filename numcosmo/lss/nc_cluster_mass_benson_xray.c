/***************************************************************************
 *            nc_cluster_mass_benson_xray.c
 *
 *  Tue July 9 17:01:24 2012
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
 * SECTION:nc_cluster_mass_benson_xray
 * @title: SZ and X ray Cluster Abundance Mass Distributions
 * @short_description: Sunyaev-Zel'dovich FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>
#include <gsl/gsl_integration.h>

G_DEFINE_TYPE (NcClusterMassBensonXRay, nc_cluster_mass_benson_xray, NC_TYPE_CLUSTER_MASS);

#define VECTOR (NCM_MODEL (msz)->params)
#define A_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_A_SZ))
#define B_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_B_SZ))
#define C_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_C_SZ))
#define D_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_D_SZ))
#define A_X    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_A_X))
#define B_X    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_B_X))
#define C_X    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_C_X))
#define D_X    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_D_X))

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
  NcClusterMassBensonXRay *msz;
  gdouble *obs_params; //*xi_params;
  gdouble z;
  gdouble lnM;
  gdouble *obs; //*xi;
} integrand_data;


guint _nc_cluster_mass_benson_xray_obs_len (NcClusterMass *clusterm) { return 2; }
guint _nc_cluster_mass_benson_xray_obs_params_len (NcClusterMass *clusterm) { return 0; }

static void
_nc_cluster_mass_benson_xray_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcClusterMassBensonXRay *mbxr = NC_CLUSTER_MASS_BENSON_XRAY (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_BENSON_XRAY (object));

  switch (prop_id)
  {
    case PROP_SIGNIFICANCE_OBS_MIN:
      mbxr->signif_obs_min = g_value_get_double (value);
      break;
	case PROP_SIGNIFICANCE_OBS_MAX:
      mbxr->signif_obs_max = g_value_get_double (value);
      break;
	case PROP_Z0:
	  mbxr->z0 = g_value_get_double (value);
	  break;
	case PROP_M0:
	  mbxr->M0 = g_value_get_double (value);
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_benson_xray_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassBensonXRay *mbrx = NC_CLUSTER_MASS_BENSON_XRAY (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_BENSON_XRAY (object));

  switch (prop_id)
  {
    case PROP_SIGNIFICANCE_OBS_MIN:
      g_value_set_double (value, mbrx->signif_obs_min);
      break;
	case PROP_SIGNIFICANCE_OBS_MAX:
      g_value_set_double (value, mbrx->signif_obs_max);
      break;
	case PROP_Z0:
      g_value_set_double (value, mbrx->z0);
      break;
	case PROP_M0:
      g_value_set_double (value, mbrx->M0);
      break;
	default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_mass_benson_xray_init (NcClusterMassBensonXRay *mbxr)
{
  mbxr->signif_obs_min = 0.0;
  mbxr->signif_obs_max = 0.0;
  mbxr->z0 = 0.0;
  mbxr->M0 = 0.0;
}

static void
_nc_cluster_mass_benson_xray_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_benson_xray_parent_class)->finalize (object);
}

static void
nc_cluster_mass_benson_xray_class_init (NcClusterMassBensonXRayClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterMassClass* parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  parent_class->obs_len = &_nc_cluster_mass_benson_xray_obs_len;
  parent_class->obs_params_len = &_nc_cluster_mass_benson_xray_obs_params_len;

  object_class->finalize = _nc_cluster_mass_benson_xray_finalize;
  object_class->set_property = &ncm_model_class_set_property;
  object_class->get_property = &ncm_model_class_get_property;


  model_class->set_property = &_nc_cluster_mass_benson_xray_set_property;
  model_class->get_property = &_nc_cluster_mass_benson_xray_get_property;

  ncm_model_class_add_params (model_class, 8, 0, PROP_SIZE);
  ncm_model_class_set_name_nick (model_class, "Benson- SZ and XRay", "Benson_SZ_XRay");

  /**
   * NcClusterMassBensonXRay:signif_obs_min:
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
   * NcClusterMassBensonXray:signif_obs_max:
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
   * NcClusterMassBensonXRay:z0:
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
   * NcClusterMassBensonXRay:M0:
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
  /*
   * SZ signal-mass scaling parameter: Asz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_A_SZ, "A_{SZ}", "Asz",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_A_SZ,
                               NCM_PARAM_TYPE_FREE);

  /*
   * SZ signal-mass scaling parameter: Bsz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_B_SZ, "B_{SZ}", "Bsz",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_B_SZ,
                               NCM_PARAM_TYPE_FIXED);

  /*
   * SZ signal-mass scaling parameter: Csz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_C_SZ, "C_{SZ}", "Csz",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_C_SZ,
                               NCM_PARAM_TYPE_FIXED);
  
 /*
   * SZ signal-mass scaling parameter: Dsz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_D_SZ, "D_{SZ}", "Dsz",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_D_SZ,
                               NCM_PARAM_TYPE_FIXED);
  /*
   * X-ray signal-mass scaling parameter: Ax.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_A_X, "A_{X}", "Ax",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_A_X,
                               NCM_PARAM_TYPE_FREE);

  /*
   * X-ray signal-mass scaling parameter: Bx.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_B_X, "B_{X}", "Bx",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_B_X,
                               NCM_PARAM_TYPE_FIXED);

  /*
   * X-ray signal-mass scaling parameter: Cx.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_C_X, "C_{X}", "Cx",
                               -2.0,  8.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_C_X,
                               NCM_PARAM_TYPE_FIXED);
  
 /*
   * X-ray signal-mass scaling parameter: Dx.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_D_X, "D_{X}", "Dx",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_D_X,
                               NCM_PARAM_TYPE_FIXED);
  
}

