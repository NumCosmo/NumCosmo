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
  PROP_AX,
  PROP_BX,
  PROP_CX,
  PROP_DX,
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
	case PROP_ASZ:
	  mbxr->Asz = g_value_get_double (value);
	  break;
	case PROP_BSZ:
	  mbxr->Bsz = g_value_get_double (value);
	  break;
	case PROP_CSZ:
	  mbxr->Csz = g_value_get_double (value);
	  break;
	case PROP_DSZ:
	  mbxr->Dsz = g_value_get_double (value);
	  break;
	case PROP_AX:
	  mbxr->Ax = g_value_get_double (value);
	  break;
	case PROP_BX:
	  mbxr->Bx = g_value_get_double (value);
	  break;
	case PROP_CX:
	  mbxr->Cx = g_value_get_double (value);
	  break;
	case PROP_DX:
	  mbxr->Dx = g_value_get_double (value);
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
	case PROP_ASZ:
      g_value_set_double (value, mbrx->Asz);
      break;
	case PROP_BSZ:
      g_value_set_double (value, mbrx->Bsz);
      break;
	case PROP_CSZ:
      g_value_set_double (value, mbrx->Csz);
      break;
	case PROP_DSZ:
      g_value_set_double (value, mbrx->Dsz);
      break;
	case PROP_AX:
      g_value_set_double (value, mbrx->Ax);
      break;
	case PROP_BX:
      g_value_set_double (value, mbrx->Bx);
      break;
	case PROP_CX:
      g_value_set_double (value, mbrx->Cx);
      break;
	case PROP_DX:
      g_value_set_double (value, mbrx->Dx);
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
  mbxr->Asz = 0.0;
  mbxr->Bsz = 0.0;
  mbxr->Csz = 0.0;
  mbxr->Dsz = 0.0;
  mbxr->Ax = 0.0;
  mbxr->Bx = 0.0;
  mbxr->Cx = 0.0;
  mbxr->Dx = 0.0;
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

  parent_class->obs_len = &_nc_cluster_mass_benson_xray_obs_len;
  parent_class->obs_params_len = &_nc_cluster_mass_benson_xray_obs_params_len;

  object_class->finalize = _nc_cluster_mass_benson_xray_finalize;
  object_class->set_property = &_nc_cluster_mass_benson_xray_set_property;
  object_class->get_property = &_nc_cluster_mass_benson_xray_get_property;
}

