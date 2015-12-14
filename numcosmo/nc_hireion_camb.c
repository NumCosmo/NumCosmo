/***************************************************************************
 *            nc_hireion_camb.c
 *
 *  Thu December 10 11:56:38 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>>
 ****************************************************************************/
/*
 * nc_hireion_camb.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hireion_camb
 * @title: NcHIReionCamb
 * @short_description: CAMB-like reionization object.
 *
 * This object implements the reionization as done in CAMB. 
 * 
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hireion_camb.h"

enum
{
  PROP_0,
  PROP_HII_HEII_REION_DELTA,
  PROP_HEIII_REION_DELTA,
  PROP_HII_HEII_REION_EXPO,
  PROP_HEIII_REIONIZED,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcHIReionCamb, nc_hireion_camb, NC_TYPE_HIREION);

#define VECTOR     (NCM_MODEL (reion)->params)
#define HII_HEII_Z (ncm_vector_get (VECTOR, NC_HIREION_CAMB_HII_HEII_Z))
#define HEIII_Z    (ncm_vector_get (VECTOR, NC_HICOSMO_CAMB_HEIII_Z))

static void
nc_hireion_camb_init (NcHIReionCamb *reion_camb)
{
  reion_camb->HII_HeII_reion_delta      = 0.0;
  reion_camb->HeIII_reion_delta         = 0.0;
  reion_camb->HII_HeII_reion_expo       = 0.0;
  reion_camb->HII_HeII_reion_delta_eff  = 0.0;
  reion_camb->HII_HeII_reion_x_pow_expo = 0.0;
  reion_camb->HEII_reionized            = FALSE;
}

static void
nc_hireion_camb_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hireion_camb_parent_class)->finalize (object);
}

static void
nc_hireion_camb_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (object);
  g_return_if_fail (NC_IS_HIREION_CAMB (object));

  switch (prop_id)
  {
    case PROP_HII_HEII_REION_DELTA:
      reion_camb->HII_HeII_reion_delta = g_value_get_double (value);
      break;
    case PROP_HEIII_REION_DELTA:
      reion_camb->HeIII_reion_delta = g_value_get_double (value);
      break;
    case PROP_HII_HEII_REION_EXPO:
      reion_camb->HII_HeII_reion_expo = g_value_get_double (value);
      break;
    case PROP_HEIII_REIONIZED:
      reion_camb->HEII_reionized = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hireion_camb_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (object);
  g_return_if_fail (NC_IS_HIREION_CAMB (object));

  switch (prop_id)
  {
    case PROP_HII_HEII_REION_DELTA:
      g_value_set_double (value, reion_camb->HII_HeII_reion_delta);
      break;
    case PROP_HEIII_REION_DELTA:
      g_value_set_double (value, reion_camb->HeIII_reion_delta);
      break;
    case PROP_HII_HEII_REION_EXPO:
      g_value_set_double (value, reion_camb->HII_HeII_reion_expo);
      break;
    case PROP_HEIII_REIONIZED:
      g_value_set_boolean (value, reion_camb->HEII_reionized);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _nc_hireion_camb_update_Xe (NcHIReion *reion, NcRecomb *recomb);

static void
nc_hireion_camb_class_init (NcHIReionCambClass *klass)
{
  GObjectClass* object_class  = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class  = NCM_MODEL_CLASS (klass);
  NcHIReionClass *reion_class = NC_HIREION_CLASS (klass);

  object_class->finalize    = nc_hireion_camb_finalize;

  model_class->set_property = &nc_hireion_camb_set_property;
  model_class->get_property = &nc_hireion_camb_get_property;
  
  ncm_model_class_set_name_nick (model_class, "Reion-CAMB", "REION_CAMB");
  ncm_model_class_add_params (model_class, NC_HIREION_CAMB_SPARAM_LEN, 0, PROP_SIZE);

  /* Set HII_HEII_Z param info */
  ncm_model_class_set_sparam (model_class, NC_HIREION_CAMB_HII_HEII_Z, "z_\\mathrm{re}", "z_re",
                               0.0, 50.0, 0.1,
                               NC_HIREION_DEFAULT_PARAMS_ABSTOL, NC_HIREION_CAMB_DEFAULT_HII_HEII_Z,
                               NCM_PARAM_TYPE_FIXED);
  /* Set HEIII_Z param info */
  ncm_model_class_set_sparam (model_class, NC_HIREION_CAMB_HEIII_Z, "z^\\mathrm{He}_\\mathrm{re}", "z_He_re",
                               0.0,  10.0, 1.0e-2,
                               NC_HIREION_DEFAULT_PARAMS_ABSTOL, NC_HIREION_CAMB_DEFAULT_HEIII_Z,
                               NCM_PARAM_TYPE_FIXED);

  g_object_class_install_property (object_class,
                                   PROP_HII_HEII_REION_DELTA,
                                   g_param_spec_double ("HII-HeII-reion-delta",
                                                        NULL,
                                                        "Window size for HII and HeII reionization",
                                                        0.0, G_MAXDOUBLE, NC_HIREION_CAMB_DEFAULT_HII_HEII_REION_DELTA,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_HEIII_REION_DELTA,
                                   g_param_spec_double ("HeIII-reion-delta",
                                                        NULL,
                                                        "Window size for HeIII reionization",
                                                        0.0, G_MAXDOUBLE, NC_HIREION_CAMB_DEFAULT_HEIII_REION_DELTA,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_HII_HEII_REION_EXPO,
                                   g_param_spec_double ("HII-HeII-reion-exponent",
                                                        NULL,
                                                        "Exponent for HII and HeII reionization transition",
                                                        0.0, G_MAXDOUBLE, NC_HIREION_CAMB_DEFAULT_HII_HEII_REION_EXPO,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HEIII_REIONIZED,
                                   g_param_spec_boolean ("HeII-reionized",
                                                         NULL,
                                                         "Whether HeIII is reionized",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  

  reion_class->update_Xe = &_nc_hireion_camb_update_Xe;
}



static void 
_nc_hireion_camb_update_Xe (NcHIReion *reion, NcRecomb *recomb)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (reion);
  if (!ncm_model_state_is_update (NCM_MODEL (reion)))
  {
    const gdouble xre = 1.0 + HII_HEII_Z;

    reion_camb->HII_HeII_reion_x_pow_expo = pow (xre, reion_camb->HII_HeII_reion_expo);
    reion_camb->HII_HeII_reion_delta_eff  = 
      reion_camb->HII_HeII_reion_expo * 
      reion_camb->HII_HeII_reion_x_pow_expo * 
      reion_camb->HII_HeII_reion_delta / xre;
    ncm_model_state_set_update (NCM_MODEL (reion));
  }

  
}

NcHIReionCamb *
nc_hireion_camb_new (void)
{
  NcHIReionCamb *reion_camb = g_object_new (NC_TYPE_HIREION_CAMB, 
                                            NULL);

  return reion_camb;
}

