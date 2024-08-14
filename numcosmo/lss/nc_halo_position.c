/***************************************************************************
 *            nc_halo_position.c
 *
 *  Thu Aug 08 08:45:20 2024
 *  Copyright  2024
 *  <code.caio@limadeoliveira.me>
 ****************************************************************************/
/*
 * nc_halo_position.c
 * Copyright (C) 2024 Caio Lima de Oliveira <code.caio@limadeoliveira.me>
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
 * SECTION:nc_halo_position
 * @title: NcHaloPosition
 * @short_description: A class to represent the center of a halo.
 * @stability: Unstable
 *
 * This class describes the center of a halo.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "lss/nc_halo_position.h"

typedef struct _NcHaloPositionPrivate
{
  gdouble ra;
  gdouble dec;
  gdouble z;
} NcHaloPositionPrivate;

struct _NcHaloPosition
{
  NcmModel parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloPosition, nc_halo_position, NCM_TYPE_MODEL)

#define VECTOR (NCM_MODEL (hc))
#define RA     (ncm_model_orig_param_get (VECTOR, NC_HALO_POSITION_RA))
#define DEC    (ncm_model_orig_param_get (VECTOR, NC_HALO_POSITION_DEC))
#define Z      (ncm_model_orig_param_get (VECTOR, NC_HALO_POSITION_Z))

static void
nc_halo_position_init (NcHaloPosition *hc)
{
}

static void
_nc_halo_position_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloPosition *hc = NC_HALO_POSITION (object);

  g_return_if_fail (NC_IS_HALO_POSITION (hc));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_position_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloPosition *hc = NC_HALO_POSITION (object);

  g_return_if_fail (NC_IS_HALO_POSITION (hc));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_position_dispose (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_position_parent_class)->dispose (object);
}

static void
_nc_halo_position_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_position_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_halo_position, NC_TYPE_HALO_POSITION)

static void
nc_halo_position_class_init (NcHaloPositionClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->set_property = &_nc_halo_position_set_property;
  object_class->get_property = &_nc_halo_position_get_property;
  object_class->dispose      = &_nc_halo_position_dispose;
  object_class->finalize     = &_nc_halo_position_finalize;

  ncm_model_class_set_name_nick (model_class, "Halo center", "Halo center");
  ncm_model_class_add_params (model_class, NC_HALO_POSITION_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcHaloPosition:ra:
   *
   * The right ascension of the halo center.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_POSITION_RA, "ra", "ra", 0.0, 2.0 * M_PI, 1.0e-2, NC_HALO_POSITION_DEFAULT_PARAMS_ABSTOL, NC_HALO_POSITION_DEFAULT_RA, NCM_PARAM_TYPE_FIXED);

  /**
   * NcHaloPosition:dec:
   *
   * The declination of the halo center.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_POSITION_DEC, "dec", "dec", -M_PI_2, M_PI_2, 1.0e-2, NC_HALO_POSITION_DEFAULT_PARAMS_ABSTOL, NC_HALO_POSITION_DEFAULT_DEC, NCM_PARAM_TYPE_FIXED);

  /**
   * NcHaloPosition:z:
   *
   * The redshift of the halo center.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_POSITION_Z, "z", "z", 0.0, 1100.0, 1.0e-3, NC_HALO_POSITION_DEFAULT_PARAMS_ABSTOL, NC_HALO_POSITION_DEFAULT_Z, NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);
}

/**
 * nc_halo_position_new:
 *
 * Creates a new halo center.
 *
 * Returns: (transfer full): A new halo center.
 *
 */
NcHaloPosition *
nc_halo_position_new (void)
{
  return g_object_new (NC_TYPE_HALO_POSITION, NULL);
}

/**
 * nc_halo_position_ref:
 * @hc: A halo center.
 *
 * Increases the reference count of @hc.
 *
 * Returns: (transfer full): The halo center.
 *
 */
NcHaloPosition *
nc_halo_position_ref (NcHaloPosition *hc)
{
  return g_object_ref (hc);
}

/**
 * nc_halo_position_free:
 * @hc: A halo center.
 *
 * Decreases the reference count of @hc.
 *
 */
void
nc_halo_position_free (NcHaloPosition *hc)
{
  g_object_unref (hc);
}

/**
 * nc_halo_position_clear:
 * @hc: A halo center.
 *
 * Clears the halo center.
 *
 */
void
nc_halo_position_clear (NcHaloPosition **hc)
{
  g_clear_object (hc);
}

