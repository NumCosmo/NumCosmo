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
#include "nc_distance.h"
#include "lss/nc_halo_position.h"
#include "math/ncm_quaternion.h"


typedef struct _NcHaloPositionPrivate
{
  NcDistance *dist;
  NcmModelCtrl *ctrl_cosmo;
  NcmQuaternion q;
} NcHaloPositionPrivate;

struct _NcHaloPosition
{
  /*< private >*/
  NcmModel parent_instance;
};

enum
{
  PROP_0,
  PROP_DIST,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloPosition, nc_halo_position, NCM_TYPE_MODEL)

#define VECTOR (NCM_MODEL (hp))
#define RA     (ncm_model_orig_param_get (VECTOR, NC_HALO_POSITION_RA))
#define DEC    (ncm_model_orig_param_get (VECTOR, NC_HALO_POSITION_DEC))
#define Z      (ncm_model_orig_param_get (VECTOR, NC_HALO_POSITION_Z))

static void
nc_halo_position_init (NcHaloPosition *hp)
{
  NcHaloPositionPrivate *self = nc_halo_position_get_instance_private (hp);

  self->dist       = NULL;
  self->ctrl_cosmo = ncm_model_ctrl_new (NULL);
}

static void
_nc_halo_position_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloPosition *hp                 = NC_HALO_POSITION (object);
  NcHaloPositionPrivate * const self = nc_halo_position_get_instance_private (hp);

  g_return_if_fail (NC_IS_HALO_POSITION (hp));

  switch (prop_id)
  {
    case PROP_DIST:
      self->dist = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_position_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloPosition *hp                 = NC_HALO_POSITION (object);
  NcHaloPositionPrivate * const self = nc_halo_position_get_instance_private (hp);

  g_return_if_fail (NC_IS_HALO_POSITION (hp));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, self->dist);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_position_dispose (GObject *object)
{
  NcHaloPosition *hp                 = NC_HALO_POSITION (object);
  NcHaloPositionPrivate * const self = nc_halo_position_get_instance_private (hp);

  if (self->dist)
    nc_distance_clear (&self->dist);

  if (self->ctrl_cosmo)
    ncm_model_ctrl_clear (&self->ctrl_cosmo);

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

  model_class->set_property = &_nc_halo_position_set_property;
  model_class->get_property = &_nc_halo_position_get_property;
  object_class->dispose     = &_nc_halo_position_dispose;
  object_class->finalize    = &_nc_halo_position_finalize;

  ncm_model_class_set_name_nick (model_class, "Halo position", "Halo position");
  ncm_model_class_add_params (model_class, NC_HALO_POSITION_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcHaloPosition:ra:
   *
   * The right ascension of the halo position.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_POSITION_RA, "ra", "ra", 0.0, 2.0 * M_PI, 1.0e-2, NC_HALO_POSITION_DEFAULT_PARAMS_ABSTOL, NC_HALO_POSITION_DEFAULT_RA, NCM_PARAM_TYPE_FIXED);

  /**
   * NcHaloPosition:dec:
   *
   * The declination of the halo position.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_POSITION_DEC, "dec", "dec", -M_PI_2, M_PI_2, 1.0e-2, NC_HALO_POSITION_DEFAULT_PARAMS_ABSTOL, NC_HALO_POSITION_DEFAULT_DEC, NCM_PARAM_TYPE_FIXED);

  /**
   * NcHaloPosition:z:
   *
   * The redshift of the halo position.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_POSITION_Z, "z", "z", 0.0, 1100.0, 1.0e-3, NC_HALO_POSITION_DEFAULT_PARAMS_ABSTOL, NC_HALO_POSITION_DEFAULT_Z, NCM_PARAM_TYPE_FIXED);

  /**
   * NcHaloPosition:dist:
   *
   * A #NcDistance object.
   *
   */
  g_object_class_install_property (object_class, PROP_DIST,
                                   g_param_spec_object ("dist", "Distance",
                                                        "A NcDistance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY));

  ncm_mset_model_register_id (model_class, "NcHaloPosition", "Halo position.", NULL, TRUE, NCM_MSET_MODEL_MAIN);

  ncm_model_class_check_params_info (model_class);
}

/**
 * nc_halo_position_new:
 *
 * Creates a new halo position.
 *
 * Returns: (transfer full): A new halo position.
 *
 */
NcHaloPosition *
nc_halo_position_new (NcDistance *dist)
{
  return g_object_new (NC_TYPE_HALO_POSITION,
                       "dist", dist,
                       NULL);
}

/**
 * nc_halo_position_ref:
 * @hp: A #NcHaloPosition.
 *
 * Increases the reference count of @hp.
 *
 * Returns: (transfer full): The halo position.
 *
 */
NcHaloPosition *
nc_halo_position_ref (NcHaloPosition *hp)
{
  return g_object_ref (hp);
}

/**
 * nc_halo_position_free:
 * @hp: A #NcHaloPosition.
 *
 * Decreases the reference count of @hp.
 *
 */
void
nc_halo_position_free (NcHaloPosition *hp)
{
  g_object_unref (hp);
}

/**
 * nc_halo_position_clear:
 * @hp: A #NcHaloPosition.
 *
 * Clears the halo position.
 *
 */
void
nc_halo_position_clear (NcHaloPosition **hp)
{
  g_clear_object (hp);
}

/**
 * nc_halo_position_prepare:
 * @hp: A #NcHaloPosition.
 * @cosmo: A #NcHICosmo.
 *
 * Prepares the #NcHaloPosition object @hp for computation.
 *
 */
void
nc_halo_position_prepare (NcHaloPosition *hp, NcHICosmo *cosmo)
{
  NcHaloPositionPrivate *self = nc_halo_position_get_instance_private (hp);

  nc_distance_prepare_if_needed (self->dist, cosmo);
}

/**
 * nc_halo_position_prepare_if_needed:
 * @hp: A #NcHaloPosition.
 * @cosmo: A #NcHICosmo.
 *
 * Prepares the #NcHaloPosition object @hp for computation if necessary.
 *
 */
void
nc_halo_position_prepare_if_needed (NcHaloPosition *hp, NcHICosmo *cosmo)
{
  NcHaloPositionPrivate *self = nc_halo_position_get_instance_private (hp);

  if (ncm_model_ctrl_update (self->ctrl_cosmo, NCM_MODEL (cosmo)))
    nc_distance_prepare_if_needed (self->dist, cosmo);
}

/**
 * nc_halo_position_get_redshift:
 * @hp: A #NcHaloPosition
 *
 * Returns the redshift of the halo position.
 */
gdouble
nc_halo_position_get_redshift (NcHaloPosition *hp)
{
  return Z;
}

/**
 * nc_halo_position_get_ra_dec:
 * @hp: A #NcHaloPosition
 * @ra: (out): The right ascension
 * @dec: (out): The declination
 *
 * Returns the right ascension and declination of the halo position.
 *
 */
void
nc_halo_position_get_ra_dec (NcHaloPosition *hp, gdouble *ra, gdouble *dec)
{
  *ra  = RA;
  *dec = DEC;
}

/**
 * nc_halo_position_polar_angles:
 * @hp: A #NcHaloPosition.
 * @ra: The right ascension.
 * @dec: The declination.
 * @theta: (out): The polar angle.
 * @phi: (out): The azimuthal angle.
 *
 * Considering a spherical coordinate system where the halo is at the north pole, this
 * function calculates the polar and azimuthal angles of a point in the sky. That is,
 * $\theta$ = @theta provides the angular separation from the halo to the point in the
 * sky, and $\phi$ = @phi provides the azimuthal angle.
 *
 */
void
nc_halo_position_polar_angles (NcHaloPosition *hp, gdouble ra, gdouble dec, gdouble *theta, gdouble *phi)
{
  NcHaloPositionPrivate *self = nc_halo_position_get_instance_private (hp);
  NcmModel *model             = NCM_MODEL (hp);
  NcmTriVec v;
  gdouble r;

  if (!ncm_model_state_is_update (model))
  {
    const gdouble ra_halo  = RA;
    const gdouble dec_halo = DEC;

    ncm_trivec_set_astro_ra_dec (&v, 1.0, ra_halo, dec_halo);
    ncm_quaternion_set_to_rotate_to_z (&self->q, &v);

    ncm_model_state_set_update (model);
  }

  ncm_trivec_set_astro_ra_dec (&v, 1.0, ra, dec);
  ncm_quaternion_rotate (&self->q, &v);
  ncm_trivec_get_spherical_coord (&v, &r, theta, phi);
}

/**
 * nc_halo_position_projected_radius:
 * @hp: A #NcHaloPosition
 * @cosmo: A #NcHICosmo
 * @theta: The angular separation
 *
 * Calculates the projected radius of the halo position.
 *
 * Returns: The projected radius in Mpc.
 */
gdouble
nc_halo_position_projected_radius (NcHaloPosition *hp, NcHICosmo *cosmo, gdouble theta)
{
  NcHaloPositionPrivate *self = nc_halo_position_get_instance_private (hp);
  const gdouble RH_Mpc        = nc_distance_hubble (self->dist, cosmo);
  const gdouble z_halo        = Z;
  const gdouble dA            = nc_distance_angular_diameter (self->dist, cosmo, z_halo);

  return ncm_util_projected_radius (theta, dA * RH_Mpc);
}

/**
 * nc_halo_position_projected_radius_from_ra_dec:
 * @hp: A #NcHaloPosition
 * @cosmo: A #NcHICosmo
 * @ra: The right ascension
 * @dec: The declination
 *
 * Calculates the projected radius of the halo position from the right ascension and declination.
 *
 * Returns: The projected radius in Mpc.
 */
gdouble
nc_halo_position_projected_radius_from_ra_dec (NcHaloPosition *hp, NcHICosmo *cosmo, gdouble ra, gdouble dec)
{
  gdouble theta, phi;

  nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);

  return nc_halo_position_projected_radius (hp, cosmo, theta);
}

