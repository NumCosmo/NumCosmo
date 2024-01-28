/***************************************************************************
 *            nc_acosmo_lowz.c
 *
 *  Fri Jan 26 15:37:04 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2024 <vitenti@uel.br>
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
 * SECTION:nc_acosmo_lowz
 * @title: NcACosmoLowz
 * @short_description: Anisotropic model for low-z computations
 *
 * Cosmological model for low-z computations assuming an anisotropic
 * but homogeneos cosmology.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_acosmo_lowz.h"
#include "math/ncm_mset.h"
#include "math/ncm_c.h"

struct _NcACosmoLowz
{
  /*< private >*/
  NcmModel parent_instance;
  gint place_holder;
};

G_DEFINE_TYPE (NcACosmoLowz, nc_acosmo_lowz, NCM_TYPE_MODEL)

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_acosmo_lowz_init (NcACosmoLowz *acosmo)
{
  /* Nothing to do */
  acosmo->place_holder = 0;
}

static void
_nc_acosmo_lowz_constructed (GObject *object)
{
  {
    NcmModel *model       = NCM_MODEL (object);
    const guint accel_len = ncm_model_vparam_len (model, NC_ACOSMO_LOWZ_ACCEL);
    const guint shear_len = ncm_model_vparam_len (model, NC_ACOSMO_LOWZ_SHEAR);

    if (accel_len != 3)
      g_error ("_nc_acosmo_lowz_constructed: acceleration vector must have three components.");

    if (shear_len != 5)
      g_error ("_nc_acosmo_lowz_constructed: shear vector must have five components.");
  }
  /* Chain up : start */
  G_OBJECT_CLASS (nc_acosmo_lowz_parent_class)->constructed (object);
}

static void
_nc_acosmo_lowz_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_acosmo_lowz_parent_class)->dispose (object);
}

static void
_nc_acosmo_lowz_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_acosmo_lowz_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_acosmo_lowz, NC_TYPE_ACOSMO_LOWZ);

static void
nc_acosmo_lowz_class_init (NcACosmoLowzClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->constructed = &_nc_acosmo_lowz_constructed;
  object_class->dispose     = &_nc_acosmo_lowz_dispose;
  object_class->finalize    = &_nc_acosmo_lowz_finalize;

  ncm_model_class_set_name_nick (model_class, "Anisotropic model for low-z", "NcACosmoLowz");

  ncm_mset_model_register_id (model_class,
                              "NcACosmoLowz",
                              "Anisotropic cosmological model for low-z",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);


  ncm_model_class_add_params (model_class,
                              NC_ACOSMO_LOWZ_SPARAM_LEN, NC_ACOSMO_LOWZ_VPARAM_LEN, PROP_SIZE);

  ncm_model_class_set_sparam (model_class, NC_ACOSMO_LOWZ_H0, "H_0", "H0",
                              31.0, 99.0, 1.0,
                              0.0, 70.0, NCM_PARAM_TYPE_FREE);

  /* Set massive neutrinos mass vector param */
  ncm_model_class_set_vparam (model_class, NC_ACOSMO_LOWZ_ACCEL, 3, "a", "a",
                              -0.4, 0.4, 0.01,
                              0.0, 0.0,
                              NCM_PARAM_TYPE_FREE);

  /* Set massive neutrinos temperature vector param */
  ncm_model_class_set_vparam (model_class, NC_ACOSMO_LOWZ_SHEAR, 5, "\\sigma", "sigma",
                              -0.4, 0.4, 0.01, 0.0, 0.0,
                              NCM_PARAM_TYPE_FREE);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
}

#define VECTOR   (NCM_MODEL (aclz))
#define MACRO_H0 (ncm_model_orig_param_get (VECTOR, NC_ACOSMO_LOWZ_H0))

/**
 * nc_acosmo_lowz_new:
 *
 * Creates a new #NcACosmoLowz.
 *
 * Returns: the newly created #NcACosmoLowz.
 */
NcACosmoLowz *
nc_acosmo_lowz_new (void)
{
  NcACosmoLowz *aclz = g_object_new (NC_TYPE_ACOSMO_LOWZ, NULL);

  return aclz;
}

/**
 * nc_acosmo_lowz_distance_modulus:
 * @aclz: a #NcACosmoLowz
 * @z: redshift $z$
 * @theta: $\theta$ angle
 * @phi: $\phi$ angle
 *
 *
 * Computes the distance modulus $\mu(\theta, \phi)$ to redshift $z$ and
 * and direction $(\theta, \phi)$.
 *
 * Returns: distance modulus $\mu$.
 */
gdouble
nc_acosmo_lowz_distance_modulus (NcACosmoLowz *aclz, const gdouble z, const gdouble theta, const gdouble phi)
{
  const gdouble H0            = MACRO_H0;
  const gdouble RH_Mpc        = ncm_c_c () / (1.0e3 * H0);
  const gdouble a_x           = ncm_model_orig_vparam_get (NCM_MODEL (aclz), NC_ACOSMO_LOWZ_ACCEL, 0);
  const gdouble a_y           = ncm_model_orig_vparam_get (NCM_MODEL (aclz), NC_ACOSMO_LOWZ_ACCEL, 1);
  const gdouble a_z           = ncm_model_orig_vparam_get (NCM_MODEL (aclz), NC_ACOSMO_LOWZ_ACCEL, 2);
  const gdouble shear_11      = ncm_model_orig_vparam_get (NCM_MODEL (aclz), NC_ACOSMO_LOWZ_SHEAR, 0);
  const gdouble shear_12      = ncm_model_orig_vparam_get (NCM_MODEL (aclz), NC_ACOSMO_LOWZ_SHEAR, 1);
  const gdouble shear_13      = ncm_model_orig_vparam_get (NCM_MODEL (aclz), NC_ACOSMO_LOWZ_SHEAR, 2);
  const gdouble shear_22      = ncm_model_orig_vparam_get (NCM_MODEL (aclz), NC_ACOSMO_LOWZ_SHEAR, 3);
  const gdouble shear_23      = ncm_model_orig_vparam_get (NCM_MODEL (aclz), NC_ACOSMO_LOWZ_SHEAR, 4);
  const gdouble shear_33      = -shear_11 - shear_22;
  const gdouble n_x           = sin (theta) * cos (phi);
  const gdouble n_y           = sin (theta) * sin (phi);
  const gdouble n_z           = cos (theta);
  const gdouble n_dot_a       = n_x * a_x + n_y * a_y + n_z * a_z;
  const gdouble shear_dot_n_n = (shear_11 * n_x * n_x + shear_22 * n_y * n_y + shear_33 * n_z * n_z
                                 + 2.0 * (shear_12 * n_x * n_y + shear_13 * n_x * n_z + shear_23 * n_y * n_z));
  const gdouble D_Mpc = RH_Mpc * z / (1.0 - n_dot_a + shear_dot_n_n);

  /* printf ("% 22.15g % 22.15g % 22.15g\n", n_x, n_y, n_z); */

  return 5.0 * log10 (D_Mpc) + 25.0;
}

