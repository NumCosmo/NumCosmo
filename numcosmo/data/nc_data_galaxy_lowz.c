/***************************************************************************
 *            nc_data_galaxy_lowz.c
 *
 *  Fri Jan 26 16:27:32 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:nc_data_galaxy_lowz
 * @title: NcDataGalaxyLowz
 * @short_description: Distance modulus data for low-z galaxies.
 *
 * Distance modulus for low-z galaxies including directions.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_galaxy_lowz.h"
#include "model/nc_acosmo_lowz.h"
#include "math/ncm_data_gauss_diag.h"
#include "math/ncm_cfg.h"
#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_POS_Z_THETA_PHI,
  PROP_SIZE,
};

struct _NcDataGalaxyLowz
{
  /*< private >*/
  NcmDataGaussDiag parent_instance;
  NcmMatrix *position_z_theta_phi;
};


G_DEFINE_TYPE (NcDataGalaxyLowz, nc_data_galaxy_lowz, NCM_TYPE_DATA_GAUSS_DIAG)

static void
nc_data_galaxy_lowz_init (NcDataGalaxyLowz *galaxy_lowz)
{
  galaxy_lowz->position_z_theta_phi = NULL;
}

static void
_nc_data_galaxy_lowz_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_galaxy_lowz_parent_class)->constructed (object);
}

static void
nc_data_galaxy_lowz_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataGalaxyLowz *galaxy_lowz = NC_DATA_GALAXY_LOWZ (object);

  g_return_if_fail (NC_IS_DATA_GALAXY_LOWZ (object));

  switch (prop_id)
  {
    case PROP_POS_Z_THETA_PHI:
      ncm_matrix_substitute (&galaxy_lowz->position_z_theta_phi, g_value_get_object (value), TRUE);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_galaxy_lowz_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataGalaxyLowz *galaxy_lowz = NC_DATA_GALAXY_LOWZ (object);

  g_return_if_fail (NC_IS_DATA_GALAXY_LOWZ (object));

  switch (prop_id)
  {
    case PROP_POS_Z_THETA_PHI:
      g_value_set_object (value, galaxy_lowz->position_z_theta_phi);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_galaxy_lowz_dispose (GObject *object)
{
  NcDataGalaxyLowz *galaxy_lowz = NC_DATA_GALAXY_LOWZ (object);

  ncm_matrix_clear (&galaxy_lowz->position_z_theta_phi);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_galaxy_lowz_parent_class)->dispose (object);
}

static void
nc_data_galaxy_lowz_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_galaxy_lowz_parent_class)->finalize (object);
}

static void _nc_data_galaxy_lowz_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_galaxy_lowz_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);
static void _nc_data_galaxy_lowz_set_size (NcmDataGaussDiag *diag, guint np);

static void
nc_data_galaxy_lowz_class_init (NcDataGalaxyLowzClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass *diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->constructed  = &_nc_data_galaxy_lowz_constructed;
  object_class->set_property = &nc_data_galaxy_lowz_set_property;
  object_class->get_property = &nc_data_galaxy_lowz_get_property;
  object_class->dispose      = &nc_data_galaxy_lowz_dispose;
  object_class->finalize     = &nc_data_galaxy_lowz_finalize;

  g_object_class_install_property (object_class,
                                   PROP_POS_Z_THETA_PHI,
                                   g_param_spec_object ("pos-z-theta-phi",
                                                        NULL,
                                                        "Position matrix (z, theta, phi)",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  data_class->prepare   = &_nc_data_galaxy_lowz_prepare;
  diag_class->mean_func = &_nc_data_galaxy_lowz_mean_func;
  diag_class->set_size  = &_nc_data_galaxy_lowz_set_size;
}

static void
_nc_data_galaxy_lowz_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataGalaxyLowz *galaxy_lowz = NC_DATA_GALAXY_LOWZ (data);
  NcACosmoLowz *aclz            = NC_ACOSMO_LOWZ (ncm_mset_peek (mset, nc_acosmo_lowz_id ()));

  g_assert_nonnull (aclz);
}

static void
_nc_data_galaxy_lowz_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataGalaxyLowz *galaxy_lowz = NC_DATA_GALAXY_LOWZ (diag);
  NcACosmoLowz *aclz            = NC_ACOSMO_LOWZ (ncm_mset_peek (mset, nc_acosmo_lowz_id ()));
  const guint np                = ncm_data_gauss_diag_get_size (diag);
  guint i;

  for (i = 0; i < np; i++)
  {
    const gdouble z     = ncm_matrix_get (galaxy_lowz->position_z_theta_phi, i, 0);
    const gdouble theta = ncm_matrix_get (galaxy_lowz->position_z_theta_phi, i, 1);
    const gdouble phi   = ncm_matrix_get (galaxy_lowz->position_z_theta_phi, i, 2);
    const gdouble mu    = nc_acosmo_lowz_distance_modulus (aclz, z, theta, phi);

    ncm_vector_set (vp, i, mu);
  }
}

static void
_nc_data_galaxy_lowz_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataGalaxyLowz *galaxy_lowz = NC_DATA_GALAXY_LOWZ (diag);
  const guint cnp               = ncm_data_gauss_diag_get_size (diag);

  if ((np == 0) || (np != cnp))
    ncm_matrix_clear (&galaxy_lowz->position_z_theta_phi);

  if ((np != 0) && (np != cnp))
    galaxy_lowz->position_z_theta_phi = ncm_matrix_new (np, 3);

  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_galaxy_lowz_parent_class)->set_size (diag, np);
}

/**
 * nc_data_galaxy_lowz_new_empty:
 *
 * Creates an empty #NcDataGalaxyLowz.
 *
 * Returns: the newly created #NcDataGalaxyLowz.
 */
NcDataGalaxyLowz *
nc_data_galaxy_lowz_new_empty (void)
{
  NcDataGalaxyLowz *galaxy_lowz = g_object_new (NC_TYPE_DATA_GALAXY_LOWZ,
                                                NULL);

  return galaxy_lowz;
}

/**
 * nc_data_galaxy_lowz_peek_position_z_theta_phi:
 * @galaxy_lowz: a #NcDataGalaxyLowz
 *
 * Returns: (transfer none): the position matrix $(z, \theta, \phi)$.
 */
NcmMatrix *
nc_data_galaxy_lowz_peek_position_z_theta_phi (NcDataGalaxyLowz *galaxy_lowz)
{
  g_return_val_if_fail (NC_IS_DATA_GALAXY_LOWZ (galaxy_lowz), NULL);

  return galaxy_lowz->position_z_theta_phi;
}

