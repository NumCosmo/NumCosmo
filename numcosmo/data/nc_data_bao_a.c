/***************************************************************************
 *            nc_data_bao_a.c
 *
 *  Thu Apr 22 15:31:32 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcDataBaoA:
 *
 * Baryon oscillation data -- acoustic scale $A$.
 *
 * The acoustic scale is defined as $$A \equiv D_V (z) \frac{\sqrt{\Omega_{m0}
 * H_0^2}}{z c},$$ where $\Omega_{m0}$ is the matter density parameter
 * [nc_hicosmo_Omega_m0()], $c$ is the speed of light [ncm_c_c()], $H_0$ is the Hubble
 * parameter [nc_hicosmo_H0()] and $D_V(z)$ is the dilation scale
 * [nc_distance_dilation_scale()]. See Section 4.5 from [Eisenstein et al.
 * (2005)][XEisenstein2005].
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_a.h"
#include "nc_enum_types.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoA, nc_data_bao_a, NCM_TYPE_DATA_GAUSS_DIAG)

static void
nc_data_bao_a_init (NcDataBaoA *bao_a)
{
  bao_a->dist = nc_distance_new (2.0);
  bao_a->x    = NULL;
}

static void
nc_data_bao_a_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (object);

  g_return_if_fail (NC_IS_DATA_BAO_A (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_bao_a_set_dist (bao_a, g_value_get_object (value));
      break;
    case PROP_Z:
      ncm_vector_substitute (&bao_a->x, g_value_get_object (value), TRUE);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_a_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (object);

  g_return_if_fail (NC_IS_DATA_BAO_A (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, bao_a->dist);
      break;
    case PROP_Z:
      g_value_set_object (value, bao_a->x);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_a_dispose (GObject *object)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (object);

  nc_distance_clear (&bao_a->dist);
  ncm_vector_clear (&bao_a->x);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_a_parent_class)->dispose (object);
}

static void
nc_data_bao_a_finalize (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_bao_a_parent_class)->finalize (object);
}

static void _nc_data_bao_a_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_bao_a_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);
static void _nc_data_bao_a_set_size (NcmDataGaussDiag *diag, guint np);

static void
nc_data_bao_a_class_init (NcDataBaoAClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass *diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->set_property = &nc_data_bao_a_set_property;
  object_class->get_property = &nc_data_bao_a_get_property;
  object_class->dispose      = &nc_data_bao_a_dispose;
  object_class->finalize     = &nc_data_bao_a_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_object ("z",
                                                        NULL,
                                                        "Data redshift",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare   = &_nc_data_bao_a_prepare;
  diag_class->mean_func = &_nc_data_bao_a_mean_func;
  diag_class->set_size  = &_nc_data_bao_a_set_size;
}

static void
_nc_data_bao_a_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (data);
  NcHICosmo *cosmo  = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  nc_distance_prepare_if_needed (bao_a->dist, cosmo);
}

static void
_nc_data_bao_a_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (diag);
  NcHICosmo *cosmo  = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const guint np    = ncm_data_gauss_diag_get_size (diag);
  guint i;

  for (i = 0; i < np; i++)
  {
    const gdouble z = ncm_vector_get (bao_a->x, i);
    const gdouble A = nc_distance_bao_A_scale (bao_a->dist, cosmo, z);

    ncm_vector_set (vp, i, A);
  }
}

/**
 * nc_data_bao_a_new_from_file:
 * @filename: file containing a serialized #NcDataBaoA.
 *
 * Creates a new #NcDataBaoA from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataBaoA.
 */
NcDataBaoA *
nc_data_bao_a_new_from_file (const gchar *filename)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_BAO_A (bao_a));

  return bao_a;
}

/**
 * nc_data_bao_a_new:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * Creates a new acustic scale data object #NcDataBaoA from @id.
 * This object requires a #NcDistance object to be set.
 *
 * Returns: a #NcDataBaoA
 */
NcDataBaoA *
nc_data_bao_a_new_from_id (NcDistance *dist, NcDataBaoId id)
{
  NcDataBaoA *bao_a;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_BAO_A_EISENSTEIN2005:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_a_eisenstein2005.obj", TRUE);
      break;
    default:
      g_error ("nc_data_bao_a_new_from_id: id %d not recognized.", id);
      break;
  }

  bao_a = nc_data_bao_a_new_from_file (filename);
  nc_data_bao_a_set_dist (bao_a, dist);
  g_free (filename);

  return bao_a;
}

static void
_nc_data_bao_a_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (diag);
  const guint cnp   = ncm_data_gauss_diag_get_size (diag);

  if ((np == 0) || (np != cnp))
    ncm_vector_clear (&bao_a->x);

  if ((np != 0) && (np != cnp))
    bao_a->x = ncm_vector_new (np);

  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_bao_a_parent_class)->set_size (diag, np);
}

/**
 * nc_data_bao_a_set_dist:
 * @bao_a: a #NcDataBaoA
 * @dist: a #NcDistance
 *
 * Sets the distance object.
 *
 */
void
nc_data_bao_a_set_dist (NcDataBaoA *bao_a, NcDistance *dist)
{
  nc_distance_clear (&bao_a->dist);
  bao_a->dist = nc_distance_ref (dist);
}

