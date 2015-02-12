/***************************************************************************
 *            nc_data_bao_a.c
 *
 *  Thu Apr 22 15:31:32 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_data_bao_a
 * @title: Baryonic Oscillation Data -- Acoustic Scale
 * @short_description: BAO acoustic scale estimator
 *
 * The acoustic scale is defined as 
 * $$ A \equiv D_V (z) \frac{\sqrt{\Omega_m H_0^2}}{z c},$$
 * where $\Omega_m$ is the matter density parameter [nc_hicosmo_Omega_m()], $c$ is the speed of light [ncm_c_c()], 
 * $H_0$ is the Hubble parameter [nc_hicosmo_H0()] and $D_V(z)$ is the dilation scale [nc_distance_dilation_scale()].
 * See Section 4.5 from <link linkend="XEisenstein2005">Eisenstein et al. (2005)</link>.
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

G_DEFINE_TYPE (NcDataBaoA, nc_data_bao_a, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_bao_a_init (NcDataBaoA *bao_a)
{
  bao_a->dist = NULL;
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
      nc_distance_clear (&bao_a->dist);
      bao_a->dist = g_value_dup_object (value);
      break;
    case PROP_Z:
      ncm_vector_set_from_variant (bao_a->x, g_value_get_variant (value));
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
      g_value_take_variant (value, ncm_vector_get_variant (bao_a->x));
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
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass* diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_variant ("z",
                                                         NULL,
                                                         "Data redshift",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare   = &_nc_data_bao_a_prepare;
  diag_class->mean_func = &_nc_data_bao_a_mean_func;
  diag_class->set_size  = &_nc_data_bao_a_set_size;
}

static void
_nc_data_bao_a_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (bao_a->dist, cosmo);
}

static void 
_nc_data_bao_a_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (diag);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  guint i;
  
  for (i = 0; i < diag->np; i++)
  {
    const gdouble z = ncm_vector_get (bao_a->x, i);
    const gdouble A = nc_distance_bao_A_scale (bao_a->dist, cosmo, z);
    ncm_vector_set (vp, i, A);
  }
}

/**
 * nc_data_bao_a_new:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * FIXME
 *
 * Returns: a #NcmData
 */
NcmData *
nc_data_bao_a_new (NcDistance *dist, NcDataBaoId id)
{
  NcmData *data = g_object_new (NC_TYPE_DATA_BAO_A,
                                "dist", dist,
                                NULL);
  nc_data_bao_a_set_sample (NC_DATA_BAO_A (data), id);
  return data;
}

static void 
_nc_data_bao_a_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataBaoA *bao_a = NC_DATA_BAO_A (diag);

  if ((np == 0) || (np != diag->np))
    ncm_vector_clear (&bao_a->x);

  if ((np != 0) && (np != diag->np))
    bao_a->x = ncm_vector_new (np);
  
  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_bao_a_parent_class)->set_size (diag, np);
}

/**
 * nc_data_bao_a_set_sample:
 * @bao_a: a #NcDataBaoA
 * @id: a #NcDataBaoId
 *
 * FIXME
 *
 */
void 
nc_data_bao_a_set_sample (NcDataBaoA *bao_a, NcDataBaoId id)
{
  NcmData *data = NCM_DATA (bao_a);
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (bao_a);
  
  g_assert (id == NC_DATA_BAO_A_EISENSTEIN2005);

  ncm_data_set_desc (data, "Eisenstein 2005, BAO Sample A");

  ncm_data_gauss_diag_set_size (diag, 1);

  ncm_vector_set (bao_a->x,    0, ncm_c_bao_eisenstein_z ());
  ncm_vector_set (diag->y,     0, ncm_c_bao_eisenstein_A ());
  ncm_vector_set (diag->sigma, 0, ncm_c_bao_eisenstein_sigma_A ());

  ncm_data_set_init (data, TRUE);
}
