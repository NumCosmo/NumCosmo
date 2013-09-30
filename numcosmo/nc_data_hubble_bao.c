/***************************************************************************
 *            nc_data_hubble_bao.c
 *
 *  Thu Apr 22 14:34:54 2010
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
 * SECTION:nc_data_hubble_bao
 * @title: Hubble Function Data from BAO
 * @short_description: Object representing Hubble Function BAO data
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_hubble_bao.h"

#include "math/ncm_data_gauss_diag.h"
#include "math/ncm_cfg.h"
#include "nc_hicosmo.h"
#include "nc_enum_types.h"

#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataHubbleBao, nc_data_hubble_bao, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_hubble_bao_init (NcDataHubbleBao *hubble_bao)
{
  hubble_bao->dist = NULL;
  hubble_bao->x    = NULL;
}

static void
_nc_data_hubble_bao_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_hubble_bao_parent_class)->constructed (object);
}

static void
nc_data_hubble_bao_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataHubbleBao *hubble_bao = NC_DATA_HUBBLE_BAO (object);
  g_return_if_fail (NC_IS_DATA_HUBBLE_BAO (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_distance_clear (&hubble_bao->dist);
      hubble_bao->dist = g_value_dup_object (value);
      break;
    case PROP_Z:
      ncm_vector_set_from_variant (hubble_bao->x, g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_hubble_bao_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataHubbleBao *hubble_bao = NC_DATA_HUBBLE_BAO (object);
  g_return_if_fail (NC_IS_DATA_HUBBLE_BAO (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, hubble_bao->dist);
      break;
    case PROP_Z:
      g_value_take_variant (value, ncm_vector_get_variant (hubble_bao->x));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_hubble_bao_dispose (GObject *object)
{
  NcDataHubbleBao *hubble_bao = NC_DATA_HUBBLE_BAO (object);

  nc_distance_clear (&hubble_bao->dist);
  ncm_vector_clear (&hubble_bao->x);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_hubble_bao_parent_class)->dispose (object);
}

static void
nc_data_hubble_bao_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_hubble_bao_parent_class)->finalize (object);
}

static void _nc_data_hubble_bao_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_hubble_bao_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);
static void _nc_data_hubble_bao_set_size (NcmDataGaussDiag *diag, guint np);

static void
nc_data_hubble_bao_class_init (NcDataHubbleBaoClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass* diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->constructed  = &_nc_data_hubble_bao_constructed;
  object_class->set_property = &nc_data_hubble_bao_set_property;
  object_class->get_property = &nc_data_hubble_bao_get_property;
  object_class->dispose      = &nc_data_hubble_bao_dispose;
  object_class->finalize     = &nc_data_hubble_bao_finalize;

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
                                                         "Data redshifts",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare   = &_nc_data_hubble_bao_prepare;
  diag_class->mean_func = &_nc_data_hubble_bao_mean_func;
  diag_class->set_size  = &_nc_data_hubble_bao_set_size;
}

static void
_nc_data_hubble_bao_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataHubbleBao *hubble_bao = NC_DATA_HUBBLE_BAO (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (hubble_bao->dist, cosmo);
}

static void 
_nc_data_hubble_bao_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataHubbleBao *hubble_bao = NC_DATA_HUBBLE_BAO (diag);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  gdouble r_zd;
  guint i;

  if (ncm_model_impl (NCM_MODEL (cosmo)) & NC_HICOSMO_IMPL_as_drag)
    r_zd = nc_hicosmo_as_drag (cosmo);
  else
  {   
    gdouble zd = nc_distance_drag_redshift (hubble_bao->dist, cosmo);
    r_zd = nc_distance_sound_horizon (hubble_bao->dist, cosmo, zd);
  }

  for (i = 0; i < diag->np; i++)
  {
    const gdouble z = ncm_vector_get (hubble_bao->x, i);
    const gdouble E = nc_hicosmo_E (cosmo, z);
    const gdouble HBAO = E * r_zd / (1.0 + z);
    ncm_vector_set (vp, i, HBAO);
  }
}

/**
 * nc_data_hubble_bao_new:
 * @dist: FIXME
 * @id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_hubble_bao_new (NcDistance *dist, NcDataHubbleBaoId id)
{
  NcmData *data = g_object_new (NC_TYPE_DATA_HUBBLE_BAO,
                                "dist", dist,
                                NULL);
  nc_data_hubble_bao_set_sample (NC_DATA_HUBBLE_BAO (data), id);
  return data;
}

static void 
_nc_data_hubble_bao_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataHubbleBao *hubble_bao = NC_DATA_HUBBLE_BAO (diag);

  if (diag->np != 0)
    g_assert (hubble_bao->x != NULL && ncm_vector_len (hubble_bao->x) == diag->np);
  
  if ((np == 0) || (np != diag->np))
    ncm_vector_clear (&hubble_bao->x);

  if ((np != 0) && (np != diag->np))
    hubble_bao->x = ncm_vector_new (np);

  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_hubble_bao_parent_class)->set_size (diag, np);  
}

/**
 * nc_data_hubble_bao_set_sample:
 * @hubble_bao: a #NcDataHubbleBao.
 * @id: FIXME
 *
 * FIXME
 *
 */
void
nc_data_hubble_bao_set_sample (NcDataHubbleBao *hubble_bao, NcDataHubbleBaoId id)
{
  NcmData *data = NCM_DATA (hubble_bao);
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble_bao);
  
  g_assert (id < NC_DATA_HUBBLE_BAO_NSAMPLES);

  ncm_data_set_desc (data, "Busca 2013 H(z)r_s(z_d)/(1+z) from BAO"); /* Busca 2013 Eq. (27) */
  ncm_data_gauss_diag_set_size (diag, 1);
  
  ncm_vector_set (hubble_bao->x,  0, 2.3);
  ncm_vector_set (diag->y,        0, 3.455724026e-2); /* (1.036 * 10^4 km / s) / c */
  ncm_vector_set (diag->sigma,    0, 0.120083074e-2); /* (0.036 * 10^4 km / s) / c */

  ncm_data_set_init (data, TRUE);
}
