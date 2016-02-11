/***************************************************************************
 *            nc_data_hubble.c
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
 * SECTION:nc_data_hubble
 * @title: NcDataHubble
 * @short_description: Hubble function data.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_hubble.h"

#include "math/ncm_data_gauss_diag.h"
#include "math/ncm_cfg.h"
#include "nc_hicosmo.h"
#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_Z,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataHubble, nc_data_hubble, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_hubble_init (NcDataHubble *hubble)
{
  hubble->x = NULL;
}

static void
_nc_data_hubble_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_hubble_parent_class)->constructed (object);
}

static void
nc_data_hubble_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataHubble *hubble = NC_DATA_HUBBLE (object);
  g_return_if_fail (NC_IS_DATA_HUBBLE (object));

  switch (prop_id)
  {
    case PROP_Z:
      ncm_vector_clear (&hubble->x);
      hubble->x = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_hubble_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataHubble *hubble = NC_DATA_HUBBLE (object);
  g_return_if_fail (NC_IS_DATA_HUBBLE (object));

  switch (prop_id)
  {
    case PROP_Z:
      g_value_set_object (value, hubble->x);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_hubble_dispose (GObject *object)
{
  NcDataHubble *hubble = NC_DATA_HUBBLE (object);

  ncm_vector_clear (&hubble->x);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_hubble_parent_class)->dispose (object);
}

static void
nc_data_hubble_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_hubble_parent_class)->finalize (object);
}

static void _nc_data_hubble_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_hubble_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);
static void _nc_data_hubble_set_size (NcmDataGaussDiag *diag, guint np);

static void
nc_data_hubble_class_init (NcDataHubbleClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass* diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->constructed  = &_nc_data_hubble_constructed;
  object_class->set_property = &nc_data_hubble_set_property;
  object_class->get_property = &nc_data_hubble_get_property;
  object_class->dispose      = &nc_data_hubble_dispose;
  object_class->finalize     = &nc_data_hubble_finalize;

  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_object ("z",
                                                        NULL,
                                                        "Data redshifts",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->prepare   = &_nc_data_hubble_prepare;
  diag_class->mean_func = &_nc_data_hubble_mean_func;
  diag_class->set_size  = &_nc_data_hubble_set_size;
}

static void
_nc_data_hubble_prepare (NcmData *data, NcmMSet *mset)
{
  NCM_UNUSED (data);
  NCM_UNUSED (mset);
  /* Nothing to do... */
}

static void 
_nc_data_hubble_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataHubble *hubble = NC_DATA_HUBBLE (diag);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  guint i;

  for (i = 0; i < diag->np; i++)
  {
    const gdouble z = ncm_vector_get (hubble->x, i);
    const gdouble H = nc_hicosmo_H (cosmo, z); 
    ncm_vector_set (vp, i, H);
  }
}

void 
_nc_data_hubble_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataHubble *hubble = NC_DATA_HUBBLE (diag);
  
  if ((np == 0) || (np != diag->np))
    ncm_vector_clear (&hubble->x);

  if ((np != 0) && (np != diag->np))
    hubble->x = ncm_vector_new (np);
  
  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_hubble_parent_class)->set_size (diag, np);  
}

/**
 * nc_data_hubble_new_empty:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataHubble *
nc_data_hubble_new_empty (void)
{
  NcDataHubble *hubble = g_object_new (NC_TYPE_DATA_HUBBLE,
                                       NULL);
  return hubble;
}

/**
 * nc_data_hubble_new_from_file:
 * @filename: file containing a serialized #NcDataHubble
 * 
 * Creates a new #NcDataHubble from @filename.
 * 
 * Returns: (transfer full): the newly created #NcDataHubble.
 */
NcDataHubble *
nc_data_hubble_new_from_file (const gchar *filename)
{
  NcDataHubble *hubble = NC_DATA_HUBBLE (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_HUBBLE (hubble));

  return hubble;
}

/**
 * nc_data_hubble_new_from_id:
 * @id: a #NcDataHubbleId
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataHubble *
nc_data_hubble_new_from_id (NcDataHubbleId id)
{
  NcDataHubble *hubble;
  gchar *filename;
  switch (id)
  {
    case NC_DATA_HUBBLE_SIMON2005:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_simon2005.obj", TRUE);
      break;
    case NC_DATA_HUBBLE_CABRE:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_cabre.obj", TRUE);
      break;
    case NC_DATA_HUBBLE_STERN2009:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_stern2009.obj", TRUE);
      break;
    case NC_DATA_HUBBLE_MORESCO2012_BC03:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_moresco2012_bc03.obj", TRUE);
      break;
    case NC_DATA_HUBBLE_MORESCO2012_MASTRO:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_moresco2012_mastro.obj", TRUE);
      break;
    case NC_DATA_HUBBLE_BUSCA2013_BAO_WMAP:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_busca2013_bao_wmap.obj", TRUE);
      break;
    case NC_DATA_HUBBLE_RIESS2008_HST:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_riess2008_hst.obj", TRUE);
      break;
    case NC_DATA_HUBBLE_ZHANG2012:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_zhang2012.obj", TRUE);
      break;
    default:
      g_error ("nc_data_hubble_new_from_id: id %d not recognized.", id);
      break;
  }

  hubble = nc_data_hubble_new_from_file (filename);
  g_free (filename);

  return hubble;
}
