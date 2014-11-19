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
 * @title: Hubble Function Data
 * @short_description: Object representing Hubble Function data
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

#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

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
      ncm_vector_set_from_variant (hubble->x, g_value_get_variant (value));
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
      g_value_take_variant (value, ncm_vector_get_variant (hubble->x));
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
                                   g_param_spec_variant ("z",
                                                         NULL,
                                                         "Data redshifts",
                                                         G_VARIANT_TYPE ("ad"), NULL,
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

#ifdef NUMCOSMO_HAVE_SQLITE3
static gchar *_nc_data_hubble_function_query[] =
{
  "Simon 2005 H(z) sample", "SELECT z,p,s FROM kinematics WHERE param='Hz_Simon2005' ORDER BY z",
  "Cabre sample", "SELECT z,p,s FROM kinematics WHERE param='H_CABRE' ORDER BY z",
  "Stern 2009 H(z) sample", "SELECT z,p,s FROM kinematics WHERE param='Hz_Stern2009' ORDER BY z",
  "Moresco 2012 H(z) BC03  sample", "SELECT z,p,s FROM kinematics WHERE param='BC03_Moresco2012' ORDER BY z",
  "Moresco 2012 H(z) MaStro sample", "SELECT z,p,s FROM kinematics WHERE param='MaStro_Moresco2012' ORDER BY z",
  "Busca 2013 H(z) BAO+WMAP sample", "SELECT z,p,s FROM kinematics WHERE param='BAO+WMAP_Busca2013' ORDER BY z",
  "Riess 2008 HST Project sample", "SELECT z,p,s FROM kinematics WHERE param='HST_Riess2008' ORDER BY z",
  "Zhang 2012 Sloan sample", "SELECT z,p,s FROM kinematics WHERE param='SloanDR7_Zhang2012' ORDER BY z",
};
#endif

/**
 * nc_data_hubble_new:
 * @id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_hubble_new (NcDataHubbleId id)
{
  NcmData *data = g_object_new (NC_TYPE_DATA_HUBBLE,
                                NULL);
  nc_data_hubble_set_sample (NC_DATA_HUBBLE (data), id);
  return data;
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
 * nc_data_hubble_set_sample:
 * @hubble: a #NcDataHubble.
 * @id: FIXME
 *
 * FIXME
 *
 */
void
nc_data_hubble_set_sample (NcDataHubble *hubble, NcDataHubbleId id)
{
#ifdef NUMCOSMO_HAVE_SQLITE3
  NcmData *data = NCM_DATA (hubble);
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  
  g_assert (id < NC_DATA_HUBBLE_NSAMPLES);

  {
    const gchar *query = _nc_data_hubble_function_query[id * 2 + 1];
    gint i, nrow, qncol, ret;
    gchar **res;
    gchar *err_str;

    sqlite3 *db = ncm_cfg_get_default_sqlite3 ();

    ncm_data_set_desc (data, _nc_data_hubble_function_query[id * 2]);    

    g_assert (db != NULL);  

    ret = sqlite3_get_table (db, query, &res, &nrow, &qncol, &err_str);
    if (ret != SQLITE_OK)
    {
      sqlite3_free_table (res);
      g_error ("nc_data_hubble_set_sample: Query error: %s", err_str);
    }

    g_assert_cmpint (nrow, >, 0);

    ncm_data_gauss_diag_set_size (diag, nrow);

    for (i = 0; i < nrow; i++)
    {
      gint j = 0;
      ncm_vector_set (hubble->x,   i, atof (res[(i + 1) * qncol + j++]));
      ncm_vector_set (diag->y,     i, atof (res[(i + 1) * qncol + j++]));  
      ncm_vector_set (diag->sigma, i, atof (res[(i + 1) * qncol + j++]));
    }

    sqlite3_free_table (res);

    ncm_data_set_init (data, TRUE);
  }
#else
  g_error (PACKAGE_NAME" compiled without support for sqlite3, Hubble data not avaliable.");
#endif
}
