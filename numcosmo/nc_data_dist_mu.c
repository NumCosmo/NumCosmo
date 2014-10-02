/***************************************************************************
 *            nc_data_distance_mu.c
 *
 *  Thu Apr 22 10:37:22 2010
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
 * SECTION:nc_data_dist_mu
 * @title: Distance Modulus Data
 * @short_description: Data samples of distance modulus
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_dist_mu.h"
#include "math/ncm_data_gauss_diag.h"
#include "math/ncm_cfg.h"
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

G_DEFINE_TYPE (NcDataDistMu, nc_data_dist_mu, NCM_TYPE_DATA_GAUSS_DIAG);

static void
nc_data_dist_mu_init (NcDataDistMu *dist_mu)
{
  dist_mu->x    = NULL;
  dist_mu->dist = NULL;
}

static void
_nc_data_dist_mu_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_dist_mu_parent_class)->constructed (object);
}

static void
nc_data_dist_mu_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataDistMu *dist_mu = NC_DATA_DIST_MU (object);
  g_return_if_fail (NC_IS_DATA_DIST_MU (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_distance_clear (&dist_mu->dist);
      dist_mu->dist = g_value_dup_object (value);
      break;
    case PROP_Z:
      ncm_vector_set_from_variant (dist_mu->x, g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_dist_mu_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataDistMu *dist_mu = NC_DATA_DIST_MU (object);
  g_return_if_fail (NC_IS_DATA_DIST_MU (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, dist_mu->dist);
      break;
    case PROP_Z:
      g_value_take_variant (value, ncm_vector_get_variant (dist_mu->x));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_dist_mu_dispose (GObject *object)
{
  NcDataDistMu *dist_mu = NC_DATA_DIST_MU (object);

  ncm_vector_clear (&dist_mu->x);
  nc_distance_clear (&dist_mu->dist);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_dist_mu_parent_class)->dispose (object);
}

static void
nc_data_dist_mu_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_dist_mu_parent_class)->finalize (object);
}

static void _nc_data_dist_mu_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_dist_mu_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);
static void _nc_data_dist_mu_set_size (NcmDataGaussDiag *diag, guint np);

static void
nc_data_dist_mu_class_init (NcDataDistMuClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass* diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->constructed  = &_nc_data_dist_mu_constructed;
  object_class->set_property = &nc_data_dist_mu_set_property;
  object_class->get_property = &nc_data_dist_mu_get_property;
  object_class->dispose      = &nc_data_dist_mu_dispose;
  object_class->finalize     = &nc_data_dist_mu_finalize;

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


  data_class->prepare   = &_nc_data_dist_mu_prepare;
  diag_class->mean_func = &_nc_data_dist_mu_mean_func;
  diag_class->set_size  = &_nc_data_dist_mu_set_size;
}

static void
_nc_data_dist_mu_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataDistMu *dist_mu = NC_DATA_DIST_MU (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (dist_mu->dist, cosmo);
}

#ifdef NUMCOSMO_HAVE_SQLITE3
static gchar *_nc_data_snia_query[] =
{
  "Gold sample 157", "SELECT z,mu-0.32 AS muc,s FROM supernovae WHERE quality='Gold-2004' ORDER BY z",
  "Gold sample 182 - removed low redshift", "SELECT z,mu-0.32 AS muc,s FROM supernovae WHERE quality='Gold' AND z >= 0.0233 ORDER BY z",
  "Gold sample 182", "SELECT z,mu-0.32 AS muc,s FROM supernovae WHERE quality='Gold' ORDER BY z",
  "ESSENCE sample",  "SELECT z,mu,s FROM supernovae WHERE quality='ESSENCE' ORDER BY z",
  "Legacy sample",   "SELECT z,mu+19.308,s FROM supernovae WHERE quality='LEGACY' ORDER BY z",
  "Union sample",    "SELECT z,mu,s FROM supernovae WHERE quality='UNION' ORDER BY z",
  "CfA3 sample",     "SELECT z,mu,s FROM supernovae WHERE quality='CfA3' ORDER BY z",
  "Union2 sample",   "SELECT z,mu,s FROM supernovae WHERE quality='Union2' ORDER BY z",
  "Union2.1 sample", "SELECT z,mu,s FROM supernovae WHERE quality='Union2.1' ORDER BY z",
  "SDSS sample by Emille", "SELECT z,mu,s FROM supernovae WHERE quality='SDSS-Emille' ORDER BY z",
};
#endif

static void 
_nc_data_dist_mu_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataDistMu *dist_mu = NC_DATA_DIST_MU (diag);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  guint i;

  for (i = 0; i < diag->np; i++)
  {
    const gdouble z  = ncm_vector_get (dist_mu->x, i);
    const gdouble mu = nc_distance_modulus (dist_mu->dist, cosmo, z);
    ncm_vector_set (vp, i, mu);
  }
}

/**
 * nc_data_dist_mu_new:
 * @dist: FIXME
 * @id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_dist_mu_new (NcDistance *dist, NcDataSNIAId id)
{
  NcmData *data = g_object_new (NC_TYPE_DATA_DIST_MU,
                                "dist", dist,
                                "w-mean", TRUE,
                                NULL);
  nc_data_dist_mu_set_sample (NC_DATA_DIST_MU (data), id);
  return data;
}

static void 
_nc_data_dist_mu_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataDistMu *dist_mu = NC_DATA_DIST_MU (diag);

  if ((np == 0) || (np != diag->np))
    ncm_vector_clear (&dist_mu->x);

  if ((np != 0) && (np != diag->np))
    dist_mu->x = ncm_vector_new (np);
  
  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_dist_mu_parent_class)->set_size (diag, np);
}

/**
 * nc_data_dist_mu_set_sample:
 * @dist_mu: a #NcDataDistMu.
 * @id: FIXME
 *
 * FIXME
 *
 */
void
nc_data_dist_mu_set_sample (NcDataDistMu *dist_mu, NcDataSNIAId id)
{
#ifdef NUMCOSMO_HAVE_SQLITE3
  NcmData *data = NCM_DATA (dist_mu);
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (dist_mu);

  g_assert (id <= NC_DATA_SNIA_SIMPLE_END && id >= NC_DATA_SNIA_SIMPLE_START);
  id -= NC_DATA_SNIA_SIMPLE_START;

  {
    const gchar *query = _nc_data_snia_query[id * 2 + 1];
    gint i, nrow, qncol, ret;
    gchar **res;
    gchar *err_str;

    sqlite3 *db = ncm_cfg_get_default_sqlite3 ();

    ncm_data_set_desc (data, _nc_data_snia_query[id * 2]);

    g_assert (db != NULL);  

    ret = sqlite3_get_table (db, query, &res, &nrow, &qncol, &err_str);
    if (ret != SQLITE_OK)
    {
      sqlite3_free_table (res);
      g_error ("nc_data_dist_mu_set_sample: Query error: %s", err_str);
    }

    ncm_data_gauss_diag_set_size (diag, nrow);

    for (i = 0; i < nrow; i++)
    {
      gint j = 0;
      ncm_vector_set (dist_mu->x,  i, atof (res[(i + 1) * qncol + j++]));
      ncm_vector_set (diag->y,     i, atof (res[(i + 1) * qncol + j++]));  
      ncm_vector_set (diag->sigma, i, atof (res[(i + 1) * qncol + j++]));
    }

    sqlite3_free_table (res);

    ncm_data_set_init (data, TRUE);
  }
#else
  g_error (PACKAGE_NAME" compiled without support for sqlite3, SN Ia data not avaliable.");
#endif

}
