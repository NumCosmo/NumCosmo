/***************************************************************************
 *            nc_data_distance_mu.c
 *
 *  Thu Apr 22 10:37:22 2010
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
 * SECTION:nc_data_dist_mu
 * @title: NcDataDistMu
 * @short_description: Distance modulus data.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_dist_mu.h"
#include "math/ncm_data_gauss_diag.h"
#include "math/ncm_cfg.h"
#include "nc_enum_types.h"

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
      nc_data_dist_mu_set_dist (dist_mu, g_value_get_object (value));
      break;
    case PROP_Z:
      ncm_vector_substitute (&dist_mu->x, g_value_get_object (value), TRUE);
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
      g_value_set_object (value, dist_mu->x);
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
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_object ("z",
                                                        NULL,
                                                        "Data redshift",
                                                        NCM_TYPE_VECTOR,
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

static void 
_nc_data_dist_mu_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataDistMu *dist_mu = NC_DATA_DIST_MU (diag);
  NcHICosmo *cosmo      = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const gdouble f_lnRH  = 5.0 * log10 (nc_hicosmo_RH_Mpc (cosmo));
  const guint np        = ncm_data_gauss_diag_get_size (diag);
  guint i;

  for (i = 0; i < np; i++)
  {
    const gdouble z   = ncm_vector_get (dist_mu->x, i);
    const gdouble dmu = nc_distance_dmodulus (dist_mu->dist, cosmo, z);
    ncm_vector_set (vp, i, dmu + f_lnRH);
  }
}

static void 
_nc_data_dist_mu_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataDistMu *dist_mu = NC_DATA_DIST_MU (diag);
  const guint cnp       = ncm_data_gauss_diag_get_size (diag);

  if ((np == 0) || (np != cnp))
    ncm_vector_clear (&dist_mu->x);

  if ((np != 0) && (np != cnp))
    dist_mu->x = ncm_vector_new (np);
  
  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_dist_mu_parent_class)->set_size (diag, np);
}

/**
 * nc_data_dist_mu_new_empty:
 * @dist: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataDistMu *
nc_data_dist_mu_new_empty (NcDistance *dist)
{
  NcDataDistMu *dist_mu = g_object_new (NC_TYPE_DATA_DIST_MU,
                                        "dist", dist,
                                        "w-mean", TRUE,
                                        NULL);
  return dist_mu;
}

/**
 * nc_data_dist_mu_new_from_file:
 * @filename: file containing a serialized #NcDataDistMu
 * 
 * Creates a new #NcDataDistMu from @filename.
 * 
 * Returns: (transfer full): the newly created #NcDataDistMu.
 */
NcDataDistMu *
nc_data_dist_mu_new_from_file (const gchar *filename)
{
  NcDataDistMu *dist_mu = NC_DATA_DIST_MU (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_DIST_MU (dist_mu));

  return dist_mu;
}

/**
 * nc_data_dist_mu_new_from_id:
 * @dist: a #NcDistance
 * @id: a #NcDataSNIAId
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataDistMu *
nc_data_dist_mu_new_from_id (NcDistance *dist, NcDataSNIAId id)
{
  NcDataDistMu *dist_mu;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_SNIA_SIMPLE_GOLD_157:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_gold_157.obj", TRUE);
      break;
    case NC_DATA_SNIA_SIMPLE_GOLD_182:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_gold_182.obj", TRUE);
      break;
    case NC_DATA_SNIA_SIMPLE_GOLD_182_FULL:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_gold_206.obj", TRUE);
      break;
    case NC_DATA_SNIA_SIMPLE_ESSENCE:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_essence.obj", TRUE);
      break;
    case NC_DATA_SNIA_SIMPLE_LEGACY:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_legacy.obj", TRUE);
      break;
    case NC_DATA_SNIA_SIMPLE_UNION:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_union.obj", TRUE);
      break;
    case NC_DATA_SNIA_SIMPLE_CfA3:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_cfa3.obj", TRUE);
      break;
    case NC_DATA_SNIA_SIMPLE_UNION2:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_union2.obj", TRUE);
      break;
    case NC_DATA_SNIA_SIMPLE_UNION2_1:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_union2_1.obj", TRUE);
      break;
    case NC_DATA_SNIA_SIMPLE_SDSS_EMILLE:
      filename = ncm_cfg_get_data_filename ("nc_data_snia_diag_sdss_emille.obj", TRUE);
      break;
    default:
      g_error ("nc_data_dist_mu_new_from_id: id %d not recognized.", id);
      break;
  }

  dist_mu = nc_data_dist_mu_new_from_file (filename);
  nc_data_dist_mu_set_dist (dist_mu, dist);
  g_free (filename);

  return dist_mu;
}

/**
 * nc_data_dist_mu_set_dist:
 * @dist_mu: a #NcDataDistMu
 * @dist: a #NcDistance
 * 
 * Sets the distance object.
 * 
 */
void 
nc_data_dist_mu_set_dist (NcDataDistMu *dist_mu, NcDistance *dist)
{
  nc_distance_clear (&dist_mu->dist);
  dist_mu->dist = nc_distance_ref (dist);
}
