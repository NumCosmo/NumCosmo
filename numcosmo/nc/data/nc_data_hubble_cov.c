/***************************************************************************
 *            nc_data_hubble_cov.c
 *
 *  Wed Oct 21 13:26:04 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_hubble_cov.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcDataHubbleCov:
 *
 * Hubble function data with a full covariance.
 *
 * Gaussian likelihood for measurements of $H(z)$ at the redshifts in
 * #NcDataHubbleCov:z, with the covariance held by the #NcmDataGaussCov parent.
 * The theory vector is nc_hicosmo_H() at each redshift. Use it instead of
 * #NcDataHubble when the measurements are correlated, for example cosmic
 * chronometers sharing a systematic error.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_hubble_cov.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/background/nc_hicosmo.h"
#include "nc_enum_types.h"

struct _NcDataHubbleCov
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  NcmVector *x;
};

enum
{
  PROP_0,
  PROP_Z,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataHubbleCov, nc_data_hubble_cov, NCM_TYPE_DATA_GAUSS_COV)

static void
nc_data_hubble_cov_init (NcDataHubbleCov *hubble_cov)
{
  hubble_cov->x = NULL;
}

static void
_nc_data_hubble_cov_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataHubbleCov *hubble_cov = NC_DATA_HUBBLE_COV (object);

  g_return_if_fail (NC_IS_DATA_HUBBLE_COV (object));

  switch (prop_id)
  {
    case PROP_Z:
      ncm_vector_substitute (&hubble_cov->x, g_value_get_object (value), TRUE);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_data_hubble_cov_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataHubbleCov *hubble_cov = NC_DATA_HUBBLE_COV (object);

  g_return_if_fail (NC_IS_DATA_HUBBLE_COV (object));

  switch (prop_id)
  {
    case PROP_Z:
      g_value_set_object (value, hubble_cov->x);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_data_hubble_cov_dispose (GObject *object)
{
  NcDataHubbleCov *hubble_cov = NC_DATA_HUBBLE_COV (object);

  ncm_vector_clear (&hubble_cov->x);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_hubble_cov_parent_class)->dispose (object);
}

static void
_nc_data_hubble_cov_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_hubble_cov_parent_class)->finalize (object);
}

static void _nc_data_hubble_cov_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);
static void _nc_data_hubble_cov_set_size (NcmDataGaussCov *gauss, guint np);

static void
nc_data_hubble_cov_class_init (NcDataHubbleCovClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &_nc_data_hubble_cov_set_property;
  object_class->get_property = &_nc_data_hubble_cov_get_property;
  object_class->dispose      = &_nc_data_hubble_cov_dispose;
  object_class->finalize     = &_nc_data_hubble_cov_finalize;

  /**
   * NcDataHubbleCov:z:
   *
   * Redshifts of the $H(z)$ measurements.
   */
  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_object ("z",
                                                        NULL,
                                                        "Data redshifts",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gauss_class->mean_func = &_nc_data_hubble_cov_mean_func;
  gauss_class->set_size  = &_nc_data_hubble_cov_set_size;
}

static void
_nc_data_hubble_cov_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataHubbleCov *hubble_cov = NC_DATA_HUBBLE_COV (gauss);
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const guint np              = ncm_data_gauss_cov_get_size (gauss);
  guint i;

  for (i = 0; i < np; i++)
  {
    const gdouble z = ncm_vector_get (hubble_cov->x, i);

    ncm_vector_set (vp, i, nc_hicosmo_H (cosmo, z));
  }
}

static void
_nc_data_hubble_cov_set_size (NcmDataGaussCov *gauss, guint np)
{
  NcDataHubbleCov *hubble_cov = NC_DATA_HUBBLE_COV (gauss);
  const guint cnp             = ncm_data_gauss_cov_get_size (gauss);

  if ((np == 0) || (np != cnp))
    ncm_vector_clear (&hubble_cov->x);

  if ((np != 0) && (np != cnp))
    hubble_cov->x = ncm_vector_new (np);

  /* Chain up : end */
  NCM_DATA_GAUSS_COV_CLASS (nc_data_hubble_cov_parent_class)->set_size (gauss, np);
}

/**
 * nc_data_hubble_cov_new_empty:
 *
 * Creates a new empty #NcDataHubbleCov.
 *
 * Returns: (transfer full): a new #NcDataHubbleCov.
 */
NcDataHubbleCov *
nc_data_hubble_cov_new_empty (void)
{
  NcDataHubbleCov *hubble_cov = g_object_new (NC_TYPE_DATA_HUBBLE_COV,
                                              NULL);

  return hubble_cov;
}

/**
 * nc_data_hubble_cov_new_from_file:
 * @filename: file containing a serialized #NcDataHubbleCov
 *
 * Creates a new #NcDataHubbleCov from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataHubbleCov.
 */
NcDataHubbleCov *
nc_data_hubble_cov_new_from_file (const gchar *filename)
{
  NcDataHubbleCov *hubble_cov = NC_DATA_HUBBLE_COV (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_HUBBLE_COV (hubble_cov));

  return hubble_cov;
}

/**
 * nc_data_hubble_cov_new_from_id:
 * @id: a #NcDataHubbleCovId
 *
 * Creates a new #NcDataHubbleCov from the data file of @id.
 *
 * Returns: (transfer full): a new #NcDataHubbleCov.
 */
NcDataHubbleCov *
nc_data_hubble_cov_new_from_id (NcDataHubbleCovId id)
{
  NcDataHubbleCov *hubble_cov;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_HUBBLE_COV_LOUBSER2025:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_cov_loubser2025.obj", TRUE);
      break;
    case NC_DATA_HUBBLE_COV_MORESCO2020:
      filename = ncm_cfg_get_data_filename ("nc_data_hubble_cov_moresco2020.obj", TRUE);
      break;
    default:                                                                 /* LCOV_EXCL_LINE */
      g_error ("nc_data_hubble_cov_new_from_id: id %d not recognized.", id); /* LCOV_EXCL_LINE */
      break;                                                                 /* LCOV_EXCL_LINE */
  }

  hubble_cov = nc_data_hubble_cov_new_from_file (filename);
  g_free (filename);

  return hubble_cov;
}

