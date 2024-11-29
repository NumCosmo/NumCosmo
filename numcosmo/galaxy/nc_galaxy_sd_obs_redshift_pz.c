/***************************************************************************
 *            nc_galaxy_sd_obs_redshift_pz.c
 *
 *  Mon Nov 25 20:28:56 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift_pz.c
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:nc_galaxy_sd_obs_redshift_pz
 * @title: NcGalaxySDObsRedshiftPz
 * @short_description: Class describing photometric redshift observations with a spline.
 * @stability: Unstable
 *
 *
 * Class describing photometric redshift observations with a spline.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "galaxy/nc_galaxy_sd_obs_redshift_pz.h"
#include "math/ncm_spline.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_stats_dist1d_spline.h"

typedef struct _NcGalaxySDObsRedshiftPzPrivate
{
} NcGalaxySDObsRedshiftPzPrivate;

struct _NcGalaxySDObsRedshiftPz
{
  NcGalaxySDObsRedshift parent_instance;
};

typedef struct _NcGalaxySDObsRedshiftPzData
{
  NcmSpline *pz;
} NcGalaxySDObsRedshiftPzData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDObsRedshiftPz, nc_galaxy_sd_obs_redshift_pz, NC_TYPE_GALAXY_SD_OBS_REDSHIFT);

static void
nc_galaxy_sd_obs_redshift_pz_init (NcGalaxySDObsRedshiftPz *gsdorpz)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_obs_redshift_pz_dispose (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_pz_parent_class)->dispose (object);
}

static void
nc_galaxy_sd_obs_redshift_pz_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_pz_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_obs_redshift_pz_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng);
static NcGalaxySDObsRedshiftIntegrand *_nc_galaxy_sd_obs_redshift_pz_integ (NcGalaxySDObsRedshift *gsdor);
static void _nc_galaxy_sd_obs_redshift_pz_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data);

static void
nc_galaxy_sd_obs_redshift_pz_class_init (NcGalaxySDObsRedshiftPzClass *klass)
{
  NcGalaxySDObsRedshiftClass *gsdor_class = NC_GALAXY_SD_OBS_REDSHIFT_CLASS (klass);
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_galaxy_sd_obs_redshift_pz_dispose;
  object_class->finalize = &nc_galaxy_sd_obs_redshift_pz_finalize;

  ncm_model_class_set_name_nick (model_class, "P(z) Observed Redshift", "GalaxySDObsRedshiftPz");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_model_class_check_params_info (model_class);

  gsdor_class->gen       = &_nc_galaxy_sd_obs_redshift_pz_gen;
  gsdor_class->integ     = &_nc_galaxy_sd_obs_redshift_pz_integ;
  gsdor_class->data_init = &_nc_galaxy_sd_obs_redshift_pz_data_init;
}

static void
_nc_galaxy_sd_obs_redshift_pz_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftPzData * const ldata = (NcGalaxySDObsRedshiftPzData *) data->ldata;
  NcmVector *xv                             = ncm_spline_get_xv (ldata->pz);
  NcmVector *yv                             = ncm_spline_get_yv (ldata->pz);
  NcmVector *m2lnyv                         = ncm_vector_new (ncm_vector_len (yv));
  NcmSpline *m2lnp;
  NcmStatsDist1d *dist;
  guint i;

  for (i = 0; i < ncm_vector_len (yv); i++)
    ncm_vector_set (m2lnyv, i, -2.0 * log (ncm_vector_get (yv, i)));

  m2lnp   = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, m2lnyv, TRUE));
  dist    = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new (m2lnp));
  data->z = ncm_stats_dist1d_gen (dist, rng);
}

struct _IntegData
{
  NcGalaxySDObsRedshiftPz *gsdorpz;
};

static gpointer
_integ_data_copy (gpointer idata)
{
  struct _IntegData *new_idata = g_new0 (struct _IntegData, 1);

  *new_idata = *(struct _IntegData *) idata;

  return new_idata;
}

static void
_integ_data_free (gpointer idata)
{
  g_free (idata);
}

static gdouble
_nc_galaxy_sd_obs_redshift_pz_integ_f (gpointer callback_data, const gdouble z, NcGalaxySDObsRedshiftData *data)
{
  NcGalaxySDObsRedshiftPzData * const ldata = (NcGalaxySDObsRedshiftPzData *) data->ldata;

  return ncm_spline_eval (ldata->pz, z);
}

static NcGalaxySDObsRedshiftIntegrand *
_nc_galaxy_sd_obs_redshift_pz_integ (NcGalaxySDObsRedshift *gsdor)
{
  NcGalaxySDObsRedshiftPz *gsdorpz      = NC_GALAXY_SD_OBS_REDSHIFT_PZ (gsdor);
  struct _IntegData *int_data           = g_new0 (struct _IntegData, 1);
  NcGalaxySDObsRedshiftIntegrand *integ = nc_galaxy_sd_obs_redshift_integrand_new (&_nc_galaxy_sd_obs_redshift_pz_integ_f,
                                                                                   &_integ_data_free,
                                                                                   &_integ_data_copy,
                                                                                   NULL,
                                                                                   int_data);

  int_data->gsdorpz = gsdorpz;

  return integ;
}

static void
_nc_galaxy_sd_obs_redshift_pz_ldata_free (gpointer ldata)
{
  NcGalaxySDObsRedshiftPzData *data = (NcGalaxySDObsRedshiftPzData *) ldata;

  ncm_spline_free (data->pz);
  g_free (ldata);
}

static void
_nc_galaxy_sd_obs_redshift_pz_ldata_read_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDObsRedshiftPzData *ldata = (NcGalaxySDObsRedshiftPzData *) data->ldata;

  ldata->pz = ncm_spline_ref (nc_galaxy_wl_obs_peek_pz (obs, i));
}

static void
_nc_galaxy_sd_obs_redshift_pz_ldata_write_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDObsRedshiftPzData *ldata = (NcGalaxySDObsRedshiftPzData *) data->ldata;

  nc_galaxy_wl_obs_set_pz (obs, i, ldata->pz);
}

static void
_nc_galaxy_sd_obs_redshift_pz_ldata_required_columns (NcGalaxySDObsRedshiftData *data, GList *columns)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_obs_redshift_pz_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  NcGalaxySDObsRedshiftPzData *ldata = g_new0 (NcGalaxySDObsRedshiftPzData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_sd_obs_redshift_pz_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_sd_obs_redshift_pz_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_sd_obs_redshift_pz_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_sd_obs_redshift_pz_ldata_required_columns;
}

/**
 * nc_galaxy_sd_obs_redshift_pz_new:
 *
 * Creates a new #NcGalaxySDObsRedshiftPz object.
 *
 * Returns: (transfer full): a new #NcGalaxySDObsRedshiftPz object.
 */
NcGalaxySDObsRedshiftPz *
nc_galaxy_sd_obs_redshift_pz_new ()
{
  NcGalaxySDObsRedshiftPz *gsdorpz = g_object_new (NC_TYPE_GALAXY_SD_OBS_REDSHIFT_PZ, NULL);

  return gsdorpz;
}

/**
 * nc_galaxy_sd_obs_redshift_pz_ref:
 * @gsdorpz: a #NcGalaxySDObsRedshiftPz
 *
 * Increases the reference count of @gsdorpz by one.
 *
 * Returns: (transfer full): @gsdorpz.
 */
NcGalaxySDObsRedshiftPz *
nc_galaxy_sd_obs_redshift_pz_ref (NcGalaxySDObsRedshiftPz *gsdorpz)
{
  return g_object_ref (gsdorpz);
}

/**
 * nc_galaxy_sd_obs_redshift_pz_free:
 * @gsdorpz: a #NcGalaxySDObsRedshiftPz
 *
 * Decreases the reference count of @gsdorpz by one.
 *
 */
void
nc_galaxy_sd_obs_redshift_pz_free (NcGalaxySDObsRedshiftPz *gsdorpz)
{
  g_object_unref (gsdorpz);
}

/**
 * nc_galaxy_sd_obs_redshift_pz_clear:
 * @gsdorpz: a #NcGalaxySDObsRedshiftPz
 *
 * Decreases the reference count of @gsdorpz by one, and sets the pointer *@gsdorpz to
 * NULL.
 *
 */
void
nc_galaxy_sd_obs_redshift_pz_clear (NcGalaxySDObsRedshiftPz **gsdorpz)
{
  g_clear_object (gsdorpz);
}

/**
 * nc_galaxy_sd_obs_redshift_pz_gen:
 * @gsdorpz: a #NcGalaxySDObsRedshiftPz
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDObsRedshiftData
 * @rng: a #NcmRNG
 *
 * Sets the required columns for the data and generates a redshift observation.
 */
void
nc_galaxy_sd_obs_redshift_pz_gen (NcGalaxySDObsRedshiftPz *gsdorpz, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftClass *klass = NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdorpz);

  klass->gen (NC_GALAXY_SD_OBS_REDSHIFT (gsdorpz), data, rng);
}

/**
 * nc_galaxy_sd_obs_redshift_pz_data_set:
 * @gsdorpz: a #NcGalaxySDObsRedshiftPz
 * @data: a #NcGalaxySDObsRedshiftPzData
 * @spline: a #NcmSpline
 *
 * Sets the data of the photometric redshift observations.
 *
 */
void
nc_galaxy_sd_obs_redshift_pz_data_set (NcGalaxySDObsRedshiftPz *gsdorpz, NcGalaxySDObsRedshiftData *data, NcmSpline *spline)
{
  NcGalaxySDObsRedshiftPzData * const ldata = (NcGalaxySDObsRedshiftPzData *) data->ldata;

  ldata->pz = spline;
}

/**
 * nc_galaxy_sd_obs_redshift_pz_data_get:
 * @gsdorpz: a #NcGalaxySDObsRedshiftPz
 * @data: a #NcGalaxySDObsRedshiftPzData
 * @spline: (out): a #NcmSpline
 *
 * Gets the data of the photometric redshift observations.
 *
 */
void
nc_galaxy_sd_obs_redshift_pz_data_get (NcGalaxySDObsRedshiftPz *gsdorpz, NcGalaxySDObsRedshiftData *data, NcmSpline **spline)
{
  NcGalaxySDObsRedshiftPzData * const ldata = (NcGalaxySDObsRedshiftPzData *) data->ldata;

  *spline = ldata->pz;
}

