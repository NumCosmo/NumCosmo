/***************************************************************************
 *            nc_galaxy_redshift_obs_gauss.c
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_obs_gauss.c
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyRedshiftObsGauss:
 *
 * Gaussian photometric-redshift observable model.
 *
 * The photometric redshift is Gaussian about the true redshift with a
 * redshift-dependent scatter $\sigma_z = \sigma_0 (1 + z)$:
 * $P(z_\mathrm{phot}|z) = \mathcal{N}(z_\mathrm{phot}; z, \sigma_z)$. The
 * per-galaxy observation $(z_\mathrm{phot}, \sigma_0)$ is carried together in
 * the #NcGalaxyRedshiftObsData, read from the "zp" and "sigma0" columns.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_obs_gauss.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyRedshiftObsGauss
{
  NcGalaxyRedshiftObs parent_instance;
};

typedef struct _NcGalaxyRedshiftObsGaussLData
{
  gdouble zp;     /* photometric-redshift point estimate */
  gdouble sigma0; /* scatter parameter: sigma_z = sigma0 (1 + z) */
} NcGalaxyRedshiftObsGaussLData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyRedshiftObsGauss, nc_galaxy_redshift_obs_gauss, NC_TYPE_GALAXY_REDSHIFT_OBS);

static void
nc_galaxy_redshift_obs_gauss_init (NcGalaxyRedshiftObsGauss *gsdreg)
{
}

static void
_nc_galaxy_redshift_obs_gauss_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_obs_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_redshift_obs_gauss_data_init (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data);
static gdouble _nc_galaxy_redshift_obs_gauss_eval (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z);
static gdouble _nc_galaxy_redshift_obs_gauss_gen (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, NcmRNG *rng);
static gdouble _nc_galaxy_redshift_obs_gauss_window_mass (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, const gdouble obs_lo, const gdouble obs_hi);
static void _nc_galaxy_redshift_obs_gauss_get_true_z_lim (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, gdouble *z_min, gdouble *z_max);

static void
nc_galaxy_redshift_obs_gauss_class_init (NcGalaxyRedshiftObsGaussClass *klass)
{
  NcGalaxyRedshiftObsClass *gsdre_class = NC_GALAXY_REDSHIFT_OBS_CLASS (klass);
  GObjectClass *object_class                = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_redshift_obs_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian photometric-redshift observable model", "GaussRedshiftObs");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_model_class_check_params_info (model_class);

  gsdre_class->data_init      = &_nc_galaxy_redshift_obs_gauss_data_init;
  gsdre_class->eval           = &_nc_galaxy_redshift_obs_gauss_eval;
  gsdre_class->gen            = &_nc_galaxy_redshift_obs_gauss_gen;
  gsdre_class->window_mass    = &_nc_galaxy_redshift_obs_gauss_window_mass;
  gsdre_class->get_true_z_lim = &_nc_galaxy_redshift_obs_gauss_get_true_z_lim;
}

static void
_nc_galaxy_redshift_obs_gauss_ldata_read_row (NcGalaxyRedshiftObsData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyRedshiftObsGaussLData *ldata = (NcGalaxyRedshiftObsGaussLData *) data->ldata;

  ldata->zp     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_ZP, i, NULL);
  ldata->sigma0 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0, i, NULL);
}

static void
_nc_galaxy_redshift_obs_gauss_ldata_write_row (NcGalaxyRedshiftObsData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyRedshiftObsGaussLData *ldata = (NcGalaxyRedshiftObsGaussLData *) data->ldata;

  nc_galaxy_wl_obs_set (obs, NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_ZP, i, ldata->zp, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0, i, ldata->sigma0, NULL);
}

static void
_nc_galaxy_redshift_obs_gauss_ldata_required_columns (NcGalaxyRedshiftObsData *data, GList **columns)
{
  *columns = g_list_append (*columns, g_strdup (NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_ZP));
  *columns = g_list_append (*columns, g_strdup (NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0));
}

static void
_nc_galaxy_redshift_obs_gauss_data_init (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data)
{
  NcGalaxyRedshiftObsGaussLData *ldata = g_new0 (NcGalaxyRedshiftObsGaussLData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_redshift_obs_gauss_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_redshift_obs_gauss_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_redshift_obs_gauss_ldata_required_columns;
}

static gdouble
_nc_galaxy_redshift_obs_gauss_eval (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z)
{
  NcGalaxyRedshiftObsGaussLData *ldata = (NcGalaxyRedshiftObsGaussLData *) data->ldata;
  const gdouble sigmaz                     = ldata->sigma0 * (1.0 + z);
  const gdouble arg                        = (ldata->zp - z) / sigmaz;

  return exp (-0.5 * arg * arg) / (sqrt (2.0 * M_PI) * sigmaz);
}

static gdouble
_nc_galaxy_redshift_obs_gauss_gen (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, NcmRNG *rng)
{
  NcGalaxyRedshiftObsGaussLData *ldata = (NcGalaxyRedshiftObsGaussLData *) data->ldata;
  const gdouble sigmaz                     = ldata->sigma0 * (1.0 + z);

  ldata->zp = ncm_rng_gaussian_gen (rng, z, sigmaz);

  return ldata->zp;
}

static gdouble
_nc_galaxy_redshift_obs_gauss_window_mass (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, const gdouble obs_lo, const gdouble obs_hi)
{
  NcGalaxyRedshiftObsGaussLData *ldata = (NcGalaxyRedshiftObsGaussLData *) data->ldata;
  const gdouble sigmaz                          = ldata->sigma0 * (1.0 + z);

  /* P(zp in [obs_lo, obs_hi] | z) for zp ~ N(z, sigmaz). */
  return ncm_util_gaussian_integral (obs_lo, obs_hi, z, sigmaz);
}

static void
_nc_galaxy_redshift_obs_gauss_get_true_z_lim (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, gdouble *z_min, gdouble *z_max)
{
  NcGalaxyRedshiftObsGaussLData *ldata = (NcGalaxyRedshiftObsGaussLData *) data->ldata;
  /* Effective support around the point estimate zp: outside +-7 sigma the
   * Gaussian kernel is negligible. sigma_z = sigma0 (1 + z), estimated at the
   * upper edge z ~ zp + 7 sigma_z so the band comfortably covers the peak. */
  const gdouble sigma_max  = ldata->sigma0 * ((1.0 + 7.0 * ldata->sigma0) * (1.0 + ldata->zp));
  const gdouble half_width = 7.0 * sigma_max;

  *z_min = ldata->zp - half_width;
  *z_max = ldata->zp + half_width;
}

/**
 * nc_galaxy_redshift_obs_gauss_new:
 *
 * Creates a new #NcGalaxyRedshiftObsGauss.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftObsGauss.
 */
NcGalaxyRedshiftObsGauss *
nc_galaxy_redshift_obs_gauss_new (void)
{
  NcGalaxyRedshiftObsGauss *gsdreg = g_object_new (NC_TYPE_GALAXY_REDSHIFT_OBS_GAUSS,
                                                       NULL);

  return gsdreg;
}

/**
 * nc_galaxy_redshift_obs_gauss_ref:
 * @gsdreg: a #NcGalaxyRedshiftObsGauss
 *
 * Increases the reference count of @gsdreg by one.
 *
 * Returns: (transfer full): @gsdreg.
 */
NcGalaxyRedshiftObsGauss *
nc_galaxy_redshift_obs_gauss_ref (NcGalaxyRedshiftObsGauss *gsdreg)
{
  return g_object_ref (gsdreg);
}

/**
 * nc_galaxy_redshift_obs_gauss_free:
 * @gsdreg: a #NcGalaxyRedshiftObsGauss
 *
 * Decreases the reference count of @gsdreg by one.
 *
 */
void
nc_galaxy_redshift_obs_gauss_free (NcGalaxyRedshiftObsGauss *gsdreg)
{
  g_object_unref (gsdreg);
}

/**
 * nc_galaxy_redshift_obs_gauss_clear:
 * @gsdreg: a #NcGalaxyRedshiftObsGauss
 *
 * Decreases the reference count of *@gsdreg by one, and sets the pointer
 * *@gsdreg to NULL.
 *
 */
void
nc_galaxy_redshift_obs_gauss_clear (NcGalaxyRedshiftObsGauss **gsdreg)
{
  g_clear_object (gsdreg);
}

/**
 * nc_galaxy_redshift_obs_gauss_data_set:
 * @gsdreg: a #NcGalaxyRedshiftObsGauss
 * @data: a #NcGalaxyRedshiftObsData
 * @zp: the photometric-redshift point estimate
 * @sigma0: the scatter parameter, sigma_z = sigma0 (1 + z)
 *
 * Sets the per-galaxy photometric-redshift observation of @data.
 *
 */
void
nc_galaxy_redshift_obs_gauss_data_set (NcGalaxyRedshiftObsGauss *gsdreg, NcGalaxyRedshiftObsData *data, const gdouble zp, const gdouble sigma0)
{
  NcGalaxyRedshiftObsGaussLData *ldata = (NcGalaxyRedshiftObsGaussLData *) data->ldata;

  ldata->zp     = zp;
  ldata->sigma0 = sigma0;
}

/**
 * nc_galaxy_redshift_obs_gauss_data_get:
 * @gsdreg: a #NcGalaxyRedshiftObsGauss
 * @data: a #NcGalaxyRedshiftObsData
 * @zp: (out): the photometric-redshift point estimate
 * @sigma0: (out): the scatter parameter
 *
 * Gets the per-galaxy photometric-redshift observation of @data.
 *
 */
void
nc_galaxy_redshift_obs_gauss_data_get (NcGalaxyRedshiftObsGauss *gsdreg, NcGalaxyRedshiftObsData *data, gdouble *zp, gdouble *sigma0)
{
  NcGalaxyRedshiftObsGaussLData *ldata = (NcGalaxyRedshiftObsGaussLData *) data->ldata;

  *zp     = ldata->zp;
  *sigma0 = ldata->sigma0;
}
