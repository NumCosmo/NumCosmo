/***************************************************************************
 *            nc_cluster_photoz_gauss.c
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_photoz_gauss
 * @title: Individual Gaussian Photoz Cluster
 * @short_description: Gaussian photometric redshift
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_randist.h>

G_DEFINE_TYPE (NcClusterPhotozGauss, nc_cluster_photoz_gauss, NC_TYPE_CLUSTER_REDSHIFT);

enum
{
  PROP_0,
  PROP_PZ_MIN,
  PROP_PZ_MAX,
};

/**
 * nc_cluster_photoz_gauss_new:
 *
 * FIXME
 *
 * Returns: A new #NcClusterRedshift.
 */
NcClusterRedshift *
nc_cluster_photoz_gauss_new ()
{
  return g_object_new (NC_TYPE_CLUSTER_PHOTOZ_GAUSS, NULL);
}

static gdouble
_nc_cluster_photoz_gauss_dist_eval (NcClusterRedshift *clusterz, gdouble z, gdouble lnM, gdouble *z_obs, gdouble *z_obs_params)
{
  const gdouble z_bias = z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_BIAS];
  const gdouble sigma = z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_SIGMA];
  const gdouble sqrt2_sigma = M_SQRT2 * sigma;
  const gdouble y1 = (z_obs[0] - z - z_bias) / sqrt2_sigma;

  return M_2_SQRTPI / M_SQRT2 * exp (- y1 * y1) / (sigma * (1.0 + gsl_sf_erf ( z / sqrt2_sigma )));
}

static gboolean
_nc_cluster_photoz_gauss_resample (NcClusterRedshift *clusterz, gdouble z, gdouble lnM, gdouble *z_obs, gdouble *z_obs_params)
{
  NcClusterPhotozGauss *pzg = NC_CLUSTER_PHOTOZ_GAUSS (clusterz);
  gsl_rng *rng = ncm_get_rng ();
  gdouble sigma_z;

  z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_BIAS] = 0.0;
  z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_SIGMA] = 0.03 * (1.0 + z);

  sigma_z = z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_SIGMA];

  do {
	z_obs[0] = z + z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_BIAS] + gsl_ran_gaussian (rng, sigma_z);
  } while (z_obs[0] < 0);

  return (z_obs[0] <= pzg->pz_max) && (z_obs[0] >= pzg->pz_min);
}

static void
_nc_cluster_photoz_gauss_integ_limits (NcClusterRedshift *clusterz, gdouble *z_obs, gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper)
{
  const gdouble mean = z_obs[0] - z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_BIAS];
  const gdouble zl = GSL_MAX (mean - 10.0 * z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_SIGMA], 0.0);
  const gdouble zu = mean + 10.0 * z_obs_params[NC_CLUSTER_PHOTOZ_GAUSS_SIGMA];

  *z_lower = zl;
  *z_upper = zu;

  return;
}

static void
_nc_cluster_photoz_gauss_z_limits (NcClusterRedshift *clusterz, gdouble *z_lower, gdouble *z_upper)
{
  NcClusterPhotozGauss *pzg = NC_CLUSTER_PHOTOZ_GAUSS (clusterz);
  const gdouble zl = GSL_MAX (pzg->pz_min - 10.0 * 0.03 * (1.0 + pzg->pz_min), 0.0);
  const gdouble zu = pzg->pz_max + 10.0 * 0.03 * (1.0 + pzg->pz_max);

  *z_lower = zl;
  *z_upper = zu;

  return;
}

guint _nc_cluster_photoz_gauss_obs_len (NcClusterRedshift *clusterz) { return 1; }
guint _nc_cluster_photoz_gauss_obs_params_len (NcClusterRedshift *clusterz) { return 2; }

static void
nc_cluster_photoz_gauss_init (NcClusterPhotozGauss *pzg)
{
  pzg->pz_max = 0.0;
  pzg->pz_min = 0.0;
}

static void
_nc_cluster_photoz_gauss_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_photoz_gauss_parent_class)->finalize (object);
}

static void
_nc_cluster_photoz_gauss_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcClusterPhotozGauss *pzg = NC_CLUSTER_PHOTOZ_GAUSS (object);
  g_return_if_fail (NC_IS_CLUSTER_PHOTOZ_GAUSS (object));

  switch (prop_id)
  {
    case PROP_PZ_MIN:
      pzg->pz_min = g_value_get_double (value);
      break;
	case PROP_PZ_MAX:
      pzg->pz_max = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_photoz_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterPhotozGauss *pzg = NC_CLUSTER_PHOTOZ_GAUSS (object);
  g_return_if_fail (NC_IS_CLUSTER_PHOTOZ_GAUSS (object));

  switch (prop_id)
  {
    case PROP_PZ_MIN:
      g_value_set_double (value, pzg->pz_min);
      break;
	case PROP_PZ_MAX:
      g_value_set_double (value, pzg->pz_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_photoz_gauss_class_init (NcClusterPhotozGaussClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterRedshiftClass* parent_class = NC_CLUSTER_REDSHIFT_CLASS (klass);

  parent_class->dist_eval      = &_nc_cluster_photoz_gauss_dist_eval;
  parent_class->resample       = &_nc_cluster_photoz_gauss_resample;
  parent_class->integ_limits   = &_nc_cluster_photoz_gauss_integ_limits;
  parent_class->z_limits       = &_nc_cluster_photoz_gauss_z_limits;
  parent_class->obs_len        = &_nc_cluster_photoz_gauss_obs_len;
  parent_class->obs_params_len = &_nc_cluster_photoz_gauss_obs_params_len;

  object_class->finalize = _nc_cluster_photoz_gauss_finalize;
  object_class->set_property = _nc_cluster_photoz_gauss_set_property;
  object_class->get_property = _nc_cluster_photoz_gauss_get_property;

  /**
   * NcClusterPhotozGauss:pz_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_PZ_MIN,
                                   g_param_spec_double ("pz-min",
                                                        NULL,
                                                        "Minimum photoz",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterPhotozGauss:pz_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_PZ_MAX,
                                   g_param_spec_double ("pz-max",
                                                        NULL,
                                                        "Maximum photoz",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}
