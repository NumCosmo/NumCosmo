/***************************************************************************
 *            nc_cluster_photoz_gauss_global.c
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
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_photoz_gauss_global
 * @title: Global Gaussian Photoz Cluster
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
//#include <gsl/gsl_randist.h>

G_DEFINE_TYPE (NcClusterPhotozGaussGlobal, nc_cluster_photoz_gauss_global, NC_TYPE_CLUSTER_PHOTOZ);

enum
{
  PROP_0,
  PROP_Z_BIAS,
  PROP_SIGMA0
};

/**
 * nc_cluster_photoz_gauss_global_new:
 * @z_bias: FIXME
 * @sigma0: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcClusterPhotoz.
 */
NcClusterPhotoz *
nc_cluster_photoz_gauss_global_new (gdouble z_bias, gdouble sigma0)
{
  return g_object_new (NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL,
                       "z-bias", z_bias,
                       "sigma0", sigma0,
                       NULL);
}

static gdouble
_nc_cluster_photoz_gauss_global_dist_eval (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble z_real)
{
  NcClusterPhotozGaussGlobal *pz_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (photo);
  const gdouble sqrt2_sigma = M_SQRT2 * pz_global->sigma0 * (1.0 + z_real);
  const gdouble y1 = (z_photo - z_real - pz_global->z_bias) / sqrt2_sigma;

  return M_2_SQRTPI / M_SQRT2 * exp (- y1 * y1) / (pz_global->sigma0 * (1.0 + z_real) * (1.0 + gsl_sf_erf (z_real/sqrt2_sigma)));
}

static void
_nc_cluster_photoz_gauss_global_resample (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_real, gdouble *result)
{
  NcClusterPhotozGaussGlobal *pz_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (photo);
  const gdouble sigma_z = pz_global->sigma0 * (1.0 + z_real);
  gsl_rng *rng = ncm_get_rng ();
  gdouble z_photo = z_real + pz_global->z_bias + gsl_ran_gaussian (rng, sigma_z);

  *result = z_photo;

  return;
}

static void
_nc_cluster_photoz_gauss_global_integ_limits (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble *z_lower, gdouble *z_upper)
{
  NcClusterPhotozGaussGlobal *pz_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (photo);
  const gdouble zl = GSL_MAX (z_photo - 10.0 * pz_global->sigma0 * (1.0 + z_photo), 0.0);
  const gdouble zu = z_photo + 10.0 * pz_global->sigma0 * (1.0 + z_photo);

  *z_lower = zl;
  *z_upper = zu;

  return;
}

/**
 * nc_cluster_photoz_gauss_global_set_z_bias:
 * @pz_global: a #NcClusterPhotozGaussGlobal.
 * @z_bias: value of #NcClusterPhotozGaussGlobal:z-bias.
 *
 * Sets the value @z_bias to the #NcClusterPhotozGaussGlobal:z-bias property.
 *
 */
void
nc_cluster_photoz_gauss_global_set_z_bias (NcClusterPhotozGaussGlobal *pz_global, gdouble z_bias)
{
  pz_global->z_bias = z_bias;
}

/**
 * nc_cluster_photoz_gauss_global_get_z_bias:
 * @pz_global: a #NcClusterPhotozGaussGlobal.
 *
 * Returns: the value of #NcClusterPhotozGaussGlobal:z-bias property.
 */
gdouble
nc_cluster_photoz_gauss_global_get_z_bias (const NcClusterPhotozGaussGlobal *pz_global)
{
  return pz_global->z_bias;
}

/**
 * nc_cluster_photoz_gauss_global_set_sigma0:
 * @pz_global: a #NcClusterPhotozGaussGlobal.
 * @sigma0: value of #NcClusterPhotozGaussGlobal:sigma0.
 *
 * Sets the value @sigma0 to the #NcClusterPhotozGaussGlobal:sigma0 property.
 *
 */
void
nc_cluster_photoz_gauss_global_set_sigma0 (NcClusterPhotozGaussGlobal *pz_global, gdouble sigma0)
{
  pz_global->sigma0 = sigma0;
}

/**
 * nc_cluster_photoz_gauss_global_get_sigma0:
 * @pz_global: a #NcClusterPhotozGaussGlobal.
 *
 * Returns: the value of #NcClusterPhotozGaussGlobal:sigma0 property.
 */
gdouble
nc_cluster_photoz_gauss_global_get_sigma0 (const NcClusterPhotozGaussGlobal *pz_global)
{
  return pz_global->sigma0;
}

static void
nc_cluster_photoz_gauss_global_init (NcClusterPhotozGaussGlobal *pz_global)
{
  pz_global->z_bias = 0.0;
  pz_global->sigma0 = 0.01;
}

static void
_nc_cluster_photoz_gauss_global_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_cluster_photoz_gauss_global_parent_class)->finalize (object);
}

static void
_nc_cluster_photoz_gauss_global_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcClusterPhotozGaussGlobal *pz_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (object);
  g_return_if_fail (NC_IS_CLUSTER_PHOTOZ_GAUSS_GLOBAL (object));

  switch (prop_id)
  {
    case PROP_Z_BIAS:
      pz_global->z_bias = g_value_get_double (value);
      break;
	case PROP_SIGMA0:
      pz_global->sigma0 = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_photoz_gauss_global_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterPhotozGaussGlobal *pz_global = NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL (object);
  g_return_if_fail (NC_IS_CLUSTER_PHOTOZ_GAUSS_GLOBAL (object));

  switch (prop_id)
  {
    case PROP_Z_BIAS:
      g_value_set_double (value, pz_global->z_bias);
      break;
	case PROP_SIGMA0:
      g_value_set_double (value, pz_global->sigma0);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_photoz_gauss_global_class_init (NcClusterPhotozGaussGlobalClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterPhotozClass* parent_class = NC_CLUSTER_PHOTOZ_CLASS (klass);

  parent_class->dist_eval = &_nc_cluster_photoz_gauss_global_dist_eval;
  parent_class->resample = &_nc_cluster_photoz_gauss_global_resample;
  parent_class->integ_limits = &_nc_cluster_photoz_gauss_global_integ_limits;

  object_class->finalize = _nc_cluster_photoz_gauss_global_finalize;
  object_class->set_property = _nc_cluster_photoz_gauss_global_set_property;
  object_class->get_property = _nc_cluster_photoz_gauss_global_get_property;

  /**
   * NcClusterPhotozGaussGlobal:z_bias:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_BIAS,
                                   g_param_spec_double ("z-bias",
                                                        NULL,
                                                        "z-bias",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterPhotozGaussGlobal:sigma0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGMA0,
                                   g_param_spec_double ("sigma0",
                                                        NULL,
                                                        "sigma0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

