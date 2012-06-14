/***************************************************************************
 *            nc_cluster_photoz_gauss_ind.c
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
 * SECTION:nc_cluster_photoz_gauss_ind
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

G_DEFINE_TYPE (NcClusterPhotozGaussInd, nc_cluster_photoz_gauss_ind, NC_TYPE_CLUSTER_PHOTOZ);

/**
 * nc_cluster_photoz_gauss_ind_new:
 *
 * FIXME
 *
 * Returns: A new #NcClusterPhotoz.
 */
NcClusterPhotoz *
nc_cluster_photoz_gauss_ind_new ()
{
  return g_object_new (NC_TYPE_CLUSTER_PHOTOZ_GAUSS_IND, NULL);
}

static gdouble
_nc_cluster_photoz_gauss_ind_dist_eval (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble z_real)
{
  const gdouble z_bias = pz_data[0];
  const gdouble sigma0 = pz_data[1];
  const gdouble sqrt2_sigma = M_SQRT2 * sigma0 * (1.0 + z_real);
  const gdouble y1 = (z_photo - z_real - z_bias) / sqrt2_sigma;

  return M_2_SQRTPI / M_SQRT2 * exp (- y1 * y1) / (sigma0 * (1.0 + z_real) * (1.0 + gsl_sf_erf (z_real/sqrt2_sigma)));
}

static void
_nc_cluster_photoz_gauss_ind_resample (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_real, gdouble *result)
{
  gsl_rng *rng = ncm_get_rng ();
  const gdouble sigma_z = pz_data[1] * (1.0 + z_real);
  gdouble z_photo = z_real + pz_data[0] + gsl_ran_gaussian (rng, sigma_z);

  *result = z_photo;

  return;
}

static void
_nc_cluster_photoz_gauss_ind_integ_limits (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble *z_lower, gdouble *z_upper)
{
  const gdouble zl = GSL_MAX (z_photo - 10.0 * pz_data[1] * (1.0 + z_photo), 0.0);
  const gdouble zu = z_photo + 10.0 * pz_data[1] * (1.0 + z_photo);

  *z_lower = zl;
  *z_upper = zu;

  return;
}

static void
nc_cluster_photoz_gauss_ind_init (NcClusterPhotozGaussInd *nc_cluster_photoz_gauss_ind)
{
  /* TODO: Add initialization code here */
}

static void
_nc_cluster_photoz_gauss_ind_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_cluster_photoz_gauss_ind_parent_class)->finalize (object);
}

static void
nc_cluster_photoz_gauss_ind_class_init (NcClusterPhotozGaussIndClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterPhotozClass* parent_class = NC_CLUSTER_PHOTOZ_CLASS (klass);

  parent_class->dist_eval = &_nc_cluster_photoz_gauss_ind_dist_eval;
  parent_class->resample = &_nc_cluster_photoz_gauss_ind_resample;
  parent_class->integ_limits = &_nc_cluster_photoz_gauss_ind_integ_limits;

  object_class->finalize = _nc_cluster_photoz_gauss_ind_finalize;
}

