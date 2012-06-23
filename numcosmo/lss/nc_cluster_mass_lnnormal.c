/***************************************************************************
 *            nc_cluster_mass_lnnormal.c
 *
 *  Fri June 22 17:48:55 2012
 *  Copyright  2012  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
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
 * SECTION:nc_cluster_mass_lnnormal
 * @title: Cluster Abundance Mass Ln Normal Distribution
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

G_DEFINE_TYPE (NcClusterMassLnnormal, nc_cluster_mass_lnnormal, NC_TYPE_CLUSTER_MASS);

static gdouble
_nc_cluster_mass_lnnormal_dist_eval (NcClusterMass *clusterm, gdouble z, gdouble lnM, gdouble *lnM_obs, gdouble *lnM_obs_params)
{
  const gdouble lnMobs = lnM_obs[0];
  const gdouble lnM_bias = lnM_obs_params[NC_CLUSTER_MASS_LNNORMAL_BIAS];
  const gdouble sigma = lnM_obs_params[NC_CLUSTER_MASS_LNNORMAL_SIGMA];
  const gdouble sqrt2_sigma = M_SQRT2 * sigma;
  const gdouble x = (lnMobs - lnM - lnM_bias) / sqrt2_sigma;

  return M_2_SQRTPI / (2.0 * M_SQRT2) * exp (- x * x) / (sigma);
}


static void
nc_cluster_mass_lnnormal_init (NcClusterMassLnnormal *nc_cluster_mass_lnnormal)
{
}

static void
nc_cluster_mass_lnnormal_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_lnnormal_parent_class)->finalize (object);
}

static void
nc_cluster_mass_lnnormal_class_init (NcClusterMassLnnormalClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterMassClass *clusterm_class = NC_CLUSTER_MASS_CLASS (klass);
  clusterm_class->dist_eval = &_nc_cluster_mass_lnnormal_dist_eval;


  object_class->finalize = nc_cluster_mass_lnnormal_finalize;
}

