/***************************************************************************
 *            nc_cluster_mass_ascaso.h
 *
 *  Thu Jan 26 18:25:11 2017
 *  Copyright  2017  Mariana Penna Lima and Begoña Ascaso
 *  <pennalima@gmail.com>, <bego.ascaso.work@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima and Begoña Ascaso 2017 <pennalima@gmail.com>, <bego.ascaso.work@gmail.com>
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

#ifndef _NC_CLUSTER_MASS_ASCASO_H_
#define _NC_CLUSTER_MASS_ASCASO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass_richness.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_ASCASO (nc_cluster_mass_ascaso_get_type ())

G_DECLARE_FINAL_TYPE (NcClusterMassAscaso, nc_cluster_mass_ascaso, NC, CLUSTER_MASS_ASCASO, NcClusterMassRichness)

/**
 * NcClusterMassAscasoSParams:
 * @NC_CLUSTER_MASS_ASCASO_MU_P0: constant term (bias) in the mean log-richness
 * @NC_CLUSTER_MASS_ASCASO_MU_P1: linear mass coefficient in the mean log-richness
 * @NC_CLUSTER_MASS_ASCASO_MU_P2: redshift evolution coefficient in the mean log-richness
 * @NC_CLUSTER_MASS_ASCASO_SIGMA_P0: constant term (bias) in the standard deviation
 * @NC_CLUSTER_MASS_ASCASO_SIGMA_P1: linear mass coefficient in the standard deviation
 * @NC_CLUSTER_MASS_ASCASO_SIGMA_P2: redshift evolution coefficient in the standard deviation
 *
 * Parameters for the Ascaso cluster mass-richness relation.
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_ASCASO_SPARAMS >*/
{
  NC_CLUSTER_MASS_ASCASO_MU_P0 = NC_CLUSTER_MASS_RICHNESS_SPARAM_LEN,
  NC_CLUSTER_MASS_ASCASO_MU_P1,
  NC_CLUSTER_MASS_ASCASO_MU_P2,
  NC_CLUSTER_MASS_ASCASO_SIGMA_P0,
  NC_CLUSTER_MASS_ASCASO_SIGMA_P1,
  NC_CLUSTER_MASS_ASCASO_SIGMA_P2,
  /* < private > */
  NC_CLUSTER_MASS_ASCASO_SPARAM_LEN, /*< skip >*/
} NcClusterMassAscasoSParams;

#define NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P0  (3.19)
#define NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P1  (2.0 / M_LN10)
#define NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P2  (-0.7 / M_LN10)
#define NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P0  (0.33)
#define NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P1  (-0.08 / M_LN10)
#define NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P2  (0.0)
#define NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_ASCASO_H_ */

