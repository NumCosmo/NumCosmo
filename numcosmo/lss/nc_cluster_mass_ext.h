/***************************************************************************
 *            nc_cluster_mass_ext.h
 *
 *  Tue Dez 02 18:25:11 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2025 <vitenti@uel.br>
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

#ifndef _NC_CLUSTER_MASS_EXT_H_
#define _NC_CLUSTER_MASS_EXT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass_richness.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_EXT (nc_cluster_mass_ext_get_type ())

G_DECLARE_FINAL_TYPE (NcClusterMassExt, nc_cluster_mass_ext, NC, CLUSTER_MASS_EXT, NcClusterMassRichness)

/**
 * NcClusterMassExtSParams:
 * @NC_CLUSTER_MASS_EXT_MU_P0: constant term (bias) in the mean log-richness
 * @NC_CLUSTER_MASS_EXT_MU_P1: linear mass coefficient in the mean log-richness
 * @NC_CLUSTER_MASS_EXT_MU_P2: quadratic mass coefficient in the mean log-richness
 * @NC_CLUSTER_MASS_EXT_MU_P3: cubic mass coefficient in the mean log-richness
 * @NC_CLUSTER_MASS_EXT_SIGMA_P0: constant term in ln(sigma)
 * @NC_CLUSTER_MASS_EXT_SIGMA_P1: linear coefficient of mu in ln(sigma)
 * @NC_CLUSTER_MASS_EXT_SIGMA_P2: quadratic coefficient of mu in ln(sigma)
 *
 * Parameters for the extended cluster mass-richness relation.
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_EXT_SPARAMS >*/
{
  NC_CLUSTER_MASS_EXT_MU_P0 = NC_CLUSTER_MASS_RICHNESS_SPARAM_LEN,
  NC_CLUSTER_MASS_EXT_MU_P1,
  NC_CLUSTER_MASS_EXT_MU_P2,
  NC_CLUSTER_MASS_EXT_MU_P3,
  NC_CLUSTER_MASS_EXT_SIGMA_P0,
  NC_CLUSTER_MASS_EXT_SIGMA_P1,
  NC_CLUSTER_MASS_EXT_SIGMA_P2,
  /* < private > */
  NC_CLUSTER_MASS_EXT_SPARAM_LEN, /*< skip >*/
} NcClusterMassExtSParams;

#define NC_CLUSTER_MASS_EXT_DEFAULT_MU_P0  (4.0)
#define NC_CLUSTER_MASS_EXT_DEFAULT_MU_P1  (1.0)
#define NC_CLUSTER_MASS_EXT_DEFAULT_MU_P2  (0.1)
#define NC_CLUSTER_MASS_EXT_DEFAULT_MU_P3  (0.01)
#define NC_CLUSTER_MASS_EXT_DEFAULT_SIGMA_P0  (-0.3)
#define NC_CLUSTER_MASS_EXT_DEFAULT_SIGMA_P1  (-0.3)
#define NC_CLUSTER_MASS_EXT_DEFAULT_SIGMA_P2  (0.1)
#define NC_CLUSTER_MASS_EXT_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_EXT_H_ */

