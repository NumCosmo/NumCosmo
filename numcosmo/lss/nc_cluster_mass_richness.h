/***************************************************************************
 *            nc_cluster_mass_richness.h
 *
 *  Tue May 27 14:05:00 2025
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

#ifndef _NC_CLUSTER_MASS_RICHNESS_H_
#define _NC_CLUSTER_MASS_RICHNESS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_RICHNESS (nc_cluster_mass_richness_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcClusterMassRichness, nc_cluster_mass_richness, NC, CLUSTER_MASS_RICHNESS, NcClusterMass)

/**
 * NcClusterMassRichnessSParams:
 * @NC_CLUSTER_MASS_RICHNESS_CUT: cut in richness
 *
 * Scalar parameters for the cluster mass-richness relation.
 *
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_RICHNESS_SPARAMS >*/
{
  NC_CLUSTER_MASS_RICHNESS_CUT,
  /* < private > */
  NC_CLUSTER_MASS_RICHNESS_SPARAM_LEN, /*< skip >*/
} NcClusterMassRichnessSParams;

#define NC_CLUSTER_MASS_RICHNESS_DEFAULT_CUT  (0.0)
#define NC_CLUSTER_MASS_RICHNESS_DEFAULT_PARAMS_ABSTOL (0.0)

struct _NcClusterMassRichnessClass
{
  /*< private >*/
  NcClusterMassClass parent_class;

  gdouble (*mu) (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
  gdouble (*sigma) (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
  void (*mu_sigma) (NcClusterMassRichness *mr, gdouble lnM, gdouble z, gdouble *mu, gdouble *sigma);

  /* Padding to allow adding up to 18 virtual functions without breaking ABI. */
  gpointer padding[15];
};

gdouble nc_cluster_mass_richness_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_richness_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
void nc_cluster_mass_richness_mu_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z, gdouble *mu, gdouble *sigma);
gdouble nc_cluster_mass_richness_get_cut (NcClusterMassRichness *mr);

void nc_cluster_mass_richness_set_sample_full_dist (NcClusterMassRichness *mr, gboolean on);
gboolean nc_cluster_mass_richness_get_sample_full_dist (NcClusterMassRichness *mr);

gdouble nc_cluster_mass_richness_get_mean (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_richness_get_std (NcClusterMassRichness *mr, gdouble lnM, gdouble z);

gdouble nc_cluster_mass_richness_compute_truncated_mean (NcClusterMassRichness *mr, gdouble lnR_mean, gdouble lnR_sigma);
gdouble nc_cluster_mass_richness_compute_truncated_std (NcClusterMassRichness *mr, gdouble lnR_mean, gdouble lnR_sigma);

gdouble nc_cluster_mass_richness_lnM0 (NcClusterMassRichness *mr);
gdouble nc_cluster_mass_richness_ln1pz0 (NcClusterMassRichness *mr);

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_RICHNESS_H_ */

