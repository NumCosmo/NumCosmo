/***************************************************************************
 *            nc_cluster_mass_selection.h
 *
 *  Thu Jan 26 18:25:11 2017
 *  Copyright  2017  Mariana Penna Lima and Begoña Selection
 *  <pennalima@gmail.com>, <bego.selection.work@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima and Begoña Selection 2017 <pennalima@gmail.com>, <bego.selection.work@gmail.com>
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

#ifndef _NC_CLUSTER_MASS_SELECTION_H_
#define _NC_CLUSTER_MASS_SELECTION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline2d_bicubic.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_SELECTION             (nc_cluster_mass_selection_get_type ())
#define NC_CLUSTER_MASS_SELECTION(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_SELECTION, NcClusterMassSelection))
#define NC_CLUSTER_MASS_SELECTION_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_SELECTION, NcClusterMassSelectionClass))
#define NC_IS_CLUSTER_MASS_SELECTION(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_SELECTION))
#define NC_IS_CLUSTER_MASS_SELECTION_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_SELECTION))
#define NC_CLUSTER_MASS_SELECTION_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_SELECTION, NcClusterMassSelectionClass))

typedef struct _NcClusterMassSelectionClass NcClusterMassSelectionClass;
typedef struct _NcClusterMassSelection NcClusterMassSelection;
typedef struct _NcClusterMassSelectionPrivate NcClusterMassSelectionPrivate;

/**
 * NcClusterMassSelectionSParams:
 * @NC_CLUSTER_MASS_SELECTION_MU_P0: bias of the mean
 * @NC_CLUSTER_MASS_SELECTION_MU_P1: slope on the mean
 * @NC_CLUSTER_MASS_SELECTION_MU_P2: redshift dependency on the mean
 * @NC_CLUSTER_MASS_SELECTION_SIGMA_P0: bias of the standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_SELECTION_SIGMA_P1: slope on the standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_SELECTION_SIGMA_P2: redshift dependency standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_SELECTION_CUT: cut in richness
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_SELECTION_SPARAMS >*/
{
  NC_CLUSTER_MASS_SELECTION_MU_P0,
  NC_CLUSTER_MASS_SELECTION_MU_P1,
  NC_CLUSTER_MASS_SELECTION_MU_P2,
  NC_CLUSTER_MASS_SELECTION_SIGMA_P0,
  NC_CLUSTER_MASS_SELECTION_SIGMA_P1,
  NC_CLUSTER_MASS_SELECTION_SIGMA_P2,
  NC_CLUSTER_MASS_SELECTION_CUT,


  /* < private > */
  NC_CLUSTER_MASS_SELECTION_SPARAM_LEN, /*< skip >*/
}   NcClusterMassSelectionSParams;

#define NC_CLUSTER_MASS_SELECTION_DEFAULT_MU_P0  (3.19)
#define NC_CLUSTER_MASS_SELECTION_DEFAULT_MU_P1  (2.0 / M_LN10)
#define NC_CLUSTER_MASS_SELECTION_DEFAULT_MU_P2  (-0.7 / M_LN10)
#define NC_CLUSTER_MASS_SELECTION_DEFAULT_SIGMA_P0  (0.33)
#define NC_CLUSTER_MASS_SELECTION_DEFAULT_SIGMA_P1  (-0.08 / M_LN10)
#define NC_CLUSTER_MASS_SELECTION_DEFAULT_SIGMA_P2  (0.0)
#define NC_CLUSTER_MASS_SELECTION_DEFAULT_CUT  (0.0)
#define NC_CLUSTER_MASS_SELECTION_DEFAULT_PARAMS_ABSTOL (0.0)

struct _NcClusterMassSelectionClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

struct _NcClusterMassSelection
{
  /*< private >*/
  NcClusterMass parent_instance;
  NcClusterMassSelectionPrivate *priv;
};

GType nc_cluster_mass_selection_get_type (void) G_GNUC_CONST;

void nc_cluster_mass_selection_set_enable_rejection (NcClusterMassSelection *selection, gboolean on);
void nc_cluster_mass_selection_set_purity(NcClusterMassSelection *selection, NcmSpline2dBicubic *purity);
void nc_cluster_mass_selection_set_completeness(NcClusterMassSelection *selection, NcmSpline2dBicubic *completeness);


gdouble nc_cluster_mass_selection_get_mean_richness (NcClusterMassSelection *selection, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_selection_get_std_richness (NcClusterMassSelection *selection, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_selection_get_cut (NcClusterMassSelection *selection, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_selection_get_mean (NcClusterMassSelection *selection, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_selection_get_std (NcClusterMassSelection *selection, gdouble lnM, gdouble z);
gboolean nc_cluster_mass_selection_get_enable_rejection (NcClusterMassSelection *selection);
NcmSpline2dBicubic * nc_cluster_mass_selection_get_purity(NcClusterMassSelection *selection);
NcmSpline2dBicubic * nc_cluster_mass_selection_get_completeness(NcClusterMassSelection *selection);

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_SELECTION_H_ */
