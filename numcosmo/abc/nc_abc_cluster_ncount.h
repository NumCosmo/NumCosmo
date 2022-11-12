/***************************************************************************
 *            nc_abc_cluster_ncount.h
 *
 *  Mon October 13 13:29:08 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_abc_cluster_ncount.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_ABC_CLUSTER_NCOUNT_H_
#define _NC_ABC_CLUSTER_NCOUNT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_abc.h>
#include <numcosmo/data/nc_data_cluster_ncount.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_histogram2d.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NC_TYPE_ABC_CLUSTER_NCOUNT             (nc_abc_cluster_ncount_get_type ())
#define NC_ABC_CLUSTER_NCOUNT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_ABC_CLUSTER_NCOUNT, NcABCClusterNCount))
#define NC_ABC_CLUSTER_NCOUNT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_ABC_CLUSTER_NCOUNT, NcABCClusterNCountClass))
#define NC_IS_ABC_CLUSTER_NCOUNT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_ABC_CLUSTER_NCOUNT))
#define NC_IS_ABC_CLUSTER_NCOUNT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_ABC_CLUSTER_NCOUNT))
#define NC_ABC_CLUSTER_NCOUNT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_ABC_CLUSTER_NCOUNT, NcABCClusterNCountClass))

typedef struct _NcABCClusterNCountClass NcABCClusterNCountClass;
typedef struct _NcABCClusterNCount NcABCClusterNCount;

struct _NcABCClusterNCountClass
{
  /*< private >*/
  NcmABCClass parent_class;

};

/**
 * NcABCClusterNCountSummary:
 * @NC_ABC_CLUSTER_NCOUNT_SUMMARY_BIN_UNIFORM: FIXME
 * @NC_ABC_CLUSTER_NCOUNT_SUMMARY_BIN_QUANTILE: FIXME
 * @NC_ABC_CLUSTER_NCOUNT_SUMMARY_BIN_NODES: FIXME
 * @NC_ABC_CLUSTER_NCOUNT_SUMMARY_GAUSS_RBF: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcABCClusterNCountSummary 
{
  NC_ABC_CLUSTER_NCOUNT_SUMMARY_BIN_UNIFORM = 0,
  NC_ABC_CLUSTER_NCOUNT_SUMMARY_BIN_QUANTILE,
  NC_ABC_CLUSTER_NCOUNT_SUMMARY_BIN_NODES,  
  NC_ABC_CLUSTER_NCOUNT_SUMMARY_GAUSS_RBF,
  /* < private > */
  NC_ABC_CLUSTER_NCOUNT_SUMMARY_NTYPES, /*< skip >*/
} NcABCClusterNCountSummary;

/**
 * NcABCClusterNCountEpsilonUpdate:
 * @NC_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE_UNIFORM: FIXME
 * @NC_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE_QUANTILE: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcABCClusterNCountEpsilonUpdate 
{
  NC_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE_UNIFORM = 0, 
  NC_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE_QUANTILE, 
  /* < private > */
  NC_ABC_CLUSTER_NCOUNT_EPSILON_UPDATE_NTYPE, /*< skip >*/
}NcABCClusterNCountEpsilonUpdate;

struct _NcABCClusterNCount
{
  /*< private >*/
  NcmABC parent_instance;
  gboolean scale_cov;
  gsl_histogram2d *data_summary;
  NcDataClusterNCount *ncount;
  gdouble data_total;
  NcmMatrix *covar;
  NcmVector *quantiles;
  NcmVector *z_nodes;
  NcmVector *lnM_nodes;
  guint z_bins;
  guint lnM_bins;
  gdouble rbf_scale;
  NcmStatsVec *z_lnM_stats;
  gdouble sigma_z, sigma_lnM, rho;
  gdouble epsilon_update;
  NcABCClusterNCountSummary s_type;
  NcABCClusterNCountEpsilonUpdate uptype;
};

GType nc_abc_cluster_ncount_get_type (void) G_GNUC_CONST;

NcABCClusterNCount *nc_abc_cluster_ncount_new (NcmMSet *mset, NcmMSetTransKern *prior, NcmDataset *dset);
void nc_abc_cluster_ncount_set_scale_cov (NcABCClusterNCount *abcnc, gboolean on);

void nc_abc_cluster_ncount_set_epsilon_update (NcABCClusterNCount *abcnc, gdouble q);
void nc_abc_cluster_ncount_set_bin_uniform (NcABCClusterNCount *abcnc, guint z_bins, guint lnM_bins);
void nc_abc_cluster_ncount_set_bin_quantile (NcABCClusterNCount *abcnc, NcmVector *quantiles);
void nc_abc_cluster_ncount_set_bin_nodes (NcABCClusterNCount *abcnc, NcmVector *z_nodes, NcmVector *lnM_nodes);

G_END_DECLS

#endif /* _NC_ABC_CLUSTER_NCOUNT_H_ */

