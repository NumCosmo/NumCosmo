/***************************************************************************
 *            nc_data_cluster_pseudo_counts.h
 *
 *  Sun Apr 5 20:23:11 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_cluster_pseudo_counts.h
 * Copyright (C) 2015 Mariana Penna Lima <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_CLUSTER_PSEUDO_COUNTS_H_
#define _NC_DATA_CLUSTER_PSEUDO_COUNTS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>
#include <numcosmo/lss/nc_cluster_pseudo_counts.h>
#include <numcosmo/math/ncm_data.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS             (nc_data_cluster_pseudo_counts_get_type ())
#define NC_DATA_CLUSTER_PSEUDO_COUNTS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS, NcDataClusterPseudoCounts))
#define NC_DATA_CLUSTER_PSEUDO_COUNTS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS, NcDataClusterPseudoCountsClass))
#define NC_IS_DATA_CLUSTER_PSEUDO_COUNTS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS))
#define NC_IS_DATA_CLUSTER_PSEUDO_COUNTS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS))
#define NC_DATA_CLUSTER_PSEUDO_COUNTS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS, NcDataClusterPseudoCountsClass))

typedef struct _NcDataClusterPseudoCountsClass NcDataClusterPseudoCountsClass;
typedef struct _NcDataClusterPseudoCounts NcDataClusterPseudoCounts;

struct _NcDataClusterPseudoCountsClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

/**
 * NcDataClusterPseudoCountsObs:
 * @NC_DATA_CLUSTER_PSEUDO_COUNTS_Z: FIXME
 * @NC_DATA_CLUSTER_PSEUDO_COUNTS_MPL: FIXME
 * @NC_DATA_CLUSTER_PSEUDO_COUNTS_MCL: FIXME
 * @NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MPL: FIXME
 * @NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MCL: FIXME
 * 
 */
typedef enum _NcDataClusterPseudoCountsObs
{
  NC_DATA_CLUSTER_PSEUDO_COUNTS_Z = 0,
  NC_DATA_CLUSTER_PSEUDO_COUNTS_MPL,
  NC_DATA_CLUSTER_PSEUDO_COUNTS_MCL,
  NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MPL,
  NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MCL, /*< private >*/
  NC_DATA_CLUSTER_PSEUDO_COUNTS_LEN, /*< skip >*/
} NcDataClusterPseudoCountsObs;

struct _NcDataClusterPseudoCounts
{
  /*< private >*/
  NcmData parent_instance;
  NcmMatrix *obs;
  guint np;
  gchar *rnd_name;
};

GType nc_data_cluster_pseudo_counts_get_type (void) G_GNUC_CONST;

NcDataClusterPseudoCounts *nc_data_cluster_pseudo_counts_new (void);
NcDataClusterPseudoCounts *nc_data_cluster_pseudo_counts_new_from_file (const gchar *filename);
NcDataClusterPseudoCounts *nc_data_cluster_pseudo_counts_ref (NcDataClusterPseudoCounts *dcpc);
void nc_data_cluster_pseudo_counts_free (NcDataClusterPseudoCounts *dcpc);
void nc_data_cluster_pseudo_counts_clear (NcDataClusterPseudoCounts **dcpc);

void nc_data_cluster_pseudo_counts_set_obs (NcDataClusterPseudoCounts *dcpc, const NcmMatrix *m);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_PSEUDO_COUNTS_H_ */
