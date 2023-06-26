/***************************************************************************
 *           nc_data_cluster_wll.h
 *
 *  Tue Jun 15 16:24:17 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_data_cluster_wll.h
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_DATA_CLUSTER_WLL_H_
#define _NC_DATA_CLUSTER_WLL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/galaxy/nc_galaxy_wl_likelihood.h>
#include <numcosmo/math/ncm_data.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_WLL             (nc_data_cluster_wll_get_type ())
#define NC_DATA_CLUSTER_WLL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_CLUSTER_WLL, NcDataClusterWLL))
#define NC_DATA_CLUSTER_WLL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_CLUSTER_WLL, NcDataClusterWLLClass))
#define NC_IS_DATA_CLUSTER_WLL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_CLUSTER_WLL))
#define NC_IS_DATA_CLUSTER_WLL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_CLUSTER_WLL))
#define NC_DATA_CLUSTER_WLL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_CLUSTER_WLL, NcDataClusterWLLClass))

typedef struct _NcDataClusterWLLClass NcDataClusterWLLClass;
typedef struct _NcDataClusterWLL NcDataClusterWLL;
typedef struct _NcDataClusterWLLPrivate NcDataClusterWLLPrivate;

struct _NcDataClusterWLLClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

/**
 * NcDataClusterWLLObs:
 * @NC_DATA_CLUSTER_WLL_ZCLUSTER: cluster redshift
 * @NC_DATA_CLUSTER_WLL_GOBS: measured reduced shear
 * @NC_DATA_CLUSTER_WLL_PZ: redshift distribution (photometric)
 *
 */
typedef enum _NcDataClusterWLLObs
{
  NC_DATA_CLUSTER_WLL_ZCLUSTER = 0,
  NC_DATA_CLUSTER_WLL_GOBS,
  NC_DATA_CLUSTER_WLL_PZ,
  /* < private > */
  NC_DATA_CLUSTER_WLL_LEN, /*< skip >*/
} NcDataClusterWLLObs;

struct _NcDataClusterWLL
{
  /*< private >*/
  NcmData parent_instance;
  NcDataClusterWLLPrivate *priv;
};

GType nc_data_cluster_wll_get_type (void) G_GNUC_CONST;

NcDataClusterWLL *nc_data_cluster_wll_new (void);
NcDataClusterWLL *nc_data_cluster_wll_new_from_file (const gchar *filename);
NcDataClusterWLL *nc_data_cluster_wll_ref (NcDataClusterWLL *dcwll);
void nc_data_cluster_wll_free (NcDataClusterWLL *dcwll);
void nc_data_cluster_wll_clear (NcDataClusterWLL **dcwll);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_WLL_H_ */

