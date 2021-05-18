/***************************************************************************
 *           nc_data_cluster_wl.h
 *
 *  Mon Jul 27 16:10:25 2020
 *  Copyright  2020  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_cluster_wl.h
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_DATA_CLUSTER_WL_H_
#define _NC_DATA_CLUSTER_WL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/lss/nc_galaxy_wl.h>
#include <numcosmo/math/ncm_data.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_WL             (nc_data_cluster_wl_get_type ())
#define NC_DATA_CLUSTER_WL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_CLUSTER_WL, NcDataClusterWL))
#define NC_DATA_CLUSTER_WL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_CLUSTER_WL, NcDataClusterWLClass))
#define NC_IS_DATA_CLUSTER_WL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_CLUSTER_WL))
#define NC_IS_DATA_CLUSTER_WL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_CLUSTER_WL))
#define NC_DATA_CLUSTER_WL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_CLUSTER_WL, NcDataClusterWLClass))

typedef struct _NcDataClusterWLClass NcDataClusterWLClass;
typedef struct _NcDataClusterWL NcDataClusterWL;
typedef struct _NcDataClusterWLPrivate NcDataClusterWLPrivate;

struct _NcDataClusterWLClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

/**
 * NcDataClusterWLObs:
 * @NC_DATA_CLUSTER_WL_ZCLUSTER: cluster redshift
 * @NC_DATA_CLUSTER_WL_GOBS: measured reduced shear
 * @NC_DATA_CLUSTER_WL_PZ: redshift distribution (photometric)
 *
 */
typedef enum _NcDataClusterWLObs
{
  NC_DATA_CLUSTER_WL_ZCLUSTER = 0,
  NC_DATA_CLUSTER_WL_GOBS,
  NC_DATA_CLUSTER_WL_PZ,
  /* < private > */
  NC_DATA_CLUSTER_WL_LEN, /*< skip >*/
} NcDataClusterWLObs;

struct _NcDataClusterWL
{
  /*< private >*/
  NcmData parent_instance;
  NcDataClusterWLPrivate *priv;
};

GType nc_data_cluster_wl_get_type (void) G_GNUC_CONST;

NcDataClusterWL *nc_data_cluster_wl_new (void);
NcDataClusterWL *nc_data_cluster_wl_new_from_file (const gchar *filename);
NcDataClusterWL *nc_data_cluster_wl_ref (NcDataClusterWL *dcwl);
void nc_data_cluster_wl_free (NcDataClusterWL *dcwl);
void nc_data_cluster_wl_clear (NcDataClusterWL **dcwl);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_WL_H_ */

