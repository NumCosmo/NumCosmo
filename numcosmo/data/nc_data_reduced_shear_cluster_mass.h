/***************************************************************************
 *           nc_data_reduced_shear_cluster_mass.h
 *
 *  Thu Mar 22 16:10:25 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_reduced_shear_cluster_mass.h
 * Copyright (C) 2018 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_H_
#define _NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_reduced_shear_cluster_mass.h>
#include <numcosmo/math/ncm_data.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_REDUCED_SHEAR_CLUSTER_MASS             (nc_data_reduced_shear_cluster_mass_get_type ())
#define NC_DATA_REDUCED_SHEAR_CLUSTER_MASS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_REDUCED_SHEAR_CLUSTER_MASS, NcDataReducedShearClusterMass))
#define NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_REDUCED_SHEAR_CLUSTER_MASS, NcDataReducedShearClusterMassClass))
#define NC_IS_DATA_REDUCED_SHEAR_CLUSTER_MASS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_REDUCED_SHEAR_CLUSTER_MASS))
#define NC_IS_DATA_REDUCED_SHEAR_CLUSTER_MASS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_REDUCED_SHEAR_CLUSTER_MASS))
#define NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_REDUCED_SHEAR_CLUSTER_MASS, NcDataReducedShearClusterMassClass))

typedef struct _NcDataReducedShearClusterMassClass NcDataReducedShearClusterMassClass;
typedef struct _NcDataReducedShearClusterMass NcDataReducedShearClusterMass;

struct _NcDataReducedShearClusterMassClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

/**
 * NcDataReducedShearClusterMassObs:
 * @NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_ZCLUSTER: cluster redshift
 * @NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_GOBS: measured reduced shear
 * @NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_PZ: (photometric) redshift distribution 
 * 
 */
typedef enum _NcDataReducedShearClusterMassObs
{
  NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_ZCLUSTER = 0,
  NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_GOBS,
  NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_PZ,
  /* < private > */
  NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_LEN, /*< skip >*/
} NcDataReducedShearClusterMassObs;

struct _NcDataReducedShearClusterMass
{
  /*< private >*/
  NcmData parent_instance;
  NcReducedShearClusterMass *rscm;
  NcmMatrix *g_obs;
  NcmMatrix *Pz;
  NcmMatrix *true_data;
  guint ngals;
  guint nzbins;
};

GType nc_data_reduced_shear_cluster_mass_get_type (void) G_GNUC_CONST;

NcDataReducedShearClusterMass *nc_data_reduced_shear_cluster_mass_new ();
NcDataReducedShearClusterMass *nc_data_reduced_shear_cluster_mass_new_from_file (const gchar *filename);
NcDataReducedShearClusterMass *nc_data_reduced_shear_cluster_mass_ref (NcDataReducedShearClusterMass *dcpc);
void nc_data_reduced_shear_cluster_mass_free (NcDataReducedShearClusterMass *dcpc);
void nc_data_reduced_shear_cluster_mass_clear (NcDataReducedShearClusterMass **dcpc);

void nc_data_reduced_shear_cluster_mass_set_ngalaxies (NcDataReducedShearClusterMass *dcpc, guint ngals);
guint nc_data_reduced_shear_cluster_mass_get_ngalaxies (NcDataReducedShearClusterMass *dcpc);
void nc_data_reduced_shear_cluster_mass_set_nzbins (NcDataReducedShearClusterMass *dcpc, guint nzbins);
guint nc_data_reduced_shear_cluster_mass_get_nzbins (NcDataReducedShearClusterMass *dcpc);
void nc_data_reduced_shear_cluster_mass_set_gobs (NcDataReducedShearClusterMass *dcpc, const NcmMatrix *m);
void nc_data_reduced_shear_cluster_mass_set_pz (NcDataReducedShearClusterMass *dcpc, const NcmMatrix *m);
void nc_data_reduced_shear_cluster_mass_set_true_data (NcDataReducedShearClusterMass *dcpc, const NcmMatrix *m);

G_END_DECLS

#endif /* _NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_H_ */
