/***************************************************************************
 *            nc_data_cluster_abundance.h
 *
 *  Tue Apr  6 01:12:58 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_DATA_NCUSTER_ABUNDANCE_H
#define _NC_DATA_NCUSTER_ABUNDANCE_H

#include <gsl/gsl_histogram2d.h>

G_BEGIN_DECLS

typedef struct _NcDataClusterAbundance NcDataClusterAbundance;

struct _NcDataClusterAbundance
{
  /*< private >*/
  NcClusterRedshift *z;
  NcClusterMass *m;
  NcmVector *lnM_true;
  NcmVector *z_true;
  NcmMatrix *z_obs;
  NcmMatrix *z_obs_params;
  NcmMatrix *lnM_obs;
  NcmMatrix *lnM_obs_params;
  gdouble area_survey;
  glong np;
  gdouble log_np_fac;
  gboolean use_true_data;
  gsl_histogram2d *completeness;
  gsl_histogram2d *purity;
  gsl_histogram2d *sd_lnM;
  gboolean fiducial;
  gulong seed;
  gchar *rnd_name;
};

/**
 * NcDataClusterAbundanceId:
 * @NC_DATA_CLUSTER_ABUNDANCE_FIT: FIXME
 * @NC_DATA_CLUSTER_ABUNDANCE_TXT: FIXME
 * @NC_DATA_CLUSTER_ABUNDANCE_SAMPLING: FIXME
 */
typedef enum _NcDataClusterAbundanceId
{
  NC_DATA_CLUSTER_ABUNDANCE_FIT,
  NC_DATA_CLUSTER_ABUNDANCE_TXT,
  NC_DATA_CLUSTER_ABUNDANCE_SAMPLING, /*< private >*/
  NC_DATA_CLUSTER_ABUNDANCE_NSAMPLES, /*< skip >*/
} NcDataClusterAbundanceId;

NcData *nc_data_cluster_abundance_binned_new (NcClusterAbundance *cad);
void nc_data_cluster_abundance_binned_init_from_text_file_gkey (NcData *data, gboolean obs, gchar *filename);
void nc_data_cluster_abundance_binned_init_from_sampling (NcData *data, NcmMSet *mset, NcmVector *nodes, gboolean obs, gdouble area_survey, gdouble lnMi, gdouble lnMf, gdouble photoz_sigma0, gdouble photoz_bias, gdouble lnM_sigma0, gdouble lnM_bias);
void nc_data_cluster_abundance_binned_save (NcData *data, gchar *filename);
NcmMSetFunc *nc_data_cluster_abundance_binned_new_function (NcClusterAbundance *cad);

NcData *nc_data_cluster_abundance_binned_lnM_z_new (NcClusterAbundance *cad);
NcData *nc_data_cluster_abundance_unbinned_new (NcClusterAbundance *cad);
void nc_data_cluster_abundance_true_data (NcData *data, gboolean use_true_data);
void nc_data_cluster_abundance_unbinned_init_from_sampling (NcData *data, NcmMSet *mset, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble area_survey);
NcData *nc_data_cluster_abundance_unbinned_bin_data (NcData *ca_unbinned, gsl_vector *nodes);
gsl_histogram2d *nc_data_cluster_abundance_hist_lnM_z (NcData *ca_unbinned, gsl_vector *lnM_nodes, gsl_vector *z_nodes);

#ifdef NUMCOSMO_HAVE_CFITSIO
void nc_cluster_abundance_catalog_save (NcData *data, gchar *filename, gboolean overwrite);
void nc_cluster_abundance_catalog_load (NcData *data, gchar *filename);
#endif /* NUMCOSMO_HAVE_CFITSIO */

G_END_DECLS

#endif /* _NC_DATA_NCUSTER_ABUNDANCE_H */
