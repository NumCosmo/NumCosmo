/***************************************************************************
 *            nc_data_cluster_ncount.h
 *
 *  Tue Apr  6 01:12:58 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_CLUSTER_NCOUNT_H_
#define _NC_DATA_CLUSTER_NCOUNT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/lss/nc_cluster_redshift.h>
#include <numcosmo/lss/nc_cluster_mass.h>
#include <numcosmo/lss/nc_cluster_abundance.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_NCOUNT             (nc_data_cluster_ncount_get_type ())
#define NC_DATA_CLUSTER_NCOUNT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_CLUSTER_NCOUNT, NcDataClusterNCount))
#define NC_DATA_CLUSTER_NCOUNT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_CLUSTER_NCOUNT, NcDataClusterNCountClass))
#define NC_IS_DATA_CLUSTER_NCOUNT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_CLUSTER_NCOUNT))
#define NC_IS_DATA_CLUSTER_NCOUNT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_CLUSTER_NCOUNT))
#define NC_DATA_CLUSTER_NCOUNT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_CLUSTER_NCOUNT, NcDataClusterNCountClass))

typedef struct _NcDataClusterNCountClass NcDataClusterNCountClass;
typedef struct _NcDataClusterNCount NcDataClusterNCount;


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

struct _NcDataClusterNCountClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

struct _NcDataClusterNCount
{
  /*< private >*/
  NcmData parent_instance;
  NcClusterAbundance *cad;
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

GType nc_data_cluster_ncount_get_type (void) G_GNUC_CONST;

NcmData *nc_data_cluster_ncount_new (NcClusterAbundance *cad);

NcDataClusterNCount *nc_data_cluster_ncount_ref (NcDataClusterNCount *ncount);
void nc_data_cluster_ncount_free (NcDataClusterNCount *ncount);

NcmData *nc_data_cluster_ncount_binned_new (NcClusterAbundance *cad);
void nc_data_cluster_ncount_binned_init_from_text_file_gkey (NcmData *data, gboolean obs, gchar *filename);
void nc_data_cluster_ncount_binned_init_from_sampling (NcmData *data, NcmMSet *mset, NcmVector *nodes, gboolean obs, gdouble area_survey, gdouble lnMi, gdouble lnMf, gdouble photoz_sigma0, gdouble photoz_bias, gdouble lnM_sigma0, gdouble lnM_bias);
void nc_data_cluster_ncount_binned_save (NcmData *data, gchar *filename);
NcmMSetFunc *nc_data_cluster_ncount_binned_new_function (NcClusterAbundance *cad);

NcmData *nc_data_cluster_ncount_binned_lnM_z_new (NcClusterAbundance *cad);
void nc_data_cluster_ncount_true_data (NcmData *data, gboolean use_true_data);
void nc_data_cluster_ncount_init_from_sampling (NcmData *data, NcmMSet *mset, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble area_survey);
NcmData *nc_data_cluster_ncount_bin_data (NcmData *data, gsl_vector *nodes);
gsl_histogram2d *nc_data_cluster_ncount_hist_lnM_z (NcmData *data, gsl_vector *lnM_nodes, gsl_vector *z_nodes);

void nc_data_cluster_ncount_print (NcmData *data, NcHICosmo *cosmo, FILE *out, gchar *header);

#ifdef NUMCOSMO_HAVE_CFITSIO
void nc_data_cluster_ncount_catalog_save (NcmData *data, gchar *filename, gboolean overwrite);
void nc_data_cluster_ncount_catalog_load (NcmData *data, gchar *filename);
#endif /* NUMCOSMO_HAVE_CFITSIO */

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_NCOUNT_H_ */

