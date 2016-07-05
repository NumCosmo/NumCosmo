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
#include <numcosmo/build_cfg.h>
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

struct _NcDataClusterNCount
{
  /*< private >*/
  NcmData parent_instance;
  NcClusterAbundance *cad;
  NcmVector *lnM_true;
  NcmVector *z_true;
  NcmMatrix *z_obs;
  NcmMatrix *z_obs_params;
  NcmMatrix *lnM_obs;
  NcmMatrix *lnM_obs_params;
  GArray *m2lnL_a;
  gdouble area_survey;
  guint np;
  guint n_z_obs;
  guint n_z_obs_params;
  guint n_M_obs;
  guint n_M_obs_params;
  gdouble log_np_fac;
  gboolean use_true_data;
  gboolean binned;
  NcmVector *z_nodes;
  NcmVector *lnM_nodes;
  gsl_histogram2d *purity;
  gsl_histogram2d *sd_lnM;
  gsl_histogram2d *z_lnM;
  gboolean fiducial;
  guint64 seed;
  gchar *rnd_name;
};

struct _NcDataClusterNCountClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

GType nc_data_cluster_ncount_get_type (void) G_GNUC_CONST;

NcDataClusterNCount *nc_data_cluster_ncount_new (NcClusterAbundance *cad);

NcDataClusterNCount *nc_data_cluster_ncount_ref (NcDataClusterNCount *ncount);
void nc_data_cluster_ncount_free (NcDataClusterNCount *ncount);
void nc_data_cluster_ncount_clear (NcDataClusterNCount **ncount);

void nc_data_cluster_ncount_set_n_z_obs (NcDataClusterNCount *ncount, guint n_z_obs);
void nc_data_cluster_ncount_set_n_z_obs_params (NcDataClusterNCount *ncount, guint n_z_obs_params);
void nc_data_cluster_ncount_set_n_M_obs (NcDataClusterNCount *ncount, guint n_M_obs);
void nc_data_cluster_ncount_set_n_M_obs_params (NcDataClusterNCount *ncount, guint n_M_obs_params);

void nc_data_cluster_ncount_set_lnM_true (NcDataClusterNCount *ncount, const NcmVector *v);
void nc_data_cluster_ncount_set_z_true (NcDataClusterNCount *ncount, const NcmVector *v);
void nc_data_cluster_ncount_set_lnM_obs (NcDataClusterNCount *ncount, const NcmMatrix *m);
void nc_data_cluster_ncount_set_lnM_obs_params (NcDataClusterNCount *ncount, const NcmMatrix *m);
void nc_data_cluster_ncount_set_z_obs (NcDataClusterNCount *ncount, const NcmMatrix *m);
void nc_data_cluster_ncount_set_z_obs_params (NcDataClusterNCount *ncount, const NcmMatrix *m);

gboolean nc_data_cluster_ncount_has_lnM_true (NcDataClusterNCount *ncount);
gboolean nc_data_cluster_ncount_has_z_true (NcDataClusterNCount *ncount);

guint nc_data_cluster_ncount_get_len (NcDataClusterNCount *ncount);
guint nc_data_cluster_ncount_lnM_obs_len (NcDataClusterNCount *ncount);
guint nc_data_cluster_ncount_lnM_obs_params_len (NcDataClusterNCount *ncount);

guint nc_data_cluster_ncount_z_obs_len (NcDataClusterNCount *ncount);
guint nc_data_cluster_ncount_z_obs_params_len (NcDataClusterNCount *ncount);

NcmVector *nc_data_cluster_ncount_get_lnM_true (NcDataClusterNCount *ncount);
NcmVector *nc_data_cluster_ncount_get_z_true (NcDataClusterNCount *ncount);

NcmMatrix *nc_data_cluster_ncount_get_lnM_obs (NcDataClusterNCount *ncount);
NcmMatrix *nc_data_cluster_ncount_get_lnM_obs_params (NcDataClusterNCount *ncount);

NcmMatrix *nc_data_cluster_ncount_get_z_obs (NcDataClusterNCount *ncount);
NcmMatrix *nc_data_cluster_ncount_get_z_obs_params (NcDataClusterNCount *ncount);

void nc_data_cluster_ncount_true_data (NcDataClusterNCount *ncount, gboolean use_true_data);
gboolean nc_data_cluster_ncount_using_true_data (NcDataClusterNCount *ncount);
void nc_data_cluster_ncount_init_from_sampling (NcDataClusterNCount *ncount, NcmMSet *mset, gdouble area_survey, NcmRNG *rng);
void nc_data_cluster_ncount_print (NcDataClusterNCount *ncount, NcHICosmo *cosmo, FILE *out, gchar *header);

void nc_data_cluster_ncount_set_bin_by_nodes (NcDataClusterNCount *ncount, NcmVector *z_nodes, NcmVector *lnM_nodes);
void nc_data_cluster_ncount_set_bin_by_minmax (NcDataClusterNCount *ncount, guint z_nbins, guint lnM_nbins);
void nc_data_cluster_ncount_set_bin_by_quantile (NcDataClusterNCount *ncount, NcmVector *z_quantiles, NcmVector *lnM_quantiles);
void nc_data_cluster_ncount_set_binned (NcDataClusterNCount *ncount, gboolean on);

#ifdef NUMCOSMO_HAVE_CFITSIO
void nc_data_cluster_ncount_catalog_save (NcDataClusterNCount *ncount, gchar *filename, gboolean overwrite);
void nc_data_cluster_ncount_catalog_load (NcDataClusterNCount *ncount, gchar *filename);
#endif /* NUMCOSMO_HAVE_CFITSIO */


G_END_DECLS

#endif /* _NC_DATA_CLUSTER_NCOUNT_H_ */

