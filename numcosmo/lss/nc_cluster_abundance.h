/***************************************************************************
 *            nc_cluster_abundance.h
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <sandro@isoftware.com.br>
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

#ifndef _NC_CLUSTER_ABUNDANCE_H_
#define _NC_CLUSTER_ABUNDANCE_H_

#include <glib-object.h>
#include <glib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_ABUNDANCE             (nc_cluster_abundance_get_type ())
#define NC_CLUSTER_ABUNDANCE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_ABUNDANCE, NcClusterAbundance))
#define NC_CLUSTER_ABUNDANCE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_ABUNDANCE, NcClusterAbundanceClass))
#define NC_IS_CLUSTER_ABUNDANCE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_ABUNDANCE))
#define NC_IS_CLUSTER_ABUNDANCE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_ABUNDANCE))
#define NC_CLUSTER_ABUNDANCE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_ABUNDANCE, NcClusterAbundanceClass))

typedef struct _NcClusterAbundanceClass NcClusterAbundanceClass;
typedef struct _NcClusterAbundance NcClusterAbundance;
typedef struct _NcClusterAbundanceDataP NcClusterAbundanceDataP;
typedef struct _NcClusterAbundanceDataBinz NcClusterAbundanceDataBinZ;
typedef struct _NcClusterAbundanceDataBinM NcClusterAbundanceDataBinM;
typedef struct _NcClusterAbundanceDataBin NcClusterAbundanceDataBin;

typedef gdouble (*NcClusterAbundanceN) (NcClusterAbundance *cad, NcHICosmo *model);
typedef gdouble (*NcClusterAbundanceIntPd2N) (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z);

#define nc_cluster_abundance_d2NdzdlnM_val(cad,cp,lnM,z) (cad)->d2NdzdlnM_val(cad,cp,lnM,z)
#define nc_cluster_abundance_dNdz_val(cad,cp,lnMl,lnMu,z) (cad)->dNdz_val(cad,cp,lnMl,lnMu,z)
#define nc_cluster_abundance_dNdlnM_val(cad,cp,lnM,zl,zu) (cad)->dNdlnM_val(cad,cp,lnM,zl,zu)
#define nc_cluster_abundance_N_val(cad,cp,lnMl,lnMu,zl,zu) (cad)->N_val(cad,cp,lnMl,lnMu,zl,zu)

struct _NcClusterAbundanceClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcClusterAbundance
{
  /*< private >*/
  GObject parent_instance;
  NcMassFunction *mfp;
  NcHaloBiasFunc *mbiasf;    /* new FIXME */
  NcClusterRedshift *z;
  NcClusterMass *m;
  NcClusterAbundanceN N;
  NcClusterAbundanceIntPd2N intp_d2N;
  gdouble norma, log_norma;
  gdouble lnMi, lnMf, zi, zf;
  gdouble completeness_factor;
  gdouble lnM_epsilon, z_epsilon;
  gboolean optimize;
  gsl_histogram2d *completeness;
  gsl_histogram2d *purity;
  gsl_histogram2d *sd_lnM;
  NcmSpline2d *dbdlnM;    /* To compute the mean bias. FIXME*/
  NcmSpline *inv_z;
  NcmSpline *inv_lnM;
  NcmSpline2d *inv_lnM_z;
  gsl_rng *rng;
  NcmModelCtrl *ctrl;
};

GType nc_cluster_abundance_get_type (void) G_GNUC_CONST;

NcClusterAbundance *nc_cluster_abundance_new (NcMassFunction *mfp, NcHaloBiasFunc *mbiasf, NcClusterRedshift *clusterz, NcClusterMass *clusterm);
NcClusterAbundance *nc_cluster_abundance_copy (NcClusterAbundance *cad);
NcClusterAbundance *nc_cluster_abundance_ref (NcClusterAbundance *cad);
void nc_cluster_abundance_free (NcClusterAbundance *cad);

void nc_cluster_abundance_prepare (NcClusterAbundance *cad, NcHICosmo *model);
void nc_cluster_abundance_prepare_inv_dNdz (NcClusterAbundance *cad, NcHICosmo *model);
void nc_cluster_abundance_prepare_inv_dNdlnM_z (NcClusterAbundance *cad, NcHICosmo *model, gdouble z);

gdouble nc_cluster_abundance_z_p_lnm_p_d2n (NcClusterAbundance *cad, NcHICosmo *model, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *z_obs, gdouble *z_obs_params);
gdouble nc_cluster_abundance_z_p_d2n (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble *z_obs, gdouble *z_obs_params);
gdouble nc_cluster_abundance_lnm_p_d2n (NcClusterAbundance *cad, NcHICosmo *model, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble z);
gdouble nc_cluster_abundance_d2n (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z);

gdouble nc_cluster_abundance_true_n (NcClusterAbundance *cad, NcHICosmo *model);
gdouble nc_cluster_abundance_n (NcClusterAbundance *cad, NcHICosmo *model);
gdouble nc_cluster_abundance_intp_d2n (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z);

void nc_cluster_abundance_set_redshift (NcClusterAbundance *cad, NcClusterRedshift *clusterz);
void nc_cluster_abundance_set_mass (NcClusterAbundance *cad, NcClusterMass *clusterm);
NcClusterRedshift *nc_cluster_abundance_get_redshift (NcClusterAbundance *cad);
NcClusterMass *nc_cluster_abundance_get_mass (NcClusterAbundance *cad);
G_INLINE_FUNC NcClusterRedshift *nc_cluster_abundance_peek_redshift (NcClusterAbundance *cad);
G_INLINE_FUNC NcClusterMass *nc_cluster_abundance_peek_mass (NcClusterAbundance *cad);

/*
void nc_cluster_abundance_bin_realization (GArray *zr, gsl_histogram **h);
void nc_cluster_abundance_realizations_save_to_file (GPtrArray *realizations, gchar *filename);
GPtrArray *nc_cluster_abundance_realizations_read_from_file (gchar *file_realization, gint n_realizations);
gdouble nc_cluster_abundance_d2NdzdlnM_photoz (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble *z_obs, gdouble *z_obs_params);
gdouble nc_cluster_abundance_d2NdzdlnM_Mobs (NcClusterAbundance *cad, NcHICosmo *model, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble z);
gdouble nc_cluster_abundance_d2NdzdlnM_photoz_Mobs (NcClusterAbundance *cad, NcHICosmo *model, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *z_obs, gdouble *z_obs_params);
gdouble nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble z);
gdouble nc_cluster_abundance_dNdz_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs_i, gdouble lnMobs_f, gdouble z);
gdouble nc_cluster_abundance_N_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs_i, gdouble lnMobs_f, gdouble z_i, gdouble z_f);
gdouble nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble zp);
gdouble nc_cluster_abundance_dNdz_photoz_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs_i, gdouble lnMobs_f, gdouble z);
gdouble nc_cluster_abundance_N_photoz_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs_i, gdouble lnMobs_f, gdouble zp_i, gdouble zp_f);
*/

void nc_bias_mean_prepare (NcClusterAbundance *cad, NcHICosmo *model);
gdouble nc_bias_mean_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z);
gdouble nc_ca_mean_bias_numerator (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z);
gdouble nc_ca_mean_bias_denominator (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z);
gdouble nc_ca_mean_bias (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z);
gdouble nc_ca_mean_bias_Mobs_numerator (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble z);
gdouble nc_ca_mean_bias_Mobs_denominator (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble z);

#define _NC_CLUSTER_ABUNDANCE_NNODES 1000
#define _NC_CLUSTER_ABUNDANCE_MIN_Z  0.0

gdouble _nc_cad_inv_dNdz_convergence_f (gdouble n, gdouble epsilon);

G_END_DECLS

#endif /* _NC_CLUSTER_ABUNDANCE_H_ */

#ifndef _NC_CLUSTER_ABUNDANCE_INLINE_H_
#define _NC_CLUSTER_ABUNDANCE_INLINE_H_

#include <glib-object.h>
#include <glib.h>

G_BEGIN_DECLS

G_INLINE_FUNC NcClusterRedshift *
nc_cluster_abundance_peek_redshift (NcClusterAbundance *cad)
{
  return cad->z;
}

G_INLINE_FUNC NcClusterMass *
nc_cluster_abundance_peek_mass (NcClusterAbundance *cad)
{
  return cad->m;
}

G_END_DECLS

#endif /* _NC_CLUSTER_ABUNDANCE_INLINE_H_ */


