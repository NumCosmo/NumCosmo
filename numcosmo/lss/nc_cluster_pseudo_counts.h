/***************************************************************************
 *            nc_cluster_pseudo_counts.h
 *
 *  Mon Mar 30 02:05:16 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2015 <pennalima@gmail.com> 
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

#ifndef _NC_CLUSTER_PSEUDO_COUNTS_H_
#define _NC_CLUSTER_PSEUDO_COUNTS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/lss/nc_halo_mass_function.h>
#include <numcosmo/lss/nc_cluster_abundance.h>
#include <numcosmo/lss/nc_cluster_mass.h>
#include <numcosmo/lss/nc_cluster_mass_plcl.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multifit_nlin.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_PSEUDO_COUNTS            (nc_cluster_pseudo_counts_get_type ())
#define NC_CLUSTER_PSEUDO_COUNTS(obj)            (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_PSEUDO_COUNTS, NcClusterPseudoCounts))
#define NC_CLUSTER_PSEUDO_COUNTS_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_PSEUDO_COUNTS, NcClusterPseudoCountsClass))
#define NC_IS_CLUSTER_PSEUDO_COUNTS(obj)         (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_PSEUDO_COUNTS))
#define NC_IS_CLUSTER_PSEUDO_COUNTS_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_PSEUDO_COUNTS))
#define NC_CLUSTER_PSEUDO_COUNTS_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_PSEUDO_COUNTS, NcClusterPseudoCountsClass))

typedef struct _NcClusterPseudoCountsClass NcClusterPseudoCountsClass;
typedef struct _NcClusterPseudoCounts NcClusterPseudoCounts;

/**
 * NcClusterPseudoCountsSParams:
 * @NC_CLUSTER_PSEUDO_COUNTS_LNMCUT: logarithm base e of the mass cut-off
 * @NC_CLUSTER_PSEUDO_COUNTS_SD_MCUT: standard deviation of mass cut-off
 * @NC_CLUSTER_PSEUDO_COUNTS_ZMIN: minimum redshift
 * @NC_CLUSTER_PSEUDO_COUNTS_DELTAZ: redshift interval size
 * 
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_PSEUDO_COUNTS_SPARAMS >*/
{
  NC_CLUSTER_PSEUDO_COUNTS_LNMCUT = 0, 
  NC_CLUSTER_PSEUDO_COUNTS_SD_MCUT, 
  NC_CLUSTER_PSEUDO_COUNTS_ZMIN, 
  NC_CLUSTER_PSEUDO_COUNTS_DELTAZ, 
  /* < private > */
  NC_CLUSTER_PSEUDO_COUNTS_SPARAM_LEN, /*< skip >*/
} NcClusterPseudoCountsSParams;

#define NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_LNMCUT  (33.0)
#define NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_SD_MCUT (0.206)
#define NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_ZMIN (0.188)
#define NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_DELTAZ (0.99)

#define NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_PARAMS_ABSTOL (0.0)

struct _NcClusterPseudoCountsClass
{
  /*< private >*/
  NcmModelClass parent_class;
};

struct _NcClusterPseudoCounts
{
  /*< private >*/
  NcmModel parent_instance;
  guint nclusters;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  gdouble *workz;
};

GType nc_cluster_pseudo_counts_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_cluster_pseudo_counts);

NcClusterPseudoCounts *nc_cluster_pseudo_counts_new (guint nclusters);
NcClusterPseudoCounts *nc_cluster_pseudo_counts_ref (NcClusterPseudoCounts *cpc);
void nc_cluster_pseudo_counts_free (NcClusterPseudoCounts *cpc);
void nc_cluster_pseudo_counts_clear (NcClusterPseudoCounts **cpc);

gdouble nc_cluster_pseudo_counts_posterior_ndetone (NcClusterPseudoCounts *cpc, NcHaloMassFunction *mfp, NcHICosmo *cosmo, NcClusterMass *clusterm, gdouble z, gdouble Mpl, gdouble Mcl, gdouble sigma_pl, gdouble sigma_cl);
gdouble nc_cluster_pseudo_counts_selection_function (NcClusterPseudoCounts *cpc, gdouble lnM, gdouble z);
gdouble nc_cluster_pseudo_counts_selection_function_lnMi (NcClusterPseudoCounts *cpc, NcHICosmo *cosmo);
gdouble nc_cluster_pseudo_counts_ndet_no_z_integral (NcClusterPseudoCounts *cpc, NcHICosmo *cosmo, gdouble z);
gdouble nc_cluster_pseudo_counts_ndet (NcClusterPseudoCounts *cpc, NcHaloMassFunction *mfp, NcHICosmo *cosmo);
gdouble nc_cluster_pseudo_counts_posterior_numerator (NcClusterPseudoCounts *cpc, NcHaloMassFunction *mfp, NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble z, const gdouble *Mobs, const gdouble *Mobs_params);
gdouble nc_cluster_pseudo_counts_mf_lognormal_integral (NcClusterPseudoCounts *cpc, NcHaloMassFunction *mfp, NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnMsz, const gdouble lnMl, const gdouble z);
gdouble nc_cluster_pseudo_counts_posterior_numerator_plcl (NcClusterPseudoCounts *cpc, NcHaloMassFunction *mfp, NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble z, const gdouble Mpl, const gdouble Mcl, const gdouble sigma_pl, const gdouble sigma_cl);

G_END_DECLS

#endif /* _NC_CLUSTER_PSEUDO_COUNTS_H_ */
