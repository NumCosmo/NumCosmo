/***************************************************************************
 *            nc_data_cluster_pseudo_counts.c
 *
 *  Sun Apr 5 20:23:19 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_cluster_pseudo_counts.c
 * Copyright (C) 2015 Mariana Penna Lima <pennalima@gmail.com>
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

/**
 * SECTION:nc_data_cluster_pseudo_counts
 * @title: NcDataClusterPseudoCounts
 * @short_description: Galaxy clusters data -- pseudo number counts likelihood.
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_pseudo_counts.h"
#include "nc_hicosmo.h"
#include "lss/nc_cluster_redshift.h"
#include "lss/nc_cluster_mass.h"
#include "lss/nc_cluster_abundance.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_CAD,
  PROP_NP,
  PROP_OBS,
  PROP_TRUE_DATA,
  PROP_M_Z_FLAT_PRIOR, 
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataClusterPseudoCounts, nc_data_cluster_pseudo_counts, NCM_TYPE_DATA);

static void
nc_data_cluster_pseudo_counts_init (NcDataClusterPseudoCounts *dcpc)
{
  dcpc->cad           = NULL;
  dcpc->obs           = NULL;
  dcpc->true_data     = NULL;
  dcpc->np            = 0;
  dcpc->M_Z_FlatPrior = FALSE;
}

static void
nc_data_cluster_pseudo_counts_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (object);
  g_return_if_fail (NC_IS_DATA_CLUSTER_PSEUDO_COUNTS (object));

  switch (prop_id)
  {
    case PROP_CAD:
      nc_data_cluster_pseudo_counts_set_cad (dcpc, g_value_get_object (value));
      break;
    case PROP_NP:
      nc_data_cluster_pseudo_counts_set_nclusters (dcpc, g_value_get_uint (value));
      break;
    case PROP_OBS:
      nc_data_cluster_pseudo_counts_set_obs (dcpc, g_value_get_object (value));
      break;
    case PROP_TRUE_DATA:
      nc_data_cluster_pseudo_counts_set_true_data (dcpc, g_value_get_object (value));
      break;
    case PROP_M_Z_FLAT_PRIOR:
      dcpc->M_Z_FlatPrior = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_pseudo_counts_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (object);
  g_return_if_fail (NC_IS_DATA_CLUSTER_PSEUDO_COUNTS (object));

  switch (prop_id)
  {
    case PROP_CAD:
      g_value_set_object (value, dcpc->cad);
      break;
    case PROP_NP:
      g_value_set_uint (value, dcpc->np);
      break;
    case PROP_OBS:
      g_value_set_object (value, dcpc->obs);
      break;
    case PROP_TRUE_DATA:
      g_value_set_object (value, dcpc->true_data);
      break;
    case PROP_M_Z_FLAT_PRIOR:
      g_value_set_boolean (value, dcpc->M_Z_FlatPrior);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_pseudo_counts_dispose (GObject *object)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (object);

  nc_cluster_abundance_clear (&dcpc->cad);
  ncm_matrix_clear (&dcpc->obs);
  ncm_matrix_clear (&dcpc->true_data);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_pseudo_counts_parent_class)->dispose (object);
}

static void
nc_data_cluster_pseudo_counts_finalize (GObject *object)
{
  /*NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_pseudo_counts_parent_class)->finalize (object);
}

static void _nc_data_cluster_pseudo_counts_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static guint _nc_data_cluster_pseudo_counts_get_len (NcmData *data) { return NC_DATA_CLUSTER_PSEUDO_COUNTS (data)->np; }
static void _nc_data_cluster_pseudo_counts_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_pseudo_counts_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);

static void
nc_data_cluster_pseudo_counts_class_init (NcDataClusterPseudoCountsClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class = NCM_DATA_CLASS (klass);

  object_class->set_property = nc_data_cluster_pseudo_counts_set_property;
  object_class->get_property = nc_data_cluster_pseudo_counts_get_property;
  object_class->dispose      = nc_data_cluster_pseudo_counts_dispose;
  object_class->finalize     = nc_data_cluster_pseudo_counts_finalize;

  g_object_class_install_property (object_class,
                                   PROP_CAD,
                                   g_param_spec_object ("cluster-abundance",
                                                        NULL,
                                                        "Cluster abundance",
                                                        NC_TYPE_CLUSTER_ABUNDANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NP,
                                   g_param_spec_uint ("np",
                                                      NULL,
                                                      "Number of clusters",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                         NULL,
                                                         "Cluster observables",
                                                         NCM_TYPE_MATRIX,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRUE_DATA,
                                   g_param_spec_object ("true-data",
                                                         NULL,
                                                         "Cluster (halo) true data (redshift and mass)",
                                                         NCM_TYPE_MATRIX,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_M_Z_FLAT_PRIOR,
                                   g_param_spec_boolean ("M-z-flat-prior",
                                                         NULL,
                                                         "Flat priors for halo mass and selection functions.",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->m2lnL_val  = &_nc_data_cluster_pseudo_counts_m2lnL_val;
  data_class->get_length = &_nc_data_cluster_pseudo_counts_get_len;
  data_class->prepare    = &_nc_data_cluster_pseudo_counts_prepare;
  data_class->resample   = &_nc_data_cluster_pseudo_counts_resample;
}

static void
_nc_data_cluster_pseudo_counts_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (data);
  NcHICosmo *cosmo                = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterMass *clusterm         = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  NcClusterPseudoCounts *cpc      = NC_CLUSTER_PSEUDO_COUNTS (ncm_mset_peek (mset, nc_cluster_pseudo_counts_id ()));

  g_assert (cosmo != NULL && clusterm != NULL && cpc != NULL);

  if (dcpc->M_Z_FlatPrior)
  {
    /* ndet = 1, i.e., no-physical case: selection function and mass function equal to one (flat prior). */
    gint i;
    for (i = 0; i < dcpc->np; i++)
    {
      const gdouble z = ncm_matrix_get (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_Z);
      const gdouble *M = ncm_matrix_ptr (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_MPL);
      const gdouble *M_params = ncm_matrix_ptr (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MPL);
      const gdouble m2lnL_i = log (nc_cluster_pseudo_counts_posterior_ndetone (cpc, dcpc->cad->mfp, cosmo, clusterm, z, M[0], M[1], M_params[0], M_params[1]));

      if (!gsl_finite (m2lnL_i))
      {
        *m2lnL += m2lnL_i;
        break;
      }
      else
      {
        *m2lnL += m2lnL_i;
      }
    }
  }
  else
  {
    gdouble Ndet = nc_cluster_pseudo_counts_ndet (cpc, dcpc->cad->mfp, cosmo);
    gdouble lnNdet = log (Ndet);
    gint i;
    *m2lnL = 0.0;

    if (Ndet < 1.0)
    {
      *m2lnL = GSL_POSINF;
      return;
    }
    else
      *m2lnL = 0.0;

    for (i = 0; i < dcpc->np; i++)
    {
      const gdouble z = ncm_matrix_get (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_Z);
      const gdouble *M = ncm_matrix_ptr (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_MPL);
      const gdouble *M_params = ncm_matrix_ptr (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MPL);
      //const gdouble m2lnL_i = log (nc_cluster_pseudo_counts_posterior_numerator (cpc, clusterm, cosmo, z, M, M_params));
      const gdouble m2lnL_i = log (nc_cluster_pseudo_counts_posterior_numerator_plcl (cpc, dcpc->cad->mfp, clusterm, cosmo, z, M[0], M[1], M_params[0], M_params[1]));

      /*printf ("%d % 20.15g % 20.15g % 20.15g\n", i, z, log (nc_hicosmo_E (cosmo, z)), m2lnL_i);*/
      if (!gsl_finite (m2lnL_i))
      {
        *m2lnL += m2lnL_i;
        break;
      }
      else
      {
        *m2lnL += m2lnL_i;
      }
    }

    *m2lnL -= dcpc->np * lnNdet;
  }
  
  *m2lnL = -2.0 * (*m2lnL);
  return;
}

static void 
_nc_data_cluster_pseudo_counts_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (data);
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ())); 
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  NcClusterPseudoCounts *cpc  = NC_CLUSTER_PSEUDO_COUNTS (ncm_mset_peek (mset, nc_cluster_pseudo_counts_id ()));

  if (dcpc->cad == NULL)
    g_error ("_nc_data_cluster_pseudo_counts_prepare: NcClusterAbundance not set, call _l");
  
  g_assert ((cosmo != NULL) && (clusterz != NULL) && (clusterm != NULL) && (cpc != NULL));

  nc_cluster_abundance_prepare_if_needed (dcpc->cad, cosmo, clusterz, clusterm);  
}

static void
_nc_data_cluster_pseudo_counts_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (data);
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ())); 
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  NcClusterPseudoCounts *cpc  = NC_CLUSTER_PSEUDO_COUNTS (ncm_mset_peek (mset, nc_cluster_pseudo_counts_id ()));
  const gdouble lnMi          = nc_cluster_pseudo_counts_selection_function_lnMi (cpc, cosmo);
  guint tries = 0;
  guint i;

  g_assert (clusterz != NULL && clusterm != NULL && cosmo != NULL && cpc != NULL);
  g_assert (NC_IS_CLUSTER_MASS_PLCL (clusterm));
  g_assert_cmpuint (nc_cluster_redshift_obs_params_len (clusterz), ==, 0);

  nc_cluster_abundance_prepare_inv_dNdz (dcpc->cad, NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())),
                                         lnMi);

  for (i = 0; i < dcpc->np;)
  {
    ncm_rng_lock (rng);
    {
      const gdouble u1       = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng->r), dcpc->cad->z_epsilon);
      const gdouble u2       = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng->r), dcpc->cad->lnM_epsilon);
      const gdouble z_true   = ncm_spline_eval (dcpc->cad->inv_z, u1);
      const gdouble lnM_true = ncm_spline2d_eval (dcpc->cad->inv_lnM_z, u2, z_true);
      gdouble *zi_obs          = ncm_matrix_ptr (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_Z);
      gdouble *lnMi_obs        = ncm_matrix_ptr (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_MPL);
      gdouble *lnMi_obs_params = ncm_matrix_ptr (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MPL);
      gdouble val_sel = nc_cluster_pseudo_counts_selection_function (cpc, lnM_true, z_true);
      gdouble sel_u   = gsl_rng_uniform_pos (rng->r);

      ncm_rng_unlock (rng);

      if (tries++ > NC_DATA_CLUSTER_PSEUDO_COUNTS_RESAMPLE_MAX_TRIES)
        g_error ("_nc_data_cluster_pseudo_counts_resample: cannot resample, selection function cannot select mass in the current mass range.");

      if (sel_u >= val_sel)
        continue;

      if ( nc_cluster_redshift_resample (clusterz, lnM_true, z_true, zi_obs, NULL, rng) &&
          nc_cluster_mass_resample (clusterm, cosmo, lnM_true, z_true, lnMi_obs, lnMi_obs_params, rng) )
      {
        ncm_matrix_set (dcpc->true_data, i, 0, z_true);      
        ncm_matrix_set (dcpc->true_data, i, 1, lnM_true);
/*
        printf("# %d            [%.8g, %.8g, %.3g, %.8g, %.8g, %.8g, %.8g]\n", i, 
               z_true, lnM_true, 
               ncm_matrix_get (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_Z),
               ncm_matrix_get (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_MPL),
               ncm_matrix_get (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_MCL),
               ncm_matrix_get (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MPL),
               ncm_matrix_get (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MCL));
*/        
      }
      else
        continue;
      i++;
      tries = 0;
    }
  }

  {
    gchar *desc = g_strdup_printf ("Pseudo Counts Resample: Sample generated with n = %u clusters.", dcpc->np);
    ncm_data_take_desc (data, desc);
  }
}

/**
 * nc_data_cluster_pseudo_counts_new:
 * @cad: a #NcClusterAbundance
 * 
 * Creates a new #NcDataClusterPseudoCounts.
 * 
 * Returns: the newly created #NcDataClusterPseudoCounts.
 */
NcDataClusterPseudoCounts *
nc_data_cluster_pseudo_counts_new (NcClusterAbundance *cad)
{
  NcDataClusterPseudoCounts *dcpc;

  dcpc = g_object_new (NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS,
                       "cluster-abundance", cad,
                       NULL);
  
  return dcpc;
}

/**
 * nc_data_cluster_pseudo_counts_new_from_file:
 * @filename: file containing a serialized #NcDataClusterPseudoCounts
 * 
 * Creates a new #NcDataClusterPseudoCounts from @filename.
 * 
 * Returns: (transfer full): the newly created #NcDataClusterPseudoCounts.
 */
NcDataClusterPseudoCounts *
nc_data_cluster_pseudo_counts_new_from_file (const gchar *filename)
{
  NcDataClusterPseudoCounts *dcpc = NC_DATA_CLUSTER_PSEUDO_COUNTS (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_CLUSTER_PSEUDO_COUNTS (dcpc));

  return dcpc;
}

/**
 * nc_data_cluster_pseudo_counts_ref:
 * @dcpc: a #NcDataClusterPseudoCounts
 *
 * Increases the reference count of @dcpc by one.
 * 
 * Returns: (transfer full): @dcpc
 */
NcDataClusterPseudoCounts *
nc_data_cluster_pseudo_counts_ref (NcDataClusterPseudoCounts *dcpc)
{
  return g_object_ref (dcpc);
}

/**
 * nc_data_cluster_pseudo_counts_free:
 * @dcpc: a #NcDataClusterPseudoCounts
 *
 * Atomically decrements the reference count of @dcpc by one. If the reference count drops to 0,
 * all memory allocated by @dcpc is released.
 * 
 */
void
nc_data_cluster_pseudo_counts_free (NcDataClusterPseudoCounts *dcpc)
{
  g_object_unref (dcpc);
}

/**
 * nc_data_cluster_pseudo_counts_clear:
 * @dcpc: a #NcDataClusterPseudoCounts
 *
 * The reference count of @dcpc is decreased and the pointer is set to NULL.
 * 
 */
void
nc_data_cluster_pseudo_counts_clear (NcDataClusterPseudoCounts **dcpc)
{
  g_clear_object (dcpc);
}

/**
 * nc_data_cluster_pseudo_counts_set_cad:
 * @dcpc: a #NcDataClusterPseudoCounts
 * @cad: a #NcClusterAbundance
 *
 * Sets the #NcClusterAbundance object @cad to be used.
 * 
 */
void 
nc_data_cluster_pseudo_counts_set_cad (NcDataClusterPseudoCounts *dcpc, NcClusterAbundance *cad)
{
  nc_cluster_abundance_clear (&dcpc->cad);
  if (cad != NULL)
  {
    dcpc->cad = nc_cluster_abundance_ref (cad);
  }
}

/**
 * nc_data_cluster_pseudo_counts_set_nclusters:
 * @dcpc: a #NcDataClusterPseudoCounts
 * @np: number of clusters
 *
 * Sets @np representing the total number of clusters.
 * 
 */
void 
nc_data_cluster_pseudo_counts_set_nclusters (NcDataClusterPseudoCounts *dcpc, guint np)
{
  if (np == dcpc->np)
    return;
  else
  {
    ncm_matrix_clear (&dcpc->obs);
    ncm_matrix_clear (&dcpc->true_data);
    dcpc->np        = np;
    if (np > 0)
    {      
      dcpc->obs       = ncm_matrix_new (dcpc->np, NC_DATA_CLUSTER_PSEUDO_COUNTS_LEN);
      dcpc->true_data = ncm_matrix_new (dcpc->np, 2);
    }
    ncm_data_set_init (NCM_DATA (dcpc), FALSE);
  }  
}

/**
 * nc_data_cluster_pseudo_counts_get_nclusters:
 * @dcpc: a #NcDataClusterPseudoCounts
 * 
 * Returns: the number of clusters 
 */
guint 
nc_data_cluster_pseudo_counts_get_obs (NcDataClusterPseudoCounts *dcpc)
{
  return dcpc->np;
}

/**
 * nc_data_cluster_pseudo_counts_set_obs:
 * @dcpc: a #NcDataClusterPseudoCounts
 * @m: a #NcmMatrix
 *
 * Sets the matrix @m representing the cluster observables, e.g., observed redshift and 
 * mass(es) and the parameters of the redshift and/or mass-observable distributions. 
 * 
 * The function nc_data_cluster_pseudo_counts_set_nclusters must be called before this one.
 *
 */
void 
nc_data_cluster_pseudo_counts_set_obs (NcDataClusterPseudoCounts *dcpc, const NcmMatrix *m)
{
  g_assert (m != NULL);
  g_assert_cmpuint (ncm_matrix_nrows (m), ==, dcpc->np);
  g_assert_cmpuint (ncm_matrix_ncols (m), ==, NC_DATA_CLUSTER_PSEUDO_COUNTS_LEN);
  
  ncm_matrix_memcpy (dcpc->obs, m);
  ncm_data_set_init (NCM_DATA (dcpc), TRUE);
}

/**
 * nc_data_cluster_pseudo_counts_set_true_data:
 * @dcpc: a #NcDataClusterPseudoCounts
 * @m: a #NcmMatrix
 *
 * Sets the matrix @m representing the cluster (halo) true data, i.e., true redshift and 
 * mass(es).
 * 
 */
void 
nc_data_cluster_pseudo_counts_set_true_data (NcDataClusterPseudoCounts *dcpc, const NcmMatrix *m)
{
  g_assert (m != NULL);
  g_assert_cmpuint (ncm_matrix_nrows (m), ==, dcpc->np);
  g_assert_cmpuint (ncm_matrix_ncols (m), ==, 2);
       
  ncm_matrix_memcpy (dcpc->true_data, m);
}


/**
 * nc_data_cluster_pseudo_counts_init_from_sampling:
 * @dcpc: a #NcDataClusterPseudoCounts
 * @mset: a #NcmMSet
 * @rng: a #NcmRNG
 * @np: number of clusters
 * 
 * FIXME
 *
 */
void
nc_data_cluster_pseudo_counts_init_from_sampling (NcDataClusterPseudoCounts *dcpc, NcmMSet *mset, NcmRNG *rng, const gint np)
{
  const gdouble M0 = 1.0e14;
  gint i;
  
  nc_data_cluster_pseudo_counts_set_nclusters (dcpc, np);

  ncm_rng_lock (rng);
  
  for (i = 0; i < np; i++)
  {
    const gdouble sd_PL = gsl_ran_gamma (rng->r, 10.0, 0.5) * M0 / 10.0;      /* This choice is suitable for the Planck catalog, logarithm (e) scale*/
    const gdouble sd_CL = (1.8 + gsl_ran_rayleigh (rng->r, 3.5)) * M0 / 10.0; /* This choice is suitable for the CLASH catalog, logarithm (e) scale*/
    
    ncm_matrix_set (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MPL, sd_PL);
    ncm_matrix_set (dcpc->obs, i, NC_DATA_CLUSTER_PSEUDO_COUNTS_SD_MCL, sd_CL);
  }
  
  ncm_rng_unlock (rng);  

  ncm_data_set_init (NCM_DATA (dcpc), TRUE);
  
  ncm_data_resample (NCM_DATA (dcpc), mset, rng);
}
