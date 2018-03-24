/***************************************************************************
 *            nc_data_reduced_shear_cluster_mass.c
 *
 *  Thu Mar 22 16:10:25 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_reduced_shear_cluster_mass.c
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

/**
 * SECTION:nc_data_reduced_shear_cluster_mass
 * @title: NcDataReducedShearClusterMass
 * @short_description: Galaxy clusters data -- pseudo number counts likelihood.
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_reduced_shear_cluster_mass.h"
#include "nc_hicosmo.h"
#include "lss/nc_cluster_redshift.h"
#include "lss/nc_cluster_mass.h"
#include "lss/nc_cluster_abundance.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_NGALS,
  PROP_NZBINS,
  PROP_TRUE_DATA, 
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataReducedShearClusterMass, nc_data_reduced_shear_cluster_mass, NCM_TYPE_DATA);

static void
nc_data_reduced_shear_cluster_mass_init (NcDataReducedShearClusterMass *drs)
{
  drs->rscm      = NULL;
  drs->g_obs      = NULL;
  drs->Pz        = NULL;
  drs->true_data = NULL;
  drs->ngals     = 0;
  drs->nzbins    = 0;
}

static void
nc_data_reduced_shear_cluster_mass_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (object);
  g_return_if_fail (NC_IS_DATA_REDUCED_SHEAR_CLUSTER_MASS (object));

  switch (prop_id)
  {
    case PROP_NGALS:
      nc_data_reduced_shear_cluster_mass_set_ngalaxies (drs, g_value_get_uint (value));
      break;
    case PROP_NZBINS:
      nc_data_reduced_shear_cluster_mass_set_nzbins (drs, g_value_get_uint (value));
      break;
    case PROP_TRUE_DATA:
      nc_data_reduced_shear_cluster_mass_set_true_data (drs, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_reduced_shear_cluster_mass_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (object);
  g_return_if_fail (NC_IS_DATA_REDUCED_SHEAR_CLUSTER_MASS (object));

  switch (prop_id)
  {
    case PROP_NGALS:
      g_value_set_uint (value, drs->ngals);
      break;
    case PROP_NZBINS:
      g_value_set_uint (value, drs->nzbins);
      break;
    case PROP_TRUE_DATA:
      g_value_set_object (value, drs->true_data);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_reduced_shear_cluster_mass_dispose (GObject *object)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (object);

  ncm_matrix_clear (&drs->g_obs);
  ncm_matrix_clear (&drs->Pz);
  ncm_matrix_clear (&drs->true_data);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_reduced_shear_cluster_mass_parent_class)->dispose (object);
}

static void
nc_data_reduced_shear_cluster_mass_finalize (GObject *object)
{
  /*NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_reduced_shear_cluster_mass_parent_class)->finalize (object);
}

static void _nc_data_reduced_shear_cluster_mass_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static guint _nc_data_reduced_shear_cluster_mass_get_len (NcmData *data) { return NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (data)->ngals; }
//static void _nc_data_reduced_shear_cluster_mass_prepare (NcmData *data, NcmMSet *mset);

static void
nc_data_reduced_shear_cluster_mass_class_init (NcDataReducedShearClusterMassClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class = NCM_DATA_CLASS (klass);

  object_class->set_property = nc_data_reduced_shear_cluster_mass_set_property;
  object_class->get_property = nc_data_reduced_shear_cluster_mass_get_property;
  object_class->dispose      = nc_data_reduced_shear_cluster_mass_dispose;
  object_class->finalize     = nc_data_reduced_shear_cluster_mass_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NGALS,
                                   g_param_spec_uint ("ngals",
                                                      NULL,
                                                      "Number of galaxies",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NZBINS,
                                   g_param_spec_object ("nzbins",
                                                         NULL,
                                                         "Number of redshift bins",
                                                         NCM_TYPE_MATRIX,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRUE_DATA,
                                   g_param_spec_object ("true-data",
                                                         NULL,
                                                         "Cluster (halo) true data (redshift and mass)",
                                                         NCM_TYPE_MATRIX,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
 
  data_class->m2lnL_val  = &_nc_data_reduced_shear_cluster_mass_m2lnL_val;
  data_class->get_length = &_nc_data_reduced_shear_cluster_mass_get_len;
  //data_class->prepare    = &_nc_data_reduced_shear_cluster_mass_prepare;
}


static void
_nc_data_reduced_shear_cluster_mass_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (data);
  NcHICosmo *cosmo                = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  
  g_assert (cosmo != NULL);

  {
    gdouble Ndet = 0.0; //nc_cluster_pseudo_counts_ndet (cpc, drs->cad->mfp, cosmo);
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

    for (i = 0; i < drs->ngals; i++)
    {
      const gdouble z = 0.0; //ncm_matrix_get (drs->gobs, i, NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_Z);
      const gdouble *M = NULL; //ncm_matrix_ptr (drs->gobs, i, NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_MPL);
      const gdouble *M_params = NULL; //ncm_matrix_ptr (drs->gobs, i, NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_SD_MPL);
      //const gdouble m2lnL_i = log (nc_cluster_pseudo_counts_posterior_numerator (cpc, clusterm, cosmo, z, M, M_params));
      const gdouble m2lnL_i = 0.0; //log (nc_cluster_pseudo_counts_posterior_numerator_plcl (cpc, drs->cad->mfp, clusterm, cosmo, z, M[0], M[1], M_params[0], M_params[1]));

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

    *m2lnL -= drs->ngals * lnNdet;
  }
  
  *m2lnL = -2.0 * (*m2lnL);
  return;
}

/*
static void 
_nc_data_reduced_shear_cluster_mass_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (data);
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ())); 
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  NcClusterPseudoCounts *cpc  = NC_CLUSTER_PSEUDO_COUNTS (ncm_mset_peek (mset, nc_cluster_pseudo_counts_id ()));

  if (drs->cad == NULL)
    g_error ("_nc_data_reduced_shear_cluster_mass_prepare: NcClusterAbundance not set, call _l");
  
  g_assert ((cosmo != NULL) && (clusterz != NULL) && (clusterm != NULL) && (cpc != NULL));

  nc_cluster_abundance_prepare_if_needed (drs->cad, cosmo, clusterz, clusterm);  
}
*/

/**
 * nc_data_reduced_shear_cluster_mass_new:
 * 
 * Creates a new #NcDataReducedShearClusterMass.
 * 
 * Returns: the newly created #NcDataReducedShearClusterMass.
 */
NcDataReducedShearClusterMass *
nc_data_reduced_shear_cluster_mass_new ()
{
  NcDataReducedShearClusterMass *drs;

  drs = g_object_new (NC_TYPE_DATA_REDUCED_SHEAR_CLUSTER_MASS,
                       NULL);
  
  return drs;
}

/**
 * nc_data_reduced_shear_cluster_mass_new_from_file:
 * @filename: file containing a serialized #NcDataReducedShearClusterMass
 * 
 * Creates a new #NcDataReducedShearClusterMass from @filename.
 * 
 * Returns: (transfer full): the newly created #NcDataReducedShearClusterMass.
 */
NcDataReducedShearClusterMass *
nc_data_reduced_shear_cluster_mass_new_from_file (const gchar *filename)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_REDUCED_SHEAR_CLUSTER_MASS (drs));

  return drs;
}

/**
 * nc_data_reduced_shear_cluster_mass_ref:
 * @drs: a #NcDataReducedShearClusterMass
 *
 * Increases the reference count of @drs by one.
 * 
 * Returns: (transfer full): @drs
 */
NcDataReducedShearClusterMass *
nc_data_reduced_shear_cluster_mass_ref (NcDataReducedShearClusterMass *drs)
{
  return g_object_ref (drs);
}

/**
 * nc_data_reduced_shear_cluster_mass_free:
 * @drs: a #NcDataReducedShearClusterMass
 *
 * Atomically decrements the reference count of @drs by one. If the reference count drops to 0,
 * all memory allocated by @drs is released.
 * 
 */
void
nc_data_reduced_shear_cluster_mass_free (NcDataReducedShearClusterMass *drs)
{
  g_object_unref (drs);
}

/**
 * nc_data_reduced_shear_cluster_mass_clear:
 * @drs: a #NcDataReducedShearClusterMass
 *
 * The reference count of @drs is decreased and the pointer is set to NULL.
 * 
 */
void
nc_data_reduced_shear_cluster_mass_clear (NcDataReducedShearClusterMass **drs)
{
  g_clear_object (drs);
}

/**
 * nc_data_reduced_shear_cluster_mass_set_ngalaxies:
 * @drs: a #NcDataReducedShearClusterMass
 * @ngals: number of galaxies
 *
 * Sets @ngals representing the total number of galaxies that belongs to a cluster.
 * 
 */
void 
nc_data_reduced_shear_cluster_mass_set_ngalaxies (NcDataReducedShearClusterMass *drs, guint ngals)
{
  if (ngals == drs->ngals)
    return;
  else
  {
    ncm_matrix_clear (&drs->g_obs);
    ncm_matrix_clear (&drs->true_data);
    drs->ngals        = ngals;
    if (ngals > 0)
    {      
      drs->g_obs       = ncm_matrix_new (drs->ngals, NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_LEN);
      drs->true_data = ncm_matrix_new (drs->ngals, 2);
    }
    ncm_data_set_init (NCM_DATA (drs), FALSE);
  }  
}

/**
 * nc_data_reduced_shear_cluster_mass_get_ngalaxies:
 * @drs: a #NcDataReducedShearClusterMass
 * 
 * Returns: the number of galaxies 
 */
guint 
nc_data_reduced_shear_cluster_mass_get_ngalaxies (NcDataReducedShearClusterMass *drs)
{
  return drs->ngals;
}

/**
 * nc_data_reduced_shear_cluster_mass_set_nzbins:
 * @drs: a #NcDataReducedShearClusterMass
 * @nzbins: number of redshift bins
 *
 * Sets @nzbins representing the number of redshift bins of the (photometric) redshift distribution.
 * 
 */
void 
nc_data_reduced_shear_cluster_mass_set_nzbins (NcDataReducedShearClusterMass *drs, guint nzbins)
{
  if (nzbins == drs->nzbins)
    return;
  else
  {
    ncm_matrix_clear (&drs->g_obs);
    ncm_matrix_clear (&drs->true_data);
    drs->nzbins        = nzbins;
    if (nzbins > 0)
    {      
      drs->Pz       = ncm_matrix_new (drs->nzbins, NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_LEN);
      drs->true_data = ncm_matrix_new (drs->nzbins, 2);
    }
    ncm_data_set_init (NCM_DATA (drs), FALSE);
  }  
}

/**
 * nc_data_reduced_shear_cluster_mass_get_nzbins:
 * @drs: a #NcDataReducedShearClusterMass
 * 
 * Returns: the number of redshift bins
 */
guint 
nc_data_reduced_shear_cluster_mass_get_nzbins (NcDataReducedShearClusterMass *drs)
{
  return drs->nzbins;
}

/**
 * nc_data_reduced_shear_cluster_mass_set_gobs:
 * @drs: a #NcDataReducedShearClusterMass
 * @m: a #NcmMatrix
 *
 * Sets the matrix @m representing the measured reduced shear of the galaxies that belong to a cluster.
 * 
 * The function nc_data_reduced_shear_cluster_mass_set_ngalaxies must be called before this one.
 *
 */
void 
nc_data_reduced_shear_cluster_mass_set_gobs (NcDataReducedShearClusterMass *drs, const NcmMatrix *m)
{
  g_assert (m != NULL);
  g_assert_cmpuint (ncm_matrix_nrows (m), ==, drs->ngals);
  g_assert_cmpuint (ncm_matrix_ncols (m), ==, NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_LEN);
  
  ncm_matrix_memcpy (drs->g_obs, m);
  ncm_data_set_init (NCM_DATA (drs), TRUE);
}

/**
 * nc_data_reduced_shear_cluster_mass_set_pz:
 * @drs: a #NcDataReducedShearClusterMass
 * @m: a #NcmMatrix
 *
 * Sets the matrix @m representing the redshift probability distribution of a galaxy.
 * 
 * The function nc_data_reduced_shear_cluster_mass_set_nzbins must be called before this one.
 *
 */
void 
nc_data_reduced_shear_cluster_mass_set_pz (NcDataReducedShearClusterMass *drs, const NcmMatrix *m)
{
  g_assert (m != NULL);
  g_assert_cmpuint (ncm_matrix_nrows (m), ==, drs->nzbins);
  g_assert_cmpuint (ncm_matrix_ncols (m), ==, NC_DATA_REDUCED_SHEAR_CLUSTER_MASS_LEN);
  
  ncm_matrix_memcpy (drs->Pz, m);
  ncm_data_set_init (NCM_DATA (drs), TRUE);
}

/**
 * nc_data_reduced_shear_cluster_mass_set_true_data:
 * @drs: a #NcDataReducedShearClusterMass
 * @m: a #NcmMatrix
 *
 * FIXME
 * Sets the matrix @m representing 
 * 
 */
void 
nc_data_reduced_shear_cluster_mass_set_true_data (NcDataReducedShearClusterMass *drs, const NcmMatrix *m)
{
  g_assert (m != NULL);
  g_assert_cmpuint (ncm_matrix_nrows (m), ==, drs->ngals);
  g_assert_cmpuint (ncm_matrix_ncols (m), ==, 2);
       
  ncm_matrix_memcpy (drs->true_data, m);
}

