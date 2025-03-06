/***************************************************************************
 *            nc_data_cluster_mass_rich.c
 *
 *  Tue Oct 24 22:22:00 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_cluster_mass_rich.c
 * Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcDataClusterMassRich:
 *
 * Cluster mass richness data object.
 *
 * #NcDataClusterMassRich is a data object that contains the information about a cluster
 * mass richness data. Its main purpose is to be used as a data object for fitting mass
 * richness proxy models.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_mass_rich.h"
#include "lss/nc_cluster_mass_ascaso.h"
#include "lss/nc_cluster_mass_lnrich_ext.h"

#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"

#include <gsl/gsl_sf_erf.h>

typedef struct _NcDataClusterMassRichPrivate
{
  NcmVector *z_cluster;
  NcmVector *lnM_cluster;
  NcmVector *lnR_cluster;
} NcDataClusterMassRichPrivate;

enum
{
  PROP_0,
  PROP_Z_CLUSTER,
  PROP_LNM_CLUSTER,
  PROP_LNR_CLUSTER,
  PROP_SIZE,
};

struct _NcDataClusterMassRich
{
  NcmData parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterMassRich, nc_data_cluster_mass_rich, NCM_TYPE_DATA);

static void
nc_data_cluster_mass_rich_init (NcDataClusterMassRich *dmr)
{
  NcDataClusterMassRichPrivate * const self = nc_data_cluster_mass_rich_get_instance_private (dmr);

  self->z_cluster   = NULL;
  self->lnM_cluster = NULL;
  self->lnR_cluster = NULL;
}

static void
nc_data_cluster_mass_rich_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterMassRich *dmr                = NC_DATA_CLUSTER_MASS_RICH (object);
  NcDataClusterMassRichPrivate * const self = nc_data_cluster_mass_rich_get_instance_private (dmr);

  g_return_if_fail (NC_IS_DATA_CLUSTER_MASS_RICH (object));

  switch (prop_id)
  {
    case PROP_Z_CLUSTER:
      ncm_vector_clear (&self->z_cluster);
      self->z_cluster = g_value_dup_object (value);
      break;
    case PROP_LNM_CLUSTER:
      ncm_vector_clear (&self->lnM_cluster);
      self->lnM_cluster = g_value_dup_object (value);
      break;
    case PROP_LNR_CLUSTER:
      ncm_vector_clear (&self->lnR_cluster);
      self->lnR_cluster = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_mass_rich_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterMassRich *dmr                = NC_DATA_CLUSTER_MASS_RICH (object);
  NcDataClusterMassRichPrivate * const self = nc_data_cluster_mass_rich_get_instance_private (dmr);

  g_return_if_fail (NC_IS_DATA_CLUSTER_MASS_RICH (object));

  switch (prop_id)
  {
    case PROP_Z_CLUSTER:
      g_value_set_object (value, self->z_cluster);
      break;
    case PROP_LNM_CLUSTER:
      g_value_set_object (value, self->lnM_cluster);
      break;
    case PROP_LNR_CLUSTER:
      g_value_set_object (value, self->lnR_cluster);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_mass_rich_dispose (GObject *object)
{
  NcDataClusterMassRich *dmr                = NC_DATA_CLUSTER_MASS_RICH (object);
  NcDataClusterMassRichPrivate * const self = nc_data_cluster_mass_rich_get_instance_private (dmr);

  ncm_vector_clear (&self->z_cluster);
  ncm_vector_clear (&self->lnM_cluster);
  ncm_vector_clear (&self->lnR_cluster);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_mass_rich_parent_class)->dispose (object);
}

static void
nc_data_cluster_mass_rich_finalize (GObject *object)
{
  /*NcDataClusterMassRich *dmr = NC_DATA_CLUSTER_MASS_RICH (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_mass_rich_parent_class)->finalize (object);
}

static guint _nc_data_cluster_mass_rich_get_length (NcmData *data);
static guint _nc_data_cluster_mass_rich_get_dof (NcmData *data);
static void _nc_data_cluster_mass_rich_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _nc_data_cluster_mass_rich_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_mass_rich_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);

static void
nc_data_cluster_mass_rich_class_init (NcDataClusterMassRichClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = nc_data_cluster_mass_rich_set_property;
  object_class->get_property = nc_data_cluster_mass_rich_get_property;
  object_class->dispose      = nc_data_cluster_mass_rich_dispose;
  object_class->finalize     = nc_data_cluster_mass_rich_finalize;

  g_object_class_install_property (object_class,
                                   PROP_Z_CLUSTER,
                                   g_param_spec_object ("z-cluster",
                                                        NULL,
                                                        "Clusters (halo) redshift array",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_CLUSTER,
                                   g_param_spec_object ("lnM-cluster",
                                                        NULL,
                                                        "Clusters (halo) ln-mass array",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_LNR_CLUSTER,
                                   g_param_spec_object ("lnR-cluster",
                                                        NULL,
                                                        "Clusters (halo) ln-richness array",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap  = FALSE;
  data_class->get_length = &_nc_data_cluster_mass_rich_get_length;
  data_class->get_dof    = &_nc_data_cluster_mass_rich_get_dof;
  data_class->m2lnL_val  = &_nc_data_cluster_mass_rich_m2lnL_val;
  data_class->prepare    = &_nc_data_cluster_mass_rich_prepare;
  data_class->resample   = &_nc_data_cluster_mass_rich_resample;
}

static guint
_nc_data_cluster_mass_rich_get_length (NcmData *data)
{
  NcDataClusterMassRich *dmr                = NC_DATA_CLUSTER_MASS_RICH (data);
  NcDataClusterMassRichPrivate * const self = nc_data_cluster_mass_rich_get_instance_private (dmr);

  return ncm_vector_len (self->z_cluster);
}

static guint
_nc_data_cluster_mass_rich_get_dof (NcmData *data)
{
  return _nc_data_cluster_mass_rich_get_length (data);
}

static void
_nc_data_cluster_mass_rich_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterMassRich *dmr                = NC_DATA_CLUSTER_MASS_RICH (data);
  NcDataClusterMassRichPrivate * const self = nc_data_cluster_mass_rich_get_instance_private (dmr);
  NcClusterMass *cluster_mass               = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  register gdouble local_m2lnL              = 0.0;

  if (NC_IS_CLUSTER_MASS_ASCASO (cluster_mass))
  {
    NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (ncm_mset_peek (mset, nc_cluster_mass_id ()));
    const guint ncluster        = ncm_vector_len (self->z_cluster);
    guint i;

    for (i = 0; i < ncluster; i++)
    {
      const gdouble z_i        = ncm_vector_get (self->z_cluster, i);
      const gdouble lnM_i      = ncm_vector_get (self->lnM_cluster, i);
      const gdouble lnR_i      = ncm_vector_get (self->lnR_cluster, i);
      const gdouble lnR_i_mean = nc_cluster_mass_ascaso_get_mean_richness (ascaso, lnM_i, z_i);
      const gdouble lnR_i_std  = nc_cluster_mass_ascaso_get_std_richness (ascaso, lnM_i, z_i);
      const gdouble lnR_cut_i  = nc_cluster_mass_ascaso_get_cut (ascaso, lnM_i, z_i);

      local_m2lnL += gsl_pow_2 ((lnR_i - lnR_i_mean) / lnR_i_std)
                     + 2.0 * log (lnR_i_std)
                     + 2.0 * gsl_sf_log_erfc ((lnR_cut_i - lnR_i_mean) / (M_SQRT2 * lnR_i_std));
    }
  }
  else if (NC_IS_CLUSTER_MASS_LNRICH_EXT (cluster_mass))
  {
    NcClusterMassLnrichExt *lnrich_ext = NC_CLUSTER_MASS_LNRICH_EXT (ncm_mset_peek (mset, nc_cluster_mass_id ()));
    const guint ncluster               = ncm_vector_len (self->z_cluster);
    guint i;

    for (i = 0; i < ncluster; i++)
    {
      const gdouble z_i        = ncm_vector_get (self->z_cluster, i);
      const gdouble lnM_i      = ncm_vector_get (self->lnM_cluster, i);
      const gdouble lnR_i      = ncm_vector_get (self->lnR_cluster, i);
      const gdouble lnR_i_mean = nc_cluster_mass_lnrich_ext_get_mean_richness (lnrich_ext, lnM_i, z_i);
      const gdouble lnR_i_std  = nc_cluster_mass_lnrich_ext_get_std_richness (lnrich_ext, lnM_i, z_i);
      const gdouble lnR_cut_i  = nc_cluster_mass_lnrich_ext_get_cut (lnrich_ext);

      local_m2lnL += gsl_pow_2 ((lnR_i - lnR_i_mean) / lnR_i_std)
                     + 2.0 * log (lnR_i_std)
                     + 2.0 * gsl_sf_log_erfc ((lnR_cut_i - lnR_i_mean) / (M_SQRT2 * lnR_i_std));
    }
  }
  else
  {
    g_error ("nc_data_cluster_mass_rich_mean_func: unsupported cluster mass model");
  }

  *m2lnL = local_m2lnL;
}

static void
_nc_data_cluster_mass_rich_prepare (NcmData *data, NcmMSet *mset)
{
  /*NcDataClusterMassRich *dmr               = NC_DATA_CLUSTER_MASS_RICH (data);*/
  NcClusterMass *cmass = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));

  /* Currently only compatible with #NcClusterMassAscaso */
  g_assert (NC_IS_CLUSTER_MASS_ASCASO (cmass) || NC_IS_CLUSTER_MASS_LNRICH_EXT (cmass));
}

/**
 * nc_data_cluster_mass_rich_new:
 *
 * Creates a new #NcDataClusterMassRich.
 *
 * Returns: the newly created #NcDataClusterMassRich.
 */
NcDataClusterMassRich *
nc_data_cluster_mass_rich_new (void)
{
  NcDataClusterMassRich *dmr = g_object_new (NC_TYPE_DATA_CLUSTER_MASS_RICH,
                                             NULL);

  return dmr;
}

/**
 * nc_data_cluster_mass_rich_ref:
 * @dmr: a #NcDataClusterMassRich
 *
 * Increases the reference count of @dmr by one.
 *
 * Returns: (transfer full): @dmr
 */
NcDataClusterMassRich *
nc_data_cluster_mass_rich_ref (NcDataClusterMassRich *dmr)
{
  return g_object_ref (dmr);
}

/**
 * nc_data_cluster_mass_rich_free:
 * @dmr: a #NcDataClusterMassRich
 *
 * Atomically decrements the reference count of @dmr by one. If the reference count drops to 0,
 * all memory allocated by @dmr is released.
 *
 */
void
nc_data_cluster_mass_rich_free (NcDataClusterMassRich *dmr)
{
  g_object_unref (dmr);
}

/**
 * nc_data_cluster_mass_rich_clear:
 * @dmr: a #NcDataClusterMassRich
 *
 * If @dmr is not %NULL, unrefs it and sets the pointer to %NULL.
 *
 */
void
nc_data_cluster_mass_rich_clear (NcDataClusterMassRich **dmr)
{
  g_clear_object (dmr);
}

/**
 * nc_data_cluster_mass_rich_set_data:
 * @dmr: a #NcDataClusterMassRich
 * @lnM: a #NcmVector
 * @z: a #NcmVector
 * @lnR: a #NcmVector
 *
 * Sets the data of @dmr.
 *
 */
void
nc_data_cluster_mass_rich_set_data (NcDataClusterMassRich *dmr, NcmVector *lnM, NcmVector *z, NcmVector *lnR)
{
  NcDataClusterMassRichPrivate * const self = nc_data_cluster_mass_rich_get_instance_private (dmr);
  const guint len                           = ncm_vector_len (lnM);

  if ((len != ncm_vector_len (z)) || (len != ncm_vector_len (lnR)))
  {
    g_error ("nc_data_cluster_mass_rich_set_data: lnM, z and lnR must have the same length");

    return;
  }

  ncm_vector_clear (&self->z_cluster);
  ncm_vector_clear (&self->lnM_cluster);
  ncm_vector_clear (&self->lnR_cluster);

  self->z_cluster   = ncm_vector_dup (z);
  self->lnM_cluster = ncm_vector_dup (lnM);
  self->lnR_cluster = ncm_vector_dup (lnR);

  ncm_data_set_init (NCM_DATA (dmr), TRUE);
}

static void
_nc_data_cluster_mass_rich_resample (NcmData  *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataClusterMassRich *dmr             = NC_DATA_CLUSTER_MASS_RICH (data);
  NcDataClusterMassRichPrivate * const self = nc_data_cluster_mass_rich_get_instance_private (dmr);

  NcHICosmo *cosmo                        = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterMass *clusterm                 = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));

  guint i;
  guint np = _nc_data_cluster_mass_rich_get_length (data);
  gdouble *lnM_obs;
  gdouble lnM;
  gdouble z;
  NcmVector *lnR_vec = ncm_vector_new(np);
  NcmVector *z_vec   = ncm_vector_new(np);
  NcmVector *lnM_vec = ncm_vector_new(np);
    
  _nc_data_cluster_mass_rich_prepare (data, mset);

  ncm_rng_lock (rng);

  for (i = 0; i < np; i++)
  { 
    lnM_obs[0] = -1;
    lnM        = ncm_vector_get(self->lnM_cluster , i);
    z          = ncm_vector_get(self->z_cluster , i);
      
    
    if (NC_IS_CLUSTER_MASS_ASCASO (clusterm))
  {

   while(nc_cluster_mass_resample (clusterm,  cosmo, lnM, z, lnM_obs, NULL, rng) == FALSE)
       {
        nc_cluster_mass_resample (clusterm,  cosmo, lnM , z, lnM_obs, NULL, rng);   
       }
         
  }

    
   else if (NC_IS_CLUSTER_MASS_LNRICH_EXT (clusterm))
  {
   NcClusterMassLnrichExt *lnrich_ext = NC_CLUSTER_MASS_LNRICH_EXT (ncm_mset_peek (mset, nc_cluster_mass_id ()));
   
   while(nc_cluster_mass_resample (clusterm,  cosmo, lnM, z, lnM_obs, NULL, rng) ==FALSE)
       {
           while(lnM_obs[0] < nc_cluster_mass_lnrich_ext_get_cut(lnrich_ext))
        { 
           nc_cluster_mass_resample (clusterm,  cosmo, lnM, z, lnM_obs, NULL, rng);
        }
       }
   
  }

    ncm_vector_set(lnR_vec , i  ,lnM_obs[0]);
    ncm_vector_set(lnM_vec , i  ,lnM);
    ncm_vector_set(z_vec , i  ,z);  
      

      

  ncm_rng_unlock (rng);


  
}
nc_data_cluster_mass_rich_set_data (dmr , lnM_vec, z_vec, lnR_vec);
}