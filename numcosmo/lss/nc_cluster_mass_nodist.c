/***************************************************************************
 *            nc_cluster_mass_nodist.c
 *
 *  Fri June 22 13:42:49 2012
 *  Copyright  2012  Mariana Penna Lima, Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>, <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
 * Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:nc_cluster_mass_nodist
 * @title: NcClusterMassNodist
 * @short_description: Cluster mass real mass distribution.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_nodist.h"
#include "math/ncm_cfg.h"

struct _NcClusterMassNodistPrivate
{
  gdouble lnM_min;
  gdouble lnM_max;
  gdouble norma;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterMassNodist, nc_cluster_mass_nodist, NC_TYPE_CLUSTER_MASS);

enum
{
  PROP_0,
  PROP_LNM_MIN,
  PROP_LNM_MAX,
  PROP_SIZE,
};

static void
nc_cluster_mass_nodist_init (NcClusterMassNodist *mn)
{
  NcClusterMassNodistPrivate * const self = mn->priv = nc_cluster_mass_nodist_get_instance_private (mn);
  
  self->lnM_min = 0.0;
  self->lnM_max = 0.0;
  self->norma   = 0.0;
}

static void
_nc_cluster_mass_nodist_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassNodist *mnodist            = NC_CLUSTER_MASS_NODIST (object);
  NcClusterMassNodistPrivate * const self = mnodist->priv;
  
  g_return_if_fail (NC_IS_CLUSTER_MASS_NODIST (object));
  
  switch (prop_id)
  {
    case PROP_LNM_MIN:
      self->lnM_min = g_value_get_double (value);
      self->norma   = 1.0 / (self->lnM_max - self->lnM_min);
      break;
    case PROP_LNM_MAX:
      self->lnM_max = g_value_get_double (value);
      self->norma   = 1.0 / (self->lnM_max - self->lnM_min);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_nodist_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassNodist *mnodist            = NC_CLUSTER_MASS_NODIST (object);
  NcClusterMassNodistPrivate * const self = mnodist->priv;
  
  g_return_if_fail (NC_IS_CLUSTER_MASS_NODIST (object));
  
  switch (prop_id)
  {
    case PROP_LNM_MIN:
      g_value_set_double (value, self->lnM_min);
      break;
    case PROP_LNM_MAX:
      g_value_set_double (value, self->lnM_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_mass_nodist_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_nodist_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_nodist_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
static gdouble _nc_cluster_mass_nodist_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_nodist_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
static gboolean _nc_cluster_mass_nodist_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_nodist_p_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_nodist_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_nodist_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnm_lower, gdouble *lnm_upper);
static gdouble _nc_cluster_mass_nodist_volume (NcClusterMass *clusterm);

static void
nc_cluster_mass_nodist_class_init (NcClusterMassNodistClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcClusterMassClass *parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);
  
  object_class->finalize =  &nc_cluster_mass_nodist_finalize;
  
  model_class->set_property = &_nc_cluster_mass_nodist_set_property;
  model_class->get_property = &_nc_cluster_mass_nodist_get_property;
  
  ncm_model_class_set_name_nick (model_class, "No mass distribution", "No_distribution");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);
  
  /**
   * NcClusterMassNodist:lnM_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_LNM_MIN,
                                   g_param_spec_double ("lnM-min",
                                                        NULL,
                                                        "Minimum mass",
                                                        11.0 * M_LN10, G_MAXDOUBLE, log (5.0) + 13.0 * M_LN10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcClusterMassNodist:lnM_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_LNM_MAX,
                                   g_param_spec_double ("lnM-max",
                                                        NULL,
                                                        "Maximum mass",
                                                        11.0 * M_LN10, G_MAXDOUBLE, 16.0 * M_LN10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  parent_class->P              = &_nc_cluster_mass_nodist_p;
  parent_class->intP           = &_nc_cluster_mass_nodist_intp;
  parent_class->intP_bin       = &_nc_cluster_mass_nodist_intp_bin;
  parent_class->resample       = &_nc_cluster_mass_nodist_resample;
  parent_class->P_limits       = &_nc_cluster_mass_nodist_p_limits;
  parent_class->P_bin_limits   = &_nc_cluster_mass_nodist_p_bin_limits;
  parent_class->N_limits       = &_nc_cluster_mass_nodist_n_limits;
  parent_class->volume         = &_nc_cluster_mass_nodist_volume;
  parent_class->obs_len        = 1;
  parent_class->obs_params_len = 0;
  
  ncm_model_class_add_impl_opts (model_class, NC_CLUSTER_MASS_N_LIMITS, NC_CLUSTER_MASS_RESAMPLE, -1);
}

static gdouble
_nc_cluster_mass_nodist_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{
  g_error ("This object don't implement p.");
  
  return GSL_NAN;
}

static gdouble
_nc_cluster_mass_nodist_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  /*NcClusterMassNodist *mnodist            = NC_CLUSTER_MASS_NODIST (clusterm);*/
  /*NcClusterMassNodistPrivate * const self = mnodist->priv;*/
  
  return 1.0;
}

static gdouble
_nc_cluster_mass_nodist_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params)
{
  /*NcClusterMassNodist *mnodist            = NC_CLUSTER_MASS_NODIST (clusterm);*/
  /*NcClusterMassNodistPrivate * const self = mnodist->priv;*/
  
  if ((lnM < lnM_obs_lower[0]) || (lnM > lnM_obs_upper[0]))
    return 0.0;
  else
    return 1.0;
}

static gboolean
_nc_cluster_mass_nodist_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassNodist *mnodist            = NC_CLUSTER_MASS_NODIST (clusterm);
  NcClusterMassNodistPrivate * const self = mnodist->priv;
  
  lnM_obs[0] = lnM;
  
  return (lnM_obs[0] <= self->lnM_max) && (lnM_obs[0] >= self->lnM_min);
}

static void
_nc_cluster_mass_nodist_p_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassNodist *mnodist            = NC_CLUSTER_MASS_NODIST (clusterm);
  NcClusterMassNodistPrivate * const self = mnodist->priv;
  
  *lnM_lower = self->lnM_min;
  *lnM_upper = self->lnM_max;
  
  return;
}

static void
_nc_cluster_mass_nodist_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  *lnM_lower = lnM_obs_lower[0];
  *lnM_upper = lnM_obs_upper[0];
  
  return;
}

static void
_nc_cluster_mass_nodist_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnm_lower, gdouble *lnm_upper)
{
  NcClusterMassNodist *mnodist            = NC_CLUSTER_MASS_NODIST (clusterm);
  NcClusterMassNodistPrivate * const self = mnodist->priv;
  
  *lnm_lower = self->lnM_min;
  *lnm_upper = self->lnM_max;
  
  return;
}

static gdouble
_nc_cluster_mass_nodist_volume (NcClusterMass *clusterm)
{
  return 1.0;
}
