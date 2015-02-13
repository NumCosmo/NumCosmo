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

G_DEFINE_TYPE (NcClusterMassNodist, nc_cluster_mass_nodist, NC_TYPE_CLUSTER_MASS);

enum
{
  PROP_0,
  PROP_LNM_MIN,
  PROP_LNM_MAX,
  PROP_SIZE,
};

static gdouble
_nc_cluster_mass_nodist_p (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params)
{
  g_error ("This object don't implement p.");
  NCM_UNUSED (clusterm);
  NCM_UNUSED (model);
  NCM_UNUSED (lnM);
  NCM_UNUSED (z);
  NCM_UNUSED (lnM_obs);
  NCM_UNUSED (lnM_obs_params);
  return GSL_NAN;
}

static gdouble
_nc_cluster_mass_nodist_intp (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z)
{
  g_error ("This object don't implement intp.");
  NCM_UNUSED (clusterm);
  NCM_UNUSED (model);
  NCM_UNUSED (lnM);
  NCM_UNUSED (z);
  return GSL_NAN;
}

static gboolean
_nc_cluster_mass_nodist_resample (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassNodist *mnodist = NC_CLUSTER_MASS_NODIST (clusterm);
  lnM_obs[0] = lnM;

  NCM_UNUSED (model);
  NCM_UNUSED (z);
  NCM_UNUSED (lnM_obs_params);

  return (lnM_obs[0] <= mnodist->lnM_max) && (lnM_obs[0] >= mnodist->lnM_min);
}

static void
_nc_cluster_mass_nodist_p_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  g_error ("This object don't implement integ_limits.");

  NCM_UNUSED (clusterm);
  NCM_UNUSED (model);
  NCM_UNUSED (lnM_obs);
  NCM_UNUSED (lnM_obs_params);
  NCM_UNUSED (lnM_lower);
  NCM_UNUSED (lnM_upper);

  return;
}

static void
_nc_cluster_mass_nodist_n_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnm_lower, gdouble *lnm_upper)
{
  NcClusterMassNodist *mn = NC_CLUSTER_MASS_NODIST (clusterm);

  NCM_UNUSED (model);

  *lnm_lower = mn->lnM_min;
  *lnm_upper = mn->lnM_max;

  return;
}

guint _nc_cluster_mass_nodist_obs_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 1; }
guint _nc_cluster_mass_nodist_obs_params_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 0; }

static void
_nc_cluster_mass_nodist_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcClusterMassNodist *mnodist = NC_CLUSTER_MASS_NODIST (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_NODIST (object));

  switch (prop_id)
  {
    case PROP_LNM_MIN:
      mnodist->lnM_min = g_value_get_double (value);
      break;
    case PROP_LNM_MAX:
      mnodist->lnM_max = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_nodist_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassNodist *mnodist = NC_CLUSTER_MASS_NODIST (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_NODIST (object));

  switch (prop_id)
  {
    case PROP_LNM_MIN:
      g_value_set_double (value, mnodist->lnM_min);
      break;
    case PROP_LNM_MAX:
      g_value_set_double (value, mnodist->lnM_max);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_mass_nodist_init (NcClusterMassNodist *mn)
{
  mn->lnM_min = 0.0;
  mn->lnM_max = 0.0;
}

static void
nc_cluster_mass_nodist_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_nodist_parent_class)->finalize (object);
}

static void
nc_cluster_mass_nodist_class_init (NcClusterMassNodistClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterMassClass* parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  parent_class->P              = &_nc_cluster_mass_nodist_p;
  parent_class->intP           = &_nc_cluster_mass_nodist_intp;
  parent_class->resample       = &_nc_cluster_mass_nodist_resample;
  parent_class->P_limits       = &_nc_cluster_mass_nodist_p_limits;
  parent_class->N_limits       = &_nc_cluster_mass_nodist_n_limits;
  parent_class->obs_len        = &_nc_cluster_mass_nodist_obs_len;
  parent_class->obs_params_len = &_nc_cluster_mass_nodist_obs_params_len;

  parent_class->impl = NC_CLUSTER_MASS_N_LIMITS | NC_CLUSTER_MASS_RESAMPLE;

  object_class->finalize     =  &nc_cluster_mass_nodist_finalize;

  model_class->set_property = &_nc_cluster_mass_nodist_set_property;
  model_class->get_property = &_nc_cluster_mass_nodist_get_property;

  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);
  ncm_model_class_set_name_nick (model_class, "No mass distribution", "No_distribution");

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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


}
