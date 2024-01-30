/***************************************************************************
 *            nc_cluster_mass_lnrich_ext.c
 *
 *  Tue Oct 31 16:15:11 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti and Cinthia Nunes de Lima
 *  <vitenti@uel.br>, <cinthia.n.lima@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti and Cinthia Nunes de Lima 2023
 * <vitenti@uel.br>, <cinthia.n.lima@uel.br>
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
 * SECTION:nc_cluster_mass_lnrich_ext
 * @title: NcClusterMassLnrichExt
 * @short_description: Extended ln-richness proxy
 *
 * Ln-richness proxy with extentended parametrization.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_lnrich_ext.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcClusterMassLnrichExtPrivate
{
  gdouble M0;
  gdouble z0;
  gdouble lnM0;
  gdouble ln1pz0;
  gdouble lnR_max;
  gdouble lnR_min;
  gboolean use_ln1pz;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterMassLnrichExt, nc_cluster_mass_lnrich_ext, NC_TYPE_CLUSTER_MASS);

#define VECTOR   (NCM_MODEL (lnrich_ext))
#define MU       (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_MU))
#define MU_M1    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_MU_M1))
#define MU_Z1    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_MU_Z1))
#define MU_M2    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_MU_M2))
#define MU_Z2    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_MU_Z2))
#define MU_MZ    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_MU_MZ))
#define SIGMA_0  (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_0))
#define SIGMA_M1 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_M1))
#define SIGMA_Z1 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_Z1))
#define SIGMA_M2 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_M2))
#define SIGMA_Z2 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_Z2))
#define SIGMA_MZ (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_MZ))
#define CUT      (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_CUT))
#define CUT_M1   (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_CUT_M1))
#define CUT_Z1   (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_LNRICH_EXT_CUT_Z1))

enum
{
  PROP_0,
  PROP_M0,
  PROP_Z0,
  PROP_LNRICHNESS_MIN,
  PROP_LNRICHNESS_MAX,
  PROP_USE_LN1PZ,
  PROP_SIZE,
};

static void
nc_cluster_mass_lnrich_ext_init (NcClusterMassLnrichExt *lnrich_ext)
{
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv = nc_cluster_mass_lnrich_ext_get_instance_private (lnrich_ext);

  self->M0        = 0.0;
  self->z0        = 0.0;
  self->lnM0      = 0.0;
  self->ln1pz0    = 0.0;
  self->lnR_min   = GSL_NEGINF;
  self->lnR_max   = GSL_POSINF;
  self->use_ln1pz = FALSE;
}

static void
_nc_cluster_mass_lnrich_ext_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassLnrichExt *lnrich_ext         = NC_CLUSTER_MASS_LNRICH_EXT (object);
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv;

  g_return_if_fail (NC_IS_CLUSTER_MASS_LNRICH_EXT (object));

  switch (prop_id)
  {
    case PROP_M0:
      self->M0   = g_value_get_double (value);
      self->lnM0 = log (self->M0);
      break;
    case PROP_Z0:
      self->z0     = g_value_get_double (value);
      self->ln1pz0 = log1p (self->z0);
      break;
    case PROP_LNRICHNESS_MIN:
      self->lnR_min = g_value_get_double (value);
      g_assert (self->lnR_min < self->lnR_max);
      break;
    case PROP_LNRICHNESS_MAX:
      self->lnR_max = g_value_get_double (value);
      g_assert (self->lnR_min < self->lnR_max);
      break;
    case PROP_USE_LN1PZ:
      self->use_ln1pz = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_lnrich_ext_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassLnrichExt *lnrich_ext         = NC_CLUSTER_MASS_LNRICH_EXT (object);
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv;

  g_return_if_fail (NC_IS_CLUSTER_MASS_LNRICH_EXT (object));

  switch (prop_id)
  {
    case PROP_M0:
      g_value_set_double (value, self->M0);
      break;
    case PROP_Z0:
      g_value_set_double (value, self->z0);
      break;
    case PROP_LNRICHNESS_MIN:
      g_value_set_double (value, self->lnR_min);
      break;
    case PROP_LNRICHNESS_MAX:
      g_value_set_double (value, self->lnR_max);
      break;
    case PROP_USE_LN1PZ:
      g_value_set_boolean (value, self->use_ln1pz);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_lnrich_ext_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_lnrich_ext_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_lnrich_ext_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
static gdouble _nc_cluster_mass_lnrich_ext_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_lnrich_ext_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
static gboolean _nc_cluster_mass_lnrich_ext_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_lnrich_ext_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_lnrich_ext_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_lnrich_ext_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);
static gdouble _nc_cluster_mass_lnrich_ext_volume (NcClusterMass *clusterm);

static void
nc_cluster_mass_lnrich_ext_class_init (NcClusterMassLnrichExtClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcClusterMassClass *parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_cluster_mass_lnrich_ext_set_property;
  model_class->get_property = &_nc_cluster_mass_lnrich_ext_get_property;
  object_class->finalize    = &_nc_cluster_mass_lnrich_ext_finalize;

  ncm_model_class_set_name_nick (model_class, "LnrichExt Ln-normal richness distribution", "LnrichExt");
  ncm_model_class_add_params (model_class, NC_CLUSTER_MASS_LNRICH_EXT_SPARAM_LEN, 0, PROP_SIZE);


  /**
   * NcClusterMassLnrichExt:M0:
   *
   * Pivot mass FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_M0,
                                   g_param_spec_double ("M0",
                                                        NULL,
                                                        "Pivot mass",
                                                        11.0 * M_LN10, G_MAXDOUBLE, 3.0e14 / 0.71,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

/**
 * NcClusterMassLnrichExt:Z0:
 *
 * Pivot redshift FIXME Set correct values (limits)
 */
  g_object_class_install_property (object_class,
                                   PROP_Z0,
                                   g_param_spec_double ("z0",
                                                        NULL,
                                                        "Pivot redshift",
                                                        0.0, G_MAXDOUBLE, 0.6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));



  /**
   * NcClusterMassLnrichExt:lnRichness_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_LNRICHNESS_MIN,
                                   g_param_spec_double ("lnRichness-min",
                                                        NULL,
                                                        "Minimum LnRichness",
                                                        0.0, G_MAXDOUBLE, M_LN10 * 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassLnrichExt:lnRichness_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_LNRICHNESS_MAX,
                                   g_param_spec_double ("lnRichness-max",
                                                        NULL,
                                                        "Maximum LnRichness",
                                                        0.0, G_MAXDOUBLE,  M_LN10 * 2.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassLnrichExt:use_ln1pz:
   *
   * Whether we are going to use ln(1+z) or z in the fits.
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_LN1PZ,
                                   g_param_spec_boolean ("use-ln1pz",
                                                         NULL,
                                                         "Whether we are going to use ln(1+z) or z in the fits.",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassLnrichExt:MU:
   *
   * Distribution's  bias in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_MU, "\\mu", "mu",
                              0.0,  6.0, 1.0e-1,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:MU_M1:
   *
   * Distribution's slope with respect to the mass in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_MU_M1, "\\mu_\\mathrm{M_1}", "muM1",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_M1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:MU_Z1:
   *
   * Distribution's slope with respect to the redshift in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_MU_Z1, "\\mu_\\mathrm{z_1}", "muZ1",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_Z1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:MU_M2:
   *
   * Distribution's quadratic slope with respect to the mass in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_MU_M2, "\\mu_\\mathrm{M_2}", "muM2",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_M2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:MU_Z2:
   *
   * Distribution's quadratic slope with respect to the redshift in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_MU_Z2, "\\mu_\\mathrm{z_2}", "muZ2",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_Z2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:MU_MZ:
   *
   * Distribution's cross term with respect to the mass and redshift in the mean.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_MU_MZ, "\\mu_\\mathrm{MZ}", "muMZ",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_MZ,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:sigma_0:
   *
   * Distribution's bias in the standard deviation, $\sigma \in [10^{-4}, 10]$.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_0, "\\sigma_0", "sigma0",
                              1.0e-4, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:sigma_M1:
   *
   * Distribution's slope with respect to the mass in the standard deviation.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_M1, "\\sigma_\\mathrm{M_1}", "sigmaM1",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_M1,
                              NCM_PARAM_TYPE_FIXED);


  /**
   * NcClusterMassLnrichExt:sigma_Z1:
   *
   * Distribution's slope with respect to the redshift in the standard deviation.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_Z1, "\\sigma_\\mathrm{z_1}", "sigmaZ1",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_Z1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:sigma_M2:
   *
   * Distribution's quadratic slope with respect to the mass in the standard deviation.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_M2, "\\sigma_\\mathrm{M_2}", "sigmaM2",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_M2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:sigma_Z2:
   *
   * Distribution's quadratic slope with respect to the redshift in the standard deviation.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_Z2, "\\sigma_\\mathrm{z_2}", "sigmaZ2",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_Z2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:sigma_MZ:
   *
   * Distribution's cross term with respect to the mass and redshift in the standard deviation.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_MZ, "\\sigma_\\mathrm{MZ}", "sigmaMZ",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_MZ,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:cut:
   *
   * Cut in richness.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_CUT, "cut", "cut",
                              0.0, 1.0e16, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_CUT,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:cut_M1:
   *
   * Cut in richness.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_CUT_M1, "\\cut_\\mathrm{M_1}", "cutM1",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_CUT_M1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassLnrichExt:cut_Z1:
   *
   * Cut in richness.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_LNRICH_EXT_CUT_Z1, "\\cut_\\mathrm{Z_1}", "cutZ1",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_CUT_Z1,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->P               = &_nc_cluster_mass_lnrich_ext_p;
  parent_class->intP            = &_nc_cluster_mass_lnrich_ext_intp;
  parent_class->intP_bin        = &_nc_cluster_mass_lnrich_ext_intp_bin;
  parent_class->resample        = &_nc_cluster_mass_lnrich_ext_resample;
  parent_class->P_limits        = &_nc_cluster_mass_lnrich_ext_p_limits;
  parent_class->P_bin_limits    = &_nc_cluster_mass_lnrich_ext_p_bin_limits;
  parent_class->N_limits        = &_nc_cluster_mass_lnrich_ext_n_limits;
  parent_class->volume          = &_nc_cluster_mass_lnrich_ext_volume;
  parent_class->P_vec_z_lnMobs  = NULL;
  parent_class->_obs_len        = 1;
  parent_class->_obs_params_len = 0;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);
}

static void
_nc_cluster_mass_lnrich_ext_lnR_sigma (NcClusterMass *clusterm, const gdouble lnM, const gdouble z, gdouble *lnR, gdouble *sigma)
{
  NcClusterMassLnrichExt *lnrich_ext         = NC_CLUSTER_MASS_LNRICH_EXT (clusterm);
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv;
  const gdouble DlnM                         = lnM - self->lnM0;
  const gdouble Vz                           = self->use_ln1pz ? log1p (z) - self->ln1pz0 : z / self-> z0;

  lnR[0] = MU + MU_M1 * DlnM + MU_Z1 * Vz + MU_M2 * DlnM * DlnM + MU_Z2 * Vz * Vz + MU_MZ * DlnM * Vz;
  sigma[0] = SIGMA_0 + SIGMA_M1 * DlnM + SIGMA_Z1 * Vz + SIGMA_M2 * DlnM * DlnM + SIGMA_Z2 * Vz * Vz + SIGMA_MZ * DlnM * Vz;   
}

static gdouble
_nc_cluster_mass_lnrich_ext_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{
  /*NcClusterMassLnrichExt *lnrich_ext    = NC_CLUSTER_MASS_LNRICH_EXT (clusterm);*/
  gdouble lnR_true, sigma;

  _nc_cluster_mass_lnrich_ext_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  {
    const gdouble x = (lnM_obs[0] - lnR_true) / sigma;

    return 1.0 / (ncm_c_sqrt_2pi () * sigma) * exp (-0.5 * x * x);
  }
}

static gdouble
_nc_cluster_mass_lnrich_ext_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcClusterMassLnrichExt *lnrich_ext         = NC_CLUSTER_MASS_LNRICH_EXT (clusterm);
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv;
  gdouble lnR_true, sigma;

  _nc_cluster_mass_lnrich_ext_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  {
    const gdouble x_min = (lnR_true - self->lnR_min) / (M_SQRT2 * sigma);
    const gdouble x_max = (lnR_true - self->lnR_max) / (M_SQRT2 * sigma);

    if (x_max > 4.0)
      return -(erfc (x_min) - erfc (x_max)) / 2.0;
    else
      return (erf (x_min) - erf (x_max)) / 2.0;
  }
}

static gdouble
_nc_cluster_mass_lnrich_ext_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params)
{
  /*NcClusterMassLnrichExt *lnrich_ext = NC_CLUSTER_MASS_LNRICH_EXT (clusterm);*/
  gdouble lnR_true, sigma;

  _nc_cluster_mass_lnrich_ext_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  {
    const gdouble x_min = (lnR_true - lnM_obs_lower[0]) / (M_SQRT2 * sigma);
    const gdouble x_max = (lnR_true - lnM_obs_upper[0]) / (M_SQRT2 * sigma);

    if (x_max > 4.0)
      return -(erfc (x_min) - erfc (x_max)) / 2.0;
    else
      return (erf (x_min) - erf (x_max)) / 2.0;
  }
}

static gboolean
_nc_cluster_mass_lnrich_ext_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassLnrichExt *lnrich_ext         = NC_CLUSTER_MASS_LNRICH_EXT (clusterm);
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv;
  gdouble lnR_true, sigma;

  _nc_cluster_mass_lnrich_ext_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

  ncm_rng_lock (rng);
  lnM_obs[0] = lnR_true + gsl_ran_gaussian (rng->r, sigma);
  ncm_rng_unlock (rng);

  return (lnM_obs[0] <= self->lnR_max) && (lnM_obs[0] >= self->lnR_min);
}

static void
_nc_cluster_mass_lnrich_ext_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassLnrichExt *lnrich_ext = NC_CLUSTER_MASS_LNRICH_EXT (clusterm);
  const gdouble mean                 = lnM_obs[0] - MU; /* - P2 * log10(1.0 + z);  FIX This!!!! What is the mean richeness? */
  const gdouble logRichnessl         = M_LN10 * log10 (1e13);
  const gdouble logRichnessu         = M_LN10 * log10 (1e15);

  *lnM_lower = logRichnessl;
  *lnM_upper = logRichnessu;

  return;
}

static void
_nc_cluster_mass_lnrich_ext_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassLnrichExt *lnrich_ext = NC_CLUSTER_MASS_LNRICH_EXT (clusterm);
  const gdouble lnMl                 = M_LN10 * log10 (1e13);
  const gdouble lnMu                 = M_LN10 * log10 (1e15);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;
}

static void
_nc_cluster_mass_lnrich_ext_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassLnrichExt *lnrich_ext = NC_CLUSTER_MASS_LNRICH_EXT (clusterm);
  const gdouble lnMl                 =  M_LN10 * log10 (1e13);
  const gdouble lnMu                 =  M_LN10 * log10 (1e15);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}

static gdouble
_nc_cluster_mass_lnrich_ext_volume (NcClusterMass *clusterm)
{
  NcClusterMassLnrichExt *lnrich_ext         = NC_CLUSTER_MASS_LNRICH_EXT (clusterm);
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv;

  return (self->lnR_max - self->lnR_min);
}

/**
 * nc_cluster_mass_lnrich_ext_get_mean_richness:
 * @lnrich_ext: a #NcClusterMassLnrichExt
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the mean of the richness distribution.
 *
 */
gdouble
nc_cluster_mass_lnrich_ext_get_mean_richness (NcClusterMassLnrichExt *lnrich_ext, gdouble lnM, gdouble z)
{
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv;
  const gdouble DlnM                         = lnM - self->lnM0;
  const gdouble Vz                           = self->use_ln1pz ? log1p (z) - self->ln1pz0 : z / self-> z0;
    
  return MU + MU_M1 * DlnM + MU_Z1 * Vz + MU_M2 * DlnM * DlnM + MU_Z2 * Vz * Vz + MU_MZ * DlnM * Vz;
}

/**
 * nc_cluster_mass_lnrich_ext_get_std_richness:
 * @lnrich_ext: a #NcClusterMassLnrichExt
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the standard deviation of the richness distribution.
 *
 */
gdouble
nc_cluster_mass_lnrich_ext_get_std_richness (NcClusterMassLnrichExt *lnrich_ext, gdouble lnM, gdouble z)
{
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv;
  const gdouble DlnM                         = lnM - self->lnM0;
  const gdouble Vz                           = self->use_ln1pz ? log1p (z) - self->ln1pz0 : z / self-> z0;

  return SIGMA_0 + ( SIGMA_M1 * DlnM ) + ( SIGMA_Z1 * Vz ) + ( SIGMA_M2 * DlnM * DlnM ) + ( SIGMA_Z2 * Vz * Vz ) + ( SIGMA_MZ * DlnM * Vz );
}

/**
 * nc_cluster_mass_lnrich_ext_get_cut:
 * @lnrich_ext: a #NcClusterMassLnrichExt
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the cut in richness.
 *
 * Returns: the cut in richness.
 */
gdouble
nc_cluster_mass_lnrich_ext_get_cut (NcClusterMassLnrichExt *lnrich_ext, gdouble lnM, gdouble z)
{
  NcClusterMassLnrichExtPrivate * const self = lnrich_ext->priv;
  const gdouble DlnM                         = lnM - self->lnM0;
  const gdouble Dln1pz                       = log1p (z) - self->ln1pz0;

  return CUT + CUT_M1 * DlnM + CUT_Z1 * Dln1pz;
}

