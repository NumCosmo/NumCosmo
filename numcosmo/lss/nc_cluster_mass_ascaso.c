/***************************************************************************
 *            nc_cluster_mass_ascaso.c
 *
 *  Thu Jan 26 18:25:11 2017
 *  Copyright  2017  Mariana Penna Lima and Begoña Ascaso
 *  <pennalima@gmail.com>, <bego.ascaso.work@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima and Begoña Ascaso 2017 <pennalima@gmail.com>
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
 * NcClusterMassAscaso:
 *
 * Cluster mass-richness distribution model based on Ascaso et al.
 *
 * This class implements the Ascaso mass-richness relation, which models the
 * mean and scatter of the log-richness distribution as linear functions
 * of mass and redshift.
 *
 * The mean log-richness is given by:
 * $$
 * \mu = \mu_{p0} + \mu_{p1} \Delta\ln M + \mu_{p2} \Delta\ln(1+z)
 * $$
 * where $\Delta\ln M = \ln M - \ln M_0$ and $\Delta\ln(1+z) = \ln(1+z) - \ln(1+z_0)$.
 *
 * The standard deviation is given by:
 * $$
 * \sigma = \sigma_{p0} + \sigma_{p1} \Delta\ln M + \sigma_{p2} \Delta\ln(1+z)
 * $$
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_ascaso.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

typedef struct _NcClusterMassAscasoPrivate
{
  guint placeholder;
} NcClusterMassAscasoPrivate;

struct _NcClusterMassAscaso
{
  /*< private >*/
  NcClusterMassRichness parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterMassAscaso, nc_cluster_mass_ascaso, NC_TYPE_CLUSTER_MASS_RICHNESS)

#define VECTOR   (NCM_MODEL (ascaso))
#define MU_P0    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_ASCASO_MU_P0))
#define MU_P1    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_ASCASO_MU_P1))
#define MU_P2    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_ASCASO_MU_P2))
#define MU_P3    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_ASCASO_MU_P3))
#define SIGMA_P0 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_ASCASO_SIGMA_P0))
#define SIGMA_P1 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_ASCASO_SIGMA_P1))
#define SIGMA_P2 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_ASCASO_SIGMA_P2))

static void
nc_cluster_mass_ascaso_init (NcClusterMassAscaso *ascaso)
{
  NcClusterMassAscasoPrivate * const self = nc_cluster_mass_ascaso_get_instance_private (ascaso);

  self->placeholder = 0;
}

static void
_nc_cluster_mass_ascaso_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_ascaso_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_ascaso_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_ascaso_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z);

static void
nc_cluster_mass_ascaso_class_init (NcClusterMassAscasoClass *klass)
{
  GObjectClass *object_class           = G_OBJECT_CLASS (klass);
  NcClusterMassRichnessClass *mr_class = NC_CLUSTER_MASS_RICHNESS_CLASS (klass);
  NcmModelClass *model_class           = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_cluster_mass_ascaso_finalize;

  ncm_model_class_set_name_nick (model_class, "Ascaso Ln-normal richness distribution", "Ascaso");
  ncm_model_class_add_params (model_class, NC_CLUSTER_MASS_ASCASO_SPARAM_LEN - NC_CLUSTER_MASS_RICHNESS_SPARAM_LEN, 0, 1);

  /**
   * NcClusterMassAscaso:MU_P0:
   *
   * Constant term (bias) in the mean log-richness relation.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_MU_P0, "mu_p0", "mup0",
                              0.0,  6.0, 1.0e-1,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassAscaso:MU_P1:
   *
   * Linear mass coefficient in the mean log-richness, i.e., the coefficient of
   * $\Delta\ln M = \ln M - \ln M_0$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_MU_P1, "mu_p1", "mup1",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassAscaso:MU_P2:
   *
   * Redshift evolution coefficient in the mean log-richness, i.e., the coefficient of
   * $\Delta\ln(1+z) = \ln(1+z) - \ln(1+z_0)$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_MU_P2, "mu_p2", "mup2",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_MU_P2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassAscaso:sigma_P0:
   *
   * Constant term (bias) in the standard deviation of the log-richness distribution.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_SIGMA_P0, "\\sigma_p0", "sigmap0",
                              1.0e-4, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassAscaso:sigma_P1:
   *
   * Linear mass coefficient in the standard deviation, i.e., the coefficient of
   * $\Delta\ln M = \ln M - \ln M_0$.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_SIGMA_P1, "\\sigma_p1", "sigmap1",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P1,
                              NCM_PARAM_TYPE_FIXED);


  /**
   * NcClusterMassAscaso:sigma_P2:
   *
   * Redshift evolution coefficient in the standard deviation, i.e., the coefficient of
   * $\Delta\ln(1+z) = \ln(1+z) - \ln(1+z_0)$.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_ASCASO_SIGMA_P2, "\\sigma_p2", "sigmap2",
                              -10.0,  10.0, 1.0e-2,
                              NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA_P2,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  /* Set virtual methods for NcClusterMassRichness */
  mr_class->mu    = &_nc_cluster_mass_ascaso_mu;
  mr_class->sigma = &_nc_cluster_mass_ascaso_sigma;
}

static gdouble
_nc_cluster_mass_ascaso_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (mr);
  const gdouble lnM0          = nc_cluster_mass_richness_lnM0 (mr);
  const gdouble ln1pz0        = nc_cluster_mass_richness_ln1pz0 (mr);
  const gdouble DlnM          = lnM - lnM0;
  const gdouble Dln1pz        = ln1pz0 - ln1pz0;

  return MU_P0 + MU_P1 * DlnM + MU_P2 * Dln1pz;
}

static gdouble
_nc_cluster_mass_ascaso_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  NcClusterMassAscaso *ascaso = NC_CLUSTER_MASS_ASCASO (mr);
  const gdouble lnM0          = nc_cluster_mass_richness_lnM0 (mr);
  const gdouble ln1pz0        = nc_cluster_mass_richness_ln1pz0 (mr);
  const gdouble DlnM          = lnM - lnM0;
  const gdouble Dln1pz        = ln1pz0 - ln1pz0;
  const gdouble sigma         = SIGMA_P0 + SIGMA_P1 * DlnM + SIGMA_P2 * Dln1pz;

  /* Add a small number to the standard deviation to avoid numerical instabilities */
  return hypot (sigma, 1.0e-5);
}

