/***************************************************************************
 *            nc_cluster_mass_ext.c
 *
 *  Tue Dez 02 18:25:11 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2025 <vitenti@uel.br>
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
 * NcClusterMassExt:
 *
 * Extended cluster mass-richness distribution model.
 *
 * This class implements an extended mass-richness relation, which models the
 * mean and scatter of the log-richness distribution as polynomial functions
 * of mass. The model is suitable for fitting over an extended domain of log-richness.
 *
 * The mean log-richness is given by:
 * $$
 * \mu = \mu_{p0} + \mu_{p1} \Delta\ln M + \mu_{p2} (\Delta\ln M)^2 + \mu_{p3} (\Delta\ln M)^3
 * $$
 * where $\Delta\ln M = \ln M - \ln M_0$.
 *
 * The standard deviation is modeled as a function of the mean:
 * $$
 * \sigma = \exp(\sigma_{p0} + \sigma_{p1} \mu + \sigma_{p2} \mu^2)
 * $$
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_ext.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

typedef struct _NcClusterMassExtPrivate
{
  guint placeholder;
} NcClusterMassExtPrivate;

struct _NcClusterMassExt
{
  /*< private >*/
  NcClusterMassRichness parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterMassExt, nc_cluster_mass_ext, NC_TYPE_CLUSTER_MASS_RICHNESS)

#define VECTOR   (NCM_MODEL (ext))
#define MU_P0    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_EXT_MU_P0))
#define MU_P1    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_EXT_MU_P1))
#define MU_P2    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_EXT_MU_P2))
#define MU_P3    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_EXT_MU_P3))
#define SIGMA_P0 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_EXT_SIGMA_P0))
#define SIGMA_P1 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_EXT_SIGMA_P1))
#define SIGMA_P2 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_EXT_SIGMA_P2))

static void
nc_cluster_mass_ext_init (NcClusterMassExt *ext)
{
  NcClusterMassExtPrivate * const self = nc_cluster_mass_ext_get_instance_private (ext);

  self->placeholder = 0;
}

static void
_nc_cluster_mass_ext_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_ext_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_ext_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_ext_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z);

static void
nc_cluster_mass_ext_class_init (NcClusterMassExtClass *klass)
{
  GObjectClass *object_class           = G_OBJECT_CLASS (klass);
  NcClusterMassRichnessClass *mr_class = NC_CLUSTER_MASS_RICHNESS_CLASS (klass);
  NcmModelClass *model_class           = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_cluster_mass_ext_finalize;

  ncm_model_class_set_name_nick (model_class, "Extended Ln-normal richness distribution", "Extended");
  ncm_model_class_add_params (model_class, NC_CLUSTER_MASS_EXT_SPARAM_LEN - NC_CLUSTER_MASS_RICHNESS_SPARAM_LEN, 0, 1);

  /**
   * NcClusterMassExt:MU_P0:
   *
   * Bias parameter in the mean of the richness-mass relation.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_EXT_MU_P0, "mu_p0", "mup0",
                              0.0,  6.0, 1.0e-1,
                              NC_CLUSTER_MASS_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_EXT_DEFAULT_MU_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassExt:MU_P1:
   *
   * Mass slope parameter in the mean of the richness-mass relation.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_EXT_MU_P1, "mu_p1", "mup1",
                              -10.0,  10.0, 1.0e-1,
                              NC_CLUSTER_MASS_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_EXT_DEFAULT_MU_P1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassExt:MU_P2:
   *
   * Quadratic mass coefficient in the mean of the richness-mass relation, i.e.,
   * the coefficient of $(\Delta\ln M)^2$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_EXT_MU_P2, "mu_p2", "mup2",
                              -10.0,  10.0, 1.0e-1,
                              NC_CLUSTER_MASS_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_EXT_DEFAULT_MU_P2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassExt:MU_P3:
   *
   * Cubic mass coefficient in the mean of the richness-mass relation, i.e.,
   * the coefficient of $(\Delta\ln M)^3$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_EXT_MU_P3, "mu_p3", "mup3",
                              -10.0,  10.0, 1.0e-1,
                              NC_CLUSTER_MASS_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_EXT_DEFAULT_MU_P3,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassExt:sigma_P0:
   *
   * Constant term in the logarithm of the standard deviation, i.e., $\ln\sigma = \sigma_{p0} + ...$.
   * This sets the baseline scatter of the distribution.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_EXT_SIGMA_P0, "\\sigma_p0", "sigmap0",
                              -10.0, 10.0, 1.0e-1,
                              NC_CLUSTER_MASS_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_EXT_DEFAULT_SIGMA_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassExt:sigma_P1:
   *
   * Linear coefficient of $\mu$ in the logarithm of the standard deviation, i.e.,
   * $\ln\sigma = ... + \sigma_{p1} \mu + ...$. Controls how scatter varies with
   * mean log-richness.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_EXT_SIGMA_P1, "\\sigma_p1", "sigmap1",
                              -10.0, 10.0, 1.0e-1,
                              NC_CLUSTER_MASS_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_EXT_DEFAULT_SIGMA_P1,
                              NCM_PARAM_TYPE_FIXED);


  /**
   * NcClusterMassExt:sigma_P2:
   *
   * Quadratic coefficient of $\mu$ in the logarithm of the standard deviation, i.e.,
   * $\ln\sigma = ... + \sigma_{p2} \mu^2$. Controls the curvature of how scatter
   * varies with mean log-richness.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_EXT_SIGMA_P2, "\\sigma_p2", "sigmap2",
                              -10.0,  10.0, 1.0e-1,
                              NC_CLUSTER_MASS_EXT_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_EXT_DEFAULT_SIGMA_P2,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  /* Set virtual methods for NcClusterMassRichness */
  mr_class->mu    = &_nc_cluster_mass_ext_mu;
  mr_class->sigma = &_nc_cluster_mass_ext_sigma;
}

static gdouble
_nc_cluster_mass_ext_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  NcClusterMassExt *ext = NC_CLUSTER_MASS_EXT (mr);
  const gdouble lnM0    = nc_cluster_mass_richness_lnM0 (mr);
  const gdouble DlnM    = lnM - lnM0;

  return MU_P0 + MU_P1 * DlnM + MU_P2 * DlnM * DlnM + MU_P3 * DlnM * DlnM * DlnM;
}

static gdouble
_nc_cluster_mass_ext_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  NcClusterMassExt *ext = NC_CLUSTER_MASS_EXT (mr);
  const gdouble mean    = _nc_cluster_mass_ext_mu (mr, lnM, z);

  return exp (SIGMA_P0 + SIGMA_P1 * mean + SIGMA_P2 * mean * mean);
}

