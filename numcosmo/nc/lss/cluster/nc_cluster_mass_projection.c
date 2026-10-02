/***************************************************************************
 *            nc_cluster_mass_projection.c
 *
 *  Thu October 02 10:00:00 2026
 *  Copyright  2026  Cinthia N. Lima
 *  <cinthia.n.lima@uel.br>
 ****************************************************************************/
/*
 * nc_cluster_mass_projection.c
 * Copyright (C) 2026 Cinthia N. Lima <cinthia.n.lima@uel.br>
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
 * NcClusterMassProjection:
 *
 * Mass-richness relation contaminated by projection effects.
 *
 * The true richness follows the same log-normal relation as
 * #NcClusterMassAscaso, with
 * $$\mu = \mu_{p0} + \mu_{p1}\Delta\ln M + \mu_{p2}\Delta\ln(1+z), \qquad
 *   \sigma = \sigma_{p0} + \sigma_{p1}\Delta\ln M + \sigma_{p2}\Delta\ln(1+z).$$
 *
 * A fraction $f_\mathrm{prj}$ of the clusters has its richness inflated by
 * unassociated structure along the line of sight, modelled as an independent
 * exponential of rate $\tau$ added to the true richness,
 * $$\lambda = \lambda_\mathrm{true} + \Delta, \qquad \Delta \sim \mathrm{Exp}(\tau).$$
 * The observed richness is then the two-component mixture
 * $$P(\lambda \mid M, z) = (1 - f_\mathrm{prj})\, f_\mathrm{LN}(\lambda \mid \mu, \sigma)
 *                        + f_\mathrm{prj}\, T(\lambda),$$
 * where $T$ is the log-normal/exponential convolution evaluated by
 * #NcClusterRichnessProjection. This covers the first two components of the
 * mixture of Costanzi et al. 2019; percolation and masking are not included.
 *
 * $T$ costs one ODE solve per $(\ln M, z)$, controlled by
 * #NcClusterMassProjection:reltol. The default $10^{-7}$ is far tighter than a
 * likelihood needs and about ten times cheaper than the tolerance
 * #NcClusterRichnessProjection uses on its own.
 *
 * Sampling never touches the ODE: nc_cluster_mass_projection_resample() draws the
 * mixture component, then the log-normal, then adds the exponential.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/cluster/nc_cluster_mass_projection.h"
#include "nc/lss/cluster/nc_cluster_richness_projection.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_rng.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcClusterMassProjectionPrivate
{
  NcClusterRichnessProjection *crp;
  gdouble reltol;
  gdouble mu;
  gdouble sigma;
  gdouble tau;
  gdouble lnl_min;
  gdouble lnl_max;
  gboolean prepared;
} NcClusterMassProjectionPrivate;

struct _NcClusterMassProjection
{
  /*< private >*/
  NcClusterMassRichness parent_instance;
};

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterMassProjection, nc_cluster_mass_projection, NC_TYPE_CLUSTER_MASS_RICHNESS)

#define VECTOR   (NCM_MODEL (mp))
#define MU_P0    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_MU_P0))
#define MU_P1    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_MU_P1))
#define MU_P2    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_MU_P2))
#define SIGMA_P0 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_SIGMA_P0))
#define SIGMA_P1 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_SIGMA_P1))
#define SIGMA_P2 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_SIGMA_P2))
#define F_PRJ    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_F_PRJ))
#define TAU      (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_PROJECTION_TAU))
#define CUT      (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_RICHNESS_CUT))

static void
nc_cluster_mass_projection_init (NcClusterMassProjection *mp)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (mp);

  self->crp      = nc_cluster_richness_projection_new ();
  self->reltol   = NC_CLUSTER_MASS_PROJECTION_DEFAULT_RELTOL;
  self->mu       = GSL_NAN;
  self->sigma    = GSL_NAN;
  self->tau      = GSL_NAN;
  self->lnl_min  = GSL_NAN;
  self->lnl_max  = GSL_NAN;
  self->prepared = FALSE;

  nc_cluster_richness_projection_set_reltol (self->crp, self->reltol);
}

static void
_nc_cluster_mass_projection_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterMassProjection *mp = NC_CLUSTER_MASS_PROJECTION (object);

  g_return_if_fail (NC_IS_CLUSTER_MASS_PROJECTION (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      nc_cluster_mass_projection_set_reltol (mp, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_cluster_mass_projection_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassProjection *mp                 = NC_CLUSTER_MASS_PROJECTION (object);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (mp);

  g_return_if_fail (NC_IS_CLUSTER_MASS_PROJECTION (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_cluster_mass_projection_dispose (GObject *object)
{
  NcClusterMassProjection *mp                 = NC_CLUSTER_MASS_PROJECTION (object);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (mp);

  nc_cluster_richness_projection_clear (&self->crp);
  self->prepared = FALSE;

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_projection_parent_class)->dispose (object);
}

static void
_nc_cluster_mass_projection_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_projection_parent_class)->finalize (object);
}

static gdouble _nc_cluster_mass_projection_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_projection_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z);

static gdouble _nc_cluster_mass_projection_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
static gdouble _nc_cluster_mass_projection_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_projection_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
static gboolean _nc_cluster_mass_projection_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_projection_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, NcmVector *res);

static void
nc_cluster_mass_projection_class_init (NcClusterMassProjectionClass *klass)
{
  GObjectClass *object_class           = G_OBJECT_CLASS (klass);
  NcClusterMassRichnessClass *mr_class = NC_CLUSTER_MASS_RICHNESS_CLASS (klass);
  NcClusterMassClass *cm_class         = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class           = NCM_MODEL_CLASS (klass);

  /* These go on model_class, not object_class: ncm_model_class_add_params()
   * installs its own GObject handlers and dispatches here for the properties
   * this class owns. */
  model_class->set_property = &_nc_cluster_mass_projection_set_property;
  model_class->get_property = &_nc_cluster_mass_projection_get_property;

  object_class->dispose  = &_nc_cluster_mass_projection_dispose;
  object_class->finalize = &_nc_cluster_mass_projection_finalize;

  ncm_model_class_set_name_nick (model_class, "Ln-normal richness distribution with projection effects", "Projection");
  ncm_model_class_add_params (model_class, NC_CLUSTER_MASS_PROJECTION_SPARAM_LEN - NC_CLUSTER_MASS_RICHNESS_SPARAM_LEN, 0, PROP_SIZE);

  /**
   * NcClusterMassProjection:MU_P0:
   *
   * Constant term (bias) in the mean log-richness relation.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_MU_P0, "mu_p0", "mup0",
                              0.0, 6.0, 1.0e-1,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:MU_P1:
   *
   * Linear mass coefficient in the mean log-richness.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_MU_P1, "mu_p1", "mup1",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:MU_P2:
   *
   * Redshift evolution coefficient in the mean log-richness.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_MU_P2, "mu_p2", "mup2",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:sigma_p0:
   *
   * Constant term (bias) in the standard deviation of the log-richness.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_SIGMA_P0, "\\sigma_p0", "sigmap0",
                              1.0e-4, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:sigma_p1:
   *
   * Linear mass coefficient in the standard deviation.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_SIGMA_P1, "\\sigma_p1", "sigmap1",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:sigma_p2:
   *
   * Redshift evolution coefficient in the standard deviation.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_SIGMA_P2, "\\sigma_p2", "sigmap2",
                              -10.0, 10.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P2,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:f_prj:
   *
   * Fraction of clusters whose richness is inflated by projection. Zero recovers
   * the pure log-normal relation of #NcClusterMassAscaso.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_F_PRJ, "f_\\mathrm{prj}", "fprj",
                              0.0, 1.0, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_F_PRJ,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:tau:
   *
   * Rate of the exponential richness boost; the mean added richness is $1/\tau$.
   * Larger $\tau$ means weaker projection.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PROJECTION_TAU, "\\tau", "tau",
                              1.0e-4, 1.0e2, 1.0e-2,
                              NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PROJECTION_DEFAULT_TAU,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassProjection:reltol:
   *
   * Relative tolerance of the ODE that evaluates the projection term.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance of the projection ODE",
                                                        GSL_DBL_EPSILON, 1.0, NC_CLUSTER_MASS_PROJECTION_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_check_params_info (model_class);

  mr_class->mu    = &_nc_cluster_mass_projection_mu;
  mr_class->sigma = &_nc_cluster_mass_projection_sigma;

  cm_class->P              = &_nc_cluster_mass_projection_p;
  cm_class->intP           = &_nc_cluster_mass_projection_intp;
  cm_class->intP_bin       = &_nc_cluster_mass_projection_intp_bin;
  cm_class->resample       = &_nc_cluster_mass_projection_resample;
  cm_class->P_vec_z_lnMobs = &_nc_cluster_mass_projection_p_vec_z_lnMobs;
}

static gdouble
_nc_cluster_mass_projection_mu (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  NcClusterMassProjection *mp = NC_CLUSTER_MASS_PROJECTION (mr);
  const gdouble lnM0          = nc_cluster_mass_richness_lnM0 (mr);
  const gdouble ln1pz0        = nc_cluster_mass_richness_ln1pz0 (mr);
  const gdouble DlnM          = lnM - lnM0;
  const gdouble Dln1pz        = log1p (z) - ln1pz0;

  return MU_P0 + MU_P1 * DlnM + MU_P2 * Dln1pz;
}

static gdouble
_nc_cluster_mass_projection_sigma (NcClusterMassRichness *mr, gdouble lnM, gdouble z)
{
  NcClusterMassProjection *mp = NC_CLUSTER_MASS_PROJECTION (mr);
  const gdouble lnM0          = nc_cluster_mass_richness_lnM0 (mr);
  const gdouble ln1pz0        = nc_cluster_mass_richness_ln1pz0 (mr);
  const gdouble DlnM          = lnM - lnM0;
  const gdouble Dln1pz        = log1p (z) - ln1pz0;
  const gdouble sigma         = SIGMA_P0 + SIGMA_P1 * DlnM + SIGMA_P2 * Dln1pz;

  /* Add a small number to the standard deviation to avoid numerical instabilities */
  return hypot (sigma, 1.0e-5);
}

/*
 * Solves the projection ODE for the current (mu, sigma, tau), reusing the previous
 * solution whenever the three are unchanged. The range follows CUT and lnR_max,
 * which are model state and may move between likelihood evaluations.
 */
static void
_nc_cluster_mass_projection_prepare (NcClusterMassProjection *mp, const gdouble mu, const gdouble sigma)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (mp);
  const gdouble tau                           = TAU;
  gdouble lnl_min                             = CUT;
  gdouble lnl_max                             = 0.0;

  g_object_get (mp, "lnRichness-max", &lnl_max, NULL);
  g_assert_cmpfloat (lnl_min, <, lnl_max);

  if ((lnl_min != self->lnl_min) || (lnl_max != self->lnl_max))
  {
    nc_cluster_richness_projection_set_lnlambda_range (self->crp, lnl_min, lnl_max);
    self->lnl_min  = lnl_min;
    self->lnl_max  = lnl_max;
    self->prepared = FALSE;
  }

  if (!self->prepared || (mu != self->mu) || (sigma != self->sigma) || (tau != self->tau))
  {
    nc_cluster_richness_projection_prepare (self->crp, mu, sigma, tau);
    self->mu       = mu;
    self->sigma    = sigma;
    self->tau      = tau;
    self->prepared = TRUE;
  }
}

static gdouble
_nc_cluster_mass_projection_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{
  NcClusterMassProjection *mp                 = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (mp);
  NcClusterMassRichness *mr                   = NC_CLUSTER_MASS_RICHNESS (clusterm);
  gdouble mu, sigma;

  if (lnM_obs[0] < CUT)
    return 0.0;

  nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &mu, &sigma);
  {
    const gdouble sigma_t = hypot (sigma, lnM_obs_params[0]);
    const gdouble x       = (lnM_obs[0] - mu) / sigma_t;
    const gdouble clean   = exp (-0.5 * x * x) / (ncm_c_sqrt_2pi () * sigma_t);
    const gdouble f_prj   = F_PRJ;

    if (f_prj <= 0.0)
      return clean;

    _nc_cluster_mass_projection_prepare (mp, mu, sigma_t);

    return (1.0 - f_prj) * clean
           + f_prj * nc_cluster_richness_projection_eval_lnlambda (self->crp, lnM_obs[0]);
  }
}

/*
 * Integral of the log-normal over [lo, hi] in ln(lambda), taking the erf argument
 * on the side where it does not saturate.
 */
static gdouble
_nc_cluster_mass_projection_clean_int (const gdouble mu, const gdouble sigma, const gdouble lo, const gdouble hi)
{
  const gdouble x_lo = (mu - lo) / (M_SQRT2 * sigma);
  const gdouble x_hi = (mu - hi) / (M_SQRT2 * sigma);

  if ((fabs (x_hi) > 4.0) || (fabs (x_lo) > 4.0))
    return -(erfc (x_lo) - erfc (x_hi)) / 2.0;

  return (erf (x_lo) - erf (x_hi)) / 2.0;
}

static gdouble
_nc_cluster_mass_projection_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcClusterMassProjection *mp                 = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (mp);
  NcClusterMassRichness *mr                   = NC_CLUSTER_MASS_RICHNESS (clusterm);
  const gdouble f_prj                         = F_PRJ;
  const gdouble lo                            = CUT;
  gdouble mu, sigma, hi, clean;

  g_object_get (mp, "lnRichness-max", &hi, NULL);

  /* Note: intP does not use the catalog sigma - it represents the selection function */
  nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &mu, &sigma);
  clean = _nc_cluster_mass_projection_clean_int (mu, sigma, lo, hi);

  if (f_prj <= 0.0)
    return clean;

  _nc_cluster_mass_projection_prepare (mp, mu, sigma);

  return (1.0 - f_prj) * clean
         + f_prj * nc_cluster_richness_projection_eval_int (self->crp, lo, hi);
}

static gdouble
_nc_cluster_mass_projection_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params)
{
  NcClusterMassProjection *mp                 = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (mp);
  NcClusterMassRichness *mr                   = NC_CLUSTER_MASS_RICHNESS (clusterm);
  const gdouble cut                           = CUT;
  const gdouble f_prj                         = F_PRJ;
  const gdouble lo                            = GSL_MAX (lnM_obs_lower[0], cut);
  const gdouble hi                            = lnM_obs_upper[0];
  gdouble mu, sigma, sigma_t, clean;

  if (hi <= lo)
    return 0.0;

  nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &mu, &sigma);
  sigma_t = hypot (sigma, lnM_obs_params[0]);
  clean   = _nc_cluster_mass_projection_clean_int (mu, sigma_t, lo, hi);

  if (f_prj <= 0.0)
    return GSL_MAX (clean, 0.0);

  _nc_cluster_mass_projection_prepare (mp, mu, sigma_t);

  {
    const gdouble bin = (1.0 - f_prj) * clean
                        + f_prj * nc_cluster_richness_projection_eval_int (self->crp, lo, hi);

    return GSL_MAX (bin, 0.0);
  }
}

static gboolean
_nc_cluster_mass_projection_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
  NcClusterMassProjection *mp = NC_CLUSTER_MASS_PROJECTION (clusterm);
  NcClusterMassRichness *mr   = NC_CLUSTER_MASS_RICHNESS (clusterm);
  const gdouble cut           = CUT;
  const gdouble f_prj         = F_PRJ;
  const gdouble tau           = TAU;
  gdouble mu, sigma, sigma_t, lnl, lnl_max;

  g_object_get (mp, "lnRichness-max", &lnl_max, NULL);

  nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &mu, &sigma);
  sigma_t = hypot (sigma, lnM_obs_params[0]);

  ncm_rng_lock (rng);

  /* The mixture is sampled directly: no ODE is involved. */
  lnl = ncm_rng_gaussian_gen (rng, mu, sigma_t);

  if ((f_prj > 0.0) && (ncm_rng_uniform01_gen (rng) < f_prj))
    lnl = log (exp (lnl) + ncm_rng_exponential_gen (rng, 1.0 / tau));

  ncm_rng_unlock (rng);

  lnM_obs[0] = lnl;

  return (lnl <= lnl_max) && (lnl >= cut);
}

static void
_nc_cluster_mass_projection_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, NcmVector *res)
{
  const guint len = ncm_vector_len (z);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble z_i          = ncm_vector_get (z, i);
    const gdouble lnM_obs_i    = ncm_matrix_get (lnM_obs, i, 0);
    const gdouble lnM_params_i = ncm_matrix_get (lnM_obs_params, i, 0);

    ncm_vector_set (res, i, _nc_cluster_mass_projection_p (clusterm, cosmo, lnM, z_i, &lnM_obs_i, &lnM_params_i));
  }
}

/**
 * nc_cluster_mass_projection_set_reltol:
 * @mp: a #NcClusterMassProjection
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance of the ODE that evaluates the projection term. The
 * cost of each evaluation falls steeply with a looser tolerance; the default
 * $10^{-7}$ already keeps the relative error of the term near $10^{-5}$.
 *
 */
void
nc_cluster_mass_projection_set_reltol (NcClusterMassProjection *mp, gdouble reltol)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (mp);

  self->reltol   = reltol;
  self->prepared = FALSE;

  nc_cluster_richness_projection_set_reltol (self->crp, reltol);
}

/**
 * nc_cluster_mass_projection_get_reltol:
 * @mp: a #NcClusterMassProjection
 *
 * Returns: the relative tolerance of the projection ODE.
 */
gdouble
nc_cluster_mass_projection_get_reltol (NcClusterMassProjection *mp)
{
  NcClusterMassProjectionPrivate * const self = nc_cluster_mass_projection_get_instance_private (mp);

  return self->reltol;
}

/**
 * nc_cluster_mass_projection_f_prj:
 * @mp: a #NcClusterMassProjection
 *
 * Returns: the fraction $f_\mathrm{prj}$ of clusters affected by projection.
 */
gdouble
nc_cluster_mass_projection_f_prj (NcClusterMassProjection *mp)
{
  return F_PRJ;
}

/**
 * nc_cluster_mass_projection_tau:
 * @mp: a #NcClusterMassProjection
 *
 * Returns: the rate $\tau$ of the exponential richness boost.
 */
gdouble
nc_cluster_mass_projection_tau (NcClusterMassProjection *mp)
{
  return TAU;
}
