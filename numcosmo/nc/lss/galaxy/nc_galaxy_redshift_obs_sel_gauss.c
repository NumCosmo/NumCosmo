/***************************************************************************
 *            nc_galaxy_redshift_obs_sel_gauss.c
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_obs_sel_gauss.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyRedshiftObsSelGauss:
 *
 * Gaussian population photo-z observable distribution.
 *
 * A #NcGalaxyRedshiftObsSel whose population observable at true
 * redshift $z$ is Gaussian, $\mathrm{obs} \sim \mathcal{N}(z, \sigma_z)$ with
 * $\sigma_z = \sigma_0 (1 + z)$ and a single population-scatter parameter
 * $\sigma_0$. Currently shares the Gaussian-integral math with the per-galaxy
 * #NcGalaxyRedshiftObsGauss, but is a distinct abstraction free to
 * become a scatter mixture as its assumptions diverge.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_obs_sel_gauss.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyRedshiftObsSelGauss
{
  NcGalaxyRedshiftObsSel parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyRedshiftObsSelGauss, nc_galaxy_redshift_obs_sel_gauss, NC_TYPE_GALAXY_REDSHIFT_OBS_SEL);

static void
nc_galaxy_redshift_obs_sel_gauss_init (NcGalaxyRedshiftObsSelGauss *gsdropg)
{
}

static void
_nc_galaxy_redshift_obs_sel_gauss_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_obs_sel_gauss_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_redshift_obs_sel_gauss_eval (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs);
static gdouble _nc_galaxy_redshift_obs_sel_gauss_window_mass (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs_lo, const gdouble obs_hi);

static void
nc_galaxy_redshift_obs_sel_gauss_class_init (NcGalaxyRedshiftObsSelGaussClass *klass)
{
  NcGalaxyRedshiftObsSelClass *gsdrop_class = NC_GALAXY_REDSHIFT_OBS_SEL_CLASS (klass);
  GObjectClass *object_class                = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_redshift_obs_sel_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian population photo-z observable", "GaussRedshiftObsSel");
  ncm_model_class_add_params (model_class, NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxyRedshiftObsSelGauss:sigma0:
   *
   * The population photo-z scatter, sigma_z = sigma0 (1 + z).
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_SIGMA0, "\\sigma_0", "sigma0",
                              1.0e-8, 1.0, 1.0e-2,
                              NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_DEFAULT_SIGMA0,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  gsdrop_class->eval        = &_nc_galaxy_redshift_obs_sel_gauss_eval;
  gsdrop_class->window_mass = &_nc_galaxy_redshift_obs_sel_gauss_window_mass;
}

#define VECTOR (NCM_MODEL (gsdrop))
#define SIGMA0 (ncm_model_orig_param_get (VECTOR, NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_SIGMA0))

static gdouble
_nc_galaxy_redshift_obs_sel_gauss_eval (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs)
{
  const gdouble sigmaz         = SIGMA0 * (1.0 + z);
  const gdouble sqrt2pi_sigmaz = M_SQRT2 * M_SQRTPI * sigmaz;
  const gdouble gauss          = exp (-0.5 * gsl_pow_2 ((obs - z) / sigmaz));

  /* N(obs | z) density for obs ~ N(z, sigmaz). */
  return gauss / sqrt2pi_sigmaz;
}

static gdouble
_nc_galaxy_redshift_obs_sel_gauss_window_mass (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs_lo, const gdouble obs_hi)
{
  const gdouble sigmaz = SIGMA0 * (1.0 + z);

  /* P(obs in [obs_lo, obs_hi] | z) for obs ~ N(z, sigmaz). */
  return ncm_util_gaussian_integral (obs_lo, obs_hi, z, sigmaz);
}

/**
 * nc_galaxy_redshift_obs_sel_gauss_new:
 *
 * Creates a new #NcGalaxyRedshiftObsSelGauss.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftObsSelGauss.
 */
NcGalaxyRedshiftObsSelGauss *
nc_galaxy_redshift_obs_sel_gauss_new (void)
{
  NcGalaxyRedshiftObsSelGauss *gsdropg = g_object_new (NC_TYPE_GALAXY_REDSHIFT_OBS_SEL_GAUSS,
                                                       NULL);

  return gsdropg;
}

/**
 * nc_galaxy_redshift_obs_sel_gauss_ref:
 * @gsdropg: a #NcGalaxyRedshiftObsSelGauss
 *
 * Increases the reference count of @gsdropg by one.
 *
 * Returns: (transfer full): @gsdropg.
 */
NcGalaxyRedshiftObsSelGauss *
nc_galaxy_redshift_obs_sel_gauss_ref (NcGalaxyRedshiftObsSelGauss *gsdropg)
{
  return g_object_ref (gsdropg);
}

/**
 * nc_galaxy_redshift_obs_sel_gauss_free:
 * @gsdropg: a #NcGalaxyRedshiftObsSelGauss
 *
 * Decreases the reference count of @gsdropg by one.
 *
 */
void
nc_galaxy_redshift_obs_sel_gauss_free (NcGalaxyRedshiftObsSelGauss *gsdropg)
{
  g_object_unref (gsdropg);
}

/**
 * nc_galaxy_redshift_obs_sel_gauss_clear:
 * @gsdropg: a #NcGalaxyRedshiftObsSelGauss
 *
 * Decreases the reference count of @gsdropg by one, and sets the pointer
 * *@gsdropg to NULL.
 *
 */
void
nc_galaxy_redshift_obs_sel_gauss_clear (NcGalaxyRedshiftObsSelGauss **gsdropg)
{
  g_clear_object (gsdropg);
}

