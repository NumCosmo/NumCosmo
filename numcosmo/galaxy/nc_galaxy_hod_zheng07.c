/***************************************************************************
 *            nc_galaxy_hod_zheng07.c
 *
 *  Sun Jun 14 12:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_hod_zheng07.c
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyHODZheng07:
 *
 * Zheng et al. (2007) halo occupation distribution.
 *
 * Concrete #NcGalaxyHOD with the five-parameter occupation of Zheng et al.
 * (2007). The mean central occupation is
 * $$\langle N_\mathrm{cen} \rangle = \frac{1}{2}\left[1 + \mathrm{erf}\left(
 * \frac{\log_{10} M - \log_{10} M_\mathrm{min}}{\sigma_{\log M}}\right)\right],$$
 * and the mean satellite occupation is
 * $$\langle N_\mathrm{sat} \rangle = \left(\frac{\max(0, M - M_0)}{M_1}
 * \right)^{\alpha},$$
 * with $M_0 = 10^{\log_{10} M_0}$ and $M_1 = 10^{\log_{10} M_1}$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "galaxy/nc_galaxy_hod_zheng07.h"
#include "galaxy/nc_galaxy_hod.h"
#include "math/ncm_model.h"
#include "math/ncm_util.h"

#include <math.h>

struct _NcGalaxyHODZheng07
{
  NcGalaxyHOD parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyHODZheng07, nc_galaxy_hod_zheng07, NC_TYPE_GALAXY_HOD)

static void
nc_galaxy_hod_zheng07_init (NcGalaxyHODZheng07 *zheng07)
{
}

static void
_nc_galaxy_hod_zheng07_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_hod_zheng07_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_hod_zheng07_mean_n_central (NcGalaxyHOD *hod, const gdouble lnM);
static gdouble _nc_galaxy_hod_zheng07_mean_n_satellite (NcGalaxyHOD *hod, const gdouble lnM);

static void
nc_galaxy_hod_zheng07_class_init (NcGalaxyHODZheng07Class *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcGalaxyHODClass *hod_class = NC_GALAXY_HOD_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_hod_zheng07_finalize;

  ncm_model_class_set_name_nick (model_class, "Zheng07 HOD", "Zheng07");
  ncm_model_class_add_params (model_class, NC_GALAXY_HOD_ZHENG07_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxyHODZheng07:logMmin:
   *
   * The $\log_{10}$ of the central cutoff mass (solar masses).
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_HOD_ZHENG07_LOG_MMIN, "\\log_{10}M_\\mathrm{min}", "logMmin",
                              8.0, 16.0, 1.0e-2,
                              NC_GALAXY_HOD_ZHENG07_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_HOD_ZHENG07_DEFAULT_LOG_MMIN,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxyHODZheng07:sigmalogM:
   *
   * The width of the central occupation transition.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_HOD_ZHENG07_SIGMA_LOG_M, "\\sigma_{\\log M}", "sigmalogM",
                              1.0e-3, 2.0, 1.0e-2,
                              NC_GALAXY_HOD_ZHENG07_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_HOD_ZHENG07_DEFAULT_SIGMA_LOG_M,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxyHODZheng07:logM0:
   *
   * The $\log_{10}$ of the satellite cutoff mass (solar masses).
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_HOD_ZHENG07_LOG_M0, "\\log_{10}M_0", "logM0",
                              8.0, 16.0, 1.0e-2,
                              NC_GALAXY_HOD_ZHENG07_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_HOD_ZHENG07_DEFAULT_LOG_M0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxyHODZheng07:logM1:
   *
   * The $\log_{10}$ of the satellite normalization mass (solar masses).
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_HOD_ZHENG07_LOG_M1, "\\log_{10}M_1", "logM1",
                              8.0, 16.0, 1.0e-2,
                              NC_GALAXY_HOD_ZHENG07_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_HOD_ZHENG07_DEFAULT_LOG_M1,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxyHODZheng07:alpha:
   *
   * The satellite power-law slope.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_HOD_ZHENG07_ALPHA, "\\alpha", "alpha",
                              0.0, 3.0, 1.0e-2,
                              NC_GALAXY_HOD_ZHENG07_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_HOD_ZHENG07_DEFAULT_ALPHA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  hod_class->mean_n_central   = &_nc_galaxy_hod_zheng07_mean_n_central;
  hod_class->mean_n_satellite = &_nc_galaxy_hod_zheng07_mean_n_satellite;
}

#define VECTOR   (NCM_MODEL (hod))
#define LOG_MMIN (ncm_model_orig_param_get (VECTOR, NC_GALAXY_HOD_ZHENG07_LOG_MMIN))
#define SIGMA_LOG_M (ncm_model_orig_param_get (VECTOR, NC_GALAXY_HOD_ZHENG07_SIGMA_LOG_M))
#define LOG_M0   (ncm_model_orig_param_get (VECTOR, NC_GALAXY_HOD_ZHENG07_LOG_M0))
#define LOG_M1   (ncm_model_orig_param_get (VECTOR, NC_GALAXY_HOD_ZHENG07_LOG_M1))
#define ALPHA    (ncm_model_orig_param_get (VECTOR, NC_GALAXY_HOD_ZHENG07_ALPHA))

static gdouble
_nc_galaxy_hod_zheng07_mean_n_central (NcGalaxyHOD *hod, const gdouble lnM)
{
  const gdouble log10_M = lnM / M_LN10;
  const gdouble arg     = (log10_M - LOG_MMIN) / SIGMA_LOG_M;

  /* 0.5 * (1 + erf(arg)) = 0.5 + N(0, arg * sqrt(2)), the unit-normal CDF. */
  return 0.5 + ncm_util_normal_gaussian_integral (0.0, arg * M_SQRT2);
}

static gdouble
_nc_galaxy_hod_zheng07_mean_n_satellite (NcGalaxyHOD *hod, const gdouble lnM)
{
  const gdouble M    = exp (lnM);
  const gdouble M0   = pow (10.0, LOG_M0);
  const gdouble M1   = pow (10.0, LOG_M1);
  const gdouble diff = M - M0;

  if (diff <= 0.0)
    return 0.0;

  return pow (diff / M1, ALPHA);
}

/**
 * nc_galaxy_hod_zheng07_new:
 *
 * Creates a new #NcGalaxyHODZheng07 with the default Zheng et al. (2007)
 * parameters.
 *
 * Returns: (transfer full): a new #NcGalaxyHODZheng07.
 */
NcGalaxyHODZheng07 *
nc_galaxy_hod_zheng07_new (void)
{
  return g_object_new (NC_TYPE_GALAXY_HOD_ZHENG07, NULL);
}

/**
 * nc_galaxy_hod_zheng07_ref:
 * @zheng07: a #NcGalaxyHODZheng07
 *
 * Increases the reference count of @zheng07 by one.
 *
 * Returns: (transfer full): @zheng07.
 */
NcGalaxyHODZheng07 *
nc_galaxy_hod_zheng07_ref (NcGalaxyHODZheng07 *zheng07)
{
  return g_object_ref (zheng07);
}

/**
 * nc_galaxy_hod_zheng07_free:
 * @zheng07: a #NcGalaxyHODZheng07
 *
 * Decreases the reference count of @zheng07 by one.
 *
 */
void
nc_galaxy_hod_zheng07_free (NcGalaxyHODZheng07 *zheng07)
{
  g_object_unref (zheng07);
}

/**
 * nc_galaxy_hod_zheng07_clear:
 * @zheng07: a #NcGalaxyHODZheng07
 *
 * If *@zheng07 is different from %NULL, decreases its reference count and sets
 * *@zheng07 to %NULL.
 *
 */
void
nc_galaxy_hod_zheng07_clear (NcGalaxyHODZheng07 **zheng07)
{
  g_clear_object (zheng07);
}
