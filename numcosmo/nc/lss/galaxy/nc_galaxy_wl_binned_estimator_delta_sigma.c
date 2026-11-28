/***************************************************************************
 *            nc_galaxy_wl_binned_estimator_delta_sigma.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_estimator_delta_sigma.c
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * NcGalaxyWLBinnedEstimatorDeltaSigma:
 *
 * Binned excess surface mass density.
 *
 * The binned excess surface mass density as in Umetsu et al. 2020 (ApJ 890, 148, Sections
 * 3.2-3.3 and 4.1), whose estimator is the one of Mandelbaum et al. 2018a (PASJ 70, S25).
 * Equation numbers refer to Umetsu et al. 2020; the original source of each step is cited
 * where it is used.
 *
 * Each galaxy $s$ of bin $i$ has the mean inverse critical density over its redshift
 * density, $\langle\Sigma^{-1}_{\mathrm{cr},s}\rangle$ (Eq. 6, zero in front of the cluster),
 * and the weight (Eq. 7)
 * $$w_s = \frac{\langle\Sigma^{-1}_{\mathrm{cr},s}\rangle^2}{\sigma_{e,s}^2 + e_{\mathrm{rms},s}^2},$$
 * with $\sigma_{e,s}$ the column "std_noise" and $e_{\mathrm{rms},s}$ the column "e_rms"
 * (or #NcGalaxyWLBinnedEstimator:e-rms). The bin has the multiplicative bias correction
 * $1 + K = \sum w_s (1 + m_s) / \sum w_s$ (Eq. 8; Mandelbaum et al. 2018a, 2018b), the
 * responsivity $\mathcal{R} = 1 - \sum w_s e_{\mathrm{rms},s}^2 / \sum w_s$ (Eq. 9;
 * Mandelbaum et al. 2005, 2018a) and the effective radius
 * $R_i = \sum w_s / \sum w_s R_s^{-1}$, the weighted harmonic mean of the radii of its
 * galaxies (Eq. 11; Okabe & Smith 2016, Sereno et al. 2017). The estimator (Eq. 4;
 * Mandelbaum et al. 2018a), with the subtraction of the weighted mean additive bias
 * (Section 3.2; Mandelbaum et al. 2018a, Miyaoka et al. 2018, Okabe et al. 2019), is
 * $$\Delta\Sigma_+(R_i) = \frac{1}{1 + K}\left[\frac{\sum w_s e_{+,s}
 *   \langle\Sigma^{-1}_{\mathrm{cr},s}\rangle^{-1}}{2\mathcal{R}\sum w_s}
 *   - \frac{\sum w_s c_{+,s} \langle\Sigma^{-1}_{\mathrm{cr},s}\rangle^{-1}}{\sum w_s}\right],$$
 * where the second term is dropped when #NcGalaxyWLBinnedEstimatorDeltaSigma:subtract-additive-bias
 * is %FALSE. Its variance is the shape noise
 * $\sigma^2(R_i) = [4\mathcal{R}^2 (1 + K)^2 \sum w_s]^{-1}$ (Eq. 14; Miyaoka et al. 2018),
 * and the covariance is diagonal: the large-scale structure (Hoekstra 2003) and intrinsic
 * (Gruen et al. 2015) terms of Eq. 16 are not included. For a
 * catalog in the ellipticity convention (%NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET), $2\mathcal{R}$
 * is replaced by one.
 *
 * The model, with the reduced-shear correction of Umetsu et al. 2014 (Eq. 22), is
 * $$f(R_i) = \frac{\Delta\Sigma(R_i)}{1 - \langle\langle\Sigma^{-1}_{\mathrm{cr}}\rangle\rangle_i
 *   \Sigma(R_i)},\qquad \langle\langle\Sigma^{-1}_{\mathrm{cr}}\rangle\rangle_i =
 *   \frac{\sum w_s \langle\Sigma^{-1}_{\mathrm{cr},s}\rangle}{\sum w_s}\quad\text{(Eq. 23)},$$
 * with $\Sigma$ and $\Delta\Sigma$ of the #NcHaloDensityProfile of the #NcmMSet. The
 * photo-z bias correction of Eq. 18 and the stacking over lenses are not part of this
 * estimator.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_wl_binned_estimator_delta_sigma.h"
#include "nc/lss/galaxy/nc_galaxy_wl_binned_estimator_private.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "nc/background/nc_hicosmo.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyWLBinnedEstimatorDeltaSigma
{
  /*< private >*/
  NcGalaxyWLBinnedEstimator parent_instance;
  gboolean subtract_additive_bias;
};

enum
{
  PROP_0,
  PROP_SUBTRACT_ADDITIVE_BIAS,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyWLBinnedEstimatorDeltaSigma, nc_galaxy_wl_binned_estimator_delta_sigma, NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR)

static void
nc_galaxy_wl_binned_estimator_delta_sigma_init (NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds)
{
  est_ds->subtract_additive_bias = TRUE;
}

static void
_nc_galaxy_wl_binned_estimator_delta_sigma_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds = NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA (object);

  switch (prop_id)
  {
    case PROP_SUBTRACT_ADDITIVE_BIAS:
      nc_galaxy_wl_binned_estimator_delta_sigma_set_subtract_additive_bias (est_ds, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_binned_estimator_delta_sigma_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds = NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA (object);

  switch (prop_id)
  {
    case PROP_SUBTRACT_ADDITIVE_BIAS:
      g_value_set_boolean (value, nc_galaxy_wl_binned_estimator_delta_sigma_get_subtract_additive_bias (est_ds));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static GStrv
_nc_galaxy_wl_binned_estimator_delta_sigma_extra_columns (NcGalaxyWLBinnedEstimator *est)
{
  GStrv cols = g_new0 (gchar *, 4);

  cols[0] = g_strdup (NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_INV_SIGMA_CRIT);
  cols[1] = g_strdup (NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_RESPONSIVITY);
  cols[2] = g_strdup (NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_MULT_BIAS);

  return cols;
}

typedef struct _NcWLBEDeltaSigmaBin
{
  gdouble w;      /* sum w */
  gdouble w_m;    /* sum w (1 + m) */
  gdouble w_erms; /* sum w e_rms^2 */
  gdouble w_invr; /* sum w / R */
  gdouble w_s;    /* sum w <Sigma_cr^-1> */
  gdouble w_et;   /* sum w e_t / <Sigma_cr^-1> */
  gdouble w_ct;   /* sum w c_t / <Sigma_cr^-1> */
} NcWLBEDeltaSigmaBin;

static gdouble
_nc_wlbe_ds_weight (const NcGalaxyWLBinnedEstimatorSources *src, const guint i)
{
  const gdouble s_i = NC_WLBE_SRC (src, inv_sigma_crit, i);

  if (s_i > 0.0)
  {
    const gdouble std_noise = NC_WLBE_SRC (src, std_noise, i);
    const gdouble e_rms     = NC_WLBE_SRC (src, e_rms, i);

    return s_i * s_i / (std_noise * std_noise + e_rms * e_rms);
  }

  return 0.0;
}

static void
_nc_galaxy_wl_binned_estimator_delta_sigma_build_bins (NcGalaxyWLBinnedEstimator *est, const NcGalaxyWLBinnedEstimatorSources *src, NcGalaxyWLObsBinned *binned, GArray *factor)
{
  NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds = NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA (est);
  NcWLBEDeltaSigmaBin *acc                    = g_new0 (NcWLBEDeltaSigmaBin, src->n_bins);
  NcmMatrix *cov                              = nc_galaxy_wl_obs_binned_peek_cov (binned);
  guint i, b;

  for (i = 0; i < src->len; i++)
  {
    const gdouble w_i = _nc_wlbe_ds_weight (src, i);
    NcWLBEDeltaSigmaBin *a;

    if (w_i == 0.0)
      continue;

    a = &acc[NC_WLBE_SRC_BIN (src, i)];

    {
      const gdouble s_i   = NC_WLBE_SRC (src, inv_sigma_crit, i);
      const gdouble e_rms = NC_WLBE_SRC (src, e_rms, i);

      a->w      += w_i;
      a->w_m    += w_i * (1.0 + NC_WLBE_SRC (src, m, i));
      a->w_erms += w_i * e_rms * e_rms;
      a->w_invr += w_i / NC_WLBE_SRC (src, radius, i);
      a->w_s    += w_i * s_i;
      a->w_et   += w_i * NC_WLBE_SRC (src, e_t, i) / s_i;
      a->w_ct   += w_i * NC_WLBE_SRC (src, c_t, i) / s_i;
    }
  }

  ncm_matrix_set_zero (cov);

  for (b = 0; b < src->n_bins; b++)
  {
    const NcWLBEDeltaSigmaBin *a = &acc[b];
    const gdouble mult           = a->w_m / a->w;
    const gdouble resp           = src->use_responsivity ? 1.0 - a->w_erms / a->w : 1.0;
    const gdouble scale          = src->shear_scale * resp;
    gdouble delta_sigma          = a->w_et / (scale * a->w);
    gdouble offset               = 0.0;

    if (est_ds->subtract_additive_bias)
      delta_sigma -= a->w_ct / a->w;
    else
      offset = a->w_ct / a->w / mult;

    delta_sigma /= mult;

    /* A bin with zero weight gets NaN entries and is dropped by the base class */
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_RADIUS, b, a->w / a->w_invr);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_VALUE, b, delta_sigma);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_OFFSET, b, offset);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_INV_SIGMA_CRIT, b, a->w_s / a->w);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_RESPONSIVITY, b, resp);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_MULT_BIAS, b, mult);

    ncm_matrix_set (cov, b, b, 1.0 / (scale * scale * mult * mult * a->w));
  }

  /*
   * Expected value: <e_t> / (shear_scale R) = (1 + m) g + c_t, so the estimator has
   * expectation sum w (1 + m) g / <Sigma_cr^-1> / (sum w (1 + K)) plus the offset.
   */
  for (i = 0; i < src->len; i++)
  {
    const gdouble w_i = _nc_wlbe_ds_weight (src, i);

    if (w_i == 0.0)
    {
      g_array_index (factor, gdouble, i) = 0.0;
    }
    else
    {
      const NcWLBEDeltaSigmaBin *a = &acc[NC_WLBE_SRC_BIN (src, i)];

      g_array_index (factor, gdouble, i) = w_i * (1.0 + NC_WLBE_SRC (src, m, i)) / (NC_WLBE_SRC (src, inv_sigma_crit, i) * a->w_m);
    }
  }

  g_free (acc);
}

static void
_nc_galaxy_wl_binned_estimator_delta_sigma_eval_model (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObsBinned *binned, NcmVector *model)
{
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *hp          = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  const guint n_bins          = nc_galaxy_wl_obs_binned_len (binned);
  gdouble z_cl;
  guint b;

  if (!nc_galaxy_wl_obs_binned_has_column (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_INV_SIGMA_CRIT))
    g_error ("nc_galaxy_wl_binned_estimator_delta_sigma: the profile has no column `%s'.",
             NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_INV_SIGMA_CRIT);

  nc_halo_position_prepare_if_needed (hp, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (smd, cosmo);
  z_cl = nc_halo_position_get_redshift (hp);

  for (b = 0; b < n_bins; b++)
  {
    const gdouble radius   = nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_OBS_BINNED_COL_RADIUS, b);
    const gdouble inv_crit = nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_INV_SIGMA_CRIT, b);
    const gdouble sigma    = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, radius, z_cl);
    const gdouble d_sigma  = nc_wl_surface_mass_density_sigma_excess (smd, dp, cosmo, radius, z_cl);

    ncm_vector_set (model, b, d_sigma / (1.0 - inv_crit * sigma));
  }
}

static void
nc_galaxy_wl_binned_estimator_delta_sigma_class_init (NcGalaxyWLBinnedEstimatorDeltaSigmaClass *klass)
{
  GObjectClass *object_class                = G_OBJECT_CLASS (klass);
  NcGalaxyWLBinnedEstimatorClass *est_class = NC_GALAXY_WL_BINNED_ESTIMATOR_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_binned_estimator_delta_sigma_set_property;
  object_class->get_property = &_nc_galaxy_wl_binned_estimator_delta_sigma_get_property;

  /**
   * NcGalaxyWLBinnedEstimatorDeltaSigma:subtract-additive-bias:
   *
   * Whether to subtract the weighted mean additive bias of each bin (Umetsu et al. 2020,
   * Section 3.2).
   */
  g_object_class_install_property (object_class,
                                   PROP_SUBTRACT_ADDITIVE_BIAS,
                                   g_param_spec_boolean ("subtract-additive-bias",
                                                         NULL,
                                                         "Whether to subtract the additive bias",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  est_class->name          = "delta_sigma";
  est_class->extra_columns = &_nc_galaxy_wl_binned_estimator_delta_sigma_extra_columns;
  est_class->build_bins    = &_nc_galaxy_wl_binned_estimator_delta_sigma_build_bins;
  est_class->eval_model    = &_nc_galaxy_wl_binned_estimator_delta_sigma_eval_model;
}

/**
 * nc_galaxy_wl_binned_estimator_delta_sigma_new:
 *
 * Creates a new #NcGalaxyWLBinnedEstimatorDeltaSigma.
 *
 * Returns: (transfer full): a new #NcGalaxyWLBinnedEstimatorDeltaSigma.
 */
NcGalaxyWLBinnedEstimatorDeltaSigma *
nc_galaxy_wl_binned_estimator_delta_sigma_new (void)
{
  return g_object_new (NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA, NULL);
}

/**
 * nc_galaxy_wl_binned_estimator_delta_sigma_ref:
 * @est_ds: a #NcGalaxyWLBinnedEstimatorDeltaSigma
 *
 * Increases the reference count of @est_ds by one.
 *
 * Returns: (transfer full): @est_ds.
 */
NcGalaxyWLBinnedEstimatorDeltaSigma *
nc_galaxy_wl_binned_estimator_delta_sigma_ref (NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds)
{
  return g_object_ref (est_ds);
}

/**
 * nc_galaxy_wl_binned_estimator_delta_sigma_free:
 * @est_ds: a #NcGalaxyWLBinnedEstimatorDeltaSigma
 *
 * Decreases the reference count of @est_ds by one.
 */
void
nc_galaxy_wl_binned_estimator_delta_sigma_free (NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds)
{
  g_object_unref (est_ds);
}

/**
 * nc_galaxy_wl_binned_estimator_delta_sigma_clear:
 * @est_ds: a #NcGalaxyWLBinnedEstimatorDeltaSigma
 *
 * Decreases the reference count of *@est_ds by one and sets *@est_ds to %NULL.
 */
void
nc_galaxy_wl_binned_estimator_delta_sigma_clear (NcGalaxyWLBinnedEstimatorDeltaSigma **est_ds)
{
  g_clear_object (est_ds);
}

/**
 * nc_galaxy_wl_binned_estimator_delta_sigma_set_subtract_additive_bias:
 * @est_ds: a #NcGalaxyWLBinnedEstimatorDeltaSigma
 * @subtract: whether to subtract the additive bias
 *
 * Sets #NcGalaxyWLBinnedEstimatorDeltaSigma:subtract-additive-bias.
 */
void
nc_galaxy_wl_binned_estimator_delta_sigma_set_subtract_additive_bias (NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds, gboolean subtract)
{
  est_ds->subtract_additive_bias = subtract;
}

/**
 * nc_galaxy_wl_binned_estimator_delta_sigma_get_subtract_additive_bias:
 * @est_ds: a #NcGalaxyWLBinnedEstimatorDeltaSigma
 *
 * Returns: #NcGalaxyWLBinnedEstimatorDeltaSigma:subtract-additive-bias.
 */
gboolean
nc_galaxy_wl_binned_estimator_delta_sigma_get_subtract_additive_bias (NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds)
{
  return est_ds->subtract_additive_bias;
}

