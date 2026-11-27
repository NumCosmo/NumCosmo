/***************************************************************************
 *            nc_galaxy_wl_binned_estimator_gt.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_estimator_gt.c
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
 * NcGalaxyWLBinnedEstimatorGt:
 *
 * Binned reduced tangential shear.
 *
 * The weighted mean reduced tangential shear in radial bins and its model with the mean
 * lensing efficiency of the sources (Schrabback et al. 2018, Eq. 12). The profile, its
 * error and its model agree with the CLMM library (make_radial_profile() with weights and
 * the "order2" reduced tangential shear), which the tests use as a reference.
 *
 * Every galaxy has the shape weight of Sheldon et al. 2004
 * $$w_i = \frac{1}{\sigma_\mathrm{pop}^2 + \sigma_{e_1,i}^2 + \sigma_{e_2,i}^2}
 *   = \frac{1}{\sigma_\mathrm{pop}^2 + 2\sigma_{e,i}^2},$$
 * with $\sigma_\mathrm{pop}$ = #NcGalaxyWLBinnedEstimatorGt:sigma-pop, the shape noise, and
 * $\sigma_{e,i}$ the per-component measurement noise, the column "std_noise". The
 * calibration follows Mandelbaum et al. 2018a (Appendix A), taken over the whole sample:
 * the estimator computes the weighted means of the multiplicative bias $\hat m$ and of the
 * additive bias $\hat c$, and the responsivity
 * $\mathcal{R} = 1 - \sum w_i e_{\mathrm{rms},i}^2 / \sum w_i$. The shear of each galaxy is
 * $$g_i = \frac{1}{1 + \hat m}\left(\frac{e_i}{2\mathcal{R}} - \hat c\right),$$
 * and its tangential component $g_{t,i}$ is averaged with the weights $w_i$ in each bin; the
 * bin radius is the weighted mean radius of its galaxies. For a catalog in the ellipticity convention
 * (%NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET), $2\mathcal{R}$ is replaced by one. The error of
 * a bin is set by #NcGalaxyWLBinnedEstimatorGt:error-type, and the covariance is diagonal.
 *
 * The model (Schrabback et al. 2018, Eq. 12; see also Seitz & Schneider 1997) is
 * $$\langle g_t\rangle = \frac{\langle\beta_s\rangle\gamma_{t,\infty}}{1 -
 *   \langle\beta_s\rangle\kappa_\infty}\left[1 + \left(\frac{\langle\beta_s^2\rangle}
 *   {\langle\beta_s\rangle^2} - 1\right)\langle\beta_s\rangle\kappa_\infty\right],$$
 * with $\gamma_{t,\infty} = \Delta\Sigma / \Sigma_\mathrm{cr,\infty}$ and
 * $\kappa_\infty = \Sigma / \Sigma_\mathrm{cr,\infty}$ at the bin radius, and the lensing
 * efficiency $\beta_s = \Sigma_\mathrm{cr,\infty} / \Sigma_\mathrm{cr}(z_s)$ (zero in front
 * of the cluster). The moments $\langle\beta_s\rangle$ and $\langle\beta_s^2\rangle$
 * are the means over the galaxies of the profile of their averages over their own redshift
 * densities, the same in every bin, that is, the moments of the stacked redshift
 * distribution of the sources (or the means weighted by $w_i$, see
 * #NcGalaxyWLBinnedEstimatorGt:beta-weighted). Sources in front of the cluster count with
 * $\beta_s = 0$. The moments of an external $N(z)$, for instance one cut at $z_l + \Delta z$,
 * can be set directly with nc_galaxy_wl_obs_binned_set_column() on the columns
 * #NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S_MEAN and
 * #NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S2_MEAN.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_wl_binned_estimator_gt.h"
#include "nc/lss/galaxy/nc_galaxy_wl_binned_estimator_private.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "nc/background/nc_hicosmo.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyWLBinnedEstimatorGt
{
  /*< private >*/
  NcGalaxyWLBinnedEstimator parent_instance;
  gdouble sigma_pop;
  NcGalaxyWLBinnedEstimatorGtError error_type;
  gboolean beta_weighted;
};

enum
{
  PROP_0,
  PROP_SIGMA_POP,
  PROP_ERROR_TYPE,
  PROP_BETA_WEIGHTED,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyWLBinnedEstimatorGt, nc_galaxy_wl_binned_estimator_gt, NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR)

static void
nc_galaxy_wl_binned_estimator_gt_init (NcGalaxyWLBinnedEstimatorGt *est_gt)
{
  est_gt->sigma_pop     = NC_GALAXY_WL_BINNED_ESTIMATOR_GT_DEFAULT_SIGMA_POP;
  est_gt->error_type    = NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_STD_ERR;
  est_gt->beta_weighted = FALSE;
}

static void
_nc_galaxy_wl_binned_estimator_gt_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLBinnedEstimatorGt *est_gt = NC_GALAXY_WL_BINNED_ESTIMATOR_GT (object);

  switch (prop_id)
  {
    case PROP_SIGMA_POP:
      nc_galaxy_wl_binned_estimator_gt_set_sigma_pop (est_gt, g_value_get_double (value));
      break;
    case PROP_ERROR_TYPE:
      nc_galaxy_wl_binned_estimator_gt_set_error_type (est_gt, g_value_get_enum (value));
      break;
    case PROP_BETA_WEIGHTED:
      nc_galaxy_wl_binned_estimator_gt_set_beta_weighted (est_gt, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_binned_estimator_gt_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLBinnedEstimatorGt *est_gt = NC_GALAXY_WL_BINNED_ESTIMATOR_GT (object);

  switch (prop_id)
  {
    case PROP_SIGMA_POP:
      g_value_set_double (value, nc_galaxy_wl_binned_estimator_gt_get_sigma_pop (est_gt));
      break;
    case PROP_ERROR_TYPE:
      g_value_set_enum (value, nc_galaxy_wl_binned_estimator_gt_get_error_type (est_gt));
      break;
    case PROP_BETA_WEIGHTED:
      g_value_set_boolean (value, nc_galaxy_wl_binned_estimator_gt_get_beta_weighted (est_gt));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static GStrv
_nc_galaxy_wl_binned_estimator_gt_extra_columns (NcGalaxyWLBinnedEstimator *est)
{
  GStrv cols = g_new0 (gchar *, 5);

  cols[0] = g_strdup (NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S_MEAN);
  cols[1] = g_strdup (NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S2_MEAN);
  cols[2] = g_strdup (NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_RESPONSIVITY);
  cols[3] = g_strdup (NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_M_MEAN);

  return cols;
}

typedef struct _NcWLBEGtBin
{
  gdouble n;    /* number of galaxies */
  gdouble w;    /* sum w */
  gdouble w2;   /* sum w^2 */
  gdouble w2v;  /* sum w^2 var */
  gdouble w_r;  /* sum w R */
  gdouble w_g;  /* sum w g_t */
  gdouble w_g2; /* sum w g_t^2 */
  gdouble w_m;  /* sum w (1 + m) */
  gdouble w_dc; /* sum w (c_t - hat c_t) */
} NcWLBEGtBin;

static void
_nc_galaxy_wl_binned_estimator_gt_build_bins (NcGalaxyWLBinnedEstimator *est, const NcGalaxyWLBinnedEstimatorSources *src, NcGalaxyWLObsBinned *binned, GArray *factor)
{
  NcGalaxyWLBinnedEstimatorGt *est_gt = NC_GALAXY_WL_BINNED_ESTIMATOR_GT (est);
  NcWLBEGtBin *acc                    = g_new0 (NcWLBEGtBin, src->n_bins);
  NcmMatrix *cov                      = nc_galaxy_wl_obs_binned_peek_cov (binned);
  const gdouble sigma_pop2            = est_gt->sigma_pop * est_gt->sigma_pop;
  gdouble sum_w = 0.0, sum_w_m = 0.0, sum_w_c1 = 0.0, sum_w_c2 = 0.0, sum_w_erms = 0.0;
  gdouble sum_beta = 0.0, sum_beta2 = 0.0, sum_beta_w = 0.0;
  gdouble m_mean, c1_mean, c2_mean, resp, scale;
  guint i, b;

#define _GT_WEIGHT(i) (1.0 / (sigma_pop2 + 2.0 * gsl_pow_2 (NC_WLBE_SRC (src, std_noise, (i)))))

  /* Calibration over the whole sample (Eqs. 10-11) */
  for (i = 0; i < src->len; i++)
  {
    const gdouble w_i   = _GT_WEIGHT (i);
    const gdouble e_rms = NC_WLBE_SRC (src, e_rms, i);

    sum_w      += w_i;
    sum_w_m    += w_i * NC_WLBE_SRC (src, m, i);
    sum_w_c1   += w_i * NC_WLBE_SRC (src, c1, i);
    sum_w_c2   += w_i * NC_WLBE_SRC (src, c2, i);
    sum_w_erms += w_i * e_rms * e_rms;

    {
      /* Means weighted by the shape weights, or plain means */
      const gdouble wb = est_gt->beta_weighted ? w_i : 1.0;

      sum_beta_w += wb;
      sum_beta   += wb * src->sigma_crit_inf * NC_WLBE_SRC (src, inv_sigma_crit, i);
      sum_beta2  += wb * src->sigma_crit_inf * src->sigma_crit_inf * NC_WLBE_SRC (src, inv_sigma_crit2, i);
    }
  }

  m_mean  = sum_w_m / sum_w;
  c1_mean = sum_w_c1 / sum_w;
  c2_mean = sum_w_c2 / sum_w;
  resp    = src->use_responsivity ? 1.0 - sum_w_erms / sum_w : 1.0;
  scale   = src->shear_scale * resp;

  for (i = 0; i < src->len; i++)
  {
    const gdouble w_i     = _GT_WEIGHT (i);
    const gdouble hat_c_t = c1_mean * NC_WLBE_SRC (src, rot_re, i) + c2_mean * NC_WLBE_SRC (src, rot_im, i);
    const gdouble g_t     = (NC_WLBE_SRC (src, e_t, i) / scale - hat_c_t) / (1.0 + m_mean);
    const gdouble e_rms   = NC_WLBE_SRC (src, e_rms, i);
    const gdouble std_n   = NC_WLBE_SRC (src, std_noise, i);
    const gdouble var_i   = (e_rms * e_rms + std_n * std_n) / gsl_pow_2 (scale * (1.0 + m_mean));
    NcWLBEGtBin *a        = &acc[NC_WLBE_SRC_BIN (src, i)];

    a->n    += 1.0;
    a->w    += w_i;
    a->w2   += w_i * w_i;
    a->w2v  += w_i * w_i * var_i;
    a->w_r  += w_i * NC_WLBE_SRC (src, radius, i);
    a->w_g  += w_i * g_t;
    a->w_g2 += w_i * g_t * g_t;
    a->w_m  += w_i * (1.0 + NC_WLBE_SRC (src, m, i));
    a->w_dc += w_i * (NC_WLBE_SRC (src, c_t, i) - hat_c_t);
  }

  ncm_matrix_set_zero (cov);

  for (b = 0; b < src->n_bins; b++)
  {
    const NcWLBEGtBin *a = &acc[b];
    const gdouble g_mean = a->w_g / a->w;
    gdouble var;

    switch (est_gt->error_type)
    {
      case NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_STD_ERR:
        /*
         * The standard error of the weighted mean: the weighted variance of the bin times
         * sum w^2 / (sum w)^2, which is 1 / n for equal weights. A single galaxy has no
         * scatter: NaN, the bin is dropped.
         */
        var = (a->n > 1.0) ? (a->w_g2 / a->w - g_mean * g_mean) * a->w2 / (a->w * a->w) : GSL_NAN;
        break;
      case NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_SHAPE_NOISE:
        var = a->w2v / (a->w * a->w);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        var = GSL_NAN;           /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }

    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_RADIUS, b, a->w_r / a->w);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_VALUE, b, g_mean);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_OFFSET, b, a->w_dc / ((1.0 + m_mean) * a->w));
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S_MEAN, b, sum_beta / sum_beta_w);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S2_MEAN, b, sum_beta2 / sum_beta_w);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_RESPONSIVITY, b, resp);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_M_MEAN, b, m_mean);

    ncm_matrix_set (cov, b, b, var);
  }

  /*
   * Expected value: <e_t> / (shear_scale R) = (1 + m_i) g_i + c_t,i, so
   * <g_t,i> = [(1 + m_i) g_i + c_t,i - hat c_t,i] / (1 + hat m).
   */
  for (i = 0; i < src->len; i++)
  {
    const NcWLBEGtBin *a = &acc[NC_WLBE_SRC_BIN (src, i)];

    g_array_index (factor, gdouble, i) = _GT_WEIGHT (i) * (1.0 + NC_WLBE_SRC (src, m, i)) / ((1.0 + m_mean) * a->w);
  }

#undef _GT_WEIGHT

  g_free (acc);
}

static void
_nc_galaxy_wl_binned_estimator_gt_eval_model (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObsBinned *binned, NcmVector *model)
{
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *hp          = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  const guint n_bins          = nc_galaxy_wl_obs_binned_len (binned);
  gdouble z_cl, sigma_crit_inf;
  guint b;

  if (!nc_galaxy_wl_obs_binned_has_column (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S_MEAN) ||
      !nc_galaxy_wl_obs_binned_has_column (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S2_MEAN))
    g_error ("nc_galaxy_wl_binned_estimator_gt: the profile has no columns `%s' and `%s'.",
             NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S_MEAN, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S2_MEAN);

  nc_halo_position_prepare_if_needed (hp, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (smd, cosmo);
  z_cl           = nc_halo_position_get_redshift (hp);
  sigma_crit_inf = nc_wl_surface_mass_density_sigma_critical_infinity (smd, cosmo, z_cl, z_cl);

  for (b = 0; b < n_bins; b++)
  {
    const gdouble radius    = nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_OBS_BINNED_COL_RADIUS, b);
    const gdouble beta_s    = nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S_MEAN, b);
    const gdouble beta_s2   = nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S2_MEAN, b);
    const gdouble kappa_inf = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, radius, z_cl) / sigma_crit_inf;
    const gdouble gamma_inf = nc_wl_surface_mass_density_sigma_excess (smd, dp, cosmo, radius, z_cl) / sigma_crit_inf;
    const gdouble bk        = beta_s * kappa_inf;

    ncm_vector_set (model, b, beta_s * gamma_inf / (1.0 - bk) * (1.0 + (beta_s2 / (beta_s * beta_s) - 1.0) * bk));
  }
}

static void
nc_galaxy_wl_binned_estimator_gt_class_init (NcGalaxyWLBinnedEstimatorGtClass *klass)
{
  GObjectClass *object_class                = G_OBJECT_CLASS (klass);
  NcGalaxyWLBinnedEstimatorClass *est_class = NC_GALAXY_WL_BINNED_ESTIMATOR_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_binned_estimator_gt_set_property;
  object_class->get_property = &_nc_galaxy_wl_binned_estimator_gt_get_property;

  /**
   * NcGalaxyWLBinnedEstimatorGt:sigma-pop:
   *
   * The shape noise of the weights $w_i = 1 / (\sigma_\mathrm{pop}^2 + 2\sigma_{e,i}^2)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGMA_POP,
                                   g_param_spec_double ("sigma-pop",
                                                        NULL,
                                                        "Intrinsic ellipticity dispersion of the weights",
                                                        0.0, 1.0, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_DEFAULT_SIGMA_POP,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLBinnedEstimatorGt:error-type:
   *
   * The error of each bin.
   */
  g_object_class_install_property (object_class,
                                   PROP_ERROR_TYPE,
                                   g_param_spec_enum ("error-type",
                                                      NULL,
                                                      "Error of each bin",
                                                      NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR,
                                                      NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_STD_ERR,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLBinnedEstimatorGt:beta-weighted:
   *
   * Whether $\langle\beta_s\rangle$ and $\langle\beta_s^2\rangle$ are means weighted by
   * the shape weights $w_i$ instead of plain
   * means over the galaxies.
   */
  g_object_class_install_property (object_class,
                                   PROP_BETA_WEIGHTED,
                                   g_param_spec_boolean ("beta-weighted",
                                                         NULL,
                                                         "Whether the lensing efficiency moments are weighted",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  est_class->name          = "gt";
  est_class->extra_columns = &_nc_galaxy_wl_binned_estimator_gt_extra_columns;
  est_class->build_bins    = &_nc_galaxy_wl_binned_estimator_gt_build_bins;
  est_class->eval_model    = &_nc_galaxy_wl_binned_estimator_gt_eval_model;
}

/**
 * nc_galaxy_wl_binned_estimator_gt_new:
 *
 * Creates a new #NcGalaxyWLBinnedEstimatorGt.
 *
 * Returns: (transfer full): a new #NcGalaxyWLBinnedEstimatorGt.
 */
NcGalaxyWLBinnedEstimatorGt *
nc_galaxy_wl_binned_estimator_gt_new (void)
{
  return g_object_new (NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR_GT, NULL);
}

/**
 * nc_galaxy_wl_binned_estimator_gt_ref:
 * @est_gt: a #NcGalaxyWLBinnedEstimatorGt
 *
 * Increases the reference count of @est_gt by one.
 *
 * Returns: (transfer full): @est_gt.
 */
NcGalaxyWLBinnedEstimatorGt *
nc_galaxy_wl_binned_estimator_gt_ref (NcGalaxyWLBinnedEstimatorGt *est_gt)
{
  return g_object_ref (est_gt);
}

/**
 * nc_galaxy_wl_binned_estimator_gt_free:
 * @est_gt: a #NcGalaxyWLBinnedEstimatorGt
 *
 * Decreases the reference count of @est_gt by one.
 */
void
nc_galaxy_wl_binned_estimator_gt_free (NcGalaxyWLBinnedEstimatorGt *est_gt)
{
  g_object_unref (est_gt);
}

/**
 * nc_galaxy_wl_binned_estimator_gt_clear:
 * @est_gt: a #NcGalaxyWLBinnedEstimatorGt
 *
 * Decreases the reference count of *@est_gt by one and sets *@est_gt to %NULL.
 */
void
nc_galaxy_wl_binned_estimator_gt_clear (NcGalaxyWLBinnedEstimatorGt **est_gt)
{
  g_clear_object (est_gt);
}

/**
 * nc_galaxy_wl_binned_estimator_gt_set_sigma_pop:
 * @est_gt: a #NcGalaxyWLBinnedEstimatorGt
 * @sigma_pop: the intrinsic ellipticity dispersion of the weights
 *
 * Sets #NcGalaxyWLBinnedEstimatorGt:sigma-pop.
 */
void
nc_galaxy_wl_binned_estimator_gt_set_sigma_pop (NcGalaxyWLBinnedEstimatorGt *est_gt, const gdouble sigma_pop)
{
  est_gt->sigma_pop = sigma_pop;
}

/**
 * nc_galaxy_wl_binned_estimator_gt_get_sigma_pop:
 * @est_gt: a #NcGalaxyWLBinnedEstimatorGt
 *
 * Returns: #NcGalaxyWLBinnedEstimatorGt:sigma-pop.
 */
gdouble
nc_galaxy_wl_binned_estimator_gt_get_sigma_pop (NcGalaxyWLBinnedEstimatorGt *est_gt)
{
  return est_gt->sigma_pop;
}

/**
 * nc_galaxy_wl_binned_estimator_gt_set_error_type:
 * @est_gt: a #NcGalaxyWLBinnedEstimatorGt
 * @error_type: a #NcGalaxyWLBinnedEstimatorGtError
 *
 * Sets #NcGalaxyWLBinnedEstimatorGt:error-type.
 */
void
nc_galaxy_wl_binned_estimator_gt_set_error_type (NcGalaxyWLBinnedEstimatorGt *est_gt, NcGalaxyWLBinnedEstimatorGtError error_type)
{
  est_gt->error_type = error_type;
}

/**
 * nc_galaxy_wl_binned_estimator_gt_get_error_type:
 * @est_gt: a #NcGalaxyWLBinnedEstimatorGt
 *
 * Returns: #NcGalaxyWLBinnedEstimatorGt:error-type.
 */
NcGalaxyWLBinnedEstimatorGtError
nc_galaxy_wl_binned_estimator_gt_get_error_type (NcGalaxyWLBinnedEstimatorGt *est_gt)
{
  return est_gt->error_type;
}

/**
 * nc_galaxy_wl_binned_estimator_gt_set_beta_weighted:
 * @est_gt: a #NcGalaxyWLBinnedEstimatorGt
 * @beta_weighted: whether the lensing efficiency moments are weighted
 *
 * Sets #NcGalaxyWLBinnedEstimatorGt:beta-weighted.
 */
void
nc_galaxy_wl_binned_estimator_gt_set_beta_weighted (NcGalaxyWLBinnedEstimatorGt *est_gt, gboolean beta_weighted)
{
  est_gt->beta_weighted = beta_weighted;
}

/**
 * nc_galaxy_wl_binned_estimator_gt_get_beta_weighted:
 * @est_gt: a #NcGalaxyWLBinnedEstimatorGt
 *
 * Returns: #NcGalaxyWLBinnedEstimatorGt:beta-weighted.
 */
gboolean
nc_galaxy_wl_binned_estimator_gt_get_beta_weighted (NcGalaxyWLBinnedEstimatorGt *est_gt)
{
  return est_gt->beta_weighted;
}

