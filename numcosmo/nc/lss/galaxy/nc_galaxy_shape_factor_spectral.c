/***************************************************************************
 *            nc_galaxy_shape_factor_spectral.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_spectral.c
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
 * NcGalaxyShapeFactorSpectral:
 *
 * Intrinsic-ellipticity marginal by adaptive Chebyshev panels in polar coordinates.
 *
 * Computes
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - f_g(\chi_I);
 *   \sigma_\mathrm{noise}^2\big)$$
 * with #NcmSpectralPolar in the polar coordinates $(r, \theta)$ of the variable of
 * integration, for any isotropic population that implements
 * nc_galaxy_shape_pop_exponent_at_origin(), nc_galaxy_shape_pop_exponent_at_one() and
 * nc_galaxy_shape_pop_eval_p_reduced(). The two end power laws of the density of $r$ are
 * the weight of the integral and the reduced density times the noise kernel is its
 * smooth factor.
 *
 * The variable of integration is $\chi_I$ for %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET. For
 * %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE it is $\epsilon_I = T^{-1}(\chi_I)$,
 * $T(\epsilon) = 2\epsilon / (1 + |\epsilon|^2)$, the trace--determinant ellipticity of
 * the same shape, in which the shear map is the composition of the trace--determinant
 * map with $T$. In it the integrand varies on the scale $1 - |g|$ near $|g| = 1$, and in
 * $\chi_I$ on the scale $(1 - |g|)^2$.
 *
 * The shear enters through its modulus. At $|g| = 1$ every point of the open disk maps
 * to $g$ and $P = N_2(\epsilon_\mathrm{obs} - g; \sigma_\mathrm{noise}^2)$. For
 * $|g| > 1$ the integral equals the one at $g / |g|^2$, which is used. The frame is then
 * rotated by the phase of $g$, so that the shear is the real number $t = |g|$.
 *
 * For $1 - t < 10^{-5}$ the integral is not computed at $t$, where the integrand varies
 * on the scale $1 - t$ next to the unit circle. The value is the sum of the marginal to
 * second order in the distance of the sheared ellipticity to its mean and of a term
 * proportional to a power of $(1 - t) / (1 + t)$, whose coefficients are set by the
 * integral at $1 - t = 10^{-5}$, and for %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET also at
 * $1 - t = 10^{-4}$. It assumes $10^{-5} |\epsilon_\mathrm{obs} - g| \ll
 * \sigma_\mathrm{noise}^2$.
 *
 * The panels are those of ncm_spectral_polar_peak_edges() for a peak located on the ray
 * through the preimage of $\epsilon_\mathrm{obs}$, at the radius that maximizes the
 * logarithm of the density of $r$ minus
 * $|\epsilon_\mathrm{obs} - f_t(\chi_I)|^2 / 2\sigma_\mathrm{noise}^2$ on a scan of that ray, with half widths $\sigma_\mathrm{noise}$ times the smaller
 * stretch factor of the inverse map there. The scan places the panels only: the integral is
 * accepted by the test of #NcmSpectralPolar with the relative tolerance
 * #NcGalaxyShapeFactorSpectral:reltol on panels that cover the whole disk.
 *
 * The class computes $\ln P$, with the smallest $|\epsilon_\mathrm{obs} - f_t(\chi_I)|^2$
 * of the scan taken out of the exponent, and nc_galaxy_shape_factor_eval_marginal()
 * returns its exponential. A panel that reaches #NcmSpectralPolar:max-depth, or a
 * non-finite $\ln P$, is an error that aborts with the point in the message.
 *
 * With #NcGalaxyShapeFactorSpectral:use-interp, nc_galaxy_shape_factor_eval_marginal() and
 * nc_galaxy_shape_factor_eval_ln_marginal() take $\ln P$ for
 * $|g| \le$ #NcGalaxyShapeFactorSpectral:interp-g-max from an approximation of $\ln P$ in
 * $g$ built for each galaxy by ncm_spectral_polar_build() from the integral, to the
 * tolerance #NcGalaxyShapeFactorSpectral:interp-abstol, and kept in the
 * #NcGalaxyShapeFactorData of the galaxy.
 *
 * A call takes its scratch space from a pool, so one object may be used for several
 * galaxies in several threads. Two threads must not evaluate the same galaxy when
 * #NcGalaxyShapeFactorSpectral:use-interp is set, since the first evaluation builds its
 * interpolant.
 *
 * See <a href="../../theory/nc/lss/galaxy/wl_shape_marginalization_spectral.html">Shape
 * Marginalization by Spectral Panels</a> for the formulas.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_spectral.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"
#include "ncm/algebra/ncm_spectral_polar.h"
#include "ncm/algebra/ncm_spectral_polar_panels.h"
#include "ncm/core/ncm_memory_pool.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* The radii of the scan of the ray: a uniform grid, the two ends, and the peak of the noise */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM (401)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_END_DECADES (11)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_OFFSETS (10)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA (2 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_END_DECADES + 4 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_OFFSETS)

/*
 * Below this distance of min (|g|, 1 / |g|) to 1 the integral is not computed at the shear;
 * the value is formed from integrals at this distance. The errors measured against the Arb
 * reference hold for this value.
 */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA (1.0e-5)

/* Steps of the golden-section search of the maximum of the scan */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_GOLDEN_STEPS (40)

struct _NcGalaxyShapeFactorSpectral
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorSpectralKernel NcGalaxyShapeFactorSpectralKernel;

/*
 * The four functions of the object's convention act on the variable of integration z:
 * chi_I for TRACE_DET and T^{-1} (chi_I) for TRACE, T (e) = 2 e / (1 + |e|^2).
 */
typedef struct _NcGalaxyShapeFactorSpectralPrivate
{
  NcmSpectralPolarF kernel; /* the integrand */

  /* z to the sheared ellipticity, at the real shear t */
  void (*map) (gdouble t, gdouble x, gdouble y, gdouble * restrict f_x, gdouble * restrict f_y);

  /* the sheared ellipticity e_obs to z, at the shear g */
  complex double (*map_inv) (complex double g, complex double e_obs);

  /* the density of r = |z| divided by its two end power laws */
  gdouble (*p_reduced) (const NcGalaxyShapeFactorSpectralKernel *k, const gdouble r);

  gboolean axis_ratio_var; /* TRUE for TRACE: z = T^{-1} (chi_I) */

  NcmSpectralPolar *polar;        /* the integrals */
  NcmSpectralPolar *polar_interp; /* the build of the interpolant in g */
  NcmMemoryPool *work_pool;       /* NcGalaxyShapeFactorSpectralWork, one per concurrent evaluation */
  gdouble reltol;
  gboolean use_interp;
  gdouble interp_g_max;
  gdouble interp_abstol;
  guint interp_min_radial_panels;
} NcGalaxyShapeFactorSpectralPrivate;

/* The scratch space of one evaluation */
typedef struct _NcGalaxyShapeFactorSpectralWork
{
  GArray *r_edges;
  GArray *theta_edges;
  GArray *r_edges_all; /* edges of a peak before those outside a region are left out */
  GArray *theta_edges_all;
} NcGalaxyShapeFactorSpectralWork;

/*
 * The interpolant of ln P in g of one galaxy and what it was built for. It is rebuilt when
 * one of these values differs from the one of the evaluation.
 */
typedef struct _NcGalaxyShapeFactorSpectralData
{
  NcmSpectralPolarPanels *spp;
  gboolean valid;
  gdouble epsilon_obs_1;
  gdouble epsilon_obs_2;
  gdouble std_noise;
  guint64 pop_pkey;
  gdouble g_max;
  gdouble abstol;
  gdouble reltol;
  guint min_radial_panels;
  gdouble cos_e; /* direction of eps_obs */
  gdouble sin_e;
  gdouble g_c; /* center of the coordinates of the interpolant */
} NcGalaxyShapeFactorSpectralData;

/* What the integrand needs, in the frame where the shear is the real number t */
struct _NcGalaxyShapeFactorSpectralKernel
{
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  gdouble t;
  gdouble e_1;
  gdouble e_2;
  gdouble c2;        /* 1 / (2 std_noise^2) */
  gdouble norm;      /* 1 / (2 pi std_noise^2) / (2 pi) */
  gdouble shift_c2;  /* c2 times the offset taken out of the exponent */
  gdouble dens_norm; /* TRACE: 2^(alpha_o + 1), alpha_o and alpha_1 the exponents of the population */
  gdouble dens_pow;  /* TRACE: alpha_o + alpha_1 + 2 */
};

/* The peak of the integrand on the ray through the preimage of eps_obs */
typedef struct _NcGalaxyShapeFactorSpectralSeed
{
  gdouble r_star;  /* radius of the maximum of h on the scan */
  gdouble theta_0; /* angle of the ray */
  gdouble J_1;     /* radial stretch factor of the inverse map at eps_obs */
  gdouble J_2;     /* tangential stretch factor */
  gdouble shift;   /* smallest squared distance of eps_obs to the sheared ellipticity on the scan */
} NcGalaxyShapeFactorSpectralSeed;

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_USE_INTERP,
  PROP_INTERP_G_MAX,
  PROP_INTERP_ABSTOL,
  PROP_INTERP_MIN_RADIAL_PANELS,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorSpectral, nc_galaxy_shape_factor_spectral, NC_TYPE_GALAXY_SHAPE_FACTOR)

static gpointer _nc_galaxy_shape_factor_spectral_work_new (gpointer userdata);
static void _nc_galaxy_shape_factor_spectral_work_free (gpointer p);

static void
nc_galaxy_shape_factor_spectral_init (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->kernel         = NULL;
  self->map            = NULL;
  self->map_inv        = NULL;
  self->p_reduced      = NULL;
  self->axis_ratio_var = FALSE;
  self->polar          = ncm_spectral_polar_new ();
  self->polar_interp   = ncm_spectral_polar_new ();
  self->work_pool      = ncm_memory_pool_new (_nc_galaxy_shape_factor_spectral_work_new, NULL,
                                              _nc_galaxy_shape_factor_spectral_work_free);
  self->reltol        = 0.0;
  self->use_interp    = FALSE;
  self->interp_g_max  = 0.0;
  self->interp_abstol = 0.0;

  self->interp_min_radial_panels = 0;
}

static void
_nc_galaxy_shape_factor_spectral_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorSpectral *gsfs = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);

  switch (prop_id)
  {
    case PROP_RELTOL:
      nc_galaxy_shape_factor_spectral_set_reltol (gsfs, g_value_get_double (value));
      break;
    case PROP_USE_INTERP:
      nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, g_value_get_boolean (value));
      break;
    case PROP_INTERP_G_MAX:
      nc_galaxy_shape_factor_spectral_set_interp_g_max (gsfs, g_value_get_double (value));
      break;
    case PROP_INTERP_ABSTOL:
      nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, g_value_get_double (value));
      break;
    case PROP_INTERP_MIN_RADIAL_PANELS:
      nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (gsfs, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_spectral_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorSpectral *gsfs = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_reltol (gsfs));
      break;
    case PROP_USE_INTERP:
      g_value_set_boolean (value, nc_galaxy_shape_factor_spectral_get_use_interp (gsfs));
      break;
    case PROP_INTERP_G_MAX:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_interp_g_max (gsfs));
      break;
    case PROP_INTERP_ABSTOL:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_interp_abstol (gsfs));
      break;
    case PROP_INTERP_MIN_RADIAL_PANELS:
      g_value_set_uint (value, nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (gsfs));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _nc_galaxy_shape_factor_spectral_kernel_trace (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval);
static void _nc_galaxy_shape_factor_spectral_kernel_trace_det (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval);
static void _nc_galaxy_shape_factor_spectral_map_trace (gdouble t, gdouble x, gdouble y, gdouble * restrict f_x, gdouble * restrict f_y);

static complex double _nc_galaxy_shape_factor_spectral_map_inv_trace (complex double g, complex double e_obs);

static gdouble _nc_galaxy_shape_factor_spectral_p_reduced_trace (const NcGalaxyShapeFactorSpectralKernel *k, const gdouble r);
static gdouble _nc_galaxy_shape_factor_spectral_p_reduced_trace_det (const NcGalaxyShapeFactorSpectralKernel *k, const gdouble r);

static void
_nc_galaxy_shape_factor_spectral_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_spectral_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);
    NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
    const NcGalaxyWLObsEllipConv ellip_conv         = nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfs));

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        self->kernel         = &_nc_galaxy_shape_factor_spectral_kernel_trace;
        self->map            = &_nc_galaxy_shape_factor_spectral_map_trace;
        self->map_inv        = &_nc_galaxy_shape_factor_spectral_map_inv_trace;
        self->p_reduced      = &_nc_galaxy_shape_factor_spectral_p_reduced_trace;
        self->axis_ratio_var = TRUE;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->kernel         = &_nc_galaxy_shape_factor_spectral_kernel_trace_det;
        self->map            = &nc_wl_ellipticity_apply_shear_trace_det_real;
        self->map_inv        = &nc_wl_ellipticity_apply_shear_inv_trace_det;
        self->p_reduced      = &_nc_galaxy_shape_factor_spectral_p_reduced_trace_det;
        self->axis_ratio_var = FALSE;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }
  }
}

static void
_nc_galaxy_shape_factor_spectral_dispose (GObject *object)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  ncm_spectral_polar_clear (&self->polar);
  ncm_spectral_polar_clear (&self->polar_interp);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_spectral_parent_class)->dispose (object);
}

static void
_nc_galaxy_shape_factor_spectral_finalize (GObject *object)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  ncm_memory_pool_free (self->work_pool, TRUE);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_spectral_parent_class)->finalize (object);
}

static void _nc_galaxy_shape_factor_spectral_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data);
static void _nc_galaxy_shape_factor_spectral_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset);
static gdouble _nc_galaxy_shape_factor_spectral_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static gdouble _nc_galaxy_shape_factor_spectral_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);

static void
nc_galaxy_shape_factor_spectral_class_init (NcGalaxyShapeFactorSpectralClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_nc_galaxy_shape_factor_spectral_constructed;
  object_class->set_property = &_nc_galaxy_shape_factor_spectral_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_spectral_get_property;
  object_class->dispose      = &_nc_galaxy_shape_factor_spectral_dispose;
  object_class->finalize     = &_nc_galaxy_shape_factor_spectral_finalize;

  /**
   * NcGalaxyShapeFactorSpectral:reltol:
   *
   * The relative tolerance of the integral, passed to ncm_spectral_polar_integrate().
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance of the integral",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:use-interp:
   *
   * Whether nc_galaxy_shape_factor_eval_marginal() and
   * nc_galaxy_shape_factor_eval_ln_marginal() take $\ln P$ from the interpolant in $g$ of
   * the galaxy for $|g| \le$ #NcGalaxyShapeFactorSpectral:interp-g-max. The interpolant is
   * built at the first such evaluation of a galaxy and again when its observed
   * ellipticity, its noise, the parameters of the population or a tolerance changes. A
   * build computes of the order of $10^3$ integrals, so it pays when the population
   * parameters are fixed.
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_INTERP,
                                   g_param_spec_boolean ("use-interp",
                                                         NULL,
                                                         "Whether ln P is interpolated in g for each galaxy",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-g-max:
   *
   * The largest $|g|$ computed from the interpolant. It is below $1/3$: the region of the
   * interpolant reaches $|g| < 3$ times this value and stays inside $|g| < 1$.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_G_MAX,
                                   g_param_spec_double ("interp-g-max",
                                                        NULL,
                                                        "Largest |g| computed from the interpolant",
                                                        1.0e-3, 0.33, 0.3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-abstol:
   *
   * The absolute tolerance on $\ln P$ of the interpolant, passed to
   * ncm_spectral_polar_build(): a panel is accepted when $\ln P$ changes by at most this
   * value between its two grids.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_ABSTOL,
                                   g_param_spec_double ("interp-abstol",
                                                        NULL,
                                                        "Absolute tolerance on ln P of the interpolant",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-min-radial-panels:
   *
   * The #NcmSpectralPolar:min-radial-panels of the build of the interpolant: the radial
   * range of the interpolant is split into at least this many equal intervals. More
   * intervals give panels of lower degree, so a faster evaluation of the interpolant, at the
   * cost of more integrals in the build.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_MIN_RADIAL_PANELS,
                                   g_param_spec_uint ("interp-min-radial-panels",
                                                      NULL,
                                                      "Minimum number of radial intervals of the interpolant",
                                                      1, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_spectral_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_spectral_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_spectral_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_spectral_eval_ln_marginal;
}

static gpointer
_nc_galaxy_shape_factor_spectral_work_new (gpointer userdata)
{
  NcGalaxyShapeFactorSpectralWork *work = g_new0 (NcGalaxyShapeFactorSpectralWork, 1);

  work->r_edges         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->theta_edges     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->r_edges_all     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->theta_edges_all = g_array_new (FALSE, FALSE, sizeof (gdouble));

  return work;
}

static void
_nc_galaxy_shape_factor_spectral_work_free (gpointer p)
{
  NcGalaxyShapeFactorSpectralWork *work = p;

  g_array_unref (work->r_edges);
  g_array_unref (work->theta_edges);
  g_array_unref (work->r_edges_all);
  g_array_unref (work->theta_edges_all);

  g_free (work);
}

static void
_nc_galaxy_shape_factor_spectral_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_factor_spectral_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void _nc_galaxy_shape_factor_spectral_ldata_free (gpointer p);

/*
 * The state per galaxy is the interpolant of ln P in g, empty until
 * NcGalaxyShapeFactorSpectral:use-interp builds it. Nothing of it is a column of the
 * catalog.
 */
static void
_nc_galaxy_shape_factor_spectral_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorSpectralData *ldata = g_new0 (NcGalaxyShapeFactorSpectralData, 1);

  ldata->spp   = NULL;
  ldata->valid = FALSE;

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_spectral_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_spectral_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_spectral_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_spectral_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_spectral_ldata_free (gpointer p)
{
  NcGalaxyShapeFactorSpectralData *ldata = p;

  ncm_spectral_polar_panels_clear (&ldata->spp);
  g_free (ldata);
}

static void
_nc_galaxy_shape_factor_spectral_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
}

static gdouble _nc_galaxy_shape_factor_spectral_eval_ln (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble std_noise, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, NcmMatrix **panels);
static void _nc_galaxy_shape_factor_spectral_interp_build (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);

static gdouble
_nc_galaxy_shape_factor_spectral_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return exp (_nc_galaxy_shape_factor_spectral_eval_ln_marginal (gsf, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

/*
 * ln P (eps_obs | g) for the galaxy of data. With NcGalaxyShapeFactorSpectral:use-interp
 * and |g| <= NcGalaxyShapeFactorSpectral:interp-g-max it is the value of the interpolant
 * of the galaxy in g, which is built first when the galaxy has none for these eps_obs,
 * std_noise, population parameters and tolerances. Otherwise it is the value of
 * _nc_galaxy_shape_factor_spectral_eval_ln(). The interpolant is evaluated at
 * (rho, psi), g = g_c + rho exp (i psi) in the frame where eps_obs is on the positive
 * real axis, with the sign of the imaginary part of g dropped: ln P is even under that
 * reflection.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  if (!self->use_interp || (g_1 * g_1 + g_2 * g_2 > gsl_pow_2 (self->interp_g_max)))
    return _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, data->pop_data, data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, NULL);

  {
    NcGalaxyShapeFactorSpectralData *ldata = data->ldata;

    if (!ldata->valid || (ldata->epsilon_obs_1 != epsilon_obs_1) || (ldata->epsilon_obs_2 != epsilon_obs_2) ||
        (ldata->std_noise != data->std_noise) || (ldata->pop_pkey != ncm_model_state_get_pkey (NCM_MODEL (pop))) ||
        (ldata->g_max != self->interp_g_max) || (ldata->abstol != self->interp_abstol) || (ldata->reltol != self->reltol) ||
        (ldata->min_radial_panels != self->interp_min_radial_panels))
      _nc_galaxy_shape_factor_spectral_interp_build (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2);

    {
      const gdouble gf_1 = g_1 * ldata->cos_e + g_2 * ldata->sin_e - ldata->g_c;
      const gdouble gf_2 = g_2 * ldata->cos_e - g_1 * ldata->sin_e;

      return ncm_spectral_polar_panels_eval (ldata->spp, hypot (gf_1, gf_2), atan2 (fabs (gf_2), gf_1));
    }
  }
}

static void _nc_galaxy_shape_factor_spectral_interp_kernel (gpointer user_data, gdouble rho, guint n, const gdouble *psi, const gdouble *cos_psi, const gdouble *sin_psi, gdouble *fval);

/* What the function approximated in g needs */
typedef struct _NcGalaxyShapeFactorSpectralInterp
{
  NcGalaxyShapeFactorSpectral *gsfs;
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  gdouble std_noise;
  gdouble e;   /* |eps_obs| */
  gdouble g_c; /* center of the coordinates */
} NcGalaxyShapeFactorSpectralInterp;

/*
 * Builds the interpolant of ln P in g of the galaxy of data, for the observed ellipticity
 * (epsilon_obs_1, epsilon_obs_2), on a region that contains the disk
 * |g| <= NcGalaxyShapeFactorSpectral:interp-g-max, and records in the data of the galaxy
 * what it was built for. In the frame where eps_obs is on the positive real axis,
 * g = g_c + rho exp (i psi), with g_c the shear that maps chi_I = 0 to eps_obs:
 * e / (1 + sqrt (1 - e^2)) for TRACE (1 for e >= 1) and e for TRACE_DET, e = |eps_obs|.
 * ln P has its narrowest structure at g_c, of width std_noise / 2 (TRACE) or std_noise
 * (TRACE_DET). The region is [0, g_max + g_c] x [0, pi] for g_c < g_max and otherwise
 * [g_c - g_max, g_c + g_max] x [pi - asin (g_max / g_c), pi]; psi < 0 is covered by the
 * reflection. The radial edges are those of ncm_spectral_polar_peak_edges() for a peak at
 * rho = 0 of that width that are inside the region, and the angular panels span pi / 2 at
 * most. ncm_spectral_polar_build() accepts a panel when the change of ln P between its two
 * grids is at most NcGalaxyShapeFactorSpectral:interp-abstol. A build that reaches the
 * maximum depth aborts with the galaxy in the message.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_build (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                               const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  NcGalaxyShapeFactorSpectralWork **ptr           = ncm_memory_pool_get (self->work_pool);
  NcGalaxyShapeFactorSpectralWork *work           = *ptr;
  const gdouble g_max                             = self->interp_g_max;
  const gdouble e                                 = hypot (epsilon_obs_1, epsilon_obs_2);
  const gdouble g_c                               = self->axis_ratio_var ? ((e < 1.0) ? e / (1.0 + sqrt (1.0 - e * e)) : 1.0) : e;
  const gdouble L                                 = (self->axis_ratio_var ? 0.5 : 1.0) * data->std_noise;
  const gdouble rho_min                           = (g_c < g_max) ? 0.0 : g_c - g_max;
  const gdouble rho_max                           = g_c + g_max;
  const gdouble psi_min                           = (g_c < g_max) ? 0.0 : M_PI - asin (g_max / g_c);
  const gdouble psi_max                           = M_PI;
  NcGalaxyShapeFactorSpectralInterp interp        = {gsfs, pop, data->pop_data, data->std_noise, e, g_c};
  GError *error                                   = NULL;
  guint i;

  ncm_spectral_polar_peak_edges (self->polar_interp, rho_max, 0.0, 0.5 * M_PI, L, M_PI, work->r_edges_all, work->theta_edges_all);

  g_array_set_size (work->r_edges, 0);
  g_array_append_val (work->r_edges, rho_min);

  for (i = 0; i < work->r_edges_all->len; i++)
  {
    const gdouble rho = g_array_index (work->r_edges_all, gdouble, i);

    if (rho > rho_min + 0.125 * data->std_noise)
      g_array_append_val (work->r_edges, rho);
  }

  g_array_set_size (work->theta_edges, 0);
  g_array_append_val (work->theta_edges, psi_min);

  if (psi_max - psi_min > 0.5 * M_PI)
  {
    const gdouble psi_mid = 0.5 * (psi_min + psi_max);

    g_array_append_val (work->theta_edges, psi_mid);
  }

  g_array_append_val (work->theta_edges, psi_max);

  if (ldata->spp == NULL)
    ldata->spp = ncm_spectral_polar_panels_new ();

  ldata->valid = FALSE;
  ncm_spectral_polar_build (self->polar_interp, &_nc_galaxy_shape_factor_spectral_interp_kernel, &interp, work->r_edges, work->theta_edges,
                            0.0, self->interp_abstol, ldata->spp, NULL, &error);

  ncm_memory_pool_return (ptr);

  if (error != NULL)
    g_error ("nc_galaxy_shape_factor_spectral: the interpolant in g: %s for the population %s, std_noise = % 22.15g, "
             "eps_obs = (% 22.15g, % 22.15g), interp-g-max = % 22.15g.",
             error->message, G_OBJECT_TYPE_NAME (pop), data->std_noise, epsilon_obs_1, epsilon_obs_2, g_max);

  ldata->valid             = TRUE;
  ldata->epsilon_obs_1     = epsilon_obs_1;
  ldata->epsilon_obs_2     = epsilon_obs_2;
  ldata->std_noise         = data->std_noise;
  ldata->pop_pkey          = ncm_model_state_get_pkey (NCM_MODEL (pop));
  ldata->g_max             = g_max;
  ldata->abstol            = self->interp_abstol;
  ldata->reltol            = self->reltol;
  ldata->min_radial_panels = self->interp_min_radial_panels;
  ldata->cos_e             = (e > 0.0) ? epsilon_obs_1 / e : 1.0;
  ldata->sin_e             = (e > 0.0) ? epsilon_obs_2 / e : 0.0;
  ldata->g_c               = g_c;
}

/*
 * The function approximated by ncm_spectral_polar_build(): ln P at the shear
 * g_c + rho exp (i psi) for the observed ellipticity (e, 0), by the integral.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_kernel (gpointer user_data, gdouble rho, guint n, const gdouble *psi,
                                                const gdouble *cos_psi, const gdouble *sin_psi, gdouble *fval)
{
  const NcGalaxyShapeFactorSpectralInterp *interp = user_data;
  guint j;

  for (j = 0; j < n; j++)
    fval[j] = _nc_galaxy_shape_factor_spectral_eval_ln (interp->gsfs, interp->pop, interp->pop_data, interp->std_noise,
                                                        interp->g_c + rho * cos_psi[j], rho * sin_psi[j], interp->e, 0.0, NULL);
}

static gdouble _nc_galaxy_shape_factor_spectral_one_minus_abs2 (const gdouble g_1, const gdouble g_2);
static gdouble _nc_galaxy_shape_factor_spectral_integrate_ln (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralKernel *k, const gdouble std_noise, const gdouble a, const gdouble b, NcGalaxyShapeFactorSpectralWork *work, NcmMatrix **panels, GError **error);
static gdouble _nc_galaxy_shape_factor_spectral_near_unit_ln (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralKernel *k, const gdouble std_noise, const gdouble a, const gdouble b, const gdouble delta, const gdouble lambda_inv, NcGalaxyShapeFactorSpectralWork *work, NcmMatrix **panels, GError **error);

/*
 * ln P (eps_obs | g) for the population pop with resolved data pop_data: the closed form
 * at |g| = 1, the integral at t = min (|g|, 1 / |g|) for
 * 1 - t >= NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA, and below it the second-order
 * term plus the term scaled from the integral at that distance from 1. With panels not
 * NULL, also the accepted panels of the integral computed, with their angles in the frame
 * of g and eps_obs, or NULL at |g| = 1, where nothing is integrated. A failed integral or
 * a non-finite ln P aborts with the point in the message.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_eval_ln (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                          const gdouble std_noise, const gdouble g_1, const gdouble g_2,
                                          const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, NcmMatrix **panels)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  const gdouble noise_var                         = gsl_pow_2 (std_noise);
  const gdouble abs_g2                            = g_1 * g_1 + g_2 * g_2;

  if (abs_g2 == 1.0)
  {
    if (panels != NULL)
      *panels = NULL;

    return -(gsl_pow_2 (epsilon_obs_1 - g_1) + gsl_pow_2 (epsilon_obs_2 - g_2)) / (2.0 * noise_var) - log (2.0 * M_PI * noise_var);
  }

  {
    /* For |g| > 1 the integral is the one at g / |g|^2, of the same phase */
    const gdouble abs_g                   = sqrt (abs_g2);
    const gdouble u                       = fabs (_nc_galaxy_shape_factor_spectral_one_minus_abs2 (g_1, g_2));
    const gdouble delta                   = (abs_g2 > 1.0) ? u / (abs_g * (1.0 + abs_g)) : u / (1.0 + abs_g);
    const gdouble lambda_inv              = u / gsl_pow_2 (1.0 + abs_g);
    const gdouble phase                   = (abs_g2 > 0.0) ? atan2 (g_2, g_1) : 0.0;
    const gdouble cos_p                   = cos (phase);
    const gdouble sin_p                   = sin (phase);
    const gdouble a                       = nc_galaxy_shape_pop_exponent_at_origin (pop);
    const gdouble b_pop                   = nc_galaxy_shape_pop_exponent_at_one (pop);
    const gdouble b                       = self->axis_ratio_var ? 2.0 * b_pop + 1.0 : b_pop;
    NcGalaxyShapeFactorSpectralWork **ptr = ncm_memory_pool_get (self->work_pool);
    NcGalaxyShapeFactorSpectralWork *work = *ptr;
    NcGalaxyShapeFactorSpectralKernel k   = {
      .pop       = pop,
      .pop_data  = pop_data,
      .t         = (abs_g2 > 1.0) ? 1.0 / abs_g : abs_g,
      .e_1       = epsilon_obs_1 * cos_p + epsilon_obs_2 * sin_p,
      .e_2       = epsilon_obs_2 * cos_p - epsilon_obs_1 * sin_p,
      .c2        = 0.5 / noise_var,
      .norm      = 1.0 / (2.0 * M_PI * noise_var) / (2.0 * M_PI),
      .shift_c2  = 0.0,
      .dens_norm = exp2 (a + 1.0),
      .dens_pow  = a + b_pop + 2.0
    };
    GError *error = NULL;
    gdouble ln_P;

    if (delta < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA)
      ln_P = _nc_galaxy_shape_factor_spectral_near_unit_ln (self, &k, std_noise, a, b, delta, lambda_inv, work, panels, &error);
    else
      ln_P = _nc_galaxy_shape_factor_spectral_integrate_ln (self, &k, std_noise, a, b, work, panels, &error);

    ncm_memory_pool_return (ptr);

    if ((error != NULL) || !gsl_finite (ln_P))
      g_error ("nc_galaxy_shape_factor_spectral: %s for the population %s, std_noise = % 22.15g, "
               "g = (% 22.15g, % 22.15g), eps_obs = (% 22.15g, % 22.15g).",
               (error != NULL) ? error->message : "non-finite ln P",
               G_OBJECT_TYPE_NAME (pop), std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2);

    if ((panels != NULL) && (phase != 0.0))
    {
      guint i;

      for (i = 0; i < ncm_matrix_nrows (*panels); i++)
      {
        ncm_matrix_addto (*panels, i, 2, phase);
        ncm_matrix_addto (*panels, i, 3, phase);
      }
    }

    return ln_P;
  }
}

/*
 * 1 - g_1^2 - g_2^2 with the rounding errors of the two squares and of their sum removed:
 * the squares are split into their rounded value and its error by fma(), the larger
 * rounded square is subtracted from 1 first, and the errors last. Next to |g| = 1 both
 * subtractions are exact, so the result keeps the digits of the inputs.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_one_minus_abs2 (const gdouble g_1, const gdouble g_2)
{
  const gdouble p   = g_1 * g_1;
  const gdouble q   = g_2 * g_2;
  const gdouble e_p = fma (g_1, g_1, -p);
  const gdouble e_q = fma (g_2, g_2, -q);
  const gdouble hi  = GSL_MAX (p, q);
  const gdouble lo  = GSL_MIN (p, q);

  return ((1.0 - hi) - lo) - (e_p + e_q);
}

static void _nc_galaxy_shape_factor_spectral_seed (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k, const gdouble std_noise, const gdouble a, const gdouble b, NcGalaxyShapeFactorSpectralSeed *seed);

/*
 * ln of the integral at the real shear k->t, by ncm_spectral_polar_integrate() on the
 * panels placed by the seed, with the offset of the seed set in k->shift_c2. NaN with
 * error set when a panel reaches the maximum depth.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_integrate_ln (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralKernel *k,
                                               const gdouble std_noise, const gdouble a, const gdouble b,
                                               NcGalaxyShapeFactorSpectralWork *work, NcmMatrix **panels, GError **error)
{
  NcGalaxyShapeFactorSpectralSeed seed;
  gdouble L, integral;

  k->shift_c2 = 0.0;
  _nc_galaxy_shape_factor_spectral_seed (self, k, std_noise, a, b, &seed);
  k->shift_c2 = seed.shift * k->c2;

  L = std_noise * GSL_MIN (seed.J_1, seed.J_2);
  ncm_spectral_polar_peak_edges (self->polar, 1.0, seed.r_star, seed.theta_0, L, L / GSL_MAX (seed.r_star, GSL_DBL_MIN),
                                 work->r_edges, work->theta_edges);

  integral = ncm_spectral_polar_integrate (self->polar, self->kernel, k, 1.0, a, b, work->r_edges, work->theta_edges,
                                           self->reltol, 0.0, NULL, panels, error);

  return log (integral) - k->shift_c2;
}

static gdouble _nc_galaxy_shape_factor_spectral_second_order_ln (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k, const gdouble std_noise, const gdouble a, const gdouble b, const gdouble delta, NcGalaxyShapeFactorSpectralWork *work, GError **error);

/*
 * ln P at the real shear t = 1 - delta, delta < delta_s =
 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA, without an integral at t:
 * P (t) = A (t) + B (t), with A the second-order term of
 * _nc_galaxy_shape_factor_spectral_second_order_ln() and
 * B (t) = lambda^-(b + 2) (k_0 + k_1 / lambda), lambda = (1 + t) / (1 - t) = 1 / lambda_inv.
 * For TRACE k_1 = 0 and k_0 is set by B = P - A at 1 - delta_s, with P the integral there.
 * For TRACE_DET k_0 and k_1 are set by B at 1 - delta_s and at 1 - 10 delta_s. Every
 * quantity is divided by P (1 - delta_s), so that none leaves the double range. panels,
 * when not NULL, receives the panels of the integral at 1 - delta_s. NaN when an integral
 * fails or when A + B is not positive.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_near_unit_ln (NcGalaxyShapeFactorSpectralPrivate * const self, NcGalaxyShapeFactorSpectralKernel *k,
                                               const gdouble std_noise, const gdouble a, const gdouble b,
                                               const gdouble delta, const gdouble lambda_inv,
                                               NcGalaxyShapeFactorSpectralWork *work, NcmMatrix **panels, GError **error)
{
  const gdouble t        = k->t;
  const gdouble t_1      = 1.0 - NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA;
  const gdouble delta_1  = 1.0 - t_1;
  const gdouble linv_1   = delta_1 / (1.0 + t_1);
  const gdouble exponent = b + 2.0;
  gdouble ln_P_1, ln_A_1, ln_A, b_1, k_0, k_1;

  k->t   = t_1;
  ln_P_1 = _nc_galaxy_shape_factor_spectral_integrate_ln (self, k, std_noise, a, b, work, panels, error);

  if (*error != NULL)
    return GSL_NAN;

  ln_A_1 = _nc_galaxy_shape_factor_spectral_second_order_ln (self, k, std_noise, a, b, delta_1, work, error);

  if (*error != NULL)
    return GSL_NAN;

  b_1 = -expm1 (ln_A_1 - ln_P_1) * pow (linv_1, -exponent);

  if (self->axis_ratio_var)
  {
    k_0 = b_1;
    k_1 = 0.0;
  }
  else
  {
    const gdouble t_2     = 1.0 - 10.0 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA;
    const gdouble delta_2 = 1.0 - t_2;
    const gdouble linv_2  = delta_2 / (1.0 + t_2);
    gdouble ln_P_2, ln_A_2, b_2;

    k->t   = t_2;
    ln_P_2 = _nc_galaxy_shape_factor_spectral_integrate_ln (self, k, std_noise, a, b, work, NULL, error);

    if (*error != NULL)
      return GSL_NAN;

    ln_A_2 = _nc_galaxy_shape_factor_spectral_second_order_ln (self, k, std_noise, a, b, delta_2, work, error);

    if (*error != NULL)
      return GSL_NAN;

    b_2 = (exp (ln_P_2 - ln_P_1) - exp (ln_A_2 - ln_P_1)) * pow (linv_2, -exponent);
    k_1 = (b_1 - b_2) / (linv_1 - linv_2);
    k_0 = b_1 - k_1 * linv_1;
  }

  k->t = t;
  ln_A = _nc_galaxy_shape_factor_spectral_second_order_ln (self, k, std_noise, a, b, delta, work, error);

  if (*error != NULL)
    return GSL_NAN;

  return ln_P_1 + log (exp (ln_A - ln_P_1) + pow (lambda_inv, exponent) * (k_0 + k_1 * lambda_inv));
}

static void _nc_galaxy_shape_factor_spectral_var_kernel (gpointer user_data, gdouble v, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval);

/* What the integrand of the variance needs */
typedef struct _NcGalaxyShapeFactorSpectralVar
{
  NcGalaxyShapeFactorSpectralPrivate *self;
  const NcGalaxyShapeFactorSpectralKernel *k;
  gdouble delta; /* 1 - t */
} NcGalaxyShapeFactorSpectralVar;

/*
 * ln of A (t) = Phi (t) + s (t) Lap Phi (t) / 2 at t = k->t = 1 - delta, the marginal to
 * second order in the distance of the sheared ellipticity to its mean. w = m_t (z) is the
 * TRACE_DET map of the variable of integration z, which has the mean t and the variance
 * s (t) per component,
 *   s (t) = (1 - t^2)^2 / 2 INT_0^1 (density of r) r^2 / (1 - t^2 r^2) dr,
 * and Phi (w) = N_2 (eps_obs - C (w)), C the identity for TRACE_DET and T for TRACE. The
 * integral is computed in v = 1 - r by ncm_spectral_polar_integrate(), where the factor
 * 1 / (delta + t v) of 1 / (1 - t^2 r^2) varies on the scale delta, with radial edges
 * placed for that scale at v = 0. Differences to 1 are formed from delta.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_second_order_ln (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k,
                                                  const gdouble std_noise, const gdouble a, const gdouble b, const gdouble delta,
                                                  NcGalaxyShapeFactorSpectralWork *work, GError **error)
{
  const gdouble t                    = k->t;
  const gdouble noise_var            = gsl_pow_2 (std_noise);
  const gdouble one_minus_t2         = delta * (1.0 + t);
  NcGalaxyShapeFactorSpectralVar var = {self, k, delta};
  gdouble s, d, lap_ratio;

  ncm_spectral_polar_peak_edges (self->polar, 1.0, 0.0, 0.0, delta, M_PI, work->r_edges, work->theta_edges);
  s = 0.5 * gsl_pow_2 (one_minus_t2) *
      ncm_spectral_polar_integrate (self->polar, &_nc_galaxy_shape_factor_spectral_var_kernel, &var, 1.0, b, a,
                                    work->r_edges, work->theta_edges, self->reltol, 0.0, NULL, NULL, error);

  if (self->axis_ratio_var)
  {
    /* C (t) = 2 t / (1 + t^2) = 1 - delta^2 / (1 + t^2) */
    const gdouble q     = 1.0 + t * t;
    const gdouble d_1   = (k->e_1 - 1.0) + delta * delta / q;
    const gdouble d_2   = k->e_2;
    const gdouble gr_1  = -4.0 * d_1 * one_minus_t2 / (q * q);
    const gdouble gr_2  = -4.0 * d_2 / q;
    const gdouble lap_d = 8.0 * gsl_pow_2 (one_minus_t2 / (q * q)) + 8.0 / (q * q) + 32.0 * t * d_1 / gsl_pow_3 (q);

    d         = d_1 * d_1 + d_2 * d_2;
    lap_ratio = -0.5 * lap_d / noise_var + 0.25 * (gr_1 * gr_1 + gr_2 * gr_2) / gsl_pow_2 (noise_var);
  }
  else
  {
    const gdouble d_1 = (k->e_1 - 1.0) + delta;
    const gdouble d_2 = k->e_2;

    d         = d_1 * d_1 + d_2 * d_2;
    lap_ratio = d / gsl_pow_2 (noise_var) - 2.0 / noise_var;
  }

  return -0.5 * d / noise_var - log (2.0 * M_PI * noise_var) + log1p (0.5 * s * lap_ratio);
}

/*
 * The integrand of the variance for ncm_spectral_polar_integrate() in v = 1 - r, with the
 * weight v^b (1 - v)^a: (reduced density of r) r^2 / ((delta + t v) (1 + t r)) / (2 pi).
 * It does not depend on the angle.
 */
static void
_nc_galaxy_shape_factor_spectral_var_kernel (gpointer user_data, gdouble v, guint n, const gdouble *theta,
                                             const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  const NcGalaxyShapeFactorSpectralVar *var = user_data;
  const gdouble t                           = var->k->t;
  const gdouble r                           = 1.0 - v;
  const gdouble val                         = var->self->p_reduced (var->k, r) * r * r / ((var->delta + t * v) * (1.0 + t * r)) / (2.0 * M_PI);
  guint j;

  for (j = 0; j < n; j++)
    fval[j] = val;
}

static gdouble _nc_galaxy_shape_factor_spectral_h (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k, const gdouble a, const gdouble b, const gdouble r, const gdouble cos_0, const gdouble sin_0, gdouble *d);
static void _nc_galaxy_shape_factor_spectral_sort (gdouble *x, const guint n);

/*
 * Locates the peak of the integrand on the ray through c, the preimage of eps_obs by the
 * map of the variable of integration, or of the point of the unit circle nearest to
 * eps_obs when |eps_obs| >= 1: the angle of the ray, the stretch factors of the inverse
 * map at that point by finite differences, the smallest d found and the radius r* that
 * maximizes h (r) = ln (density of r) - d / 2 std_noise^2, d the squared distance of
 * eps_obs to the sheared ellipticity. h is scanned
 * on a uniform grid of radii merged with radii next to both ends and at |c| plus and minus
 * multiples of the widths std_noise J_1 and std_noise J_2, and r* is refined by golden
 * section between the neighbors of the best radius.
 */
static void
_nc_galaxy_shape_factor_spectral_seed (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k,
                                       const gdouble std_noise, const gdouble a, const gdouble b,
                                       NcGalaxyShapeFactorSpectralSeed *seed)
{
  static const gdouble offsets[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_OFFSETS] = {0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 25.0};
  const complex double eps                                                   = k->e_1 + I * k->e_2;
  const gdouble abs_eps                                                      = cabs (eps);
  const complex double eps_0                                                 = (abs_eps < 1.0) ? eps : eps / abs_eps * (1.0 - 1.0e-12);
  const complex double c                                                     = self->map_inv (k->t, eps_0);
  const gdouble abs_c                                                        = cabs (c);
  const gdouble cos_0                                                        = (abs_c > 0.0) ? creal (c) / abs_c : 1.0;
  const gdouble sin_0                                                        = (abs_c > 0.0) ? cimag (c) / abs_c : 0.0;
  const gdouble dh                                                           = 1.0e-6;
  const guint n_uniform                                                      = NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM;
  gdouble extra[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  gdouble best_h   = GSL_NEGINF;
  gdouble r_before = 0.0;
  gdouble r_after  = 0.0;
  gdouble r_prev   = 0.0;
  gdouble shift    = GSL_POSINF;
  gboolean pending = FALSE;
  guint n_extra    = 0;
  guint iu         = 0;
  guint ie         = 0;
  guint i;

  seed->theta_0 = carg (c);

  if (cabs (eps_0) > 0.0)
  {
    seed->J_1 = cabs (self->map_inv (k->t, eps_0 * (1.0 - dh)) - c) / (dh * cabs (eps_0));
    seed->J_2 = cabs (self->map_inv (k->t, eps_0 * cexp (I * dh)) - c) / (dh * cabs (eps_0));
  }
  else
  {
    seed->J_1 = seed->J_2 = cabs (self->map_inv (k->t, dh) - c) / dh;
  }

  for (i = 2; i < 2 + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_END_DECADES; i++)
  {
    extra[n_extra++] = pow (10.0, -(gdouble) i);
    extra[n_extra++] = 1.0 - pow (10.0, -(gdouble) i);
  }

  for (i = 0; i < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_OFFSETS; i++)
  {
    extra[n_extra++] = abs_c + std_noise * seed->J_1 * offsets[i];
    extra[n_extra++] = abs_c - std_noise * seed->J_1 * offsets[i];
    extra[n_extra++] = abs_c + std_noise * seed->J_2 * offsets[i];
    extra[n_extra++] = abs_c - std_noise * seed->J_2 * offsets[i];
  }

  _nc_galaxy_shape_factor_spectral_sort (extra, n_extra);

  /* The uniform radii and the extra ones, in ascending order, each radius once */
  while ((iu < n_uniform) || (ie < n_extra))
  {
    const gdouble r_u = (iu < n_uniform) ? 1.0e-4 + (1.0 - 2.0e-4) * iu / (n_uniform - 1.0) : GSL_POSINF;
    const gdouble r_e = (ie < n_extra) ? extra[ie] : GSL_POSINF;
    const gdouble r   = GSL_MIN (r_u, r_e);
    gdouble h, d;

    if (r_u <= r_e)
      iu++;
    else
      ie++;

    if (!(r > 0.0) || !(r < 1.0) || (r == r_prev))
      continue;

    h     = _nc_galaxy_shape_factor_spectral_h (self, k, a, b, r, cos_0, sin_0, &d);
    shift = GSL_MIN (shift, d);

    if (pending)
    {
      r_after = r;
      pending = FALSE;
    }

    if (h > best_h)
    {
      best_h   = h;
      r_before = (r_prev > 0.0) ? r_prev : r;
      r_after  = r;
      pending  = TRUE;
    }

    r_prev = r;
  }

  /* The maximum of h on [r_before, r_after] by golden section */
  {
    const gdouble gr = 0.5 * (sqrt (5.0) - 1.0);
    gdouble lo       = r_before;
    gdouble hi       = r_after;
    gdouble x1       = hi - gr * (hi - lo);
    gdouble x2       = lo + gr * (hi - lo);
    gdouble d;
    gdouble h1 = _nc_galaxy_shape_factor_spectral_h (self, k, a, b, x1, cos_0, sin_0, &d);
    gdouble h2 = _nc_galaxy_shape_factor_spectral_h (self, k, a, b, x2, cos_0, sin_0, &d);

    for (i = 0; i < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_GOLDEN_STEPS; i++)
    {
      if (h1 > h2)
      {
        hi = x2;
        x2 = x1;
        h2 = h1;
        x1 = hi - gr * (hi - lo);
        h1 = _nc_galaxy_shape_factor_spectral_h (self, k, a, b, x1, cos_0, sin_0, &d);
      }
      else
      {
        lo = x1;
        x1 = x2;
        h1 = h2;
        x2 = lo + gr * (hi - lo);
        h2 = _nc_galaxy_shape_factor_spectral_h (self, k, a, b, x2, cos_0, sin_0, &d);
      }
    }

    seed->r_star = 0.5 * (lo + hi);
  }

  seed->shift = shift;
}

/*
 * h (r) = ln (density of r) - d / (2 std_noise^2) at the radius r, 0 < r < 1, of the ray of
 * direction (cos_0, sin_0), with d the squared distance of eps_obs to the sheared
 * ellipticity, returned in *d. The logarithm of the density is the two end power laws,
 * of exponents a and b, plus the logarithm of the reduced density.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_h (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k,
                                    const gdouble a, const gdouble b, const gdouble r, const gdouble cos_0, const gdouble sin_0,
                                    gdouble *d)
{
  gdouble f_x, f_y;

  self->map (k->t, r * cos_0, r * sin_0, &f_x, &f_y);
  *d = gsl_pow_2 (k->e_1 - f_x) + gsl_pow_2 (k->e_2 - f_y);

  return a * log (r) + b * log1p (-r) + log (self->p_reduced (k, r)) - *d * k->c2;
}

/* Sorts x[0..n-1] in ascending order, by insertion */
static void
_nc_galaxy_shape_factor_spectral_sort (gdouble *x, const guint n)
{
  guint i;

  for (i = 1; i < n; i++)
  {
    const gdouble v = x[i];
    guint j         = i;

    while ((j > 0) && (x[j - 1] > v))
    {
      x[j] = x[j - 1];
      j--;
    }

    x[j] = v;
  }
}

/*
 * The integrand of ncm_spectral_polar_integrate() for the TRACE convention, in the polar
 * coordinates of e_I = T^{-1} (chi_I): G (r, theta) = (reduced density of r = |e_I|)
 * N_2 (eps_obs - T (m_t (e_I))) / (2 pi), with m_t the TRACE_DET map and the offset of the
 * kernel taken out of the exponent of N_2.
 */
static void
_nc_galaxy_shape_factor_spectral_kernel_trace (gpointer user_data, gdouble r, guint n, const gdouble *theta,
                                               const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  const NcGalaxyShapeFactorSpectralKernel *k = user_data;
  const gdouble pref                         = _nc_galaxy_shape_factor_spectral_p_reduced_trace (k, r) * k->norm;
  guint j;

  for (j = 0; j < n; j++)
  {
    gdouble f_x, f_y;

    _nc_galaxy_shape_factor_spectral_map_trace (k->t, r * cos_theta[j], r * sin_theta[j], &f_x, &f_y);
    fval[j] = pref * exp (k->shift_c2 - (gsl_pow_2 (k->e_1 - f_x) + gsl_pow_2 (k->e_2 - f_y)) * k->c2);
  }
}

/*
 * The integrand for the TRACE_DET convention, in the polar coordinates of chi_I:
 * G (r, theta) = (P_pop (r) / (r^a (1 - r)^b)) N_2 (eps_obs - f_t (chi_I)) / (2 pi).
 */
static void
_nc_galaxy_shape_factor_spectral_kernel_trace_det (gpointer user_data, gdouble r, guint n, const gdouble *theta,
                                                   const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  const NcGalaxyShapeFactorSpectralKernel *k = user_data;
  const gdouble pref                         = _nc_galaxy_shape_factor_spectral_p_reduced_trace_det (k, r) * k->norm;
  guint j;

  for (j = 0; j < n; j++)
  {
    gdouble f_x, f_y;

    nc_wl_ellipticity_apply_shear_trace_det_real (k->t, r * cos_theta[j], r * sin_theta[j], &f_x, &f_y);
    fval[j] = pref * exp (k->shift_c2 - (gsl_pow_2 (k->e_1 - f_x) + gsl_pow_2 (k->e_2 - f_y)) * k->c2);
  }
}

/*
 * The TRACE map at the real shear t in the variable e_I = x + i y: T (m_t (e_I)), with
 * m_t the TRACE_DET map. It equals f_t (T (e_I)), f_t the TRACE map.
 */
static void
_nc_galaxy_shape_factor_spectral_map_trace (gdouble t, gdouble x, gdouble y, gdouble * restrict f_x, gdouble * restrict f_y)
{
  gdouble m_x, m_y;

  nc_wl_ellipticity_apply_shear_trace_det_real (t, x, y, &m_x, &m_y);
  {
    const gdouble T_fac = 2.0 / (1.0 + m_x * m_x + m_y * m_y);

    *f_x = T_fac * m_x;
    *f_y = T_fac * m_y;
  }
}

/* The inverse of the map above at the shear g: m_g^{-1} (T^{-1} (e_obs)), for |e_obs| <= 1 */
static complex double
_nc_galaxy_shape_factor_spectral_map_inv_trace (complex double g, complex double e_obs)
{
  const gdouble abs_e2 = gsl_pow_2 (creal (e_obs)) + gsl_pow_2 (cimag (e_obs));

  return nc_wl_ellipticity_apply_shear_inv_trace_det (g, e_obs / (1.0 + sqrt (1.0 - abs_e2)));
}

/*
 * The density of r = |e_I| divided by r^alpha_o (1 - r)^(2 alpha_1 + 1), for a population
 * that is a density of |chi_I| = 2 r / (1 + r^2) with the exponents alpha_o and alpha_1:
 * 2^(alpha_o + 1) (1 + r) (1 + r^2)^-(alpha_o + alpha_1 + 2) times the reduced density of
 * the population at |chi_I|.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_p_reduced_trace (const NcGalaxyShapeFactorSpectralKernel *k, const gdouble r)
{
  const gdouble q = 1.0 + r * r;

  return k->dens_norm * (1.0 + r) * pow (q, -k->dens_pow) * nc_galaxy_shape_pop_eval_p_reduced (k->pop, k->pop_data, 2.0 * r / q);
}

/* The reduced density of the population at r = |chi_I| */
static gdouble
_nc_galaxy_shape_factor_spectral_p_reduced_trace_det (const NcGalaxyShapeFactorSpectralKernel *k, const gdouble r)
{
  return nc_galaxy_shape_pop_eval_p_reduced (k->pop, k->pop_data, r);
}

/**
 * nc_galaxy_shape_factor_spectral_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorSpectral with the ellipticity convention @ellip_conv.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorSpectral.
 */
NcGalaxyShapeFactorSpectral *
nc_galaxy_shape_factor_spectral_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_SPECTRAL,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_spectral_ref:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Increases the reference count of @gsfs by one.
 *
 * Returns: (transfer full): @gsfs.
 */
NcGalaxyShapeFactorSpectral *
nc_galaxy_shape_factor_spectral_ref (NcGalaxyShapeFactorSpectral *gsfs)
{
  return g_object_ref (gsfs);
}

/**
 * nc_galaxy_shape_factor_spectral_free:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Decreases the reference count of @gsfs by one.
 */
void
nc_galaxy_shape_factor_spectral_free (NcGalaxyShapeFactorSpectral *gsfs)
{
  g_object_unref (gsfs);
}

/**
 * nc_galaxy_shape_factor_spectral_clear:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Decreases the reference count of *@gsfs by one and sets *@gsfs to %NULL.
 */
void
nc_galaxy_shape_factor_spectral_clear (NcGalaxyShapeFactorSpectral **gsfs)
{
  g_clear_object (gsfs);
}

/**
 * nc_galaxy_shape_factor_spectral_set_reltol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @reltol: relative tolerance of the integral
 *
 * Sets #NcGalaxyShapeFactorSpectral:reltol.
 */
void
nc_galaxy_shape_factor_spectral_set_reltol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble reltol)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->reltol = reltol;
}

/**
 * nc_galaxy_shape_factor_spectral_get_reltol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:reltol.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_reltol (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->reltol;
}

/**
 * nc_galaxy_shape_factor_spectral_set_use_interp:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @use_interp: whether $\ln P$ is interpolated in $g$ for each galaxy
 *
 * Sets #NcGalaxyShapeFactorSpectral:use-interp.
 */
void
nc_galaxy_shape_factor_spectral_set_use_interp (NcGalaxyShapeFactorSpectral *gsfs, const gboolean use_interp)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->use_interp = use_interp;
}

/**
 * nc_galaxy_shape_factor_spectral_get_use_interp:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:use-interp.
 */
gboolean
nc_galaxy_shape_factor_spectral_get_use_interp (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->use_interp;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_g_max:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @g_max: largest $|g|$ computed from the interpolant, in $(0, 1/3)$
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-g-max.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_g_max (NcGalaxyShapeFactorSpectral *gsfs, const gdouble g_max)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  g_assert_cmpfloat (g_max, >, 0.0);
  g_assert_cmpfloat (3.0 * g_max, <, 1.0);

  self->interp_g_max = g_max;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_g_max:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-g-max.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_interp_g_max (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_g_max;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_abstol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @abstol: absolute tolerance on $\ln P$ of the interpolant
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-abstol.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_abstol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble abstol)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  g_assert_cmpfloat (abstol, >, 0.0);

  self->interp_abstol = abstol;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_abstol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-abstol.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_interp_abstol (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_abstol;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @n: minimum number of radial intervals of the interpolant
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-min-radial-panels, which must be at least 1.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (NcGalaxyShapeFactorSpectral *gsfs, const guint n)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  ncm_spectral_polar_set_min_radial_panels (self->polar_interp, n);

  self->interp_min_radial_panels = n;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-min-radial-panels.
 */
guint
nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_min_radial_panels;
}

/**
 * nc_galaxy_shape_factor_spectral_eval_ln_direct:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @pop: a #NcGalaxyShapePop
 * @pop_data: a resolved #NcGalaxyShapePopData of @pop
 * @g_1: reduced shear, real component
 * @g_2: reduced shear, imaginary component
 * @epsilon_obs_1: observed ellipticity, real component
 * @epsilon_obs_2: observed ellipticity, imaginary component
 * @std_noise: measurement noise standard deviation
 *
 * Computes $\ln P(\epsilon_\mathrm{obs} \mid g)$ from the population and the values of
 * one galaxy, without a #NcGalaxyShapeFactorData. @pop_data must have been prepared with
 * nc_galaxy_shape_pop_prepare() since the parameters of @pop last changed.
 * nc_galaxy_shape_factor_eval_ln_marginal() computes the same value for the galaxy of a
 * #NcGalaxyShapeFactorData.
 *
 * Returns: $\ln P(\epsilon_\mathrm{obs} \mid g)$.
 */
gdouble
nc_galaxy_shape_factor_spectral_eval_ln_direct (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                                const gdouble g_1, const gdouble g_2,
                                                const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                                const gdouble std_noise)
{
  return _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, pop_data, std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, NULL);
}

/**
 * nc_galaxy_shape_factor_spectral_eval_panels:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @pop: a #NcGalaxyShapePop
 * @pop_data: a resolved #NcGalaxyShapePopData of @pop
 * @g_1: reduced shear, real component
 * @g_2: reduced shear, imaginary component
 * @epsilon_obs_1: observed ellipticity, real component
 * @epsilon_obs_2: observed ellipticity, imaginary component
 * @std_noise: measurement noise standard deviation
 * @ln_marginal: (out) (optional): $\ln P(\epsilon_\mathrm{obs} \mid g)$
 *
 * Computes $\ln P(\epsilon_\mathrm{obs} \mid g)$ as
 * nc_galaxy_shape_factor_spectral_eval_ln_direct() and returns the panels the integral
 * was accepted on, for inspection: one row per panel with the columns of
 * ncm_spectral_polar_integrate(), $r_a$, $r_b$, $\theta_a$, $\theta_b$, the two levels,
 * the value and the change of the panel. $(r, \theta)$ are the polar coordinates of the
 * variable of integration in the frame of @g_1, @g_2 and of the observed ellipticity:
 * $\chi_I$ for %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, and for
 * %NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE the variable with $|\chi_I| = 2 r / (1 + r^2)$ and
 * the angle of $\chi_I$. For $1 - \min(|g|, 1 / |g|) < 10^{-5}$ the panels are those of
 * the integral at $1 - 10^{-5}$. For $|g| > 1$ the
 * panels are those of the integral at $g / |g|^2$, and at $|g| = 1$ there is no integral
 * and the function returns %NULL.
 *
 * Returns: (transfer full) (nullable): the accepted panels.
 */
NcmMatrix *
nc_galaxy_shape_factor_spectral_eval_panels (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                             const gdouble g_1, const gdouble g_2,
                                             const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                             const gdouble std_noise, gdouble *ln_marginal)
{
  NcmMatrix *panels  = NULL;
  const gdouble ln_P = _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, pop_data, std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, &panels);

  if (ln_marginal != NULL)
    *ln_marginal = ln_P;

  return panels;
}

