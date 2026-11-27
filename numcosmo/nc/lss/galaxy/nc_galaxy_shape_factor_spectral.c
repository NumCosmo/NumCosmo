/***************************************************************************
 *            nc_galaxy_shape_factor_spectral.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_spectral.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * A galaxy may hold an interpolant: an approximation of $\ln P$ in $g$, an
 * #NcmSpectralTensor built by #NcmSpectralTensorBuilder from the integral, to the tolerance
 * #NcGalaxyShapeFactorSpectral:interp-abstol, kept in the #NcGalaxyShapeFactorData of the
 * galaxy. nc_galaxy_shape_factor_eval_marginal() and
 * nc_galaxy_shape_factor_eval_ln_marginal() take $\ln P$ from it for every shear in its
 * range, $|g| \le$ #NcGalaxyShapeFactorSpectral:interp-g-max or less with
 * #NcGalaxyShapeFactorSpectral:restrict-range, and from the integral otherwise. The
 * interpolant is built only by nc_galaxy_shape_factor_spectral_build_interp(), which
 * nc_data_cluster_wl_factor_write_shape_interp() calls for every galaxy of a catalog, or
 * read from a catalog that carries one; the evaluations never build it. Each free parameter
 * of the population is one more variable of the approximation, over the bounds of the
 * parameter (#NcGalaxyShapeFactorSpectral:interp-pop-level-min,
 * #NcGalaxyShapeFactorSpectral:interp-pop-level-cap): a change of a free parameter inside
 * its bounds uses the same approximation, and a value outside them aborts. With
 * #NcGalaxyShapeFactorSpectral:interp-pop-one-integral the populations of a call of the
 * builder share one integral per shear.
 *
 * The interpolants can be stored in the catalog: nc_galaxy_shape_factor_data_write_row()
 * writes them (#NcGalaxyWLObs:shape-interp), or
 * nc_data_cluster_wl_factor_write_shape_interp_file() writes them to a file the catalog
 * refers to (#NcGalaxyWLObsShapeInterp), and
 * nc_galaxy_shape_factor_spectral_set_obs_source() records the shape factor and the
 * population that built them. Reading the rows of such a catalog loads the interpolants,
 * and an interpolant that is present is always used. The first use compares the record with
 * the fit (the convention, reltol, the interp-* settings, the population parameters that
 * are not variables of the stored interpolants, which must be fixed at the recorded values
 * within a relative tolerance of $10^{-15}$, and the free ones, which must be variables with
 * the same bounds; a variable fixed in the fit holds any value in its range) and aborts on
 * any difference. An
 * interpolant, loaded or built, that no longer holds for an evaluation (another observed
 * ellipticity, noise, setting or population) aborts. A change of the rows of the catalog, as
 * a resample, drops its interpolants.
 *
 * A call takes its scratch space from a pool, so one object may be used for several
 * galaxies in several threads. Two threads must not evaluate the same galaxy that holds an
 * interpolant, since an evaluation records the state of the population in it.
 *
 * See <a href="../../theory/nc/lss/galaxy/wl_shape_marginalization_spectral.html">Shape
 * Marginalization by Spectral Panels</a> for the formulas.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_spectral.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_moments_private.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/algebra/ncm_spectral_polar.h"
#include "ncm/algebra/ncm_spectral_tensor.h"
#include "ncm/algebra/ncm_spectral_tensor_builder.h"
#include "ncm/algebra/ncm_spectral_polar_private.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_vmath_private.h"

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
 * The functions of the convention act on the variable of integration z: chi_I for
 * TRACE_DET and T^{-1} (chi_I) for TRACE, T (e) = 2 e / (1 + |e|^2). Those of the integrals
 * take the convention from NcGalaxyShapeFactorSpectralKernel:trace, once per call.
 */
typedef struct _NcGalaxyShapeFactorSpectralPrivate
{
  /* the sheared ellipticity e_obs to z, at the shear g */
  complex double (*map_inv) (complex double g, complex double e_obs);

  gboolean axis_ratio_var; /* TRUE for TRACE: z = T^{-1} (chi_I) */

  NcmSpectralPolar *polar;                  /* the integrals */
  NcmSpectralTensorBuilder *builder_interp; /* the build of the interpolant in g */
  NcmMemoryPool *work_pool;                 /* NcGalaxyShapeFactorSpectralWork, one per concurrent evaluation */
  NcmMemoryPool *line_pool;                 /* NcGalaxyShapeFactorSpectralLineWork, one per concurrent eval_line */
  gdouble line_abstol;
  gdouble reltol;
  gdouble interp_g_max;
  gdouble interp_abstol;
  guint interp_min_radial_panels;
  gdouble interp_trim_frac;
  guint interp_pop_level_min;
  guint interp_pop_level_cap;
  gboolean interp_pop_one_integral;
  gdouble interp_compress_frac;

  /*
   * restrict-range: whether each interpolant covers |g| only up to the bound of its galaxy
   * (the bound fields of #NcGalaxyShapeFactorData)
   */
  gboolean restrict_range;
} NcGalaxyShapeFactorSpectralPrivate;

/* The scratch space of one evaluation */
typedef struct _NcGalaxyShapeFactorSpectralWork
{
  GArray *r_edges;
  GArray *theta_edges;
  GArray *r_edges_all; /* radial edges of the interpolant before those below its range are left out */
  GArray *var_r;       /* the radii of a call of the radial factor of the variance */
  GArray *g_arr;       /* the shears of a call of eval_nodes that the interpolant takes, two values per shear */
  GArray *x_arr;       /* their points of the interpolant without the population variables, two values per point */
  GArray *x_pop;       /* the population variables of those points, the same for all of them */
  GArray *f_arr;       /* their values */
  GArray *j_arr;       /* guint, the index of each of those points among the shears of the call */
  GArray *lower;       /* the domain of the interpolant, one value per variable */
  GArray *upper;
  GArray *level_min; /* guint, the levels of its build, one value per variable */
  GArray *level_cap;
  GArray *series;        /* NcmSpectralPolarRadialSeries, the radial series of a reference integral */
  GArray *series_coeffs; /* and their coefficients */
} NcGalaxyShapeFactorSpectralWork;

/* The scratch space of one eval_line: the 1D expansion and its coefficients */
typedef struct _NcGalaxyShapeFactorSpectralLineWork
{
  NcmSpectral *spectral;
  GArray *coeffs;
} NcGalaxyShapeFactorSpectralLineWork;

/*
 * Relative tolerance of the comparison of a loaded interpolant's population parameters and
 * settings with those of the fit: a few ulps, which a decimal round trip of the record can
 * change; a parameter difference this small changes ln P far below interp-abstol.
 */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SOURCE_RELTOL (1.0e-15)

/* Ratio between consecutive radial edge distances of the interpolant in g */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INTERP_EDGE_RATIO (4.0)

/*
 * The table of ln P along the line of the shears of a galaxy (eval_line): the first grid has
 * 2^LINE_LEVEL_MIN + 1 points and the finest at most 2^LINE_LEVEL_MAX + 1, and every shear of
 * the line has |g| <= LINE_G_MAX.
 */
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MIN (1)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MAX (7)
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_G_MAX (0.9)

/*
 * The interpolant of ln P in g of one galaxy and what it was built for. It is rebuilt when
 * one of these values differs from the one of the evaluation, unless it was loaded from the
 * catalog (loaded): then src_gsfs and src_pop are the record that built it, which
 * data_prepare() compares with the fit, and a loaded interpolant is never rebuilt. With free
 * population parameters, st has one variable for each after (rho, psi): pop_pkey is then the
 * state of pop last found inside the table, and a new state is compared with pop_idx,
 * pop_bounds and the fixed entries of pop_params.
 */
typedef struct _NcGalaxyShapeFactorSpectralData
{
  NcmSpectralTensor *st;
  gboolean valid;
  gboolean loaded;
  NcGalaxyShapeFactorSpectral *src_gsfs;
  NcGalaxyShapePop *src_pop;
  gdouble epsilon_obs_1;
  gdouble epsilon_obs_2;
  gdouble std_noise;
  guint64 pop_pkey;
  GArray *pop_idx;    /* guint, the free parameters of the population, the variables 2, 3, ... of st */
  GArray *pop_bounds; /* gdouble, the lower and upper bound of each, the range of its variable */
  GArray *pop_params; /* gdouble, every parameter of the population; the free ones are not compared */
  GArray *x;          /* gdouble, a point of st: rho, psi and the free parameters */
  guint pop_level_min;
  guint pop_level_cap;
  gboolean pop_one_integral;
  gdouble compress_frac;
  gdouble g_max;
  gdouble abstol;
  gdouble reltol;
  guint min_radial_panels;
  gdouble trim_frac;
  gdouble g_target; /* the |g| range data_prepare() asks of the interpolant; 0: interp-g-max */
  gdouble cos_e;    /* direction of eps_obs */
  gdouble sin_e;
  gdouble g_c; /* center of the coordinates of the interpolant */
} NcGalaxyShapeFactorSpectralData;

/*
 * The functions that a loop of the integrals calls with the convention as a constant: inlined
 * always, so that each copy of the loop holds a single map and vectorizes
 */
#ifdef __GNUC__
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE static inline __attribute__ ((always_inline))
#else
#define NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE static inline
#endif

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
  gboolean trace;    /* the convention: TRUE for TRACE (NcGalaxyShapeFactorSpectralPrivate:axis_ratio_var) */
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
  PROP_INTERP_G_MAX,
  PROP_INTERP_ABSTOL,
  PROP_INTERP_MIN_RADIAL_PANELS,
  PROP_INTERP_TRIM_FRAC,
  PROP_INTERP_POP_LEVEL_MIN,
  PROP_INTERP_POP_LEVEL_CAP,
  PROP_INTERP_POP_ONE_INTEGRAL,
  PROP_INTERP_COMPRESS_FRAC,
  PROP_RESTRICT_RANGE,
  PROP_LINE_ABSTOL,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorSpectral, nc_galaxy_shape_factor_spectral, NC_TYPE_GALAXY_SHAPE_FACTOR)

static gpointer _nc_galaxy_shape_factor_spectral_work_new (gpointer userdata);
static void _nc_galaxy_shape_factor_spectral_work_free (gpointer p);
static gpointer _nc_galaxy_shape_factor_spectral_line_work_new (gpointer userdata);
static void _nc_galaxy_shape_factor_spectral_line_work_free (gpointer p);

static void
nc_galaxy_shape_factor_spectral_init (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->map_inv        = NULL;
  self->axis_ratio_var = FALSE;
  self->polar          = ncm_spectral_polar_new ();
  self->builder_interp = ncm_spectral_tensor_builder_new ();
  self->work_pool      = ncm_memory_pool_new (_nc_galaxy_shape_factor_spectral_work_new, NULL,
                                              _nc_galaxy_shape_factor_spectral_work_free);
  self->reltol        = 0.0;
  self->interp_g_max  = 0.0;
  self->interp_abstol = 0.0;

  self->interp_min_radial_panels = 0;
  self->interp_trim_frac         = 0.0;
  self->interp_pop_level_min     = 0;
  self->interp_pop_level_cap     = 0;
  self->interp_pop_one_integral  = FALSE;
  self->interp_compress_frac     = 0.0;
  self->restrict_range           = FALSE;

  self->line_pool = ncm_memory_pool_new (_nc_galaxy_shape_factor_spectral_line_work_new, NULL,
                                         _nc_galaxy_shape_factor_spectral_line_work_free);
  self->line_abstol = 0.0;
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
    case PROP_INTERP_G_MAX:
      nc_galaxy_shape_factor_spectral_set_interp_g_max (gsfs, g_value_get_double (value));
      break;
    case PROP_INTERP_ABSTOL:
      nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, g_value_get_double (value));
      break;
    case PROP_INTERP_MIN_RADIAL_PANELS:
      nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (gsfs, g_value_get_uint (value));
      break;
    case PROP_INTERP_TRIM_FRAC:
      nc_galaxy_shape_factor_spectral_set_interp_trim_frac (gsfs, g_value_get_double (value));
      break;
    case PROP_INTERP_POP_LEVEL_MIN:
      nc_galaxy_shape_factor_spectral_set_interp_pop_level_min (gsfs, g_value_get_uint (value));
      break;
    case PROP_INTERP_POP_LEVEL_CAP:
      nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs, g_value_get_uint (value));
      break;
    case PROP_INTERP_POP_ONE_INTEGRAL:
      nc_galaxy_shape_factor_spectral_set_interp_pop_one_integral (gsfs, g_value_get_boolean (value));
      break;
    case PROP_INTERP_COMPRESS_FRAC:
      nc_galaxy_shape_factor_spectral_set_interp_compress_frac (gsfs, g_value_get_double (value));
      break;
    case PROP_RESTRICT_RANGE:
      nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs, g_value_get_boolean (value));
      break;
    case PROP_LINE_ABSTOL:
      nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, g_value_get_double (value));
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
    case PROP_INTERP_G_MAX:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_interp_g_max (gsfs));
      break;
    case PROP_INTERP_ABSTOL:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_interp_abstol (gsfs));
      break;
    case PROP_INTERP_MIN_RADIAL_PANELS:
      g_value_set_uint (value, nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (gsfs));
      break;
    case PROP_INTERP_TRIM_FRAC:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_interp_trim_frac (gsfs));
      break;
    case PROP_INTERP_POP_LEVEL_MIN:
      g_value_set_uint (value, nc_galaxy_shape_factor_spectral_get_interp_pop_level_min (gsfs));
      break;
    case PROP_INTERP_POP_LEVEL_CAP:
      g_value_set_uint (value, nc_galaxy_shape_factor_spectral_get_interp_pop_level_cap (gsfs));
      break;
    case PROP_INTERP_POP_ONE_INTEGRAL:
      g_value_set_boolean (value, nc_galaxy_shape_factor_spectral_get_interp_pop_one_integral (gsfs));
      break;
    case PROP_INTERP_COMPRESS_FRAC:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_interp_compress_frac (gsfs));
      break;
    case PROP_RESTRICT_RANGE:
      g_value_set_boolean (value, nc_galaxy_shape_factor_spectral_get_restrict_range (gsfs));
      break;
    case PROP_LINE_ABSTOL:
      g_value_set_double (value, nc_galaxy_shape_factor_spectral_get_line_abstol (gsfs));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static complex double _nc_galaxy_shape_factor_spectral_map_inv_trace (complex double g, complex double e_obs);

static void _nc_galaxy_shape_factor_spectral_map_trace (gdouble t, gdouble x, gdouble y, gdouble * restrict f_x, gdouble * restrict f_y);
static void _nc_galaxy_shape_factor_spectral_kernel (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval);
static void _nc_galaxy_shape_factor_spectral_radial (gpointer user_data, guint n, const gdouble *r, gdouble *rho);
static void _nc_galaxy_shape_factor_spectral_density (const NcGalaxyShapeFactorSpectralKernel *k, guint n, const gdouble *r, gdouble *out);
static void _nc_galaxy_shape_factor_spectral_dist2 (const NcGalaxyShapeFactorSpectralKernel *k, guint n, const gdouble *r, gdouble cos_0, gdouble sin_0, gdouble *d);

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
        self->map_inv        = &_nc_galaxy_shape_factor_spectral_map_inv_trace;
        self->axis_ratio_var = TRUE;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->map_inv        = &nc_wl_ellipticity_apply_shear_inv_trace_det;
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
  ncm_spectral_tensor_builder_clear (&self->builder_interp);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_spectral_parent_class)->dispose (object);
}

static void
_nc_galaxy_shape_factor_spectral_finalize (GObject *object)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (object);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  ncm_memory_pool_free (self->work_pool, TRUE);
  ncm_memory_pool_free (self->line_pool, TRUE);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_spectral_parent_class)->finalize (object);
}

static void _nc_galaxy_shape_factor_spectral_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data);
static void _nc_galaxy_shape_factor_spectral_data_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max);
static gdouble _nc_galaxy_shape_factor_spectral_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static gboolean _nc_galaxy_shape_factor_spectral_eval_line (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n, const gdouble *t, const gdouble u_1, const gdouble u_2, const gdouble b_1, const gdouble b_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gboolean want_log, gdouble *out);
static void _nc_galaxy_shape_factor_spectral_eval_nodes (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n, const gdouble *g_1, const gdouble *g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gboolean want_log, gdouble *out);
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
   * ncm_spectral_tensor_builder_build(): a node becomes a panel node when $\ln P$ changes
   * by at most this value between its two grids.
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
   * The number $n \ge 1$ of equal intervals of the radial range of the interpolant: their
   * inner edges are radial edges of the interpolant, so the range is split into at least
   * $n$ intervals. More
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

  /**
   * NcGalaxyShapeFactorSpectral:interp-trim-frac:
   *
   * The #NcmSpectralTensorBuilder:trim-frac of the build of the interpolant: the fraction $f$
   * of #NcGalaxyShapeFactorSpectral:interp-abstol that the stored coefficients may give away.
   * The removed coefficients change $\ln P$ by at most $f$ times
   * #NcGalaxyShapeFactorSpectral:interp-abstol, $10^{-9}$ at the defaults. With 0 every
   * coefficient is stored.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_TRIM_FRAC,
                                   g_param_spec_double ("interp-trim-frac",
                                                        NULL,
                                                        "Fraction of interp-abstol given to the removed coefficients",
                                                        0.0, 1.0, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-pop-level-min:
   *
   * The starting level $l$ of the variables of the interpolant that are free parameters of the
   * population: the first grid has $2^l + 1$ points in each of them.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_POP_LEVEL_MIN,
                                   g_param_spec_uint ("interp-pop-level-min",
                                                      NULL,
                                                      "Starting level of the population variables of the interpolant",
                                                      1, 10, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-pop-level-cap:
   *
   * The largest level $l$ of the variables of the interpolant that are free parameters of the
   * population: a panel has at most $2^l + 1$ points in each of them, and a panel that needs
   * more is split. A build with a value below
   * #NcGalaxyShapeFactorSpectral:interp-pop-level-min aborts.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_POP_LEVEL_CAP,
                                   g_param_spec_uint ("interp-pop-level-cap",
                                                      NULL,
                                                      "Largest level of the population variables of the interpolant",
                                                      1, 10, 5,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-pop-one-integral:
   *
   * Whether the build of an interpolant with population variables computes, at each shear,
   * one integral for all the populations of a call of its builder. The integral is the one
   * of the population at the center of those populations, and the values of the others
   * follow from the radial series of its accepted panels and the factors of each
   * population, with no further value of the noise kernel. It applies to a population whose
   * reduced density has the separated form of nc_galaxy_shape_pop_reduced_factor_params(),
   * with each free parameter the parameter of one of the two factors. Otherwise, and when
   * this is %FALSE, each population takes its own integral. The reference integral resolves
   * the integrand of the population at the center, so the bounds of the free parameters
   * should contain only populations whose densities do not differ by many orders of
   * magnitude.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_POP_ONE_INTEGRAL,
                                   g_param_spec_boolean ("interp-pop-one-integral",
                                                         NULL,
                                                         "Whether the populations of a call share one integral per shear",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:interp-compress-frac:
   *
   * The fraction $f$ of #NcGalaxyShapeFactorSpectral:interp-abstol given to the compression
   * of an interpolant with population variables: after the build,
   * ncm_spectral_tensor_compress() with tolerance $f$ times
   * #NcGalaxyShapeFactorSpectral:interp-abstol stores each panel in the format
   * %NCM_SPECTRAL_TENSOR_PANEL_FORMAT_TUCKER when it holds fewer numbers. The compression
   * changes $\ln P$ by at most that tolerance, $10^{-7}$ at the defaults; the largest change
   * measured is about 0.15 of it. An interpolant evaluated at one point of the
   * population per call reduces such a panel once per call
   * (ncm_spectral_tensor_eval_slice_array()). With 0 no panel is compressed; an interpolant
   * without population variables is never compressed.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP_COMPRESS_FRAC,
                                   g_param_spec_double ("interp-compress-frac",
                                                        NULL,
                                                        "Fraction of interp-abstol given to the compression of the interpolant",
                                                        0.0, 1.0, 0.01,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:restrict-range:
   *
   * Whether the interpolant of a galaxy covers only the shears it can reach over the box of
   * the fit. nc_galaxy_shape_factor_spectral_build_interp() bounds $|g|$ at the corner of the
   * box (the smallest projected radius the box of the halo
   * center allows, the upper bound of log10MDelta, the largest source redshift of the galaxy),
   * $|1 + m| |g_t| + |c|$, and the interpolant covers $|g|$ up to the smaller of that bound and
   * #NcGalaxyShapeFactorSpectral:interp-g-max. The bound holds when only log10MDelta, the halo
   * center and the parameters of the shape population are free; otherwise, and when the models
   * of the bound are not in the fit, the interpolant covers interp-g-max. A shear above the
   * range of the interpolant takes the integral. An interpolant loaded from a catalog keeps the
   * range it was built with.
   */
  g_object_class_install_property (object_class,
                                   PROP_RESTRICT_RANGE,
                                   g_param_spec_boolean ("restrict-range",
                                                         NULL,
                                                         "Whether the interpolant covers only the shears a galaxy reaches",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorSpectral:line-abstol:
   *
   * The absolute tolerance on $\ln P$ of the table along the line of the shears of the
   * redshift nodes of a galaxy, used for a galaxy that holds no interpolant. With
   * a positive value, nc_galaxy_shape_factor_eval_shear_line() takes $\ln P$ at those shears
   * from a Chebyshev series in the position on the line, expanded from integrals on nested
   * Chebyshev-Lobatto grids and accepted when its coefficients change by at most this value
   * in $\ell_1$ norm between two grids. The table is made only when it needs at most half as
   * many integrals as there are nodes and every shear of the line has $|g| \le 0.9$; when it
   * is not made, or does not reach the tolerance, each node is integrated. With 0 (the
   * default) each node is integrated.
   */
  g_object_class_install_property (object_class,
                                   PROP_LINE_ABSTOL,
                                   g_param_spec_double ("line-abstol",
                                                        NULL,
                                                        "Absolute tolerance on ln P of the table along the line of the shears",
                                                        0.0, 1.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_spectral_data_init;
  gsf_class->data_prepare     = &_nc_galaxy_shape_factor_spectral_data_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_spectral_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_spectral_eval_ln_marginal;
  gsf_class->eval_nodes       = &_nc_galaxy_shape_factor_spectral_eval_nodes;
  gsf_class->eval_line        = &_nc_galaxy_shape_factor_spectral_eval_line;
}

static gpointer
_nc_galaxy_shape_factor_spectral_work_new (gpointer userdata)
{
  NcGalaxyShapeFactorSpectralWork *work = g_new0 (NcGalaxyShapeFactorSpectralWork, 1);

  work->r_edges     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->theta_edges = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->r_edges_all = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->var_r       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->g_arr       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->x_arr       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->x_pop       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->f_arr       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->j_arr       = g_array_new (FALSE, FALSE, sizeof (guint));
  work->lower       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->upper       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->level_min   = g_array_new (FALSE, FALSE, sizeof (guint));
  work->level_cap   = g_array_new (FALSE, FALSE, sizeof (guint));

  work->series        = g_array_new (FALSE, FALSE, sizeof (NcmSpectralPolarRadialSeries));
  work->series_coeffs = g_array_new (FALSE, FALSE, sizeof (gdouble));

  return work;
}

static void
_nc_galaxy_shape_factor_spectral_work_free (gpointer p)
{
  NcGalaxyShapeFactorSpectralWork *work = p;

  g_array_unref (work->r_edges);
  g_array_unref (work->theta_edges);
  g_array_unref (work->r_edges_all);
  g_array_unref (work->var_r);
  g_array_unref (work->g_arr);
  g_array_unref (work->x_arr);
  g_array_unref (work->x_pop);
  g_array_unref (work->f_arr);
  g_array_unref (work->j_arr);
  g_array_unref (work->lower);
  g_array_unref (work->upper);
  g_array_unref (work->level_min);
  g_array_unref (work->level_cap);
  g_array_unref (work->series);
  g_array_unref (work->series_coeffs);

  g_free (work);
}

/*
 * Takes the interpolant of row i of obs, when the catalog has one for its current rows, as a
 * loaded interpolant (shared with obs, never rebuilt in place). Otherwise the galaxy has no
 * valid interpolant.
 */
static void
_nc_galaxy_shape_factor_spectral_ldata_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorSpectralData *ldata = data->ldata;
  NcmSpectralTensor *st                  = nc_galaxy_wl_obs_peek_shape_interp (obs, i);

  if (ldata->loaded)
  {
    ncm_spectral_tensor_clear (&ldata->st);
    g_clear_object (&ldata->src_gsfs);
    nc_galaxy_shape_pop_clear (&ldata->src_pop);
    ldata->loaded = FALSE;
  }

  ldata->valid = FALSE;

  if (st != NULL)
  {
    NcmObjDictStr *source = nc_galaxy_wl_obs_peek_shape_interp_source (obs);
    GObject *src_gsfs     = ncm_obj_dict_str_peek (source, "shape-factor");
    GObject *src_pop      = ncm_obj_dict_str_peek (source, "shape-pop");

    if (!NC_IS_GALAXY_SHAPE_FACTOR_SPECTRAL (src_gsfs) || !NC_IS_GALAXY_SHAPE_POP (src_pop))
      g_error ("nc_galaxy_shape_factor_spectral: the catalog has interpolants, but its record "
               "lacks the NcGalaxyShapeFactorSpectral `shape-factor' or the NcGalaxyShapePop `shape-pop'.");

    ncm_spectral_tensor_clear (&ldata->st);
    ldata->st       = ncm_spectral_tensor_ref (st);
    ldata->src_gsfs = g_object_ref (NC_GALAXY_SHAPE_FACTOR_SPECTRAL (src_gsfs));
    ldata->src_pop  = nc_galaxy_shape_pop_ref (NC_GALAXY_SHAPE_POP (src_pop));
    ldata->loaded   = TRUE;
  }
}

/*
 * Writes the interpolant of the galaxy to row i of obs when it is valid for the observed
 * ellipticity and noise of the galaxy. nc_galaxy_shape_factor_spectral_set_obs_source()
 * makes the written interpolants current.
 */
static void
_nc_galaxy_shape_factor_spectral_ldata_write_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorSpectralData *ldata = data->ldata;

  if (ldata->valid && (ldata->epsilon_obs_1 == data->epsilon_obs_1) && (ldata->epsilon_obs_2 == data->epsilon_obs_2) &&
      (ldata->std_noise == data->std_noise))
    nc_galaxy_wl_obs_set_shape_interp (obs, i, ldata->st);
}

static void
_nc_galaxy_shape_factor_spectral_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void _nc_galaxy_shape_factor_spectral_ldata_free (gpointer p);

/*
 * The state per galaxy is the interpolant of ln P in g, empty until
 * nc_galaxy_shape_factor_spectral_build_interp() builds it or a row of a catalog that carries
 * one is read. Nothing of it is a column of the catalog.
 */
static void
_nc_galaxy_shape_factor_spectral_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorSpectralData *ldata = g_new0 (NcGalaxyShapeFactorSpectralData, 1);

  ldata->st       = NULL;
  ldata->valid    = FALSE;
  ldata->loaded   = FALSE;
  ldata->src_gsfs = NULL;
  ldata->src_pop  = NULL;

  ldata->pop_idx    = g_array_new (FALSE, FALSE, sizeof (guint));
  ldata->pop_bounds = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ldata->pop_params = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ldata->x          = g_array_new (FALSE, TRUE, sizeof (gdouble));

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_spectral_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_spectral_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_spectral_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_spectral_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_spectral_ldata_free (gpointer p)
{
  NcGalaxyShapeFactorSpectralData *ldata = p;

  ncm_spectral_tensor_clear (&ldata->st);
  g_clear_object (&ldata->src_gsfs);
  nc_galaxy_shape_pop_clear (&ldata->src_pop);
  g_array_unref (ldata->pop_idx);
  g_array_unref (ldata->pop_bounds);
  g_array_unref (ldata->pop_params);
  g_array_unref (ldata->x);
  g_free (ldata);
}

/*
 * The |g| range data_prepare() asks of the interpolant of a galaxy with restrict-range:
 * the bound of |g| = |(1 + m) g_t e^{2 i phi} + c| at the corner of the prior box, |1 +
 * m| |g_t| + |c|, with the convergence and g_t at the smallest projected radius the box
 * of the halo center allows, the source redshift z_max and the mass at its upper bound
 * (the bound fields of @data), when it is below interp-g-max. At fixed radius, source
 * redshift and concentration, the convergence and the shear grow with the mass and the
 * source redshift, so with kappa < 1 at the corner |g_t| is largest there. 0 (the whole
 * interp-g-max) when the bound is not ready, when the bound does not hold, when the
 * radius can reach the halo itself, without a finite source redshift behind the lens,
 * when kappa >= 1 at the corner, and when the bound is 0.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_reach (NcGalaxyShapeFactorSpectralPrivate *self, NcGalaxyShapeFactorData *data, const gdouble z_max)
{
  const gdouble R_min = data->bound_R_min;
  const gdouble kappa = data->bound_kappa;
  const gdouble g_t   = data->bound_g_t;
  gdouble g_cal;

  if (!self->restrict_range || !data->bound_ready || !data->bound_holds)
    return 0.0;

  if (!(R_min > 0.0) || !gsl_finite (z_max) || !(z_max > data->bound_z_cl))
    return 0.0;

  if (!gsl_finite (kappa) || !gsl_finite (g_t) || !(kappa < 1.0))
    return 0.0;

  g_cal = fabs (1.0 + data->m) * g_t + hypot (data->c1, data->c2);

  return ((g_cal > 0.0) && (g_cal < self->interp_g_max)) ? g_cal : 0.0;
}

static gdouble _nc_galaxy_shape_factor_spectral_eval_ln (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble std_noise, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, NcmMatrix **panels);
static gboolean _nc_galaxy_shape_factor_spectral_interp_valid (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static void _nc_galaxy_shape_factor_spectral_interp_build (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static void _nc_galaxy_shape_factor_spectral_interp_record (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
static void _nc_galaxy_shape_factor_spectral_check_loaded (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data);
static gdouble _nc_galaxy_shape_factor_spectral_interp_center (NcGalaxyShapeFactorSpectralPrivate *self, const gdouble e);
static void _nc_galaxy_shape_factor_spectral_pop_vars (NcGalaxyShapePop *pop, NcGalaxyShapeFactorSpectralData *ldata);
static void _nc_galaxy_shape_factor_spectral_pop_point (NcGalaxyShapePop *pop, NcGalaxyShapeFactorSpectralData *ldata);

/* The |g| range of the interpolant of a galaxy: the range data_prepare() set, or interp-g-max */
static inline gdouble
_nc_galaxy_shape_factor_spectral_interp_range (NcGalaxyShapeFactorSpectralPrivate *self, const NcGalaxyShapeFactorSpectralData *ldata)
{
  return (ldata->g_target > 0.0) ? ldata->g_target : self->interp_g_max;
}

/* Whether the galaxy of data holds an interpolant: one built for it or one loaded from its row */
static inline gboolean
_nc_galaxy_shape_factor_spectral_has_interp (const NcGalaxyShapeFactorSpectralData *ldata)
{
  return ldata->loaded || ldata->valid;
}

/*
 * The interpolant the galaxy of data holds, checked for an evaluation at eps_obs with the
 * population pop: valid as recorded; or loaded and not yet compared with the fit, which
 * _nc_galaxy_shape_factor_spectral_check_loaded() then does; otherwise the call aborts, since
 * an interpolant is never rebuilt during an evaluation.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_check (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                               const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectralData *ldata = data->ldata;

  if (_nc_galaxy_shape_factor_spectral_interp_valid (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2))
    return;

  if (ldata->loaded && !ldata->valid)
  {
    _nc_galaxy_shape_factor_spectral_check_loaded (gsfs, pop, data);

    if (_nc_galaxy_shape_factor_spectral_interp_valid (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2))
      return;
  }

  g_error ("nc_galaxy_shape_factor_spectral: the interpolant of the galaxy does not hold for this evaluation "
           "(eps_obs = (% 22.15g, % 22.15g), std_noise = % 22.15g): it was built for another observed ellipticity, noise, "
           "setting or population. Interpolants are built by nc_galaxy_shape_factor_spectral_build_interp() and never "
           "during an evaluation.", epsilon_obs_1, epsilon_obs_2, data->std_noise);
}

/*
 * Checks the interpolant of the galaxy of data, when it holds one, for the current
 * population (_nc_galaxy_shape_factor_spectral_interp_check()): the first use of an
 * interpolant loaded from the catalog compares it with the fit here. Nothing is built.
 */
static void
_nc_galaxy_shape_factor_spectral_data_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max)
{
  NcGalaxyShapeFactorSpectral *gsfs = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf);

  if (_nc_galaxy_shape_factor_spectral_has_interp (data->ldata))
  {
    NcGalaxyShapePop *pop = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

    _nc_galaxy_shape_factor_spectral_interp_check (gsfs, pop, data, data->epsilon_obs_1, data->epsilon_obs_2);
  }
}

static gdouble
_nc_galaxy_shape_factor_spectral_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return exp (_nc_galaxy_shape_factor_spectral_eval_ln_marginal (gsf, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

/*
 * ln P (eps_obs | g) for the galaxy of data. When the galaxy holds an interpolant, checked by
 * _nc_galaxy_shape_factor_spectral_interp_check(), and |g| is in its range, it is the value
 * of the interpolant. Otherwise it is the value of
 * _nc_galaxy_shape_factor_spectral_eval_ln(). The interpolant is evaluated at
 * (rho, psi), g = g_c + rho exp (i psi) in the frame where eps_obs is on the positive
 * real axis, with the sign of the imaginary part of g dropped: ln P is even under that
 * reflection.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectral *gsfs      = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf);
  NcGalaxyShapeFactorSpectralData *ldata = data->ldata;

  if (!_nc_galaxy_shape_factor_spectral_has_interp (ldata))
    return _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, data->pop_data, data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, NULL);

  _nc_galaxy_shape_factor_spectral_interp_check (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2);

  if (g_1 * g_1 + g_2 * g_2 > gsl_pow_2 (ldata->g_max))
    return _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, data->pop_data, data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2, NULL);

  {
    {
      const gdouble gf_1 = g_1 * ldata->cos_e + g_2 * ldata->sin_e - ldata->g_c;
      const gdouble gf_2 = g_2 * ldata->cos_e - g_1 * ldata->sin_e;

      gdouble *x = &g_array_index (ldata->x, gdouble, 0);

      x[0] = hypot (gf_1, gf_2);
      x[1] = atan2 (fabs (gf_2), gf_1);

      return ncm_spectral_tensor_eval (ldata->st, ldata->x->len, x);
    }
  }
}

/*
 * eval_nodes of the class: for a galaxy that holds an interpolant, checked once per call, the
 * shears with |g| in its range are evaluated in one call of
 * ncm_spectral_tensor_eval_slice_array(), with the population variables fixed.
 * Their points are those of _nc_galaxy_shape_factor_spectral_eval_ln_marginal(), computed in
 * one loop that vectorizes: rho by sqrt in place of hypot (rho < 2 here) and psi by the
 * vector atan2. The other shears, and every shear of a galaxy without an interpolant, take
 * the integral. The values of the interpolant agree with those of
 * the evaluation one shear at a time to rounding, not bit by bit.
 */
static void
_nc_galaxy_shape_factor_spectral_eval_nodes (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n,
                                             const gdouble *g_1, const gdouble *g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                             const gboolean want_log, gdouble *out)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  NcGalaxyShapeFactorSpectralWork **ptr;
  NcGalaxyShapeFactorSpectralWork *work;
  guint m_in = 0;
  gdouble *g_buf;
  guint *j_buf;
  gdouble g_max2;
  guint j, q;

  if (!_nc_galaxy_shape_factor_spectral_has_interp (ldata))
  {
    for (j = 0; j < n; j++)
    {
      const gdouble lnP = _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, data->pop_data, data->std_noise, g_1[j], g_2[j],
                                                                    epsilon_obs_1, epsilon_obs_2, NULL);

      out[j] = want_log ? lnP : exp (lnP);
    }

    return;
  }

  _nc_galaxy_shape_factor_spectral_interp_check (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2);
  g_max2 = gsl_pow_2 (ldata->g_max);

  ptr  = ncm_memory_pool_get (self->work_pool);
  work = *ptr;
  g_array_set_size (work->g_arr, 2 * n);
  g_array_set_size (work->j_arr, n);
  g_buf = &g_array_index (work->g_arr, gdouble, 0);
  j_buf = &g_array_index (work->j_arr, guint, 0);

  for (j = 0; j < n; j++)
  {
    if (g_1[j] * g_1[j] + g_2[j] * g_2[j] > g_max2)
    {
      out[j] = _nc_galaxy_shape_factor_spectral_eval_ln (gsfs, pop, data->pop_data, data->std_noise, g_1[j], g_2[j],
                                                         epsilon_obs_1, epsilon_obs_2, NULL);
      continue;
    }

    g_buf[2 * m_in]     = g_1[j];
    g_buf[2 * m_in + 1] = g_2[j];
    j_buf[m_in]         = j;
    m_in++;
  }

  if (m_in > 0)
  {
    const guint m       = m_in;
    const gdouble cos_e = ldata->cos_e;
    const gdouble sin_e = ldata->sin_e;
    const gdouble g_c   = ldata->g_c;
    const guint d       = ldata->x->len;
    const gdouble *g    = &g_array_index (work->g_arr, gdouble, 0);
    gdouble *x;

    g_array_set_size (work->x_arr, 2 * m);
    x = &g_array_index (work->x_arr, gdouble, 0);

    #pragma omp simd

    for (q = 0; q < m; q++)
    {
      const gdouble gf_1 = g[2 * q] * cos_e + g[2 * q + 1] * sin_e - g_c;
      const gdouble gf_2 = g[2 * q + 1] * cos_e - g[2 * q] * sin_e;

      x[2 * q]     = sqrt (gf_1 * gf_1 + gf_2 * gf_2);
      x[2 * q + 1] = atan2 (fabs (gf_2), gf_1);
    }

    g_array_set_size (work->x_pop, d - 2);

    if (d > 2)
      memcpy (&g_array_index (work->x_pop, gdouble, 0), &g_array_index (ldata->x, gdouble, 2), sizeof (gdouble) * (d - 2));

    ncm_spectral_tensor_eval_slice_array (ldata->st, work->x_pop, work->x_arr, &work->f_arr);

    for (q = 0; q < m; q++)
      out[j_buf[q]] = g_array_index (work->f_arr, gdouble, q);
  }

  ncm_memory_pool_return (ptr);

  if (!want_log)
    for (j = 0; j < n; j++)
      out[j] = exp (out[j]);
}

/* What the function tabulated along a line needs */
typedef struct _NcGalaxyShapeFactorSpectralLine
{
  NcGalaxyShapeFactorSpectral *gsfs;
  NcGalaxyShapePop *pop;
  NcGalaxyShapeFactorData *data;
  gdouble u_1;
  gdouble u_2;
  gdouble b_1;
  gdouble b_2;
  gdouble epsilon_obs_1;
  gdouble epsilon_obs_2;
} NcGalaxyShapeFactorSpectralLine;

/* ln P at the position t of the line, by the integral */
static gdouble
_nc_galaxy_shape_factor_spectral_line_ln_P (gpointer user_data, gdouble t)
{
  const NcGalaxyShapeFactorSpectralLine *line = user_data;

  return _nc_galaxy_shape_factor_spectral_eval_ln (line->gsfs, line->pop, line->data->pop_data, line->data->std_noise,
                                                   line->b_1 + t * line->u_1, line->b_2 + t * line->u_2,
                                                   line->epsilon_obs_1, line->epsilon_obs_2, NULL);
}

/*
 * eval_line of the class: for a galaxy without an interpolant and with a positive
 * line-abstol, the Chebyshev series of ln P at b + t u in t on [min t, max t], expanded
 * from integrals on nested grids from level NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MIN
 * up to the largest level whose grid has at most n / 2 points (and at most
 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MAX), and accepted at line-abstol. FALSE, with
 * out unchanged, for a galaxy with an interpolant, with line-abstol 0, for fewer than 6
 * shears or a single position, when an end of the line has |g| >
 * NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_G_MAX, or when the series does not reach the
 * tolerance. Nothing is kept between calls.
 */
static gboolean
_nc_galaxy_shape_factor_spectral_eval_line (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n,
                                            const gdouble *t, const gdouble u_1, const gdouble u_2, const gdouble b_1, const gdouble b_2,
                                            const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gboolean want_log, gdouble *out)
{
  NcGalaxyShapeFactorSpectral *gsfs               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (gsf);
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  guint level_cap                                 = 0;
  gdouble t_min, t_max;
  guint j;

  if (_nc_galaxy_shape_factor_spectral_has_interp (data->ldata) || (self->line_abstol == 0.0))
    return FALSE;

  while ((level_cap < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MAX) && ((1u << (level_cap + 1)) + 1 <= n / 2))
    level_cap++;

  if (level_cap <= NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MIN)
    return FALSE;

  t_min = t[0];
  t_max = t[0];

  for (j = 1; j < n; j++)
  {
    t_min = GSL_MIN (t_min, t[j]);
    t_max = GSL_MAX (t_max, t[j]);
  }

  /* |b + t u| is largest at an end of the segment */
  if (!(t_max > t_min) ||
      (hypot (b_1 + t_min * u_1, b_2 + t_min * u_2) > NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_G_MAX) ||
      (hypot (b_1 + t_max * u_1, b_2 + t_max * u_2) > NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_G_MAX))
    return FALSE;

  {
    NcGalaxyShapeFactorSpectralLine line      = {gsfs, pop, data, u_1, u_2, b_1, b_2, epsilon_obs_1, epsilon_obs_2};
    NcGalaxyShapeFactorSpectralLineWork **ptr = ncm_memory_pool_get (self->line_pool);
    NcGalaxyShapeFactorSpectralLineWork *lw   = *ptr;
    gboolean converged                        = FALSE;

    ncm_spectral_compute_chebyshev_coeffs_adaptive_try (lw->spectral, &_nc_galaxy_shape_factor_spectral_line_ln_P, t_min, t_max,
                                                        NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MIN, level_cap,
                                                        0.0, self->line_abstol, &lw->coeffs, &line, &converged);

    if (converged)
    {
      const gdouble mid  = 0.5 * (t_max + t_min);
      const gdouble ihlf = 2.0 / (t_max - t_min);

      for (j = 0; j < n; j++)
      {
        const gdouble lnP = ncm_spectral_chebyshev_eval (lw->coeffs, (t[j] - mid) * ihlf);

        out[j] = want_log ? lnP : exp (lnP);
      }
    }

    ncm_memory_pool_return (ptr);

    return converged;
  }
}

static gpointer
_nc_galaxy_shape_factor_spectral_line_work_new (gpointer userdata)
{
  NcGalaxyShapeFactorSpectralLineWork *lw = g_new0 (NcGalaxyShapeFactorSpectralLineWork, 1);

  lw->spectral = ncm_spectral_new_with_max_level (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_LINE_LEVEL_MAX);
  lw->coeffs   = g_array_new (FALSE, FALSE, sizeof (gdouble));

  return lw;
}

static void
_nc_galaxy_shape_factor_spectral_line_work_free (gpointer p)
{
  NcGalaxyShapeFactorSpectralLineWork *lw = p;

  ncm_spectral_free (lw->spectral);
  g_array_unref (lw->coeffs);
  g_free (lw);
}

static gboolean _nc_galaxy_shape_factor_spectral_pop_inside (NcGalaxyShapePop *pop, NcGalaxyShapeFactorSpectralData *ldata);

/*
 * Whether the galaxy of data has an interpolant built for these eps_obs, its std_noise,
 * the current state of pop and the current settings of the class. Without population
 * variables the state of pop is the one of the build. With them, a new state of pop is
 * accepted when _nc_galaxy_shape_factor_spectral_pop_inside() finds it inside the table.
 */
static gboolean
_nc_galaxy_shape_factor_spectral_interp_valid (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                               const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  const guint n_pop                               = ldata->pop_idx->len;

  if (!(ldata->valid && (ldata->epsilon_obs_1 == epsilon_obs_1) && (ldata->epsilon_obs_2 == epsilon_obs_2) &&
        (ldata->std_noise == data->std_noise) && (ldata->abstol == self->interp_abstol) &&
        (ldata->reltol == self->reltol) &&
        (ldata->min_radial_panels == self->interp_min_radial_panels) && (ldata->trim_frac == self->interp_trim_frac) &&
        ((n_pop == 0) || ((ldata->pop_level_min == self->interp_pop_level_min) && (ldata->pop_level_cap == self->interp_pop_level_cap) &&
                          (ldata->pop_one_integral == self->interp_pop_one_integral) &&
                          (ldata->compress_frac == self->interp_compress_frac)))))
    return FALSE;

  if (ldata->pop_pkey == ncm_model_state_get_pkey (NCM_MODEL (pop)))
    return TRUE;

  return (n_pop > 0) && _nc_galaxy_shape_factor_spectral_pop_inside (pop, ldata);
}

/*
 * Whether the table of a galaxy with population variables holds for the current state of
 * pop: every free parameter is a variable of the table with the bounds recorded, and every
 * parameter that is not a variable is fixed at the value recorded. A variable may also be
 * fixed. Then the values of the variables become the population part of the point of the
 * table (ldata->x) and the state is recorded. A value outside the range of its variable
 * aborts.
 */
static gboolean
_nc_galaxy_shape_factor_spectral_pop_inside (NcGalaxyShapePop *pop, NcGalaxyShapeFactorSpectralData *ldata)
{
  NcmModel *model   = NCM_MODEL (pop);
  const guint len   = ncm_model_len (model);
  const guint n_pop = ldata->pop_idx->len;
  guint k, v = 0;

  if (ldata->pop_params->len != len)
    return FALSE;

  for (k = 0; k < len; k++)
  {
    const gboolean free = (ncm_model_param_get_ftype (model, k) == NCM_PARAM_TYPE_FREE);

    if ((v < n_pop) && (g_array_index (ldata->pop_idx, guint, v) == k))
    {
      if (free &&
          ((g_array_index (ldata->pop_bounds, gdouble, 2 * v) != ncm_model_param_get_lower_bound (model, k)) ||
           (g_array_index (ldata->pop_bounds, gdouble, 2 * v + 1) != ncm_model_param_get_upper_bound (model, k))))
        return FALSE;

      v++;
    }
    else if (free || (ncm_model_param_get (model, k) != g_array_index (ldata->pop_params, gdouble, k)))
    {
      return FALSE;
    }
  }

  if (v != n_pop)
    return FALSE;

  _nc_galaxy_shape_factor_spectral_pop_point (pop, ldata);

  return TRUE;
}

/*
 * The population variables of the table for the current state of pop: the indices of its
 * free parameters (ldata->pop_idx) and their bounds (ldata->pop_bounds), and the length of
 * the point of the table (ldata->x).
 */
static void
_nc_galaxy_shape_factor_spectral_pop_vars (NcGalaxyShapePop *pop, NcGalaxyShapeFactorSpectralData *ldata)
{
  NcmModel *model = NCM_MODEL (pop);
  const guint len = ncm_model_len (model);
  guint k;

  g_array_set_size (ldata->pop_idx, 0);
  g_array_set_size (ldata->pop_bounds, 0);

  for (k = 0; k < len; k++)
  {
    if (ncm_model_param_get_ftype (model, k) == NCM_PARAM_TYPE_FREE)
    {
      const gdouble lb = ncm_model_param_get_lower_bound (model, k);
      const gdouble ub = ncm_model_param_get_upper_bound (model, k);

      g_array_append_val (ldata->pop_idx, k);
      g_array_append_val (ldata->pop_bounds, lb);
      g_array_append_val (ldata->pop_bounds, ub);
    }
  }

  g_array_set_size (ldata->x, 2 + ldata->pop_idx->len);
}

/*
 * Records the state of pop in ldata: its pkey, every parameter value, and the values of the
 * population variables in the point of the table. A free value outside the range of its
 * variable aborts: the table holds only inside.
 */
static void
_nc_galaxy_shape_factor_spectral_pop_point (NcGalaxyShapePop *pop, NcGalaxyShapeFactorSpectralData *ldata)
{
  NcmModel *model = NCM_MODEL (pop);
  const guint len = ncm_model_len (model);
  guint k, v;

  g_array_set_size (ldata->pop_params, len);

  for (k = 0; k < len; k++)
    g_array_index (ldata->pop_params, gdouble, k) = ncm_model_param_get (model, k);

  for (v = 0; v < ldata->pop_idx->len; v++)
  {
    const guint idx   = g_array_index (ldata->pop_idx, guint, v);
    const gdouble lb  = g_array_index (ldata->pop_bounds, gdouble, 2 * v);
    const gdouble ub  = g_array_index (ldata->pop_bounds, gdouble, 2 * v + 1);
    const gdouble val = g_array_index (ldata->pop_params, gdouble, idx);

    if (!((val >= lb) && (val <= ub)))
      g_error ("nc_galaxy_shape_factor_spectral: the population parameter %s = % 22.15g is outside the range "
               "[% 22.15g, % 22.15g] of the interpolant.", ncm_model_param_name (model, idx), val, lb, ub);

    g_array_index (ldata->x, gdouble, 2 + v) = val;
  }

  ldata->pop_pkey = ncm_model_state_get_pkey (model);
}

static void _nc_galaxy_shape_factor_spectral_interp_kernel (guint n_vars, const guint *n, const gdouble * const *x, gdouble *fval, gpointer user_data);
static void _nc_galaxy_shape_factor_spectral_interp_r_edges (const gdouble R, const gdouble L, const guint n, GArray *r_edges);
static gint _nc_galaxy_shape_factor_spectral_cmp_double (gconstpointer a, gconstpointer b);

static void _nc_galaxy_shape_factor_spectral_interp_pop_kernel (guint n_vars, const guint *n, const gdouble * const *x, gdouble *fval, gpointer user_data);

/*
 * What the function approximated in g needs. With population variables, pop is a copy of the
 * population of the fit and pop_data its own data, whose free parameters the kernel sets.
 */
typedef struct _NcGalaxyShapeFactorSpectralInterp
{
  NcGalaxyShapeFactorSpectral *gsfs;
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  gdouble std_noise;
  gdouble e;             /* |eps_obs| */
  gdouble g_c;           /* center of the coordinates */
  const GArray *pop_idx; /* guint, the parameters of pop that are variables 2, 3, ... */

  /*
   * interp-pop-one-integral, for the populations of one call of the kernel: the Ka values of
   * the parameter of the factor at the origin times the Kb values of the one at 1. side_pops
   * holds one copy of pop per value, the Ka of the origin first, with side_datas their data.
   */
  GPtrArray *side_pops;
  GPtrArray *side_datas;
  GArray *grid_a;    /* gdouble, the exponent at the origin of each of the Ka */
  GArray *grid_bpop; /* gdouble, the exponent at 1 of the population of each of the Kb */
  GArray *grid_b;    /* gdouble, the exponent at 1 of the variable of integration of each of the Kb */
  GArray *grid_ln_N; /* gdouble, ln of the constant of each pair, Ka x Kb */
  GArray *grid_res;  /* gdouble, the integrals of a shear, Ka x Kb */
  NcmSpectralPolarWeightGrid *wg;
  guint Ka;
  guint Kb;
  gint var_o; /* the population variable that is the parameter of the origin, or -1 */
  gint var_1; /* the one of the parameter at 1, or -1 */
} NcGalaxyShapeFactorSpectralInterp;

/*
 * Builds the interpolant of ln P in g of the galaxy of data, for the observed ellipticity
 * (epsilon_obs_1, epsilon_obs_2), on a region that contains the disk
 * |g| <= NcGalaxyShapeFactorSpectral:interp-g-max, and records in the data of the galaxy
 * what it was built for. In the frame where eps_obs is on the positive real axis,
 * g = g_c + rho exp (i psi), with g_c the shear that maps chi_I = 0 to eps_obs:
 * e / (1 + sqrt (1 - e^2)) for TRACE (1 for e >= 1) and e for TRACE_DET, e = |eps_obs|
 * (_nc_galaxy_shape_factor_spectral_interp_center()). ln P has its narrowest structure at
 * g_c, of width std_noise / 2 (TRACE) or std_noise (TRACE_DET). The region is
 * [0, g_max + g_c] x [0, pi] for g_c < g_max and otherwise
 * [g_c - g_max, g_c + g_max] x [pi - asin (g_max / g_c), pi]; psi < 0 is covered by the
 * reflection. The radial edges are those of _nc_galaxy_shape_factor_spectral_interp_r_edges()
 * for that width that are inside the region, and the angular panels span pi / 2 at most. Each
 * free parameter of pop is one more variable, over its bounds, with the levels
 * interp-pop-level-min and interp-pop-level-cap; its values are computed on copies of pop, so
 * that the builds of several galaxies may run in parallel. ncm_spectral_tensor_builder_build()
 * makes a node a panel node when the change of ln P between its two grids is at most
 * NcGalaxyShapeFactorSpectral:interp-abstol. A build that reaches the maximum depth aborts
 * with the galaxy in the message.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_build (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                               const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  NcGalaxyShapeFactorSpectralWork **ptr           = ncm_memory_pool_get (self->work_pool);
  NcGalaxyShapeFactorSpectralWork *work           = *ptr;
  const gdouble g_max                             = _nc_galaxy_shape_factor_spectral_interp_range (self, ldata);
  const gdouble e                                 = hypot (epsilon_obs_1, epsilon_obs_2);
  const gdouble g_c                               = _nc_galaxy_shape_factor_spectral_interp_center (self, e);
  const gdouble L                                 = (self->axis_ratio_var ? 0.5 : 1.0) * data->std_noise;
  const gdouble rho_min                           = (g_c < g_max) ? 0.0 : g_c - g_max;
  const gdouble rho_max                           = g_c + g_max;
  const gdouble psi_min                           = (g_c < g_max) ? 0.0 : M_PI - asin (g_max / g_c);
  const gdouble psi_max                           = M_PI;
  NcGalaxyShapeFactorSpectralInterp interp        = {
    gsfs, pop, data->pop_data, data->std_noise, e, g_c, ldata->pop_idx,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, -1, -1
  };
  GError *error = NULL;
  guint i, n_pop, n_vars;

  _nc_galaxy_shape_factor_spectral_interp_r_edges (rho_max, L, self->interp_min_radial_panels, work->r_edges_all);

  g_array_set_size (work->r_edges, 0);
  g_array_append_val (work->r_edges, rho_min);

  for (i = 0; i < work->r_edges_all->len; i++)
  {
    const gdouble rho = g_array_index (work->r_edges_all, gdouble, i);

    if (rho > rho_min + 0.125 * data->std_noise)
      g_array_append_val (work->r_edges, rho);
  }

  /* A range narrower than std_noise / 8 (restrict-range) keeps its upper end */
  if (work->r_edges->len == 1)
    g_array_append_val (work->r_edges, rho_max);

  g_array_set_size (work->theta_edges, 0);
  g_array_append_val (work->theta_edges, psi_min);

  if (psi_max - psi_min > 0.5 * M_PI)
  {
    const gdouble psi_mid = 0.5 * (psi_min + psi_max);

    g_array_append_val (work->theta_edges, psi_mid);
  }

  g_array_append_val (work->theta_edges, psi_max);

  ldata->valid = FALSE;
  _nc_galaxy_shape_factor_spectral_pop_vars (pop, ldata);
  _nc_galaxy_shape_factor_spectral_pop_point (pop, ldata);
  n_pop  = ldata->pop_idx->len;
  n_vars = 2 + n_pop;

  if ((n_pop > 0) && (self->interp_pop_level_min > self->interp_pop_level_cap))
    g_error ("nc_galaxy_shape_factor_spectral: interp-pop-level-min (%u) is above interp-pop-level-cap (%u).",
             self->interp_pop_level_min, self->interp_pop_level_cap);

  if ((ldata->st != NULL) && (ncm_spectral_tensor_get_n_vars (ldata->st) != n_vars))
    ncm_spectral_tensor_clear (&ldata->st);

  if (ldata->st == NULL)
    ldata->st = ncm_spectral_tensor_new (n_vars);

  g_array_set_size (work->lower, n_vars);
  g_array_set_size (work->upper, n_vars);
  g_array_set_size (work->level_min, n_vars);
  g_array_set_size (work->level_cap, n_vars);

  g_array_index (work->lower, gdouble, 0)   = rho_min;
  g_array_index (work->lower, gdouble, 1)   = psi_min;
  g_array_index (work->upper, gdouble, 0)   = g_array_index (work->r_edges, gdouble, work->r_edges->len - 1);
  g_array_index (work->upper, gdouble, 1)   = psi_max;
  g_array_index (work->level_min, guint, 0) = ncm_spectral_tensor_builder_get_level_min (self->builder_interp);
  g_array_index (work->level_min, guint, 1) = ncm_spectral_tensor_builder_get_level_min (self->builder_interp);
  g_array_index (work->level_cap, guint, 0) = ncm_spectral_tensor_builder_get_level_cap (self->builder_interp);
  g_array_index (work->level_cap, guint, 1) = ncm_spectral_tensor_builder_get_level_cap (self->builder_interp);

  for (i = 0; i < n_pop; i++)
  {
    g_array_index (work->lower, gdouble, 2 + i)   = g_array_index (ldata->pop_bounds, gdouble, 2 * i);
    g_array_index (work->upper, gdouble, 2 + i)   = g_array_index (ldata->pop_bounds, gdouble, 2 * i + 1);
    g_array_index (work->level_min, guint, 2 + i) = self->interp_pop_level_min;
    g_array_index (work->level_cap, guint, 2 + i) = self->interp_pop_level_cap;
  }

  ncm_spectral_tensor_set_domain (ldata->st, n_vars, &g_array_index (work->lower, gdouble, 0), n_vars, &g_array_index (work->upper, gdouble, 0));
  ncm_spectral_tensor_split_empty (ldata->st, 0, work->r_edges->len, &g_array_index (work->r_edges, gdouble, 0));
  ncm_spectral_tensor_split_empty (ldata->st, 1, work->theta_edges->len, &g_array_index (work->theta_edges, gdouble, 0));

  if (n_pop == 0)
  {
    ncm_spectral_tensor_builder_build (self->builder_interp, &_nc_galaxy_shape_factor_spectral_interp_kernel, ldata->st,
                                       &g_array_index (work->level_min, guint, 0), &g_array_index (work->level_cap, guint, 0),
                                       0.0, self->interp_abstol, &interp, NULL, &error);
  }
  else
  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

    interp.pop        = NC_GALAXY_SHAPE_POP (ncm_model_dup (NCM_MODEL (pop), ser));
    interp.pop_data   = nc_galaxy_shape_pop_data_new (interp.pop);
    interp.side_pops  = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_shape_pop_free);
    interp.side_datas = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_shape_pop_data_unref);
    interp.grid_a     = g_array_new (FALSE, FALSE, sizeof (gdouble));
    interp.grid_bpop  = g_array_new (FALSE, FALSE, sizeof (gdouble));
    interp.grid_b     = g_array_new (FALSE, FALSE, sizeof (gdouble));
    interp.grid_ln_N  = g_array_new (FALSE, FALSE, sizeof (gdouble));
    interp.grid_res   = g_array_new (FALSE, FALSE, sizeof (gdouble));
    ncm_serialize_free (ser);

    ncm_spectral_tensor_builder_build (self->builder_interp, &_nc_galaxy_shape_factor_spectral_interp_pop_kernel, ldata->st,
                                       &g_array_index (work->level_min, guint, 0), &g_array_index (work->level_cap, guint, 0),
                                       0.0, self->interp_abstol, &interp, NULL, &error);

    if ((error == NULL) && (self->interp_compress_frac > 0.0))
      ncm_spectral_tensor_compress (ldata->st, self->interp_compress_frac * self->interp_abstol);

    g_clear_pointer (&interp.wg, _ncm_spectral_polar_weight_grid_free);
    g_ptr_array_unref (interp.side_pops);
    g_ptr_array_unref (interp.side_datas);
    g_array_unref (interp.grid_a);
    g_array_unref (interp.grid_bpop);
    g_array_unref (interp.grid_b);
    g_array_unref (interp.grid_ln_N);
    g_array_unref (interp.grid_res);
    nc_galaxy_shape_pop_data_unref (interp.pop_data);
    nc_galaxy_shape_pop_free (interp.pop);
  }

  ncm_memory_pool_return (ptr);

  if (error != NULL)
    g_error ("nc_galaxy_shape_factor_spectral: the interpolant in g: %s for the population %s, std_noise = % 22.15g, "
             "eps_obs = (% 22.15g, % 22.15g), interp-g-max = % 22.15g.",
             error->message, G_OBJECT_TYPE_NAME (pop), data->std_noise, epsilon_obs_1, epsilon_obs_2, g_max);

  _nc_galaxy_shape_factor_spectral_interp_record (gsfs, pop, data, epsilon_obs_1, epsilon_obs_2);
}

/*
 * The first radial grid of the interpolant in g on [0, R], in ascending order, for the
 * narrowest structure of ln P, of width L, at rho = 0; ncm_spectral_tensor_builder_build()
 * refines it. The edges are 0 and R, the inner edges k R / n of n equal intervals, and the
 * distances L w^k below R, w = NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INTERP_EDGE_RATIO. An edge
 * at most L / 4 above the edge kept before it is left out, except R.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_r_edges (const gdouble R, const gdouble L, const guint n, GArray *r_edges)
{
  const gdouble org = 0.0;
  gdouble w;
  guint i, m;

  g_array_set_size (r_edges, 0);
  g_array_append_val (r_edges, org);
  g_array_append_val (r_edges, R);

  for (i = 1; i < n; i++)
  {
    const gdouble x = R * i / n;

    g_array_append_val (r_edges, x);
  }

  for (w = L; w < 2.0 * R; w *= NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INTERP_EDGE_RATIO)
  {
    if (w < R)
      g_array_append_val (r_edges, w);
  }

  g_array_sort (r_edges, &_nc_galaxy_shape_factor_spectral_cmp_double);

  for (i = 1, m = 1; i < r_edges->len; i++)
  {
    const gdouble x = g_array_index (r_edges, gdouble, i);

    if ((x - g_array_index (r_edges, gdouble, m - 1) > 0.25 * L) || (i + 1 == r_edges->len))
      g_array_index (r_edges, gdouble, m++) = x;
  }

  g_array_set_size (r_edges, m);
}

/* Orders two doubles for g_array_sort() */
static gint
_nc_galaxy_shape_factor_spectral_cmp_double (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x > y) - (x < y);
}

/*
 * Marks the interpolant of the galaxy of data as valid and records what it is valid for:
 * these eps_obs, the noise of the galaxy, the current state of pop and the current settings
 * of the class, the frame of the interpolant (the direction of eps_obs and the center), and
 * its range in |g|, the upper end of its radial variable minus the center, which is the
 * range it was built for, also when it was loaded from a catalog. The population variables
 * (ldata->pop_idx and ldata->pop_bounds) are set before.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_record (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  const gdouble e                                 = hypot (epsilon_obs_1, epsilon_obs_2);

  ldata->valid             = TRUE;
  ldata->epsilon_obs_1     = epsilon_obs_1;
  ldata->epsilon_obs_2     = epsilon_obs_2;
  ldata->std_noise         = data->std_noise;
  ldata->pop_level_min     = self->interp_pop_level_min;
  ldata->pop_level_cap     = self->interp_pop_level_cap;
  ldata->pop_one_integral  = self->interp_pop_one_integral;
  ldata->compress_frac     = self->interp_compress_frac;
  ldata->abstol            = self->interp_abstol;
  ldata->reltol            = self->reltol;
  ldata->min_radial_panels = self->interp_min_radial_panels;
  ldata->trim_frac         = self->interp_trim_frac;
  ldata->cos_e             = (e > 0.0) ? epsilon_obs_1 / e : 1.0;
  ldata->sin_e             = (e > 0.0) ? epsilon_obs_2 / e : 0.0;
  ldata->g_c               = _nc_galaxy_shape_factor_spectral_interp_center (self, e);

  {
    gdouble rho_lo, rho_hi;

    ncm_spectral_tensor_get_range (ldata->st, 0, &rho_lo, &rho_hi);
    ldata->g_max = rho_hi - ldata->g_c;
  }

  _nc_galaxy_shape_factor_spectral_pop_point (pop, ldata);
}

/*
 * Compares the record of the loaded interpolant of the galaxy of data with the fit: the
 * convention, reltol and the interp-* settings of the class, the type of the population, and
 * its parameters against the population variables of the interpolant, the parameters free
 * in the record. A free parameter of the fit must be a variable whose range is its bounds; a
 * variable fixed in the fit holds any value in its range; a parameter that is not a variable
 * must be fixed at the value of the record. When all agree the interpolant is recorded as
 * valid for the current state of pop, with the bounds of the fit as the ranges of its free
 * variables and the ranges of the table for the fixed ones; otherwise the call aborts with
 * every difference in the message.
 */
static void
_nc_galaxy_shape_factor_spectral_check_loaded (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  NcGalaxyShapeFactorSpectralPrivate * const src  = nc_galaxy_shape_factor_spectral_get_instance_private (ldata->src_gsfs);
  const gdouble tol                               = NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SOURCE_RELTOL;
  GString *diff                                   = g_string_new (NULL);
  guint k;

  if (nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfs)) != nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (ldata->src_gsfs)))
    g_string_append (diff, " the ellipticity convention;");

  if (ncm_cmp (self->reltol, src->reltol, tol, 0.0) != 0)
    g_string_append_printf (diff, " reltol (% 22.15g, built with % 22.15g);", self->reltol, src->reltol);

  if (ncm_cmp (self->interp_g_max, src->interp_g_max, tol, 0.0) != 0)
    g_string_append_printf (diff, " interp-g-max (% 22.15g, built with % 22.15g);", self->interp_g_max, src->interp_g_max);

  if (ncm_cmp (self->interp_abstol, src->interp_abstol, tol, 0.0) != 0)
    g_string_append_printf (diff, " interp-abstol (% 22.15g, built with % 22.15g);", self->interp_abstol, src->interp_abstol);

  if (self->interp_min_radial_panels != src->interp_min_radial_panels)
    g_string_append_printf (diff, " interp-min-radial-panels (%u, built with %u);", self->interp_min_radial_panels, src->interp_min_radial_panels);

  if (ncm_cmp (self->interp_trim_frac, src->interp_trim_frac, tol, 0.0) != 0)
    g_string_append_printf (diff, " interp-trim-frac (% 22.15g, built with % 22.15g);", self->interp_trim_frac, src->interp_trim_frac);

  if ((G_OBJECT_TYPE (pop) != G_OBJECT_TYPE (ldata->src_pop)) || (ncm_model_len (NCM_MODEL (pop)) != ncm_model_len (NCM_MODEL (ldata->src_pop))))
  {
    g_string_append_printf (diff, " the population (%s, built with %s);", G_OBJECT_TYPE_NAME (pop), G_OBJECT_TYPE_NAME (ldata->src_pop));
  }
  else
  {
    NcmModel *model     = NCM_MODEL (pop);
    NcmModel *src_model = NCM_MODEL (ldata->src_pop);
    guint n_src_free    = 0;
    guint v             = 0;

    for (k = 0; k < ncm_model_len (src_model); k++)
      n_src_free += (ncm_model_param_get_ftype (src_model, k) == NCM_PARAM_TYPE_FREE) ? 1 : 0;

    g_array_set_size (ldata->pop_idx, 0);
    g_array_set_size (ldata->pop_bounds, 0);
    g_array_set_size (ldata->x, 2 + n_src_free);

    if (ncm_spectral_tensor_get_n_vars (ldata->st) != 2 + n_src_free)
      g_string_append_printf (diff, " the record (%u free population parameters, the interpolants have %u population variables);",
                              n_src_free, ncm_spectral_tensor_get_n_vars (ldata->st) - 2);

    for (k = 0; (k < ncm_model_len (model)) && (ncm_spectral_tensor_get_n_vars (ldata->st) == 2 + n_src_free); k++)
    {
      const gboolean free = (ncm_model_param_get_ftype (model, k) == NCM_PARAM_TYPE_FREE);

      if (ncm_model_param_get_ftype (src_model, k) == NCM_PARAM_TYPE_FREE)
      {
        gdouble a, b;

        ncm_spectral_tensor_get_range (ldata->st, 2 + v, &a, &b);

        if (free)
        {
          if ((ncm_cmp (a, ncm_model_param_get_lower_bound (model, k), tol, 0.0) != 0) ||
              (ncm_cmp (b, ncm_model_param_get_upper_bound (model, k), tol, 0.0) != 0))
            g_string_append_printf (diff, " the bounds of %s ([% 22.15g, % 22.15g], built for [% 22.15g, % 22.15g]);",
                                    ncm_model_param_name (model, k), ncm_model_param_get_lower_bound (model, k),
                                    ncm_model_param_get_upper_bound (model, k), a, b);

          a = ncm_model_param_get_lower_bound (model, k);
          b = ncm_model_param_get_upper_bound (model, k);
        }

        g_array_append_val (ldata->pop_idx, k);
        g_array_append_val (ldata->pop_bounds, a);
        g_array_append_val (ldata->pop_bounds, b);
        v++;
      }
      else if (free)
      {
        g_string_append_printf (diff, " the population parameter %s is free, not a variable of the interpolants;", ncm_model_param_name (model, k));
      }
      else if (ncm_cmp (ncm_model_param_get (model, k), ncm_model_param_get (src_model, k), tol, 0.0) != 0)
      {
        g_string_append_printf (diff, " the population parameter %s (% 22.15g, built with % 22.15g);", ncm_model_param_name (model, k),
                                ncm_model_param_get (model, k), ncm_model_param_get (src_model, k));
      }
    }

    if (ldata->pop_idx->len > 0)
    {
      if (self->interp_pop_level_min != src->interp_pop_level_min)
        g_string_append_printf (diff, " interp-pop-level-min (%u, built with %u);", self->interp_pop_level_min, src->interp_pop_level_min);

      if (self->interp_pop_level_cap != src->interp_pop_level_cap)
        g_string_append_printf (diff, " interp-pop-level-cap (%u, built with %u);", self->interp_pop_level_cap, src->interp_pop_level_cap);

      if (self->interp_pop_one_integral != src->interp_pop_one_integral)
        g_string_append_printf (diff, " interp-pop-one-integral (%d, built with %d);", self->interp_pop_one_integral, src->interp_pop_one_integral);

      if (ncm_cmp (self->interp_compress_frac, src->interp_compress_frac, tol, 0.0) != 0)
        g_string_append_printf (diff, " interp-compress-frac (% 22.15g, built with % 22.15g);", self->interp_compress_frac, src->interp_compress_frac);
    }
  }

  if (diff->len > 0)
    g_error ("nc_galaxy_shape_factor_spectral: the interpolants of the catalog were built for other settings:%s "
             "prepare the catalog again for this fit.", diff->str);

  g_string_free (diff, TRUE);

  _nc_galaxy_shape_factor_spectral_interp_record (gsfs, pop, data, data->epsilon_obs_1, data->epsilon_obs_2);
}

/*
 * The center g_c of the coordinates of the interpolant for |eps_obs| = e: the shear that maps
 * chi_I = 0 to eps_obs, e / (1 + sqrt (1 - e^2)) for TRACE (1 for e >= 1) and e for
 * TRACE_DET.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_interp_center (NcGalaxyShapeFactorSpectralPrivate *self, const gdouble e)
{
  if (!self->axis_ratio_var)
    return e;

  return (e < 1.0) ? e / (1.0 + sqrt (1.0 - e * e)) : 1.0;
}

/*
 * The function approximated by ncm_spectral_tensor_builder_build(): ln P at the shear
 * g_c + rho exp (i psi) for the observed ellipticity (e, 0), by the integral, on the grid
 * of the radii x[0] and the angles x[1], the angle running fastest.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_kernel (guint n_vars, const guint *n, const gdouble * const *x, gdouble *fval, gpointer user_data)
{
  const NcGalaxyShapeFactorSpectralInterp *interp = user_data;
  guint i, j;

  for (j = 0; j < n[1]; j++)
  {
    const gdouble cos_psi = cos (x[1][j]);
    const gdouble sin_psi = sin (x[1][j]);

    for (i = 0; i < n[0]; i++)
      fval[i * n[1] + j] = _nc_galaxy_shape_factor_spectral_eval_ln (interp->gsfs, interp->pop, interp->pop_data, interp->std_noise,
                                                                     interp->g_c + x[0][i] * cos_psi, x[0][i] * sin_psi,
                                                                     interp->e, 0.0, NULL);
  }
}

static gboolean _nc_galaxy_shape_factor_spectral_interp_pop_grid_setup (NcGalaxyShapeFactorSpectralInterp *interp, guint n_vars, const guint *n, const gdouble * const *x);
static void _nc_galaxy_shape_factor_spectral_interp_pop_grid_shear (NcGalaxyShapeFactorSpectralInterp *interp, const gdouble g_1, const gdouble g_2);

/*
 * The function approximated with population variables: ln P at the shear
 * g_c + rho exp (i psi) for the observed ellipticity (e, 0) and the population whose free
 * parameters are the values of the variables 2, 3, ..., on the tensor product of the points
 * of every variable, the last running fastest. With interp-pop-one-integral, for a
 * population with the separated form whose free parameters are those of its two factors
 * (_nc_galaxy_shape_factor_spectral_interp_pop_grid_setup()), each shear takes one integral
 * for all the populations of the call (_nc_galaxy_shape_factor_spectral_interp_pop_grid_shear()).
 * Otherwise each population takes its own integral, its parameters set on the copy
 * interp->pop once per population point.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_pop_kernel (guint n_vars, const guint *n, const gdouble * const *x, gdouble *fval, gpointer user_data)
{
  NcGalaxyShapeFactorSpectralInterp *interp = user_data;
  const guint n_pop                         = n_vars - 2;
  guint n_p                                 = 1;
  guint i, j, p, v;

  for (v = 0; v < n_pop; v++)
    n_p *= n[2 + v];

  if (_nc_galaxy_shape_factor_spectral_interp_pop_grid_setup (interp, n_vars, n, x))
  {
    const gdouble *res   = &g_array_index (interp->grid_res, gdouble, 0);
    const guint Ka       = interp->Ka;
    const guint Kb       = interp->Kb;
    const gboolean o_1st = (interp->var_o < interp->var_1) || (interp->var_1 < 0);

    for (j = 0; j < n[1]; j++)
    {
      const gdouble cos_psi = cos (x[1][j]);
      const gdouble sin_psi = sin (x[1][j]);

      for (i = 0; i < n[0]; i++)
      {
        gdouble *out = &fval[(gsize) (i * n[1] + j) * n_p];
        guint ja, kb;

        _nc_galaxy_shape_factor_spectral_interp_pop_grid_shear (interp, interp->g_c + x[0][i] * cos_psi, x[0][i] * sin_psi);

        /* The variable of lower index runs slowest in fval */
        for (ja = 0; ja < Ka; ja++)
          for (kb = 0; kb < Kb; kb++)
            out[o_1st ? ja * Kb + kb : kb * Ka + ja] = res[ja * Kb + kb];
      }
    }

    return;
  }

  for (p = 0; p < n_p; p++)
  {
    guint r = p;

    for (v = n_pop; v-- > 0;)
    {
      ncm_model_param_set (NCM_MODEL (interp->pop), g_array_index (interp->pop_idx, guint, v), x[2 + v][r % n[2 + v]]);
      r /= n[2 + v];
    }

    nc_galaxy_shape_pop_prepare (interp->pop, interp->pop_data);

    for (j = 0; j < n[1]; j++)
    {
      const gdouble cos_psi = cos (x[1][j]);
      const gdouble sin_psi = sin (x[1][j]);

      for (i = 0; i < n[0]; i++)
        fval[(i * n[1] + j) * n_p + p] = _nc_galaxy_shape_factor_spectral_eval_ln (interp->gsfs, interp->pop, interp->pop_data, interp->std_noise,
                                                                                   interp->g_c + x[0][i] * cos_psi, x[0][i] * sin_psi,
                                                                                   interp->e, 0.0, NULL);
    }
  }
}

/*
 * Prepares interp for the populations of one call of the kernel, for
 * interp-pop-one-integral. FALSE, with nothing prepared, when the property is off, when the
 * population has no separated form, or when a population variable is not the parameter of
 * one of its two factors. Otherwise: one copy of the population per value of each of the
 * two parameters, with its exponent; the constant of each pair; the reference population,
 * at the center of the values of each parameter, in interp->pop; and the weight grid of the
 * exponents. A parameter that is not a variable keeps the value of the fit.
 */
static gboolean
_nc_galaxy_shape_factor_spectral_interp_pop_grid_setup (NcGalaxyShapeFactorSpectralInterp *interp, guint n_vars, const guint *n, const gdouble * const *x)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (interp->gsfs);
  const guint n_pop                               = n_vars - 2;
  gint param_o                                    = -1;
  gint param_1                                    = -1;
  gint var_o                                      = -1;
  gint var_1                                      = -1;
  guint Ka, Kb, k, l, v;

  if (!self->interp_pop_one_integral || !nc_galaxy_shape_pop_reduced_factor_params (interp->pop, &param_o, &param_1))
    return FALSE;

  for (v = 0; v < n_pop; v++)
  {
    const gint idx = g_array_index (interp->pop_idx, guint, v);

    if (idx == param_o)
      var_o = v;
    else if (idx == param_1)
      var_1 = v;
    else
      return FALSE;
  }

  Ka = (var_o >= 0) ? n[2 + var_o] : 1;
  Kb = (var_1 >= 0) ? n[2 + var_1] : 1;

  interp->Ka    = Ka;
  interp->Kb    = Kb;
  interp->var_o = var_o;
  interp->var_1 = var_1;

  while (interp->side_pops->len < Ka + Kb)
  {
    NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
    NcGalaxyShapePop *pop = NC_GALAXY_SHAPE_POP (ncm_model_dup (NCM_MODEL (interp->pop), ser));

    g_ptr_array_add (interp->side_pops, pop);
    g_ptr_array_add (interp->side_datas, nc_galaxy_shape_pop_data_new (pop));
    ncm_serialize_free (ser);
  }

  g_array_set_size (interp->grid_a, Ka);
  g_array_set_size (interp->grid_bpop, Kb);
  g_array_set_size (interp->grid_b, Kb);
  g_array_set_size (interp->grid_ln_N, Ka * Kb);
  g_array_set_size (interp->grid_res, Ka * Kb);

  for (k = 0; k < Ka; k++)
  {
    NcGalaxyShapePop *pop = g_ptr_array_index (interp->side_pops, k);

    if (var_o >= 0)
      ncm_model_param_set (NCM_MODEL (pop), param_o, x[2 + var_o][k]);

    nc_galaxy_shape_pop_prepare (pop, g_ptr_array_index (interp->side_datas, k));
    g_array_index (interp->grid_a, gdouble, k) = nc_galaxy_shape_pop_exponent_at_origin (pop);
  }

  for (l = 0; l < Kb; l++)
  {
    NcGalaxyShapePop *pop = g_ptr_array_index (interp->side_pops, Ka + l);
    gdouble b_pop;

    if (var_1 >= 0)
      ncm_model_param_set (NCM_MODEL (pop), param_1, x[2 + var_1][l]);

    nc_galaxy_shape_pop_prepare (pop, g_ptr_array_index (interp->side_datas, Ka + l));
    b_pop = nc_galaxy_shape_pop_exponent_at_one (pop);

    g_array_index (interp->grid_bpop, gdouble, l) = b_pop;
    g_array_index (interp->grid_b, gdouble, l)    = self->axis_ratio_var ? 2.0 * b_pop + 1.0 : b_pop;
  }

  /* The constant of each pair, then the reference at the center of the values */
  {
    gdouble lo_o = 0.0, hi_o = 0.0, lo_1 = 0.0, hi_1 = 0.0;

    for (k = 0; k < Ka; k++)
    {
      if (var_o >= 0)
      {
        const gdouble val = x[2 + var_o][k];

        ncm_model_param_set (NCM_MODEL (interp->pop), param_o, val);
        lo_o = (k == 0) ? val : MIN (lo_o, val);
        hi_o = (k == 0) ? val : MAX (hi_o, val);
      }

      for (l = 0; l < Kb; l++)
      {
        gdouble f_o, f_1;

        if (var_1 >= 0)
        {
          const gdouble val = x[2 + var_1][l];

          ncm_model_param_set (NCM_MODEL (interp->pop), param_1, val);
          lo_1 = (l == 0) ? val : MIN (lo_1, val);
          hi_1 = (l == 0) ? val : MAX (hi_1, val);
        }

        nc_galaxy_shape_pop_prepare (interp->pop, interp->pop_data);
        g_array_index (interp->grid_ln_N, gdouble, k * Kb + l) =
          log (nc_galaxy_shape_pop_eval_p_reduced_factors (interp->pop, interp->pop_data, 0.0, &f_o, &f_1));
      }
    }

    if (var_o >= 0)
      ncm_model_param_set (NCM_MODEL (interp->pop), param_o, 0.5 * (lo_o + hi_o));

    if (var_1 >= 0)
      ncm_model_param_set (NCM_MODEL (interp->pop), param_1, 0.5 * (lo_1 + hi_1));

    nc_galaxy_shape_pop_prepare (interp->pop, interp->pop_data);
  }

  g_clear_pointer (&interp->wg, _ncm_spectral_polar_weight_grid_free);
  interp->wg = _ncm_spectral_polar_weight_grid_new (self->polar, 1.0, Ka, &g_array_index (interp->grid_a, gdouble, 0),
                                                    Kb, &g_array_index (interp->grid_b, gdouble, 0));

  return TRUE;
}

/*
 * The factor of the population of index k among the Ka of the origin, at the radii of the
 * variable of integration: f_o at chi (r), and for TRACE, chi = 2 r / (1 + r^2), times the
 * factor of the convention of the exponent a at the origin, 2^(a + 1) (1 + r^2)^-a.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_pop_fa (guint k, guint n, const gdouble *r, gdouble *f, gpointer user_data)
{
  NcGalaxyShapeFactorSpectralInterp *interp       = user_data;
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (interp->gsfs);
  NcGalaxyShapePop *pop                           = g_ptr_array_index (interp->side_pops, k);
  NcGalaxyShapePopData *pop_data                  = g_ptr_array_index (interp->side_datas, k);
  const gdouble a                                 = g_array_index (interp->grid_a, gdouble, k);
  const gdouble c                                 = exp2 (a + 1.0);
  guint i;

  for (i = 0; i < n; i++)
  {
    gdouble f_1;

    nc_galaxy_shape_pop_eval_p_reduced_factors (pop, pop_data, self->axis_ratio_var ? 2.0 * r[i] / (1.0 + r[i] * r[i]) : r[i], &f[i], &f_1);
  }

  if (self->axis_ratio_var)
  {
    #pragma omp simd

    for (i = 0; i < n; i++)
      f[i] *= c * exp (-a * log (1.0 + r[i] * r[i]));
  }
}

/*
 * The factor of the population of index l among the Kb of the end at 1: f_1 at chi (r), and
 * for TRACE times the factor of the convention of the exponent b of the population at 1 and
 * of the change of variable, (1 + r) (1 + r^2)^-(b + 2).
 */
static void
_nc_galaxy_shape_factor_spectral_interp_pop_fb (guint l, guint n, const gdouble *r, gdouble *f, gpointer user_data)
{
  NcGalaxyShapeFactorSpectralInterp *interp       = user_data;
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (interp->gsfs);
  NcGalaxyShapePop *pop                           = g_ptr_array_index (interp->side_pops, interp->Ka + l);
  NcGalaxyShapePopData *pop_data                  = g_ptr_array_index (interp->side_datas, interp->Ka + l);
  const gdouble b_pop                             = g_array_index (interp->grid_bpop, gdouble, l);
  guint i;

  for (i = 0; i < n; i++)
  {
    gdouble f_o;

    nc_galaxy_shape_pop_eval_p_reduced_factors (pop, pop_data, self->axis_ratio_var ? 2.0 * r[i] / (1.0 + r[i] * r[i]) : r[i], &f_o, &f[i]);
  }

  if (self->axis_ratio_var)
  {
    #pragma omp simd

    for (i = 0; i < n; i++)
      f[i] *= (1.0 + r[i]) * exp (-(b_pop + 2.0) * log (1.0 + r[i] * r[i]));
  }
}

static gdouble _nc_galaxy_shape_factor_spectral_one_minus_abs2 (const gdouble g_1, const gdouble g_2);
static void _nc_galaxy_shape_factor_spectral_seed (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k, const gdouble std_noise, const gdouble a, const gdouble b, NcGalaxyShapeFactorSpectralSeed *seed);

/*
 * ln P at the shear g, |g| < 1, for the populations that
 * _nc_galaxy_shape_factor_spectral_interp_pop_grid_setup() prepared, in interp->grid_res,
 * the value of the pair (k, l) at k Kb + l: the integral of the reference population, in
 * interp->pop, as _nc_galaxy_shape_factor_spectral_eval_ln() computes it, with the radial
 * series of its accepted panels; then the integrals of the weights of the populations
 * against those series, by the weight grid. The seed, the edges and the offset are those of
 * the reference. A shear within NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA of |g| = 1,
 * which the region of an interpolant does not reach, a failed integral and a non-finite
 * ln P abort with the point in the message.
 */
static void
_nc_galaxy_shape_factor_spectral_interp_pop_grid_shear (NcGalaxyShapeFactorSpectralInterp *interp, const gdouble g_1, const gdouble g_2)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (interp->gsfs);
  const gdouble noise_var                         = gsl_pow_2 (interp->std_noise);
  const gdouble abs_g2                            = g_1 * g_1 + g_2 * g_2;
  const gdouble abs_g                             = sqrt (abs_g2);
  const gdouble u                                 = fabs (_nc_galaxy_shape_factor_spectral_one_minus_abs2 (g_1, g_2));
  const gdouble delta                             = u / (1.0 + abs_g);
  const gdouble phase                             = (abs_g2 > 0.0) ? atan2 (g_2, g_1) : 0.0;
  const gdouble a_ref                             = nc_galaxy_shape_pop_exponent_at_origin (interp->pop);
  const gdouble b_pop                             = nc_galaxy_shape_pop_exponent_at_one (interp->pop);
  const gdouble b_ref                             = self->axis_ratio_var ? 2.0 * b_pop + 1.0 : b_pop;
  NcGalaxyShapeFactorSpectralWork **ptr           = ncm_memory_pool_get (self->work_pool);
  NcGalaxyShapeFactorSpectralWork *work           = *ptr;
  NcGalaxyShapeFactorSpectralKernel k             = {
    .pop       = interp->pop,
    .pop_data  = interp->pop_data,
    .t         = abs_g,
    .e_1       = interp->e * cos (phase),
    .e_2       = -interp->e * sin (phase),
    .c2        = 0.5 / noise_var,
    .norm      = 1.0 / (2.0 * M_PI * noise_var) / (2.0 * M_PI),
    .shift_c2  = 0.0,
    .dens_norm = exp2 (a_ref + 1.0),
    .dens_pow  = a_ref + b_pop + 2.0,
    .trace     = self->axis_ratio_var
  };
  const guint Ka      = interp->Ka;
  const guint Kb      = interp->Kb;
  gdouble *res        = &g_array_index (interp->grid_res, gdouble, 0);
  const gdouble *ln_N = &g_array_index (interp->grid_ln_N, gdouble, 0);
  NcGalaxyShapeFactorSpectralSeed seed;
  GError *error = NULL;
  gdouble L, ln_norm;
  guint q;

  if (!(abs_g2 < 1.0) || (delta < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA))
    g_error ("nc_galaxy_shape_factor_spectral: the interpolant asks for the shear g = (% 22.15g, % 22.15g), "
             "which is not inside |g| < 1 by at least %g.", g_1, g_2, NC_GALAXY_SHAPE_FACTOR_SPECTRAL_NEAR_UNIT_DELTA);

  _nc_galaxy_shape_factor_spectral_seed (self, &k, interp->std_noise, a_ref, b_ref, &seed);
  k.shift_c2 = seed.shift * k.c2;

  L = interp->std_noise * GSL_MIN (seed.J_1, seed.J_2);
  ncm_spectral_polar_peak_edges (self->polar, 1.0, seed.r_star, seed.theta_0, L, L / GSL_MAX (seed.r_star, GSL_DBL_MIN),
                                 work->r_edges, work->theta_edges);

  _ncm_spectral_polar_integrate_radial_series (self->polar, &_nc_galaxy_shape_factor_spectral_kernel, &_nc_galaxy_shape_factor_spectral_radial,
                                               &k, 1.0, a_ref, b_ref, work->r_edges, work->theta_edges, self->reltol, 0.0, NULL,
                                               work->series, work->series_coeffs, &error);

  if (error != NULL)
    g_error ("nc_galaxy_shape_factor_spectral: %s for the population %s, std_noise = % 22.15g, "
             "g = (% 22.15g, % 22.15g), eps_obs = (% 22.15g, % 22.15g).",
             error->message, G_OBJECT_TYPE_NAME (interp->pop), interp->std_noise, g_1, g_2, interp->e, 0.0);

  _ncm_spectral_polar_weight_grid_integrate (interp->wg, work->series, work->series_coeffs,
                                             &_nc_galaxy_shape_factor_spectral_interp_pop_fa, &_nc_galaxy_shape_factor_spectral_interp_pop_fb,
                                             interp, res);

  ncm_memory_pool_return (ptr);

  ln_norm = log (k.norm) - k.shift_c2;

  for (q = 0; q < Ka * Kb; q++)
  {
    res[q] = log (res[q]) + ln_N[q] + ln_norm;

    if (!gsl_finite (res[q]))
      g_error ("nc_galaxy_shape_factor_spectral: non-finite ln P for the population %s of index %u of a call of the interpolant, "
               "std_noise = % 22.15g, g = (% 22.15g, % 22.15g), eps_obs = (% 22.15g, % 22.15g).",
               G_OBJECT_TYPE_NAME (interp->pop), q, interp->std_noise, g_1, g_2, interp->e, 0.0);
  }
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
      .dens_pow  = a + b_pop + 2.0,
      .trace     = self->axis_ratio_var
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
 * ln of the integral at the real shear k->t, by _ncm_spectral_polar_integrate_radial() on the
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

  integral = _ncm_spectral_polar_integrate_radial (self->polar, &_nc_galaxy_shape_factor_spectral_kernel, &_nc_galaxy_shape_factor_spectral_radial,
                                                   k, 1.0, a, b, work->r_edges, work->theta_edges, self->reltol, 0.0, NULL, panels, error);

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
static void _nc_galaxy_shape_factor_spectral_var_radial (gpointer user_data, guint n, const gdouble *v, gdouble *rho);

/* What the integrand of the variance needs */
typedef struct _NcGalaxyShapeFactorSpectralVar
{
  const NcGalaxyShapeFactorSpectralKernel *k;
  gdouble delta; /* 1 - t */
  gdouble *r;    /* the radii r = 1 - v of a call of the radial factor, in the workspace */
} NcGalaxyShapeFactorSpectralVar;

/*
 * ln of A (t) = Phi (t) + s (t) Lap Phi (t) / 2 at t = k->t = 1 - delta, the marginal to
 * second order in the distance of the sheared ellipticity to its mean. w = m_t (z) is the
 * TRACE_DET map of the variable of integration z, which has the mean t and the variance
 * s (t) per component,
 *   s (t) = (1 - t^2)^2 / 2 INT_0^1 (density of r) r^2 / (1 - t^2 r^2) dr,
 * and Phi (w) = N_2 (eps_obs - C (w)), C the identity for TRACE_DET and T for TRACE. The
 * integral is computed in v = 1 - r by _ncm_spectral_polar_integrate_radial(), where the factor
 * 1 / (delta + t v) of 1 / (1 - t^2 r^2) varies on the scale delta, with radial edges
 * placed for that scale at v = 0. The integrand depends only on r, so all of it is the
 * radial factor, computed once per radius. Differences to 1 are formed from delta.
 */
static gdouble
_nc_galaxy_shape_factor_spectral_second_order_ln (NcGalaxyShapeFactorSpectralPrivate * const self, const NcGalaxyShapeFactorSpectralKernel *k,
                                                  const gdouble std_noise, const gdouble a, const gdouble b, const gdouble delta,
                                                  NcGalaxyShapeFactorSpectralWork *work, GError **error)
{
  const gdouble t                    = k->t;
  const gdouble noise_var            = gsl_pow_2 (std_noise);
  const gdouble one_minus_t2         = delta * (1.0 + t);
  NcGalaxyShapeFactorSpectralVar var = {k, delta, NULL};
  gdouble s, d, lap_ratio;

  /* A call of the radial factor has at most 2^level-cap + 1 radii */
  g_array_set_size (work->var_r, (1u << ncm_spectral_polar_get_level_cap (self->polar)) + 1);
  var.r = &g_array_index (work->var_r, gdouble, 0);

  ncm_spectral_polar_peak_edges (self->polar, 1.0, 0.0, 0.0, delta, M_PI, work->r_edges, work->theta_edges);
  s = 0.5 * gsl_pow_2 (one_minus_t2) *
      _ncm_spectral_polar_integrate_radial (self->polar, &_nc_galaxy_shape_factor_spectral_var_kernel, &_nc_galaxy_shape_factor_spectral_var_radial,
                                            &var, 1.0, b, a, work->r_edges, work->theta_edges, self->reltol, 0.0, NULL, NULL, error);

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
 * The integrand of the variance for _ncm_spectral_polar_integrate_radial() in v = 1 - r, with
 * the weight v^b (1 - v)^a, is (reduced density of r) r^2 / ((delta + t v) (1 + t r)) / (2 pi):
 * it does not depend on the angle, and all of it is the radial factor
 * (_nc_galaxy_shape_factor_spectral_var_radial()). The rest is 1.
 */
static void
_nc_galaxy_shape_factor_spectral_var_kernel (gpointer user_data, gdouble v, guint n, const gdouble *theta,
                                             const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  guint j;

  for (j = 0; j < n; j++)
    fval[j] = 1.0;
}

/*
 * The integrand of the variance at the n values v[i] of v = 1 - r, in rho[i]. All of it
 * is a radial factor because the population is isotropic: its density depends on r alone.
 */
static void
_nc_galaxy_shape_factor_spectral_var_radial (gpointer user_data, guint n, const gdouble *v, gdouble *rho)
{
  const NcGalaxyShapeFactorSpectralVar *var = user_data;
  const gdouble t                           = var->k->t;
  const gdouble delta                       = var->delta;
  gdouble *r                                = var->r;
  guint i;

  for (i = 0; i < n; i++)
    r[i] = 1.0 - v[i];

  _nc_galaxy_shape_factor_spectral_density (var->k, n, r, rho);

  for (i = 0; i < n; i++)
    rho[i] *= r[i] * r[i] / ((delta + t * v[i]) * (1.0 + t * r[i])) / (2.0 * M_PI);
}

static void _nc_galaxy_shape_factor_spectral_h (const NcGalaxyShapeFactorSpectralKernel *k, const gdouble a, const gdouble b, const guint n, const gdouble *r, const gdouble cos_0, const gdouble sin_0, gdouble *h, gdouble *d, gdouble *x);
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
  gdouble scan_r[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  gdouble scan_h[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  gdouble scan_d[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  gdouble scan_x[NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_UNIFORM + NC_GALAXY_SHAPE_FACTOR_SPECTRAL_SCAN_EXTRA];
  guint n_scan     = 0;
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

    if (r_u <= r_e)
      iu++;
    else
      ie++;

    if (!(r > 0.0) || !(r < 1.0) || ((n_scan > 0) && (r == scan_r[n_scan - 1])))
      continue;

    scan_r[n_scan++] = r;
  }

  _nc_galaxy_shape_factor_spectral_h (k, a, b, n_scan, scan_r, cos_0, sin_0, scan_h, scan_d, scan_x);

  for (i = 0; i < n_scan; i++)
  {
    const gdouble r = scan_r[i];
    const gdouble h = scan_h[i];

    shift = GSL_MIN (shift, scan_d[i]);

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
    gdouble d, xs, h1, h2;

    _nc_galaxy_shape_factor_spectral_h (k, a, b, 1, &x1, cos_0, sin_0, &h1, &d, &xs);
    _nc_galaxy_shape_factor_spectral_h (k, a, b, 1, &x2, cos_0, sin_0, &h2, &d, &xs);

    for (i = 0; i < NC_GALAXY_SHAPE_FACTOR_SPECTRAL_GOLDEN_STEPS; i++)
    {
      if (h1 > h2)
      {
        hi = x2;
        x2 = x1;
        h2 = h1;
        x1 = hi - gr * (hi - lo);
        _nc_galaxy_shape_factor_spectral_h (k, a, b, 1, &x1, cos_0, sin_0, &h1, &d, &xs);
      }
      else
      {
        lo = x1;
        x1 = x2;
        h1 = h2;
        x2 = lo + gr * (hi - lo);
        _nc_galaxy_shape_factor_spectral_h (k, a, b, 1, &x2, cos_0, sin_0, &h2, &d, &xs);
      }
    }

    seed->r_star = 0.5 * (lo + hi);
  }

  seed->shift = shift;
}

/* log1p (x) for x > -1, in a loop that vectorizes: Kahan's log (u) x / (u - 1), u = 1 + x */
static inline gdouble
_nc_galaxy_shape_factor_spectral_log1p (const gdouble x)
{
  const gdouble u = 1.0 + x;

  return ncm_vmath_select (u == 1.0, x, log (u) * (x / (u - 1.0)));
}

/*
 * h (r) = ln (density of r) - d / (2 std_noise^2) at the n radii r[i], 0 < r < 1, of the
 * ray of direction (cos_0, sin_0), in h[i], with d the squared distance of eps_obs to the
 * sheared ellipticity, in d[i]. The logarithm of the density is the two end power laws,
 * of exponents a and b, plus the logarithm of the density of the convention divided by
 * them (in x, scratch). Both pieces of the convention are computed at the n points at once.
 */
static void
_nc_galaxy_shape_factor_spectral_h (const NcGalaxyShapeFactorSpectralKernel *k, const gdouble a, const gdouble b, const guint n,
                                    const gdouble *r, const gdouble cos_0, const gdouble sin_0, gdouble *h, gdouble *d, gdouble *x)
{
  const gdouble c2 = k->c2;
  guint i;

  _nc_galaxy_shape_factor_spectral_dist2 (k, n, r, cos_0, sin_0, d);
  _nc_galaxy_shape_factor_spectral_density (k, n, r, x);

  #pragma omp simd

  for (i = 0; i < n; i++)
    h[i] = a * log (r[i]) + b * _nc_galaxy_shape_factor_spectral_log1p (-r[i]) + log (x[i]) - d[i] * c2;
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
 * The map of the convention at the real shear t, in the variable of integration z = x + i y:
 * T (m_t (z)) for TRACE (_nc_galaxy_shape_factor_spectral_map_trace()) and m_t (z) for
 * TRACE_DET, m_t the TRACE_DET map. Its callers pass trace as a constant, so that each loop
 * that calls it holds a single map once inlined, and vectorizes.
 */
NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE void
_nc_galaxy_shape_factor_spectral_map (const gboolean trace, const gdouble t, const gdouble x, const gdouble y,
                                      gdouble * restrict f_x, gdouble * restrict f_y)
{
  if (trace)
    _nc_galaxy_shape_factor_spectral_map_trace (t, x, y, f_x, f_y);
  else
    nc_wl_ellipticity_apply_shear_trace_det_real (t, x, y, f_x, f_y);
}

/* The loop of _nc_galaxy_shape_factor_spectral_kernel() for one convention */
NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE void
_nc_galaxy_shape_factor_spectral_kernel_loop (const NcGalaxyShapeFactorSpectralKernel *k, const gboolean trace, const gdouble r, const guint n,
                                              const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  const gdouble t        = k->t;
  const gdouble e_1      = k->e_1;
  const gdouble e_2      = k->e_2;
  const gdouble shift_c2 = k->shift_c2;
  const gdouble c2       = k->c2;
  guint j;

  /* The fields of k in locals: fval might alias them, which keeps the loop scalar */
  #pragma omp simd

  for (j = 0; j < n; j++)
  {
    gdouble f_x, f_y;

    _nc_galaxy_shape_factor_spectral_map (trace, t, r * cos_theta[j], r * sin_theta[j], &f_x, &f_y);
    fval[j] = exp (shift_c2 - (gsl_pow_2 (e_1 - f_x) + gsl_pow_2 (e_2 - f_y)) * c2);
  }
}

/*
 * The integrand of _ncm_spectral_polar_integrate_radial(), in the polar coordinates of the
 * variable of integration z: G (r, theta) = (reduced density of r = |z|)
 * N_2 (eps_obs - F_t (z)) / (2 pi), F_t the map of the convention
 * (_nc_galaxy_shape_factor_spectral_map()), with the offset of the kernel taken out of the
 * exponent of N_2. The factor that depends only on r, the reduced density times norm, is
 * _nc_galaxy_shape_factor_spectral_radial(): here only the exponential.
 */
static void
_nc_galaxy_shape_factor_spectral_kernel (gpointer user_data, gdouble r, guint n, const gdouble *theta,
                                         const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  const NcGalaxyShapeFactorSpectralKernel *k = user_data;

  if (k->trace)
    _nc_galaxy_shape_factor_spectral_kernel_loop (k, TRUE, r, n, cos_theta, sin_theta, fval);
  else
    _nc_galaxy_shape_factor_spectral_kernel_loop (k, FALSE, r, n, cos_theta, sin_theta, fval);
}

/*
 * The factor of the integrand that depends only on r, given to
 * _ncm_spectral_polar_integrate_radial() so that it is computed once per radius: norm times
 * the density of the convention (_nc_galaxy_shape_factor_spectral_density()).
 * It relies on an isotropic population, whose density depends on r alone.
 */
static void
_nc_galaxy_shape_factor_spectral_radial (gpointer user_data, guint n, const gdouble *r, gdouble *rho)
{
  const NcGalaxyShapeFactorSpectralKernel *k = user_data;
  guint i;

  _nc_galaxy_shape_factor_spectral_density (k, n, r, rho);

  for (i = 0; i < n; i++)
    rho[i] *= k->norm;
}

/*
 * The TRACE map at the real shear t in the variable e_I = x + i y: T (m_t (e_I)), with
 * m_t the TRACE_DET map. It equals f_t (T (e_I)), f_t the TRACE map. With
 * m_t (e_I) = N / D, N = (x + t) + i y and D = (1 + t x) + i t y, and T (m) = 2 m / (1 + |m|^2),
 * it is 2 N conj (D) / (|D|^2 + |N|^2): one division instead of the two of m_t followed by
 * T, the slowest operations of the loop of the kernel.
 */
static void
_nc_galaxy_shape_factor_spectral_map_trace (gdouble t, gdouble x, gdouble y, gdouble * restrict f_x, gdouble * restrict f_y)
{
  const gdouble num_x = x + t;
  const gdouble den_x = 1.0 + t * x;
  const gdouble den_y = t * y;
  const gdouble fac   = 2.0 / ((num_x * num_x + y * y) + (den_x * den_x + den_y * den_y));

  *f_x = (num_x * den_x + y * den_y) * fac;
  *f_y = (y * den_x - num_x * den_y) * fac;
}

/* The inverse of the map above at the shear g: m_g^{-1} (T^{-1} (e_obs)), for |e_obs| <= 1 */
static complex double
_nc_galaxy_shape_factor_spectral_map_inv_trace (complex double g, complex double e_obs)
{
  const gdouble abs_e2 = gsl_pow_2 (creal (e_obs)) + gsl_pow_2 (cimag (e_obs));

  return nc_wl_ellipticity_apply_shear_inv_trace_det (g, e_obs / (1.0 + sqrt (1.0 - abs_e2)));
}

/*
 * The density of r = |z| divided by its two end power laws, at the n radii r[i] in out[i].
 * TRACE_DET: the reduced density of the population at r = |chi_I|. TRACE: the density of
 * r = |e_I| divided by r^alpha_o (1 - r)^(2 alpha_1 + 1), for a population that is a
 * density of |chi_I| = 2 r / (1 + r^2) with the exponents alpha_o and alpha_1,
 * 2^(alpha_o + 1) (1 + r) (1 + r^2)^-(alpha_o + alpha_1 + 2) times the reduced density of
 * the population at |chi_I|, the factor of the convention applied in a loop that
 * vectorizes.
 */
static void
_nc_galaxy_shape_factor_spectral_density (const NcGalaxyShapeFactorSpectralKernel *k, guint n, const gdouble *r, gdouble *out)
{
  guint i;

  if (k->trace)
  {
    const gdouble c  = k->dens_norm;
    const gdouble mp = -k->dens_pow;

    for (i = 0; i < n; i++)
      out[i] = nc_galaxy_shape_pop_eval_p_reduced (k->pop, k->pop_data, 2.0 * r[i] / (1.0 + r[i] * r[i]));

    #pragma omp simd

    for (i = 0; i < n; i++)
      out[i] *= c * (1.0 + r[i]) * pow (1.0 + r[i] * r[i], mp);
  }
  else
  {
    for (i = 0; i < n; i++)
      out[i] = nc_galaxy_shape_pop_eval_p_reduced (k->pop, k->pop_data, r[i]);
  }
}

/* The loop of _nc_galaxy_shape_factor_spectral_dist2() for one convention */
NC_GALAXY_SHAPE_FACTOR_SPECTRAL_INLINE void
_nc_galaxy_shape_factor_spectral_dist2_loop (const NcGalaxyShapeFactorSpectralKernel *k, const gboolean trace, const guint n, const gdouble *r,
                                             const gdouble cos_0, const gdouble sin_0, gdouble *d)
{
  const gdouble t   = k->t;
  const gdouble e_1 = k->e_1;
  const gdouble e_2 = k->e_2;
  guint i;

  #pragma omp simd

  for (i = 0; i < n; i++)
  {
    gdouble f_x, f_y;

    _nc_galaxy_shape_factor_spectral_map (trace, t, r[i] * cos_0, r[i] * sin_0, &f_x, &f_y);
    d[i] = gsl_pow_2 (e_1 - f_x) + gsl_pow_2 (e_2 - f_y);
  }
}

/*
 * The squared distance of eps_obs to the map of the convention of z = r (cos_0, sin_0)
 * (_nc_galaxy_shape_factor_spectral_map()), at the n radii r[i] in d[i]
 */
static void
_nc_galaxy_shape_factor_spectral_dist2 (const NcGalaxyShapeFactorSpectralKernel *k, guint n, const gdouble *r,
                                        gdouble cos_0, gdouble sin_0, gdouble *d)
{
  if (k->trace)
    _nc_galaxy_shape_factor_spectral_dist2_loop (k, TRUE, n, r, cos_0, sin_0, d);
  else
    _nc_galaxy_shape_factor_spectral_dist2_loop (k, FALSE, n, r, cos_0, sin_0, d);
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
 * nc_galaxy_shape_factor_spectral_data_take_interp:
 * @data: a #NcGalaxyShapeFactorData of a #NcGalaxyShapeFactorSpectral
 *
 * To store the interpolant of a galaxy without keeping it in memory, takes the interpolant
 * that nc_galaxy_shape_factor_spectral_build_interp() built for the row of @data and leaves
 * the galaxy without one: until the next build, the galaxy is evaluated by the integrals.
 *
 * Returns: (transfer full) (nullable): the interpolant, or %NULL when the galaxy holds none
 * built for its row.
 */
NcmSpectralTensor *
nc_galaxy_shape_factor_spectral_data_take_interp (NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorSpectralData *ldata = data->ldata;
  NcmSpectralTensor *st;

  if (ldata->loaded || !ldata->valid || (ldata->epsilon_obs_1 != data->epsilon_obs_1) || (ldata->epsilon_obs_2 != data->epsilon_obs_2) ||
      (ldata->std_noise != data->std_noise))
    return NULL;

  st           = ldata->st;
  ldata->st    = NULL;
  ldata->valid = FALSE;

  return st;
}

/**
 * nc_galaxy_shape_factor_spectral_set_obs_source:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @pop: the #NcGalaxyShapePop that built the interpolants
 * @obs: a #NcGalaxyWLObs
 *
 * Sets the record of the interpolants written into @obs by
 * nc_galaxy_shape_factor_data_write_row() (nc_galaxy_wl_obs_set_shape_interp_source()):
 * copies of @gsfs under `shape-factor` and of @pop with its parameter values under
 * `shape-pop`. Call it after the last row is written. A fit that reads the rows compares its
 * shape factor and population with the record in nc_galaxy_shape_factor_data_prepare().
 */
void
nc_galaxy_shape_factor_spectral_set_obs_source (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyWLObs *obs)
{
  NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmObjDictStr *source = ncm_obj_dict_str_new ();
  GObject *gsfs_dup     = ncm_serialize_dup_obj (ser, G_OBJECT (gsfs));
  GObject *pop_dup;

  ncm_serialize_reset (ser, FALSE);
  pop_dup = ncm_serialize_dup_obj (ser, G_OBJECT (pop));

  ncm_obj_dict_str_add (source, "shape-factor", gsfs_dup);
  ncm_obj_dict_str_add (source, "shape-pop", pop_dup);
  nc_galaxy_wl_obs_set_shape_interp_source (obs, source);

  g_object_unref (gsfs_dup);
  g_object_unref (pop_dup);
  ncm_obj_dict_str_unref (source);
  ncm_serialize_free (ser);
}

/**
 * nc_galaxy_shape_factor_spectral_build_interp:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @mset: a #NcmMSet
 * @data: a #NcGalaxyShapeFactorData
 * @z_max: the largest source redshift of the galaxy
 *
 * Builds the interpolant of $\ln P$ in $g$ of the galaxy of @data for its observed
 * ellipticity and noise, with the population of @mset and the interp-* settings of @gsfs,
 * and records what it holds for. It replaces the interpolant the galaxy held, built or
 * loaded; a loaded one stays unchanged in the catalog. With
 * #NcGalaxyShapeFactorSpectral:restrict-range the range is bounded with @z_max (see that
 * property). nc_galaxy_shape_factor_prepare() must have run on @mset and the population of
 * @data must be resolved (nc_galaxy_shape_factor_prepare_data_array()). Builds of different
 * galaxies may run in parallel. A build that reaches the maximum depth aborts with the
 * galaxy in the message.
 */
void
nc_galaxy_shape_factor_spectral_build_interp (NcGalaxyShapeFactorSpectral *gsfs, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);
  NcGalaxyShapeFactorSpectralData *ldata          = data->ldata;
  NcGalaxyShapePop *pop                           = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

  if (pop == NULL)
    g_error ("nc_galaxy_shape_factor_spectral_build_interp: the model set has no NcGalaxyShapePop.");

  if (ldata->loaded)
  {
    ncm_spectral_tensor_clear (&ldata->st);
    g_clear_object (&ldata->src_gsfs);
    nc_galaxy_shape_pop_clear (&ldata->src_pop);
    ldata->loaded = FALSE;
  }

  ldata->g_target = _nc_galaxy_shape_factor_spectral_reach (self, data, z_max);
  _nc_galaxy_shape_factor_spectral_interp_build (gsfs, pop, data, data->epsilon_obs_1, data->epsilon_obs_2);
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

  if (n < 1)
    g_error ("nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels: interp-min-radial-panels must be at least 1.");

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
 * nc_galaxy_shape_factor_spectral_set_interp_trim_frac:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @f: fraction of interp-abstol given to the removed coefficients, in $[0, 1)$
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-trim-frac. A value outside $[0, 1)$ aborts.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_trim_frac (NcGalaxyShapeFactorSpectral *gsfs, const gdouble f)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  ncm_spectral_tensor_builder_set_trim_frac (self->builder_interp, f);

  self->interp_trim_frac = f;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_trim_frac:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-trim-frac.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_interp_trim_frac (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_trim_frac;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_compress_frac:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @f: fraction of interp-abstol given to the compression, in $[0, 1]$
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-compress-frac. A value outside $[0, 1]$ aborts.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_compress_frac (NcGalaxyShapeFactorSpectral *gsfs, const gdouble f)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  if (!((f >= 0.0) && (f <= 1.0)))
    g_error ("nc_galaxy_shape_factor_spectral_set_interp_compress_frac: %g is outside [0, 1].", f);

  self->interp_compress_frac = f;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_compress_frac:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-compress-frac.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_interp_compress_frac (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_compress_frac;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_pop_level_min:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @level_min: the starting level of the population variables, from 1 to 10
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-pop-level-min.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_pop_level_min (NcGalaxyShapeFactorSpectral *gsfs, const guint level_min)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->interp_pop_level_min = level_min;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_pop_level_min:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-pop-level-min.
 */
guint
nc_galaxy_shape_factor_spectral_get_interp_pop_level_min (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_pop_level_min;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @level_cap: the largest level of the population variables, from 1 to 10
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-pop-level-cap.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (NcGalaxyShapeFactorSpectral *gsfs, const guint level_cap)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->interp_pop_level_cap = level_cap;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_pop_level_cap:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-pop-level-cap.
 */
guint
nc_galaxy_shape_factor_spectral_get_interp_pop_level_cap (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_pop_level_cap;
}

/**
 * nc_galaxy_shape_factor_spectral_set_interp_pop_one_integral:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @one_integral: whether the populations of a call share one integral per shear
 *
 * Sets #NcGalaxyShapeFactorSpectral:interp-pop-one-integral.
 */
void
nc_galaxy_shape_factor_spectral_set_interp_pop_one_integral (NcGalaxyShapeFactorSpectral *gsfs, const gboolean one_integral)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->interp_pop_one_integral = one_integral;
}

/**
 * nc_galaxy_shape_factor_spectral_get_interp_pop_one_integral:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:interp-pop-one-integral.
 */
gboolean
nc_galaxy_shape_factor_spectral_get_interp_pop_one_integral (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->interp_pop_one_integral;
}

/**
 * nc_galaxy_shape_factor_spectral_set_restrict_range:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @restrict_range: whether the interpolant covers only the shears a galaxy reaches
 *
 * Sets #NcGalaxyShapeFactorSpectral:restrict-range.
 */
void
nc_galaxy_shape_factor_spectral_set_restrict_range (NcGalaxyShapeFactorSpectral *gsfs, const gboolean restrict_range)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  self->restrict_range = restrict_range;
}

/**
 * nc_galaxy_shape_factor_spectral_get_restrict_range:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:restrict-range.
 */
gboolean
nc_galaxy_shape_factor_spectral_get_restrict_range (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->restrict_range;
}

/**
 * nc_galaxy_shape_factor_spectral_set_line_abstol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 * @abstol: absolute tolerance on $\ln P$ of the table along a line, in $[0, 1]$
 *
 * Sets #NcGalaxyShapeFactorSpectral:line-abstol. A value outside $[0, 1]$ aborts.
 */
void
nc_galaxy_shape_factor_spectral_set_line_abstol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble abstol)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  if (!((abstol >= 0.0) && (abstol <= 1.0)))
    g_error ("nc_galaxy_shape_factor_spectral_set_line_abstol: %g is outside [0, 1].", abstol);

  self->line_abstol = abstol;
}

/**
 * nc_galaxy_shape_factor_spectral_get_line_abstol:
 * @gsfs: a #NcGalaxyShapeFactorSpectral
 *
 * Returns: #NcGalaxyShapeFactorSpectral:line-abstol.
 */
gdouble
nc_galaxy_shape_factor_spectral_get_line_abstol (NcGalaxyShapeFactorSpectral *gsfs)
{
  NcGalaxyShapeFactorSpectralPrivate * const self = nc_galaxy_shape_factor_spectral_get_instance_private (gsfs);

  return self->line_abstol;
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

