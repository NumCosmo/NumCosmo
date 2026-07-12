/***************************************************************************
 *            nc_galaxy_shape_factor_knots.c
 *
 *  Fri Jul 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_knots.c
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
 * NcGalaxyShapeFactorKnots:
 *
 * Fixed-node cubature evaluation of the intrinsic-ellipticity marginal.
 *
 * A playground sibling of #NcGalaxyShapeFactorQuad: it evaluates the exact
 * same integral,
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - f_g(\chi_I);
 *   \sigma_\mathrm{noise}^2\big),$$
 * over the same plane-substituted, re-centered box (see #NcGalaxyShapeFactorQuad's
 * documentation for the $(u,v)$ substitution and the re-centering rationale,
 * both reused verbatim here), but with a FIXED set of $(u,v)$ nodes and
 * weights instead of Cuhre's adaptive subdivision. The knots are generated
 * once from #NcGalaxyShapeFactorKnots:bound, #NcGalaxyShapeFactorKnots:n and
 * #NcGalaxyShapeFactorKnots:method (whenever any of the three changes), then
 * simply translated to the per-galaxy re-centering point $(u_0,v_0)$ and
 * summed at every evaluation -- no refinement, no error estimate.
 *
 * This exists to let different fixed-knot schemes be tried and inspected
 * directly against the Cuhre reference (accuracy, node count needed, cost)
 * without touching the validated Quad backend. #NcGalaxyShapeFactorKnotsMethod
 * has two schemes: @NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN, a plain
 * tensor-product Gauss-Legendre grid laid out once on $[-B,B]^2$ and simply
 * translated to $(u_0,v_0)$ at every evaluation (no notion of scale); and
 * @NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN_ADAPTIVE, which keeps the
 * same $n\times n$ node COUNT but re-centers on the joint mode of
 * $P_\mathrm{pop}(\chi_I)\,N_2(\epsilon_\mathrm{obs}-f_g(\chi_I))$ (found by
 * nc_galaxy_shape_intrinsic_mode_find(), a nested search exploiting that
 * $P_\mathrm{pop}$ has no angular dependence, see that file's documentation)
 * and rescales the box to a few widths of the local curvature there, capped
 * at $B$. This matters because the fixed scheme's uniform $[-B,B]^2$ grid
 * (generous enough to cover the whole disc for a broad population) puts
 * almost none of its $n^2$ nodes near an isolated, narrow, off-center peak --
 * the adaptive scheme spends the same node budget where the integrand
 * actually lives instead. Adding another scheme is a matter of a new enum
 * value plus a matching case in `_nc_galaxy_shape_factor_knots_regenerate()`
 * and `_nc_galaxy_shape_factor_knots_eval()`.
 * nc_galaxy_shape_factor_knots_peek_nodes() and
 * nc_galaxy_shape_factor_knots_peek_weights() expose the generated 1D
 * nodes/weights directly for inspection (e.g. from Python, as numpy array
 * views) without needing to evaluate anything; for the adaptive scheme these
 * are the REFERENCE nodes on $[-1,1]$ before the per-evaluation re-centering
 * and rescaling (there is no single fixed $(u,v)$ layout to expose).
 *
 * Being untested against sub-percent accuracy claims, this class is NOT a
 * recommended likelihood backend -- use #NcGalaxyShapeFactorQuad for that.
 * It is meant to be poked at.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_knots.h"
#include "nc/lss/galaxy/nc_galaxy_shape_intrinsic_mode.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"
#include "ncm/algebra/ncm_vector.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define NC_GALAXY_SHAPE_FACTOR_KNOTS_DEFAULT_BOUND  8.0
#define NC_GALAXY_SHAPE_FACTOR_KNOTS_DEFAULT_N      32
#define NC_GALAXY_SHAPE_FACTOR_KNOTS_DEFAULT_METHOD NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN

struct _NcGalaxyShapeFactorKnots
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorKnotsPrivate
{
  /* Resolved once at construction, see #NcGalaxyShapeFactorQuad's private
   * struct for why. */
  complex double (*apply_shear) (complex double g, complex double chi);
  complex double (*apply_shear_inv) (complex double g, complex double eps_obs);

  gdouble bound;
  guint n;
  NcGalaxyShapeFactorKnotsMethod method;

  /* 1D nodes/weights on [-bound, bound]; the 2D rule is their tensor
   * product, regenerated whenever bound/n/method change. */
  NcmVector *nodes;
  NcmVector *weights;
} NcGalaxyShapeFactorKnotsPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorKnots, nc_galaxy_shape_factor_knots, NC_TYPE_GALAXY_SHAPE_FACTOR)

enum
{
  PROP_0,
  PROP_BOUND,
  PROP_N,
  PROP_METHOD,
  PROP_LEN,
};

/* Knot generation ------------------------------------------------------- */

static void
_nc_galaxy_shape_factor_knots_gl_cartesian (const gdouble bound, const guint n, gdouble *nodes, gdouble *weights)
{
  gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc (n);
  guint i;

  for (i = 0; i < n; i++)
    gsl_integration_glfixed_point (-bound, bound, i, &nodes[i], &weights[i], table);

  gsl_integration_glfixed_table_free (table);
}

static void
_nc_galaxy_shape_factor_knots_regenerate (NcGalaxyShapeFactorKnots *gsfk)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);
  gdouble *nodes                               = g_new (gdouble, self->n);
  gdouble *weights                             = g_new (gdouble, self->n);

  switch (self->method)
  {
    case NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN:
      _nc_galaxy_shape_factor_knots_gl_cartesian (self->bound, self->n, nodes, weights);
      break;
    case NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN_ADAPTIVE:
      /* Reference rule on [-1,1]: rescaled per evaluation (see
       * _nc_galaxy_shape_factor_knots_eval()), so bound plays no role here
       * beyond capping that per-evaluation scale. */
      _nc_galaxy_shape_factor_knots_gl_cartesian (1.0, self->n, nodes, weights);
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }

  ncm_vector_clear (&self->nodes);
  ncm_vector_clear (&self->weights);
  self->nodes   = ncm_vector_new_data_malloc (nodes, self->n, 1);
  self->weights = ncm_vector_new_data_malloc (weights, self->n, 1);
}

/* Integrand (shared physics with NcGalaxyShapeFactorQuad) --------------- */

static gdouble
_nc_galaxy_shape_factor_knots_integrand (const gdouble u, const gdouble v, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, complex double (*apply_shear) (complex double, complex double), const complex double g, const complex double eps_obs, const gdouble std_noise)
{
  const gdouble rho2            = u * u + v * v;
  const gdouble s2              = 1.0 + rho2;
  const gdouble jac             = 1.0 / gsl_pow_2 (s2);
  const complex double chi_i    = (u + I * v) / sqrt (s2);
  const gdouble p_pop           = nc_galaxy_shape_pop_eval_p_rho2 (pop, pop_data, rho2) / M_PI;
  const complex double eps_true = apply_shear (g, chi_i);
  const complex double delta    = eps_obs - eps_true;
  const gdouble noise_var       = gsl_pow_2 (std_noise);
  const gdouble noise           = exp (-(gsl_pow_2 (creal (delta)) + gsl_pow_2 (cimag (delta))) / (2.0 * noise_var)) / (2.0 * M_PI * noise_var);

  return p_pop * noise * jac;
}

/* Local curvature of the (u,v)-space log-integrand at the re-centering point,
 * used only by the adaptive scheme to pick a per-galaxy box half-width: the
 * SMALLER of the two axis curvatures sets the slower-decaying direction, so
 * scaling by it (not the larger) keeps that direction adequately covered
 * too. Falls back to bound_cap whenever the curvature estimate is degenerate
 * (e.g. a genuinely flat/broad population, or numerical underflow at the
 * center itself), which recovers the classic fixed-box behaviour. */
static gdouble
_nc_galaxy_shape_factor_knots_adaptive_scale (const gdouble u0, const gdouble v0,
                                              NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                              complex double (*apply_shear) (complex double, complex double),
                                              const complex double g, const complex double eps_obs,
                                              const gdouble std_noise, const gdouble bound_cap)
{
  const gdouble h       = 1.0e-3;
  const gdouble k_sigma = 8.0;
  gdouble f00, fpu, fmu, fpv, fmv, huu, hvv, curv_min;

  #define LNI(uu, vv)                                                                                                             \
          ({                                                                                                                      \
    const gdouble _val = _nc_galaxy_shape_factor_knots_integrand ((uu), (vv), pop, pop_data, apply_shear, g, eps_obs, std_noise); \
    ((_val > 0.0) && isfinite (_val)) ? log (_val) : -G_MAXDOUBLE;                                                                \
  })
  f00 = LNI (u0, v0);
  fpu = LNI (u0 + h, v0);
  fmu = LNI (u0 - h, v0);
  fpv = LNI (u0, v0 + h);
  fmv = LNI (u0, v0 - h);
  #undef LNI

  huu      = (fpu - 2.0 * f00 + fmu) / gsl_pow_2 (h);
  hvv      = (fpv - 2.0 * f00 + fmv) / gsl_pow_2 (h);
  curv_min = MIN (-huu, -hvv);

  if (!(curv_min > 0.0) || !isfinite (curv_min))
    return bound_cap;

  return CLAMP (k_sigma / sqrt (curv_min), bound_cap * 1.0e-3, bound_cap);
}

static void
_nc_galaxy_shape_factor_knots_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_factor_knots_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_knots_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  /* No persistent per-galaxy state: each evaluation is a self-contained
   * weighted sum over the precomputed knots. */
  data->ldata                  = NULL;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_knots_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_knots_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_knots_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_knots_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
}

static gdouble
_nc_galaxy_shape_factor_knots_eval (NcGalaxyShapeFactorKnots *gsfk, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);
  const complex double g                       = g_1 + I * g_2;
  const complex double eps_obs                 = epsilon_obs_1 + I * epsilon_obs_2;
  gdouble u0, v0, scale, result;
  guint i, j;

  if (self->method == NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN_ADAPTIVE)
  {
    /* Re-center on the (u,v) pre-image of the TRUE joint mode of
     * P_pop(chi_I) N_2(eps_obs-f_g(chi_I)) (not just the naive noiseless
     * inverse map): chi_I lives in exactly the same disc this class already
     * integrates over via the (u,v) compactification, so no extra forward/
     * inverse shear map is needed to place it, unlike Quad's chi_L-space
     * third hint. */
    NcGalaxyShapeIntrinsicMode mode;
    complex double chi_i_mode;
    gdouble chi_abs, chi_abs_clamped, map_scale;
    complex double chi_centered;

    nc_galaxy_shape_intrinsic_mode_find (self->apply_shear, self->apply_shear_inv,
                                         pop, data->pop_data, g, eps_obs, data->std_noise, &mode);
    chi_i_mode      = mode.rho * cexp (I * mode.theta);
    chi_abs         = cabs (chi_i_mode);
    chi_abs_clamped = MIN (chi_abs, 0.999);
    chi_centered    = (chi_abs > 0.0) ? chi_i_mode * (chi_abs_clamped / chi_abs) : chi_i_mode;
    map_scale       = 1.0 / sqrt (1.0 - gsl_pow_2 (chi_abs_clamped));
    u0              = creal (chi_centered) * map_scale;
    v0              = cimag (chi_centered) * map_scale;

    scale = _nc_galaxy_shape_factor_knots_adaptive_scale (u0, v0, pop, data->pop_data, self->apply_shear,
                                                          g, eps_obs, data->std_noise, self->bound);
  }
  else
  {
    /* Re-center on the (u,v) pre-image of mu = f_g^{-1}(eps_obs), the naive
     * noiseless intrinsic estimate (ignores the population entirely). */
    const complex double mu        = self->apply_shear_inv (g, eps_obs);
    const gdouble mu_abs           = cabs (mu);
    const gdouble mu_abs_clamped   = MIN (mu_abs, 0.999);
    const complex double mu_center = (mu_abs > 0.0) ? mu * (mu_abs_clamped / mu_abs) : mu;
    const gdouble map_scale        = 1.0 / sqrt (1.0 - gsl_pow_2 (mu_abs_clamped));

    u0    = creal (mu_center) * map_scale;
    v0    = cimag (mu_center) * map_scale;
    scale = 1.0; /* nodes already laid out on [-bound,bound] */
  }

  result = 0.0;

  for (i = 0; i < self->n; i++)
  {
    const gdouble u  = u0 + scale * ncm_vector_fast_get (self->nodes, i);
    const gdouble wu = scale * ncm_vector_fast_get (self->weights, i);

    for (j = 0; j < self->n; j++)
    {
      const gdouble v  = v0 + scale * ncm_vector_fast_get (self->nodes, j);
      const gdouble wv = scale * ncm_vector_fast_get (self->weights, j);

      result += wu * wv * _nc_galaxy_shape_factor_knots_integrand (u, v, pop, data->pop_data, self->apply_shear, g, eps_obs, data->std_noise);
    }
  }

  return result;
}

static gdouble
_nc_galaxy_shape_factor_knots_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return _nc_galaxy_shape_factor_knots_eval (NC_GALAXY_SHAPE_FACTOR_KNOTS (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

static gdouble
_nc_galaxy_shape_factor_knots_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return log (_nc_galaxy_shape_factor_knots_eval (NC_GALAXY_SHAPE_FACTOR_KNOTS (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

/* GObject boilerplate --------------------------------------------------- */

static void
nc_galaxy_shape_factor_knots_init (NcGalaxyShapeFactorKnots *gsfk)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  self->apply_shear     = NULL;
  self->apply_shear_inv = NULL;
  self->bound           = NC_GALAXY_SHAPE_FACTOR_KNOTS_DEFAULT_BOUND;
  self->n               = NC_GALAXY_SHAPE_FACTOR_KNOTS_DEFAULT_N;
  self->method          = NC_GALAXY_SHAPE_FACTOR_KNOTS_DEFAULT_METHOD;
  self->nodes           = NULL;
  self->weights         = NULL;
}

static void
_nc_galaxy_shape_factor_knots_dispose (GObject *object)
{
  NcGalaxyShapeFactorKnots *gsfk               = NC_GALAXY_SHAPE_FACTOR_KNOTS (object);
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  ncm_vector_clear (&self->nodes);
  ncm_vector_clear (&self->weights);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_knots_parent_class)->dispose (object);
}

static void
_nc_galaxy_shape_factor_knots_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_knots_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorKnots *gsfk               = NC_GALAXY_SHAPE_FACTOR_KNOTS (object);
    NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);
    const NcGalaxyWLObsEllipConv ellip_conv      = nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfk));

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        self->apply_shear     = &nc_wl_ellipticity_apply_shear_trace_c;
        self->apply_shear_inv = &nc_wl_ellipticity_apply_shear_inv_trace_c;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->apply_shear     = &nc_wl_ellipticity_apply_shear_trace_det_c;
        self->apply_shear_inv = &nc_wl_ellipticity_apply_shear_inv_trace_det_c;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }

    _nc_galaxy_shape_factor_knots_regenerate (gsfk);
  }
}

static void
_nc_galaxy_shape_factor_knots_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorKnots *gsfk = NC_GALAXY_SHAPE_FACTOR_KNOTS (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_FACTOR_KNOTS (gsfk));

  switch (prop_id)
  {
    case PROP_BOUND:
      nc_galaxy_shape_factor_knots_set_bound (gsfk, g_value_get_double (value));
      break;
    case PROP_N:
      nc_galaxy_shape_factor_knots_set_n (gsfk, g_value_get_uint (value));
      break;
    case PROP_METHOD:
      nc_galaxy_shape_factor_knots_set_method (gsfk, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_knots_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorKnots *gsfk = NC_GALAXY_SHAPE_FACTOR_KNOTS (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_FACTOR_KNOTS (gsfk));

  switch (prop_id)
  {
    case PROP_BOUND:
      g_value_set_double (value, nc_galaxy_shape_factor_knots_get_bound (gsfk));
      break;
    case PROP_N:
      g_value_set_uint (value, nc_galaxy_shape_factor_knots_get_n (gsfk));
      break;
    case PROP_METHOD:
      g_value_set_enum (value, nc_galaxy_shape_factor_knots_get_method (gsfk));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_knots_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_knots_parent_class)->finalize (object);
}

static void
nc_galaxy_shape_factor_knots_class_init (NcGalaxyShapeFactorKnotsClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_nc_galaxy_shape_factor_knots_constructed;
  object_class->set_property = &_nc_galaxy_shape_factor_knots_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_knots_get_property;
  object_class->dispose      = &_nc_galaxy_shape_factor_knots_dispose;
  object_class->finalize     = &_nc_galaxy_shape_factor_knots_finalize;

  /**
   * NcGalaxyShapeFactorKnots:bound:
   *
   * Half-width $B$ of the box $[u_0-B,u_0+B]\times[v_0-B,v_0+B]$ the fixed
   * knots are laid out over (re-centered on $(u_0,v_0)$ at every evaluation,
   * see the class documentation). Same role and same caveats as
   * #NcGalaxyShapeFactorQuad:bound.
   */
  g_object_class_install_property (object_class,
                                   PROP_BOUND,
                                   g_param_spec_double ("bound",
                                                        NULL,
                                                        "Plane-integration box half-width",
                                                        0.0, G_MAXDOUBLE, NC_GALAXY_SHAPE_FACTOR_KNOTS_DEFAULT_BOUND,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorKnots:n:
   *
   * Number of knots per axis; the 2D rule is the $n\times n$ tensor product.
   */
  g_object_class_install_property (object_class,
                                   PROP_N,
                                   g_param_spec_uint ("n",
                                                      NULL,
                                                      "Number of knots per axis",
                                                      1, G_MAXUINT, NC_GALAXY_SHAPE_FACTOR_KNOTS_DEFAULT_N,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorKnots:method:
   *
   * Fixed-knot generation scheme, see #NcGalaxyShapeFactorKnotsMethod.
   */
  g_object_class_install_property (object_class,
                                   PROP_METHOD,
                                   g_param_spec_enum ("method",
                                                      NULL,
                                                      "Fixed-knot generation scheme",
                                                      NC_TYPE_GALAXY_SHAPE_FACTOR_KNOTS_METHOD,
                                                      NC_GALAXY_SHAPE_FACTOR_KNOTS_DEFAULT_METHOD,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_knots_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_knots_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_knots_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_knots_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_knots_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorKnots.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorKnots.
 */
NcGalaxyShapeFactorKnots *
nc_galaxy_shape_factor_knots_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_KNOTS,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_knots_ref:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 *
 * Increases the reference count of @gsfk by one.
 *
 * Returns: (transfer full): @gsfk.
 */
NcGalaxyShapeFactorKnots *
nc_galaxy_shape_factor_knots_ref (NcGalaxyShapeFactorKnots *gsfk)
{
  return g_object_ref (gsfk);
}

/**
 * nc_galaxy_shape_factor_knots_free:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 *
 * Decreases the reference count of @gsfk by one.
 *
 */
void
nc_galaxy_shape_factor_knots_free (NcGalaxyShapeFactorKnots *gsfk)
{
  g_object_unref (gsfk);
}

/**
 * nc_galaxy_shape_factor_knots_clear:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 *
 * Decreases the reference count of *@gsfk by one, and sets the pointer
 * *@gsfk to NULL.
 *
 */
void
nc_galaxy_shape_factor_knots_clear (NcGalaxyShapeFactorKnots **gsfk)
{
  g_clear_object (gsfk);
}

/**
 * nc_galaxy_shape_factor_knots_set_bound:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 * @bound: the new box half-width $B$
 *
 * Sets the half-width of the $[-B,B]^2$ knot-layout box and regenerates the
 * knots.
 *
 */
void
nc_galaxy_shape_factor_knots_set_bound (NcGalaxyShapeFactorKnots *gsfk, const gdouble bound)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  g_assert_cmpfloat (bound, >, 0.0);

  self->bound = bound;
  _nc_galaxy_shape_factor_knots_regenerate (gsfk);
}

/**
 * nc_galaxy_shape_factor_knots_get_bound:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 *
 * Returns: the box half-width $B$.
 */
gdouble
nc_galaxy_shape_factor_knots_get_bound (NcGalaxyShapeFactorKnots *gsfk)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  return self->bound;
}

/**
 * nc_galaxy_shape_factor_knots_set_n:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 * @n: the new number of knots per axis
 *
 * Sets the number of knots per axis and regenerates the knots.
 *
 */
void
nc_galaxy_shape_factor_knots_set_n (NcGalaxyShapeFactorKnots *gsfk, const guint n)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  g_assert_cmpuint (n, >, 0);

  self->n = n;
  _nc_galaxy_shape_factor_knots_regenerate (gsfk);
}

/**
 * nc_galaxy_shape_factor_knots_get_n:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 *
 * Returns: the number of knots per axis.
 */
guint
nc_galaxy_shape_factor_knots_get_n (NcGalaxyShapeFactorKnots *gsfk)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  return self->n;
}

/**
 * nc_galaxy_shape_factor_knots_set_method:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 * @method: the new #NcGalaxyShapeFactorKnotsMethod
 *
 * Sets the fixed-knot generation scheme and regenerates the knots.
 *
 */
void
nc_galaxy_shape_factor_knots_set_method (NcGalaxyShapeFactorKnots *gsfk, const NcGalaxyShapeFactorKnotsMethod method)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  g_assert_cmpuint (method, <, NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_LEN);

  self->method = method;
  _nc_galaxy_shape_factor_knots_regenerate (gsfk);
}

/**
 * nc_galaxy_shape_factor_knots_get_method:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 *
 * Returns: the #NcGalaxyShapeFactorKnotsMethod.
 */
NcGalaxyShapeFactorKnotsMethod
nc_galaxy_shape_factor_knots_get_method (NcGalaxyShapeFactorKnots *gsfk)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  return self->method;
}

/**
 * nc_galaxy_shape_factor_knots_peek_nodes:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 *
 * Exposes the current 1D nodes on $[-B,B]$ for direct inspection (the 2D
 * rule is their tensor product). Invalidated by any call that changes
 * #NcGalaxyShapeFactorKnots:bound, #NcGalaxyShapeFactorKnots:n or
 * #NcGalaxyShapeFactorKnots:method.
 *
 * Returns: (transfer none): the 1D nodes.
 */
NcmVector *
nc_galaxy_shape_factor_knots_peek_nodes (NcGalaxyShapeFactorKnots *gsfk)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  return self->nodes;
}

/**
 * nc_galaxy_shape_factor_knots_peek_weights:
 * @gsfk: a #NcGalaxyShapeFactorKnots
 *
 * Exposes the current 1D weights matching nc_galaxy_shape_factor_knots_peek_nodes().
 * Invalidated by any call that changes #NcGalaxyShapeFactorKnots:bound,
 * #NcGalaxyShapeFactorKnots:n or #NcGalaxyShapeFactorKnots:method.
 *
 * Returns: (transfer none): the 1D weights.
 */
NcmVector *
nc_galaxy_shape_factor_knots_peek_weights (NcGalaxyShapeFactorKnots *gsfk)
{
  NcGalaxyShapeFactorKnotsPrivate * const self = nc_galaxy_shape_factor_knots_get_instance_private (gsfk);

  return self->weights;
}

