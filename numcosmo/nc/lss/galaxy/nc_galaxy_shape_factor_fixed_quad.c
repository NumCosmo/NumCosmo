/***************************************************************************
 *            nc_galaxy_shape_factor_fixed_quad.c
 *
 *  Thu Jul 9 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_fixed_quad.c
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
 * NcGalaxyShapeFactorFixedQuad:
 *
 * Fixed-node lens-domain quadrature evaluation of the intrinsic-ellipticity
 * marginal.
 *
 * Like #NcGalaxyShapeFactorQuad, evaluates
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_L|<1} \mathrm{d}^2\chi_L\,
 *   P_\mathrm{pop}\big(f_g^{-1}(\chi_L)\big)\,
 *   \left|\det J_{f_g^{-1}}(\chi_L)\right|\,
 *   N_2\big(\epsilon_\mathrm{obs} - \chi_L; \sigma_\mathrm{noise}^2\big)
 * $$
 * exactly (no series truncation in $g$), with a FIXED node count over the
 * INTERSECTION of two discs instead of Quad's adaptive Divonne cubature over
 * a generic box: the noise kernel is supported near $\epsilon_\mathrm{obs}$
 * (radius $\sim n_\sigma\sigma_\mathrm{noise}$), $P_\mathrm{pop}\circ f_g^{-1}$
 * only over the unit disc, and their overlap is a two-circle "lens" in
 * general, a plain disc when one contains the other.
 *
 * The noise kernel lives in $\chi_L$-space and does not depend on $g$, so
 * the whole quadrature domain (node positions, weights, and the
 * noise-kernel value at each node) is $g$-INDEPENDENT: it depends only on
 * $(R,\phi,\sigma_\mathrm{noise})=(\lvert\epsilon_\mathrm{obs}\rvert,
 * \arg\epsilon_\mathrm{obs},\sigma_\mathrm{noise})$, cached per galaxy and
 * reused across every $g$ a fit tries. Validated end to end against an
 * independent scipy oracle and against #NcGalaxyShapeFactorQuad itself (see
 * the test suite and
 * <a href="../../theory/wl_shape_marginalization_fixed_quad.html">Fixed-Node
 * Lens-Domain Quadrature</a>).
 *
 * Every term summed is manifestly non-negative (quadrature weights,
 * population density, $\lvert\det J\rvert$, and the two-arc domain's
 * Jacobian are all non-negative), so unlike #NcGalaxyShapeFactorSeriesLensed
 * there is no truncated polynomial that can cross zero: this class stays
 * accurate at any physical $g$, real or complex, through $\lvert g\rvert=0.99$.
 *
 * Works for ANY population (not just Gaussian): each node evaluates
 * nc_galaxy_shape_pop_eval_p() at $r_i=\lvert\chi_I\rvert$ (this class's own
 * per-node $x_i=\lvert\chi_I\rvert^2$, computed by the shear-map kernels
 * below, is sqrt()'d once before the population call -- see
 * nc_galaxy_shape_pop.h's own eval_p() contract), then converts the
 * returned r-marginal density to the 2D area density this class's own
 * quadrature needs via the explicit $P_\mathrm{pop}(r_i)/(2\pi r_i)$
 * factor, so unlike #NcGalaxyShapeFactorSeriesLensed there is no
 * Gaussian-only guard.
 *
 * Limitation: a fixed grid cannot resolve a population much narrower than
 * its node spacing ($\sigma_\mathrm{pop}\lesssim0.05$, or a sharply
 * concentrated Beta population); use #NcGalaxyShapeFactorQuad for narrower
 * or more exotic populations. Production only uses Gaussian populations
 * with $\sigma_\mathrm{pop}\in(0.2,0.4)$, comfortably inside this class's
 * validated regime. See docs/theory/wl_shape_factor_history.md for why an
 * adaptive alternative was tried and rejected for the narrow-population
 * case, and for the design history of this class more generally.
 *
 * Cost cliff: at $\sigma_\mathrm{noise}$ where the noise disk is comparable
 * in size to the unit disc (roughly $\sigma_\mathrm{noise}\in(0.05,0.2)$,
 * see docs/theory/wl_shape_marginalization_fixed_quad.qmd's own "cost
 * cliff" section), nearly every galaxy lands in the genuine-lens branch
 * ($\mathtt{n\_lens}^2$ nodes, 1681 at the default 41) rather than the
 * cheaper contained branches -- expensive, though this project's actual
 * production regime ($\sigma_\mathrm{noise}\sim0.3$) is unaffected either
 * side of it. #NcGalaxyShapeFactorFixedQuad:auto-lens-nodes (default
 * %FALSE, opt-in) calibrates a per-galaxy lens-branch node count instead
 * of always using the configured #NcGalaxyShapeFactorFixedQuad:n-lens,
 * cutting cost in that regime (~2x fewer nodes typical, more in the
 * expensive middle) with no change to shipped behavior unless explicitly
 * enabled -- see _calibrate_n_lens()'s own docs for the calibration
 * strategy.
 *
 * Repeated calls at many different g for the SAME galaxy (e.g. a z-integral
 * over source-redshift quadrature nodes, or many fit/MCMC iterations)
 * recompute this whole per-node sum from scratch every time.
 * #NcGalaxyShapeFactorFixedQuad:use-marginal-spline (default %FALSE, opt-in)
 * instead caches $\ln P(\epsilon_\mathrm{obs}\mid g)$ as a bivariate
 * function of $g$ over a square
 * $[-g_\mathrm{max},g_\mathrm{max}]^2$ (#NcGalaxyShapeFactorFixedQuad:spline-g-max,
 * a caller-chosen box matching the shear range actually explored, NOT this
 * class' own validated $\lvert g\rvert<0.99$ regime), built lazily once per
 * galaxy and reused across every subsequent g inside that square until the
 * domain cache rebuilds or the population's own parameters change; g
 * outside the square always falls back to the exact direct computation.
 * Two build strategies, chosen automatically per population (see
 * _build_g_spline()'s own docs): this codebase's existing "autoknots" 2D
 * spline machinery (ncm_spline2d_set_function(), the same mechanism
 * nc_halo_mass_function.c uses) for populations whose area density stays
 * bounded at $r=0$; a plain fixed knot grid
 * (_build_g_spline_fixed_knots(), never adaptive, cannot abort) for
 * populations that diverge there -- e.g. alpha<2 Beta populations (the same
 * divergent-density regime #NcGalaxyShapeFactorQuad's own tests already
 * flag near $g\sim0.18$: see this class' own test suite), which the
 * codebase's own users care about most and for which the adaptive path
 * would abort via g_error. The fixed-grid path trades some interpolation
 * accuracy (bounded, not unbounded, worst case -- see its own docs) for a
 * real, if smaller, cache speedup instead of forcing always-direct
 * evaluation.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_fixed_quad.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"
#include "ncm/spline/ncm_spline2d_bicubic.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Window half-width, in units of std_noise, used to size the noise disk
 * (R2 = NSIGMA*std_noise): 8 sigma leaves a Gaussian tail of exp(-32) ~
 * 1e-14. Same convention as NcGalaxyShapeFactorSeriesLensed's own
 * NC_GALAXY_SHAPE_FACTOR_SERIES_LENSED_WINDOW_NSIGMA. */
#define NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_NSIGMA (8.0)

/* The genuine-lens (two-circle partial overlap) branch's (u,v) parametrization
 * degenerates at d=R1+R2: both two-arc half-angles (alpha, beta in
 * _regen_lens) go to zero there, collapsing the grid onto a single point
 * regardless of node count -- a parametrization degeneracy, not a resolution
 * one. _regen_domain() therefore grows the noise disk's EFFECTIVE radius
 * (R2_eff, distinct from the fixed NSIGMA window used for branch selection)
 * so the lens branch always has at least NSIGMA_TAIL sigma of noise-kernel
 * tail depth inside the unit disc, keeping d=R1+R2_eff unreachable for any
 * std_noise>0. See docs/theory/wl_shape_factor_history.md for why this
 * (rather than falling back to full-disc quadrature) is the correct fix. */
#define NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_NSIGMA_TAIL (NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_NSIGMA)

/* Defensive floor against double-precision underflow (every branch always
 * has a nonempty domain -- see _regen_domain()), not a divergence guard:
 * every summed term is already manifestly non-negative (see the class docs
 * above). Same role and value as NcGalaxyShapeFactorSeriesLensed's own
 * MIN_MARGINAL constant. */
#define NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_MIN_MARGINAL (1.0e-300)

/* Grid resolution for _build_g_spline_fixed_knots()'s fallback when the
 * population is unsafe for the adaptive autoknots build (e.g. Beta with
 * alpha<2). Unlike the adaptive path's node count, this is a plain, safe
 * dial: raising it can only improve accuracy (more, closer-together fixed
 * samples of an always-finite function) at a predictable O(N^2) build-cost
 * increase -- there is no discontinuity-detection/abort logic in a
 * fixed-knot build to trip, however sharp the sampled function gets. 33
 * verified (dev session, alpha=1.55/beta=1.62, spline_g_max=0.3): ~4.5ms
 * build, ~12x per-call speedup over always-direct, median interpolation
 * error ~3e-4 in ln(marginal), with a bounded (not exploding) worst-case
 * error around 0.1-0.15 in ln(marginal) at isolated points near wherever a
 * domain node's shear map lands close to chi_I=0 for that g. */
#define NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_UNSAFE_SPLINE_N_KNOTS (33)

/* ===========================================================================
 * GObject boilerplate
 * ===========================================================================
 */

struct _NcGalaxyShapeFactorFixedQuad
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorFixedQuadPrivate
{
  /* Resolved once at construction, same pattern as Quad/VarAdd. Uses
   * the direct (non-log) Jacobian: this class sums linearly, not in
   * log-space, so exp(lndet_jac(...)) would just be a wasted log+exp
   * round-trip -- see nc_wl_ellipticity_det_jac_trace_det()'s docs. */
  complex double (*apply_shear_inv) (complex double g, complex double chi_L);

  gdouble (*det_jac) (complex double g, complex double chi_L);

  /* Resolved once at construction alongside apply_shear_inv/det_jac above.
   * The per-galaxy hot loop in _nc_galaxy_shape_factor_fixed_quad_marginal()
   * branches on this ONCE per call (not per node) to pick between two
   * duplicated loop bodies that call nc_wl_ellipticity.h's direct,
   * non-pointer _trace/_trace_det kernels -- the per-node
   * apply_shear_inv/det_jac indirection above was blocking inlining and
   * vectorization there. apply_shear_inv/det_jac stay as-is for the
   * (non-hot, calibration-only) _lens_quad_at_n(). */
  NcGalaxyWLObsEllipConv ellip_conv;

  guint n_radial;
  guint n_angular;
  guint n_lens; /* forced odd, see constructed() */
  guint n_max;  /* max (n_radial*n_angular, n_lens*n_lens): ldata buffer size */

  /* Opt-in per-galaxy lens-branch node-count calibration -- see
   * _calibrate_n_lens()'s own docs. Off by default: zero behavior change
   * unless explicitly enabled. CONSTRUCT_ONLY (like n-radial/n-angular/
   * n-lens above), not mutable mid-run -- this class's per-galaxy cache
   * (ldata->cache_valid) is invalidated only by epsilon_obs/std_noise
   * changes, so a mid-run change to either of these would silently reuse
   * a stale grid built under the old setting; simplest fix is to not allow
   * that in the first place, matching this class's own existing
   * convention for n-lens itself. */
  gboolean auto_lens_nodes;
  gdouble lens_node_reltol;

  /* Opt-in per-galaxy cache of the marginal itself as a function of g (see
   * the g-spline block below _regen_domain()): CONSTRUCT_ONLY like every
   * other node-count knob on this class, off by default. */
  gboolean use_marginal_spline;
  gdouble spline_g_max;
  gdouble spline_rel_err;
} NcGalaxyShapeFactorFixedQuadPrivate;

enum
{
  PROP_0,
  PROP_N_RADIAL,
  PROP_N_ANGULAR,
  PROP_N_LENS,
  PROP_AUTO_LENS_NODES,
  PROP_LENS_NODE_RELTOL,
  PROP_USE_MARGINAL_SPLINE,
  PROP_SPLINE_G_MAX,
  PROP_SPLINE_REL_ERR,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorFixedQuad, nc_galaxy_shape_factor_fixed_quad, NC_TYPE_GALAXY_SHAPE_FACTOR)

static void
nc_galaxy_shape_factor_fixed_quad_init (NcGalaxyShapeFactorFixedQuad *gsffq)
{
  NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (gsffq);

  self->apply_shear_inv  = NULL;
  self->det_jac          = NULL;
  self->ellip_conv       = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE;
  self->n_radial         = 15;
  self->n_angular        = 15;
  self->n_lens           = 41;
  self->n_max            = 0;
  self->auto_lens_nodes  = FALSE;
  self->lens_node_reltol = 1.0e-4;

  self->use_marginal_spline = FALSE;
  self->spline_g_max        = 0.3;
  self->spline_rel_err      = 1.0e-4;
}

static void
_nc_galaxy_shape_factor_fixed_quad_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorFixedQuad *gsffq              = NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD (object);
  NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (gsffq);

  switch (prop_id)
  {
    case PROP_N_RADIAL:
      self->n_radial = g_value_get_uint (value);
      break;
    case PROP_N_ANGULAR:
      self->n_angular = g_value_get_uint (value);
      break;
    case PROP_N_LENS:
      self->n_lens = g_value_get_uint (value);
      break;
    case PROP_AUTO_LENS_NODES:
      self->auto_lens_nodes = g_value_get_boolean (value);
      break;
    case PROP_LENS_NODE_RELTOL:
      self->lens_node_reltol = g_value_get_double (value);
      break;
    case PROP_USE_MARGINAL_SPLINE:
      self->use_marginal_spline = g_value_get_boolean (value);
      break;
    case PROP_SPLINE_G_MAX:
      self->spline_g_max = g_value_get_double (value);
      break;
    case PROP_SPLINE_REL_ERR:
      self->spline_rel_err = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_fixed_quad_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorFixedQuad *gsffq              = NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD (object);
  NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (gsffq);

  switch (prop_id)
  {
    case PROP_N_RADIAL:
      g_value_set_uint (value, self->n_radial);
      break;
    case PROP_N_ANGULAR:
      g_value_set_uint (value, self->n_angular);
      break;
    case PROP_N_LENS:
      g_value_set_uint (value, self->n_lens);
      break;
    case PROP_AUTO_LENS_NODES:
      g_value_set_boolean (value, self->auto_lens_nodes);
      break;
    case PROP_LENS_NODE_RELTOL:
      g_value_set_double (value, self->lens_node_reltol);
      break;
    case PROP_USE_MARGINAL_SPLINE:
      g_value_set_boolean (value, self->use_marginal_spline);
      break;
    case PROP_SPLINE_G_MAX:
      g_value_set_double (value, self->spline_g_max);
      break;
    case PROP_SPLINE_REL_ERR:
      g_value_set_double (value, self->spline_rel_err);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_fixed_quad_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_fixed_quad_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactor *gsf                         = NC_GALAXY_SHAPE_FACTOR (object);
    NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD (object));
    const NcGalaxyWLObsEllipConv ellip_conv          = nc_galaxy_shape_factor_get_ellip_conv (gsf);

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        self->apply_shear_inv = &nc_wl_ellipticity_apply_shear_inv_trace;
        self->det_jac         = &nc_wl_ellipticity_det_jac_trace;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->apply_shear_inv = &nc_wl_ellipticity_apply_shear_inv_trace_det;
        self->det_jac         = &nc_wl_ellipticity_det_jac_trace_det;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }

    self->ellip_conv = ellip_conv;

    /* Force n_lens odd: guarantees a Gauss-Legendre node lands exactly on
     * the u=0.5 symmetry line, which always passes through the noise-disk
     * center in the genuine-lens branch (see _regen_lens below). */
    self->n_lens |= 1;

    self->n_max = MAX (self->n_radial * self->n_angular, self->n_lens * self->n_lens);
  }
}

static void
_nc_galaxy_shape_factor_fixed_quad_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_fixed_quad_parent_class)->finalize (object);
}

/* ===========================================================================
 * Per-galaxy cache: domain nodes/weights, g-independent.
 * ===========================================================================
 */

typedef struct _NcGalaxyShapeFactorFixedQuadData
{
  gboolean cache_valid;
  gdouble cached_epsilon_obs_1;
  gdouble cached_epsilon_obs_2;
  gdouble cached_std_noise;
  guint n_used;
  complex double *chi_L; /* size n_max */
  gdouble *eff_weight;   /* size n_max, = quadrature_weight * noise_value */
  gdouble *jac;          /* size n_max: |det J_{f_g^-1}| at each node, recomputed every marginal() call */
  GArray *x_arr;         /* size n_max: x_i = |chi_I(chi_L_i,g)|^2 from the kernels, sqrt()'d in place into r_i before feeding eval_p_array() -- see _direct_marginal_at_g()'s own comment */
  GArray *p_arr;         /* size n_max: eval_p_array()'s output (P_pop(r_i)), reused across every g */

  /* Opt-in marginal-as-function-of-g cache (see :use-marginal-spline).
   * Invalidated by _regen_domain() itself (same epoch as chi_L/eff_weight
   * above) plus an independent check against pop's own live pkey -- see
   * _build_g_spline()'s docs. */
  gboolean g_spline_valid;
  guint64 g_spline_pop_pkey;
  NcmSpline2d *g_spline; /* ln(marginal) over [-spline_g_max,spline_g_max]^2, autoknots-built */

  /* Set by _regen_domain(): FALSE whenever any domain node has |chi_L|>1
   * (genuine-lens branch only), which would make the g-spline's adaptive
   * build walk into a real den(g,chi_L)=0 singularity -- see that
   * function's own docs. Overrides :use-marginal-spline for this domain
   * epoch: direct evaluation is always used instead, regardless of the
   * property. */
  gboolean g_spline_safe;

  /* Set by _build_g_spline(): TRUE whenever a spline (of EITHER kind -- see
   * that function's own docs) was actually built for the current pop-pkey
   * epoch. Only FALSE in the (expected-rare) case where the population's
   * own eval_p() isn't even well-defined near r=0 (non-finite or
   * non-positive), which no spline construction can work around -- direct
   * evaluation is used instead. Checked once per pop-pkey epoch,
   * independent of domain. */
  gboolean g_spline_built;
} NcGalaxyShapeFactorFixedQuadData;

static void
_nc_galaxy_shape_factor_fixed_quad_ldata_destroy (gpointer p)
{
  NcGalaxyShapeFactorFixedQuadData *ldata = (NcGalaxyShapeFactorFixedQuadData *) p;

  ncm_spline2d_clear (&ldata->g_spline);

  g_free (ldata->chi_L);
  g_free (ldata->eff_weight);
  g_free (ldata->jac);
  g_array_unref (ldata->x_arr);
  g_array_unref (ldata->p_arr);
  g_free (ldata);
}

static void
_nc_galaxy_shape_factor_fixed_quad_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_factor_fixed_quad_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_fixed_quad_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD (gsf));
  NcGalaxyShapeFactorFixedQuadData *ldata          = g_new0 (NcGalaxyShapeFactorFixedQuadData, 1);

  ldata->cache_valid = FALSE;
  ldata->n_used      = 0;
  ldata->chi_L       = g_new (complex double, self->n_max);
  ldata->eff_weight  = g_new (gdouble, self->n_max);
  ldata->jac         = g_new (gdouble, self->n_max);
  ldata->x_arr       = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), self->n_max);
  ldata->p_arr       = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), self->n_max);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_fixed_quad_ldata_destroy;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_fixed_quad_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_fixed_quad_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_fixed_quad_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_fixed_quad_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  /* No population-capability guard: works for any population, since each
   * node evaluates nc_galaxy_shape_pop_eval_p() directly (see the
   * class docs). */
}

static inline gdouble
_noise_val (complex double delta, gdouble sig2)
{
  const gdouble d2 = gsl_pow_2 (creal (delta)) + gsl_pow_2 (cimag (delta));

  return exp (-d2 / (2.0 * sig2)) / (2.0 * M_PI * sig2);
}

/* Branch 1: noise disk (radius R2) contained in the unit disc, centered at
 * eps_obs -- this project's production regime (std_noise~0.3). Equally-spaced
 * (not Gauss-Legendre) theta nodes, offset by phi=arg(eps_obs) so the grid is
 * a rigid rotation of the same canonical pattern for every galaxy (matching
 * _regen_lens()'s own local-frame convention below): the noise kernel depends
 * on r alone here (radially symmetric about eps_obs), so only the smooth,
 * 2pi-periodic population/Jacobian factor varies with theta, for which
 * equally-spaced sampling is spectrally accurate regardless of the phi
 * offset: n_angular=8-10 already reaches the float floor for that factor,
 * comfortably covered by the default 15 -- but ONLY for smooth (e.g.
 * Gaussian) populations; a population near-singular at chi_I=0 breaks the
 * "smooth" assumption, and without the phi offset here, swapping which side
 * of an exact rotational-covariance identity carries a rotation (e.g. this
 * class' own eval_marginal(g,eps_obs) vs eval_marginal(g*e^{ia},
 * eps_obs*e^{ia})) would silently sample a DIFFERENT, unrelated point of an
 * already-coarse grid instead of the same grid rigidly rotated -- see
 * docs/theory/wl_shape_factor_history.md. */
static void
_regen_noise_contained (NcGalaxyShapeFactorFixedQuadPrivate * const self, complex double eps_obs, gdouble R2, gdouble phi, gdouble sig2,
                        complex double *chi_L, gdouble *eff_weight, guint *n_used)
{
  gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc (self->n_radial);
  const gdouble w_theta                = 2.0 * M_PI / self->n_angular;
  guint i, j, idx = 0;

  for (i = 0; i < self->n_radial; i++)
  {
    gdouble r, wr;

    gsl_integration_glfixed_point (0.0, R2, i, &r, &wr, table);

    for (j = 0; j < self->n_angular; j++)
    {
      const gdouble theta       = phi + 2.0 * M_PI * j / self->n_angular;
      const complex double disp = r * cexp (I * theta);

      chi_L[idx]      = eps_obs + disp;
      eff_weight[idx] = wr * r * w_theta * _noise_val (-disp, sig2);
      idx++;
    }
  }

  gsl_integration_glfixed_table_free (table);
  *n_used = idx;
}

/* Branch 2: full-disc quadrature, centered at the origin -- used both when
 * the unit disc is fully inside the noise disk (large std_noise) AND when
 * the noise disk misses the unit disc entirely (see _regen_domain()'s docs
 * for why both reduce to the same computation). The noise kernel is NOT
 * radially symmetric about the origin (only about eps_obs), so plain
 * Gauss-Legendre in theta -- same phi=arg(eps_obs) offset as Branch 1 above
 * (added to the raw GL node, harmless for a 2pi-periodic integrand: a
 * translated Gauss-Legendre rule over a full period integrates the
 * rotated integrand with the same accuracy, just at rotated node positions)
 * and for the same reason: with NSIGMA=8, this branch is hit by any galaxy
 * with std_noise>~0.125 (R2>R1) -- a real, common part of this project's
 * per-galaxy std_noise distribution (observed range ~0.025-0.47 in a real
 * production catalog, dev session 2026-07-29), not just an edge case -- and
 * without the phi offset its absolute (non-rotating) grid made this class
 * NOT rotation-covariant for narrow/singular populations there: eval_marginal
 * on a rotated config swung by up to 96% of its mean at the class default
 * n_angular=15 before this fix. */
static void
_regen_unit_contained (NcGalaxyShapeFactorFixedQuadPrivate * const self, complex double eps_obs, gdouble phi, gdouble sig2,
                       complex double *chi_L, gdouble *eff_weight, guint *n_used)
{
  gsl_integration_glfixed_table *table_r     = gsl_integration_glfixed_table_alloc (self->n_radial);
  gsl_integration_glfixed_table *table_theta = gsl_integration_glfixed_table_alloc (self->n_angular);
  guint i, j, idx = 0;

  for (i = 0; i < self->n_radial; i++)
  {
    gdouble r, wr;

    gsl_integration_glfixed_point (0.0, 1.0, i, &r, &wr, table_r);

    for (j = 0; j < self->n_angular; j++)
    {
      gdouble theta_raw, theta, wtheta;
      complex double chi;

      gsl_integration_glfixed_point (0.0, 2.0 * M_PI, j, &theta_raw, &wtheta, table_theta);
      theta           = theta_raw + phi;
      chi             = r * cexp (I * theta);
      chi_L[idx]      = chi;
      eff_weight[idx] = wr * r * wtheta * _noise_val (eps_obs - chi, sig2);
      idx++;
    }
  }

  gsl_integration_glfixed_table_free (table_r);
  gsl_integration_glfixed_table_free (table_theta);
  *n_used = idx;
}

/* Branch 3: genuine two-circle partial overlap ("lens"). Two-arc Coons-patch
 * blend in the LOCAL frame (real axis along the line joining the disc
 * centers, i.e. along eps_obs), rotated by phi at the end. See
 * dev-notes/wl_fixed_quad_lens_domain_prototype.py's lens_nodes() for the
 * reference Python implementation.
 *
 * @n_lens: nodes per axis for THIS call -- normally self->n_lens, but
 * _calibrate_n_lens() also calls this at smaller trial values, and the
 * production call site passes a calibrated value when auto-lens-nodes is
 * on (see _regen_domain). Always odd; callers are responsible for that
 * (self->n_lens is forced odd once in constructed(), and
 * _calibrate_n_lens() only ever tries odd values itself). */
static void
_regen_lens (NcGalaxyShapeFactorFixedQuadPrivate * const self, complex double eps_obs, gdouble R2, gdouble d, gdouble phi, gdouble sig2,
             guint n_lens, complex double *chi_L, gdouble *eff_weight, guint *n_used)
{
  const gdouble R1                     = 1.0;
  const gdouble x0                     = (gsl_pow_2 (d) + gsl_pow_2 (R1) - gsl_pow_2 (R2)) / (2.0 * d);
  const gdouble alpha                  = acos (CLAMP (x0 / R1, -1.0, 1.0));
  const gdouble beta                   = acos (CLAMP ((d - x0) / R2, -1.0, 1.0));
  const complex double phase           = cexp (I * phi);
  gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc (n_lens);
  guint i, j, idx = 0;

  for (i = 0; i < n_lens; i++)
  {
    gdouble u, wu;
    gdouble theta1, theta2, x1, y1, x2, y2, dx1_du, dy1_du, dx2_du, dy2_du;

    gsl_integration_glfixed_point (0.0, 1.0, i, &u, &wu, table);

    theta1 = (2.0 * u - 1.0) * alpha;
    theta2 = (2.0 * u - 1.0) * beta;
    x1     = R1 * cos (theta1);
    y1     = R1 * sin (theta1);
    x2     = d - R2 * cos (theta2);
    y2     = R2 * sin (theta2);

    dx1_du = -R1 *sin (theta1) * (2.0 * alpha);

    dy1_du = R1 * cos (theta1) * (2.0 * alpha);
    dx2_du = R2 * sin (theta2) * (2.0 * beta);
    dy2_du = R2 * cos (theta2) * (2.0 * beta);

    for (j = 0; j < n_lens; j++)
    {
      gdouble v, wv;
      gdouble x, y, dx_du, dy_du, dx_dv, dy_dv, jac;
      complex double chi;

      gsl_integration_glfixed_point (0.0, 1.0, j, &v, &wv, table);

      x     = (1.0 - v) * x1 + v * x2;
      y     = (1.0 - v) * y1 + v * y2;
      dx_du = (1.0 - v) * dx1_du + v * dx2_du;
      dy_du = (1.0 - v) * dy1_du + v * dy2_du;
      dx_dv = x2 - x1;
      dy_dv = y2 - y1;
      jac   = fabs (dx_du * dy_dv - dx_dv * dy_du);

      chi             = (x + I * y) * phase; /* rotate local frame into place */
      chi_L[idx]      = chi;
      eff_weight[idx] = wu * wv * jac * _noise_val (eps_obs - chi, sig2);
      idx++;
    }
  }

  gsl_integration_glfixed_table_free (table);
  *n_used = idx;
}

/* Fixed, documented calibration shear for _calibrate_n_lens() -- a
 * representative moderate shear, not tied to any particular galaxy's own
 * g. Validated (dev-notes/wl_fixed_quad_lens_nlens_calibration.py) that
 * calibrating the (g-independent) domain's node count at this one g value
 * generalizes across the full g range a fit actually explores; see
 * docs/theory/wl_shape_factor_history.md. */
#define NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_LENS_CALIB_G (0.15)

/* One-off evaluation of the lens-branch marginal at a trial node count @n,
 * used only by _calibrate_n_lens() below -- never on the per-galaxy
 * production path (which reuses ldata's own persistent buffers via
 * _regen_lens directly). Allocates its own small scratch buffers: this
 * runs at most a handful of times per galaxy (calibration only, not per
 * node/per-g), nowhere near the per-node-per-eval scale the SeriesLensed
 * malloc-churn fix addressed. */
static gdouble
_lens_quad_at_n (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                 complex double eps_obs, gdouble R2, gdouble d, gdouble phi, gdouble sig2, guint n, complex double g)
{
  complex double *chi_L = g_new (complex double, n * n);
  gdouble *eff_weight   = g_new (gdouble, n * n);
  gdouble result        = 0.0;
  guint n_used          = 0;
  guint i;

  _regen_lens (self, eps_obs, R2, d, phi, sig2, n, chi_L, eff_weight, &n_used);

  for (i = 0; i < n_used; i++)
  {
    const complex double chi_i = self->apply_shear_inv (g, chi_L[i]);
    const gdouble x_i          = gsl_pow_2 (creal (chi_i)) + gsl_pow_2 (cimag (chi_i));
    const gdouble r_i          = sqrt (x_i);
    const gdouble p_pop        = nc_galaxy_shape_pop_eval_p (pop, data->pop_data, r_i) / (2.0 * M_PI * r_i);
    const gdouble jac          = self->det_jac (g, chi_L[i]);

    result += eff_weight[i] * p_pop * jac;
  }

  g_free (chi_L);
  g_free (eff_weight);

  return result;
}

/* Calibrates the minimal odd lens-branch node count (per axis, capped at
 * self->n_lens) whose marginal (at the fixed calibration shear above)
 * matches a self-consistent, always-more-accurate reference to
 * self->lens_node_reltol -- same "self-consistent high-resolution
 * reference, no external oracle, geometric bracket then bisection"
 * strategy as ncm_integral_fixed_calibrate() (numcosmo/ncm/integration/
 * ncm_integrate.c), adapted here to this branch's 2D (u,v) grid instead of
 * that function's 1D fixed rule -- the domain shapes genuinely differ, so
 * this is a purpose-built calibration, not a call into the existing one.
 * Reference resolution is 2*self->n_lens+1 (always odd, always stricter
 * than the shipped default), so a calibrated result can never be trusted
 * less than this class's own baseline accuracy. Never returns more than
 * self->n_lens (the search is capped there, matching this property's role
 * as both "the exact count when auto is off" and "the upper bound when
 * auto is on"). */
static guint
_calibrate_n_lens (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                   complex double eps_obs, gdouble R2, gdouble d, gdouble phi, gdouble sig2)
{
  const complex double g_calib = NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_LENS_CALIB_G;
  const guint n_ref            = 2 * self->n_lens + 1;
  const gdouble I_ref          = _lens_quad_at_n (self, pop, data, eps_obs, R2, d, phi, sig2, n_ref, g_calib);
  const gdouble denom          = (I_ref != 0.0) ? fabs (I_ref) : 1.0;
  guint n                      = 5;
  guint last_fail              = 0;
  guint hi                     = 0;

  if (self->n_lens <= 5)
    return self->n_lens;

  while (TRUE)
  {
    const gboolean at_ceiling = (n >= self->n_lens);
    const guint n_try         = (at_ceiling ? self->n_lens : n) | 1;
    const gdouble I_n         = _lens_quad_at_n (self, pop, data, eps_obs, R2, d, phi, sig2, n_try, g_calib);
    const gdouble err         = fabs (I_n - I_ref) / denom;

    if (err < self->lens_node_reltol)
    {
      hi = n_try;
      break;
    }

    last_fail = n_try;

    if (at_ceiling)
      return self->n_lens;  /* no passing config below the cap */

    n = (guint) ceil (n * 1.5);
  }

  {
    guint lo = (last_fail > 0) ? last_fail : 3;

    while (hi - lo > 2)
    {
      guint mid = ((lo + hi) / 2) | 1;

      if (mid <= lo)
        mid = lo + 2;

      {
        const gdouble I_mid   = _lens_quad_at_n (self, pop, data, eps_obs, R2, d, phi, sig2, mid, g_calib);
        const gdouble err_mid = fabs (I_mid - I_ref) / denom;

        if (err_mid < self->lens_node_reltol)
          hi = mid;
        else
          lo = mid;
      }
    }
  }

  return hi;
}

/* R/phi (and the hypot/atan2 needed to get them) are only needed here, on
 * the path that rebuilds the domain -- the cache-validity check in
 * _marginal() below compares the raw epsilon_obs_1/epsilon_obs_2 doubles
 * directly instead, so a fit that holds a galaxy's observed ellipticity
 * fixed across many g values never pays for hypot/atan2. */
static void
_regen_domain (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, NcGalaxyShapeFactorFixedQuadData *ldata,
               gdouble epsilon_obs_1, gdouble epsilon_obs_2, gdouble std_noise)
{
  const gdouble R              = hypot (epsilon_obs_1, epsilon_obs_2);
  const gdouble phi            = atan2 (epsilon_obs_2, epsilon_obs_1);
  const gdouble R1             = 1.0;
  const gdouble R2             = NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_NSIGMA * std_noise;
  const gdouble d              = R;
  const gdouble sig2           = gsl_pow_2 (std_noise);
  const complex double eps_obs = R * cexp (I * phi);

  if ((R2 <= R1) && (d <= R1 - R2))
  {
    /* Noise disk fully inside the unit disc: this project's production
     * regime (std_noise~0.3). */
    _regen_noise_contained (self, eps_obs, R2, phi, sig2, ldata->chi_L, ldata->eff_weight, &ldata->n_used);
  }
  else if (d + R1 <= R2)
  {
    /* Unit disc already fully inside the (fixed-window) noise disk: full
     * disc quadrature's original, validated use case. */
    _regen_unit_contained (self, eps_obs, phi, sig2, ldata->chi_L, ldata->eff_weight, &ldata->n_used);
  }
  else
  {
    /* Genuine partial overlap. R2_eff (not the fixed NSIGMA window used for
     * branch selection above) grows the noise disk's effective radius so
     * that R1+R2_eff > d always holds for std_noise>0, keeping the lens
     * branch's d=R1+R2 parametrization degeneracy unreachable (see
     * NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_NSIGMA_TAIL's docs and
     * docs/theory/wl_shape_factor_history.md). Guard: if this growth would
     * itself make R2_eff fully contain the unit disc (only possible for
     * std_noise comparable to R1/NSIGMA_TAIL), fall through to the
     * full-disc computation instead. */
    const gdouble R2_eff = fmax (R2, (d - R1) + NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_NSIGMA_TAIL * std_noise);

    if (d + R1 <= R2_eff)
    {
      _regen_unit_contained (self, eps_obs, phi, sig2, ldata->chi_L, ldata->eff_weight, &ldata->n_used);
    }
    else
    {
      const guint n_lens = self->auto_lens_nodes ?
                           _calibrate_n_lens (self, pop, data, eps_obs, R2_eff, d, phi, sig2) :
                           self->n_lens;

      _regen_lens (self, eps_obs, R2_eff, d, phi, sig2, n_lens, ldata->chi_L, ldata->eff_weight, &ldata->n_used);
    }
  }

  ldata->cached_epsilon_obs_1 = epsilon_obs_1;
  ldata->cached_epsilon_obs_2 = epsilon_obs_2;
  ldata->cached_std_noise     = std_noise;
  ldata->cache_valid          = TRUE;

  /* The domain just moved (new chi_L/eff_weight), so any previously built
   * g-spline (see :use-marginal-spline) no longer applies -- this is the
   * single hook that keeps it correct without duplicating this function's
   * own epsilon_obs/std_noise change detection. */
  ldata->g_spline_valid = FALSE;

  /* Safety gate for :use-marginal-spline: any domain node with |chi_L|>1
   * (only possible in the genuine-lens branch above -- noise-contained and
   * unit-contained both place nodes with |chi_L|<=1 by construction) makes
   * the per-node shear-inversion Jacobian's denominator den(g,chi_L) vanish
   * on a real locus in g-space: solving den=0 for g at fixed chi_L gives
   * the circle |g-chi_L|^2=|chi_L|^2-1 (TRACE) or the single point
   * g=1/conj(chi_L) (TRACE_DET) -- both require |chi_L|>1 to have any real
   * solution at all. If that locus falls near the cached box, the adaptive
   * autoknots build in _build_g_spline() will correctly detect it cannot
   * resolve a genuine unbounded singularity there and abort (verified
   * against a real production catalog, reproduced identically with
   * OMP_NUM_THREADS=1, ruling out a threading race) -- the direct path
   * never hits this because it evaluates once per g and returns a
   * big-but-finite number instead of adaptively bisecting toward the
   * singularity. Scanning nodes here (not analytically checking whether the
   * locus actually reaches the box) is deliberately conservative: cheap
   * relative to the domain build itself, branch-agnostic, and safe by
   * construction rather than by geometry that would need re-deriving per
   * ellip-conv. */
  {
    gboolean node_beyond_unit_disc = FALSE;
    guint i;

    for (i = 0; i < ldata->n_used; i++)
    {
      if (creal (ldata->chi_L[i] * conj (ldata->chi_L[i])) > 1.0)
      {
        node_beyond_unit_disc = TRUE;
        break;
      }
    }

    ldata->g_spline_safe = !node_beyond_unit_disc;
  }
}

/* Raw (un-floored) marginal at a single g, given an already-valid domain
 * cache (ldata->chi_L/eff_weight, see _regen_domain()). Extracted out of
 * _nc_galaxy_shape_factor_fixed_quad_marginal() so both that function's
 * direct path AND _build_g_spline()'s per-sample-point evaluations below
 * share one implementation -- behavior is unchanged from before this
 * extraction. Mutates ldata's own scratch buffers (x_arr/jac/p_arr), so
 * repeated calls (as _build_g_spline() makes many of, one per (g_1,g_2)
 * grid/slice sample) are safe but not reentrant -- same constraint the
 * un-extracted code always had. */
static gdouble
_direct_marginal_at_g (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                       NcGalaxyShapeFactorFixedQuadData *ldata, const complex double g)
{
  gdouble result = 0.0;
  guint i;

  /* Pass 1: every x_i is known before eval_p() is ever called (the whole
   * point of a FIXED node count), so batch them through eval_p_array()
   * instead of n_used one-at-a-time vfunc calls -- see nc_galaxy_shape_pop.h's
   * eval_p_array doc comment. jac is cheap (unlike p_pop) and stays a plain
   * per-node computation, stashed here rather than recomputed in pass 2.
   * ldata->x_arr holds x_i=|chi_I|^2 from the kernels below (the Möbius
   * algebra naturally produces x, not r), then gets sqrt()'d in place into
   * r_i before the population call, since nc_galaxy_shape_pop_eval_p_array()
   * takes r now (see the eval_p/eval_p_rho2 contract collapse) -- nothing
   * else reads x_arr's contents in between, so the in-place conversion is
   * safe. */
  g_array_set_size (ldata->x_arr, ldata->n_used);

  {
    gdouble * const x_data = (gdouble *) ldata->x_arr->data;

    /* Branch on ellip_conv ONCE per call (it is fixed for the object's
     * lifetime), not per node: each duplicated loop below calls
     * nc_wl_ellipticity.h's direct, non-pointer fused *_kernel() (x_i and
     * jac from a single g_conj/abs_g2/den, no separate apply_shear_inv +
     * det_jac calls -- see its own docs), letting the compiler inline/
     * vectorize freely -- the per-node self->apply_shear_inv/self->det_jac
     * indirection was blocking that. TRACE additionally prepares its
     * g-only terms once before the loop (measured ~12% faster, see
     * nc_wl_ellipticity_trace_kernel_prepare()'s docs); TRACE_DET has no
     * equivalent win and stays a single call per node (see
     * nc_wl_ellipticity_trace_det_kernel()'s docs). */
    if (self->ellip_conv == NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE)
    {
      NcWLEllipticityTraceKernelPrep prep;

      nc_wl_ellipticity_trace_kernel_prepare (g, &prep);

      for (i = 0; i < ldata->n_used; i++)
        nc_wl_ellipticity_trace_kernel_apply (&prep, ldata->chi_L[i], &x_data[i], &ldata->jac[i]);
    }
    else
    {
      for (i = 0; i < ldata->n_used; i++)
        nc_wl_ellipticity_trace_det_kernel (g, ldata->chi_L[i], &x_data[i], &ldata->jac[i]);
    }

    for (i = 0; i < ldata->n_used; i++)
      x_data[i] = sqrt (x_data[i]);
  }

  nc_galaxy_shape_pop_eval_p_array (pop, data->pop_data, ldata->x_arr, &ldata->p_arr);

  {
    const gdouble * const p_data = (const gdouble *) ldata->p_arr->data;
    const gdouble * const r_data = (const gdouble *) ldata->x_arr->data;

    /* p_data[i]=P_pop(r_i); the 2D area density these terms need is
     * P_pop(r_i)/(2*pi*r_i) -- the same physical division that used to be
     * folded into the old eval_p(x)'s own exponent, now explicit here
     * (see the class docs' own note on this being an accepted, deliberate
     * minor precision cost in exchange for the interface simplification). */
    for (i = 0; i < ldata->n_used; i++)
      result += ldata->eff_weight[i] * (p_data[i] / (2.0 * M_PI * r_data[i])) * ldata->jac[i];
  }

  return result;
}

/* Shared mutable args behind _build_g_spline()'s Fx/Fy gsl_function slices:
 * the SAME struct backs both (one field toggling which axis is the free
 * variable and which is held fixed), following ncm_spline2d_set_function()'s
 * own contract (see its docs) and the exact idiom
 * _nc_halo_mass_function_generate_2Dspline_knots() (numcosmo/nc/lss/halo/
 * nc_halo_mass_function.c) already uses for a different 2D function. */
typedef struct
{
  NcGalaxyShapeFactorFixedQuadPrivate *self;
  NcGalaxyShapePop *pop;
  NcGalaxyShapeFactorData *data;
  NcGalaxyShapeFactorFixedQuadData *ldata;
  gboolean slice_is_g1; /* TRUE: F(g_1) at fixed g_2=0; FALSE: F(g_2) at fixed g_1=0 */
} GSplineSliceArgs;

static gdouble
_g_spline_slice_func (gdouble x, gpointer p)
{
  GSplineSliceArgs *a    = (GSplineSliceArgs *) p;
  const complex double g = a->slice_is_g1 ? x : x * I;
  const gdouble v        = _direct_marginal_at_g (a->self, a->pop, a->data, a->ldata, g);

  return log (fmax (v, NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_MIN_MARGINAL));
}

/* Builds (or rebuilds) ldata's g-spline: ln(marginal), cached as a genuine
 * bivariate function of (g_1,g_2) over the square
 * [-spline_g_max,spline_g_max]^2 -- NOT a reduction to fewer variables (an
 * earlier design tried reducing to (|g|,arg(g)-arg(epsilon_obs)) via a
 * per-radius real-DFT-in-angle decomposition, exploiting this class' own
 * exact rotation-covariance -- see nc_wl_ellipticity_apply_shear_inv_trace()/
 * _trace_det()'s docs for that covariance property, still true and still the
 * reason this cache needs no knowledge of the per-galaxy calibration terms
 * the caller folds into g -- but a small, fixed harmonic
 * count could not track this function's actual angular sharpness: verified
 * numerically to blow up by orders of magnitude away from the build nodes,
 * worst for exactly the alpha<2 Beta populations this class' users care
 * about most; see the branch history for the failed prototype). Working in
 * ln-space here is what actually matters: the marginal itself spans many
 * orders of magnitude over this square (verified: ratios of 1e4-1e10
 * between nearby points are common), which defeats any polynomial/spline
 * fit in linear space regardless of node placement.
 *
 * Node placement uses this codebase's existing "autoknots" machinery
 * instead of a fixed grid: ncm_spline2d_set_function() (numcosmo/ncm/
 * spline/ncm_spline2d.c), the same mechanism
 * _nc_halo_mass_function_generate_2Dspline_knots() uses for an unrelated 2D
 * function. It adaptively places knots along each axis from two 1D slices
 * (_g_spline_slice_func() above, sliced at the OTHER coordinate fixed to 0 --
 * a natural choice here since it's always inside the square and contains no
 * special structure that would make it a bad representative) to
 * spline_rel_err, then this function fills in the TRUE 2D ln(marginal) at
 * the resulting tensor grid (peek_xv/yv/zm) before calling
 * ncm_spline2d_prepare() -- so accuracy away from those two representative
 * slices depends on how well 1D-slice-calibrated knots generalize across
 * the other axis, which is NOT guaranteed in general (verified: usually
 * very good, but a population with a genuinely sharp -- not divergent, just
 * sharp -- feature away from both slices could still miss the target
 * rel_err by orders of magnitude off them; this path is only ever reached
 * by populations that already passed the r=0 area-density check below, so
 * a truly divergent population never exercises it at all -- see
 * _build_g_spline_fixed_knots() for that case instead). Node count is
 * therefore data-dependent, not a user-chosen constant: anywhere from tens
 * to (rarely) tens of thousands of builds depending on how sharply the
 * population's own density varies. */

/* Fallback for populations unsafe for the adaptive autoknots build above
 * (e.g. Beta with alpha<2): a plain, FIXED, uniform knot grid over
 * [-spline_g_max,spline_g_max]^2, filled by the same _direct_marginal_at_g()
 * ground truth and set via ncm_spline2d_set()+prepare() instead of
 * ncm_spline2d_set_function() -- i.e. no adaptive refinement, hence no
 * discontinuity-detection/abort logic to trip no matter how sharp
 * ln(marginal) actually gets near a domain node's own chi_I=0 crossing (see
 * _build_g_spline()'s own docs for why that crossing is unavoidable for
 * these populations). This trades the adaptive path's near-exact accuracy
 * for a strictly bounded, non-crashing worst case: verified (dev session)
 * a median ln(marginal) interpolation error of order 1e-4 with an isolated
 * worst case around 0.1-0.15 at the one or two grid cells straddling a
 * chi_I=0 crossing, in exchange for a real (~12x per-call, dev-session
 * measurement) cache speedup that always-direct evaluation cannot offer at
 * all. See NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_UNSAFE_SPLINE_N_KNOTS's own
 * docs for why raising the knot count is always safe here. */
static void
_build_g_spline_fixed_knots (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                             NcGalaxyShapeFactorFixedQuadData *ldata)
{
  const gdouble g_max = self->spline_g_max;
  const guint n_knots = NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_UNSAFE_SPLINE_N_KNOTS;
  NcmVector *xv       = ncm_vector_new (n_knots);
  NcmVector *yv       = ncm_vector_new (n_knots);
  NcmMatrix *zm       = ncm_matrix_new (n_knots, n_knots);
  guint i, j;

  for (i = 0; i < n_knots; i++)
  {
    const gdouble v = -g_max + (2.0 * g_max) * i / (n_knots - 1.0);

    ncm_vector_set (xv, i, v);
    ncm_vector_set (yv, i, v);
  }

  for (i = 0; i < n_knots; i++)
  {
    const gdouble g_2 = ncm_vector_get (yv, i);

    for (j = 0; j < n_knots; j++)
    {
      const gdouble g_1      = ncm_vector_get (xv, j);
      const complex double g = g_1 + I * g_2;
      const gdouble v        = _direct_marginal_at_g (self, pop, data, ldata, g);

      ncm_matrix_set (zm, i, j, log (fmax (v, NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_MIN_MARGINAL)));
    }
  }

  ncm_spline2d_clear (&ldata->g_spline);
  ldata->g_spline = NCM_SPLINE2D (ncm_spline2d_bicubic_notaknot_new ());
  ncm_spline2d_set (ldata->g_spline, xv, yv, zm, TRUE);
  ncm_spline2d_prepare (ldata->g_spline);

  ncm_vector_free (xv);
  ncm_vector_free (yv);
  ncm_matrix_free (zm);
}

static void
_build_g_spline (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                 NcGalaxyShapeFactorFixedQuadData *ldata)
{
  const gdouble g_max     = self->spline_g_max;
  GSplineSliceArgs args_x = { self, pop, data, ldata, TRUE };
  GSplineSliceArgs args_y = { self, pop, data, ldata, FALSE };
  gsl_function Fx         = { &_g_spline_slice_func, &args_x };
  gsl_function Fy         = { &_g_spline_slice_func, &args_y };
  NcmVector *xv, *yv;
  NcmMatrix *zm;
  guint i, j;

  ldata->g_spline_pop_pkey = ncm_model_state_get_pkey (NCM_MODEL (pop));
  ldata->g_spline_valid    = TRUE;

  /* Population-divergence check: the AREA DENSITY this class actually
   * needs, P_pop(r)/(2*pi*r) (see the class doc and _direct_marginal_at_g()
   * below), diverges at r=0 whenever eval_p(r) vanishes slower than r
   * itself as r->0 -- e.g. NcGalaxyShapePopBeta with alpha<2, whose
   * eval_p(r) ~ r^(alpha-1) (see _nc_galaxy_shape_pop_beta_eval_p()).
   * Testing eval_p AT r=0 directly (as an earlier version of this guard
   * did, before the eval_p(x)/eval_p_rho2(rho2) contract collapse made
   * eval_p() itself r-native) can no longer distinguish this: eval_p(0.0)
   * collapses to exactly 0.0 in floating point for ANY alpha>1, no matter
   * how close to 1 -- the vanishing ORDER is lost, not just its value.
   * Instead, probe the local power-law exponent via two small,
   * well-separated r: eval_p(r) ~ C*r^p near 0 for every population
   * implemented here (Gauss: p=1 exactly, the safe boundary; Beta: p =
   * alpha-1), so the area density diverges iff the finite log-log slope
   * between the two samples is below 1 -- scale-invariant, no new vfunc
   * needed, and reduces to the same alpha>=2 safety boundary as the old
   * x-space test. This matters here specifically -- not just as the
   * already-documented direct-quadrature accuracy caveat -- because EVERY
   * domain node chi_L has some g where the implied intrinsic ellipticity
   * chi_I(g,chi_L) is exactly 0 (solving apply_shear_inv(g,chi_L)=0 always
   * has a root), so a population divergent at r=0 turns each of the
   * (hundreds to thousands of) domain nodes into its own genuine
   * unbounded spike somewhere in g-space. A box of any reasonable size is
   * essentially guaranteed to contain at least one, and the adaptive
   * autoknots build below cannot resolve a true singularity and would
   * abort via g_error (verified against a real production catalog crash --
   * see git history) -- so this population instead gets
   * _build_g_spline_fixed_knots()'s fixed grid, which cannot abort no
   * matter how sharp the sampled function gets (see its own docs for the
   * accuracy/speed tradeoff that buys). Checked once per pop-pkey epoch
   * (population-only, not domain-dependent), not per galaxy, so this never
   * touches the per-node chi_L scan below unless the population is
   * actually safe for the ADAPTIVE path specifically. */
  {
    const gdouble r1             = 1.0e-6;
    const gdouble r2             = 1.0e-3;
    const gdouble p1             = nc_galaxy_shape_pop_eval_p (pop, data->pop_data, r1);
    const gdouble p2             = nc_galaxy_shape_pop_eval_p (pop, data->pop_data, r2);
    const gboolean well_defined  = gsl_finite (p1) && gsl_finite (p2) && (p1 > 0.0) && (p2 > 0.0);
    const gboolean adaptive_safe = well_defined && (log (p2 / p1) / log (r2 / r1) >= 1.0 - 1.0e-4);

    if (!well_defined)
    {
      /* Can't even probe the local behavior near r=0 (e.g. eval_p itself
       * misbehaves there) -- no spline construction can work around that;
       * disable caching entirely for this pkey epoch. */
      ldata->g_spline_built = FALSE;

      return;
    }

    ldata->g_spline_built = TRUE;

    if (!adaptive_safe)
    {
      _build_g_spline_fixed_knots (self, pop, data, ldata);

      return;
    }
  }

  ncm_spline2d_clear (&ldata->g_spline);
  ldata->g_spline = NCM_SPLINE2D (ncm_spline2d_bicubic_notaknot_new ());

  ncm_spline2d_set_function (ldata->g_spline, NCM_SPLINE_FUNCTION_SPLINE, &Fx, &Fy,
                             -g_max, g_max, -g_max, g_max, self->spline_rel_err);

  xv = ncm_spline2d_peek_xv (ldata->g_spline);
  yv = ncm_spline2d_peek_yv (ldata->g_spline);
  zm = ncm_spline2d_peek_zm (ldata->g_spline);

  for (i = 0; i < ncm_vector_len (yv); i++)
  {
    const gdouble g_2 = ncm_vector_get (yv, i);

    for (j = 0; j < ncm_vector_len (xv); j++)
    {
      const gdouble g_1      = ncm_vector_get (xv, j);
      const complex double g = g_1 + I * g_2;
      const gdouble v        = _direct_marginal_at_g (self, pop, data, ldata, g);

      ncm_matrix_set (zm, i, j, log (fmax (v, NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_MIN_MARGINAL)));
    }
  }

  ncm_spline2d_prepare (ldata->g_spline);
}

/* Dispatch for :use-marginal-spline: falls back to the exact direct
 * computation outside the cached square (see :spline-g-max's docs),
 * rebuilds the spline whenever it's not valid or pop's own live pkey moved
 * on since it was built (_build_g_spline()'s own docs explain why nothing
 * else needs checking here), and otherwise evaluates the cached
 * ln(marginal) surface directly. */
static gdouble
_eval_marginal_spline (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                       NcGalaxyShapeFactorFixedQuadData *ldata, const complex double g)
{
  const gdouble g_1 = creal (g);
  const gdouble g_2 = cimag (g);

  if ((fabs (g_1) > self->spline_g_max) || (fabs (g_2) > self->spline_g_max))
    return _direct_marginal_at_g (self, pop, data, ldata, g);

  {
    const guint64 pop_pkey = ncm_model_state_get_pkey (NCM_MODEL (pop));

    if (!ldata->g_spline_valid || (ldata->g_spline_pop_pkey != pop_pkey))
      _build_g_spline (self, pop, data, ldata);
  }

  /* eval_p() isn't even well-defined near r=0 for this pkey epoch -- see
   * _build_g_spline()'s own docs: it returned early without building
   * anything (this is separate from, and rarer than, the adaptive-vs-fixed-
   * knots choice _build_g_spline() otherwise makes silently). */
  if (!ldata->g_spline_built)
    return _direct_marginal_at_g (self, pop, data, ldata, g);

  return exp (ncm_spline2d_eval (ldata->g_spline, g_1, g_2));
}

static gdouble
_nc_galaxy_shape_factor_fixed_quad_marginal (NcGalaxyShapeFactorFixedQuad *gsffq, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                             const gdouble g_1, const gdouble g_2,
                                             const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (gsffq);
  NcGalaxyShapeFactorFixedQuadData *ldata          = (NcGalaxyShapeFactorFixedQuadData *) data->ldata;
  const complex double g                           = g_1 + I * g_2;
  gdouble result;

  /* Compares the raw epsilon_obs_1/epsilon_obs_2 directly -- R/phi (and the
   * hypot/atan2 to compute them) are only needed inside _regen_domain, on
   * the rare path that actually rebuilds the domain. See that function's
   * docs. */
  if (!ldata->cache_valid || (ldata->cached_epsilon_obs_1 != epsilon_obs_1) ||
      (ldata->cached_epsilon_obs_2 != epsilon_obs_2) || (ldata->cached_std_noise != data->std_noise))
    _regen_domain (self, pop, data, ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise);

  result = (self->use_marginal_spline && ldata->g_spline_safe) ?
           _eval_marginal_spline (self, pop, data, ldata, g) :
           _direct_marginal_at_g (self, pop, data, ldata, g);

  /* Purely the empty-domain / underflow floor -- see the class docs. Also
   * catches the rare small-negative dip a truncated-harmonic reconstruction
   * of a strictly-positive function can produce near sharp features. */
  if (!(result > 0.0))
    return NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_MIN_MARGINAL;

  return result;
}

static gdouble
_nc_galaxy_shape_factor_fixed_quad_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return _nc_galaxy_shape_factor_fixed_quad_marginal (NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

static gdouble
_nc_galaxy_shape_factor_fixed_quad_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return log (_nc_galaxy_shape_factor_fixed_quad_marginal (NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

static void
nc_galaxy_shape_factor_fixed_quad_class_init (NcGalaxyShapeFactorFixedQuadClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_shape_factor_fixed_quad_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_fixed_quad_get_property;
  object_class->constructed  = &_nc_galaxy_shape_factor_fixed_quad_constructed;
  object_class->finalize     = &_nc_galaxy_shape_factor_fixed_quad_finalize;

  /**
   * NcGalaxyShapeFactorFixedQuad:n-radial:
   *
   * Number of fixed Gauss-Legendre nodes in the radial direction (branches
   * where one disc is contained in the other). Default 15.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_RADIAL,
                                   g_param_spec_uint ("n-radial",
                                                      "Number of radial nodes",
                                                      "Number of fixed Gauss-Legendre nodes in the radial direction",
                                                      1, G_MAXUINT, 15,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorFixedQuad:n-angular:
   *
   * Number of angular nodes (equally-spaced when the noise disk is
   * contained in the unit disc, Gauss-Legendre otherwise). Default 15.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_ANGULAR,
                                   g_param_spec_uint ("n-angular",
                                                      "Number of angular nodes",
                                                      "Number of angular quadrature nodes",
                                                      1, G_MAXUINT, 15,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorFixedQuad:n-lens:
   *
   * Number of fixed Gauss-Legendre nodes per axis in the genuine two-circle
   * "lens" (partial-overlap) branch; always rounded up to the next odd
   * number (see the class docs for why). Default 41.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_LENS,
                                   g_param_spec_uint ("n-lens",
                                                      "Number of lens-branch nodes",
                                                      "Number of fixed Gauss-Legendre nodes per axis in the genuine-lens branch",
                                                      3, G_MAXUINT, 41,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorFixedQuad:auto-lens-nodes:
   *
   * When %TRUE, calibrates a per-galaxy lens-branch node count (capped at
   * #NcGalaxyShapeFactorFixedQuad:n-lens) instead of always using
   * #NcGalaxyShapeFactorFixedQuad:n-lens for every galaxy -- see
   * _calibrate_n_lens()'s docs. Default %FALSE (zero behavior change unless
   * explicitly enabled).
   */
  g_object_class_install_property (object_class,
                                   PROP_AUTO_LENS_NODES,
                                   g_param_spec_boolean ("auto-lens-nodes",
                                                         "Auto lens-branch nodes",
                                                         "Calibrate a per-galaxy lens-branch node count instead of always using n-lens",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorFixedQuad:lens-node-reltol:
   *
   * Target relative tolerance for #NcGalaxyShapeFactorFixedQuad:auto-lens-nodes's
   * calibration. Default 1e-4.
   */
  g_object_class_install_property (object_class,
                                   PROP_LENS_NODE_RELTOL,
                                   g_param_spec_double ("lens-node-reltol",
                                                        "Lens-branch node calibration reltol",
                                                        "Target relative tolerance for auto-lens-nodes' calibration",
                                                        0.0, 1.0, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorFixedQuad:use-marginal-spline:
   *
   * When %TRUE, caches the marginal itself as a function of g, per galaxy:
   * built once (lazily, on first use), over the square
   * $[-g_\mathrm{max},g_\mathrm{max}]^2$ set by
   * #NcGalaxyShapeFactorFixedQuad:spline-g-max, and reused across every
   * subsequent g request for that galaxy inside that square until the
   * domain cache rebuilds (epsilon_obs/std_noise change) or the
   * population's parameters change -- see _build_g_spline()'s own docs. A
   * request for g outside that square always falls back to the exact
   * direct computation (never extrapolated), so enabling this can only
   * trade cache-region accuracy (bounded by
   * #NcGalaxyShapeFactorFixedQuad:spline-rel-err) for speed, never
   * correctness outside it. Default %FALSE (zero behavior change unless
   * explicitly enabled).
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_MARGINAL_SPLINE,
                                   g_param_spec_boolean ("use-marginal-spline",
                                                         "Use marginal spline",
                                                         "Cache the marginal as a function of g instead of recomputing it every call",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorFixedQuad:spline-g-max:
   *
   * Half-side of the square $[-g_\mathrm{max},g_\mathrm{max}]^2$
   * #NcGalaxyShapeFactorFixedQuad:use-marginal-spline's cache covers --
   * NOT this class' own validated $\lvert g\rvert<0.99$ regime: a much
   * smaller box matching the shear range a given analysis actually
   * explores keeps the autoknots build (see _build_g_spline()'s docs)
   * cheap; g outside the box still gets the exact direct computation.
   * Default 0.3.
   */
  g_object_class_install_property (object_class,
                                   PROP_SPLINE_G_MAX,
                                   g_param_spec_double ("spline-g-max",
                                                        "g-spline cached box half-side",
                                                        "Half-side of the square use-marginal-spline's cache covers",
                                                        0.0, 1.0, 0.3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorFixedQuad:spline-rel-err:
   *
   * Target relative error for #NcGalaxyShapeFactorFixedQuad:use-marginal-spline's
   * autoknots build (passed directly to ncm_spline2d_set_function()).
   * Default 1e-4.
   */
  g_object_class_install_property (object_class,
                                   PROP_SPLINE_REL_ERR,
                                   g_param_spec_double ("spline-rel-err",
                                                        "g-spline target relative error",
                                                        "Target relative error for use-marginal-spline's autoknots build",
                                                        0.0, 1.0, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_fixed_quad_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_fixed_quad_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_fixed_quad_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_fixed_quad_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_fixed_quad_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorFixedQuad.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorFixedQuad.
 */
NcGalaxyShapeFactorFixedQuad *
nc_galaxy_shape_factor_fixed_quad_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_FIXED_QUAD,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_fixed_quad_ref:
 * @gsffq: a #NcGalaxyShapeFactorFixedQuad
 *
 * Increases the reference count of @gsffq by one.
 *
 * Returns: (transfer full): @gsffq.
 */
NcGalaxyShapeFactorFixedQuad *
nc_galaxy_shape_factor_fixed_quad_ref (NcGalaxyShapeFactorFixedQuad *gsffq)
{
  return g_object_ref (gsffq);
}

/**
 * nc_galaxy_shape_factor_fixed_quad_free:
 * @gsffq: a #NcGalaxyShapeFactorFixedQuad
 *
 * Decreases the reference count of @gsffq by one.
 *
 */
void
nc_galaxy_shape_factor_fixed_quad_free (NcGalaxyShapeFactorFixedQuad *gsffq)
{
  g_object_unref (gsffq);
}

/**
 * nc_galaxy_shape_factor_fixed_quad_clear:
 * @gsffq: a #NcGalaxyShapeFactorFixedQuad
 *
 * Decreases the reference count of *@gsffq by one, and sets the pointer
 * *@gsffq to NULL.
 *
 */
void
nc_galaxy_shape_factor_fixed_quad_clear (NcGalaxyShapeFactorFixedQuad **gsffq)
{
  g_clear_object (gsffq);
}

