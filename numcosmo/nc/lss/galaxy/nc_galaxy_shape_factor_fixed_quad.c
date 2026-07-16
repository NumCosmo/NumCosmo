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
 * nc_galaxy_shape_pop_eval_p() directly at $x_i=\lvert\chi_I\rvert^2$, the
 * same convention #NcGalaxyShapeFactorQuad uses (not eval_p_rho2(), whose
 * $\rho^2$ argument is the pre-compactification radius of the $(u,v)$-plane
 * substitution this class does not use), so unlike
 * #NcGalaxyShapeFactorSeriesLensed there is no Gaussian-only guard.
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
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_fixed_quad.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"

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
} NcGalaxyShapeFactorFixedQuadPrivate;

enum
{
  PROP_0,
  PROP_N_RADIAL,
  PROP_N_ANGULAR,
  PROP_N_LENS,
  PROP_AUTO_LENS_NODES,
  PROP_LENS_NODE_RELTOL,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorFixedQuad, nc_galaxy_shape_factor_fixed_quad, NC_TYPE_GALAXY_SHAPE_FACTOR)

static void
nc_galaxy_shape_factor_fixed_quad_init (NcGalaxyShapeFactorFixedQuad *gsffq)
{
  NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (gsffq);

  self->apply_shear_inv  = NULL;
  self->det_jac          = NULL;
  self->n_radial         = 15;
  self->n_angular        = 15;
  self->n_lens           = 41;
  self->n_max            = 0;
  self->auto_lens_nodes  = FALSE;
  self->lens_node_reltol = 1.0e-4;
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
} NcGalaxyShapeFactorFixedQuadData;

static void
_nc_galaxy_shape_factor_fixed_quad_ldata_destroy (gpointer p)
{
  NcGalaxyShapeFactorFixedQuadData *ldata = (NcGalaxyShapeFactorFixedQuadData *) p;

  g_free (ldata->chi_L);
  g_free (ldata->eff_weight);
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
 * (not Gauss-Legendre) theta nodes: the noise kernel depends on r alone here
 * (radially symmetric about eps_obs), so only the smooth, 2pi-periodic
 * population/Jacobian factor varies with theta, for which equally-spaced
 * sampling is spectrally accurate: n_angular=8-10 already reaches the float
 * floor, comfortably covered by the default 15. */
static void
_regen_noise_contained (NcGalaxyShapeFactorFixedQuadPrivate * const self, complex double eps_obs, gdouble R2, gdouble sig2,
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
      const gdouble theta       = 2.0 * M_PI * j / self->n_angular;
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
 * Gauss-Legendre in theta. */
static void
_regen_unit_contained (NcGalaxyShapeFactorFixedQuadPrivate * const self, complex double eps_obs, gdouble sig2,
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
      gdouble theta, wtheta;
      complex double chi;

      gsl_integration_glfixed_point (0.0, 2.0 * M_PI, j, &theta, &wtheta, table_theta);
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
    const gdouble p_pop        = nc_galaxy_shape_pop_eval_p (pop, data->pop_data, x_i) / M_PI;
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
    _regen_noise_contained (self, eps_obs, R2, sig2, ldata->chi_L, ldata->eff_weight, &ldata->n_used);
  }
  else if (d + R1 <= R2)
  {
    /* Unit disc already fully inside the (fixed-window) noise disk: full
     * disc quadrature's original, validated use case. */
    _regen_unit_contained (self, eps_obs, sig2, ldata->chi_L, ldata->eff_weight, &ldata->n_used);
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
      _regen_unit_contained (self, eps_obs, sig2, ldata->chi_L, ldata->eff_weight, &ldata->n_used);
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
  guint i;

  /* Compares the raw epsilon_obs_1/epsilon_obs_2 directly -- R/phi (and the
   * hypot/atan2 to compute them) are only needed inside _regen_domain, on
   * the rare path that actually rebuilds the domain. See that function's
   * docs. */
  if (!ldata->cache_valid || (ldata->cached_epsilon_obs_1 != epsilon_obs_1) ||
      (ldata->cached_epsilon_obs_2 != epsilon_obs_2) || (ldata->cached_std_noise != data->std_noise))
    _regen_domain (self, pop, data, ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise);

  result = 0.0;

  for (i = 0; i < ldata->n_used; i++)
  {
    const complex double chi_L = ldata->chi_L[i];
    const complex double chi_i = self->apply_shear_inv (g, chi_L);
    const gdouble x_i          = gsl_pow_2 (creal (chi_i)) + gsl_pow_2 (cimag (chi_i));
    const gdouble p_pop        = nc_galaxy_shape_pop_eval_p (pop, data->pop_data, x_i) / M_PI;
    const gdouble jac          = self->det_jac (g, chi_L);

    result += ldata->eff_weight[i] * p_pop * jac;
  }

  /* Purely the empty-domain / underflow floor -- see the class docs. */
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

