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
 * Fixed-node quadrature evaluation of the intrinsic-ellipticity marginal.
 *
 * Evaluates
 * $$
 * P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - f_g(\chi_I);
 *   \sigma_\mathrm{noise}^2\big)
 * $$
 * exactly (no series truncation in $g$), like #NcGalaxyShapeFactorQuad, via one of two
 * fixed quadrature schemes chosen from $\epsilon_\mathrm{obs}$/$\sigma_\mathrm{noise}$
 * alone, never $g$:
 *
 * - **Two-panel $\psi$** (default). Quad's $\psi=f_h^{-1}(\chi_L)$
 *   reparametrization ($h$: $f_h(0)=\epsilon_\mathrm{obs}$), radially split
 *   into $[0,R_\sigma)$/$[R_\sigma,1)$ Gauss-Legendre panels (see
 *   _exact_r_sigma()), with Quad's puncture correction.
 * - **Native $\chi_I$-polar**. Used when $1+\lvert\epsilon_\mathrm{obs}\rvert\leq
 *   N_\sigma\sigma_\mathrm{noise}$ ($N_\sigma=8$, see _use_chi_i_native()): integrates
 *   $\chi_I$'s own polar coordinates directly, no reparametrization, no Jacobian, no
 *   puncture correction (see _marginal_chi_i_native()).
 *
 * Both are exact; the switch only picks which grid resolves the integrand.
 * #NcGalaxyShapeFactorFixedQuad:n-radial/:n-angular size whichever grid is chosen
 * (two-panel: nodes PER PANEL; native: nodes in the one grid). The per-galaxy domain
 * depends only on $\epsilon_\mathrm{obs}$/$\sigma_\mathrm{noise}$, cached and reused
 * across every $g$.
 *
 * Works for any population; a fixed grid cannot resolve one narrower than
 * its node spacing ($\sigma_\mathrm{pop}\lesssim0.05$) -- use Quad instead.
 *
 * #NcGalaxyShapeFactorFixedQuad:use-marginal-spline additionally caches
 * $\ln P(\epsilon_\mathrm{obs}\mid g)$ over
 * $[-g_\mathrm{max},g_\mathrm{max}]^2$
 * (#NcGalaxyShapeFactorFixedQuad:spline-g-max), built lazily via autoknots
 * 2D splines for populations bounded at $r=0$, or a fixed knot grid
 * (_build_g_spline_fixed_knots()) for populations that diverge there.
 *
 * See docs/theory/wl_shape_factor_history.md for the design rationale.
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
#include <stdio.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_EXPM1_SAFE_BOUND (1.0e-1)
#define NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_NSIGMA (8.0)
#define NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_R_SIGMA_MAX (0.98)
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
  /* det_jac is the direct (non-log) Jacobian: this class sums linearly.
   * apply_shear builds the two-panel domain (chi_L=f_h(psi)) and the
   * native domain's per-g chi_L=apply_shear(g,chi_I); shear_at_origin
   * computes h. No apply_shear_inv pointer: _marginal_two_panel()'s hot
   * loop calls nc_wl_ellipticity.h's fused, non-pointer _trace/_trace_det
   * kernels directly, and _exact_r_sigma() gets its inverse from
   * apply_shear(-h,.) (Mobius inverse identity). */
  complex double (*apply_shear) (complex double g, complex double chi_i);
  complex double (*shear_at_origin) (complex double target);

  gdouble (*det_jac) (complex double g, complex double chi_L);

  NcGalaxyWLObsEllipConv ellip_conv;

  guint n_radial;
  guint n_angular;
  guint n_max; /* 2*n_radial*n_angular: two-panel uses all of it, native chi_I half */

  gboolean use_marginal_spline;
  gdouble spline_g_max;
  gdouble spline_rel_err;
} NcGalaxyShapeFactorFixedQuadPrivate;

enum
{
  PROP_0,
  PROP_N_RADIAL,
  PROP_N_ANGULAR,
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

  self->apply_shear     = NULL;
  self->shear_at_origin = NULL;
  self->det_jac         = NULL;
  self->ellip_conv      = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE;
  self->n_radial        = 0;
  self->n_angular       = 0;
  self->n_max           = 0;

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
        self->apply_shear     = &nc_wl_ellipticity_apply_shear_trace;
        self->shear_at_origin = &nc_wl_ellipticity_shear_at_origin_trace;
        self->det_jac         = &nc_wl_ellipticity_det_jac_trace;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->apply_shear     = &nc_wl_ellipticity_apply_shear_trace_det;
        self->shear_at_origin = &nc_wl_ellipticity_shear_at_origin_trace_det;
        self->det_jac         = &nc_wl_ellipticity_det_jac_trace_det;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }

    self->ellip_conv = ellip_conv;
    self->n_max      = 2 * self->n_radial * self->n_angular;
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

typedef enum
{
  NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_TWO_PANEL,
  NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_CHI_I_NATIVE,
} NcGalaxyShapeFactorFixedQuadMethod;

typedef struct _NcGalaxyShapeFactorFixedQuadData
{
  gboolean cache_valid;
  gdouble cached_epsilon_obs_1;
  gdouble cached_epsilon_obs_2;
  gdouble cached_std_noise;
  NcGalaxyShapeFactorFixedQuadMethod cached_method;

  guint n_used;
  complex double *base; /* two-panel: chi_L=f_h(psi_i); native chi_I: chi_I_i itself */
  gdouble *weight;      /* two-panel: quadrature_weight*|det J_{f_h}|; native chi_I: quadrature_weight/(2*pi) */
  gdouble *jac;         /* two-panel only: |det J_{f_g^-1}| per node, recomputed every call */
  GArray *x_arr;        /* two-panel: r_i, recomputed every call; native chi_I: r_i, fixed, filled once */
  GArray *p_arr;        /* eval_p_array()'s output, P_pop(r_i) */

  /* :use-marginal-spline cache. Every cached node, either branch, has a
   * strictly finite image under any g, so there is no singular-node
   * safety gate to track. */
  gboolean g_spline_valid;
  guint64 g_spline_pop_pkey;
  NcmSpline2d *g_spline; /* ln(marginal) over [-spline_g_max,spline_g_max]^2 */

  /* TRUE iff a spline was actually built for the current pop-pkey epoch;
   * FALSE only when eval_p() isn't well-defined near r=0. */
  gboolean g_spline_built;
} NcGalaxyShapeFactorFixedQuadData;

static void
_nc_galaxy_shape_factor_fixed_quad_ldata_destroy (gpointer p)
{
  NcGalaxyShapeFactorFixedQuadData *ldata = (NcGalaxyShapeFactorFixedQuadData *) p;

  ncm_spline2d_clear (&ldata->g_spline);

  g_free (ldata->base);
  g_free (ldata->weight);
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
  ldata->base        = g_new (complex double, self->n_max);
  ldata->weight      = g_new (gdouble, self->n_max);
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
}

static inline gdouble
_noise_val (complex double delta, gdouble sig2)
{
  const gdouble d2 = gsl_pow_2 (creal (delta)) + gsl_pow_2 (cimag (delta));

  return exp (-d2 / (2.0 * sig2)) / (2.0 * M_PI * sig2);
}

/* Controls whether we should integrate over chi_I instead of the transformed variable
 * $\psi$.
 *
 * TRUE iff the whole unit chi_L-disc lies within NSIGMA*std_noise of eps_obs -- the
 * farthest disc point from eps_obs is the diametrically opposite boundary point, at
 * distance 1+|eps_obs|.
 *
 * Also TRUE whenever |eps_obs|>=1. A noise-corrupted observed distortion/ellipticity
 * CAN legitimately land outside the unit disc (the noise model is not itself
 * disc-bounded).
 */
static inline gboolean
_use_chi_i_native (complex double eps_obs, gdouble std_noise)
{
  if (cabs (eps_obs) >= 1.0)
    return TRUE;

  return (1.0 + cabs (eps_obs)) <= NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_NSIGMA * std_noise;
}

/* Exact center/radius of the circle through 3 points: Mobius maps send circles to
 * circles exactly, so circumscribing 3 mapped boundary points recovers the true image
 * circle losslessly (see _exact_r_sigma()). */
static void
_circumcircle (complex double z0, complex double z1, complex double z2, complex double *center, gdouble *radius)
{
  const gdouble x1 = creal (z0), y1 = cimag (z0);
  const gdouble x2 = creal (z1), y2 = cimag (z1);
  const gdouble x3 = creal (z2), y3 = cimag (z2);
  const gdouble d  = 2.0 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
  const gdouble ux = ((x1 * x1 + y1 * y1) * (y2 - y3) + (x2 * x2 + y2 * y2) * (y3 - y1) + (x3 * x3 + y3 * y3) * (y1 - y2)) / d;
  const gdouble uy = ((x1 * x1 + y1 * y1) * (x3 - x2) + (x2 * x2 + y2 * y2) * (x1 - x3) + (x3 * x3 + y3 * y3) * (x2 - x1)) / d;

  *center = ux + I * uy;
  *radius = cabs (*center - z0);
}

/* Radius (centered at psi=0) containing the exact psi-preimage of the physical noise
 * disc {chi_L : |chi_L-eps_obs|<=NSIGMA*std_noise}: three boundary points mapped
 * through f_h^{-1}=apply_shear(-h,.) and circumscribed. Capped at R_SIGMA_MAX so the
 * outer panel keeps positive width.
 */
static gdouble
_exact_r_sigma (NcGalaxyShapeFactorFixedQuadPrivate * const self, complex double eps_obs, gdouble std_noise, complex double h)
{
  static const gdouble thetas[3] = { 0.0, 2.0 * M_PI / 3.0, 4.0 * M_PI / 3.0 };
  const gdouble target           = NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_NSIGMA * std_noise;
  const complex double neg_h     = -h;
  complex double psi[3];
  complex double center;
  gdouble radius;
  guint k;

  for (k = 0; k < 3; k++)
  {
    const complex double chi_L = eps_obs + target * cexp (I * thetas[k]);

    psi[k] = self->apply_shear (neg_h, chi_L);
  }

  _circumcircle (psi[0], psi[1], psi[2], &center, &radius);

  return fmin (NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_R_SIGMA_MAX, cabs (center) + radius);
}

/* Builds the two-panel psi-polar grid, g-independent: base[i]=f_h(psi_i),
 * weight[i]=quadrature_weight*|det J_{f_h}(psi_i)|. Radial coordinate uses two
 * Gauss-Legendre panels, [0,R_sigma) and [R_sigma,1); angular is equally-spaced,
 * offset by phi=arg(eps_obs) for rotation covariance. Every node has rho<1 strictly,
 * so every chi_L has |chi_L|<1 strictly (f_h is a disc automorphism).
 */
static void
_regen_domain_two_panel (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapeFactorFixedQuadData *ldata,
                         gdouble epsilon_obs_1, gdouble epsilon_obs_2, gdouble std_noise)
{
  const complex double eps_obs             = epsilon_obs_1 + I * epsilon_obs_2;
  const complex double h                   = self->shear_at_origin (eps_obs);
  const complex double neg_h               = -h;
  const gdouble phi                        = carg (eps_obs);
  const gdouble r_sigma                    = _exact_r_sigma (self, eps_obs, std_noise, h);
  const gdouble panel_lo[2]                = { 0.0, r_sigma };
  const gdouble panel_hi[2]                = { r_sigma, 1.0 };
  const guint n_panel_nodes[2]             = { self->n_radial, 15 };
  gsl_integration_glfixed_table *table1    = gsl_integration_glfixed_table_alloc (n_panel_nodes[0]);
  gsl_integration_glfixed_table *table2    = gsl_integration_glfixed_table_alloc (n_panel_nodes[1]);
  gsl_integration_glfixed_table *tables[2] = { table1, table2 };
  const gdouble w_theta                    = 2.0 * M_PI / self->n_angular;
  guint idx                                = 0;
  guint panel, i, j;

  for (panel = 0; panel < 2; panel++)
  {
    for (i = 0; i < n_panel_nodes[panel]; i++)
    {
      gdouble rho, wr;

      gsl_integration_glfixed_point (panel_lo[panel], panel_hi[panel], i, &rho, &wr, tables[panel]);

      for (j = 0; j < self->n_angular; j++)
      {
        const gdouble theta        = phi + w_theta * j;
        const complex double psi   = rho * cexp (I * theta);
        const complex double chi_L = self->apply_shear (h, psi);
        const gdouble jac_psi      = self->det_jac (neg_h, psi);

        ldata->base[idx]   = chi_L;
        ldata->weight[idx] = wr * rho * w_theta * jac_psi;
        idx++;
      }
    }
  }

  gsl_integration_glfixed_table_free (table1);
  gsl_integration_glfixed_table_free (table2);

  ldata->n_used               = idx;
  ldata->cached_method        = NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_TWO_PANEL;
  ldata->cached_epsilon_obs_1 = epsilon_obs_1;
  ldata->cached_epsilon_obs_2 = epsilon_obs_2;
  ldata->cached_std_noise     = std_noise;
  ldata->cache_valid          = TRUE;
  ldata->g_spline_valid       = FALSE;
}

/* Builds the native chi_I-polar grid: fixed nodes chi_I_i=rho_i*e^{i theta_i},
 * weight[i]=quadrature_weight/(2*pi) -- NO rho_i factor, since the polar measure's own
 * rho and the population's area density P_pop(rho)/(2*pi*rho) share a rho that cancels
 * analytically (this grid is centered on chi_I=0 itself). r_i is filled into x_arr
 * once, since it is g-independent (no inverse map in this branch).
 */
static void
_regen_domain_chi_i_native (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapeFactorFixedQuadData *ldata,
                            gdouble epsilon_obs_1, gdouble epsilon_obs_2, gdouble std_noise)
{
  const complex double eps_obs         = epsilon_obs_1 + I * epsilon_obs_2;
  const gdouble phi                    = carg (eps_obs);
  gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc (self->n_radial);
  const gdouble w_theta                = 2.0 * M_PI / self->n_angular;
  guint idx                            = 0;
  guint i, j;

  g_array_set_size (ldata->x_arr, self->n_radial * self->n_angular);

  {
    gdouble * const r_data = (gdouble *) ldata->x_arr->data;

    for (i = 0; i < self->n_radial; i++)
    {
      gdouble rho, wr;

      gsl_integration_glfixed_point (0.0, 1.0, i, &rho, &wr, table);

      for (j = 0; j < self->n_angular; j++)
      {
        const gdouble theta = phi + w_theta * j;

        ldata->base[idx]   = rho * cexp (I * theta);
        ldata->weight[idx] = wr * w_theta / (2.0 * M_PI);
        r_data[idx]        = rho;
        idx++;
      }
    }
  }

  gsl_integration_glfixed_table_free (table);

  ldata->n_used               = idx;
  ldata->cached_method        = NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_CHI_I_NATIVE;
  ldata->cached_epsilon_obs_1 = epsilon_obs_1;
  ldata->cached_epsilon_obs_2 = epsilon_obs_2;
  ldata->cached_std_noise     = std_noise;
  ldata->cache_valid          = TRUE;
  ldata->g_spline_valid       = FALSE;
}

/* Auto-switch domain build, used by eval_marginal()/eval_ln_marginal().
 * nc_galaxy_shape_factor_fixed_quad_eval_two_panel()/_eval_chi_i_native() call the two
 * branch-specific builders directly instead, forcing a branch regardless of this
 * switch.
 */
static void
_regen_domain (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapeFactorFixedQuadData *ldata,
               gdouble epsilon_obs_1, gdouble epsilon_obs_2, gdouble std_noise)
{
  const complex double eps_obs = epsilon_obs_1 + I * epsilon_obs_2;

  if (_use_chi_i_native (eps_obs, std_noise))
    _regen_domain_chi_i_native (self, ldata, epsilon_obs_1, epsilon_obs_2, std_noise);
  else
    _regen_domain_two_panel (self, ldata, epsilon_obs_1, epsilon_obs_2, std_noise);
}

static inline gboolean
_domain_matches_method (NcGalaxyShapeFactorFixedQuadData *ldata, gdouble epsilon_obs_1, gdouble epsilon_obs_2, gdouble std_noise,
                        NcGalaxyShapeFactorFixedQuadMethod method)
{
  return ldata->cache_valid && (ldata->cached_epsilon_obs_1 == epsilon_obs_1) &&
         (ldata->cached_epsilon_obs_2 == epsilon_obs_2) && (ldata->cached_std_noise == std_noise) &&
         (ldata->cached_method == method);
}

/* Also checks cached_method, not just eps_obs/std_noise: a domain last built by
 * eval_two_panel()/eval_chi_i_native() forcing the OTHER branch must still be detected
 * as stale. */
static inline gboolean
_domain_matches_auto (NcGalaxyShapeFactorFixedQuadData *ldata, gdouble epsilon_obs_1, gdouble epsilon_obs_2, gdouble std_noise)
{
  const complex double eps_obs                    = epsilon_obs_1 + I * epsilon_obs_2;
  const NcGalaxyShapeFactorFixedQuadMethod method = _use_chi_i_native (eps_obs, std_noise) ?
                                                    NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_CHI_I_NATIVE :
                                                    NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_TWO_PANEL;

  return _domain_matches_method (ldata, epsilon_obs_1, epsilon_obs_2, std_noise, method);
}

/* Raw (un-floored) two-panel marginal, given a valid TWO_PANEL domain. Sums the
 * puncture-correction bracket [N(eps_obs-chi_L)-N0], N0=N(eps_obs-chi_L0),
 * chi_L0=f_g(0), then adds N0 back once via integral(P_pop)=1. Mutates ldata's scratch
 * buffers (x_arr/jac/p_arr): safe for repeated calls, not reentrant.
 */
static gdouble
_marginal_two_panel (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                     NcGalaxyShapeFactorFixedQuadData *ldata, const complex double g)
{
  const complex double eps_obs = ldata->cached_epsilon_obs_1 + I * ldata->cached_epsilon_obs_2;
  const gdouble sig2           = gsl_pow_2 (ldata->cached_std_noise);
  gdouble alpha_o;
  guint i;

  /* x_i is known for every node before eval_p() is called, so batch through
   * eval_p_array() instead of one-at-a-time vfunc calls. x_arr holds x_i=|chi_I|^2
   * from the kernels, sqrt()'d in place into r_i. */
  g_array_set_size (ldata->x_arr, ldata->n_used);

  {
    gdouble * const x_data = (gdouble *) ldata->x_arr->data;

    /* Branch on ellip_conv once per call, not per node: each loop calls
     * nc_wl_ellipticity.h's fused, non-pointer kernel directly. */
    if (self->ellip_conv == NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE)
    {
      NcWLEllipticityTraceKernelPrep prep;

      nc_wl_ellipticity_trace_kernel_prepare (g, &prep);

      for (i = 0; i < ldata->n_used; i++)
        nc_wl_ellipticity_trace_kernel_apply (&prep, ldata->base[i], &x_data[i], &ldata->jac[i]);
    }
    else
    {
      for (i = 0; i < ldata->n_used; i++)
        nc_wl_ellipticity_trace_det_kernel (g, ldata->base[i], &x_data[i], &ldata->jac[i]);
    }

    for (i = 0; i < ldata->n_used; i++)
      x_data[i] = sqrt (x_data[i]);
  }

  nc_galaxy_shape_pop_eval_p_array (pop, data->pop_data, ldata->x_arr, &ldata->p_arr);
  alpha_o = nc_galaxy_shape_pop_exponent_at_origin (pop);

  /* Non-divergent case or too narrow noise disk */
  if ((alpha_o >= 1.0) || (ldata->cached_std_noise < 0.05))
  {
    const gdouble * const p_data = (const gdouble *) ldata->p_arr->data;
    const gdouble * const r_data = (const gdouble *) ldata->x_arr->data;
    gdouble sum                  = 0.0;

    for (i = 0; i < ldata->n_used; i++)
    {
      const complex double delta_i = eps_obs - ldata->base[i];
      const gdouble p_2d           = p_data[i] / (2.0 * M_PI * r_data[i]);
      const gdouble noise_i        = _noise_val (delta_i, sig2);

      sum += ldata->weight[i] * ldata->jac[i] * p_2d * noise_i;
    }

    return sum;
  }
  else
  {
    const gdouble * const p_data = (const gdouble *) ldata->p_arr->data;
    const gdouble * const r_data = (const gdouble *) ldata->x_arr->data;
    const complex double chi_L0  = self->apply_shear (g, 0.0);
    const gdouble d0             = gsl_pow_2 (creal (eps_obs - chi_L0)) + gsl_pow_2 (cimag (eps_obs - chi_L0));
    const gdouble N0             = _noise_val (eps_obs - chi_L0, sig2);
    gdouble corr_sum             = 0.0;

    for (i = 0; i < ldata->n_used; i++)
    {
      const complex double delta_i = eps_obs - ldata->base[i];
      const gdouble d_i            = gsl_pow_2 (creal (delta_i)) + gsl_pow_2 (cimag (delta_i));
      const gdouble delta_ratio    = (d0 - d_i) / (2.0 * sig2);
      const gdouble p_2d           = p_data[i] / (2.0 * M_PI * r_data[i]);
      gdouble bracket;

      /* expm1 avoids cancellation near delta_ratio=0; falls back to a direct,
       * individually-bounded difference where expm1 would overflow instead.
       */
      if ((fabs (delta_ratio) <= NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_EXPM1_SAFE_BOUND))
      {
        bracket = N0 * expm1 (delta_ratio);
      }
      else
      {
        const gdouble noise_i = _noise_val (delta_i, sig2);

        bracket = noise_i - N0;
      }

      corr_sum += ldata->weight[i] * ldata->jac[i] * p_2d * bracket;
    }

    return fabs (N0 + corr_sum);
  }
}

/* Verbose, node-by-node re-evaluation of _marginal_two_panel(), for
 * _direct_marginal_at_g()'s failure branch only: dumps every quantity that
 * feeds the sum (including which expm1/direct branch each node took) to
 * stderr right before the g_error() abort. Never called on the hot path.
 */
static gdouble
_marginal_two_panel_debug (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                           NcGalaxyShapeFactorFixedQuadData *ldata, const complex double g)
{
  const complex double eps_obs = ldata->cached_epsilon_obs_1 + I * ldata->cached_epsilon_obs_2;
  const gdouble sig2           = gsl_pow_2 (ldata->cached_std_noise);
  const complex double chi_L0  = self->apply_shear (g, 0.0);
  const gdouble d0             = gsl_pow_2 (creal (eps_obs - chi_L0)) + gsl_pow_2 (cimag (eps_obs - chi_L0));
  const gdouble N0             = _noise_val (eps_obs - chi_L0, sig2);
  gdouble corr_sum             = 0.0;
  guint i;

  ncm_model_params_print_all (NCM_MODEL (pop), stdout);

  fprintf (stderr, "---- _marginal_two_panel_debug: eps_obs=(% .15g,% .15g) sig2=% .15g g=(% .15g,% .15g) n_used=%u\n",
           creal (eps_obs), cimag (eps_obs), sig2, creal (g), cimag (g), ldata->n_used);
  fprintf (stderr, "     chi_L0=(% .15g,% .15g) d0=% .15g N0=% .15g SNR=% .15g\n", creal (chi_L0), cimag (chi_L0), d0, N0, cabs (eps_obs - chi_L0) / ldata->cached_std_noise);

  g_array_set_size (ldata->x_arr, ldata->n_used);

  {
    gdouble * const x_data = (gdouble *) ldata->x_arr->data;

    if (self->ellip_conv == NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE)
    {
      NcWLEllipticityTraceKernelPrep prep;

      nc_wl_ellipticity_trace_kernel_prepare (g, &prep);

      for (i = 0; i < ldata->n_used; i++)
        nc_wl_ellipticity_trace_kernel_apply (&prep, ldata->base[i], &x_data[i], &ldata->jac[i]);
    }
    else
    {
      for (i = 0; i < ldata->n_used; i++)
        nc_wl_ellipticity_trace_det_kernel (g, ldata->base[i], &x_data[i], &ldata->jac[i]);
    }

    for (i = 0; i < ldata->n_used; i++)
      x_data[i] = sqrt (x_data[i]);
  }

  nc_galaxy_shape_pop_eval_p_array (pop, data->pop_data, ldata->x_arr, &ldata->p_arr);

  {
    const gdouble * const p_data = (const gdouble *) ldata->p_arr->data;
    const gdouble * const r_data = (const gdouble *) ldata->x_arr->data;

    for (i = 0; i < ldata->n_used; i++)
    {
      const complex double delta_i = eps_obs - ldata->base[i];
      const gdouble d_i            = gsl_pow_2 (creal (delta_i)) + gsl_pow_2 (cimag (delta_i));
      const gdouble delta_ratio    = (d0 - d_i) / (2.0 * sig2);
      const gdouble p_2d           = p_data[i] / (2.0 * M_PI * r_data[i]);
      const gboolean via_expm1     = (fabs (delta_ratio) <= NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_EXPM1_SAFE_BOUND);
      gdouble bracket;
      gdouble contrib;


      if (via_expm1)
      {
        bracket = N0 * expm1 (delta_ratio);
      }
      else
      {
        const gdouble noise_i = _noise_val (delta_i, sig2);

        bracket = noise_i - N0;
      }

      contrib = ldata->weight[i] * ldata->jac[i] * p_2d * bracket;

      fprintf (stderr,
               "  [%4u] chi_L=(% .15g,% .15g) r=% .15g weight=% .15g jac=% .15g p=% .15g p_2d=% .15g d_i=% .15g "
               "delta_ratio=% .15g via=%s bracket=% .15g contrib=% .15g running_corr_sum=% .15g\n",
               i,
               creal (ldata->base[i]),
               cimag (ldata->base[i]),
               r_data[i],
               ldata->weight[i],
               ldata->jac[i],
               p_data[i],
               p_2d,
               d_i,
               delta_ratio,
               via_expm1 ? "expm1" : "direct",
               bracket,
               contrib,
               corr_sum + contrib);

      corr_sum += contrib;
    }
  }

  fprintf (stderr, "---- _marginal_two_panel_debug: N0=% .15g corr_sum=% .15g result=% .15g\n", N0, corr_sum, N0 + corr_sum);

  return N0 + corr_sum;
}

/* Raw (un-floored) native-chi_I marginal, given a valid CHI_I_NATIVE domain. Plain
 * forward evaluation, no bracket: this grid is centered on chi_I=0, so there is no
 * unrelated singular point to protect against. r_i is g-independent, already filled by
 * _regen_domain_chi_i_native().
 */
static gdouble
_marginal_chi_i_native (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                        NcGalaxyShapeFactorFixedQuadData *ldata, const complex double g)
{
  const complex double eps_obs = ldata->cached_epsilon_obs_1 + I * ldata->cached_epsilon_obs_2;
  const gdouble sig2           = gsl_pow_2 (ldata->cached_std_noise);
  gdouble total                = 0.0;
  guint i;

  nc_galaxy_shape_pop_eval_p_array (pop, data->pop_data, ldata->x_arr, &ldata->p_arr);

  {
    const gdouble * const p_data = (const gdouble *) ldata->p_arr->data;

    for (i = 0; i < ldata->n_used; i++)
    {
      const complex double chi_L = self->apply_shear (g, ldata->base[i]);
      const gdouble noise_i      = _noise_val (eps_obs - chi_L, sig2);

      total += ldata->weight[i] * p_data[i] * noise_i;
    }
  }

  return total;
}

/* Verbose, node-by-node re-evaluation of _marginal_chi_i_native(), for
 * _direct_marginal_at_g()'s failure branch only: dumps every quantity that
 * feeds the sum to stderr right before the g_error() abort. Never called on
 * the hot path.
 */
static gdouble
_marginal_chi_i_native_debug (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                              NcGalaxyShapeFactorFixedQuadData *ldata, const complex double g)
{
  const complex double eps_obs = ldata->cached_epsilon_obs_1 + I * ldata->cached_epsilon_obs_2;
  const gdouble sig2           = gsl_pow_2 (ldata->cached_std_noise);
  gdouble total                = 0.0;
  guint i;

  fprintf (stderr, "---- _marginal_chi_i_native_debug: eps_obs=(% .15g,% .15g) sig2=% .15g g=(% .15g,% .15g) n_used=%u\n",
           creal (eps_obs), cimag (eps_obs), sig2, creal (g), cimag (g), ldata->n_used);

  nc_galaxy_shape_pop_eval_p_array (pop, data->pop_data, ldata->x_arr, &ldata->p_arr);

  {
    const gdouble * const p_data = (const gdouble *) ldata->p_arr->data;
    const gdouble * const r_data = (const gdouble *) ldata->x_arr->data;

    for (i = 0; i < ldata->n_used; i++)
    {
      const complex double chi_L = self->apply_shear (g, ldata->base[i]);
      const gdouble noise_i      = _noise_val (eps_obs - chi_L, sig2);
      const gdouble contrib      = ldata->weight[i] * p_data[i] * noise_i;

      fprintf (stderr,
               "  [%4u] chi_I=(% .15g,% .15g) r=% .15g weight=% .15g p=% .15g chi_L=(% .15g,% .15g) noise=% .15g "
               "contrib=% .15g running_total=% .15g\n",
               i, creal (ldata->base[i]), cimag (ldata->base[i]), r_data[i], ldata->weight[i], p_data[i],
               creal (chi_L), cimag (chi_L), noise_i, contrib, total + contrib);

      total += contrib;
    }
  }

  fprintf (stderr, "---- _marginal_chi_i_native_debug: total=% .15g\n", total);

  return total;
}

/* Ground truth for use-marginal-spline's builders below. A non-finite result is always
 * a bug (invalid input the domain build should have rejected, or a genuine defect) --
 * never a deep-tail underflow, which shows up as a finite, small/negative value
 * instead (see the MIN_MARGINAL floor at this function's callers).
 */
static gdouble
_direct_marginal_at_g (NcGalaxyShapeFactorFixedQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                       NcGalaxyShapeFactorFixedQuadData *ldata, const complex double g)
{
  gdouble result;

  switch (ldata->cached_method)
  {
    case NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_CHI_I_NATIVE:
      result = _marginal_chi_i_native (self, pop, data, ldata, g);
      break;

    case NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_TWO_PANEL:
    default:
      result = _marginal_two_panel (self, pop, data, ldata, g);
      break;
  }

  if (!isfinite (result))
  {
    switch (ldata->cached_method)
    {
      case NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_CHI_I_NATIVE:
        _marginal_chi_i_native_debug (self, pop, data, ldata, g);
        break;

      case NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_TWO_PANEL:
      default:
        _marginal_two_panel_debug (self, pop, data, ldata, g);
        break;
    }

    g_error ("nc_galaxy_shape_factor_fixed_quad: non-finite marginal at g=(% .6g,% .6g), "
             "eps_obs=(% .6g,% .6g), std_noise=% .6g, method=%s, result=% .6g.",
             creal (g), cimag (g), ldata->cached_epsilon_obs_1, ldata->cached_epsilon_obs_2, ldata->cached_std_noise,
             (ldata->cached_method == NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_CHI_I_NATIVE) ? "chi_i_native" : "two_panel",
             result);
  }

  return result > 0.0 ? result : 1.0e-300;
}

/* Shared args behind _build_g_spline()'s Fx/Fy gsl_function slices,
 * following ncm_spline2d_set_function()'s own contract. */
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

  return log (v);
}

/* Builds ldata's g-spline: ln(marginal) as a bivariate function of
 * (g_1,g_2) over [-spline_g_max,spline_g_max]^2. Ln-space is required --
 * the marginal spans many orders of magnitude over this square. Node
 * placement uses this codebase's "autoknots" 2D spline machinery
 * (ncm_spline2d_set_function()): knots placed adaptively from two 1D
 * slices (_g_spline_slice_func(), each at the other coordinate fixed to
 * 0) to spline_rel_err, then the true 2D surface is filled in at the
 * resulting tensor grid. Node count is data-dependent.
 *
 * The r=0 divergence check below only matters for the TWO_PANEL branch --
 * the native chi_I branch's per-g function is always smooth by
 * construction. The check is population-only, so a chi_I-native galaxy
 * with a divergent population is still (conservatively) routed to the
 * fixed-knots path. */
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

      ncm_matrix_set (zm, i, j, log (v));
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

  /* The area density the two-panel branch needs, P_pop(r)/(2*pi*r),
   * diverges at r=0 whenever eval_p(r) vanishes slower than r itself
   * (e.g. Beta with alpha<2, eval_p(r)~r^(alpha-1)). eval_p(0.0) cannot
   * detect this directly (it collapses to exactly 0.0 for any alpha>1),
   * so probe the local power-law exponent via two small, well-separated
   * r instead: the area density diverges iff the log-log slope between
   * them is below 1.
   *
   * Every two-panel domain node has some g mapping it to chi_I=0, so a
   * divergent population turns each node into its own spike in g-space;
   * the adaptive autoknots build cannot resolve that and aborts via
   * g_error, so this population gets the fixed-knots grid instead, which
   * cannot abort. */
  {
    const gdouble r1             = 1.0e-6;
    const gdouble r2             = 1.0e-3;
    const gdouble p1             = nc_galaxy_shape_pop_eval_p (pop, data->pop_data, r1);
    const gdouble p2             = nc_galaxy_shape_pop_eval_p (pop, data->pop_data, r2);
    const gboolean well_defined  = gsl_finite (p1) && gsl_finite (p2) && (p1 > 0.0) && (p2 > 0.0);
    const gboolean adaptive_safe = well_defined && (log (p2 / p1) / log (r2 / r1) >= 1.0 - 1.0e-4);

    if (!well_defined)
    {
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

      ncm_matrix_set (zm, i, j, log (v));
    }
  }

  ncm_spline2d_prepare (ldata->g_spline);
}

/* Dispatch for :use-marginal-spline: falls back to direct computation
 * outside the cached box, rebuilds the spline if invalid or stale. */
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

  if (!_domain_matches_auto (ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise))
    _regen_domain (self, ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise);

  result = self->use_marginal_spline ?
           _eval_marginal_spline (self, pop, data, ldata, g) :
           _direct_marginal_at_g (self, pop, data, ldata, g);

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
   * Number of fixed Gauss-Legendre radial nodes: PER PANEL for the
   * two-panel psi branch, or in the single grid for the native chi_I
   * branch. Default 21.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_RADIAL,
                                   g_param_spec_uint ("n-radial",
                                                      "Number of radial nodes",
                                                      "Number of fixed Gauss-Legendre nodes in the radial direction",
                                                      1, G_MAXUINT, 21,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorFixedQuad:n-angular:
   *
   * Number of equally-spaced angular nodes of whichever grid is chosen.
   * Default 21.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_ANGULAR,
                                   g_param_spec_uint ("n-angular",
                                                      "Number of angular nodes",
                                                      "Number of angular quadrature nodes",
                                                      1, G_MAXUINT, 21,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapeFactorFixedQuad:use-marginal-spline:
   *
   * When %TRUE, caches the marginal as a function of g, per galaxy, over
   * $[-g_\mathrm{max},g_\mathrm{max}]^2$
   * (#NcGalaxyShapeFactorFixedQuad:spline-g-max); g outside the box always
   * falls back to the exact direct computation. Default %FALSE.
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
   * Half-side of the square #NcGalaxyShapeFactorFixedQuad:use-marginal-spline's
   * cache covers. Default 0.3.
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
   * autoknots build. Default 1e-4.
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

/**
 * nc_galaxy_shape_factor_fixed_quad_eval_two_panel:
 * @gsffq: a #NcGalaxyShapeFactorFixedQuad
 * @pop: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: reduced shear, real component
 * @g_2: reduced shear, imaginary component
 * @epsilon_obs_1: observed ellipticity/distortion, real component
 * @epsilon_obs_2: observed ellipticity/distortion, imaginary component
 *
 * Evaluates the marginal using the two-panel psi branch, regardless of
 * what the class' own switch (_use_chi_i_native()) would pick. Not used
 * by any production path -- eval_marginal()/eval_ln_marginal() always use
 * the switch; this and
 * nc_galaxy_shape_factor_fixed_quad_eval_chi_i_native() exist to test each
 * branch in isolation.
 *
 * Returns: $P(\epsilon_\mathrm{obs} \mid g)$, via the two-panel psi branch.
 */
gdouble
nc_galaxy_shape_factor_fixed_quad_eval_two_panel (NcGalaxyShapeFactorFixedQuad *gsffq, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                  const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (gsffq);
  NcGalaxyShapeFactorFixedQuadData *ldata          = (NcGalaxyShapeFactorFixedQuadData *) data->ldata;
  const complex double g                           = g_1 + I * g_2;

  if (!_domain_matches_method (ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise, NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_TWO_PANEL))
    _regen_domain_two_panel (self, ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise);

  return _direct_marginal_at_g (self, pop, data, ldata, g);
}

/**
 * nc_galaxy_shape_factor_fixed_quad_eval_chi_i_native:
 * @gsffq: a #NcGalaxyShapeFactorFixedQuad
 * @pop: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: reduced shear, real component
 * @g_2: reduced shear, imaginary component
 * @epsilon_obs_1: observed ellipticity/distortion, real component
 * @epsilon_obs_2: observed ellipticity/distortion, imaginary component
 *
 * Evaluates the marginal using the native chi_I branch, regardless of
 * what the class' own switch would pick. See
 * nc_galaxy_shape_factor_fixed_quad_eval_two_panel()'s docs.
 *
 * Returns: $P(\epsilon_\mathrm{obs} \mid g)$, via the native chi_I branch.
 */
gdouble
nc_galaxy_shape_factor_fixed_quad_eval_chi_i_native (NcGalaxyShapeFactorFixedQuad *gsffq, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                     const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (gsffq);
  NcGalaxyShapeFactorFixedQuadData *ldata          = (NcGalaxyShapeFactorFixedQuadData *) data->ldata;
  const complex double g                           = g_1 + I * g_2;

  if (!_domain_matches_method (ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise, NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_CHI_I_NATIVE))
    _regen_domain_chi_i_native (self, ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise);

  return _direct_marginal_at_g (self, pop, data, ldata, g);
}

/**
 * nc_galaxy_shape_factor_fixed_quad_peek_domain:
 * @gsffq: a #NcGalaxyShapeFactorFixedQuad
 * @pop: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapeFactorData
 * @g_1: reduced shear, real component
 * @g_2: reduced shear, imaginary component
 * @weights: (out) (transfer full): per-node quadrature weight times the
 * noise-kernel value at @g_1/@g_2 -- NOT the node's full contribution to
 * the marginal (that also needs the population density and, for
 * two-panel, the inverse-map Jacobian and the puncture-correction
 * baseline), so summing @weights does not reproduce eval_marginal()'s
 * result
 *
 * Diagnostic accessor exposing @data's own per-galaxy quadrature domain,
 * mapped to chi_L-space at @g_1/@g_2. Rebuilds the domain first if
 * @data's epsilon_obs/std_noise changed (same cache check as
 * eval_marginal()'s own, i.e. the auto-switch, never a forced branch).
 * Not used by any production evaluation path -- added for `numcosmo
 * inspect galaxy-shape-integrand` and similar diagnostics.
 *
 * The two-panel branch's nodes are already positioned in chi_L-space
 * (g-independent), so @g_1/@g_2 only affect the returned weight there; the
 * native chi_I branch's nodes are fixed in chi_I-space, so @g_1/@g_2 are
 * needed to place them in chi_L-space at all.
 *
 * Returns: (transfer full): an n_used x 2 #NcmMatrix of chi_L node
 * positions, column 0 = Re(chi_L), column 1 = Im(chi_L).
 */
NcmMatrix *
nc_galaxy_shape_factor_fixed_quad_peek_domain (NcGalaxyShapeFactorFixedQuad *gsffq, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                               const gdouble g_1, const gdouble g_2, NcmVector **weights)
{
  NcGalaxyShapeFactorFixedQuadPrivate * const self = nc_galaxy_shape_factor_fixed_quad_get_instance_private (gsffq);
  NcGalaxyShapeFactorFixedQuadData *ldata          = (NcGalaxyShapeFactorFixedQuadData *) data->ldata;
  const gdouble epsilon_obs_1                      = data->epsilon_obs_1;
  const gdouble epsilon_obs_2                      = data->epsilon_obs_2;
  const complex double g                           = g_1 + I * g_2;
  gboolean chi_i_native;
  NcmMatrix *chi_L_mat;
  complex double eps_obs;
  gdouble sig2;
  guint i;

  if (!_domain_matches_auto (ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise))
    _regen_domain (self, ldata, epsilon_obs_1, epsilon_obs_2, data->std_noise);

  eps_obs      = ldata->cached_epsilon_obs_1 + I * ldata->cached_epsilon_obs_2;
  sig2         = gsl_pow_2 (ldata->cached_std_noise);
  chi_i_native = (ldata->cached_method == NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_METHOD_CHI_I_NATIVE);

  chi_L_mat = ncm_matrix_new (ldata->n_used, 2);
  *weights  = ncm_vector_new (ldata->n_used);

  for (i = 0; i < ldata->n_used; i++)
  {
    const complex double chi_L = chi_i_native ? self->apply_shear (g, ldata->base[i]) : ldata->base[i];
    const gdouble noise_i      = _noise_val (eps_obs - chi_L, sig2);

    ncm_matrix_set (chi_L_mat, i, 0, creal (chi_L));
    ncm_matrix_set (chi_L_mat, i, 1, cimag (chi_L));
    ncm_vector_set (*weights, i, ldata->weight[i] * noise_i);
  }

  return chi_L_mat;
}

