/***************************************************************************
 *            nc_xcor_kernel.c
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcXcorKernel:
 *
 * Base object for the kernels of projected observables used in cross-correlations.
 *
 * The projected field and its kernel are linked by
 * \begin{equation}
 * $A(\hat{\mathbf{n}}) = \int_0^\infty dz \ W^A(z) \ \delta(\chi(z)\hat{\mathbf{n}}, z)$
 * \end{equation}
 * where $\delta$ is the matter density field.
 *
 * Kernels also implement the noise power spectrum.
 *
 * See <a href="../../theory/ncm/specfunc/sbessel_projection.html">UltraLevin: Non-Limber
 * Angular Power Spectra</a> for the pipeline a kernel drives: the adaptive
 * $k$ domain, the closure interpolating $W_\ell(k)$, and the error estimate built
 * from its interpolation error.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/powspec/ncm_powspec.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"
#include "ncm/specfunc/ncm_sbessel_ode_solver.h"
#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "ncm/stats/ncm_function_sample_set.h"
#include "nc/background/nc_distance.h"
#include "nc/xcor/nc_xcor_kernel.h"
#include "ncm/model/ncm_model_ctrl.h"
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/core/ncm_memory_pool.h"
#include "nc/xcor/nc_xcor_kernel_component.h"
#include "nc/xcor/nc_xcor.h"
#include "nc_enum_types.h"
#include "nc/xcor/kernel/nc_xcor_kernel_private.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_POWSPEC,
  PROP_INTEGRATOR,
  PROP_LMAX,
  PROP_L_LIMBER,
  PROP_ADAPTIVE_EPSILON,
  PROP_ADAPTIVE_BOUNDARY_TRIES,
  PROP_RELTOL,
  PROP_PEAK_EPSILON,
  PROP_MAX_BORDER_EXPANSIONS,
  PROP_MAX_ITER,
  PROP_EXPANSION_FACTOR,
  PROP_TRACK_CLOSURE_ERROR,
  PROP_PANEL_ORDER_CAP,
  PROP_PANELS_PER_EFOLD,
  PROP_PANEL_LEVEL_MIN,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcXcorKernel, nc_xcor_kernel, NCM_TYPE_MODEL)
G_DEFINE_BOXED_TYPE (NcXcorKinetic, nc_xcor_kinetic, nc_xcor_kinetic_copy, nc_xcor_kinetic_free)

static void
nc_xcor_kernel_init (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->dist                     = NULL;
  self->ps                       = NULL;
  self->sbi                      = NULL;
  self->cosmo_ctrl               = ncm_model_ctrl_new (NULL);
  self->prepared_pkey            = 0;
  self->outdated                 = TRUE;
  self->lmax                     = 0;
  self->l_limber                 = 0;
  self->adaptive_epsilon         = 0.0;
  self->adaptive_boundary_tries  = 0;
  self->reltol                   = 0.0;
  self->peak_epsilon             = 0.0;
  self->max_border_expansions    = 0;
  self->max_iter                 = 0;
  self->expansion_factor         = 0.0;
  self->panel_order_cap          = 0;
  self->panels_per_efold         = 0.0;
  self->panel_level_min          = 0;
  self->track_closure_error      = FALSE;
  self->tolerance_balance_warned = FALSE;
  self->constructed              = FALSE;
}

static void
_nc_xcor_kernel_dispose (GObject *object)
{
  NcXcorKernel *xclk        = NC_XCOR_KERNEL (object);
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  nc_distance_clear (&self->dist);
  ncm_powspec_clear (&self->ps);
  ncm_sbessel_integrator_clear (&self->sbi);
  ncm_model_ctrl_clear (&self->cosmo_ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_parent_class)->dispose (object);
}

static void
_nc_xcor_kernel_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_parent_class)->finalize (object);
}

static void
_nc_xcor_kernel_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_parent_class)->constructed (object);
  {
    NcXcorKernel *xclk        = NC_XCOR_KERNEL (object);
    NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

    if (self->dist == NULL)
      g_error ("nc_xcor_kernel_constructed: dist property was not set. "
               "The 'dist' property must be provided at construction time.");

    if (self->ps == NULL)
      g_error ("nc_xcor_kernel_constructed: powspec property was not set. "
               "The 'powspec' property must be provided at construction time.");

    if ((self->l_limber != 0) && (self->sbi == NULL))
      g_error ("nc_xcor_kernel_constructed: l_limber property is set to %d but "
               "integrator property was not set. "
               "The 'integrator' property must be provided at construction time "
               "to use the non-Limber method.",
               self->l_limber);

    nc_distance_compute_inv_comoving (self->dist, TRUE);
    nc_distance_require_zf (self->dist, 1.0e10);

    self->constructed = TRUE;
  }
}

static void
_nc_xcor_kernel_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernel *xclk        = NC_XCOR_KERNEL (object);
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_return_if_fail (NC_IS_XCOR_KERNEL (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_distance_clear (&self->dist);
      self->dist = g_value_dup_object (value);
      break;
    case PROP_POWSPEC:
      ncm_powspec_clear (&self->ps);
      self->ps = g_value_dup_object (value);
      break;
    case PROP_INTEGRATOR:
      ncm_sbessel_integrator_clear (&self->sbi);
      self->sbi = g_value_dup_object (value);
      break;
    case PROP_LMAX:
      nc_xcor_kernel_set_lmax (xclk, g_value_get_uint (value));
      break;
    case PROP_L_LIMBER:
      nc_xcor_kernel_set_l_limber (xclk, g_value_get_int (value));
      break;
    case PROP_ADAPTIVE_EPSILON:
      nc_xcor_kernel_set_adaptive_epsilon (xclk, g_value_get_double (value));
      break;
    case PROP_ADAPTIVE_BOUNDARY_TRIES:
      nc_xcor_kernel_set_adaptive_boundary_tries (xclk, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      nc_xcor_kernel_set_reltol (xclk, g_value_get_double (value));
      break;
    case PROP_PEAK_EPSILON:
      nc_xcor_kernel_set_peak_epsilon (xclk, g_value_get_double (value));
      break;
    case PROP_MAX_BORDER_EXPANSIONS:
      nc_xcor_kernel_set_max_border_expansions (xclk, g_value_get_uint (value));
      break;
    case PROP_MAX_ITER:
      nc_xcor_kernel_set_max_iter (xclk, g_value_get_uint (value));
      break;
    case PROP_EXPANSION_FACTOR:
      nc_xcor_kernel_set_expansion_factor (xclk, g_value_get_double (value));
      break;
    case PROP_TRACK_CLOSURE_ERROR:
      nc_xcor_kernel_set_track_closure_error (xclk, g_value_get_boolean (value));
      break;
    case PROP_PANEL_ORDER_CAP:
      nc_xcor_kernel_set_panel_order_cap (xclk, g_value_get_uint (value));
      break;
    case PROP_PANELS_PER_EFOLD:
      nc_xcor_kernel_set_panels_per_efold (xclk, g_value_get_double (value));
      break;
    case PROP_PANEL_LEVEL_MIN:
      nc_xcor_kernel_set_panel_level_min (xclk, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernel *xclk        = NC_XCOR_KERNEL (object);
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_return_if_fail (NC_IS_XCOR_KERNEL (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, self->dist);
      break;
    case PROP_POWSPEC:
      g_value_set_object (value, self->ps);
      break;
    case PROP_INTEGRATOR:
      g_value_set_object (value, self->sbi);
      break;
    case PROP_LMAX:
      g_value_set_uint (value, nc_xcor_kernel_get_lmax (xclk));
      break;
    case PROP_L_LIMBER:
      g_value_set_int (value, nc_xcor_kernel_get_l_limber (xclk));
      break;
    case PROP_ADAPTIVE_EPSILON:
      g_value_set_double (value, nc_xcor_kernel_get_adaptive_epsilon (xclk));
      break;
    case PROP_ADAPTIVE_BOUNDARY_TRIES:
      g_value_set_uint (value, nc_xcor_kernel_get_adaptive_boundary_tries (xclk));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_xcor_kernel_get_reltol (xclk));
      break;
    case PROP_PEAK_EPSILON:
      g_value_set_double (value, nc_xcor_kernel_get_peak_epsilon (xclk));
      break;
    case PROP_MAX_BORDER_EXPANSIONS:
      g_value_set_uint (value, nc_xcor_kernel_get_max_border_expansions (xclk));
      break;
    case PROP_MAX_ITER:
      g_value_set_uint (value, nc_xcor_kernel_get_max_iter (xclk));
      break;
    case PROP_EXPANSION_FACTOR:
      g_value_set_double (value, nc_xcor_kernel_get_expansion_factor (xclk));
      break;
    case PROP_TRACK_CLOSURE_ERROR:
      g_value_set_boolean (value, nc_xcor_kernel_get_track_closure_error (xclk));
      break;
    case PROP_PANEL_ORDER_CAP:
      g_value_set_uint (value, nc_xcor_kernel_get_panel_order_cap (xclk));
      break;
    case PROP_PANELS_PER_EFOLD:
      g_value_set_double (value, nc_xcor_kernel_get_panels_per_efold (xclk));
      break;
    case PROP_PANEL_LEVEL_MIN:
      g_value_set_uint (value, nc_xcor_kernel_get_panel_level_min (xclk));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

NCM_MSET_MODEL_REGISTER_ID (nc_xcor_kernel, NC_TYPE_XCOR_KERNEL);

/* LCOV_EXCL_START */

static void
_nc_xcor_kernel_get_z_range_not_implemented (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  g_error ("nc_xcor_kernel_get_z_range: get_z_range virtual method not implemented for %s",
           G_OBJECT_TYPE_NAME (xclk));
}

/* LCOV_EXCL_STOP */

static void
nc_xcor_kernel_class_init (NcXcorKernelClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_xcor_kernel_set_property;
  model_class->get_property = &_nc_xcor_kernel_get_property;
  object_class->constructed = &_nc_xcor_kernel_constructed;
  object_class->dispose     = &_nc_xcor_kernel_dispose;
  object_class->finalize    = &_nc_xcor_kernel_finalize;

  ncm_model_class_set_name_nick (model_class, "Cross-correlation Kernels", "xcor-kernel");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));

  /**
   * NcXcorKernel:dist:
   *
   * Distance object used to compute the comoving distance $\chi(z)$ and its inverse
   * $z(\chi)$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:powspec:
   *
   * Power spectrum object used to compute the cross-correlation.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_POWSPEC,
                                   g_param_spec_object ("powspec",
                                                        NULL,
                                                        "Power spectrum object",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:integrator:
   *
   * Spherical Bessel integrator object used to compute the non-Limber $W_\ell(k)$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_INTEGRATOR,
                                   g_param_spec_object ("integrator",
                                                        NULL,
                                                        "Spherical Bessel integrator object",
                                                        NCM_TYPE_SBESSEL_INTEGRATOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:lmax:
   *
   * Maximum multipole $\ell$ to compute $W_\ell(k)$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "Maximum multipole",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:l-limber:
   *
   * Limber approximation threshold $\ell_\mathrm{limber}$. The Limber approximation is
   * used for $\ell \ge \ell_\mathrm{limber}$, and the non-Limber method is used for
   * $\ell < \ell_\mathrm{limber}$. A value of 0 means that the Limber approximation is
   * used for all multipoles, and a value of -1 means that the non-Limber method is
   * used for all multipoles.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_L_LIMBER,
                                   g_param_spec_int ("l-limber",
                                                     NULL,
                                                     "Limber approximation threshold (-1: never, 0: always, N>0: use for l>=N)",
                                                     -1, G_MAXINT, 0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:adaptive-epsilon:
   *
   * Threshold $\epsilon$ below which $W_\ell(k)$ is considered negligible. With
   * $\Vert W(k) \Vert_2$ the norm over the multipoles of the block and
   * $W_{\max}$ its largest value over the $k$ already evaluated, a value of
   * $k$ satisfies the test when $\Vert W(k) \Vert_2 < \epsilon\, W_{\max}$.
   *
   * The test enters the construction of a closure at two points.
   *
   * - The domain. ncm_function_sample_set_expand_domain() extends each end of
   *   the $k$ range by the factor $1 \pm$ #NcXcorKernel:expansion-factor per
   *   step. An end is fixed after #NcXcorKernel:adaptive-boundary-tries
   *   consecutive values of $k$ satisfy the test, or when it reaches the $k$
   *   range of the components.
   * - The components. In the non-Limber case the test is applied to each
   *   component separately, with its own norm in place of $\Vert W(k)
   *   \Vert_2$. After the same number of consecutive values of $k$ satisfy
   *   it, the component is removed from the sum beyond that $k$. Components whose
   *   $\chi$ supports touch or overlap are tested on the norm of their sum and
   *   removed together. In the Chebyshev closure that $k$ is a panel edge.
   */
  g_object_class_install_property (object_class,
                                   PROP_ADAPTIVE_EPSILON,
                                   g_param_spec_double ("adaptive-epsilon",
                                                        NULL,
                                                        "Convergence threshold for adaptive k-range determination",
                                                        0.0, 1.0, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:adaptive-boundary-tries:
   *
   * Number of consecutive values of $k$ below the #NcXcorKernel:adaptive-epsilon
   * threshold before stopping the extension of the $k$ range or removing a component
   * from the sum.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ADAPTIVE_BOUNDARY_TRIES,
                                   g_param_spec_uint ("adaptive-boundary-tries",
                                                      NULL,
                                                      "Number of consecutive boundary points below threshold before stopping extension",
                                                      1, G_MAXUINT, 5,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:reltol:
   *
   * Relative tolerance of the closure of $W_\ell(k)$. The absolute tolerance $a$ is
   * #NcXcorKernel:peak-epsilon times the smallest of the block's peaks. The
   * norm $\Vert \cdot \Vert_2$ is taken over the multipoles of the block.
   *
   * - Spline closure. ncm_function_sample_set_adaptive_midpoint() accepts an interval
   *   when, at its midpoint, $\Vert W - \tilde W \Vert_2 \le \mathrm{reltol}\, \Vert W
   *   \Vert_2 + a$, with $\tilde W$ the spline. The two terms are added, so the larger
   *   of them determines the accuracy reached.
   * - Chebyshev closure. A panel is accepted when, for every multipole, the
   *   coefficient vectors $c_N$ and $c_{2N}$ at consecutive orders satisfy $\Vert
   *   c_{2N} - c_N \Vert_2 < \max (\mathrm{reltol}\, \Vert c_{2N} \Vert_2, a)$.
   *
   * Each computed value of $W_\ell(k)$ carries the relative error of the integrator,
   * and a closure is not more accurate than the values it interpolates. The
   * construction of a closure therefore aborts with g_error() when the larger of
   * @reltol and #NcXcorKernel:peak-epsilon is below the integrator's relative
   * tolerance. The spline closure warns once when @reltol and
   * #NcXcorKernel:peak-epsilon differ by more than two orders of magnitude.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance for adaptive midpoint refinement",
                                                        0.0, 1.0, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:peak-epsilon:
   *
   * Absolute tolerance of the closure of $W_\ell(k)$, as a fraction of the block's
   * smallest peak: $a = $ `peak-epsilon` $\times \min_\ell \max_k \vert W_\ell(k)
   * \vert$. The smallest peak is used so that the multipole of smallest amplitude in
   * the block is resolved to a fraction of its own peak. How $a$ enters each closure
   * is given under #NcXcorKernel:reltol.
   *
   * In the regions of $k$ where $a$ exceeds the relative tolerance, the closure is
   * accepted with an absolute error of order $a$ in $W_\ell$, independent of the value
   * of $W_\ell$ there. The quantity is a tolerance on $W_\ell(k)$, not on $C_\ell$:
   * the integrand of $C_\ell$ is $k^2 W_\ell^i W_\ell^j$, so the error in a product of
   * two closures is of order $a^2$ in those regions.
   *
   * nc_xcor_kernel_set_peak_epsilon() warns for values below
   * %NC_XCOR_KERNEL_MIN_USEFUL_PEAK_EPSILON.
   */
  g_object_class_install_property (object_class,
                                   PROP_PEAK_EPSILON,
                                   g_param_spec_double ("peak-epsilon",
                                                        NULL,
                                                        "Peak-relative floor of the adaptive refinement of the k-space closure",
                                                        GSL_DBL_MIN, 1.0, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:max-border-expansions:
   *
   * Maximum number of border expansion iterations. Each iteration extends the $k$
   * range by the factor $1 \pm$ #NcXcorKernel:expansion-factor.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_BORDER_EXPANSIONS,
                                   g_param_spec_uint ("max-border-expansions",
                                                      NULL,
                                                      "Maximum number of border expansion iterations",
                                                      1, G_MAXUINT, 500,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:max-iter:
   *
   * Maximum number of adaptive midpoint refinement iterations per interval.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_ITER,
                                   g_param_spec_uint ("max-iter",
                                                      NULL,
                                                      "Maximum number of adaptive midpoint refinement iterations",
                                                      1, G_MAXUINT, 10000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:track-closure-error:
   *
   * Whether the closure records its interpolation error $\delta W_\ell$ on each
   * interval on which it is a single polynomial: a knot interval of the spline
   * closure, a panel of the Chebyshev closure. nc_xcor_compute_full() propagates
   * this record to the error estimate of $C_\ell$. Without the record it uses
   * #NcXcorKernel:reltol and #NcXcorKernel:peak-epsilon in its place.
   *
   * Each closure records the quantity its acceptance test measured.
   *
   * - Spline closure. On the interval between two consecutive knots,
   *   $\delta W_\ell$ is the largest value of $\vert W_\ell - \tilde W_\ell \vert$
   *   measured at a midpoint of that interval during refinement; see
   *   ncm_function_sample_set_get_residuals().
   * - Chebyshev closure. On a panel accepted at order $N$ after the test at
   *   order $N/2$, $\delta W_\ell = \sum_{j \ge N/2} \vert c_j \vert$, the sum
   *   over the coefficients the order $N/2$ did not carry. The sum bounds the
   *   difference between the two expansions in the supremum norm.
   *
   * The record is one double per interval per multipole of the block.
   */
  g_object_class_install_property (object_class,
                                   PROP_TRACK_CLOSURE_ERROR,
                                   g_param_spec_boolean ("track-closure-error",
                                                         NULL,
                                                         "Whether to record the interpolation error the closure achieved",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:panel-order-cap:
   *
   * Highest level $k$ of the Chebyshev expansion tried on one panel, with $N = 2^k +
   * 1$ nodes and coefficients at level $k$. Zero selects the default, $k = 5$.
   *
   * A panel is expanded from #NcXcorKernel:panel-level-min upward, doubling the order
   * until the acceptance test of #NcXcorKernel:reltol passes. When the test fails at
   * the cap, or when the coefficients at the level below it predict that it will fail
   * there (see ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap()), the panel
   * is bisected and each half is expanded again from the minimum level. The nodes of
   * the parent panel are not nodes of its halves, so the values of $W_\ell$ computed
   * for a failed panel are not reused. A panel narrower than $10^{-6}$ of its upper
   * edge is kept at the cap without further bisection.
   */
  g_object_class_install_property (object_class,
                                   PROP_PANEL_ORDER_CAP,
                                   g_param_spec_uint ("panel-order-cap",
                                                      NULL,
                                                      "Highest Chebyshev order tried per panel before bisecting",
                                                      0, 12, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:panels-per-efold:
   *
   * Number of panels per e-fold of $u = k/\nu$ in the initial grid of the Limber
   * Chebyshev closure. The initial panel edges are the band edges of the components
   * together with a geometric grid in $u$ of this density over the closure's domain; a
   * grid point closer than a quarter of the grid spacing to another edge is dropped.
   * Each initial panel is then expanded and bisected where it does not converge, see
   * #NcXcorKernel:panel-order-cap. Zero disables the grid, and each segment between
   * band edges starts as one panel. A non-Limber closure takes no grid: its
   * $W_\ell(k)$ oscillates, and its panels start from the whole domain, cut at the
   * component boundaries.
   */
  g_object_class_install_property (object_class,
                                   PROP_PANELS_PER_EFOLD,
                                   g_param_spec_double ("panels-per-efold",
                                                        NULL,
                                                        "Panels per e-fold of k in the initial Chebyshev grid",
                                                        0.0, 1000.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:panel-level-min:
   *
   * Level $k$ at which the Chebyshev expansion of a panel starts, with $N = 2^k + 1$
   * nodes; the order is doubled from there until the acceptance test of
   * #NcXcorKernel:reltol passes or #NcXcorKernel:panel-order-cap is reached. The first
   * test compares the levels $k$ and $k + 1$. A low level accepts a panel on which
   * $W_\ell(k)$ is near zero or near linear with few nodes; the check of each accepted
   * panel against the values computed during the domain expansion guards against a
   * coarse first grid that misses a feature.
   */
  g_object_class_install_property (object_class,
                                   PROP_PANEL_LEVEL_MIN,
                                   g_param_spec_uint ("panel-level-min",
                                                      NULL,
                                                      "Level at which the Chebyshev expansion of a panel starts",
                                                      0, 12, 3,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:expansion-factor:
   *
   * Factor by which the $k$ range is extended at each iteration of the domain
   * expansion. The $k$ range is extended by the factor $1 \pm$
   * #NcXcorKernel:expansion-factor at each end of the range. The expansion stops when
   * #NcXcorKernel:adaptive-boundary-tries consecutive values of $k$ satisfy the
   * #NcXcorKernel:adaptive-epsilon test, or when the $k$ range reaches the $k$ range
   * of the components.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_EXPANSION_FACTOR,
                                   g_param_spec_double ("expansion-factor",
                                                        NULL,
                                                        "Expansion factor for domain extension",
                                                        0.0, 1.0, 0.2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_mset_model_register_id (model_class, "NcXcorKernel", "Cross-correlation Kernels",
                              NULL, TRUE, NCM_MSET_MODEL_MAIN);

  klass->get_z_range = &_nc_xcor_kernel_get_z_range_not_implemented;
}

/*
 * Returns the private data of @xclk, for the other files of the kernel
 * (nc_xcor_kernel_private.h).
 */
NcXcorKernelPrivate *
_nc_xcor_kernel_get_private (NcXcorKernel *xclk)
{
  return nc_xcor_kernel_get_instance_private (xclk);
}

/**
 * nc_xcor_kernel_ref:
 * @xclk: a #NcXcorKernel
 *
 * Increases the reference count of @xclk by one.
 *
 * Returns: (transfer full): @xclk.
 */
NcXcorKernel *
nc_xcor_kernel_ref (NcXcorKernel *xclk)
{
  return g_object_ref (xclk);
}

/**
 * nc_xcor_kernel_free:
 * @xclk: a #NcXcorKernel
 *
 * Decreases the reference count of @xclk by one. If the reference count
 * reaches zero, the object is freed.
 *
 */
void
nc_xcor_kernel_free (NcXcorKernel *xclk)
{
  g_object_unref (xclk);
}

/**
 * nc_xcor_kernel_clear:
 * @xclk: a #NcXcorKernel
 *
 * Atomically decrements the reference count of @xclk by one.
 * If the reference count drops to zero, all memory allocated by @xclk is
 * released. @xclk is set to NULL after being freed.
 *
 */
void
nc_xcor_kernel_clear (NcXcorKernel **xclk)
{
  g_clear_object (xclk);
}

/**
 * nc_xcor_kinetic_copy:
 * @xck: a #NcXcorKinetic
 *
 * Creates a copy of @xck.
 *
 * Returns: (transfer full): a new #NcXcorKinetic copy of @xck.
 */
NcXcorKinetic *
nc_xcor_kinetic_copy (NcXcorKinetic *xck)
{
  NcXcorKinetic *xck_copy = g_new (NcXcorKinetic, 1);

  xck_copy[0] = xck[0];

  return xck_copy;
}

/**
 * nc_xcor_kinetic_free:
 * @xck: a #NcXcorKinetic
 *
 * Frees @xck.
 *
 */
void
nc_xcor_kinetic_free (NcXcorKinetic *xck)
{
  g_free (xck);
}

/**
 * nc_xcor_kernel_obs_len: (virtual obs_len)
 * @xclk: a #NcXcorKernel
 *
 * Gets the number of observables required by this kernel.
 *
 * Returns: the number of observables
 */
guint
nc_xcor_kernel_obs_len (NcXcorKernel *xclk)
{
  return NC_XCOR_KERNEL_GET_CLASS (xclk)->obs_len (xclk);
}

/**
 * nc_xcor_kernel_obs_params_len: (virtual obs_params_len)
 * @xclk: a #NcXcorKernel
 *
 * Gets the number of parameters needed to describe the observables
 * for this kernel (e.g., measurement uncertainties, systematic parameters).
 *
 * Returns: the number of observable parameters
 */
guint
nc_xcor_kernel_obs_params_len (NcXcorKernel *xclk)
{
  return NC_XCOR_KERNEL_GET_CLASS (xclk)->obs_params_len (xclk);
}

/**
 * nc_xcor_kernel_get_z_range: (virtual get_z_range)
 * @xclk: a #NcXcorKernel
 * @zmin: (out): minimum redshift
 * @zmax: (out): maximum redshift
 * @zmid: (out) (allow-none): mid redshift
 *
 * Get the redshift range of the kernel. This is a virtual method that
 * must be implemented by subclasses.
 *
 */
void
nc_xcor_kernel_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  NC_XCOR_KERNEL_GET_CLASS (xclk)->get_z_range (xclk, zmin, zmax, zmid);
}

/**
 * nc_xcor_kernel_peek_dist:
 * @xclk: a #NcXcorKernel
 *
 * Peeks the distance object from the kernel. This method is intended
 * for use by subclass implementations.
 *
 * Returns: (transfer none): the distance object.
 */
NcDistance *
nc_xcor_kernel_peek_dist (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->dist;
}

/**
 * nc_xcor_kernel_peek_powspec:
 * @xclk: a #NcXcorKernel
 *
 * Peeks the power spectrum object from the kernel. This method is intended
 * for use by subclass implementations.
 *
 * Returns: (transfer none): the power spectrum object.
 */
NcmPowspec *
nc_xcor_kernel_peek_powspec (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->ps;
}

/**
 * nc_xcor_kernel_peek_integrator:
 * @xclk: a #NcXcorKernel
 *
 * Peeks the spherical Bessel integrator object from the kernel. This method is
 * intended for use by subclass implementations. Returns NULL if no integrator is set.
 *
 * Returns: (transfer none) (nullable): the spherical Bessel integrator object or NULL.
 */
NcmSBesselIntegrator *
nc_xcor_kernel_peek_integrator (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->sbi;
}

/**
 * nc_xcor_kernel_get_k_range:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @l: multipole
 * @kmin: (out): minimum wavenumber
 * @kmax: (out): maximum wavenumber
 *
 * Gets the valid k range for the kernel at multipole @l.
 * Uses the component-based implementation.
 */
void
nc_xcor_kernel_get_k_range (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax)
{
  NcXcorKernelClass *klass = NC_XCOR_KERNEL_GET_CLASS (xclk);
  GPtrArray *comp_list     = klass->get_component_list (xclk);
  gdouble global_kmin      = 0.0;
  gdouble global_kmax      = G_MAXDOUBLE;
  const gdouble nu         = l + 0.5;
  guint i;

  if ((comp_list == NULL) || (comp_list->len == 0))
  {
    if (comp_list != NULL)
      g_ptr_array_unref (comp_list);

    g_error ("nc_xcor_kernel_get_k_range: kernel %s returned empty component list",
             G_OBJECT_TYPE_NAME (xclk));

    return;
  }

  for (i = 0; i < comp_list->len; i++)
  {
    NcXcorKernelComponent *comp = g_ptr_array_index (comp_list, i);
    gdouble chi_min, chi_max, k_min, k_max;

    nc_xcor_kernel_component_get_limits (comp, cosmo, &chi_min, &chi_max, &k_min, &k_max);

    {
      const gdouble k_min_limb = nu / chi_max;
      const gdouble k_max_limb = nu / chi_min;

      k_min = GSL_MAX (k_min, k_min_limb);
      k_max = GSL_MIN (k_max, k_max_limb);
    }

    global_kmin = GSL_MAX (global_kmin, k_min);
    global_kmax = GSL_MIN (global_kmax, k_max);
  }

  g_ptr_array_unref (comp_list);

  *kmin = global_kmin;
  *kmax = global_kmax;
}

/**
 * nc_xcor_kernel_get_eval:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @l: multipole
 * @closure_type: how to represent $W_\ell(k)$, see #NcXcor:closure-type
 *
 * Gets an evaluation function for the kernel at multipole @l.
 * Convenience wrapper around nc_xcor_kernel_get_eval_vectorized() for a single multipole.
 *
 * Returns: (transfer full): the evaluation function for the kernel.
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_get_eval (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, NcXcorKernelClosure closure_type)
{
  return nc_xcor_kernel_get_eval_vectorized (xclk, cosmo, l, l, closure_type);
}

/**
 * nc_xcor_kernel_get_eval_vectorized:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 * @closure_type: how to represent $W_\ell(k)$, see #NcXcor:closure-type
 *
 * Gets a vectorized evaluation function for the kernel over a range of multipoles.
 * The returned integrand will have len = lmax - lmin + 1, and will evaluate all
 * multipoles in the range [lmin, lmax] simultaneously.
 *
 * Uses the base class implementation which checks the l-limber property:
 *
 * - If lmin >= l_limber (or l_limber == 0), uses component-based Limber approximation
 * - If l_limber < 0, use the non-Limber method
 * - Otherwise falls back to single-l get_eval for lmin
 *
 * Returns: (transfer full): the vectorized evaluation function for the kernel.
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_get_eval_vectorized (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcXcorKernelClosure closure_type)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return nc_xcor_kernel_get_eval_vectorized_full (xclk, cosmo, lmin, lmax, self->sbi, closure_type);
}

/**
 * nc_xcor_kernel_get_eval_vectorized_full:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 * @sbi: (nullable): the #NcmSBesselIntegrator to use, or %NULL for @xclk's own
 * @closure_type: how to represent $W_\ell(k)$, see #NcXcor:closure-type
 *
 * Same as nc_xcor_kernel_get_eval_vectorized(), but integrates with @sbi
 * instead of the kernel's `integrator` property.
 *
 * A #NcmSBesselIntegratorLevin holds reusable state tied to one multipole
 * range, so a caller computing several ell blocks does better keeping one
 * integrator per block and passing it in here than letting every kernel carry
 * its own. Passing the integrator rather than storing it also keeps @xclk free
 * of per-call state, so one kernel can be evaluated for several blocks at once
 * as long as each gets its own @sbi.
 *
 * @sbi is unused when the block is taken under Limber (see
 * #NcXcorKernel:l-limber), where no spherical Bessel integral is performed.
 *
 * Returns: (transfer full): the kernel integrand over [@lmin, @lmax]
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_get_eval_vectorized_full (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi, NcXcorKernelClosure closure_type)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  if ((self->l_limber == 0) || ((self->l_limber > 0) && (lmin >= self->l_limber)))
    return _nc_xcor_kernel_build_limber_integrand (xclk, cosmo, lmin, lmax, closure_type);
  else
    return _nc_xcor_kernel_build_non_limber_integrand (xclk, cosmo, lmin, lmax, sbi, closure_type);
}

/**
 * nc_xcor_kernel_get_lmax:
 * @xclk: a #NcXcorKernel
 *
 * Gets the maximum multipole for the kernel.
 *
 * Returns: the maximum multipole
 */
guint
nc_xcor_kernel_get_lmax (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->lmax;
}

/**
 * nc_xcor_kernel_set_lmax:
 * @xclk: a #NcXcorKernel
 * @lmax: the maximum multipole
 *
 * Sets the maximum multipole for the kernel.
 *
 */
void
nc_xcor_kernel_set_lmax (NcXcorKernel *xclk, guint lmax)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->lmax = lmax;
}

/**
 * nc_xcor_kernel_get_l_limber:
 * @xclk: a #NcXcorKernel
 *
 * Gets the Limber approximation threshold for the kernel.
 * Returns -1 for never using Limber, 0 for always using Limber,
 * or N > 0 to use Limber for l >= N.
 *
 * Returns: the Limber threshold
 */
gint
nc_xcor_kernel_get_l_limber (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->l_limber;
}

/**
 * nc_xcor_kernel_set_l_limber:
 * @xclk: a #NcXcorKernel
 * @l_limber: the Limber threshold (-1: never, 0: always, N>0: use for l>=N)
 *
 * Sets the Limber approximation threshold for the kernel.
 *
 */
void
nc_xcor_kernel_set_l_limber (NcXcorKernel *xclk, gint l_limber)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  if ((self->constructed) && (l_limber != 0) && (self->sbi == NULL))
    g_error ("nc_xcor_kernel_set_l_limber: cannot set l_limber to %d "
             "for kernel %s because no integrator is set. "
             "The 'integrator' property must be provided to use the non-Limber method.",
             l_limber, G_OBJECT_TYPE_NAME (xclk));

  self->l_limber = l_limber;
}

/**
 * nc_xcor_kernel_get_adaptive_epsilon:
 * @xclk: a #NcXcorKernel
 *
 * Gets the convergence threshold for adaptive k-range determination in the
 * non-Limber integrand. The algorithm stops extending the k range when all
 * component contributions drop below epsilon times the maximum kernel value.
 *
 * Returns: the adaptive epsilon value
 */
gdouble
nc_xcor_kernel_get_adaptive_epsilon (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->adaptive_epsilon;
}

/**
 * nc_xcor_kernel_set_adaptive_epsilon:
 * @xclk: a #NcXcorKernel
 * @adaptive_epsilon: the convergence threshold (must be > 0)
 *
 * Sets the convergence threshold for adaptive k-range determination in the
 * non-Limber integrand. Typical values range from 1e-4 to 1e-8, with smaller
 * values providing more accurate integration at the cost of more evaluations.
 *
 */
void
nc_xcor_kernel_set_adaptive_epsilon (NcXcorKernel *xclk, gdouble adaptive_epsilon)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert (adaptive_epsilon > 0.0);
  self->adaptive_epsilon = adaptive_epsilon;
}

/**
 * nc_xcor_kernel_get_adaptive_boundary_tries:
 * @xclk: a #NcXcorKernel
 *
 * Gets the number of consecutive boundary points that must be below the
 * convergence threshold before stopping boundary extension. This helps
 * avoid false positives where a single low point prematurely stops the
 * adaptive k-range determination.
 *
 * Returns: the number of required consecutive tries
 */
guint
nc_xcor_kernel_get_adaptive_boundary_tries (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->adaptive_boundary_tries;
}

/**
 * nc_xcor_kernel_set_adaptive_boundary_tries:
 * @xclk: a #NcXcorKernel
 * @adaptive_boundary_tries: the number of consecutive tries (must be >= 1)
 *
 * Sets the number of consecutive boundary points that must be below the
 * convergence threshold before stopping boundary extension. Higher values
 * provide more robust convergence detection at the cost of additional
 * function evaluations. Typical values range from 2 to 5.
 *
 */
void
nc_xcor_kernel_set_adaptive_boundary_tries (NcXcorKernel *xclk, guint adaptive_boundary_tries)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert (adaptive_boundary_tries >= 1);
  self->adaptive_boundary_tries = adaptive_boundary_tries;
}

/**
 * nc_xcor_kernel_get_reltol:
 * @xclk: a #NcXcorKernel
 *
 * Gets the relative tolerance used for adaptive midpoint refinement in the
 * non-Limber integrand construction.
 *
 * Returns: the relative tolerance value
 */
gdouble
nc_xcor_kernel_get_reltol (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->reltol;
}

/**
 * nc_xcor_kernel_set_reltol:
 * @xclk: a #NcXcorKernel
 * @reltol: the relative tolerance (must be > 0)
 *
 * Sets the relative tolerance for adaptive midpoint refinement. Smaller values
 * provide more accurate spline interpolation at the cost of more knots.
 * Typical values range from 1e-4 to 1e-8.
 *
 */
void
nc_xcor_kernel_set_reltol (NcXcorKernel *xclk, gdouble reltol)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert (reltol > 0.0);
  self->reltol = reltol;
}

/**
 * nc_xcor_kernel_get_peak_epsilon:
 * @xclk: a #NcXcorKernel
 *
 *
 */
gdouble
nc_xcor_kernel_get_peak_epsilon (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->peak_epsilon;
}

/**
 * nc_xcor_kernel_set_peak_epsilon:
 * @xclk: a #NcXcorKernel
 * @peak_epsilon: the absolute tolerance as a fraction of the peak (must be > 0)
 *
 * Sets #NcXcorKernel:peak-epsilon. Values below
 * %NC_XCOR_KERNEL_MIN_USEFUL_PEAK_EPSILON are accepted with a warning.
 */
void
nc_xcor_kernel_set_peak_epsilon (NcXcorKernel *xclk, gdouble peak_epsilon)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert_cmpfloat (peak_epsilon, >, 0.0);

  if (peak_epsilon < NC_XCOR_KERNEL_MIN_USEFUL_PEAK_EPSILON)
    g_warning ("nc_xcor_kernel_set_peak_epsilon: %.3e is below the useful floor of %.0e. "
               "This floor is measured against the peak of W(k), but the C_l integrand "
               "is k^2 W_a W_b, so it enters squared: %.3e here is %.3e on the integrand, "
               "past what the outer integral carries. It cannot improve the result and can "
               "cost orders of magnitude in spline knots.",
               peak_epsilon, NC_XCOR_KERNEL_MIN_USEFUL_PEAK_EPSILON,
               peak_epsilon, peak_epsilon * peak_epsilon);

  self->peak_epsilon = peak_epsilon;
}

/**
 * nc_xcor_kernel_get_max_border_expansions:
 * @xclk: a #NcXcorKernel
 *
 * Gets the maximum number of border expansion iterations allowed during domain
 * extension in the non-Limber integrand construction.
 *
 * Returns: the maximum number of border expansions
 */
guint
nc_xcor_kernel_get_max_border_expansions (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->max_border_expansions;
}

/**
 * nc_xcor_kernel_set_max_border_expansions:
 * @xclk: a #NcXcorKernel
 * @max_border_expansions: the maximum number of expansions (must be >= 1)
 *
 * Sets the maximum number of border expansion iterations. Higher values allow
 * the domain to extend further when needed, at the cost of potentially more
 * function evaluations. Typical values range from 1000 to 10000.
 *
 */
void
nc_xcor_kernel_set_max_border_expansions (NcXcorKernel *xclk, guint max_border_expansions)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert (max_border_expansions >= 1);
  self->max_border_expansions = max_border_expansions;
}

/**
 * nc_xcor_kernel_get_max_iter:
 * @xclk: a #NcXcorKernel
 *
 * Gets the maximum number of adaptive midpoint refinement iterations allowed
 * in the non-Limber integrand construction.
 *
 * Returns: the maximum number of iterations
 */
guint
nc_xcor_kernel_get_max_iter (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->max_iter;
}

/**
 * nc_xcor_kernel_set_max_iter:
 * @xclk: a #NcXcorKernel
 * @max_iter: the maximum number of iterations (must be >= 1)
 *
 * Sets the maximum number of adaptive midpoint refinement iterations. Higher values
 * allow for more refinement passes when needed. Typical values range from 1000 to
 * 100000.
 *
 */
void
nc_xcor_kernel_set_max_iter (NcXcorKernel *xclk, guint max_iter)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert (max_iter >= 1);
  self->max_iter = max_iter;
}

/**
 * nc_xcor_kernel_get_expansion_factor:
 * @xclk: a #NcXcorKernel
 *
 * Gets the expansion factor used for domain extension in the non-Limber integrand
 * construction. This determines how much the domain is extended in each iteration.
 *
 * Returns: the expansion factor
 */
gdouble
nc_xcor_kernel_get_expansion_factor (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->expansion_factor;
}

/**
 * nc_xcor_kernel_set_expansion_factor:
 * @xclk: a #NcXcorKernel
 * @expansion_factor: the expansion factor (must be > 0 and < 1)
 *
 * Sets the expansion factor for domain extension. Larger values result in more
 * aggressive expansion. Typical values range from 0.1 to 0.5.
 *
 */
void
nc_xcor_kernel_set_expansion_factor (NcXcorKernel *xclk, gdouble expansion_factor)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert (expansion_factor > 0.0 && expansion_factor < 1.0);
  self->expansion_factor = expansion_factor;
}

/**
 * nc_xcor_kernel_get_panel_order_cap:
 * @xclk: a #NcXcorKernel
 *
 * Returns: the panel order cap, or 0 for the default. See
 * #NcXcorKernel:panel-order-cap.
 */
guint
nc_xcor_kernel_get_panel_order_cap (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->panel_order_cap;
}

/**
 * nc_xcor_kernel_set_panel_order_cap:
 * @xclk: a #NcXcorKernel
 * @panel_order_cap: the cap, or 0 for the default
 *
 * Sets #NcXcorKernel:panel-order-cap. Read when a closure is built, so one already
 * built keeps the panels it was built with.
 *
 */
void
nc_xcor_kernel_set_panel_order_cap (NcXcorKernel *xclk, guint panel_order_cap)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->panel_order_cap = panel_order_cap;
}

/**
 * nc_xcor_kernel_get_panels_per_efold:
 * @xclk: a #NcXcorKernel
 *
 * Returns: the number of panels per e-fold of $k$ in the initial Chebyshev grid, or
 * 0.0 for no grid. See #NcXcorKernel:panels-per-efold.
 */
gdouble
nc_xcor_kernel_get_panels_per_efold (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->panels_per_efold;
}

/**
 * nc_xcor_kernel_set_panels_per_efold:
 * @xclk: a #NcXcorKernel
 * @panels_per_efold: panels per e-fold of $k$, or 0.0 for no grid
 *
 * Sets #NcXcorKernel:panels-per-efold. Read when a closure is built, so one already
 * built keeps the panels it was built with.
 */
void
nc_xcor_kernel_set_panels_per_efold (NcXcorKernel *xclk, gdouble panels_per_efold)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert_cmpfloat (panels_per_efold, >=, 0.0);
  self->panels_per_efold = panels_per_efold;
}

/**
 * nc_xcor_kernel_get_panel_level_min:
 * @xclk: a #NcXcorKernel
 *
 * Returns: the level at which the Chebyshev expansion of a panel starts. See
 * #NcXcorKernel:panel-level-min.
 */
guint
nc_xcor_kernel_get_panel_level_min (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->panel_level_min;
}

/**
 * nc_xcor_kernel_set_panel_level_min:
 * @xclk: a #NcXcorKernel
 * @panel_level_min: the starting level
 *
 * Sets #NcXcorKernel:panel-level-min. Read when a closure is built, so one already
 * built keeps the panels it was built with. Aborts when the level is not below
 * #NcXcorKernel:panel-order-cap at build time.
 */
void
nc_xcor_kernel_set_panel_level_min (NcXcorKernel *xclk, guint panel_level_min)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->panel_level_min = panel_level_min;
}

/**
 * nc_xcor_kernel_get_track_closure_error:
 * @xclk: a #NcXcorKernel
 *
 * Returns: whether the closure records its interpolation error. See
 * #NcXcorKernel:track-closure-error.
 */
gboolean
nc_xcor_kernel_get_track_closure_error (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->track_closure_error;
}

/**
 * nc_xcor_kernel_set_track_closure_error:
 * @xclk: a #NcXcorKernel
 * @track_closure_error: whether to record the interpolation error
 *
 * Sets #NcXcorKernel:track-closure-error. It is read when a closure is built, so a
 * closure already built keeps whatever it was built with.
 *
 */
void
nc_xcor_kernel_set_track_closure_error (NcXcorKernel *xclk, gboolean track_closure_error)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->track_closure_error = track_closure_error;
}

/**
 * nc_xcor_kernel_eval_limber_z: (virtual eval_limber_z)
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @z: a #gdouble
 * @xck: a #NcXcorKinetic
 * @l: a #gint
 *
 * Evaluates the Limber kernel at redshift @z for multipole @l. The kinetic quantities
 * (comoving distance and Hubble parameter) are provided in @xck. Returns zero if @z is
 * outside the kernel's redshift range.
 *
 * Returns: the kernel value $W(z,\ell)$
 */
gdouble
nc_xcor_kernel_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  return NC_XCOR_KERNEL_GET_CLASS (xclk)->eval_limber_z (xclk, cosmo, z, xck, l);
}

/**
 * nc_xcor_kernel_eval_limber_z_prefactor:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @l: a #gint
 *
 * Evaluates the Limber approximation redshift-dependent prefactor for multipole @l.
 *
 * Returns: the Limber redshift prefactor.
 */
gdouble
nc_xcor_kernel_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  return NC_XCOR_KERNEL_GET_CLASS (xclk)->eval_limber_z_prefactor (xclk, cosmo, l);
}

/**
 * nc_xcor_kernel_eval_limber_z_full:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @z: a #gdouble
 * @dist: a #NcDistance
 * @l: a #gint
 *
 * Evaluates the Limber kernel at redshift @z for multipole @l, including the
 * normalization factor. This function computes the kinetic quantities internally using
 * @dist and applies the kernel's constant factor.
 *
 * Returns: the normalized kernel value $c \times W(z,\ell)$
 */
gdouble
nc_xcor_kernel_eval_limber_z_full (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, NcDistance *dist, gint l)
{
  const gdouble chi_z     = nc_distance_comoving (dist, cosmo, z); /* in units of Hubble radius */
  const gdouble E_z       = nc_hicosmo_E (cosmo, z);
  const NcXcorKinetic xck = { chi_z, E_z, z };
  const gdouble prefactor = nc_xcor_kernel_eval_limber_z_prefactor (xclk, cosmo, l);

  return NC_XCOR_KERNEL_GET_CLASS (xclk)->eval_limber_z (xclk, cosmo, z, &xck, l) * prefactor;
}

/**
 * nc_xcor_kernel_add_noise: (virtual add_noise)
 * @xclk: a #NcXcorKernel
 * @vp1: a #NcmVector
 * @vp2: a #NcmVector
 * @lmin: a #guint
 *
 * vp2 = vp1 + noise spectrum
 *
 */
void
nc_xcor_kernel_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NC_XCOR_KERNEL_GET_CLASS (xclk)->add_noise (xclk, vp1, vp2, lmin);
}

/**
 * nc_xcor_kernel_prepare: (virtual prepare)
 * @xclk: a #NcXcorKernel
 * @cosmo: a NcHICosmo
 *
 * Prepares the kernel for evaluation with the given cosmological model. This may
 * involve precomputing quantities that depend on @cosmo but not on redshift or
 * multipole.
 *
 */
void
nc_xcor_kernel_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  NC_XCOR_KERNEL_GET_CLASS (xclk)->prepare (xclk, cosmo);

  ncm_model_ctrl_update (self->cosmo_ctrl, NCM_MODEL (cosmo));
  self->prepared_pkey = ncm_model_state_get_pkey (NCM_MODEL (xclk));
  self->outdated      = FALSE;
}

/**
 * nc_xcor_kernel_prepare_if_needed:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 *
 * Calls nc_xcor_kernel_prepare() only when something the prepared state depends on has
 * changed since the last preparation: the cosmology (tracked through a #NcmModelCtrl),
 * the kernel's own parameters (its #NcmModel pkey), or a change announced with
 * nc_xcor_kernel_mark_outdated(). Repeated solves at one cosmology then pay nothing
 * here, which is what a sampler that keeps its kernels and its #NcXcorSolver across
 * steps relies on.
 */
void
nc_xcor_kernel_prepare_if_needed (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);
  const gboolean cosmo_up   = ncm_model_ctrl_update (self->cosmo_ctrl, NCM_MODEL (cosmo));
  const gboolean self_up    = ncm_model_state_get_pkey (NCM_MODEL (xclk)) != self->prepared_pkey;

  if (cosmo_up || self_up || self->outdated)
    nc_xcor_kernel_prepare (xclk, cosmo);
}

/**
 * nc_xcor_kernel_mark_outdated:
 * @xclk: a #NcXcorKernel
 *
 * Announces that data the kernel is built on changed outside its parameters (a
 * replaced window table, say), so the next nc_xcor_kernel_prepare_if_needed() prepares
 * again even at an unchanged cosmology.
 */
void
nc_xcor_kernel_mark_outdated (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->outdated = TRUE;
}

/**
 * nc_xcor_kernel_get_component_list: (virtual get_component_list)
 * @xclk: a #NcXcorKernel
 *
 * Gets the list of components that make up this kernel.
 *
 * Returns: (transfer container) (element-type NcXcorKernelComponent): a #GPtrArray of
 * #NcXcorKernelComponent
 */
GPtrArray *
nc_xcor_kernel_get_component_list (NcXcorKernel *xclk)
{
  return NC_XCOR_KERNEL_GET_CLASS (xclk)->get_component_list (xclk);
}

static void
_nc_xcor_kernel_log_all_models_go (GType model_type, guint n)
{
  guint nc, i, j;
  GType *models = g_type_children (model_type, &nc);

  for (i = 0; i < nc; i++)
  {
    guint ncc;
    GType *model_sc = g_type_children (models[i], &ncc);

    g_message ("#  ");

    for (j = 0; j < n; j++)
      g_message (" ");

    g_message ("%s\n", g_type_name (models[i]));

    if (ncc)
      _nc_xcor_kernel_log_all_models_go (models[i], n + 2);

    g_free (model_sc);
  }

  g_free (models);
}

/**
 * nc_xcor_kernel_log_all_models:
 *
 * Logs all registered #NcXcorLimberKernel subclasses to the message log. This is
 * useful for debugging and discovering available kernel implementations.
 *
 */
void
nc_xcor_kernel_log_all_models (void)
{
  g_message ("# Registered NcXcorKernel:%s are:\n",
             g_type_name (NC_TYPE_XCOR_KERNEL));
  _nc_xcor_kernel_log_all_models_go (NC_TYPE_XCOR_KERNEL, 0);
}

