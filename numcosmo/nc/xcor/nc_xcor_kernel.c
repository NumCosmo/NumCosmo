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
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/core/ncm_memory_pool.h"
#include "nc/xcor/nc_xcor_kernel_component.h"
#include "nc/xcor/nc_xcor.h"
#include "nc_enum_types.h"

/* #define DEBUG */

typedef struct _NcXcorKernelPrivate
{
  /*< private >*/
  NcmModel parent_instance;
  NcDistance *dist;
  NcmPowspec *ps;
  NcmSBesselIntegrator *sbi;
  guint lmax;
  gint l_limber;
  gdouble adaptive_epsilon;
  guint adaptive_boundary_tries;
  gdouble reltol;
  gdouble scaled_abstol;
  guint max_border_expansions;
  guint max_iter;
  gdouble expansion_factor;
  guint panel_order_cap;
  gboolean track_fit_residual;
  gboolean tolerance_balance_warned;
  gboolean constructed;
} NcXcorKernelPrivate;

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
  PROP_SCALED_ABSTOL,
  PROP_MAX_BORDER_EXPANSIONS,
  PROP_MAX_ITER,
  PROP_EXPANSION_FACTOR,
  PROP_TRACK_FIT_RESIDUAL,
  PROP_PANEL_ORDER_CAP,
  PROP_SIZE,
};


G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcXcorKernel, nc_xcor_kernel, NCM_TYPE_MODEL)
G_DEFINE_BOXED_TYPE (NcXcorKinetic, nc_xcor_kinetic, nc_xcor_kinetic_copy, nc_xcor_kinetic_free)
G_DEFINE_BOXED_TYPE (NcXcorKernelIntegrand, nc_xcor_kernel_integrand, nc_xcor_kernel_integrand_ref, nc_xcor_kernel_integrand_unref)

static void
nc_xcor_kernel_init (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->dist                     = NULL;
  self->ps                       = NULL;
  self->sbi                      = NULL;
  self->lmax                     = 0;
  self->l_limber                 = 0;
  self->adaptive_epsilon         = 0.0;
  self->adaptive_boundary_tries  = 0;
  self->reltol                   = 0.0;
  self->scaled_abstol            = 0.0;
  self->max_border_expansions    = 0;
  self->max_iter                 = 0;
  self->expansion_factor         = 0.0;
  self->panel_order_cap          = 0;
  self->track_fit_residual       = FALSE;
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
    case PROP_SCALED_ABSTOL:
      nc_xcor_kernel_set_scaled_abstol (xclk, g_value_get_double (value));
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
    case PROP_TRACK_FIT_RESIDUAL:
      nc_xcor_kernel_set_track_fit_residual (xclk, g_value_get_boolean (value));
      break;
    case PROP_PANEL_ORDER_CAP:
      nc_xcor_kernel_set_panel_order_cap (xclk, g_value_get_uint (value));
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
    case PROP_SCALED_ABSTOL:
      g_value_set_double (value, nc_xcor_kernel_get_scaled_abstol (xclk));
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
    case PROP_TRACK_FIT_RESIDUAL:
      g_value_set_boolean (value, nc_xcor_kernel_get_track_fit_residual (xclk));
      break;
    case PROP_PANEL_ORDER_CAP:
      g_value_set_uint (value, nc_xcor_kernel_get_panel_order_cap (xclk));
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

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_POWSPEC,
                                   g_param_spec_object ("powspec",
                                                        NULL,
                                                        "Power spectrum object",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_INTEGRATOR,
                                   g_param_spec_object ("integrator",
                                                        NULL,
                                                        "Spherical Bessel integrator object",
                                                        NCM_TYPE_SBESSEL_INTEGRATOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "Maximum multipole",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

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
   * Convergence threshold for the adaptive $k$-range determination: sampling stops
   * extending outwards once $\Vert W\Vert$ falls below this fraction of its running
   * maximum, for #NcXcorKernel:adaptive-boundary-tries consecutive steps.
   *
   * The criterion is based on the amplitude of $W$ relative to its peak. For a tail
   * decaying as $k^{-p}$, the resulting range scales approximately as
   * $(1/\epsilon)^{1/p}$ times the peak location. Consequently, tightening $\epsilon$
   * can substantially increase the sampled range without a corresponding increase in
   * the contribution to the $C_\ell$ integral.
   *
   * This is particularly relevant for kernels with compact support in $z$, such as
   * top-hat redshift bins, whose transforms exhibit slowly decaying oscillatory tails.
   * Smooth kernels (galaxy, weak lensing, CMB lensing, tSZ) have amplitude and
   * integral contributions that fall off together, making the criterion a more direct
   * indication of the relevant integration range.
   */
  g_object_class_install_property (object_class,
                                   PROP_ADAPTIVE_EPSILON,
                                   g_param_spec_double ("adaptive-epsilon",
                                                        NULL,
                                                        "Convergence threshold for adaptive k-range determination",
                                                        0.0, 1.0, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

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
   * Relative tolerance for the adaptive refinement of the $W_i(k)$ spline.
   *
   * The tolerance is applied only where #NcXcorKernel:scaled-abstol does not bind. For
   * typical settings, the scaled absolute tolerance controls most of the refinement,
   * so this parameter primarily affects regions where the spline amplitude is
   * sufficiently large for relative accuracy to be relevant. See
   * #NcXcorKernel:scaled-abstol for the complementary criterion.
   *
   * ## The two tolerances gate each other -- move them together
   *
   * Refinement accepts an interval when
   * $\Vert f - \tilde f \Vert_2 \le \mathrm{reltol} \Vert f \Vert_2 + a \Vert f
   * \Vert_2^\mathrm{max}$, a **sum**. Whichever term is larger decides where
   * refinement stops, so tightening the other one alone changes nothing at all.
   * Measured on a Gaussian kernel, accuracy gained over the 1e-4/1e-4 defaults:
   *
   * | | a = 1e-4 | a = 1e-5 | a = 1e-6 |
   * |---|---|---|---|
   * | reltol 1e-4 | 1x | 1.0x | 1.0x |
   * | reltol 1e-6 | 1.2x | 15x | 25x |
   * | reltol 1e-8 | 1.2x | 21x | 280x |
   *
   * The first column and the first row are flat: one knob alone buys nothing,
   * whichever one it is. Tightening reltol is cheap -- a hundredfold costs about
   * 2% in knots -- but cheap and inert, and the knots are still paid for. A
   * kernel built with the two more than two orders apart says so on stderr.
   *
   * The defaults are equal for that reason: equal terms means neither is wasted.
   * They are a *cheap* balanced pair, not an accurate one. Moving both to
   * reltol 1e-6 with @a 1e-5 -- what #NcXcorSSCSij sets for itself -- is worth
   * 15-46x in accuracy for 1.6x the solve on smooth kernels and 3.0x on cluster
   * top-hats, whose fit is floor-driven and gains nothing from reltol at all.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance for adaptive midpoint refinement",
                                                        0.0, 1.0, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:scaled-abstol:
   *
   * Absolute tolerance for the adaptive refinement of the $W_i(k)$ spline, expressed
   * as a fraction of the peak sampled integrand. An interval is accepted once its
   * estimated interpolation error falls below `scaled-abstol` $\times \max\vert
   * F\vert$.
   *
   * This criterion sets the absolute accuracy of the spline where it binds. Tightening
   * #NcXcorKernel:reltol has no effect in those regions. Since the tolerance is
   * absolute, the corresponding relative error can become large where the resulting
   * $C_\ell$ is small; these contributions are typically dominated by cancellation and
   * have little impact on the total signal.
   *
   * Tightening this tolerance can substantially increase the number of spline knots,
   * particularly for kernels with slowly decaying oscillatory tails (see
   * #NcXcorKernel:adaptive-epsilon). The resulting increase in resolution can also
   * make the outer $k$ integral more difficult to evaluate, since the spline may
   * resolve oscillations that contribute negligibly to the integral.
   *
   * The useful precision is ultimately limited by the radial integration used to
   * compute $W_i(k)$: far below its peak that integral is dominated by cancellation,
   * so refining the spline beyond that level does not add reliable information.
   *
   * **Do not set this below $10^{-6}$.** The floor is measured against the peak of
   * $W_i(k)$, but the quantity actually integrated is $k^2 W_i W_j$, so the floor
   * enters *squared*: the $10^{-4}$ default is already $10^{-8}$ on the integrand,
   * which is about what the outer $k$ integral carries, and $10^{-6}$ here is
   * $10^{-12}$ there. Asking for $10^{-8}$ means $10^{-16}$ on the integrand, below
   * double precision, so it cannot improve the answer and measurably does not --
   * while costing up to two orders of magnitude in knots on a compactly supported
   * window. What a tighter floor would sharpen is the tail-times-tail part of the
   * product; a spectrum dominated by that is one whose two kernels barely overlap,
   * where the signal is negligible to begin with. Values below $10^{-6}$ emit a
   * warning from nc_xcor_kernel_set_scaled_abstol().
   */
  g_object_class_install_property (object_class,
                                   PROP_SCALED_ABSTOL,
                                   g_param_spec_double ("scaled-abstol",
                                                        NULL,
                                                        "Absolute tolerance scaled by the maximum kernel value for adaptive midpoint refinement",
                                                        GSL_DBL_MIN, 1.0, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MAX_BORDER_EXPANSIONS,
                                   g_param_spec_uint ("max-border-expansions",
                                                      NULL,
                                                      "Maximum number of border expansion iterations",
                                                      1, G_MAXUINT, 500,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MAX_ITER,
                                   g_param_spec_uint ("max-iter",
                                                      NULL,
                                                      "Maximum number of adaptive midpoint refinement iterations",
                                                      1, G_MAXUINT, 10000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:track-fit-residual:
   *
   * Whether the closure records the residual its fit actually achieved on each
   * knot interval, which is what nc_xcor_compute_full() turns into an error
   * estimate. On by default: without it the estimate has only
   * #NcXcorKernel:reltol and #NcXcorKernel:scaled-abstol to work from -- the
   * tolerances the fit was asked for, which it beats by 12 to 3100 times
   * depending on the kernel, so the resulting bound tracks the pair's
   * cancellation rather than its accuracy.
   *
   * The record costs one double per knot per multipole in the block, about
   * what the closure's own spline data costs, and #NcXcorSolver holds one
   * closure per kernel per $\ell$ block. Turn it off to get that memory back
   * from a run that never asks for an error.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_TRACK_FIT_RESIDUAL,
                                   g_param_spec_boolean ("track-fit-residual",
                                                         NULL,
                                                         "Whether to record the residual the closure fit achieved",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernel:panel-order-cap:
   *
   * Highest Chebyshev order tried on one panel before it is bisected, as
   * $N = 2^\\mathrm{cap} + 1$ coefficients.
   *
   * **This is a heuristic, and the default is fitted rather than derived.** The
   * cap trades waste against panel count: a panel that fails to converge at it
   * discards its whole grid before splitting, so a high cap wastes more per
   * failure while a low one fails more often. Neither side has a closed form,
   * so 5 comes from a sweep on two kernel families, solve time:
   *
   * | cap | $N \\le$ | galaxy + weak lensing | cluster top-hat |
   * |---|---|---|---|
   * | 5 | 33 | 1.97 s | 2.00 s |
   * | 6 | 65 | 2.64 s | 2.04 s |
   * | 7 | 129 | 3.73 s | 2.19 s |
   * | 8 | 257 | 5.74 s | 2.56 s |
   *
   * Uniformly best at 5 there, and accuracy did not move with it. But those are
   * two kernel families at one multipole range on one machine, and the optimum
   * depends on how a window's phase is distributed across its domain -- which
   * is a property of the kernel, not of the library. **A caller with a
   * different kernel should sweep this rather than assume 5 transfers**, and it
   * is a property rather than a compile-time constant so that they can.
   *
   * Zero selects the default.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PANEL_ORDER_CAP,
                                   g_param_spec_uint ("panel-order-cap",
                                                      NULL,
                                                      "Highest Chebyshev order tried per panel before bisecting",
                                                      0, 12, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

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
 * Spline-based integrand infrastructure (used by both Limber and non-Limber)
 */

/*
 * Spline-based integrand data (used by both Limber and non-Limber)
 */
typedef struct _SplineIntegrandData
{
  NcHICosmo *cosmo;
  gdouble RH_Mpc;
  gint lmin;
  guint len;
  NcmSplineVec *spline_vec;
  NcmVector *eval_result;
  gdouble k_min;
  gdouble k_max;
  gdouble *k_min_comp;
  gdouble *k_max_comp;
} SplineIntegrandData;

static NcmVector *
_spline_integrand_get_knots (gpointer data)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  /* Every component of the spline vec is sampled on the same abscissa, so the first
   * component's knots are the knots of the whole block. */
  return ncm_spline_peek_xv (ncm_spline_vec_peek_spline (sid->spline_vec, 0));
}

static void
_spline_integrand_eval (gpointer data, gdouble k, gdouble *W)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;
  guint i;

  ncm_spline_vec_eval (sid->spline_vec, k, sid->eval_result);

  for (i = 0; i < sid->len; i++)
  {
    W[i] = ncm_vector_get (sid->eval_result, i);
  }
}

/*
 * The components of a block share one abscissa, so one index lookup serves the
 * whole run -- the same saving ncm_spline_vec_eval() makes over the whole
 * vector, restricted to the components the caller is integrating.
 */
static void
_spline_integrand_eval_comps (gpointer data, gdouble k, guint offset, guint len, gdouble *W)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;
  NcmSpline *spline_0      = ncm_spline_vec_peek_spline (sid->spline_vec, 0);
  const gsize idx          = ncm_spline_get_index (spline_0, k);
  guint i;

  for (i = 0; i < len; i++)
  {
    NcmSpline *spline_i = ncm_spline_vec_peek_spline (sid->spline_vec, offset + i);

    W[offset + i] = ncm_spline_eval_idx (spline_i, k, idx);
  }
}

static void
_spline_integrand_get_range (gpointer data, gdouble *kmin, gdouble *kmax)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  *kmin = sid->k_min;
  *kmax = sid->k_max;
}

/*
 * Under the Limber approximation a multipole's window is supported only on its
 * own band in k, so the shared domain of an ell block carries one step per
 * multipole. Reporting the per-component support lets the outer integral put
 * each step on an integration limit, where it is not a discontinuity of the
 * integrand at all.
 */
static void
_spline_integrand_get_range_comp (gpointer data, guint i, gdouble *kmin, gdouble *kmax)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  *kmin = sid->k_min_comp[i];
  *kmax = sid->k_max_comp[i];
}

static void
_spline_integrand_data_free (gpointer data)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  nc_hicosmo_clear (&sid->cosmo);
  ncm_spline_vec_clear (&sid->spline_vec);
  ncm_vector_clear (&sid->eval_result);

  g_free (sid->k_min_comp);
  g_free (sid->k_max_comp);
  g_free (data);
}

/*
 * Chebyshev-based integrand data. Mirrors SplineIntegrandData -- same domain,
 * same per-component ranges, same sampled function -- and differs only in
 * carrying coefficients instead of a spline.
 */
static gpointer
_nc_xcor_spectral_alloc (gpointer userdata)
{
  return ncm_spectral_new ();
}

static void
_nc_xcor_spectral_free (gpointer p)
{
  ncm_spectral_free (NCM_SPECTRAL (p));
}

/*
 * A borrowed NcmSpectral. It carries sampling buffers, coefficient scratch and
 * cached FFTW plans, none of which can be shared: one kernel is evaluated
 * concurrently for different ell blocks, and one closure is restricted
 * concurrently for different pairs. Pooled rather than created per call so the
 * plans stay warm -- planning dominates a cold expansion.
 *
 * Return it with ncm_memory_pool_return().
 */
static NcmSpectral **
_nc_xcor_spectral_get (void)
{
  G_LOCK_DEFINE_STATIC (create_lock);

  static NcmMemoryPool *mp = NULL;

  G_LOCK (create_lock);

  if (mp == NULL)
    mp = ncm_memory_pool_new (_nc_xcor_spectral_alloc, NULL, _nc_xcor_spectral_free);

  G_UNLOCK (create_lock);

  return (NcmSpectral **) ncm_memory_pool_get (mp);
}

typedef struct _ChebPanel
{
  gdouble a;
  gdouble b;
  NcmMatrix *coeffs; /* len x N, one row per multipole */
  guint N;
} ChebPanel;

typedef struct _ChebIntegrandData
{
  NcHICosmo *cosmo;
  gdouble RH_Mpc;
  gint lmin;
  guint len;
  GArray *panels; /* ChebPanel, ascending and contiguous */
  GArray *edges;  /* gdouble, panels->len + 1 entries, for the lookup */
  gdouble k_min;
  gdouble k_max;
  gdouble *k_min_comp;
  gdouble *k_max_comp;
} ChebIntegrandData;

/*
 * Which panel holds @k. Panels are contiguous and ascending, so this is a
 * bisection over the edges -- the same lookup a spline does over its knots,
 * against far fewer of them.
 */
static const ChebPanel *
_cheb_integrand_find_panel (const ChebIntegrandData *cid, const gdouble k)
{
  const gdouble *edges = (const gdouble *) cid->edges->data;
  guint lo             = 0;
  guint hi             = cid->panels->len - 1;

  while (lo < hi)
  {
    const guint mid = (lo + hi + 1) / 2;

    if (k >= edges[mid])
      lo = mid;
    else
      hi = mid - 1;
  }

  return &g_array_index (cid->panels, ChebPanel, lo);
}

/*
 * Clenshaw, one component. O(N) against a spline's O(1), which is the price of
 * the representation for callers that evaluate pointwise; the exact method does
 * not evaluate at all, it works on the coefficients.
 */
static gdouble
_cheb_panel_eval_one (const ChebPanel *panel, guint comp, const gdouble t)
{
  const gdouble two_t = 2.0 * t;
  gdouble b_1         = 0.0;
  gdouble b_2         = 0.0;
  gint n;

  for (n = (gint) panel->N - 1; n >= 1; n--)
  {
    const gdouble b_0 = two_t * b_1 - b_2 + ncm_matrix_get (panel->coeffs, comp, n);

    b_2 = b_1;
    b_1 = b_0;
  }

  return t * b_1 - b_2 + ncm_matrix_get (panel->coeffs, comp, 0);
}

static void
_cheb_integrand_eval (gpointer data, gdouble k, gdouble *W)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;
  const ChebPanel *panel = _cheb_integrand_find_panel (cid, k);
  const gdouble t        = ncm_spectral_x_to_t (panel->a, panel->b, k);
  guint i;

  for (i = 0; i < cid->len; i++)
    W[i] = _cheb_panel_eval_one (panel, i, t);
}

static void
_cheb_integrand_eval_comps (gpointer data, gdouble k, guint offset, guint len, gdouble *W)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;
  const ChebPanel *panel = _cheb_integrand_find_panel (cid, k);
  const gdouble t        = ncm_spectral_x_to_t (panel->a, panel->b, k);
  guint i;

  for (i = 0; i < len; i++)
    W[offset + i] = _cheb_panel_eval_one (panel, offset + i, t);
}

static void
_cheb_integrand_get_range (gpointer data, gdouble *kmin, gdouble *kmax)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;

  *kmin = cid->k_min;
  *kmax = cid->k_max;
}

static void
_cheb_integrand_get_range_comp (gpointer data, guint i, gdouble *kmin, gdouble *kmax)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;

  *kmin = cid->k_min_comp[i];
  *kmax = cid->k_max_comp[i];
}

/*
 * Reports the expansion only when there is a single panel, since that is the
 * case a caller working on coefficients can use directly. With several panels
 * the closure is still exactly evaluable through eval(), which is what the
 * quadratures use.
 */
static gboolean
_cheb_integrand_get_spectral (gpointer data, NcmMatrix **coeffs, gdouble *k_min, gdouble *k_max)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;

  if (cid->panels->len != 1)
    return FALSE;

  {
    const ChebPanel *panel = &g_array_index (cid->panels, ChebPanel, 0);

    *coeffs = panel->coeffs;
    *k_min  = panel->a;
    *k_max  = panel->b;
  }

  return TRUE;
}

static guint
_cheb_integrand_get_panels (gpointer data)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;

  return cid->panels->len;
}

static void
_cheb_integrand_peek_panel (gpointer data, guint i, NcmMatrix **coeffs, gdouble *a, gdouble *b)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;
  const ChebPanel *panel = &g_array_index (cid->panels, ChebPanel, i);

  *coeffs = panel->coeffs;
  *a      = panel->a;
  *b      = panel->b;
}

/*
 * Coefficients on [a, b], which must lie inside one panel. A polynomial
 * restricted to a subinterval is still a polynomial of the same degree, so this
 * is an exact change of basis and not a refit -- and it costs O(N^2) in the
 * panel's coefficient count, bounded by #NcXcorKernel:panel-order-cap, against
 * N radial solves to sample the subinterval afresh.
 */
static gboolean
_cheb_integrand_restrict (gpointer data, gdouble a, gdouble b, NcmMatrix **coeffs)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;
  const ChebPanel *panel = _cheb_integrand_find_panel (cid, 0.5 * (a + b));
  GArray *row            = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), panel->N);
  NcmSpectral **spectral = _nc_xcor_spectral_get ();
  GArray *out            = NULL;
  guint c, i;

  if ((a < panel->a) || (b > panel->b))
  {
    g_array_unref (row);
    ncm_memory_pool_return (spectral);

    return FALSE;
  }

  ncm_matrix_clear (coeffs);
  *coeffs = ncm_matrix_new (cid->len, panel->N);

  g_array_set_size (row, panel->N);

  for (c = 0; c < cid->len; c++)
  {
    for (i = 0; i < panel->N; i++)
      g_array_index (row, gdouble, i) = ncm_matrix_get (panel->coeffs, c, i);

    ncm_spectral_chebyshev_rebase (*spectral, row, panel->N,
                                   panel->a, panel->b, a, b, &out);

    for (i = 0; i < panel->N; i++)
      ncm_matrix_set (*coeffs, c, i, g_array_index (out, gdouble, i));
  }

  g_array_unref (row);
  g_clear_pointer (&out, g_array_unref);
  ncm_memory_pool_return (spectral);

  return TRUE;
}

static void
_cheb_integrand_data_free (gpointer data)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;
  guint i;

  nc_hicosmo_clear (&cid->cosmo);

  for (i = 0; i < cid->panels->len; i++)
    ncm_matrix_clear (&g_array_index (cid->panels, ChebPanel, i).coeffs);

  g_array_unref (cid->panels);
  g_array_unref (cid->edges);

  g_free (cid->k_min_comp);
  g_free (cid->k_max_comp);
  g_free (data);
}

typedef struct _ComponentParams
{
  NcXcorKernelComponent *comp;
  NcHICosmo *cosmo;
} ComponentParams;

gdouble
_nc_xcor_kernel_component_kernel_integ (gpointer params, gdouble x, gdouble k)
{
  const ComponentParams *nlcp = (const ComponentParams *) params;
  const gdouble kernel        = nc_xcor_kernel_component_eval_kernel (nlcp->comp, nlcp->cosmo, x, k);

  return kernel / (k * sqrt (k));
}

#define MAX_ELL_BLOCK NC_XCOR_KERNEL_MAX_ELL_BLOCK
#define MAX_COMP_BLOCK 6

typedef struct _ComponentState
{
  NcXcorKernelComponent *comp;
  guint comp_idx;
  gdouble xi_min;
  gdouble xi_max;
  gdouble k_min_hard;
  gdouble k_max_hard;
  gdouble last_k_left;
  gdouble last_k_right;
  gdouble last_values_left[MAX_ELL_BLOCK];
  gdouble last_values_right[MAX_ELL_BLOCK];
  guint left_boundary_found;
  guint right_boundary_found;
  gdouble k_min_limber_ell[MAX_ELL_BLOCK]; /* Per-ell minimum k for Limber */
  gdouble k_max_limber_ell[MAX_ELL_BLOCK]; /* Per-ell maximum k for Limber */
  ComponentParams params;
} ComponentState;

typedef struct _ComponentStates
{
  ComponentState states[MAX_COMP_BLOCK];
  NcXcorKernel *xclk;
  NcmSBesselIntegrator *sbi; /* Integrator for this call; not owned */
  gdouble k_min_hard;
  gdouble k_max_hard;
  gdouble l2_norm;
  const guint n_comp;
  const guint lmin;
  const guint n_l;
  const gdouble epsilon;
  const guint adaptive_boundary_tries;
  const gboolean is_limber; /* k_min_limber_ell/k_max_limber_ell are set only then */
} ComponentStates;

static void
_component_state_init (ComponentState *state, NcXcorKernelComponent *comp, guint comp_idx, NcHICosmo *cosmo, guint n_l)
{
  state->comp                 = nc_xcor_kernel_component_ref (comp);
  state->comp_idx             = comp_idx;
  state->last_k_left          = G_MAXDOUBLE;
  state->last_k_right         = 0.0;
  state->left_boundary_found  = 0;
  state->right_boundary_found = 0;
  state->params.comp          = comp;
  state->params.cosmo         = cosmo;

  nc_xcor_kernel_component_get_limits (comp, cosmo,
                                       &state->xi_min, &state->xi_max,
                                       &state->k_min_hard, &state->k_max_hard);
}

gboolean
_is_new_k (gdouble k, gdouble *arr, guint n)
{
  guint m;

  for (m = 0; m < n; m++)
    if (gsl_fcmp (k, arr[m], 1e-3) == 0)
      return FALSE;

  return TRUE;
}

static GPtrArray *
_nc_xcor_kernel_validate_component_list (NcXcorKernel *xclk, guint n_l)
{
  NcXcorKernelClass *klass = NC_XCOR_KERNEL_GET_CLASS (xclk);
  GPtrArray *comp_list     = klass->get_component_list (xclk);

  if ((comp_list == NULL) || (comp_list->len == 0))
  {
    if (comp_list != NULL)
      g_ptr_array_unref (comp_list);

    g_error ("_nc_xcor_kernel_validate_component_list: kernel %s returned empty component list",
             G_OBJECT_TYPE_NAME (xclk));

    return NULL;
  }

  /* Hard errors rather than g_assert(): both bound writes into fixed-size
   * stack arrays (ComponentStates::last_values_*, kernel_out[][]), and asserts
   * compile out under -Dnumcosmo_assert=false, which would turn an
   * out-of-range block into a stack overflow instead of a clean abort. */
  if (n_l > MAX_ELL_BLOCK)
    g_error ("_nc_xcor_kernel_validate_component_list: kernel %s asked for %u multipoles "
             "in a single block, but at most %d fit (NC_XCOR_KERNEL_MAX_ELL_BLOCK). "
             "Split the range into blocks.",
             G_OBJECT_TYPE_NAME (xclk), n_l, MAX_ELL_BLOCK);

  if (comp_list->len > MAX_COMP_BLOCK)
    g_error ("_nc_xcor_kernel_validate_component_list: kernel %s has %u components, "
             "but at most %d fit in a single block.",
             G_OBJECT_TYPE_NAME (xclk), comp_list->len, MAX_COMP_BLOCK);

  return comp_list; /* Caller must unref */
}

static ComponentStates
_component_states_init_non_limber (NcXcorKernel *xclk, gint lmin, guint n_l,
                                   GPtrArray *comp_list, NcHICosmo *cosmo,
                                   NcmSBesselIntegrator *sbi)
{
  NcXcorKernelPrivate *self   = nc_xcor_kernel_get_instance_private (xclk);
  const guint n_comp          = comp_list->len;
  ComponentStates comp_states = {
    .xclk                    = xclk,
    .sbi                     = sbi,
    .k_min_hard              = 0.0,
    .k_max_hard              = G_MAXDOUBLE,
    .l2_norm                 = 0.0,
    .n_comp                  = n_comp,
    .lmin                    = lmin,
    .n_l                     = n_l,
    .epsilon                 = self->adaptive_epsilon,
    .adaptive_boundary_tries = self->adaptive_boundary_tries,
    .is_limber               = FALSE
  };
  guint i;

  /* Initialize each component state */
  for (i = 0; i < n_comp; i++)
  {
    ComponentState *state = &comp_states.states[i];

    _component_state_init (state, g_ptr_array_index (comp_list, i), i, cosmo, n_l);

    /* Compute global hard limits as intersection of all component limits */
    comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, state->k_min_hard);
    comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, state->k_max_hard);
  }

  g_assert_cmpfloat (comp_states.k_min_hard, <, comp_states.k_max_hard);

  return comp_states;
}

/*
 * Limber peak value of one component at (k, l): the stationary-phase
 * approximation of int K(chi, k) j_l^(d)(k chi) dchi, including the
 * component prefactor. For d > 0 the derivative is expanded in j_l and
 * j_{l+1} through the downward recurrences and each term is
 * peak-approximated at its own nu / k. The j_{l+1} peak sits at
 * (nu + 1) / k, slightly deeper than the j_l one; when it falls beyond the
 * component's support the window there is zero and the term is dropped.
 */
static gdouble
_component_limber_eval (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi_max, gdouble k, gint l)
{
  const guint deriv       = nc_xcor_kernel_component_get_bessel_deriv (comp);
  const gdouble nu        = l + 0.5;
  const gdouble prefactor = nc_xcor_kernel_component_eval_prefactor (comp, cosmo, k, l);
  const gdouble peak_l    = sqrt (M_PI / (2.0 * nu)) * nc_xcor_kernel_component_eval_kernel (comp, cosmo, nu / k, k);
  gdouble val;

  if (deriv == 0)
  {
    val = peak_l;
  }
  else
  {
    const gdouble nup      = nu + 1.0;
    const gdouble xi_p     = nup / k;
    const gdouble peak_lp1 = (xi_p <= xi_max) ?
                             sqrt (M_PI / (2.0 * nup)) * nc_xcor_kernel_component_eval_kernel (comp, cosmo, xi_p, k) :
                             0.0;

    if (deriv == 1) /* j_l' = (l/y) j_l - j_{l+1} */
      val = (l / nu) * peak_l - peak_lp1;
    else /* j_l'' = (l (l-1)/y^2 - 1) j_l + (2/y) j_{l+1} */
      val = -(2.0 * l + 0.25) / (nu * nu) * peak_l + (2.0 / nup) * peak_lp1;
  }

  return prefactor * val / k;
}

static ComponentStates
_component_states_init_limber (NcXcorKernel *xclk, gint lmin, guint n_l,
                               GPtrArray *comp_list, NcHICosmo *cosmo)
{
  NcXcorKernelPrivate *self   = nc_xcor_kernel_get_instance_private (xclk);
  const guint n_comp          = comp_list->len;
  ComponentStates comp_states = {
    .xclk                    = xclk,
    .k_min_hard              = 0.0,
    .k_max_hard              = G_MAXDOUBLE,
    .l2_norm                 = 0.0,
    .n_comp                  = n_comp,
    .lmin                    = lmin,
    .n_l                     = n_l,
    .epsilon                 = self->adaptive_epsilon,
    .adaptive_boundary_tries = self->adaptive_boundary_tries,
    .is_limber               = TRUE
  };
  gdouble k_min_union = G_MAXDOUBLE; /* For union of Limber bounds */
  gdouble k_max_union = 0.0;         /* For union of Limber bounds */
  guint i, j;

  /* Initialize each component state and compute hard limit intersection */
  for (i = 0; i < n_comp; i++)
  {
    ComponentState *state = &comp_states.states[i];

    _component_state_init (state, g_ptr_array_index (comp_list, i), i, cosmo, n_l);

    /* Compute intersection of component hard limits */
    comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, state->k_min_hard);
    comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, state->k_max_hard);

    /* Compute per-ell Limber constraints for this component */
    for (j = 0; j < n_l; j++)
    {
      const gint l_j     = lmin + j;
      const gdouble nu_j = l_j + 0.5;

      /* Limber constraint: xi = nu/k must be in [xi_min, xi_max]
       * Therefore: k must be in [nu/xi_max, nu/xi_min]
       */
      gdouble k_min_limber_j = nu_j / state->xi_max;
      gdouble k_max_limber_j = nu_j / state->xi_min;

      state->k_min_limber_ell[j] = k_min_limber_j;
      state->k_max_limber_ell[j] = k_max_limber_j;

      /* Update union of all (component, ell) Limber bounds */
      k_min_union = GSL_MIN (k_min_union, state->k_min_limber_ell[j]);
      k_max_union = GSL_MAX (k_max_union, state->k_max_limber_ell[j]);
    }
  }

  /* Override global bounds with union of Limber bounds (within hard limits) */
  comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, k_min_union);
  comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, k_max_union);

  g_assert_cmpfloat (comp_states.k_min_hard, <, comp_states.k_max_hard);

  /* Compute boundary values for extrapolation */
  for (i = 0; i < n_comp; i++)
  {
    ComponentState *state = &comp_states.states[i];

    for (j = 0; j < n_l; j++)
    {
      const gint l_j        = lmin + j;
      const gdouble k_min_j = state->k_min_limber_ell[j];
      const gdouble k_max_j = state->k_max_limber_ell[j];

      /* Evaluate at the two boundaries */
      state->last_values_left[j]  = _component_limber_eval (state->comp, cosmo, state->xi_max, k_min_j, l_j);
      state->last_values_right[j] = _component_limber_eval (state->comp, cosmo, state->xi_max, k_max_j, l_j);
    }
  }

  return comp_states;
}

static void
_component_states_compute_k_seeds (ComponentStates *comp_states, GArray *k_seeds)
{
  gdouble log_k_center_sum = 0.0;
  gdouble k_comp_scales[MAX_COMP_BLOCK];
  gdouble k_min_soft = G_MAXDOUBLE;
  gdouble k_max_soft = 0.0;
  gdouble k_center;
  guint i, j;

  /* Clear the array and prepare for new values */
  g_array_set_size (k_seeds, 0);

  /* Compute soft limits and component scales */
  for (i = 0; i < comp_states->n_comp; i++)
  {
    ComponentState *state = &comp_states->states[i];
    gdouble ln_k_scale    = 0.0;
    gdouble n_k           = 0.0;

    for (j = 0; j < comp_states->n_l; j++)
    {
      const gdouble nu       = comp_states->lmin + j + 0.5;
      const gdouble k_max_ij = nc_xcor_kernel_component_eval_k_max (state->comp, nu);
      gdouble k_upper_ij     = k_max_ij * 1.01;
      gdouble k_lower_ij     = k_max_ij * 0.99;

      if ((k_lower_ij > comp_states->k_min_hard) && (k_upper_ij < comp_states->k_max_hard))
      {
        ln_k_scale += log (k_max_ij);
        n_k        += 1.0;
        k_max_soft  = GSL_MAX (k_max_soft, k_upper_ij);
        k_min_soft  = GSL_MIN (k_min_soft, k_lower_ij);
      }
    }

    if (n_k > 0.0)
    {
      log_k_center_sum += ln_k_scale / n_k;
      k_comp_scales[i]  = exp (ln_k_scale / n_k);
    }
    else
    {
      /* If no valid k_max was found for this component, use the geometric mean of hard limits */
      k_comp_scales[i]  = sqrt (comp_states->k_min_hard * comp_states->k_max_hard);
      log_k_center_sum += log (k_comp_scales[i]);
      k_max_soft        = GSL_MAX (k_max_soft, k_comp_scales[i] * (1.0 + 1.0e-5));
      k_min_soft        = GSL_MIN (k_min_soft, k_comp_scales[i] * (1.0 - 1.0e-5));
    }
  }

  k_center = exp (log_k_center_sum / comp_states->n_comp);

  if ((k_center < k_min_soft) || (k_center > k_max_soft))
    k_center = (k_min_soft + k_max_soft) / 2.0;

  g_assert_cmpfloat (k_min_soft, <, k_max_soft);
  g_assert_cmpfloat (comp_states->k_min_hard, <=, k_min_soft);
  g_assert_cmpfloat (k_max_soft, <=, comp_states->k_max_hard);
  g_assert_cmpfloat (k_min_soft, <, k_center);
  g_assert_cmpfloat (k_center, <, k_max_soft);

  /* Add k_min_soft */
  g_array_append_val (k_seeds, k_min_soft);

  /* Add k_center if unique */
  if (_is_new_k (k_center, (gdouble *) k_seeds->data, k_seeds->len))
    g_array_append_val (k_seeds, k_center);

  /* Add unique component scales */
  for (i = 0; i < comp_states->n_comp; i++)
  {
    if (_is_new_k (k_comp_scales[i], (gdouble *) k_seeds->data, k_seeds->len))
      g_array_append_val (k_seeds, k_comp_scales[i]);
  }

  /* Add k_max_soft */
  if (_is_new_k (k_max_soft, (gdouble *) k_seeds->data, k_seeds->len))
    g_array_append_val (k_seeds, k_max_soft);
}

#define DECAY_RATE 1.0e10

static void
_component_states_compute_non_limber (const gdouble k, NcmVector *y, gpointer user_data)
{
  ComponentStates *comp_states = (ComponentStates *) user_data;
  gdouble kernel_out[MAX_COMP_BLOCK][MAX_ELL_BLOCK];
  gdouble l2_norm = 0.0;
  guint ci, i;

  /* Compute kernel for each component */
  for (ci = 0; ci < comp_states->n_comp; ci++)
  {
    NcmVector *integ_result              = ncm_vector_new_data_static (kernel_out[ci], comp_states->n_l, 1);
    ComponentState *state                = &comp_states->states[ci];
    const gboolean right_boundary_found  = state->right_boundary_found >= comp_states->adaptive_boundary_tries;
    const gboolean left_boundary_found   = state->left_boundary_found >= comp_states->adaptive_boundary_tries;
    const gboolean within_left_boundary  = !left_boundary_found || (k >= state->last_k_left);
    const gboolean within_right_boundary = !right_boundary_found || (k <= state->last_k_right);

    if (within_left_boundary && within_right_boundary)
    {
      gdouble component_l2_norm = 0.0;
      gboolean below_epsilon    = FALSE;

      /* Exact integration within boundaries */
      ncm_sbessel_integrator_integrate_deriv (
        comp_states->sbi, _nc_xcor_kernel_component_kernel_integ,
        state->xi_min, state->xi_max, k,
        nc_xcor_kernel_component_get_bessel_deriv (state->comp),
        integ_result, &state->params
      );

      for (i = 0; i < comp_states->n_l; i++)
      {
        const gdouble prefactor = nc_xcor_kernel_component_eval_prefactor (
          state->comp, state->params.cosmo, k, comp_states->lmin + i
        );

        kernel_out[ci][i] *= prefactor * k * sqrt (k);
        component_l2_norm += kernel_out[ci][i] * kernel_out[ci][i];
      }

      below_epsilon = component_l2_norm < gsl_pow_2 (comp_states->epsilon * comp_states->l2_norm);

      /* Update boundary tracking */
      if (k > state->last_k_right)
      {
        state->last_k_right = k;

        for (i = 0; i < comp_states->n_l; i++)
          state->last_values_right[i] = kernel_out[ci][i];

        if (below_epsilon)
          state->right_boundary_found++;
        else
          state->right_boundary_found = 0;
      }
      else if (k < state->last_k_left)
      {
        state->last_k_left = k;

        for (i = 0; i < comp_states->n_l; i++)
          state->last_values_left[i] = kernel_out[ci][i];

        if (below_epsilon)
          state->left_boundary_found++;
        else
          state->left_boundary_found = 0;
      }
    }
    else
    {
      /* Exponential tail extrapolation beyond boundaries */
      if (!within_right_boundary)
      {
        for (i = 0; i < comp_states->n_l; i++)
        {
          const gdouble val              = state->last_values_right[i];
          const gdouble delta_k          = k - state->last_k_right;
          const gdouble decay_rate       = DECAY_RATE;
          const gdouble val_extrapolated = val * exp (-decay_rate * delta_k / state->last_k_right);

          kernel_out[ci][i] = val_extrapolated;
        }
      }
      else
      {
        for (i = 0; i < comp_states->n_l; i++)
        {
          const gdouble val              = state->last_values_left[i];
          const gdouble delta_k          = state->last_k_left - k;
          const gdouble decay_rate       = DECAY_RATE;
          const gdouble val_extrapolated = val * exp (-decay_rate * delta_k / state->last_k_left);

          kernel_out[ci][i] = val_extrapolated;
        }
      }
    }

    ncm_vector_clear (&integ_result);
  }

  /* Sum contributions from all components and compute total L2 norm */
  g_assert_cmpuint (ncm_vector_len (y), ==, comp_states->n_l);

  for (i = 0; i < comp_states->n_l; i++)
  {
    gdouble sum = 0.0;

    for (ci = 0; ci < comp_states->n_comp; ci++)
      sum += kernel_out[ci][i];

    ncm_vector_set (y, i, sum);
    l2_norm += sum * sum;
  }

  l2_norm = sqrt (l2_norm);

  /* Update reference L2 norm for convergence testing */
  if (l2_norm > comp_states->l2_norm)
    comp_states->l2_norm = l2_norm;
}

static void
_component_states_compute_limber (const gdouble k, NcmVector *y, gpointer user_data)
{
  ComponentStates *comp_states = (ComponentStates *) user_data;
  gdouble kernel_out[MAX_COMP_BLOCK][MAX_ELL_BLOCK];
  gdouble l2_norm = 0.0;
  guint ci, i;

  /* Compute kernel for each component using Limber approximation */
  for (ci = 0; ci < comp_states->n_comp; ci++)
  {
    ComponentState *state = &comp_states->states[ci];

    /* Evaluate with per-ell boundary checking and extrapolation */
    for (i = 0; i < comp_states->n_l; i++)
    {
      const gint l                = comp_states->lmin + i;
      const gboolean within_range = (k >= state->k_min_limber_ell[i]) && (k <= state->k_max_limber_ell[i]);

      if (within_range)
      {
        /* Normal Limber evaluation within valid range */
        kernel_out[ci][i] = _component_limber_eval (state->comp, state->params.cosmo, state->xi_max, k, l);
      }
      else
      {
        /* Exponential extrapolation outside valid Limber range */
        if (k < state->k_min_limber_ell[i])
        {
          /* Left extrapolation */
          const gdouble delta_k    = state->k_min_limber_ell[i] - k;
          const gdouble decay_rate = DECAY_RATE;
          const gdouble decay      = exp (-gsl_pow_2 (decay_rate * delta_k / state->k_min_limber_ell[i]));

          kernel_out[ci][i] = state->last_values_left[i] * decay;
        }
        else /* k > state->k_max_limber_ell[i] */
        {
          /* Right extrapolation */
          const gdouble delta_k    = k - state->k_max_limber_ell[i];
          const gdouble decay_rate = DECAY_RATE;
          const gdouble decay      = exp (-gsl_pow_2 (decay_rate * delta_k / state->k_max_limber_ell[i]));

          kernel_out[ci][i] = state->last_values_right[i] * decay;
        }
      }
    }
  }

  /* Sum contributions from all components and compute total L2 norm */
  for (i = 0; i < comp_states->n_l; i++)
  {
    gdouble sum = 0.0;

    for (ci = 0; ci < comp_states->n_comp; ci++)
      sum += kernel_out[ci][i];

    ncm_vector_set (y, i, sum);
    l2_norm += sum * sum;
  }

  l2_norm = sqrt (l2_norm);

  /* Update reference L2 norm for overall convergence testing */
  if (l2_norm > comp_states->l2_norm)
    comp_states->l2_norm = l2_norm;
}

/*
 * The refinement criterion is reltol * ||f||_2 + a * ||f||_2^max, a *sum*, so
 * the larger of the two terms sets what refinement stops at and the smaller one
 * is inert. Setting one far tighter than the other therefore buys nothing while
 * still being paid for in knots, and nothing else reports it.
 *
 * The two terms are not compared directly -- one is scaled to the block's norm
 * and the other to its peak, a factor of order sqrt(n_l) apart -- so the
 * threshold is deliberately loose at two orders, firing only where the
 * imbalance cannot be anything else. Once per kernel: closures are built per
 * ell block, and the tolerances do not change between them.
 */
static void
_nc_xcor_kernel_check_tolerance_balance (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);
  const gdouble ratio       = self->reltol / self->scaled_abstol;

  if (self->tolerance_balance_warned)
    return;

  if ((ratio > 1.0e2) || (ratio < 1.0e-2))
  {
    const gboolean reltol_inert = (ratio < 1.0e-2);

    self->tolerance_balance_warned = TRUE;

    g_warning ("_nc_xcor_kernel_check_tolerance_balance: %s has reltol %.3e and "
               "scaled-abstol %.3e, %.0f orders apart. The refinement criterion adds "
               "the two, so the looser one decides where refinement stops and %s is "
               "inert -- tightening it alone cannot improve the result, and it is "
               "still paid for in spline knots. Move them together.",
               G_OBJECT_TYPE_NAME (xclk), self->reltol, self->scaled_abstol,
               fabs (log10 (ratio)), reltol_inert ? "reltol" : "scaled-abstol");
  }
}

/*
 * NcmSpectralFBatch takes user data first; the samplers here take x first, as
 * NcmFunctionSampleSetFunc does. One adapter rather than changing either.
 */
typedef struct _ChebSampler
{
  void (*compute_func) (const gdouble, NcmVector *, gpointer);

  gpointer comp_states;
} ChebSampler;

static void
_cheb_sampler_call (gpointer user_data, gdouble k, NcmVector *y)
{
  ChebSampler *sampler = (ChebSampler *) user_data;

  sampler->compute_func (k, y, sampler->comp_states);
}

/* Default panel-order cap for Chebyshev closure fits. */
#define NC_XCOR_KERNEL_CHEB_PANEL_K_CAP (5)
#define NC_XCOR_KERNEL_CHEB_MIN_PANEL_FRAC (1.0e-6)

/*
 * Expands on [a, b], bisecting where the capped order does not converge.
 * Panels are appended in ascending order, so the result is contiguous.
 */
static void
_nc_xcor_kernel_cheb_split (NcmSpectral *spectral, gpointer sampler, guint n_l,
                            gdouble a, gdouble b, gdouble reltol, gdouble abstol,
                            guint k_cap, GArray *panels)
{
  NcmMatrix *coeffs = NULL;
  const guint k_ord = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (
    spectral, _cheb_sampler_call, n_l, a, b, 3,
    k_cap, reltol, abstol, FALSE, &coeffs, sampler);

  if (k_ord > 0)
  {
    ChebPanel panel = { a, b, coeffs, (1u << k_ord) + 1u };

    g_array_append_val (panels, panel);

    return;
  }

  {
    const gdouble mid = 0.5 * (a + b);

    /* A panel that will not converge however far it is split is not a
     * resolution problem; refusing to bisect past a fraction of the domain
     * turns a hang into a diagnosable expansion. */
    if ((b - a) < NC_XCOR_KERNEL_CHEB_MIN_PANEL_FRAC * b)
    {
      const guint k_forced = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (
        spectral, _cheb_sampler_call, n_l, a, b, 3,
        k_cap, reltol, abstol, TRUE, &coeffs, sampler);
      ChebPanel panel = { a, b, coeffs, (1u << k_forced) + 1u };

      g_array_append_val (panels, panel);

      return;
    }

    _nc_xcor_kernel_cheb_split (spectral, sampler, n_l, a, mid, reltol, abstol, k_cap, panels);
    _nc_xcor_kernel_cheb_split (spectral, sampler, n_l, mid, b, reltol, abstol, k_cap, panels);
  }
}

/*
 * Builds the closure as a Chebyshev series rather than a refined spline.
 *
 * The domain is found exactly as the spline path finds it -- the seeds and
 * ncm_function_sample_set_expand_domain() are shared, since where W is
 * negligible is a property of the kernel and not of the representation. What
 * changes is everything after: instead of bisecting until a fit criterion is
 * met, the whole ell block is expanded on one Chebyshev-Lobatto grid, doubling
 * the order until every multipole's coefficients converge.
 *
 * The samples that the domain expansion took are discarded, which the spline
 * path also does -- it keeps their abscissa but re-tests everything.
 */
static NcXcorKernelIntegrand *
_nc_xcor_kernel_build_cheb_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax,
                                      ComponentStates *comp_states,
                                      void (*compute_func) (const gdouble, NcmVector *, gpointer),
                                      const gdouble reltol, const gdouble abs_reltol)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);
  ChebIntegrandData *cid    = g_new0 (ChebIntegrandData, 1);
  const guint n_l           = lmax - lmin + 1;

  {
    NcmFunctionSampleSet *fss = ncm_function_sample_set_new (n_l);
    GArray *k_seeds           = g_array_new (FALSE, FALSE, sizeof (gdouble));
    ChebSampler sampler       = { compute_func, comp_states };
    gdouble abstol;
    guint i;

    _component_states_compute_k_seeds (comp_states, k_seeds);

    for (i = 0; i < k_seeds->len; i++)
    {
      const gdouble k_seed = g_array_index (k_seeds, gdouble, i);

      ncm_function_sample_set_add_old_func (fss, k_seed, compute_func, comp_states);
    }

    g_array_unref (k_seeds);

    ncm_function_sample_set_expand_domain (
      fss,
      compute_func,
      comp_states->k_min_hard,
      comp_states->k_max_hard,
      self->expansion_factor,
      comp_states->epsilon,
      self->max_border_expansions,
      comp_states->adaptive_boundary_tries,
      comp_states
    );

    cid->k_min = ncm_function_sample_set_get_x_min (fss);
    cid->k_max = ncm_function_sample_set_get_x_max (fss);

    /* Same meaning the spline path gives it: a floor scaled to the smallest of
     * the block's peaks, so a sub-dominant multipole is not held to a
     * tolerance relative to its neighbours. */
    abstol = ncm_function_sample_set_get_absmaxF_min (fss) * abs_reltol;

    ncm_function_sample_set_clear (&fss);

    cid->panels = g_array_new (FALSE, FALSE, sizeof (ChebPanel));
    cid->edges  = g_array_new (FALSE, FALSE, sizeof (gdouble));

    {
      NcmSpectral **spectral = _nc_xcor_spectral_get ();

      _nc_xcor_kernel_cheb_split (*spectral, &sampler, n_l,
                                  cid->k_min, cid->k_max, reltol, abstol,
                                  (self->panel_order_cap == 0) ?
                                  NC_XCOR_KERNEL_CHEB_PANEL_K_CAP : self->panel_order_cap,
                                  cid->panels);

      ncm_memory_pool_return (spectral);
    }

    {
      const gdouble first = g_array_index (cid->panels, ChebPanel, 0).a;

      g_array_append_val (cid->edges, first);

      for (i = 0; i < cid->panels->len; i++)
        g_array_append_val (cid->edges, g_array_index (cid->panels, ChebPanel, i).b);
    }

    cid->k_min_comp = g_new (gdouble, n_l);
    cid->k_max_comp = g_new (gdouble, n_l);

    for (i = 0; i < n_l; i++)
    {
      gdouble k_min_i = cid->k_min;
      gdouble k_max_i = cid->k_max;

      if (comp_states->is_limber)
      {
        gdouble band_min = G_MAXDOUBLE;
        gdouble band_max = 0.0;
        guint ci;

        for (ci = 0; ci < comp_states->n_comp; ci++)
        {
          band_min = GSL_MIN (band_min, comp_states->states[ci].k_min_limber_ell[i]);
          band_max = GSL_MAX (band_max, comp_states->states[ci].k_max_limber_ell[i]);
        }

        k_min_i = GSL_MAX (k_min_i, band_min);
        k_max_i = GSL_MIN (k_max_i, band_max);

        if (k_min_i >= k_max_i)
          k_min_i = k_max_i = cid->k_min;
      }

      cid->k_min_comp[i] = k_min_i;
      cid->k_max_comp[i] = k_max_i;
    }

    cid->lmin   = lmin;
    cid->len    = n_l;
    cid->RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);
    cid->cosmo  = nc_hicosmo_ref (cosmo);

    {
      NcXcorKernelIntegrand *integrand = nc_xcor_kernel_integrand_new (n_l,
                                                                       _cheb_integrand_eval,
                                                                       _cheb_integrand_get_range,
                                                                       cid,
                                                                       _cheb_integrand_data_free);

      nc_xcor_kernel_integrand_set_get_range_comp (integrand, _cheb_integrand_get_range_comp);
      nc_xcor_kernel_integrand_set_eval_comps (integrand, _cheb_integrand_eval_comps);
      nc_xcor_kernel_integrand_set_get_spectral (integrand, _cheb_integrand_get_spectral);
      nc_xcor_kernel_integrand_set_panel_accessors (integrand,
                                                    _cheb_integrand_get_panels,
                                                    _cheb_integrand_peek_panel);
      nc_xcor_kernel_integrand_set_restrict (integrand, _cheb_integrand_restrict);
      nc_xcor_kernel_integrand_set_tolerances (integrand, reltol, abs_reltol);

      return integrand;
    }
  }
}

static NcXcorKernelIntegrand *
_nc_xcor_kernel_build_spline_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax,
                                        ComponentStates *comp_states,
                                        void (*compute_func) (const gdouble, NcmVector *, gpointer),
                                        const gdouble reltol, const gdouble abs_reltol)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);
  SplineIntegrandData *sid  = g_new0 (SplineIntegrandData, 1);
  const guint n_l           = lmax - lmin + 1;

  {
    NcmFunctionSampleSet *fss = ncm_function_sample_set_new (n_l);
    NcmSpline *spline         = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    GArray *k_seeds           = g_array_new (FALSE, FALSE, sizeof (gdouble));
    NcmMatrix *residuals      = NULL;
    guint i;

    _nc_xcor_kernel_check_tolerance_balance (xclk);
    ncm_function_sample_set_set_track_residual (fss, self->track_fit_residual);

    /* Compute k-seeds for initial sampling. Local to the call, not kernel
     * state: one kernel may be evaluated concurrently for different ell
     * blocks, each with its own integrator. */
    _component_states_compute_k_seeds (comp_states, k_seeds);

    /* Add all k-seeds as initial sampling points */
    for (i = 0; i < k_seeds->len; i++)
    {
      const gdouble k_seed = g_array_index (k_seeds, gdouble, i);

      ncm_function_sample_set_add_old_func (fss, k_seed, compute_func, comp_states);
    }

    g_array_unref (k_seeds);

    {
      /* Domain expansion */
      ncm_function_sample_set_expand_domain (
        fss,
        compute_func,
        comp_states->k_min_hard,
        comp_states->k_max_hard,
        self->expansion_factor,
        comp_states->epsilon,
        self->max_border_expansions,
        comp_states->adaptive_boundary_tries,
        comp_states
      );
    }
    ncm_function_sample_set_mark_all_old (fss);
    ncm_function_sample_set_reset_interval_ok (fss);
    {
      const gdouble max_absF_total = ncm_function_sample_set_get_absmaxF_min (fss);

      ncm_function_sample_set_adaptive_midpoint (
        fss, compute_func,
        reltol, max_absF_total * abs_reltol, self->max_iter, 1,
        spline, comp_states
      );
    }

    residuals       = ncm_function_sample_set_get_residuals (fss);
    sid->spline_vec = ncm_function_sample_set_to_spline_vec (fss, spline);
    sid->k_min      = ncm_function_sample_set_get_x_min (fss);
    sid->k_max      = ncm_function_sample_set_get_x_max (fss);
    sid->k_min_comp = g_new (gdouble, n_l);
    sid->k_max_comp = g_new (gdouble, n_l);

    /* Per-multipole support within the block's shared domain. Only the Limber
     * branch confines a multipole to a band of its own; outside it the window
     * is zero, so the band edge falling inside the shared domain is a step.
     * The band is taken over all components, since the window is their sum. */
    for (i = 0; i < n_l; i++)
    {
      gdouble k_min_i = sid->k_min;
      gdouble k_max_i = sid->k_max;

      if (comp_states->is_limber)
      {
        gdouble band_min = G_MAXDOUBLE;
        gdouble band_max = 0.0;
        guint ci;

        for (ci = 0; ci < comp_states->n_comp; ci++)
        {
          band_min = GSL_MIN (band_min, comp_states->states[ci].k_min_limber_ell[i]);
          band_max = GSL_MAX (band_max, comp_states->states[ci].k_max_limber_ell[i]);
        }

        k_min_i = GSL_MAX (k_min_i, band_min);
        k_max_i = GSL_MIN (k_max_i, band_max);

        /* A band disjoint from the fitted domain leaves nothing to integrate;
         * report the empty range as the domain's lower edge rather than an
         * inverted one. */
        if (k_min_i >= k_max_i)
          k_min_i = k_max_i = sid->k_min;
      }

      sid->k_min_comp[i] = k_min_i;
      sid->k_max_comp[i] = k_max_i;
    }

    sid->lmin        = lmin;
    sid->len         = n_l;
    sid->RH_Mpc      = nc_hicosmo_RH_Mpc (cosmo);
    sid->cosmo       = nc_hicosmo_ref (cosmo);
    sid->eval_result = ncm_vector_new (n_l);

    ncm_function_sample_set_clear (&fss);
    ncm_spline_free (spline);

    {
      NcXcorKernelIntegrand *integrand = nc_xcor_kernel_integrand_new (n_l,
                                                                       _spline_integrand_eval,
                                                                       _spline_integrand_get_range,
                                                                       sid,
                                                                       _spline_integrand_data_free);

      nc_xcor_kernel_integrand_set_get_knots (integrand, _spline_integrand_get_knots);
      nc_xcor_kernel_integrand_set_get_range_comp (integrand, _spline_integrand_get_range_comp);
      nc_xcor_kernel_integrand_set_eval_comps (integrand, _spline_integrand_eval_comps);

      nc_xcor_kernel_integrand_set_tolerances (integrand, reltol, abs_reltol);
      nc_xcor_kernel_integrand_set_residuals (integrand, residuals);
      ncm_matrix_clear (&residuals);

      return integrand;
    }
  }
}

static NcXcorKernelIntegrand *
_nc_xcor_kernel_build_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);
  const guint n_l           = lmax - lmin + 1;
  GPtrArray *comp_list      = _nc_xcor_kernel_validate_component_list (xclk, n_l);

  if (comp_list == NULL)
    return NULL;

  /* Initialize with Limber-specific limits */
  {
    ComponentStates comp_states = _component_states_init_limber (xclk, lmin, n_l, comp_list, cosmo);

    g_ptr_array_unref (comp_list);

    /* Always the spline here, whatever #NcXcor:closure-type asks for.
     * Under Limber a multipole's window is supported only on its own band in k
     * and is zero outside it, so the block's shared domain carries one step per
     * multipole -- see _spline_integrand_get_range_comp(). A Chebyshev series
     * converges on this kernel because W_l(k) is entire in k, and a step is
     * not: the expansion would never converge and the panel splitter would
     * bisect until it gave up. The Limber closure is also the cheap one, so
     * there is nothing to win by trying. */
    return _nc_xcor_kernel_build_spline_integrand (xclk, cosmo, lmin, lmax,
                                                   &comp_states,
                                                   _component_states_compute_limber,
                                                   self->reltol, self->scaled_abstol);
  }
}

static NcXcorKernelIntegrand *
_nc_xcor_kernel_build_non_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi, NcXcorKernelClosure closure_type)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);
  const guint n_l           = lmax - lmin + 1;
  GPtrArray *comp_list      = _nc_xcor_kernel_validate_component_list (xclk, n_l);

  if (comp_list == NULL)
    return NULL;

  if (sbi == NULL)
  {
    g_ptr_array_unref (comp_list);
    g_error ("_nc_xcor_kernel_build_non_limber_integrand: no integrator for kernel %s. "
             "Either set the 'integrator' property or pass one to "
             "nc_xcor_kernel_get_eval_vectorized_full().",
             G_OBJECT_TYPE_NAME (xclk));

    return NULL;
  }

  ncm_sbessel_integrator_set_ell_range (sbi, lmin, lmax);

  /* Initialize with standard (non-Limber) limits */
  {
    ComponentStates comp_states = _component_states_init_non_limber (xclk, lmin, n_l, comp_list, cosmo, sbi);

    g_ptr_array_unref (comp_list);

    if (closure_type == NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV)
      return _nc_xcor_kernel_build_cheb_integrand (xclk, cosmo, lmin, lmax,
                                                   &comp_states,
                                                   _component_states_compute_non_limber,
                                                   self->reltol, self->scaled_abstol);

    return _nc_xcor_kernel_build_spline_integrand (xclk, cosmo, lmin, lmax,
                                                   &comp_states,
                                                   _component_states_compute_non_limber,
                                                   self->reltol, self->scaled_abstol);
  }
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
 * nc_xcor_kernel_integrand_new:
 * @len: number of components in the integrand
 * @eval: (scope async): function to evaluate the integrand
 * @get_range: (scope async): function to get the k range
 * @data: (nullable): user data to pass to @eval and @get_range
 * @data_free: (nullable): function to free @data
 *
 * Creates a new #NcXcorKernelIntegrand with reference count of 1.
 *
 * Returns: (transfer full): a new #NcXcorKernelIntegrand
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_integrand_new (guint len, void (*eval) (gpointer, gdouble, gdouble *), void (*get_range) (gpointer, gdouble *, gdouble *), gpointer data, GDestroyNotify data_free)
{
  NcXcorKernelIntegrand *integrand = g_new (NcXcorKernelIntegrand, 1);

  integrand->refcount       = 1;
  integrand->len            = len;
  integrand->eval_func      = eval;
  integrand->get_range_func = get_range;
  integrand->data           = data;
  integrand->data_free      = data_free;
  integrand->get_knots_func = NULL;

  integrand->get_range_comp_func = NULL;
  integrand->eval_comps_func     = NULL;
  integrand->get_spectral_func   = NULL;
  integrand->get_panels_func     = NULL;
  integrand->peek_panel_func     = NULL;
  integrand->restrict_func       = NULL;

  integrand->residuals     = NULL;
  integrand->reltol        = 0.0;
  integrand->scaled_abstol = 0.0;

  return integrand;
}

/**
 * nc_xcor_kernel_integrand_set_get_knots: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @get_knots: (scope async): function returning @integrand's knots
 *
 * Declares @integrand as spline-backed, by installing the accessor returning
 * the knots its components are represented on. Left unset by
 * nc_xcor_kernel_integrand_new(), so integrands that are not spline-backed
 * report no knots.
 *
 */
void
nc_xcor_kernel_integrand_set_get_knots (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandGetKnots get_knots)
{
  integrand->get_knots_func = get_knots;
}

/**
 * nc_xcor_kernel_integrand_set_get_spectral: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @get_spectral: (scope async): function reporting a spectral representation
 *
 * Installs the accessor reporting @integrand's Chebyshev expansion. Left unset
 * by nc_xcor_kernel_integrand_new(), in which case @integrand has none and
 * nc_xcor_kernel_integrand_peek_spectral() returns %FALSE.
 *
 */
void
nc_xcor_kernel_integrand_set_get_spectral (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandGetSpectral get_spectral)
{
  integrand->get_spectral_func = get_spectral;
}

/**
 * nc_xcor_kernel_integrand_peek_spectral:
 * @integrand: a #NcXcorKernelIntegrand
 * @coeffs: (out) (transfer none): the coefficient matrix, one row per component
 * @k_min: (out): lower end of the expansion interval
 * @k_max: (out): upper end of the expansion interval
 *
 * Peeks @integrand's Chebyshev expansion, when it has one.
 *
 * A pair of integrands that both report one, over the same interval, can have
 * their outer integral evaluated on the coefficients rather than by quadrature:
 * a product of Chebyshev series is a Chebyshev series, and its integral is a
 * fixed weighted sum of the coefficients.
 *
 * Returns: %TRUE when @integrand carries an expansion
 */
gboolean
nc_xcor_kernel_integrand_peek_spectral (NcXcorKernelIntegrand *integrand, NcmMatrix **coeffs, gdouble *k_min, gdouble *k_max)
{
  if (integrand->get_spectral_func == NULL)
    return FALSE;

  return integrand->get_spectral_func (integrand->data, coeffs, k_min, k_max);
}

/**
 * nc_xcor_kernel_integrand_set_panel_accessors: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @get_panels: (scope async): function reporting the panel count
 * @peek_panel: (scope async): function reporting one panel
 *
 * Installs the accessors enumerating @integrand's panels, for a spectral
 * representation split into more than one.
 *
 */
void
nc_xcor_kernel_integrand_set_panel_accessors (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandGetPanels get_panels, NcXcorKernelIntegrandPeekPanel peek_panel)
{
  integrand->get_panels_func = get_panels;
  integrand->peek_panel_func = peek_panel;
}

/**
 * nc_xcor_kernel_integrand_set_restrict: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @restrict_func: (scope async): function restricting a panel to a subinterval
 *
 * Installs the accessor producing coefficients on a subinterval of a panel.
 *
 */
void
nc_xcor_kernel_integrand_set_restrict (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandRestrict restrict_func)
{
  integrand->restrict_func = restrict_func;
}

/**
 * nc_xcor_kernel_integrand_restrict:
 * @integrand: a #NcXcorKernelIntegrand
 * @a: lower edge of the target interval
 * @b: upper edge of the target interval
 * @coeffs: (out) (transfer full): coefficients on [@a, @b], one row per component
 *
 * Produces @integrand's coefficients on [@a, @b], which has to lie inside a
 * single panel.
 *
 * This is what lets a pair of spectral closures be integrated on the common
 * refinement of their panel edges: on each merged panel both are polynomials
 * over the same interval, so the product is exact and needs no quadrature.
 * Restricting is a change of basis rather than a refit, so it costs arithmetic
 * at panel order rather than fresh radial solves.
 *
 * Returns: %TRUE when @integrand could produce them
 */
gboolean
nc_xcor_kernel_integrand_restrict (NcXcorKernelIntegrand *integrand, gdouble a, gdouble b, NcmMatrix **coeffs)
{
  if (integrand->restrict_func == NULL)
    return FALSE;

  return integrand->restrict_func (integrand->data, a, b, coeffs);
}

/**
 * nc_xcor_kernel_integrand_get_n_panels:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Returns: how many panels @integrand is split into, or 0 when it carries no
 * spectral representation
 */
guint
nc_xcor_kernel_integrand_get_n_panels (NcXcorKernelIntegrand *integrand)
{
  if (integrand->get_panels_func == NULL)
    return 0;

  return integrand->get_panels_func (integrand->data);
}

/**
 * nc_xcor_kernel_integrand_peek_panel:
 * @integrand: a #NcXcorKernelIntegrand
 * @i: panel index, below nc_xcor_kernel_integrand_get_n_panels()
 * @coeffs: (out) (transfer none): the panel's coefficients, one row per component
 * @a: (out): the panel's lower edge
 * @b: (out): the panel's upper edge
 *
 * Peeks one panel. Panels are contiguous and ascending, so panel @i ends where
 * panel @i + 1 begins.
 *
 */
void
nc_xcor_kernel_integrand_peek_panel (NcXcorKernelIntegrand *integrand, guint i, NcmMatrix **coeffs, gdouble *a, gdouble *b)
{
  g_assert (integrand->peek_panel_func != NULL);
  g_assert_cmpuint (i, <, nc_xcor_kernel_integrand_get_n_panels (integrand));

  integrand->peek_panel_func (integrand->data, i, coeffs, a, b);
}

/**
 * nc_xcor_kernel_integrand_set_get_range_comp: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @get_range_comp: (scope async): function returning one component's k range
 *
 * Installs the accessor returning the k range a single component of @integrand
 * is supported on. Left unset by nc_xcor_kernel_integrand_new(), in which case
 * every component reports the whole range.
 *
 */
void
nc_xcor_kernel_integrand_set_get_range_comp (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandGetRangeComp get_range_comp)
{
  integrand->get_range_comp_func = get_range_comp;
}

/**
 * nc_xcor_kernel_integrand_set_eval_comps: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @eval_comps: (scope async): function evaluating a run of components
 *
 * Installs the accessor evaluating a contiguous run of @integrand's
 * components, for callers that integrate the run on its own. Left unset by
 * nc_xcor_kernel_integrand_new(), in which case a run is served by evaluating
 * every component.
 *
 */
void
nc_xcor_kernel_integrand_set_eval_comps (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandEvalComps eval_comps)
{
  integrand->eval_comps_func = eval_comps;
}

/**
 * nc_xcor_kernel_integrand_peek_knots:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Peeks the knots @integrand's components are represented on, shared by every
 * component (multipole) it carries, or %NULL when @integrand is not
 * spline-backed.
 *
 * These knots are what makes the outer $k$ integral exactly integrable, and
 * are why %NC_XCOR_METHOD_KERNEL_EXACT needs no tolerance. Each component is a
 * cubic spline in $k$, so on any interval over which both members of a pair
 * are a single cubic piece, the product $k^2 W_i(k) W_j(k)$ entering $C_\ell$
 * is a polynomial of degree $8$ and a $5$-node Gauss-Legendre rule integrates
 * it exactly.
 *
 * Two kernels are sampled independently and so do not share an abscissa: the
 * intervals with that property are the panels of the *common refinement* of
 * the two knot sets. Merging two sorted knot vectors is all the coupling the
 * argument needs -- sampling the kernels jointly onto one shared abscissa
 * would also work but costs about twice as much to produce, and forces every
 * pair to carry every kernel's knots.
 *
 * Returns: (transfer none) (nullable): the knot vector, or %NULL.
 */

/**
 * nc_xcor_kernel_integrand_set_tolerances:
 * @integrand: a #NcXcorKernelIntegrand
 * @reltol: the relative half of the fit criterion
 * @scaled_abstol: its floor, as a fraction of the fitted function's own peak
 *
 * Records the criterion @integrand was fitted to, in the two parts it actually
 * has. ncm_function_sample_set_refine() accepts a point when
 *
 * $$ \Vert f - \tilde f \Vert_2 \le \mathrm{reltol} \Vert f \Vert_2 + a \Vert f \Vert_2^\mathrm{max} $$
 *
 * so the two govern different regions: the floor is *added*, not maxed, and is
 * scaled to the peak, which leaves @reltol biting only where the function is
 * within a few orders of that peak. Keeping them apart matters downstream --
 * they reach a $C_\ell$ differently, the floor through the product of two
 * closures. See nc_xcor_compute_full().
 *
 * A spline-backed integrand is a fit, not the function, and the quadratures
 * that consume it have no other way to learn how good a fit.
 * %NC_XCOR_METHOD_KERNEL_EXACT integrates it *exactly*, so this is the whole of
 * its error budget.
 *
 */
void
nc_xcor_kernel_integrand_set_tolerances (NcXcorKernelIntegrand *integrand, gdouble reltol, gdouble scaled_abstol)
{
  g_return_if_fail (integrand != NULL);
  g_return_if_fail (reltol >= 0.0);
  g_return_if_fail (scaled_abstol >= 0.0);

  integrand->reltol        = reltol;
  integrand->scaled_abstol = scaled_abstol;
}

/**
 * nc_xcor_kernel_integrand_get_reltol:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Returns: the relative half of the fit criterion, or 0.0 when exact or
 * unknown. See nc_xcor_kernel_integrand_set_tolerances().
 */
gdouble
nc_xcor_kernel_integrand_get_reltol (NcXcorKernelIntegrand *integrand)
{
  g_return_val_if_fail (integrand != NULL, 0.0);

  return integrand->reltol;
}

/**
 * nc_xcor_kernel_integrand_get_scaled_abstol:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Returns: the floor of the fit criterion as a fraction of the fitted
 * function's own peak, or 0.0 when there was none. See
 * nc_xcor_kernel_integrand_set_tolerances().
 */
gdouble
nc_xcor_kernel_integrand_get_scaled_abstol (NcXcorKernelIntegrand *integrand)
{
  g_return_val_if_fail (integrand != NULL, 0.0);

  return integrand->scaled_abstol;
}

/**
 * nc_xcor_kernel_integrand_set_residuals:
 * @integrand: a #NcXcorKernelIntegrand
 * @residuals: (nullable): the achieved fit residuals, or %NULL
 *
 * Records the residual the fit *achieved* on each knot interval, one row per
 * knot of nc_xcor_kernel_integrand_peek_knots() and one column per component,
 * as produced by ncm_function_sample_set_get_residuals().
 *
 * This is the sharper half of the pair with
 * nc_xcor_kernel_integrand_set_tolerances(): the tolerances say what the fit
 * was *asked* for, and refinement beats its own request by orders -- by 12 to
 * 3100 depending on the kernel, which is enough for an error built from the
 * tolerances alone to be unable to tell a well-fitted pair from a badly
 * fitted one. Where these residuals are present %NC_XCOR_METHOD_KERNEL_EXACT
 * uses them and falls back to the tolerances only on intervals that carry no
 * record (NaN). See nc_xcor_compute_full().
 *
 */
void
nc_xcor_kernel_integrand_set_residuals (NcXcorKernelIntegrand *integrand, NcmMatrix *residuals)
{
  g_return_if_fail (integrand != NULL);

  ncm_matrix_clear (&integrand->residuals);

  if (residuals != NULL)
    integrand->residuals = ncm_matrix_ref (residuals);
}

/**
 * nc_xcor_kernel_integrand_peek_residuals:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Peeks the achieved fit residuals, or %NULL when the closure was built
 * without residual tracking. See
 * nc_xcor_kernel_integrand_set_residuals().
 *
 * Returns: (transfer none) (nullable): the residual matrix, or %NULL
 */
NcmMatrix *
nc_xcor_kernel_integrand_peek_residuals (NcXcorKernelIntegrand *integrand)
{
  g_return_val_if_fail (integrand != NULL, NULL);

  return integrand->residuals;
}

NcmVector *
nc_xcor_kernel_integrand_peek_knots (NcXcorKernelIntegrand *integrand)
{
  if (integrand->get_knots_func == NULL)
    return NULL;

  return integrand->get_knots_func (integrand->data);
}

/**
 * nc_xcor_kernel_integrand_ref:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Increases the reference count of @integrand by one atomically.
 *
 * Returns: (transfer full): @integrand
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_integrand_ref (NcXcorKernelIntegrand *integrand)
{
  g_atomic_int_inc (&integrand->refcount);

  return integrand;
}

/**
 * nc_xcor_kernel_integrand_unref:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Decreases the reference count of @integrand by one atomically.
 * When the reference count reaches zero, frees @integrand and its
 * associated data using the free function provided at creation time
 * (if any).
 */
void
nc_xcor_kernel_integrand_unref (NcXcorKernelIntegrand *integrand)
{
  if (g_atomic_int_dec_and_test (&integrand->refcount))
  {
    if (integrand->data_free != NULL)
      integrand->data_free (integrand->data);

    ncm_matrix_clear (&integrand->residuals);

    g_free (integrand);
  }
}

/**
 * nc_xcor_kernel_integrand_clear:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * If *@integrand is not %NULL, decreases its reference count and
 * sets the pointer to %NULL.
 */
void
nc_xcor_kernel_integrand_clear (NcXcorKernelIntegrand **integrand)
{
  if (*integrand != NULL)
  {
    nc_xcor_kernel_integrand_unref (*integrand);
    *integrand = NULL;
  }
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
    gdouble xi_min, xi_max, k_min, k_max;

    nc_xcor_kernel_component_get_limits (comp, cosmo, &xi_min, &xi_max, &k_min, &k_max);

    {
      const gdouble k_min_limb = nu / xi_max;
      const gdouble k_max_limb = nu / xi_min;

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
 * @closure_type: how to represent the sampled window, see #NcXcor:closure-type
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
 * @closure_type: how to represent the sampled window, see #NcXcor:closure-type
 *
 * Gets a vectorized evaluation function for the kernel over a range of multipoles.
 * The returned integrand will have len = lmax - lmin + 1, and will evaluate all
 * multipoles in the range [lmin, lmax] simultaneously.
 *
 * Uses the base class implementation which checks the l-limber property:
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
 * @closure_type: how to represent the sampled window, see #NcXcor:closure-type
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
 * @sbi is unused below the kernel's l-limber threshold, where no spherical
 * Bessel integral is performed. @closure_type is likewise unused there: a
 * Limber window carries a step per multipole and only the spline closure
 * represents that.
 *
 * Returns: (transfer full): the kernel integrand over [@lmin, @lmax]
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_get_eval_vectorized_full (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi, NcXcorKernelClosure closure_type)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  if ((self->l_limber == 0) || ((self->l_limber > 0) && (lmin >= self->l_limber)))
    return _nc_xcor_kernel_build_limber_integrand (xclk, cosmo, lmin, lmax);
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
 * provide more accurate spline interpolation at the cost of more sample points.
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
 * nc_xcor_kernel_get_scaled_abstol:
 * @xclk: a #NcXcorKernel
 *
 *
 */
gdouble
nc_xcor_kernel_get_scaled_abstol (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->scaled_abstol;
}

/**
 * nc_xcor_kernel_set_scaled_abstol:
 * @xclk: a #NcXcorKernel
 * @scaled_abstol: the absolute minimum, as a fraction of the peak (must be > 0)
 *
 * Sets the absolute minimum threshold for adaptive midpoint refinement. This parameter
 * helps prevent excessive refinement in cases where the kernel has very low amplitude,
 * by providing a floor below which the refinement will stop regardless of the relative
 * tolerance.
 *
 * Values below %NC_XCOR_KERNEL_MIN_USEFUL_SCALED_ABSTOL are accepted but warned about:
 * the floor enters the $C_\ell$ integrand squared, so they ask for accuracy the outer
 * integral cannot carry and pay for it in spline knots. See
 * #NcXcorKernel:scaled-abstol.
 */
void
nc_xcor_kernel_set_scaled_abstol (NcXcorKernel *xclk, gdouble scaled_abstol)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert_cmpfloat (scaled_abstol, >, 0.0);

  if (scaled_abstol < NC_XCOR_KERNEL_MIN_USEFUL_SCALED_ABSTOL)
    g_warning ("nc_xcor_kernel_set_scaled_abstol: %.3e is below the useful floor of %.0e. "
               "This tolerance is measured against the peak of W(k), but the C_l integrand "
               "is k^2 W_a W_b, so it enters squared: %.3e here is %.3e on the integrand, "
               "past what the outer integral carries. It cannot improve the result and can "
               "cost orders of magnitude in spline knots.",
               scaled_abstol, NC_XCOR_KERNEL_MIN_USEFUL_SCALED_ABSTOL,
               scaled_abstol, scaled_abstol * scaled_abstol);

  self->scaled_abstol = scaled_abstol;
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
 * Sets the maximum number of adaptive midpoint refinement iterations. Higher
 * values allow for more refinement passes when needed. Typical values range
 * from 1000 to 100000.
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
 * Gets the expansion factor used for domain extension in the non-Limber
 * integrand construction. This determines how much the domain is extended
 * in each iteration.
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
 * Sets the expansion factor for domain extension. Larger values result in
 * more aggressive expansion. Typical values range from 0.1 to 0.5.
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
 * Sets #NcXcorKernel:panel-order-cap. Read when a closure is built, so one
 * already built keeps the panels it was built with.
 *
 */
void
nc_xcor_kernel_set_panel_order_cap (NcXcorKernel *xclk, guint panel_order_cap)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->panel_order_cap = panel_order_cap;
}

/**
 * nc_xcor_kernel_get_track_fit_residual:
 * @xclk: a #NcXcorKernel
 *
 * Returns: whether the closure records the residual its fit achieved. See
 * #NcXcorKernel:track-fit-residual.
 */
gboolean
nc_xcor_kernel_get_track_fit_residual (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return self->track_fit_residual;
}

/**
 * nc_xcor_kernel_set_track_fit_residual:
 * @xclk: a #NcXcorKernel
 * @track_fit_residual: whether to record the achieved residual
 *
 * Sets #NcXcorKernel:track-fit-residual. It is read when a closure is built,
 * so a closure already built keeps whatever it was built with.
 *
 */
void
nc_xcor_kernel_set_track_fit_residual (NcXcorKernel *xclk, gboolean track_fit_residual)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->track_fit_residual = track_fit_residual;
}

/**
 * nc_xcor_kernel_eval_limber_z: (virtual eval_limber_z)
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @z: a #gdouble
 * @xck: a #NcXcorKinetic
 * @l: a #gint
 *
 * Evaluates the Limber kernel at redshift @z for multipole @l.
 * The kinetic quantities (comoving distance and Hubble parameter) are
 * provided in @xck. Returns zero if @z is outside the kernel's redshift range.
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
 * Evaluates the Limber kernel at redshift @z for multipole @l, including
 * the normalization factor. This function computes the kinetic quantities
 * internally using @dist and applies the kernel's constant factor.
 *
 * Returns: the normalized kernel value $c \times W(z,\ell)$
 */
gdouble
nc_xcor_kernel_eval_limber_z_full (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, NcDistance *dist, gint l)
{
  const gdouble xi_z      = nc_distance_comoving (dist, cosmo, z); /* in units of Hubble radius */
  const gdouble E_z       = nc_hicosmo_E (cosmo, z);
  const NcXcorKinetic xck = { xi_z, E_z };
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
 * Prepares the kernel for evaluation with the given cosmological model.
 * This may involve precomputing quantities that depend on @cosmo but not
 * on redshift or multipole.
 *
 */
void
nc_xcor_kernel_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  return NC_XCOR_KERNEL_GET_CLASS (xclk)->prepare (xclk, cosmo);
}

/**
 * nc_xcor_kernel_get_component_list: (virtual get_component_list)
 * @xclk: a #NcXcorKernel
 *
 * Gets the list of components that make up this kernel.
 *
 * Returns: (transfer container) (element-type NcXcorKernelComponent): a #GPtrArray of #NcXcorKernelComponent
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
 * Logs all registered #NcXcorLimberKernel subclasses to the message log.
 * This is useful for debugging and discovering available kernel implementations.
 *
 */
void
nc_xcor_kernel_log_all_models (void)
{
  g_message ("# Registered NcXcorKernel:%s are:\n",
             g_type_name (NC_TYPE_XCOR_KERNEL));
  _nc_xcor_kernel_log_all_models_go (NC_TYPE_XCOR_KERNEL, 0);
}

/**
 * nc_xcor_kernel_integrand_get_range:
 * @integrand: a #NcXcorKernelIntegrand
 * @k_min: (out): minimum k value
 * @k_max: (out): maximum k value
 *
 * Gets the valid k range for this integrand.
 */
/**
 * nc_xcor_kernel_integrand_get_range_comp:
 * @integrand: a #NcXcorKernelIntegrand
 * @i: component index
 * @k_min: (out): minimum k value
 * @k_max: (out): maximum k value
 *
 * Gets the k range component @i is supported on, which can be a part of the
 * range nc_xcor_kernel_integrand_get_range() reports for the whole integrand:
 * a block of multipoles shares one domain, and under the Limber approximation
 * each of them vanishes outside its own band within it. Integrating a
 * component over its own range keeps that band edge on an integration limit
 * instead of leaving a step inside the interval.
 *
 * Falls back to the whole range for integrands that do not distinguish their
 * components.
 */
/**
 * nc_xcor_kernel_integrand_eval_comps: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @k: wavenumber
 * @offset: index of the first component to evaluate
 * @len: number of components to evaluate
 * @W: (array) (out caller-allocates): full-length array to store results in
 *
 * Evaluates components [@offset, @offset + @len) at wavenumber @k, writing
 * them at their own indices in @W. Integrands that can only evaluate every
 * component at once do so, filling the whole of @W; either way the entries
 * the caller asked for are valid.
 */
/**
 * nc_xcor_kernel_integrand_eval: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @k: wavenumber
 * @W: (array) (out caller-allocates): array of length @len to store results
 *
 * Evaluates the integrand at wavenumber @k, storing @len results in @W.
 */
/**
 * nc_xcor_kernel_integrand_get_len:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Gets the number of components in the integrand.
 *
 * Returns: the number of components
 */
/**
 * nc_xcor_kernel_integrand_eval_array:
 * @integrand: a #NcXcorKernelIntegrand
 * @k: wavenumber
 *
 * Evaluates the integrand at wavenumber @k and returns the results
 * in a newly allocated #GArray. This is a convenience wrapper around
 * nc_xcor_kernel_integrand_eval() that handles array allocation.
 *
 * Returns: (transfer full) (element-type gdouble): a #GArray containing @len #gdouble values
 */

