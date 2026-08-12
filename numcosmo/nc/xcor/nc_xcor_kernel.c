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
  PROP_SIZE,
};


G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcXcorKernel, nc_xcor_kernel, NCM_TYPE_MODEL)
G_DEFINE_BOXED_TYPE (NcXcorKinetic, nc_xcor_kinetic, nc_xcor_kinetic_copy, nc_xcor_kinetic_free)
G_DEFINE_BOXED_TYPE (NcXcorKernelIntegrand, nc_xcor_kernel_integrand, nc_xcor_kernel_integrand_ref, nc_xcor_kernel_integrand_unref)

static void
nc_xcor_kernel_init (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  self->dist                    = NULL;
  self->ps                      = NULL;
  self->sbi                     = NULL;
  self->lmax                    = 0;
  self->l_limber                = 0;
  self->adaptive_epsilon        = 0.0;
  self->adaptive_boundary_tries = 0;
  self->reltol                  = 0.0;
  self->scaled_abstol           = 0.0;
  self->max_border_expansions   = 0;
  self->max_iter                = 0;
  self->expansion_factor        = 0.0;
  self->constructed             = FALSE;
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
   * Convergence threshold for the adaptive $k$-range determination: sampling
   * stops extending outwards once $\Vert W\Vert$ falls below this fraction of
   * its running maximum, for #NcXcorKernel:adaptive-boundary-tries consecutive
   * steps.
   *
   * **Tightening this is not free, and not monotonic in usefulness.** The test
   * is on an amplitude relative to a peak, so on a tail decaying as $k^{-p}$ it
   * truncates around $(1/\epsilon)^{1/p}$ times the peak location. Kernels with
   * compact support in $z$ -- a top-hat redshift bin, for instance -- have a
   * transform that *rings*, decaying as a slow power law rather than
   * exponentially, so a smaller $\epsilon$ buys a great deal of extra range in
   * a region that is purely oscillatory. Measured on a top-hat bin at
   * $\ell = 0$: 71% of the knots lay beyond the point at which 99.9% of the
   * $C_\ell$ integral had already accumulated.
   *
   * Smooth kernels (galaxy, weak lensing, CMB lensing, tSZ) do not show this:
   * their amplitude and their contribution to the integral fall off together,
   * so the criterion cuts where it should.
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
   * Only active where #NcXcorKernel:scaled-abstol does not bind, which for
   * realistic settings is most of the domain -- so this is usually *not* the
   * knob that governs accuracy. See #NcXcorKernel:scaled-abstol before
   * adjusting either.
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
   * Absolute floor for the adaptive refinement of the $W_i(k)$ spline, in units
   * of the peak of the sampled integrand: the refinement stops subdividing an
   * interval once its error falls below
   * `scaled-abstol` $\times \max\vert F\vert$.
   *
   * **This, not #NcXcorKernel:reltol, is what limits the accuracy of the
   * resulting $C_\ell$.** While this floor binds, tightening `reltol` changes
   * nothing at all. Because the floor is absolute, it produces a roughly
   * constant absolute error on every $C_\ell$ of a block, so entries that are
   * themselves small carry a correspondingly large *relative* error --
   * generally on exactly the pairs that contribute least.
   *
   * **Do not tighten it reflexively.** Two effects bite, in order:
   *
   * 1. The knot count grows steeply, since refinement chases the oscillatory
   *    tail described under #NcXcorKernel:adaptive-epsilon. Measured on seven
   *    top-hat bins at $\ell = 32$: 493 knots per kernel at `1e-4`, 2367 at
   *    `1e-6`, 4573 at `1e-7`, ~8759 at `1e-8`.
   * 2. The *outer* $k$ integral then has a genuinely oscillatory integrand to
   *    handle over several decades, and adaptive quadrature stops converging:
   *    in the same test %NC_XCOR_METHOD_KERNEL_CUBATURE succeeds through `1e-7`
   *    and aborts at `1e-8`. The abort is a true report about the integrand, and
   *    the integrand was made difficult by the refinement, not by the physics.
   *
   * There is also a hard floor below which the request is meaningless. The
   * radial integral is computed to an absolute tolerance of
   * %NC_XCOR_KERNEL_INTEG_ABSTOL_FRAC times its own running maximum, so where
   * $\vert W\vert$ has fallen far below its peak its *relative* accuracy is
   * already poor -- the value there is a small residual of a strongly
   * cancelling oscillatory integral. Asking the spline to resolve that region
   * to a tighter relative precision than the radial integral can deliver
   * subdivides against round-off, buying knots and no information.
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

  /* Per-component support. Components of a joint (multi-kernel) integrand
   * share one knot set but not one domain, and NcmSpline does not range-check
   * on evaluation, so each component's own support has to be carried
   * alongside. See nc_xcor_kernel_integrand_get_component_range(). */
  gdouble *supp_k_min;
  gdouble *supp_k_max;
} SplineIntegrandData;

static NcmVector *
_spline_integrand_get_knots (gpointer data)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  /* Every component of the spline vec is sampled on the same abscissa, so the
   * first component's knots are the knots of the whole block. */
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

static void
_spline_integrand_get_range (gpointer data, gdouble *kmin, gdouble *kmax)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  *kmin = sid->k_min;
  *kmax = sid->k_max;
}

static void
_spline_integrand_data_free (gpointer data)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  nc_hicosmo_clear (&sid->cosmo);
  ncm_spline_vec_clear (&sid->spline_vec);
  ncm_vector_clear (&sid->eval_result);
  g_clear_pointer (&sid->supp_k_min, g_free);
  g_clear_pointer (&sid->supp_k_max, g_free);

  g_free (data);
}

/*
 * Records, for each component, the outermost knots at which it still carries
 * signal: the first and last knot where |y| >= epsilon * max|y| over the
 * component. Outside that the component contributed nothing to the sampling
 * and its spline holds only whatever the cubic does between negligible
 * values -- which, weighted by k^2 in the outer integral, is not negligible at
 * all. This is the same criterion the per-kernel adaptive boundary search
 * converges to, applied after the fact to the sampled data so that it does not
 * depend on the order in which the sampler explored.
 */
static void
_spline_integrand_compute_supports (SplineIntegrandData *sid, const gdouble epsilon)
{
  const guint n_knots = ncm_spline_vec_get_nknots (sid->spline_vec);
  NcmVector *xv       = ncm_spline_peek_xv (ncm_spline_vec_peek_spline (sid->spline_vec, 0));
  guint c;

  sid->supp_k_min = g_new (gdouble, sid->len);
  sid->supp_k_max = g_new (gdouble, sid->len);

  for (c = 0; c < sid->len; c++)
  {
    NcmVector *yv    = ncm_spline_peek_yv (ncm_spline_vec_peek_spline (sid->spline_vec, c));
    gdouble max_absy = 0.0;
    gdouble cut;
    guint i;
    gint lo = -1;
    gint hi = -1;

    for (i = 0; i < n_knots; i++)
      max_absy = GSL_MAX (max_absy, fabs (ncm_vector_get (yv, i)));

    cut = epsilon * max_absy;

    for (i = 0; i < n_knots; i++)
    {
      if (fabs (ncm_vector_get (yv, i)) >= cut)
      {
        if (lo < 0)
          lo = i;

        hi = i;
      }
    }

    if (lo < 0) /* identically zero component: empty support */
    {
      sid->supp_k_min[c] = sid->k_max;
      sid->supp_k_max[c] = sid->k_min;
    }
    else
    {
      /* Widen by one knot on each side: the criterion is evaluated at knots,
       * and the last panel before the cut still carries signal. */
      const guint lo_u = (lo > 0) ? (guint) (lo - 1) : 0;
      const guint hi_u = ((guint) hi + 1 < n_knots) ? (guint) (hi + 1) : n_knots - 1;

      sid->supp_k_min[c] = ncm_vector_get (xv, lo_u);
      sid->supp_k_max[c] = ncm_vector_get (xv, hi_u);
    }
  }
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

/* Fraction of a component's largest spherical-Bessel integral, over all k,
 * below which a further integral cannot affect the k-spline built from them. */
#define NC_XCOR_KERNEL_INTEG_ABSTOL_FRAC 1.0e-16

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
  gdouble integ_max;                       /* Running max |integrator output| over k, pre-prefactor */
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
  state->integ_max            = 0.0;
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

  g_assert_cmpuint (n_l, <=, MAX_ELL_BLOCK);
  g_assert_cmpuint (comp_list->len, <=, MAX_COMP_BLOCK);

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
    .adaptive_boundary_tries = self->adaptive_boundary_tries
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
    .adaptive_boundary_tries = self->adaptive_boundary_tries
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
      const gdouble nu_j    = l_j + 0.5;
      const gdouble k_min_j = state->k_min_limber_ell[j];
      const gdouble k_max_j = state->k_max_limber_ell[j];

      /* Evaluate at left boundary */
      {
        const gdouble xi_left         = nu_j / k_min_j;
        const gdouble limber_k_left   = 1.0 / k_min_j;
        const gdouble prefactor_limb  = sqrt (M_PI / (2.0 * nu_j));
        const gdouble kernel_val_left = nc_xcor_kernel_component_eval_kernel (state->comp, cosmo, xi_left, k_min_j);
        const gdouble prefactor_left  = nc_xcor_kernel_component_eval_prefactor (state->comp, cosmo, k_min_j, l_j);

        state->last_values_left[j] = prefactor_limb * prefactor_left * limber_k_left * kernel_val_left;
      }

      /* Evaluate at right boundary */
      {
        const gdouble xi_right         = nu_j / k_max_j;
        const gdouble limber_k_right   = 1.0 / k_max_j;
        const gdouble prefactor_limb   = sqrt (M_PI / (2.0 * nu_j));
        const gdouble kernel_val_right = nc_xcor_kernel_component_eval_kernel (state->comp, cosmo, xi_right, k_max_j);
        const gdouble prefactor_right  = nc_xcor_kernel_component_eval_prefactor (state->comp, cosmo, k_max_j, l_j);

        state->last_values_right[j] = prefactor_limb * prefactor_right * limber_k_right * kernel_val_right;
      }
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
      {
        /* The integrator cannot know the scale its result feeds into: this
         * component's kernel over all k. Below a k-independent floor tied to
         * that scale the result cannot move the k-spline built from it, and
         * chasing relative accuracy there is both wasted and, where the
         * integrand has an unresolvable feature, unattainable. */
        ncm_sbessel_integrator_set_abstol (comp_states->sbi, NC_XCOR_KERNEL_INTEG_ABSTOL_FRAC * state->integ_max);

        ncm_sbessel_integrator_integrate (
          comp_states->sbi, _nc_xcor_kernel_component_kernel_integ,
          state->xi_min, state->xi_max, k, integ_result, &state->params
        );
      }

      for (i = 0; i < comp_states->n_l; i++)
        state->integ_max = GSL_MAX (state->integ_max, fabs (kernel_out[ci][i]));

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
      const gdouble nu            = l + 0.5;
      const gboolean within_range = (k >= state->k_min_limber_ell[i]) && (k <= state->k_max_limber_ell[i]);

      if (within_range)
      {
        /* Normal Limber evaluation within valid range */
        const gdouble xi               = nu / k;
        const gdouble limber_k         = 1.0 / k;
        const gdouble prefactor_limber = sqrt (M_PI / (2.0 * nu));
        const gdouble kernel_val       = nc_xcor_kernel_component_eval_kernel (state->comp, state->params.cosmo, xi, k);
        const gdouble prefactor        = nc_xcor_kernel_component_eval_prefactor (state->comp, state->params.cosmo, k, l);

        kernel_out[ci][i] = prefactor_limber * prefactor * limber_k * kernel_val;
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
 * Joint (multi-kernel) sampling.
 *
 * Every kernel keeps its own ComponentStates -- the per-kernel adaptive
 * boundary tracking, integrator abstol and extrapolation logic are unchanged.
 * What is shared is the *abscissa*: one NcmFunctionSampleSet carrying
 * n_kernels * n_l components, so a single knot set represents every kernel of
 * the block. That is what lets all kernel pairs be integrated on common
 * panels (and hence be assembled as one matrix product) instead of on the
 * union of two independently-adapted knot sets, which measurement shows is
 * essentially the *sum* of the two: independently adapted splines resolve the
 * same oscillation at nearby but distinct abscissas, so unioning them buys
 * nothing.
 */
typedef struct _JointComputeData
{
  ComponentStates *per_kernel;
  guint n_kernels;
  guint n_l;
} JointComputeData;

static void
_joint_states_compute_non_limber (const gdouble k, NcmVector *y, gpointer user_data)
{
  JointComputeData *joint = (JointComputeData *) user_data;
  guint ik;

  g_assert_cmpuint (ncm_vector_len (y), ==, joint->n_kernels * joint->n_l);

  for (ik = 0; ik < joint->n_kernels; ik++)
  {
    NcmVector *sub = ncm_vector_get_subvector (y, ik * joint->n_l, joint->n_l);

    _component_states_compute_non_limber (k, sub, &joint->per_kernel[ik]);

    ncm_vector_free (sub);
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
    guint i;

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

    sid->spline_vec  = ncm_function_sample_set_to_spline_vec (fss, spline);
    sid->k_min       = ncm_function_sample_set_get_x_min (fss);
    sid->k_max       = ncm_function_sample_set_get_x_max (fss);
    sid->lmin        = lmin;
    sid->len         = n_l;
    sid->RH_Mpc      = nc_hicosmo_RH_Mpc (cosmo);
    sid->cosmo       = nc_hicosmo_ref (cosmo);
    sid->eval_result = ncm_vector_new (n_l);

    _spline_integrand_compute_supports (sid, comp_states->epsilon);

    ncm_function_sample_set_clear (&fss);
    ncm_spline_free (spline);

    {
      NcXcorKernelIntegrand *integrand = nc_xcor_kernel_integrand_new (n_l,
                                                                       _spline_integrand_eval,
                                                                       _spline_integrand_get_range,
                                                                       sid,
                                                                       _spline_integrand_data_free);

      nc_xcor_kernel_integrand_set_get_knots (integrand, _spline_integrand_get_knots);

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

    return _nc_xcor_kernel_build_spline_integrand (xclk, cosmo, lmin, lmax,
                                                   &comp_states,
                                                   _component_states_compute_limber,
                                                   self->reltol, self->scaled_abstol);
  }
}

static NcXcorKernelIntegrand *
_nc_xcor_kernel_build_non_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi)
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

    return _nc_xcor_kernel_build_spline_integrand (xclk, cosmo, lmin, lmax,
                                                   &comp_states,
                                                   _component_states_compute_non_limber,
                                                   self->reltol, self->scaled_abstol);
  }
}

/**
 * nc_xcor_kernel_get_eval_vectorized_joint:
 * @kernels: (element-type NcXcorKernel): the kernels to sample together
 * @cosmo: a #NcHICosmo
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 * @sbi: a #NcmSBesselIntegrator
 *
 * Builds one integrand carrying every kernel of @kernels over
 * [@lmin, @lmax], represented on a *single shared knot set*. The returned
 * integrand has `len = n_kernels * n_l` components, ordered kernel-major:
 * component `ik * n_l + il` is kernel `ik` at multipole `lmin + il`.
 *
 * Sampling jointly is what makes the outer $k$ integral exactly and cheaply
 * computable for every kernel pair at once: each component is a cubic spline
 * on the *same* abscissa, so one $5$-node Gauss-Legendre rule per knot panel
 * integrates $k^2 W_i W_j$ exactly for all pairs simultaneously (see
 * nc_xcor_kernel_integrand_peek_knots()), and the pair assembly reduces to a
 * single matrix product. Building each kernel separately instead forces every
 * pair onto the union of two distinct knot sets, which costs as much as the
 * sum of the two.
 *
 * The joint refinement uses the strictest tolerances of the kernels involved,
 * so no kernel is sampled more loosely than it would have been on its own.
 *
 * Returns: (transfer full): the joint integrand.
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_get_eval_vectorized_joint (GPtrArray *kernels, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi)
{
  const guint n_kernels = (kernels != NULL) ? kernels->len : 0;
  const guint n_l       = lmax - lmin + 1;
  JointComputeData joint;
  SplineIntegrandData *sid;
  gdouble reltol              = GSL_POSINF;
  gdouble scaled_abstol       = GSL_POSINF;
  gdouble expansion_factor    = GSL_POSINF;
  guint max_border_expansions = 0;
  guint max_iter              = 0;
  gdouble k_min_hard          = GSL_POSINF;
  gdouble k_max_hard          = 0.0;
  guint ik;

  g_assert_cmpuint (n_kernels, >, 0);
  g_assert_cmpint (lmax, >=, lmin);

  if (sbi == NULL)
    g_error ("nc_xcor_kernel_get_eval_vectorized_joint: no integrator given.");

  if (n_l > NC_XCOR_KERNEL_MAX_ELL_BLOCK)
    g_error ("nc_xcor_kernel_get_eval_vectorized_joint: ell block %u exceeds the "
             "maximum %d.", n_l, NC_XCOR_KERNEL_MAX_ELL_BLOCK);

  ncm_sbessel_integrator_set_ell_range (sbi, lmin, lmax);

  joint.per_kernel = g_new0 (ComponentStates, n_kernels);
  joint.n_kernels  = n_kernels;
  joint.n_l        = n_l;

  for (ik = 0; ik < n_kernels; ik++)
  {
    NcXcorKernel *xclk        = NC_XCOR_KERNEL (g_ptr_array_index (kernels, ik));
    NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);
    GPtrArray *comp_list      = _nc_xcor_kernel_validate_component_list (xclk, n_l);

    if (comp_list == NULL)
      g_error ("nc_xcor_kernel_get_eval_vectorized_joint: kernel %s has no components.",
               G_OBJECT_TYPE_NAME (xclk));

    {
      /* Struct assignment is unavailable: ComponentStates carries const
       * members, so it can only be initialized, never assigned. The
       * destination is freshly allocated storage, not a const object. */
      ComponentStates cs = _component_states_init_non_limber (xclk, lmin, n_l, comp_list, cosmo, sbi);

      memcpy (&joint.per_kernel[ik], &cs, sizeof (ComponentStates));
    }

    g_ptr_array_unref (comp_list);

    /* Strictest tolerance across the group: joint sampling must never relax
     * what any single kernel would have demanded of itself. */
    reltol                = GSL_MIN (reltol, self->reltol);
    scaled_abstol         = GSL_MIN (scaled_abstol, self->scaled_abstol);
    expansion_factor      = GSL_MIN (expansion_factor, self->expansion_factor);
    max_border_expansions = GSL_MAX (max_border_expansions, self->max_border_expansions);
    max_iter              = GSL_MAX (max_iter, self->max_iter);

    /* The shared domain is the union: a component contributing nothing over
     * part of it is represented there anyway, at negligible cost, rather than
     * truncating the domain of the kernels that do contribute. */
    k_min_hard = GSL_MIN (k_min_hard, joint.per_kernel[ik].k_min_hard);
    k_max_hard = GSL_MAX (k_max_hard, joint.per_kernel[ik].k_max_hard);
  }

  sid = g_new0 (SplineIntegrandData, 1);

  {
    NcmFunctionSampleSet *fss = ncm_function_sample_set_new (n_kernels * n_l);
    NcmSpline *spline         = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    GArray *k_seeds           = g_array_new (FALSE, FALSE, sizeof (gdouble));
    guint i;

    /* Seed with every kernel's own scales, so no kernel's structure is missed
     * by the initial sampling. */
    for (ik = 0; ik < n_kernels; ik++)
    {
      GArray *kernel_seeds = g_array_new (FALSE, FALSE, sizeof (gdouble));

      _component_states_compute_k_seeds (&joint.per_kernel[ik], kernel_seeds);

      for (i = 0; i < kernel_seeds->len; i++)
      {
        const gdouble k_seed = g_array_index (kernel_seeds, gdouble, i);

        if (_is_new_k (k_seed, (gdouble *) k_seeds->data, k_seeds->len))
          g_array_append_val (k_seeds, k_seed);
      }

      g_array_unref (kernel_seeds);
    }

    for (i = 0; i < k_seeds->len; i++)
    {
      const gdouble k_seed = g_array_index (k_seeds, gdouble, i);

      ncm_function_sample_set_add_old_func (fss, k_seed, _joint_states_compute_non_limber, &joint);
    }

    g_array_unref (k_seeds);

    ncm_function_sample_set_expand_domain (
      fss,
      _joint_states_compute_non_limber,
      k_min_hard,
      k_max_hard,
      expansion_factor,
      joint.per_kernel[0].epsilon,
      max_border_expansions,
      joint.per_kernel[0].adaptive_boundary_tries,
      &joint
    );

    ncm_function_sample_set_mark_all_old (fss);
    ncm_function_sample_set_reset_interval_ok (fss);

    {
      /* get_absmaxF_min() is the minimum over components, so pooling kernels
       * tightens the floor rather than loosening it: the weakest component of
       * the group sets it. Accuracy is safe by construction; the cost is that
       * strong components are refined more than they need. */
      const gdouble max_absF_total = ncm_function_sample_set_get_absmaxF_min (fss);

      ncm_function_sample_set_adaptive_midpoint (
        fss, _joint_states_compute_non_limber,
        reltol, max_absF_total * scaled_abstol, max_iter, 1,
        spline, &joint
      );
    }

    sid->spline_vec  = ncm_function_sample_set_to_spline_vec (fss, spline);
    sid->k_min       = ncm_function_sample_set_get_x_min (fss);
    sid->k_max       = ncm_function_sample_set_get_x_max (fss);
    sid->lmin        = lmin;
    sid->len         = n_kernels * n_l;
    sid->RH_Mpc      = nc_hicosmo_RH_Mpc (cosmo);
    sid->cosmo       = nc_hicosmo_ref (cosmo);
    sid->eval_result = ncm_vector_new (n_kernels * n_l);

    _spline_integrand_compute_supports (sid, joint.per_kernel[0].epsilon);

    ncm_function_sample_set_clear (&fss);
    ncm_spline_free (spline);
  }

  g_free (joint.per_kernel);

  {
    NcXcorKernelIntegrand *integrand = nc_xcor_kernel_integrand_new (n_kernels * n_l,
                                                                     _spline_integrand_eval,
                                                                     _spline_integrand_get_range,
                                                                     sid,
                                                                     _spline_integrand_data_free);

    nc_xcor_kernel_integrand_set_get_knots (integrand, _spline_integrand_get_knots);

    return integrand;
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
 * nc_xcor_kernel_integrand_peek_knots:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Peeks the knots @integrand's components are represented on, shared by every
 * component (multipole) it carries, or %NULL when @integrand is not
 * spline-backed.
 *
 * These knots are what makes the outer $k$ integral exactly integrable: each
 * component is a cubic spline in $k$, so on every knot panel the product
 * $k^2 W_i(k) W_j(k)$ entering $C_\ell$ is a polynomial of degree $8$, which a
 * $5$-node Gauss-Legendre rule integrates exactly.
 *
 * Returns: (transfer none) (nullable): the knot vector, or %NULL.
 */
NcmVector *
nc_xcor_kernel_integrand_peek_knots (NcXcorKernelIntegrand *integrand)
{
  if (integrand->get_knots_func == NULL)
    return NULL;

  return integrand->get_knots_func (integrand->data);
}

/**
 * nc_xcor_kernel_integrand_get_component_range:
 * @integrand: a #NcXcorKernelIntegrand
 * @component: component index, in [0, nc_xcor_kernel_integrand_get_len())
 * @k_min: (out): lowest wavenumber at which @component carries signal
 * @k_max: (out): highest wavenumber at which @component carries signal
 *
 * Gets the support of a single component, which for a joint (multi-kernel)
 * integrand is narrower than the shared domain reported by
 * nc_xcor_kernel_integrand_get_range().
 *
 * #NcmSpline deliberately does not range-check on evaluation, so a component
 * evaluated outside its own support returns extrapolation, not a small number:
 * measured on a mid-redshift top-hat bin, the true value at $\tilde k = 2703$
 * is $2\times10^{-5}$ while the spline gives $14.9$. Weighted by $k^2$ in the
 * outer integral that dominates the result. Callers must therefore restrict a
 * pair $(i,j)$ to
 * $[\max(k^{\rm min}_i, k^{\rm min}_j), \min(k^{\rm max}_i, k^{\rm max}_j)]$.
 * Since those bounds are knots of the shared set, the panels remain spline
 * panels and the quadrature remains exact.
 *
 * When @integrand is not spline-backed the shared range is reported for every
 * component.
 *
 * Returns: %TRUE if @component has non-empty support.
 */
gboolean
nc_xcor_kernel_integrand_get_component_range (NcXcorKernelIntegrand *integrand, guint component, gdouble *k_min, gdouble *k_max)
{
  SplineIntegrandData *sid;

  g_assert_cmpuint (component, <, integrand->len);

  if (integrand->get_knots_func == NULL)
  {
    integrand->get_range_func (integrand->data, k_min, k_max);

    return TRUE;
  }

  sid = (SplineIntegrandData *) integrand->data;

  *k_min = sid->supp_k_min[component];
  *k_max = sid->supp_k_max[component];

  return (*k_min) < (*k_max);
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
 *
 * Gets an evaluation function for the kernel at multipole @l.
 * Convenience wrapper around nc_xcor_kernel_get_eval_vectorized() for a single multipole.
 *
 * Returns: (transfer full): the evaluation function for the kernel.
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_get_eval (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  return nc_xcor_kernel_get_eval_vectorized (xclk, cosmo, l, l);
}

/**
 * nc_xcor_kernel_get_eval_vectorized:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @lmin: minimum multipole
 * @lmax: maximum multipole
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
nc_xcor_kernel_get_eval_vectorized (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  return nc_xcor_kernel_get_eval_vectorized_full (xclk, cosmo, lmin, lmax, self->sbi);
}

/**
 * nc_xcor_kernel_get_eval_vectorized_full:
 * @xclk: a #NcXcorKernel
 * @cosmo: a #NcHICosmo
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 * @sbi: (nullable): the #NcmSBesselIntegrator to use, or %NULL for @xclk's own
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
 * Bessel integral is performed.
 *
 * Returns: (transfer full): the kernel integrand over [@lmin, @lmax]
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_get_eval_vectorized_full (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  if ((self->l_limber == 0) || ((self->l_limber > 0) && (lmin >= self->l_limber)))
    return _nc_xcor_kernel_build_limber_integrand (xclk, cosmo, lmin, lmax);
  else
    return _nc_xcor_kernel_build_non_limber_integrand (xclk, cosmo, lmin, lmax, sbi);
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
 * @scaled_abstol: the absolute minimum (must be > 0)
 *
 * Sets the absolute minimum threshold for adaptive midpoint refinement. This parameter
 * helps prevent excessive refinement in cases where the kernel has very low amplitude,
 * by providing a floor below which the refinement will stop regardless of the relative
 * tolerance.
 */
void
nc_xcor_kernel_set_scaled_abstol (NcXcorKernel *xclk, gdouble scaled_abstol)
{
  NcXcorKernelPrivate *self = nc_xcor_kernel_get_instance_private (xclk);

  g_assert_cmpfloat (scaled_abstol, >, 0.0);

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

