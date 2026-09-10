/***************************************************************************
 *            ncm_sbessel_integrator_levin.c
 *
 *  Sat January 25 00:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sbessel_integrator_levin.c
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcmSBesselIntegratorLevin:
 *
 * Levin-Bessel method for spherical Bessel function integration.
 *
 * Uses a Levin-type method for low multipoles and vector cubature for high
 * multipoles.
 *
 * The integral is written in the dimensionless variable $x = k \chi$, so that
 * $\int K(\chi, k) j_\ell(k \chi) \mathrm{d}\chi = \int F(x) j_\ell(x) \mathrm{d}x$ with
 * $F(x) = K(x / k, k) / k$. A single panel set therefore supports every $k$.
 *
 * For low ell values, the contribution of a panel $[a, b]$ is obtained by solving
 * $x^2 w''(x) + 2 x w'(x) + (x^2 - \ell(\ell+1)) w(x) = F(x)$ with boundary conditions
 * $w(a) = w(b) = 0$. Since $j_\ell$ solves the homogeneous equation, the combination
 * $x^2 (w' j_\ell - j_\ell' w)$ has derivative $F j_\ell$, so the panel integral is
 * the boundary term $b^2 w'(b) j_\ell(b) - a^2 w'(a) j_\ell(a)$.
 *
 * In practice the equation is solved for $u = x w$, which removes the
 * first-derivative term and gives $x^2 u'' + (x^2 - \ell(\ell+1)) u = x F(x)$, so
 * the forcing handed to #NcmSBesselOdeSolver is the weighted $x F(x)$. Since $w$
 * vanishes at the endpoints, $u'(a) = a w'(a)$ and $u'(b) = b w'(b)$, and the
 * panel contribution evaluated is $b j_\ell(b) u'(b) - a j_\ell(a) u'(a)$.
 *
 * For high multipoles, vector cubature evaluates the integrand and all requested
 * spherical Bessel functions together.
 *
 * See <a href="../../theory/sbessel_projection.html">UltraLevin: Non-Limber
 * Angular Power Spectra</a> for the derivation, the fixed panel grid in $x$,
 * and the conjugate-point condition on a panel's span.
 *
 * ## Accuracy limit from panel placement
 *
 * Panel edges come from the fixed $x$-knot grid
 * (#NcmSBesselIntegratorLevin:x-knots-min, :x-knots-max, :n-knots), which is what
 * lets one panel set serve every $k$. They therefore fall where that grid says
 * rather than where the integrand would prefer, and adjacent panels can nearly
 * cancel: for a $4\sigma$-truncated Gaussian at $\ell = 50$, $k = 8247$, two
 * panels contribute $-1.71\times10^{-4}$ and $+1.75\times10^{-4}$ against a total
 * of $8.4\times10^{-9}$.
 *
 * The relative error there reaches $7\times10^{-5}$, against $\sim10^{-9}$ for
 * #NcmSBesselIntegratorGL on the same integrand, so it is not the conditioning of
 * the integral. It is not a tolerance either: changing
 * #NcmSBesselIntegratorLevin:cheb-reltol from $10^{-8}$ to $10^{-14}$ has no
 * effect. Scaling the whole knot grid moves the error by four orders in either
 * direction, and no offset is good for every $k$.
 *
 * The error is bounded in absolute terms. It appears only where the integral is
 * $10^{-7}$ to $10^{-9}$ of its own peak over $k$, in the deep oscillatory tail
 * $x \gg \ell$, and the worst absolute error measured is $2\times10^{-11}$ of
 * that peak; near the peak the same scan gives $4\times10^{-11}$ relative. Every
 * consumer in the library reaches this class through #NcXcorKernel, whose
 * #NcXcorKernel:scaled-abstol floors the $k$-spline at $10^{-4}$ of the peak,
 * $10^{-5}$ for #NcXcorSSCSij, and refuses to go below $10^{-6}$, leaving the
 * panel error at least four orders under a floor that is applied anyway.
 *
 * Use #NcmSBesselIntegratorGL where a small value has to be accurate in its own
 * right rather than as part of a larger integral.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "ncm/specfunc/ncm_sbessel_ode_solver.h"
#include "ncm/specfunc/ncm_sf_sbessel.h"
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/algebra/ncm_lapack.h"
#include "ncm/core/ncm_c.h"
#include "ncm/integration/ncm_integral_nd.h"
#include "ncm/core/ncm_dtuple.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_integration.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

/*
 * Per-panel diagnostic record. Off by default and zero cost when off.
 *
 * The integral is assembled as a sum of per-panel boundary terms that cancel
 * heavily, so the attainable relative accuracy is bounded by the cancellation
 * ratio sum|contrib| / |total| times the machine epsilon. That ratio cannot be
 * inferred from the result; recording the individual contributions is the only
 * way to measure it.
 */
typedef struct _NcmSBesselIntegratorLevinPanelRec
{
  gdouble a;
  gdouble b;
  gint ell;
  gdouble contrib;
} NcmSBesselIntegratorLevinPanelRec;

struct _NcmSBesselIntegratorLevin
{
  /*< private >*/
  NcmSBesselIntegrator parent_instance;
  guint max_order;
  gdouble reltol;
  guint cheb_min_order;
  gdouble cheb_reltol;
  NcmSBesselOdeSolver *ode_solver;
  NcmSBesselOdeOperator *ode_operator;
  NcmSFSBesselArray *sba; /* Allocation tracking */
  guint alloc_max_order;
  guint alloc_ell_min;
  guint alloc_ell_max;
  gboolean constructed;
  guint deriv; /* Bessel-derivative order of the running integrate call */
  /* Pre-allocated working arrays */
  GArray *cheb_coeffs;
  GArray *edge_cheb_coeffs;
  GArray *gegen_coeffs;
  GArray *deriv_gegen_coeffs;
  GArray *rhs;
  GArray *values_result;
  gdouble *j_array_a;
  gdouble *j_array_b;
  GArray *endpoints_result;
  gdouble *jl_arr;
  gboolean record_panels; /* Diagnostic panel recording (off by default) */

  /* Per-panel constraint rule. A panel whose oscillation count, measured from the
   * highest turning point in the block, exceeds this uses the tau constraint;
   * every other panel keeps Dirichlet data. Zero disables the rule and the
   * solver-wide setting governs. */
  gdouble tau_constraint_min_osc;
  guint n_constraint_fallbacks;   /* Panels the guard sent back to Dirichlet (diagnostic) */
  guint n_panel_solves;           /* Panel solves since the last constraint reset */
  guint n_tau_solves;             /* Of those, solved and kept with the tau constraint */
  guint n_locked_eligible_solves; /* Dirichlet solves on a rule-eligible panel whose forcing would have passed the order check */
  gdouble tau_constraint_guard;   /* Factor above the forcing's scale at which a tau solve is rejected */
  gdouble tau_constraint_order_fraction;
  gboolean dead_edge_cells;
  GArray *panel_records; /* NcmSBesselIntegratorLevinPanelRec, valid when recording */
  /* Knots-based paneling */
  gdouble x_knots_min;
  gdouble x_knots_max;
  guint n_knots;
  guint ell_cache_max;                        /* Maximum ell for precomputed j_l at knots */
  GArray *base_knots;                         /* Log-spaced knots from the construct properties */
  GArray *knots;                              /* Working grid: base_knots plus the turning-point knot */
  gdouble turning_knot_margin;                /* Turning-point knot at margin * nu_max; 0 disables */
  gdouble cur_turning_knot;                   /* Turning-point knot in the working grid, 0 if none */
  guint alloc_n_knots;                        /* Working knot count the operators were built for */
  NcmMatrix *jl_base;                         /* j_l at the base knots, from the Bessel unit's shared store */
  gdouble *jl_split_row;                      /* j_l at the turning-point knot, this instance's own row */
  gint split_pos;                             /* Index of that knot in the working grid, -1 when there is none */
  GPtrArray *operators;                       /* Operators for each panel between consecutive knots */
  GPtrArray *dirichlet_operators;             /* Dirichlet twin of a panel operator, built on that panel's first fallback */
  GHashTable *edge_operators;                 /* Fixed-cell operators (half-octave ladder) used by moving edge panels */
  NcmSBesselOdeOperator *ode_operator_temp_a; /* Temporary operator for [a, smallest_knot > a] */
  NcmSBesselOdeOperator *ode_operator_temp_b; /* Temporary operator for [largest_knot < b, b] */
  gboolean ode_operator_temp_a_valid;         /* True when temp_a matches the cached panel */
  gboolean ode_operator_temp_b_valid;         /* True when temp_b matches the cached panel */
  gdouble ode_operator_temp_a_a;
  gdouble ode_operator_temp_a_b;
  gdouble ode_operator_temp_b_a;
  gdouble ode_operator_temp_b_b;
  guint ode_operator_temp_a_ell_min;
  guint ode_operator_temp_a_ell_max;
  guint ode_operator_temp_b_ell_min;
  guint ode_operator_temp_b_ell_max;
};

enum
{
  PROP_0,
  PROP_MAX_ORDER,
  PROP_RELTOL,
  PROP_CHEB_MIN_ORDER,
  PROP_CHEB_RELTOL,
  PROP_X_KNOTS_MIN,
  PROP_X_KNOTS_MAX,
  PROP_N_KNOTS,
  PROP_ELL_CACHE_MAX,
  PROP_TAU_CONSTRAINT_MIN_OSC,
  PROP_TAU_CONSTRAINT_GUARD,
  PROP_TAU_CONSTRAINT_ORDER_FRACTION,
  PROP_TURNING_KNOT_MARGIN,
  PROP_DEAD_EDGE_CELLS,
};

static void _ncm_sbessel_integrator_levin_prepare_knots_array (NcmSBesselIntegratorLevin *sbilv);
static void _ncm_sbessel_integrator_levin_prepare_ell_cache (NcmSBesselIntegratorLevin *sbilv);
static gboolean _ncm_sbessel_integrator_levin_rebuild_working_knots (NcmSBesselIntegratorLevin *sbilv, guint ell_max);
static gdouble _ncm_sbessel_integrator_levin_turning_knot (NcmSBesselIntegratorLevin *sbilv, guint ell_max);
static void _ncm_sbessel_integrator_levin_prepare_jl_knots (NcmSBesselIntegratorLevin *sbilv);

/* j_l for one knot of the working grid: a shared base row, or this instance's own row
 * for the inserted turning-point knot. */
static inline const gdouble *
_ncm_sbessel_integrator_levin_jl_row (NcmSBesselIntegratorLevin *sbilv, gint idx)
{
  if (sbilv->split_pos < 0)
    return ncm_matrix_const_ptr (sbilv->jl_base, idx, 0);

  if (idx < sbilv->split_pos)
    return ncm_matrix_const_ptr (sbilv->jl_base, idx, 0);

  if (idx == sbilv->split_pos)
    return sbilv->jl_split_row;

  return ncm_matrix_const_ptr (sbilv->jl_base, idx - 1, 0);
}

static void _ncm_sbessel_operator_unref_null_ok (gpointer op);
static void _ncm_sbessel_integrator_levin_ensure_prepared (NcmSBesselIntegratorLevin *sbilv, guint max_order, guint ell_min, guint ell_max);
static void _ncm_sbessel_integrator_levin_apply_constraint (NcmSBesselIntegratorLevin *sbilv, NcmSBesselOdeOperator *op, gdouble a, gdouble b, guint ell_max);
static gboolean _ncm_sbessel_integrator_levin_tau_constraint_order_ok (NcmSBesselIntegratorLevin *sbilv, NcmSBesselOdeOperator *op, gdouble a, gdouble b, guint ell_max);
static gdouble _ncm_sbessel_integrator_levin_osc (gdouble a, gdouble b, guint ell_max);
static gboolean _ncm_sbessel_integrator_levin_tau_constraint_blew_up (NcmSBesselIntegratorLevin *sbilv, NcmSBesselOdeOperator *op, GArray *cheb, gdouble a, gdouble b, guint ell_min, guint ell_max);
static void _ncm_sbessel_integrator_levin_prepare_knots_operators (NcmSBesselIntegratorLevin *sbilv, guint ell_min, guint ell_max);
static void _ncm_sbessel_integrator_levin_compute_rhs (NcmSBesselIntegratorLevin *sbilv, NcmSpectral *spectral, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, gpointer user_data);
static void _ncm_sbessel_integrator_levin_set_ell_range (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max);
static void _ncm_sbessel_integrator_levin_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, NcmVector *result, gpointer user_data);
static void _ncm_sbessel_integrator_levin_integrate_deriv (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, guint deriv, NcmVector *result, gpointer user_data);

G_DEFINE_TYPE (NcmSBesselIntegratorLevin, ncm_sbessel_integrator_levin, NCM_TYPE_SBESSEL_INTEGRATOR)

static void
ncm_sbessel_integrator_levin_init (NcmSBesselIntegratorLevin *sbilv)
{
  sbilv->max_order                     = 0;
  sbilv->reltol                        = 0.0;
  sbilv->cheb_min_order                = 0;
  sbilv->cheb_reltol                   = 0.0;
  sbilv->tau_constraint_min_osc        = 0.0;
  sbilv->n_constraint_fallbacks        = 0;
  sbilv->n_panel_solves                = 0;
  sbilv->n_tau_solves                  = 0;
  sbilv->n_locked_eligible_solves      = 0;
  sbilv->tau_constraint_guard          = 0.0;
  sbilv->tau_constraint_order_fraction = 0.0;
  sbilv->dead_edge_cells               = TRUE;
  sbilv->ode_solver                    = ncm_sbessel_ode_solver_new ();
  sbilv->ode_operator                  = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, 2, 2);
  sbilv->sba                           = ncm_sf_sbessel_array_new ();
  sbilv->alloc_max_order               = 0;
  sbilv->alloc_ell_min                 = -1;
  sbilv->alloc_ell_max                 = -1;
  sbilv->deriv                         = 0;
  sbilv->cheb_coeffs                   = NULL;
  sbilv->edge_cheb_coeffs              = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sbilv->gegen_coeffs                  = NULL;
  sbilv->deriv_gegen_coeffs            = NULL;
  sbilv->rhs                           = NULL;
  sbilv->values_result                 = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sbilv->j_array_a                     = NULL;
  sbilv->j_array_b                     = NULL;
  sbilv->endpoints_result              = NULL;
  sbilv->jl_arr                        = NULL;
  sbilv->constructed                   = FALSE;
  sbilv->record_panels                 = FALSE;
  sbilv->panel_records                 = g_array_new (FALSE, FALSE, sizeof (NcmSBesselIntegratorLevinPanelRec));
  /* Knots-based paneling */
  sbilv->x_knots_min         = 0.0;
  sbilv->x_knots_max         = 0.0;
  sbilv->n_knots             = 0;
  sbilv->ell_cache_max       = 0;
  sbilv->base_knots          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sbilv->knots               = g_array_new (FALSE, FALSE, sizeof (gdouble));
  sbilv->turning_knot_margin = 0.0;
  sbilv->cur_turning_knot    = -1.0;
  sbilv->alloc_n_knots       = 0;
  sbilv->jl_base             = NULL;
  sbilv->jl_split_row        = NULL;
  sbilv->split_pos           = -1;
  sbilv->operators           = NULL;
  sbilv->dirichlet_operators = NULL;
  sbilv->edge_operators      = g_hash_table_new_full (g_int64_hash, g_int64_equal, g_free,
                                                      (GDestroyNotify) ncm_sbessel_ode_operator_unref);
  sbilv->ode_operator_temp_a       = NULL;
  sbilv->ode_operator_temp_b       = NULL;
  sbilv->ode_operator_temp_a_valid = FALSE;
  sbilv->ode_operator_temp_b_valid = FALSE;
}

static void
_ncm_sbessel_integrator_levin_dispose (GObject *object)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  ncm_sbessel_ode_solver_clear (&sbilv->ode_solver);
  ncm_sbessel_ode_operator_clear (&sbilv->ode_operator);
  ncm_sf_sbessel_array_clear (&sbilv->sba);
  g_clear_pointer (&sbilv->panel_records, g_array_unref);
  g_clear_pointer (&sbilv->cheb_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->edge_cheb_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->gegen_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->deriv_gegen_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->rhs, g_array_unref);
  g_clear_pointer (&sbilv->endpoints_result, g_array_unref);

  /* Clear knots-based paneling resources */
  g_clear_pointer (&sbilv->base_knots, g_array_unref);
  g_clear_pointer (&sbilv->knots, g_array_unref);

  ncm_matrix_clear (&sbilv->jl_base);

  if (sbilv->jl_split_row != NULL)
  {
    g_free (sbilv->jl_split_row);
    sbilv->jl_split_row = NULL;
  }

  if (sbilv->operators != NULL)
  {
    g_ptr_array_unref (sbilv->operators);
    sbilv->operators = NULL;
  }

  if (sbilv->dirichlet_operators != NULL)
  {
    g_ptr_array_unref (sbilv->dirichlet_operators);
    sbilv->dirichlet_operators = NULL;
  }

  ncm_sbessel_ode_operator_clear (&sbilv->ode_operator_temp_a);
  ncm_sbessel_ode_operator_clear (&sbilv->ode_operator_temp_b);
  g_clear_pointer (&sbilv->edge_operators, g_hash_table_unref);

  g_clear_pointer (&sbilv->values_result, g_array_unref);

  if (sbilv->j_array_a != NULL)
  {
    g_free (sbilv->j_array_a);
    sbilv->j_array_a = NULL;
  }

  if (sbilv->j_array_b != NULL)
  {
    g_free (sbilv->j_array_b);
    sbilv->j_array_b = NULL;
  }

  if (sbilv->jl_arr != NULL)
  {
    g_free (sbilv->jl_arr);
    sbilv->jl_arr = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_levin_parent_class)->dispose (object);
}

static void
_ncm_sbessel_integrator_levin_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_levin_parent_class)->finalize (object);
}

static void
_ncm_sbessel_integrator_levin_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_sbessel_integrator_levin_parent_class)->constructed (object);

  {
    NcmSBesselIntegrator *sbi        = NCM_SBESSEL_INTEGRATOR (object);
    NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);
    guint ell_min, ell_max;

    /* Prepare knots array and precompute spherical Bessel functions at construction time
     * since knots parameters and ell_cache_max are CONSTRUCT_ONLY */
    _ncm_sbessel_integrator_levin_prepare_knots_array (sbilv);
    _ncm_sbessel_integrator_levin_prepare_ell_cache (sbilv);
    ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
    _ncm_sbessel_integrator_levin_rebuild_working_knots (sbilv, ell_max);
    _ncm_sbessel_integrator_levin_prepare_jl_knots (sbilv);

    /* Mark as constructed and trigger set_ell_range to initialize operators */
    sbilv->constructed = TRUE;
    ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
    ncm_sbessel_integrator_set_ell_range (sbi, ell_min, ell_max);
  }
}

static void
_ncm_sbessel_integrator_levin_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_LEVIN (object));

  switch (prop_id)
  {
    case PROP_MAX_ORDER:
      ncm_sbessel_integrator_levin_set_max_order (sbilv, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      ncm_sbessel_integrator_levin_set_reltol (sbilv, g_value_get_double (value));
      break;
    case PROP_CHEB_MIN_ORDER:
      ncm_sbessel_integrator_levin_set_cheb_min_order (sbilv, g_value_get_uint (value));
      break;
    case PROP_CHEB_RELTOL:
      ncm_sbessel_integrator_levin_set_cheb_reltol (sbilv, g_value_get_double (value));
      break;
    case PROP_TAU_CONSTRAINT_MIN_OSC:
      ncm_sbessel_integrator_levin_set_tau_constraint_min_osc (sbilv, g_value_get_double (value));
      break;
    case PROP_TAU_CONSTRAINT_GUARD:
      ncm_sbessel_integrator_levin_set_tau_constraint_guard (sbilv, g_value_get_double (value));
      break;
    case PROP_TAU_CONSTRAINT_ORDER_FRACTION:
      ncm_sbessel_integrator_levin_set_tau_constraint_order_fraction (sbilv, g_value_get_double (value));
      break;
    case PROP_TURNING_KNOT_MARGIN:
      ncm_sbessel_integrator_levin_set_turning_knot_margin (sbilv, g_value_get_double (value));
      break;
    case PROP_DEAD_EDGE_CELLS:
      ncm_sbessel_integrator_levin_set_dead_edge_cells (sbilv, g_value_get_boolean (value));
      break;
    case PROP_X_KNOTS_MIN:
      sbilv->x_knots_min = g_value_get_double (value);
      break;
    case PROP_X_KNOTS_MAX:
      sbilv->x_knots_max = g_value_get_double (value);
      break;
    case PROP_N_KNOTS:
      sbilv->n_knots = g_value_get_uint (value);
      break;
    case PROP_ELL_CACHE_MAX:
      sbilv->ell_cache_max = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_integrator_levin_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (object);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR_LEVIN (object));

  switch (prop_id)
  {
    case PROP_MAX_ORDER:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_max_order (sbilv));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_reltol (sbilv));
      break;
    case PROP_CHEB_MIN_ORDER:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_cheb_min_order (sbilv));
      break;
    case PROP_CHEB_RELTOL:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_cheb_reltol (sbilv));
      break;
    case PROP_TAU_CONSTRAINT_MIN_OSC:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_tau_constraint_min_osc (sbilv));
      break;
    case PROP_TAU_CONSTRAINT_GUARD:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_tau_constraint_guard (sbilv));
      break;
    case PROP_TAU_CONSTRAINT_ORDER_FRACTION:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_tau_constraint_order_fraction (sbilv));
      break;
    case PROP_TURNING_KNOT_MARGIN:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_turning_knot_margin (sbilv));
      break;
    case PROP_DEAD_EDGE_CELLS:
      g_value_set_boolean (value, ncm_sbessel_integrator_levin_get_dead_edge_cells (sbilv));
      break;
    case PROP_X_KNOTS_MIN:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_x_knots_min (sbilv));
      break;
    case PROP_X_KNOTS_MAX:
      g_value_set_double (value, ncm_sbessel_integrator_levin_get_x_knots_max (sbilv));
      break;
    case PROP_N_KNOTS:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_n_knots (sbilv));
      break;
    case PROP_ELL_CACHE_MAX:
      g_value_set_uint (value, ncm_sbessel_integrator_levin_get_ell_cache_max (sbilv));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_sbessel_integrator_levin_class_init (NcmSBesselIntegratorLevinClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmSBesselIntegratorClass *parent_class = NCM_SBESSEL_INTEGRATOR_CLASS (klass);

  object_class->set_property = &_ncm_sbessel_integrator_levin_set_property;
  object_class->get_property = &_ncm_sbessel_integrator_levin_get_property;
  object_class->constructed  = &_ncm_sbessel_integrator_levin_constructed;
  object_class->dispose      = &_ncm_sbessel_integrator_levin_dispose;
  object_class->finalize     = &_ncm_sbessel_integrator_levin_finalize;

  /**
   * NcmSBesselIntegratorLevin:max-order:
   *
   * Maximum order of Chebyshev decomposition for the Levin method. Higher order may
   * give better accuracy but is more expensive. Default is 2^14.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_ORDER,
                                   g_param_spec_uint ("max-order",
                                                      NULL,
                                                      "Maximum Chebyshev order",
                                                      2, G_MAXUINT, 1 << 14,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:reltol:
   *
   * Relative tolerance of the ODE solve. Bounds the result together with
   * #NcmSBesselIntegratorLevin:cheb-reltol; the looser of the two wins.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "ODE solve relative tolerance",
                                                        0.0, 1.0, 1.0e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:cheb-min-order:
   *
   * Minimum order of Chebyshev decomposition used when computing the RHS for the
   * Levin method.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHEB_MIN_ORDER,
                                   g_param_spec_uint ("cheb-min-order",
                                                      NULL,
                                                      "Minimum Chebyshev order for RHS",
                                                      1, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:cheb-reltol:
   *
   * Relative tolerance of the integrand's Chebyshev fit, which forms the RHS.
   * Bounds the result together with #NcmSBesselIntegratorLevin:reltol; the
   * looser of the two wins.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHEB_RELTOL,
                                   g_param_spec_double ("cheb-reltol",
                                                        NULL,
                                                        "Integrand Chebyshev fit relative tolerance",
                                                        0.0, 1.0, 1.0e-8,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:tau-constraint-min-osc:
   *
   * A panel solve uses the tau constraint when the panel's oscillation count,
   * $\frac{2}{\pi}\,(b - \max(a, \nu_{\max}))$ with $\nu_{\max}$ the turning point of
   * the highest multipole in the block, exceeds this value; otherwise the panel takes
   * the pinned constraint, with its pins at the peak of the homogeneous spectrum, and
   * Dirichlet data only when the panel is too short to leave room for them. Zero
   * disables the rule, leaving every panel on the Dirichlet data it is created with.
   * The tau constraint is valid only where the homogeneous solutions are
   * unrepresentable at the working order, which is what the count measures. That condition is $N_{\min} \gtrsim 1.5\,n_F$, from the floor
   * $1.1\,n_F$ against the order check's $0.75\,N_{\min}$, so the count a kernel needs
   * follows its forcing order: 15 to 30 for the analytic windows ($n_F$ of 9 to 18) and
   * 60 to 200 for tabulated kernels, whose sampling noise gives them an algebraic
   * coefficient tail ($n_F$ of 40 to 135). The default of 50 is the lowest value the
   * gaussian truth table still passes: at 45 one panel at $\ell = 228$ takes the tau
   * constraint and loses the k reach the table asks for, while the Arb window table is
   * unchanged from 200 down to 5. Lower is cheaper, 126.3 MiB of resident memory at 50
   * against 136.8 at 200 on the lensing blocks, since a tau operator stores
   * $1.1\,n_F$ columns of rotations against a Dirichlet one's $N_{\min}$.
   */
  g_object_class_install_property (object_class,
                                   PROP_TAU_CONSTRAINT_MIN_OSC,
                                   g_param_spec_double ("tau-constraint-min-osc",
                                                        NULL,
                                                        "Panel oscillation count above which the tau constraint is used (0 disables)",
                                                        0.0, G_MAXDOUBLE, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_MIN_OSC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:tau-constraint-guard:
   *
   * A tau-constraint solve is rejected, and the panel redone with Dirichlet data, when
   * its largest Chebyshev coefficient exceeds this factor times the bound
   * $b\,\sum_k|c_k| / (a^2 - \nu_\ell^2)$ on the smooth member, judged per multipole
   * of the block. A valid solve sits at a ratio of order one and a contaminated one at
   * $10^{11}$ and above; the default of 10 is the smallest value tested that catches
   * every marginal case in the Arb truth tables with no false rejection on valid
   * panels. Larger values admit mild contamination that costs digits silently.
   */
  g_object_class_install_property (object_class,
                                   PROP_TAU_CONSTRAINT_GUARD,
                                   g_param_spec_double ("tau-constraint-guard",
                                                        NULL,
                                                        "Coefficient-norm factor above the forcing scale that rejects a tau-constraint solve",
                                                        1.0, G_MAXDOUBLE, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_GUARD,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:tau-constraint-order-fraction:
   *
   * Largest working order of a tau-constraint solve as a fraction of the panel's
   * oscillation count $N_\min = 2\,\mathrm{span}/\pi$. A forcing whose fit needs more
   * coefficients than this is solved with Dirichlet data instead.
   */
  g_object_class_install_property (object_class,
                                   PROP_TAU_CONSTRAINT_ORDER_FRACTION,
                                   g_param_spec_double ("tau-constraint-order-fraction",
                                                        NULL,
                                                        "Tau-constraint working order as a fraction of N_min",
                                                        0.0, 1.0, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_ORDER_FRACTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_STRINGS));

  /**
   * NcmSBesselIntegratorLevin:turning-knot-margin:
   *
   * One knot is added to the panel grid of each block, at this multiple of the turning
   * point $\nu_{\max} = \sqrt{\ell_{\max}(\ell_{\max}+1)}$ of the block's highest
   * multipole, whenever that lands strictly inside the grid. Zero disables it.
   *
   * The bound on the smooth member that the constraint guard tests against carries
   * $\min|x^2 - \nu^2|$ in its denominator, so on a panel straddling the turning point no
   * finite bound exists and the tau constraint cannot be certified there. When the kernel's
   * reach in $x$ ends before the next base knot, that straddling panel is the whole
   * oscillatory region and the tau constraint is then unusable at every $k$: measured on the
   * analytic windows, zero tau solves for $\ell \gtrsim 500$ on a decade-aligned grid,
   * against 1.5 to 5.2 times faster with the knot in place and no fallbacks left at any
   * multipole.
   *
   * The default of 1.05 clears the turning point by five percent, enough to make the bound
   * finite without moving the knot far into the oscillatory region.
   */
  g_object_class_install_property (object_class,
                                   PROP_TURNING_KNOT_MARGIN,
                                   g_param_spec_double ("turning-knot-margin",
                                                        NULL,
                                                        "Turning-point knot at this multiple of nu_max (0 disables)",
                                                        0.0, G_MAXDOUBLE, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TURNING_KNOT_MARGIN,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:dead-edge-cells:
   *
   * Whether a moving end piece whose forcing has died at the junction with the
   * rest of its cell is solved on the cached cell operator, with the true
   * forcing inside the piece and zero outside. Saves the per-$k$ factorization of
   * the moving operator at the cost of keeping those cells' factorizations.
   */
  g_object_class_install_property (object_class,
                                   PROP_DEAD_EDGE_CELLS,
                                   g_param_spec_boolean ("dead-edge-cells",
                                                         NULL,
                                                         "Solve dead-junction end pieces on cached cells",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_STRINGS));

  /**
   * NcmSBesselIntegratorLevin:x-knots-min:
   *
   * Minimum value for knots in log-spaced grid. Set to 0 to disable knots-based
   * paneling. This property can only be set during construction.
   */
  g_object_class_install_property (object_class,
                                   PROP_X_KNOTS_MIN,
                                   g_param_spec_double ("x-knots-min",
                                                        NULL,
                                                        "Minimum knot value",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:x-knots-max:
   *
   * Maximum value for knots in log-spaced grid. Set to 0 to disable knots-based
   * paneling. This property can only be set during construction.
   */
  g_object_class_install_property (object_class,
                                   PROP_X_KNOTS_MAX,
                                   g_param_spec_double ("x-knots-max",
                                                        NULL,
                                                        "Maximum knot value",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:n-knots:
   *
   * Number of knots in the log-spaced grid. The knots will be equally spaced in log
   * space between x-knots-min and x-knots-max. Set to 0 to disable knots-based
   * paneling. This property can only be set during construction.
   *
   * This is the base grid, shared by every multipole block. What a block is solved on is
   * that grid plus at most the one knot of
   * #NcmSBesselIntegratorLevin:turning-knot-margin.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_KNOTS,
                                   g_param_spec_uint ("n-knots",
                                                      NULL,
                                                      "Number of knots",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegratorLevin:ell-cache-max:
   *
   * Maximum ell value for precomputed spherical Bessel functions at knots. The
   * integrator will precompute j_ell(knot) for all knots and all ell from 0 to
   * ell-cache-max. This enables fast lookup during integration when the requested
   * ell values are within the cached range. For ell values beyond ell-cache-max,
   * the integrator will compute spherical Bessel functions on-the-fly. This property
   * can only be set during construction.
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL_CACHE_MAX,
                                   g_param_spec_uint ("ell-cache-max",
                                                      NULL,
                                                      "Maximum ell for cache",
                                                      0, G_MAXUINT, 500,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->set_ell_range   = &_ncm_sbessel_integrator_levin_set_ell_range;
  parent_class->integrate       = &_ncm_sbessel_integrator_levin_integrate;
  parent_class->integrate_deriv = &_ncm_sbessel_integrator_levin_integrate_deriv;
}

static void
_ncm_sbessel_integrator_levin_ensure_prepared (NcmSBesselIntegratorLevin *sbilv, guint max_order, guint ell_min, guint ell_max)
{
  if (sbilv->alloc_max_order == max_order)
    return;

  /* Free existing allocations */
  g_clear_pointer (&sbilv->cheb_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->gegen_coeffs, g_array_unref);
  g_clear_pointer (&sbilv->rhs, g_array_unref);
  g_clear_pointer (&sbilv->endpoints_result, g_array_unref);

  /* Allocate arrays for spectral coefficients */
  sbilv->cheb_coeffs  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), max_order);
  sbilv->gegen_coeffs = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), max_order);
  sbilv->rhs          = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), max_order + 2);

  /* Allocate result matrix for batched endpoint computation (max block size is 8) */
  sbilv->endpoints_result = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 8 * 3);

  sbilv->alloc_max_order = max_order;
}

/**
 * _ncm_sbessel_integrator_levin_prepare_knots_array:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Creates the log-spaced knots array from x_knots_min to x_knots_max.
 * This only needs to be done when knot parameters change, not when ell range changes.
 */
static void
_ncm_sbessel_integrator_levin_prepare_knots_array (NcmSBesselIntegratorLevin *sbilv)
{
  if (sbilv->n_knots < 2)
    return;

  g_assert_cmpfloat (sbilv->x_knots_min, >, 0.0);
  g_assert_cmpfloat (sbilv->x_knots_max, >, sbilv->x_knots_min);
  g_assert_cmpuint (sbilv->n_knots, <, G_MAXUINT / 2); /* Prevent overflow in log spacing calculation */
  g_array_set_size (sbilv->base_knots, sbilv->n_knots);

  {
    const gdouble ln_x_min = log (sbilv->x_knots_min);
    const gdouble ln_x_max = log (sbilv->x_knots_max);
    const gdouble L        = ln_x_max - ln_x_min;
    const gdouble dL       = L / (sbilv->n_knots - 1.0);
    const gdouble expm1_dL = expm1 (dL);
    gdouble x0             = exp (ln_x_min);
    guint i;

    g_array_index (sbilv->base_knots, gdouble, 0) = x0;

    for (i = 1; i < sbilv->n_knots; i++)
    {
      const gdouble dy = x0 * expm1_dL;
      const gdouble x  = x0 + dy;

      g_array_index (sbilv->base_knots, gdouble, i) = x;
      x0                                            = x;
    }
  }
}

/*
 * Where to split the panel that straddles the turning point, if splitting it pays.
 *
 * Only the panel holding nu_max is blocked by the guard: its left edge is at or below the
 * turning point, so min|x^2 - nu^2| vanishes inside it and the bound on the smooth member
 * is infinite. Every panel above starts clear of nu and is certifiable already. So the
 * knot has to land inside that one panel -- putting it at margin * nu_max unclamped can
 * overshoot into the next panel, splitting a healthy one and leaving a sliver behind.
 *
 * Splitting pays only when the piece above the split would itself qualify under the
 * per-panel rule, so the same oscillation count that gates the rule gates the insertion.
 * That is what keeps low multipoles untouched: there the straddling panel is narrow in
 * oscillations and Dirichlet either way, and an extra panel would be one more forcing fit
 * for nothing.
 *
 * Returns: the split point, or 0.0 when no knot should be added.
 */
static gdouble
_ncm_sbessel_integrator_levin_turning_knot (NcmSBesselIntegratorLevin *sbilv, guint ell_max)
{
  const guint n_base   = sbilv->base_knots->len;
  const gdouble nu_max = sqrt (ell_max * (ell_max + 1.0));
  gdouble x_lo, x_hi, x_split;
  guint j;

  if ((sbilv->turning_knot_margin <= 0.0) || (sbilv->tau_constraint_min_osc <= 0.0) || (n_base < 2))
    return 0.0;

  /* The panel [x_lo, x_hi] holding nu_max. */
  for (j = 0; j + 1 < n_base; j++)
  {
    x_lo = g_array_index (sbilv->base_knots, gdouble, j);
    x_hi = g_array_index (sbilv->base_knots, gdouble, j + 1);

    if ((x_lo <= nu_max) && (nu_max < x_hi))
      break;
  }

  if (j + 1 >= n_base)
    return 0.0;

  x_split = sbilv->turning_knot_margin * nu_max;

  /* Clamped strictly inside the straddling panel, and not so close to either edge that
   * the split leaves a panel of negligible width. */
  if (x_split >= x_hi)
    x_split = sqrt (nu_max * x_hi);

  if ((x_split <= nu_max) ||
      (x_split - nu_max <= 1.0e-3 * nu_max) ||
      (x_hi - x_split <= 1.0e-3 * x_hi))
    return 0.0;

  /* Does the piece above the split earn its own panel? */
  if (2.0 * (x_hi - x_split) / M_PI <= sbilv->tau_constraint_min_osc)
    return 0.0;

  return x_split;
}

/*
 * Working grid for one block: the base grid plus, when it falls strictly inside it, one
 * knot just above the turning point of the block's highest multipole.
 *
 * Without it the panel holding the turning point straddles it, so no finite bound on the
 * smooth member exists there -- the bound carries min|x^2 - nu^2| in its denominator --
 * and the guard of _ncm_sbessel_integrator_levin_tau_constraint_blew_up() has to refuse the
 * tau constraint on that panel. When the kernel's reach in x stops before the next base knot
 * that panel is the whole oscillatory region, and the tau constraint is then unusable at
 * every k: measured on the analytic windows, zero tau solves for ell >= 1000 on a
 * decade-aligned grid. One knot lifts the left edge clear of the turning point and the
 * oscillatory region gets a panel the guard can certify.
 *
 * Returns: whether the working knot count changed.
 */
static gboolean
_ncm_sbessel_integrator_levin_rebuild_working_knots (NcmSBesselIntegratorLevin *sbilv, guint ell_max)
{
  const guint n_base    = sbilv->base_knots->len;
  const gdouble x_split = (n_base < 2) ? 0.0 : _ncm_sbessel_integrator_levin_turning_knot (sbilv, ell_max);
  guint i;

  /* The grid is settled by the base knots and this one value, so nothing changes unless
   * the value does. Comparing lengths is not enough: two blocks can each insert a knot,
   * at different x, and leave the count alone. */
  if ((sbilv->knots->len > 0) && (x_split == sbilv->cur_turning_knot))
    return FALSE;

  sbilv->cur_turning_knot = x_split;
  g_array_set_size (sbilv->knots, 0);

  if (n_base < 2)
    return TRUE;

  {
    if (x_split > 0.0)
    {
      gboolean inserted = FALSE;

      for (i = 0; i < n_base; i++)
      {
        const gdouble x = g_array_index (sbilv->base_knots, gdouble, i);

        if (!inserted && (x > x_split))
        {
          g_array_append_val (sbilv->knots, x_split);
          inserted = TRUE;
        }

        g_array_append_val (sbilv->knots, x);
      }

      return TRUE;
    }
  }

  for (i = 0; i < n_base; i++)
  {
    const gdouble x = g_array_index (sbilv->base_knots, gdouble, i);

    g_array_append_val (sbilv->knots, x);
  }

  return TRUE;
}

/**
 * _ncm_sbessel_integrator_levin_prepare_ell_cache:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * Precomputes spherical Bessel functions j_ell at all knots for ell from 0 to ell_cache_max.
 * This is only needed when ell_cache_max or knot parameters change, not when ell range changes.
 */
static void
_ncm_sbessel_integrator_levin_prepare_ell_cache (NcmSBesselIntegratorLevin *sbilv)
{
  const guint n_ell = sbilv->ell_cache_max + 1;
  guint i;

  g_assert_cmpuint (sbilv->ell_cache_max, <, G_MAXUINT / 2);

  /* Allocate arrays for spherical Bessel functions */
  sbilv->j_array_a = g_new0 (gdouble, n_ell);
  sbilv->j_array_b = g_new0 (gdouble, n_ell);
  sbilv->jl_arr    = g_new (gdouble, n_ell);

  if (sbilv->base_knots->len == 0)
    return;

  /* The base knots never move, so their rows come from the shared store of the Bessel
   * unit; only the turning-point knot's row belongs to this instance. */
  sbilv->jl_base      = ncm_sf_sbessel_array_ref_table (sbilv->sba, sbilv->base_knots, sbilv->ell_cache_max);
  sbilv->jl_split_row = g_new (gdouble, n_ell);

  (void) i;
}

/*
 * Fills j_ell at every knot of the working grid, for ell up to ell_cache_max. Rebuilt
 * whenever the working grid changes, which the turning-point knot makes per block.
 */
static void
_ncm_sbessel_integrator_levin_prepare_jl_knots (NcmSBesselIntegratorLevin *sbilv)
{
  const guint n_base = sbilv->base_knots->len;
  const guint n_work = sbilv->knots->len;
  guint pos;

  if ((n_work == 0) || (sbilv->jl_base == NULL))
    return;

  if (n_work == n_base)
  {
    sbilv->split_pos = -1;

    return;
  }

  g_assert_cmpuint (n_work, ==, n_base + 1);

  /* The working grid is the base grid with one knot inserted; find where, and evaluate
   * that row alone. Every other row is a shared base row. */
  for (pos = 0; pos < n_base; pos++)
  {
    if (g_array_index (sbilv->knots, gdouble, pos) != g_array_index (sbilv->base_knots, gdouble, pos))
      break;
  }

  sbilv->split_pos = (gint) pos;
  ncm_sf_sbessel_array_eval (sbilv->sba, sbilv->ell_cache_max,
                             g_array_index (sbilv->knots, gdouble, pos),
                             sbilv->jl_split_row);
}

/**
 * _ncm_sbessel_integrator_levin_prepare_knots_operators:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 *
 * Prepares ODE operators for the knots-based paneling system:
 *
 * - Pre-allocated ODE operators for each panel between consecutive knots
 * - Two temporary operators for edge panels [a, smallest_knot > a] and [largest_knot < b, b]
 *
 * Uses ncm_sbessel_ode_solver_reconfigure_operator() to update operators in place when the
 * multipole range or the tolerance changes.
 */
static void
_ncm_sbessel_integrator_levin_prepare_knots_operators (NcmSBesselIntegratorLevin *sbilv, guint ell_min, guint ell_max)
{
  if (sbilv->n_knots < 2)
  {
    return;
  }
  else
  {
    const gboolean need_create = (sbilv->operators == NULL) || (sbilv->alloc_n_knots != sbilv->knots->len);

    /* Reconfiguring reads the working grid, so a moved knot is covered by need_reset;
     * only a changed count needs the array rebuilt. */
    const gboolean need_reset = !need_create && ((sbilv->alloc_ell_min != ell_min) || (sbilv->alloc_ell_max != ell_max));
    guint i;

    if (!need_create && !need_reset)
      return;

    if (need_create)
    {
      g_clear_pointer (&sbilv->operators, g_ptr_array_unref);
      sbilv->operators = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_sbessel_ode_operator_unref);
      g_clear_pointer (&sbilv->dirichlet_operators, g_ptr_array_unref);
      sbilv->dirichlet_operators = g_ptr_array_new_with_free_func (_ncm_sbessel_operator_unref_null_ok);
      g_ptr_array_set_size (sbilv->dirichlet_operators, sbilv->knots->len - 1);
      g_hash_table_remove_all (sbilv->edge_operators);

      for (i = 0; i < sbilv->knots->len - 1; i++)
      {
        const gdouble x_a = g_array_index (sbilv->knots, gdouble, i);
        const gdouble x_b = g_array_index (sbilv->knots, gdouble, i + 1);
        NcmSBesselOdeOperator *op;

        op = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, x_a, x_b, ell_min, ell_max);
        _ncm_sbessel_integrator_levin_apply_constraint (sbilv, op, x_a, x_b, ell_max);
        g_ptr_array_add (sbilv->operators, op);
      }

      ncm_sbessel_ode_operator_clear (&sbilv->ode_operator_temp_a);
      ncm_sbessel_ode_operator_clear (&sbilv->ode_operator_temp_b);

      sbilv->ode_operator_temp_a       = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, ell_min, ell_max);
      sbilv->ode_operator_temp_b       = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, 0.0, 1.0, ell_min, ell_max);
      sbilv->ode_operator_temp_a_valid = FALSE;
      sbilv->ode_operator_temp_b_valid = FALSE;
    }
    else if (need_reset)
    {
      g_hash_table_remove_all (sbilv->edge_operators);

      for (i = 0; i < sbilv->operators->len; i++)
      {
        NcmSBesselOdeOperator *op = g_ptr_array_index (sbilv->operators, i);
        const gdouble x_a         = g_array_index (sbilv->knots, gdouble, i);
        const gdouble x_b         = g_array_index (sbilv->knots, gdouble, i + 1);

        ncm_sbessel_ode_solver_reconfigure_operator (sbilv->ode_solver, op, x_a, x_b, ell_min, ell_max);
        _ncm_sbessel_integrator_levin_apply_constraint (sbilv, op, x_a, x_b, ell_max);
      }

      /* Reset temporary operators */
      ncm_sbessel_ode_solver_reconfigure_operator (sbilv->ode_solver, sbilv->ode_operator_temp_a, 0.0, 1.0, ell_min, ell_max);
      ncm_sbessel_ode_solver_reconfigure_operator (sbilv->ode_solver, sbilv->ode_operator_temp_b, 0.0, 1.0, ell_min, ell_max);
      sbilv->ode_operator_temp_a_valid = FALSE;
      sbilv->ode_operator_temp_b_valid = FALSE;

      /* The twins hold the old multipole range. Dropping them is cheaper than
       * reconfiguring twins no fallback may ask for again. */
      g_ptr_array_set_size (sbilv->dirichlet_operators, 0);
      g_ptr_array_set_size (sbilv->dirichlet_operators, sbilv->knots->len - 1);
    }

    sbilv->alloc_ell_min = ell_min;
    sbilv->alloc_ell_max = ell_max;
    sbilv->alloc_n_knots = sbilv->knots->len;
  }
}

/* Wrapper function that transforms K(chi, k) -> f(x) = K(x/k, k)/k */
typedef struct _NcmSBesselIntegratorLevinWrapper
{
  NcmSBesselIntegratorF K;
  gdouble k;
  gpointer user_data;
} NcmSBesselIntegratorLevinWrapper;

static gdouble
_ncm_sbessel_integrator_levin_wrapper_func (gpointer data, gdouble x)
{
  NcmSBesselIntegratorLevinWrapper *wrapper = (NcmSBesselIntegratorLevinWrapper *) data;

  const gdouble chi   = x / wrapper->k;
  const gdouble K_val = wrapper->K (wrapper->user_data, chi, wrapper->k);

  return chi * K_val;
}

/* Same change of variable without the chi weight: samples F(x) = K(x/k, k)/k. */
static gdouble
_ncm_sbessel_integrator_levin_wrapper_func_plain (gpointer data, gdouble x)
{
  NcmSBesselIntegratorLevinWrapper *wrapper = (NcmSBesselIntegratorLevinWrapper *) data;

  const gdouble chi   = x / wrapper->k;
  const gdouble K_val = wrapper->K (wrapper->user_data, chi, wrapper->k);

  return K_val / wrapper->k;
}

static NcmSpectralF
_ncm_sbessel_integrator_levin_peek_wrapper_func (NcmSBesselIntegratorLevin *sbilv)
{
  if (sbilv->deriv == 0)
    return &_ncm_sbessel_integrator_levin_wrapper_func;

  return &_ncm_sbessel_integrator_levin_wrapper_func_plain;
}

static void
_ncm_sbessel_integrator_levin_build_rhs_from_gegen (NcmSBesselIntegratorLevin *sbilv)
{
  /* Two boundary rows make three coefficients the solver's safe minimum. */
  if (sbilv->gegen_coeffs->len < 3)
  {
    const guint old_len = sbilv->gegen_coeffs->len;
    guint i;

    g_array_set_size (sbilv->gegen_coeffs, 3);

    for (i = old_len; i < 3; i++)
      g_array_index (sbilv->gegen_coeffs, gdouble, i) = 0.0;
  }

  g_array_set_size (sbilv->rhs, sbilv->gegen_coeffs->len + 2);
  {
    gdouble *rhs_data                = (gdouble *) sbilv->rhs->data;
    const gdouble *gegen_coeffs_data = (gdouble *) sbilv->gegen_coeffs->data;

    rhs_data[0] = 0.0;
    rhs_data[1] = 0.0;
    memcpy (&rhs_data[2], gegen_coeffs_data, sbilv->gegen_coeffs->len * sizeof (gdouble));
  }
}

/*
 * Builds the ODE forcing from the Chebyshev fit in cheb_coeffs, expressed on
 * [a, b]. For deriv == 0 the fit is of x F(x) and the forcing is its C^(2)
 * representation. For deriv > 0 the fit is of F(x) itself and the forcing is
 * x F'(x) or x F''(x): the derivative is read off in the C^(2) basis (a
 * banded, respectively diagonal, map of the same coefficients), scaled by the
 * chain-rule factor of the interval, and multiplied by x(t) = mid + half t.
 * The deriv == 1 forcing carries a minus sign, matching the single
 * integration by parts int F j' = [F j] - int F' j.
 */
static void
_ncm_sbessel_integrator_levin_build_rhs (NcmSBesselIntegratorLevin *sbilv, gdouble a, gdouble b)
{
  switch (sbilv->deriv)
  {
    case 0:
      ncm_spectral_chebT_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->gegen_coeffs);
      break;
    case 1:
    case 2:
    {
      const gdouble half  = 0.5 * (b - a);
      const gdouble mid   = 0.5 * (a + b);
      const gdouble dt_dy = 2.0 / (b - a);
      const gdouble scale = (sbilv->deriv == 1) ? -dt_dy : dt_dy * dt_dy;
      guint i;

      if (sbilv->deriv == 1)
        ncm_spectral_chebT_deriv_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->deriv_gegen_coeffs);
      else
        ncm_spectral_chebT_deriv2_to_gegenbauer_alpha2 (sbilv->cheb_coeffs, &sbilv->deriv_gegen_coeffs);

      for (i = 0; i < sbilv->deriv_gegen_coeffs->len; i++)
        g_array_index (sbilv->deriv_gegen_coeffs, gdouble, i) *= scale;

      ncm_spectral_gegenbauer_alpha2_xmul (sbilv->deriv_gegen_coeffs, half, mid, &sbilv->gegen_coeffs);
      break;
    }
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }

  _ncm_sbessel_integrator_levin_build_rhs_from_gegen (sbilv);
}

/**
 * _ncm_sbessel_integrator_levin_compute_rhs:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @spectral: spectral methods object
 * @F: integrand function K(chi, k)
 * @a: lower integration bound in x-space
 * @b: upper integration bound in x-space
 * @k: wave number parameter
 * @user_data: user data for integrand
 *
 * Computes the RHS for the Levin ODE by:
 *
 * 1. Computing Chebyshev coefficients for f(x) = K(x/k, k)/k
 * 2. Converting to Gegenbauer C^(2) basis
 * 3. Setting up RHS with homogeneous boundary conditions
 */
static void
_ncm_sbessel_integrator_levin_compute_rhs (NcmSBesselIntegratorLevin *sbilv,
                                           NcmSpectral *spectral,
                                           NcmSBesselIntegratorF F,
                                           gdouble a, gdouble b,
                                           gdouble k,
                                           gpointer user_data)
{
  NcmSBesselIntegratorLevinWrapper wrapper = {F, k, user_data};

  ncm_spectral_compute_chebyshev_coeffs_adaptive_full (spectral, _ncm_sbessel_integrator_levin_peek_wrapper_func (sbilv),
                                                       a, b, sbilv->cheb_min_order, sbilv->cheb_reltol, 0.0,
                                                       &sbilv->cheb_coeffs, &wrapper);

  _ncm_sbessel_integrator_levin_build_rhs (sbilv, a, b);
}

/**
 * _ncm_sbessel_integrator_levin_get_panel_resources:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @a_p_idx: knot index for left endpoint (-1 if not on a knot)
 * @b_p_idx: knot index for right endpoint (-1 if not on a knot)
 * @a_p: lower bound of panel
 * @b_p: upper bound of panel
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 * @j_a_p_out: (out): pointer to j_ell array at a_p
 * @j_b_p_out: (out): pointer to j_ell array at b_p
 *
 * Gets panel resources: j_ell arrays and appropriate ODE operator.
 * If a_p_idx == -1, computes j_a_p and configures temp_a if the panel changed.
 * If b_p_idx == -1, computes j_b_p and configures temp_b if the panel changed.
 * Otherwise, uses cached j_ell values and operators[a_p_idx] without reset.
 *
 * Returns: (transfer none): the appropriate ODE operator for this panel
 */
static NcmSBesselOdeOperator *
_ncm_sbessel_integrator_levin_get_panel_resources (NcmSBesselIntegratorLevin *sbilv,
                                                   gint a_p_idx, gint b_p_idx,
                                                   gdouble a_p, gdouble b_p,
                                                   guint ell_min, guint ell_max,
                                                   const gdouble **j_a_p_out,
                                                   const gdouble **j_b_p_out)
{
  NcmSBesselOdeOperator *op;

  /* Get j_a_p: compute at a moving endpoint or use a knot cache. */
  if (a_p_idx < 0)
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, a_p, sbilv->j_array_a);
    *j_a_p_out = sbilv->j_array_a;
  }
  else
  {
    *j_a_p_out = _ncm_sbessel_integrator_levin_jl_row (sbilv, a_p_idx);
  }

  /* Get j_b_p: compute at a moving endpoint or use a knot cache. */
  if (b_p_idx < 0)
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, b_p, sbilv->j_array_b);
    *j_b_p_out = sbilv->j_array_b;
  }
  else
  {
    *j_b_p_out = _ncm_sbessel_integrator_levin_jl_row (sbilv, b_p_idx);
  }

  /* Reuse an edge factorization until its interval or ell range changes. */
  if (a_p_idx < 0)
  {
    op = sbilv->ode_operator_temp_a;

    if (!sbilv->ode_operator_temp_a_valid ||
        (sbilv->ode_operator_temp_a_a != a_p) ||
        (sbilv->ode_operator_temp_a_b != b_p) ||
        (sbilv->ode_operator_temp_a_ell_min != ell_min) ||
        (sbilv->ode_operator_temp_a_ell_max != ell_max))
    {
      ncm_sbessel_ode_solver_reconfigure_operator (sbilv->ode_solver, op, a_p, b_p, ell_min, ell_max);
      _ncm_sbessel_integrator_levin_apply_constraint (sbilv, op, a_p, b_p, ell_max);
      sbilv->ode_operator_temp_a_a       = a_p;
      sbilv->ode_operator_temp_a_b       = b_p;
      sbilv->ode_operator_temp_a_ell_min = ell_min;
      sbilv->ode_operator_temp_a_ell_max = ell_max;
      sbilv->ode_operator_temp_a_valid   = TRUE;
    }
  }
  else if (b_p_idx < 0)
  {
    op = sbilv->ode_operator_temp_b;

    if (!sbilv->ode_operator_temp_b_valid ||
        (sbilv->ode_operator_temp_b_a != a_p) ||
        (sbilv->ode_operator_temp_b_b != b_p) ||
        (sbilv->ode_operator_temp_b_ell_min != ell_min) ||
        (sbilv->ode_operator_temp_b_ell_max != ell_max))
    {
      ncm_sbessel_ode_solver_reconfigure_operator (sbilv->ode_solver, op, a_p, b_p, ell_min, ell_max);
      _ncm_sbessel_integrator_levin_apply_constraint (sbilv, op, a_p, b_p, ell_max);
      sbilv->ode_operator_temp_b_a       = a_p;
      sbilv->ode_operator_temp_b_b       = b_p;
      sbilv->ode_operator_temp_b_ell_min = ell_min;
      sbilv->ode_operator_temp_b_ell_max = ell_max;
      sbilv->ode_operator_temp_b_valid   = TRUE;
    }
  }
  else
  {
    op = g_ptr_array_index (sbilv->operators, a_p_idx);
  }

  return op;
}

/* Largest extension growth an edge cell may show and still be used. */
#define NCM_SBESSEL_LEVIN_EDGE_GROWTH_MAX 1.0e4

/*
 * A panel's contribution to ell's result is
 * b_p * j_ell(b_p) * u'(b_p) - a_p * j_ell(a_p) * u'(a_p). When j_ell has
 * underflowed to zero at both endpoints for every ell in the batch, the
 * contribution is exactly zero whatever the solve would produce, so the panel
 * can be skipped without building or solving anything.
 */
static gboolean
_ncm_sbessel_integrator_levin_panel_is_null (const gdouble *j_a_p, const gdouble *j_b_p,
                                             gdouble a_p, gdouble b_p,
                                             guint ell_min, guint ell_max)
{
  guint ell;

  for (ell = ell_min; ell <= ell_max; ell++)
  {
    if ((b_p * fabs (j_b_p[ell]) != 0.0) || (a_p * fabs (j_a_p[ell]) != 0.0))
      return FALSE;
  }

  return TRUE;
}

/*
 * Endpoint values of the panel fit needed by the integration-by-parts
 * boundary terms of the derivative-weighted integrals. The fit in
 * cheb_coeffs holds F(x) on [fit_a, fit_b]; the derivative values are
 * only used (and only computed) for deriv == 2.
 */
typedef struct _NcmSBesselIntegratorLevinBoundary
{
  gdouble F_a;
  gdouble F_b;
  gdouble dF_a;
  gdouble dF_b;
} NcmSBesselIntegratorLevinBoundary;

static void
_ncm_sbessel_integrator_levin_boundary_data (NcmSBesselIntegratorLevin *sbilv,
                                             gdouble fit_a, gdouble fit_b,
                                             gdouble x_a, gdouble x_b,
                                             NcmSBesselIntegratorLevinBoundary *bd)
{
  bd->F_a  = ncm_spectral_chebyshev_eval_x (sbilv->cheb_coeffs, fit_a, fit_b, x_a);
  bd->F_b  = ncm_spectral_chebyshev_eval_x (sbilv->cheb_coeffs, fit_a, fit_b, x_b);
  bd->dF_a = 0.0;
  bd->dF_b = 0.0;

  if (sbilv->deriv == 2)
  {
    bd->dF_a = ncm_spectral_chebyshev_deriv_x (sbilv->cheb_coeffs, fit_a, fit_b, x_a);
    bd->dF_b = ncm_spectral_chebyshev_deriv_x (sbilv->cheb_coeffs, fit_a, fit_b, x_b);
  }
}

/*
 * Integration-by-parts boundary contribution for one multipole:
 * [F j_ell]_{x_a}^{x_b} for deriv == 1 and
 * [F j_ell' - F' j_ell]_{x_a}^{x_b} for deriv == 2.
 */
static inline gdouble
_ncm_sbessel_integrator_levin_boundary_contrib (NcmSBesselIntegratorLevin *sbilv,
                                                const NcmSBesselIntegratorLevinBoundary *bd,
                                                guint ell,
                                                gdouble x_a, gdouble x_b,
                                                const gdouble *j_a, const gdouble *j_b)
{
  if (sbilv->deriv == 1)
  {
    return bd->F_b * j_b[ell] - bd->F_a * j_a[ell];
  }
  else
  {
    const gdouble jp_a = ncm_sf_sbessel_jl_deriv_from_array (ell, x_a, j_a);
    const gdouble jp_b = ncm_sf_sbessel_jl_deriv_from_array (ell, x_b, j_b);

    return (bd->F_b * jp_b - bd->dF_b * j_b[ell]) - (bd->F_a * jp_a - bd->dF_a * j_a[ell]);
  }
}

/* GPtrArray calls its element free function on every slot it drops, and the twin array
 * keeps a slot per panel whether or not that panel ever fell back. */
static void
_ncm_sbessel_operator_unref_null_ok (gpointer op)
{
  if (op != NULL)
    ncm_sbessel_ode_operator_unref (op);
}

/*
 * The Dirichlet twin of the cached panel operator @panel_idx, built the first time that
 * panel falls back. A fallback must not be recorded on the panel operator, or the
 * constraint would depend on which k arrived first; and flipping the operator's own
 * constraint discards its rotations, so a panel that falls back at every k would
 * refactorize at every k. The twin keeps its own factorization, leaving both the tau
 * attempt and the Dirichlet solve on stored rotations.
 */
static NcmSBesselOdeOperator *
_ncm_sbessel_integrator_levin_dirichlet_twin (NcmSBesselIntegratorLevin *sbilv, gint panel_idx,
                                              gdouble a_p, gdouble b_p,
                                              guint ell_min, guint ell_max)
{
  NcmSBesselOdeOperator *twin = g_ptr_array_index (sbilv->dirichlet_operators, panel_idx);

  if (twin == NULL)
  {
    twin = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, a_p, b_p, ell_min, ell_max);
    ncm_sbessel_ode_operator_set_constraint (twin, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
    g_ptr_array_index (sbilv->dirichlet_operators, panel_idx) = twin;
  }

  return twin;
}

/* Solve the current RHS and add its boundary terms to result_data. */
static void
_ncm_sbessel_integrator_levin_solve_rhs_and_accumulate (NcmSBesselIntegratorLevin *sbilv,
                                                        NcmSBesselOdeOperator *operator,
                                                        gint panel_idx,
                                                        gdouble a_p, gdouble b_p,
                                                        const gdouble *j_a_p,
                                                        const gdouble *j_b_p,
                                                        guint ell_min, guint ell_max,
                                                        gdouble *result_data)
{
  NcmSBesselIntegratorLevinBoundary bd = {0.0, 0.0, 0.0, 0.0};
  guint ell;

  if (sbilv->deriv > 0)
    _ncm_sbessel_integrator_levin_boundary_data (sbilv, a_p, b_p, a_p, b_p, &bd);

  /* Dirichlet data makes u vanish at both ends, so the u term of W drops and only the
   * derivatives are needed. The tau constraint leaves u nonzero there, and the general
   * functional -- the same one the extended cells below use -- is required. */
  const NcmSBesselOdeConstraint constraint = ncm_sbessel_ode_operator_get_constraint (operator);
  const gboolean is_tau                    = constraint == NCM_SBESSEL_ODE_CONSTRAINT_TAU;

  /* Only Dirichlet data makes u vanish at the ends; every other constraint needs the
   * general functional, with the u term kept. The order check and the guard judge the
   * tau constraint alone. */
  gboolean use_values = constraint != NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET;

  if (is_tau && !_ncm_sbessel_integrator_levin_tau_constraint_order_ok (sbilv, operator, a_p, b_p, ell_max))
  {
    sbilv->n_constraint_fallbacks++;
    use_values = FALSE;
  }

  if (use_values)
  {
    ncm_sbessel_ode_operator_solve_values (operator, sbilv->rhs, a_p, b_p, &sbilv->values_result);

    if (is_tau &&
        _ncm_sbessel_integrator_levin_tau_constraint_blew_up (sbilv, operator, sbilv->cheb_coeffs, a_p, b_p, ell_min, ell_max))
    {
      /* The tau solve admitted homogeneous content: this panel is not in the regime
       * where that constraint is valid, for this right-hand side. */
      sbilv->n_constraint_fallbacks++;
      use_values = FALSE;
    }
  }

  if (!use_values)
  {
    /* A fallback holds for this right-hand side alone: the operator keeps the constraint
     * the rule gave it, so the choice is a function of the panel, the block and the
     * right-hand side in hand. Recording it instead, as this did before, made the
     * result depend on which k reached the panel first. */
    if (constraint != NCM_SBESSEL_ODE_CONSTRAINT_TAU)
    {
      ncm_sbessel_ode_operator_solve_endpoints (operator, sbilv->rhs, &sbilv->endpoints_result);
    }
    else if (panel_idx >= 0)
    {
      NcmSBesselOdeOperator *twin =
        _ncm_sbessel_integrator_levin_dirichlet_twin (sbilv, panel_idx, a_p, b_p, ell_min, ell_max);

      ncm_sbessel_ode_operator_solve_endpoints (twin, sbilv->rhs, &sbilv->endpoints_result);
    }
    else
    {
      /* A moving operator is reconfigured on every k in any case, so there is no stored
       * factorization for the flip to discard. */
      ncm_sbessel_ode_operator_set_constraint (operator, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
      ncm_sbessel_ode_operator_solve_endpoints (operator, sbilv->rhs, &sbilv->endpoints_result);
      ncm_sbessel_ode_operator_set_constraint (operator, NCM_SBESSEL_ODE_CONSTRAINT_TAU);
    }
  }

  sbilv->n_panel_solves++;

  if (use_values && is_tau)
    sbilv->n_tau_solves++;
  else if ((sbilv->tau_constraint_min_osc > 0.0) &&
           (_ncm_sbessel_integrator_levin_osc (a_p, b_p, ell_max) > sbilv->tau_constraint_min_osc) &&
           _ncm_sbessel_integrator_levin_tau_constraint_order_ok (sbilv, operator, a_p, b_p, ell_max))
    sbilv->n_locked_eligible_solves++;

  for (ell = ell_min; ell <= ell_max; ell++)
  {
    const gint ell_idx  = ell - ell_min;
    const gdouble j_l_a = j_a_p[ell];
    const gdouble j_l_b = j_b_p[ell];
    gdouble contrib;

    if (use_values)
    {
      const gdouble *values = &g_array_index (sbilv->values_result, gdouble, 4 * ell_idx);
      const gdouble xj_p_a  = ncm_sf_sbessel_xjl_deriv_from_array (ell, a_p, j_a_p);
      const gdouble xj_p_b  = ncm_sf_sbessel_xjl_deriv_from_array (ell, b_p, j_b_p);

      contrib = (b_p * j_l_b * values[3] - xj_p_b * values[2]) -
                (a_p * j_l_a * values[1] - xj_p_a * values[0]);
    }
    else
    {
      const gdouble u_prime_a = g_array_index (sbilv->endpoints_result, gdouble, ell_idx * 3 + 0);
      const gdouble u_prime_b = g_array_index (sbilv->endpoints_result, gdouble, ell_idx * 3 + 1);

      contrib = b_p * j_l_b * u_prime_b - a_p * j_l_a * u_prime_a;
    }

    if (sbilv->deriv > 0)
      contrib += _ncm_sbessel_integrator_levin_boundary_contrib (sbilv, &bd, ell, a_p, b_p, j_a_p, j_b_p);

    result_data[ell_idx] += contrib;

    if (G_UNLIKELY (sbilv->record_panels))
    {
      const NcmSBesselIntegratorLevinPanelRec rec = {a_p, b_p, (gint) ell, contrib};

      g_array_append_val (sbilv->panel_records, rec);
    }
  }
}

static void
_ncm_sbessel_integrator_levin_solve_and_accumulate (NcmSBesselIntegratorLevin *sbilv,
                                                    NcmSpectral *spectral,
                                                    NcmSBesselOdeOperator *operator,
                                                    gint panel_idx,
                                                    NcmSBesselIntegratorF F,
                                                    gdouble a_p, gdouble b_p,
                                                    const gdouble *j_a_p,
                                                    const gdouble *j_b_p,
                                                    gdouble k,
                                                    guint ell_min, guint ell_max,
                                                    gdouble *result_data,
                                                    gpointer user_data)
{
  if (_ncm_sbessel_integrator_levin_panel_is_null (j_a_p, j_b_p, a_p, b_p,
                                                   ((sbilv->deriv > 0) && (ell_min > 0)) ? ell_min - 1 : ell_min,
                                                   ell_max))
    return;

  _ncm_sbessel_integrator_levin_compute_rhs (sbilv, spectral, F, a_p, b_p, k, user_data);
  _ncm_sbessel_integrator_levin_solve_rhs_and_accumulate (sbilv, operator, panel_idx,
                                                          a_p, b_p, j_a_p, j_b_p,
                                                          ell_min, ell_max, result_data);
}

static NcmSBesselOdeOperator *
_ncm_sbessel_integrator_levin_get_edge_operator (NcmSBesselIntegratorLevin *sbilv,
                                                 guint panel_idx, gboolean right_edge,
                                                 gdouble integral_a, gdouble integral_b,
                                                 guint ell_min, guint ell_max,
                                                 gdouble *panel_a, gdouble *panel_b)
{
  const gdouble coarse_a = g_array_index (sbilv->knots, gdouble, panel_idx);
  const gdouble coarse_b = g_array_index (sbilv->knots, gdouble, panel_idx + 1);
  const gdouble span     = integral_b - integral_a;
  gdouble width          = coarse_b - coarse_a;
  guint level            = 0;

  g_assert_cmpfloat (span, >, 0.0);

  /* Select the smallest cell of a half-octave ladder that contains the edge, so
   * the piece fills more than 1/sqrt(2) of it. The fill sets how far the piece
   * fit is extrapolated when it is rebased onto the cell: on an octave ladder a
   * piece can fill as little as half the cell, and a fit of 17 to 32 coefficients
   * extrapolated over the other half exceeds the growth limit by orders of
   * magnitude, sending about eight percent of the edge pieces to the moving
   * operator and a fresh factorisation at every k. Above a fill of 0.7 the
   * rejections are a few per thousand. A quarter-octave ladder removes them
   * entirely but fails the certified top-hat at ell = 2 by 3e-8 relative. */
  while ((width * M_SQRT1_2 >= span) && (level < 104))
  {
    width *= M_SQRT1_2;
    level++;
  }

  if (right_edge)
  {
    *panel_a = coarse_a;
    *panel_b = coarse_a + width;
  }
  else
  {
    *panel_a = coarse_b - width;
    *panel_b = coarse_b;
  }

  if (level == 0)
  {
    return g_ptr_array_index (sbilv->operators, panel_idx);
  }
  else
  {
    const guint64 key_value   = ((guint64) right_edge << 63) | ((guint64) panel_idx << 16) | level;
    NcmSBesselOdeOperator *op = g_hash_table_lookup (sbilv->edge_operators, &key_value);

    if (op == NULL)
    {
      guint64 *key = g_new (guint64, 1);

      *key = key_value;
      op   = ncm_sbessel_ode_solver_create_operator (sbilv->ode_solver, *panel_a, *panel_b, ell_min, ell_max);
      _ncm_sbessel_integrator_levin_apply_constraint (sbilv, op, *panel_a, *panel_b, ell_max);
      g_hash_table_insert (sbilv->edge_operators, key, op);
    }

    return op;
  }
}

/*
 * Largest value the extension may reach, in units of the fit's own scale.
 *
 * W_a and W_b are evaluated at the true bounds but built from a solution that
 * lives on the whole cell, so an extension larger than the fit by a factor G
 * costs G times the double-precision noise in their difference.  Keeping
 * G * GSL_DBL_EPSILON under the requested accuracy is what stops a poorly
 * conditioned cell from consuming the whole error budget; the constant caps
 * G for the loose tolerances, where cancellation is not the binding limit.
 */
static gdouble
_ncm_sbessel_integrator_levin_edge_growth_limit (NcmSBesselIntegratorLevin *sbilv,
                                                 gdouble                   reference_scale)
{
  const gdouble growth = MIN (NCM_SBESSEL_LEVIN_EDGE_GROWTH_MAX, sbilv->reltol / GSL_DBL_EPSILON);

  return growth * MAX (reference_scale, G_MINDOUBLE);
}

/* Rebase the edge fit onto the fixed cell, rejecting an ill-conditioned one. */
static gboolean
_ncm_sbessel_integrator_levin_transform_edge_coeffs (NcmSBesselIntegratorLevin *sbilv,
                                                     NcmSpectral *spectral,
                                                     gdouble panel_a, gdouble panel_b,
                                                     gdouble integral_a, gdouble integral_b,
                                                     gdouble reference_scale,
                                                     guint effective_len)
{
  /* Two boundary rows make three coefficients the solver's safe minimum. */
  const guint output_len = MAX (effective_len, 3u);
  const gdouble norm     = ncm_spectral_chebyshev_rebase (spectral, sbilv->edge_cheb_coeffs,
                                                          effective_len,
                                                          integral_a, integral_b,
                                                          panel_a, panel_b,
                                                          &sbilv->cheb_coeffs);

  if (sbilv->cheb_coeffs->len < output_len)
  {
    const guint padded_from = sbilv->cheb_coeffs->len;

    g_array_set_size (sbilv->cheb_coeffs, output_len);
    memset (&g_array_index (sbilv->cheb_coeffs, gdouble, padded_from), 0,
            (output_len - padded_from) * sizeof (gdouble));
  }

  return norm <= _ncm_sbessel_integrator_levin_edge_growth_limit (sbilv, reference_scale);
}

static void
_ncm_sbessel_integrator_levin_prepare_extended_rhs_fallback (NcmSBesselIntegratorLevin *sbilv,
                                                             gdouble integral_a, gdouble integral_b)
{
  GArray *tmp = sbilv->cheb_coeffs;

  sbilv->cheb_coeffs      = sbilv->edge_cheb_coeffs;
  sbilv->edge_cheb_coeffs = tmp;
  _ncm_sbessel_integrator_levin_build_rhs (sbilv, integral_a, integral_b);
}

/*
 * Dead-junction extension.
 *
 * W(integral_b) - W(integral_a) does not depend on the forcing outside the piece,
 * so when the forcing is already negligible where the piece meets the rest of the
 * cell, the cell is fitted directly with the true forcing inside the piece and
 * zero outside. The junction is then as smooth as the forcing itself, nothing is
 * extrapolated, and the cached operator is replayed; the piece fit and the rebase
 * are skipped altogether. Two probe evaluations per junction decide.
 */
typedef struct _NcmSBesselIntegratorLevinTruncated
{
  NcmSpectralF inner;
  gpointer inner_data;
  gdouble piece_a;
  gdouble piece_b;
} NcmSBesselIntegratorLevinTruncated;

static gdouble
_ncm_sbessel_integrator_levin_truncated_func (gpointer data, gdouble x)
{
  const NcmSBesselIntegratorLevinTruncated *trunc = (const NcmSBesselIntegratorLevinTruncated *) data;

  if ((x >= trunc->piece_a) && (x <= trunc->piece_b))
    return trunc->inner (trunc->inner_data, x);

  return 0.0;
}

/* Doublings above the piece's own order allowed to the dead-edge fit before giving up. */
#define NCM_SBESSEL_LEVIN_DEAD_EDGE_MAX_LEVEL 12

#define NCM_SBESSEL_LEVIN_DEAD_PROBE_FRACTION (1.0 / 64.0)

/* A cached Dirichlet edge cell stores ~2 span / pi columns; above this many it is
 * served by the moving operator unless the tau constraint applies. */
#define NCM_SBESSEL_LEVIN_EDGE_CACHE_MAX_OSC (200.0)

static gboolean
_ncm_sbessel_integrator_levin_prepare_dead_edge_rhs (NcmSBesselIntegratorLevin *sbilv,
                                                     NcmSpectral *spectral,
                                                     NcmSBesselIntegratorF F,
                                                     gdouble panel_a, gdouble panel_b,
                                                     gdouble integral_a, gdouble integral_b,
                                                     gdouble k, gpointer user_data)
{
  NcmSBesselIntegratorLevinWrapper wrapper = {F, k, user_data};
  NcmSBesselIntegratorLevinTruncated trunc;
  NcmSpectralF func  = _ncm_sbessel_integrator_levin_peek_wrapper_func (sbilv);
  const gdouble span = integral_b - integral_a;
  gdouble scale      = 0.0;
  gdouble dead_level;
  guint i;

  /* Scale of the forcing on the piece from a few interior samples: the dead
   * level is the fit's own discard level relative to it. */
  for (i = 0; i < 4; i++)
  {
    const gdouble probe = fabs (func (&wrapper, integral_a + (0.125 + 0.25 * i) * span));

    scale = GSL_MAX (scale, probe);
  }

  dead_level = 1.0e-4 * sbilv->cheb_reltol * scale;

  if (integral_b < panel_b)
  {
    const gdouble delta = (panel_b - integral_b) * NCM_SBESSEL_LEVIN_DEAD_PROBE_FRACTION;

    if ((fabs (func (&wrapper, integral_b)) > dead_level) ||
        (fabs (func (&wrapper, integral_b - delta)) > dead_level))
      return FALSE;
  }

  if (integral_a > panel_a)
  {
    const gdouble delta = (integral_a - panel_a) * NCM_SBESSEL_LEVIN_DEAD_PROBE_FRACTION;

    if ((fabs (func (&wrapper, integral_a)) > dead_level) ||
        (fabs (func (&wrapper, integral_a + delta)) > dead_level))
      return FALSE;
  }

  /* Truncate only on a side that has an extension: the Chebyshev node at the cell
   * end is computed as mid -/+ half-width and can land an ulp outside the piece,
   * and a zero there would be a jump of the forcing's full size. */
  trunc.inner      = func;
  trunc.inner_data = &wrapper;
  trunc.piece_a    = (integral_a > panel_a) ? integral_a : -G_MAXDOUBLE;
  trunc.piece_b    = (integral_b < panel_b) ? integral_b : G_MAXDOUBLE;

  {
    gboolean converged = FALSE;

    ncm_spectral_compute_chebyshev_coeffs_adaptive_try (spectral, _ncm_sbessel_integrator_levin_truncated_func,
                                                        panel_a, panel_b,
                                                        sbilv->cheb_min_order, NCM_SBESSEL_LEVIN_DEAD_EDGE_MAX_LEVEL,
                                                        sbilv->cheb_reltol, 0.0,
                                                        &sbilv->cheb_coeffs, &trunc, &converged);

    if (!converged)
      return FALSE;
  }

  /* Drop the roundoff tail of the cell fit, as the continued path does for the
   * piece fit: it sets the stored column count of the cached operator. */
  {
    gdouble cell_scale = 0.0, discarded = 0.0;
    guint len = sbilv->cheb_coeffs->len;

    for (i = 0; i < len; i++)
      cell_scale += fabs (g_array_index (sbilv->cheb_coeffs, gdouble, i));

    while (len > 3)
    {
      const gdouble next = discarded + fabs (g_array_index (sbilv->cheb_coeffs, gdouble, len - 1));

      if (next > 1.0e-4 * sbilv->cheb_reltol * cell_scale)
        break;

      discarded = next;
      len--;
    }

    g_array_set_size (sbilv->cheb_coeffs, len);
  }

  _ncm_sbessel_integrator_levin_build_rhs (sbilv, panel_a, panel_b);

  return TRUE;
}

static gboolean
_ncm_sbessel_integrator_levin_prepare_extended_rhs (NcmSBesselIntegratorLevin *sbilv,
                                                    NcmSpectral *spectral,
                                                    NcmSBesselIntegratorF F,
                                                    gdouble panel_a, gdouble panel_b,
                                                    gdouble integral_a, gdouble integral_b,
                                                    gdouble k, gpointer user_data)
{
  NcmSBesselIntegratorLevinWrapper wrapper = {F, k, user_data};
  gdouble reference_scale                  = 0.0;
  gdouble discarded_scale                  = 0.0;
  gdouble discard_limit;
  guint effective_len;
  guint i;

  if (sbilv->dead_edge_cells &&
      _ncm_sbessel_integrator_levin_prepare_dead_edge_rhs (sbilv, spectral, F, panel_a, panel_b,
                                                           integral_a, integral_b, k, user_data))
    return TRUE;

  /* Fit only on the caller's interval.  The resulting polynomial supplies a
   * smooth extension over the fixed cell, so callbacks are never evaluated
   * outside their advertised integration domain. */
  ncm_spectral_compute_chebyshev_coeffs_adaptive_full (spectral, _ncm_sbessel_integrator_levin_peek_wrapper_func (sbilv),
                                                       integral_a, integral_b,
                                                       sbilv->cheb_min_order, sbilv->cheb_reltol, 0.0,
                                                       &sbilv->edge_cheb_coeffs, &wrapper);

  for (i = 0; i < sbilv->edge_cheb_coeffs->len; i++)
    reference_scale += fabs (g_array_index (sbilv->edge_cheb_coeffs, gdouble, i));

  /* Remove only roundoff-level tail coefficients before extrapolation.  Even
   * a 1e-16 coefficient can grow enormously under T_n(alpha t + beta). */
  discard_limit = 1.0e-4 * sbilv->cheb_reltol * reference_scale;
  effective_len = sbilv->edge_cheb_coeffs->len;

  while (effective_len > 1)
  {
    const gdouble next_scale = discarded_scale + fabs (g_array_index (sbilv->edge_cheb_coeffs,
                                                                      gdouble, effective_len - 1));

    if (next_scale > discard_limit)
      break;

    discarded_scale = next_scale;
    effective_len--;
  }

  /* The affine coefficient transform is quadratic in the polynomial degree. */
  if (effective_len > 513)
  {
    _ncm_sbessel_integrator_levin_prepare_extended_rhs_fallback (sbilv, integral_a, integral_b);

    return FALSE;
  }

  /* Extrapolated high-order noise can grow exponentially outside [-1, 1].
   * Reject such a cell and use the moving-panel solver instead. */
  if (!_ncm_sbessel_integrator_levin_transform_edge_coeffs (sbilv, spectral,
                                                            panel_a, panel_b, integral_a, integral_b,
                                                            reference_scale, effective_len))
  {
    _ncm_sbessel_integrator_levin_prepare_extended_rhs_fallback (sbilv, integral_a, integral_b);

    return FALSE;
  }

  _ncm_sbessel_integrator_levin_build_rhs (sbilv, panel_a, panel_b);

  return TRUE;
}

/*
 * Let v_ell(x) = x j_ell(x) and W = v_ell u' - v_ell' u.  The forced ODE
 * gives W' = K(x / k, k) j_ell(x) / k.  We may therefore solve on a larger,
 * fixed cell with any smooth forcing extension: W(integral_b) - W(integral_a)
 * depends only on the forcing between the true bounds.  This is what makes the
 * fixed operator factorization reusable while a or b moves.
 */
static gboolean
_ncm_sbessel_integrator_levin_integrate_extended_panel (NcmSBesselIntegratorLevin *sbilv,
                                                        NcmSpectral *spectral,
                                                        NcmSBesselOdeOperator *operator,
                                                        NcmSBesselIntegratorF F,
                                                        gdouble panel_a, gdouble panel_b,
                                                        gdouble integral_a, gdouble integral_b,
                                                        const gdouble *j_integral_a, const gdouble *j_integral_b,
                                                        gdouble k, guint ell_min, guint ell_max,
                                                        gdouble *result_data, gpointer user_data)
{
  const gdouble *j_a;
  const gdouble *j_b;
  guint ell;

  if (j_integral_a != NULL)
  {
    j_a = j_integral_a;
  }
  else
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, integral_a, sbilv->j_array_a);
    j_a = sbilv->j_array_a;
  }

  if (j_integral_b != NULL)
  {
    j_b = j_integral_b;
  }
  else
  {
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, integral_b, sbilv->j_array_b);
    j_b = sbilv->j_array_b;
  }

  if (_ncm_sbessel_integrator_levin_panel_is_null (j_a, j_b, integral_a, integral_b,
                                                   ((sbilv->deriv > 0) && (ell_min > 0)) ? ell_min - 1 : ell_min,
                                                   ell_max))
    return TRUE;

  /* A cached edge operator must stay cheap: a Dirichlet factorization on a deep
   * cell costs its span in stored columns, and there are tens of cells per panel.
   * Deep cells are served only with the tau constraint; a shallow cell may keep a
   * Dirichlet operator. Everything else goes to the moving operator, which reuses
   * its storage. */
  {
    const gboolean op_tau = ncm_sbessel_ode_operator_get_constraint (operator) == NCM_SBESSEL_ODE_CONSTRAINT_TAU;
    const gboolean deep   = _ncm_sbessel_integrator_levin_osc (panel_a, panel_b, ell_max) > NCM_SBESSEL_LEVIN_EDGE_CACHE_MAX_OSC;

    if (!op_tau && deep)
    {
      _ncm_sbessel_integrator_levin_compute_rhs (sbilv, spectral, F, integral_a, integral_b, k, user_data);

      return FALSE;
    }
  }

  if (!_ncm_sbessel_integrator_levin_prepare_extended_rhs (sbilv, spectral, F,
                                                           panel_a, panel_b, integral_a, integral_b,
                                                           k, user_data))
    return FALSE;

  if ((ncm_sbessel_ode_operator_get_constraint (operator) == NCM_SBESSEL_ODE_CONSTRAINT_TAU) &&
      !_ncm_sbessel_integrator_levin_tau_constraint_order_ok (sbilv, operator, panel_a, panel_b, ell_max))
  {
    /* Not sticky on a cached cell: the piece goes to the moving operator this time. */
    sbilv->n_constraint_fallbacks++;
    _ncm_sbessel_integrator_levin_compute_rhs (sbilv, spectral, F, integral_a, integral_b, k, user_data);

    return FALSE;
  }

  ncm_sbessel_ode_operator_solve_values (operator, sbilv->rhs,
                                         integral_a, integral_b, &sbilv->values_result);

  if ((ncm_sbessel_ode_operator_get_constraint (operator) == NCM_SBESSEL_ODE_CONSTRAINT_TAU) &&
      _ncm_sbessel_integrator_levin_tau_constraint_blew_up (sbilv, operator, sbilv->cheb_coeffs, panel_a, panel_b, ell_min, ell_max))
  {
    sbilv->n_constraint_fallbacks++;
    _ncm_sbessel_integrator_levin_compute_rhs (sbilv, spectral, F, integral_a, integral_b, k, user_data);

    return FALSE;
  }

  sbilv->n_panel_solves++;

  if (ncm_sbessel_ode_operator_get_constraint (operator) == NCM_SBESSEL_ODE_CONSTRAINT_TAU)
    sbilv->n_tau_solves++;
  else if ((sbilv->tau_constraint_min_osc > 0.0) &&
           (_ncm_sbessel_integrator_levin_osc (panel_a, panel_b, ell_max) > sbilv->tau_constraint_min_osc) &&
           _ncm_sbessel_integrator_levin_tau_constraint_order_ok (sbilv, operator, panel_a, panel_b, ell_max))
    sbilv->n_locked_eligible_solves++;

  {
    NcmSBesselIntegratorLevinBoundary bd = {0.0, 0.0, 0.0, 0.0};

    if (sbilv->deriv > 0)
      _ncm_sbessel_integrator_levin_boundary_data (sbilv, panel_a, panel_b, integral_a, integral_b, &bd);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      const guint ell_idx   = ell - ell_min;
      const gdouble *values = &g_array_index (sbilv->values_result, gdouble, 4 * ell_idx);
      const gdouble u_a     = values[0];
      const gdouble du_a    = values[1];
      const gdouble u_b     = values[2];
      const gdouble du_b    = values[3];
      const gdouble xj_p_a  = ncm_sf_sbessel_xjl_deriv_from_array (ell, integral_a, j_a);
      const gdouble xj_p_b  = ncm_sf_sbessel_xjl_deriv_from_array (ell, integral_b, j_b);
      gdouble W_a, W_b, contrib;

      W_a     = integral_a * j_a[ell] * du_a - xj_p_a * u_a;
      W_b     = integral_b * j_b[ell] * du_b - xj_p_b * u_b;
      contrib = W_b - W_a;

      if (sbilv->deriv > 0)
        contrib += _ncm_sbessel_integrator_levin_boundary_contrib (sbilv, &bd, ell, integral_a, integral_b, j_a, j_b);

      result_data[ell_idx] += contrib;

      if (G_UNLIKELY (sbilv->record_panels))
      {
        const NcmSBesselIntegratorLevinPanelRec rec = {integral_a, integral_b, (gint) ell, contrib};

        g_array_append_val (sbilv->panel_records, rec);
      }
    }
  }

  return TRUE;
}

/**
 * _ncm_sbessel_integrator_levin_integrate_panel:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @a_p_idx: knot index for left endpoint (-1 if not on a knot)
 * @b_p_idx: knot index for right endpoint (-1 if not on a knot)
 * @a_p: lower bound of panel
 * @b_p: upper bound of panel
 * @spectral: spectral methods object
 * @F: integrand function
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 * @result_data: array to accumulate results
 * @user_data: user data for integrand
 *
 * High-level wrapper that integrates a single panel by:
 *
 * 1. Acquiring panel resources (j_ell arrays and operator)
 * 2. Solving the Levin ODE and accumulating results
 *
 * This orchestrator function provides a clean interface for panel integration
 * while keeping the underlying implementation modular for testing and reuse.
 */
static void
_ncm_sbessel_integrator_levin_integrate_panel (NcmSBesselIntegratorLevin *sbilv,
                                               gint a_p_idx, gint b_p_idx,
                                               gdouble a_p, gdouble b_p,
                                               NcmSpectral *spectral,
                                               NcmSBesselIntegratorF F,
                                               gdouble k,
                                               guint ell_min, guint ell_max,
                                               gdouble *result_data,
                                               gpointer user_data,
                                               gboolean rhs_ready)
{
  const gdouble *j_a_p, *j_b_p;
  NcmSBesselOdeOperator *op;

  op = _ncm_sbessel_integrator_levin_get_panel_resources (sbilv, a_p_idx, b_p_idx,
                                                          a_p, b_p, ell_min, ell_max,
                                                          &j_a_p, &j_b_p);

  /* Only a panel sitting between two knots has a cached operator, and so a twin. */
  const gint panel_idx = ((a_p_idx >= 0) && (b_p_idx >= 0)) ? a_p_idx : -1;

  if (rhs_ready)
    _ncm_sbessel_integrator_levin_solve_rhs_and_accumulate (sbilv, op, panel_idx,
                                                            a_p, b_p, j_a_p, j_b_p,
                                                            ell_min, ell_max, result_data);
  else
    _ncm_sbessel_integrator_levin_solve_and_accumulate (sbilv, spectral, op, panel_idx,
                                                        F, a_p, b_p, j_a_p, j_b_p, k,
                                                        ell_min, ell_max, result_data, user_data);
}

static void
_ncm_sbessel_integrator_levin_set_ell_range (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);

  /* Chain up : start */
  NCM_SBESSEL_INTEGRATOR_CLASS (ncm_sbessel_integrator_levin_parent_class)->set_ell_range (sbi, ell_min, ell_max);

  /* Skip operator preparation during construction - resources aren't ready yet */
  if (!sbilv->constructed)
    return;

  if ((ell_min != sbilv->alloc_ell_min) || (ell_max != sbilv->alloc_ell_max))
  {
    if (ell_max > sbilv->ell_cache_max)
      g_error ("Requested ell_max (%u) exceeds ell_cache_max (%u). "
               "Increase ell_cache_max to enable caching for the requested range.",
               ell_max, sbilv->ell_cache_max);

    _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, ell_min, ell_max);

    if (_ncm_sbessel_integrator_levin_rebuild_working_knots (sbilv, ell_max))
      _ncm_sbessel_integrator_levin_prepare_jl_knots (sbilv);

    _ncm_sbessel_integrator_levin_prepare_knots_operators (sbilv, ell_min, ell_max);
  }
}

static guint
_ncm_sbessel_integrator_levin_get_ell_threshold (NcmSBesselIntegratorLevin *sbilv, gdouble a, gdouble b)
{
  return 1000000.0;

  /* The threshold is based on the upper bound */
  return (guint) floor (b);
}

static void
_ncm_sbessel_integrator_levin_integrate_direct (NcmSBesselIntegratorLevin *sbilv,
                                                const guint ell_min, guint ell_max,
                                                NcmSBesselIntegratorF F, const gdouble a, const gdouble b, gdouble k,
                                                NcmVector *result, gpointer user_data)
{
  const gdouble x_min                      = k * a; /* Transform to x-space */
  const gdouble x_max                      = k * b;
  const guint N                            = GSL_MAX (256, 4 * ell_max);
  const gdouble dy                         = (x_max - x_min) / N;
  gdouble * restrict result_ptr            = ncm_vector_data (result);
  NcmSBesselIntegratorLevinWrapper wrapper = {F, k, user_data};
  guint i, ell;

  /* The derivative-weighted integrals are implemented by the Levin path only. */
  if (sbilv->deriv > 0)                                                         /* LCOV_EXCL_LINE */
    g_error ("ncm_sbessel_integrator_levin: the direct cubature path does not " /* LCOV_EXCL_LINE */
             "support Bessel-derivative weights; only the Levin path does.");  /* LCOV_EXCL_LINE */

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);
  /* Initialize direct results to zero */
  memset (result_ptr, 0, sizeof (gdouble) * (ell_max - ell_min + 1));
  ell_max = GSL_MIN (ell_max, ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, x_max));

  /* First term */
  {
    const gdouble fa = _ncm_sbessel_integrator_levin_wrapper_func (&wrapper, x_min) / x_min;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, x_min, sbilv->jl_arr);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      result_ptr[ell - ell_min] = fa * sbilv->jl_arr[ell];
    }
  }

  /* Interior terms with alternating weights 4 and 2 */
  for (i = 1; i < N; i++)
  {
    const gdouble x      = x_min + i * dy;
    const gdouble weight = (i % 2 == 1) ? 4.0 : 2.0;
    const gdouble fy     = _ncm_sbessel_integrator_levin_wrapper_func (&wrapper, x) / x;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, x, sbilv->jl_arr);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      result_ptr[ell - ell_min] += weight * fy * sbilv->jl_arr[ell];
    }
  }

  /* Last term */
  {
    const gdouble fb = _ncm_sbessel_integrator_levin_wrapper_func (&wrapper, x_max) / x_max;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, x_max, sbilv->jl_arr);

    for (ell = ell_min; ell <= ell_max; ell++)
    {
      result_ptr[ell - ell_min] += fb * sbilv->jl_arr[ell];
    }
  }

  /* Apply Simpson's rule factor */
  for (ell = ell_min; ell <= ell_max; ell++)
    result_ptr[ell - ell_min] *= dy / 3.0;
}

static void
_ncm_sbessel_integrator_levin_integrate_levin (NcmSBesselIntegratorLevin *sbilv,
                                               const guint ell_min, const guint ell_max,
                                               NcmSBesselIntegratorF F, const gdouble a, const gdouble b, gdouble k,
                                               NcmVector *result, gpointer user_data)
{
  NcmSBesselOdeSolver *solver = sbilv->ode_solver;
  NcmSpectral *spectral       = ncm_sbessel_ode_solver_peek_spectral (solver);
  gdouble *result_data        = ncm_vector_data (result);
  const gdouble x_min         = k * a; /* Transform to x-space */
  const gdouble x_max         = k * b;
  gint first_knot_idx         = -1;
  gint last_knot_idx          = -1;
  const guint n_ell           = ell_max - ell_min + 1;

  /* Lowest order read from the Bessel array in a panel: the deriv > 0
   * boundary terms use j_{ell-1}, so the panel-skip gate must match
   * _ncm_sbessel_integrator_levin_panel_is_null and look one order down. */
  const guint ell_gate = ((sbilv->deriv > 0) && (ell_min > 0)) ? ell_min - 1 : ell_min;
  gdouble first_knot, last_knot;
  guint i;

  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);

  /* Initialize Levin results to zero */
  memset (result_data, 0, sizeof (gdouble) * n_ell);

  /* Find knots within [x_min, x_max] in x-space */
  for (i = 0; i < sbilv->knots->len; i++)
  {
    const gdouble knot = g_array_index (sbilv->knots, gdouble, i);

    if ((knot >= x_min) && (knot <= x_max))
    {
      if (first_knot_idx == -1)
        first_knot_idx = i;

      last_knot_idx = i;
    }
  }

  /* Use knots-based paneling if configured, otherwise single panel mode */
  if ((sbilv->knots->len > 1) && (first_knot_idx != -1) && (last_knot_idx != -1))
  {
    first_knot = g_array_index (sbilv->knots, gdouble, first_knot_idx);
    last_knot  = g_array_index (sbilv->knots, gdouble, last_knot_idx);

    /* Handle edge panel at the start [x_min, first_knot > x_min] if needed */
    if (x_min < first_knot)
    {
      const gdouble a_p = x_min;
      const gdouble b_p = g_array_index (sbilv->knots, gdouble, first_knot_idx);

      if (ell_gate <= ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p))
      {
        if (first_knot_idx > 0)
        {
          const guint panel_idx = first_knot_idx - 1;
          gdouble panel_a, panel_b;
          NcmSBesselOdeOperator *op;

          op = _ncm_sbessel_integrator_levin_get_edge_operator (sbilv, panel_idx, FALSE,
                                                                a_p, b_p, ell_min, ell_max,
                                                                &panel_a, &panel_b);

          if (!_ncm_sbessel_integrator_levin_integrate_extended_panel (sbilv, spectral, op, F,
                                                                       panel_a, panel_b, a_p, b_p,
                                                                       NULL,
                                                                       _ncm_sbessel_integrator_levin_jl_row (sbilv, first_knot_idx),
                                                                       k, ell_min, ell_max, result_data, user_data))
            _ncm_sbessel_integrator_levin_integrate_panel (sbilv, -1, first_knot_idx,
                                                           a_p, b_p, spectral, F, k,
                                                           ell_min, ell_max, result_data, user_data, TRUE);
        }
        else
        {
          _ncm_sbessel_integrator_levin_integrate_panel (sbilv, -1, first_knot_idx,
                                                         a_p, b_p, spectral, F, k,
                                                         ell_min, ell_max, result_data, user_data, FALSE);
        }
      }
    }

    for (i = first_knot_idx; i < (guint) last_knot_idx; i++)
    {
      const gdouble a_p = g_array_index (sbilv->knots, gdouble, i);
      const gdouble b_p = g_array_index (sbilv->knots, gdouble, i + 1);

      if (ell_gate <= ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p))
        _ncm_sbessel_integrator_levin_integrate_panel (sbilv, i, i + 1,
                                                       a_p, b_p, spectral, F, k,
                                                       ell_min, ell_max, result_data, user_data, FALSE);
    }

    if (x_max > last_knot)
    {
      const gdouble a_p = g_array_index (sbilv->knots, gdouble, last_knot_idx);
      const gdouble b_p = x_max;

      if (ell_gate <= ncm_sf_sbessel_array_eval_ell_cutoff (sbilv->sba, b_p))
      {
        if ((guint) last_knot_idx + 1 < sbilv->knots->len)
        {
          const guint panel_idx = last_knot_idx;
          gdouble panel_a, panel_b;
          NcmSBesselOdeOperator *op;

          op = _ncm_sbessel_integrator_levin_get_edge_operator (sbilv, panel_idx, TRUE,
                                                                a_p, b_p, ell_min, ell_max,
                                                                &panel_a, &panel_b);

          if (!_ncm_sbessel_integrator_levin_integrate_extended_panel (sbilv, spectral, op, F,
                                                                       panel_a, panel_b, a_p, b_p,
                                                                       _ncm_sbessel_integrator_levin_jl_row (sbilv, panel_idx),
                                                                       NULL,
                                                                       k, ell_min, ell_max, result_data, user_data))
            _ncm_sbessel_integrator_levin_integrate_panel (sbilv, last_knot_idx, -1,
                                                           a_p, b_p, spectral, F, k,
                                                           ell_min, ell_max, result_data, user_data, TRUE);
        }
        else
        {
          _ncm_sbessel_integrator_levin_integrate_panel (sbilv, last_knot_idx, -1,
                                                         a_p, b_p, spectral, F, k,
                                                         ell_min, ell_max, result_data, user_data, FALSE);
        }
      }
    }
  }
  else
  {
    /* No paneling: integrate over full range [x_min, x_max]
     * Note: For single panel mode without knots, we use ode_operator directly
     * rather than temp operators. This is handled by passing -1 for both indices
     * but requires special handling in get_panel_resources. */
    const gdouble *j_a_p, *j_b_p;
    NcmSBesselOdeOperator *op;

    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, x_min, sbilv->j_array_a);
    ncm_sf_sbessel_array_eval (sbilv->sba, ell_max, x_max, sbilv->j_array_b);
    j_a_p = sbilv->j_array_a;
    j_b_p = sbilv->j_array_b;

    op = sbilv->ode_operator;
    ncm_sbessel_ode_solver_reconfigure_operator (sbilv->ode_solver, op, x_min, x_max, ell_min, ell_max);
    _ncm_sbessel_integrator_levin_apply_constraint (sbilv, op, x_min, x_max, ell_max);

    _ncm_sbessel_integrator_levin_solve_and_accumulate (sbilv, spectral, op, -1,
                                                        F, x_min, x_max, j_a_p, j_b_p, k,
                                                        ell_min, ell_max, result_data, user_data);
  }
}

static void
_ncm_sbessel_integrator_levin_integrate_full (NcmSBesselIntegrator *sbi,
                                              NcmSBesselIntegratorF F,
                                              gdouble a, gdouble b,
                                              gdouble k,
                                              guint deriv,
                                              NcmVector *result,
                                              gpointer user_data)
{
  NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);
  const gdouble x_min              = k * a; /* Transform to x-space */
  const gdouble x_max              = k * b;
  guint ell_min, ell_max;
  guint n_ell, ell_threshold;

  sbilv->deriv = deriv;

  ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
  n_ell         = ell_max - ell_min + 1;
  ell_threshold = _ncm_sbessel_integrator_levin_get_ell_threshold (sbilv, x_min, x_max);

  if (G_UNLIKELY (sbilv->record_panels))
    g_array_set_size (sbilv->panel_records, 0);

  g_assert_cmpuint (ncm_vector_len (result), >=, n_ell);

  /* Ensure resources are allocated */
  _ncm_sbessel_integrator_levin_ensure_prepared (sbilv, sbilv->max_order, ell_min, ell_max);

  /* Use direct cubature integration for high ell values */
  if (ell_threshold <= ell_max)
    _ncm_sbessel_integrator_levin_integrate_direct (sbilv, ell_min, ell_max, F, a, b, k, result, user_data);
  else
    _ncm_sbessel_integrator_levin_integrate_levin (sbilv, ell_min, ell_max, F, a, b, k, result, user_data);
}

static void
_ncm_sbessel_integrator_levin_integrate (NcmSBesselIntegrator *sbi,
                                         NcmSBesselIntegratorF F,
                                         gdouble a, gdouble b,
                                         gdouble k,
                                         NcmVector *result,
                                         gpointer user_data)
{
  _ncm_sbessel_integrator_levin_integrate_full (sbi, F, a, b, k, 0, result, user_data);
}

static void
_ncm_sbessel_integrator_levin_integrate_deriv (NcmSBesselIntegrator *sbi,
                                               NcmSBesselIntegratorF F,
                                               gdouble a, gdouble b,
                                               gdouble k,
                                               guint deriv,
                                               NcmVector *result,
                                               gpointer user_data)
{
  _ncm_sbessel_integrator_levin_integrate_full (sbi, F, a, b, k, deriv, result, user_data);
}

/**
 * ncm_sbessel_integrator_levin_set_record_panels:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @record: whether to record per-panel diagnostics
 *
 * Enables or disables recording of the individual panel contributions.
 *
 * The integral is assembled as a sum of per-panel boundary terms which cancel
 * heavily, so the attainable relative accuracy is bounded by the cancellation
 * ratio $\sum_p |I_p| / |\sum_p I_p|$ times the machine epsilon. That ratio cannot
 * be recovered from the result alone. With recording enabled the contributions are
 * kept and can be read back with
 * ncm_sbessel_integrator_levin_get_panel_contrib().
 *
 * Recording is off by default and costs nothing when off. The records are cleared
 * at the start of every ncm_sbessel_integrator_integrate() call, so they always
 * describe the most recent one.
 */
void
ncm_sbessel_integrator_levin_set_record_panels (NcmSBesselIntegratorLevin *sbilv, gboolean record)
{
  sbilv->record_panels = record;

  if (!record)
    g_array_set_size (sbilv->panel_records, 0);
}

/**
 * ncm_sbessel_integrator_levin_get_record_panels:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets whether per-panel diagnostics are being recorded.
 *
 * Returns: TRUE if recording is enabled.
 */
gboolean
ncm_sbessel_integrator_levin_get_record_panels (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->record_panels;
}

/**
 * ncm_sbessel_integrator_levin_get_n_panel_records:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the number of panel records from the most recent integration. One record is
 * produced per (panel, multipole) pair.
 *
 * Returns: the number of records.
 */
guint
ncm_sbessel_integrator_levin_get_n_panel_records (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->panel_records->len;
}

#define _NCM_SBILV_REC(sbilv, i) (&g_array_index ((sbilv)->panel_records, NcmSBesselIntegratorLevinPanelRec, (i)))

/**
 * ncm_sbessel_integrator_levin_get_panel_a:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @i: record index
 *
 * Gets the lower bound, in $x = kx$, of the panel of record @i.
 *
 * Returns: the panel lower bound.
 */
gdouble
ncm_sbessel_integrator_levin_get_panel_a (NcmSBesselIntegratorLevin *sbilv, guint i)
{
  g_assert_cmpuint (i, <, sbilv->panel_records->len);

  return _NCM_SBILV_REC (sbilv, i)->a;
}

/**
 * ncm_sbessel_integrator_levin_get_panel_b:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @i: record index
 *
 * Gets the upper bound, in $x = kx$, of the panel of record @i.
 *
 * Returns: the panel upper bound.
 */
gdouble
ncm_sbessel_integrator_levin_get_panel_b (NcmSBesselIntegratorLevin *sbilv, guint i)
{
  g_assert_cmpuint (i, <, sbilv->panel_records->len);

  return _NCM_SBILV_REC (sbilv, i)->b;
}

/**
 * ncm_sbessel_integrator_levin_get_panel_ell:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @i: record index
 *
 * Gets the multipole of record @i.
 *
 * Returns: the multipole.
 */
gint
ncm_sbessel_integrator_levin_get_panel_ell (NcmSBesselIntegratorLevin *sbilv, guint i)
{
  g_assert_cmpuint (i, <, sbilv->panel_records->len);

  return _NCM_SBILV_REC (sbilv, i)->ell;
}

/**
 * ncm_sbessel_integrator_levin_get_panel_contrib:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @i: record index
 *
 * Gets the contribution of record @i to the total integral. The sum of the
 * contributions over all records with the same multipole is that multipole's result.
 *
 * Returns: the panel contribution.
 */
gdouble
ncm_sbessel_integrator_levin_get_panel_contrib (NcmSBesselIntegratorLevin *sbilv, guint i)
{
  g_assert_cmpuint (i, <, sbilv->panel_records->len);

  return _NCM_SBILV_REC (sbilv, i)->contrib;
}

/**
 * ncm_sbessel_integrator_levin_new:
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 *
 * Creates a new #NcmSBesselIntegratorLevin with default parameters:
 *
 * - x_knots_min = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_X_KNOTS_MIN
 * - x_knots_max = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_X_KNOTS_MAX
 * - n_knots = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_N_KNOTS
 * - ell_cache_max = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_ELL_CACHE_MAX
 * - reltol = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_RELTOL
 * - cheb_min_order = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_MIN_ORDER
 * - cheb_reltol = %NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_RELTOL
 *
 * Returns: (transfer full): a new #NcmSBesselIntegratorLevin
 */
NcmSBesselIntegratorLevin *
ncm_sbessel_integrator_levin_new (guint ell_min, guint ell_max)
{
  return ncm_sbessel_integrator_levin_new_full (ell_min, ell_max,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_X_KNOTS_MIN,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_X_KNOTS_MAX,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_N_KNOTS,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_ELL_CACHE_MAX,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_RELTOL,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_MIN_ORDER,
                                                NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_RELTOL);
}

/**
 * ncm_sbessel_integrator_levin_new_full:
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 * @x_knots_min: minimum value for knots in log-spaced grid (set to 0 to disable knots-based paneling)
 * @x_knots_max: maximum value for knots in log-spaced grid (set to 0 to disable knots-based paneling)
 * @n_knots: number of knots in the log-spaced grid (set to 0 to disable knots-based paneling)
 * @ell_cache_max: maximum ell value for precomputed spherical Bessel functions at knots
 * @reltol: relative tolerance for integration
 * @cheb_min_order: minimum order of Chebyshev decomposition for RHS computation
 * @cheb_reltol: relative tolerance for Chebyshev decomposition of integrand
 *
 * Creates a new #NcmSBesselIntegratorLevin with optional knots-based paneling. To
 * disable knots-based paneling and use single panel mode, set @x_knots_min,
 * @x_knots_max, or @n_knots to 0.
 *
 * The @ell_cache_max parameter controls the maximum ell value for which spherical
 * Bessel functions will be precomputed at all knots. For ell values beyond this,
 * spherical Bessel functions will be computed on-the-fly during integration.
 *
 * Returns: (transfer full): a new #NcmSBesselIntegratorLevin
 */
NcmSBesselIntegratorLevin *
ncm_sbessel_integrator_levin_new_full (guint ell_min, guint ell_max, gdouble x_knots_min, gdouble x_knots_max, guint n_knots, guint ell_cache_max, gdouble reltol, guint cheb_min_order, gdouble cheb_reltol)
{
  NcmDTuple2 ell_range             = NCM_DTUPLE2_STATIC_INIT ((gdouble) ell_min, (gdouble) ell_max);
  NcmSBesselIntegratorLevin *sbilv = g_object_new (NCM_TYPE_SBESSEL_INTEGRATOR_LEVIN,
                                                   "ell-range", &ell_range,
                                                   "x-knots-min", x_knots_min,
                                                   "x-knots-max", x_knots_max,
                                                   "n-knots", n_knots,
                                                   "ell-cache-max", ell_cache_max,
                                                   "reltol", reltol,
                                                   "cheb-min-order", cheb_min_order,
                                                   "cheb-reltol", cheb_reltol,
                                                   NULL);

  return sbilv;
}

/**
 * ncm_sbessel_integrator_levin_ref:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Increases the reference count of @sbilv by one.
 *
 * Returns: (transfer full): @sbilv
 */
NcmSBesselIntegratorLevin *
ncm_sbessel_integrator_levin_ref (NcmSBesselIntegratorLevin *sbilv)
{
  return g_object_ref (sbilv);
}

/**
 * ncm_sbessel_integrator_levin_free:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Decreases the reference count of @sbilv by one.
 */
void
ncm_sbessel_integrator_levin_free (NcmSBesselIntegratorLevin *sbilv)
{
  g_object_unref (sbilv);
}

/**
 * ncm_sbessel_integrator_levin_clear:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * If @sbilv is different from NULL, decreases the reference count of
 * @sbilv by one and sets @sbilv to NULL.
 */
void
ncm_sbessel_integrator_levin_clear (NcmSBesselIntegratorLevin **sbilv)
{
  g_clear_object (sbilv);
}

/**
 * ncm_sbessel_integrator_levin_set_max_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @max_order: maximum order
 *
 * Sets the maximum order of Clenshaw-Curtis quadrature.
 */
void
ncm_sbessel_integrator_levin_set_max_order (NcmSBesselIntegratorLevin *sbilv, guint max_order)
{
  sbilv->max_order = max_order;
}

/**
 * ncm_sbessel_integrator_levin_get_max_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the maximum order.
 *
 * Returns: the maximum order
 */
guint
ncm_sbessel_integrator_levin_get_max_order (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->max_order;
}

/**
 * ncm_sbessel_integrator_levin_set_reltol:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @reltol: relative tolerance
 *
 * Sets the ODE solve tolerance. An operator's accuracy is fixed when it is
 * built, so this discards every panel operator, cached edge cell and
 * moving-edge temporary and rebuilds them. Setting the current value is a
 * no-op.
 *
 * The rebuild is expensive and invalidates factorizations in use: prefer
 * ncm_sbessel_integrator_levin_new_full() or the `reltol` property.
 *
 * The integrand fit has a separate tolerance,
 * ncm_sbessel_integrator_levin_set_cheb_reltol(). The looser of the two
 * bounds the result.
 */

/*
 * Rebuilds the constraint of every cached panel operator after a constraint setting
 * changed. The operators are reconfigured in place through the set_ell_range path,
 * which is the only place the preparation guard lives: a plain integrate() call
 * never re-prepares, so invalidating the allocated range alone would leave the
 * operators untouched until the next range change. Before construction nothing
 * exists yet and the constructed handler prepares with the setting already in place.
 */
static void
_ncm_sbessel_integrator_levin_reprepare (NcmSBesselIntegratorLevin *sbilv)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (sbilv);
  guint ell_min, ell_max;

  g_hash_table_remove_all (sbilv->edge_operators);
  sbilv->ode_operator_temp_a_valid = FALSE;
  sbilv->ode_operator_temp_b_valid = FALSE;

  if (!sbilv->constructed)
    return;

  sbilv->alloc_ell_min = G_MAXUINT;
  sbilv->alloc_ell_max = G_MAXUINT;
  ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
  ncm_sbessel_integrator_set_ell_range (sbi, ell_min, ell_max);
}

/*
 * Returns TRUE when the last tau-constraint solve on @op admitted homogeneous content.
 * The smooth member is u_p ~ yF / (x^2 - nu^2), so its coefficients are bounded by
 * max|yF| / min(x^2 - nu^2); @cheb holds the Chebyshev coefficients of F on the
 * panel, so b * sum|c_k| bounds max|yF|.
 */
static gboolean
_ncm_sbessel_integrator_levin_tau_constraint_blew_up (NcmSBesselIntegratorLevin *sbilv, NcmSBesselOdeOperator *op, GArray *cheb, gdouble a, gdouble b, guint ell_min, guint ell_max)
{
  gdouble sum_c = 0.0;
  guint i;

  for (i = 0; i < cheb->len; i++)
    sum_c += fabs (g_array_index (cheb, gdouble, i));

  /* Each member of the block is judged against its own turning point: the bound on
   * the smooth member is max|yF| / (a^2 - nu_ell^2), and using the block's largest
   * nu for every member loosens it by (a^2 - nu_ell^2) / (a^2 - nu_max^2) for the low
   * ones, which on a short panel is an order of magnitude. A member whose turning
   * point lies inside or above the panel has no smooth member and is rejected. */
  for (i = 0; i <= ell_max - ell_min; i++)
  {
    const gdouble ell   = (gdouble) (ell_min + i);
    const gdouble denom = a * a - ell * (ell + 1.0);

    if (denom <= 0.0)
      return TRUE;

    {
      const gdouble max_coeff = ncm_sbessel_ode_operator_get_last_max_coeff (op, i);

      if (!isfinite (max_coeff) || (max_coeff > sbilv->tau_constraint_guard * b * sum_c / denom))
        return TRUE;
    }
  }

  return FALSE;
}

/**
 * ncm_sbessel_integrator_levin_set_turning_knot_margin:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @margin: multiple of $\nu_{\max}$, or 0 to add no knot
 *
 * Sets #NcmSBesselIntegratorLevin:turning-knot-margin. The panel grid and every cached
 * operator are rebuilt, since the grid itself is what changes.
 *
 */
void
ncm_sbessel_integrator_levin_set_turning_knot_margin (NcmSBesselIntegratorLevin *sbilv, gdouble margin)
{
  g_assert_cmpfloat (margin, >=, 0.0);

  if (sbilv->turning_knot_margin == margin)
    return;

  sbilv->turning_knot_margin = margin;

  if (sbilv->constructed)
    _ncm_sbessel_integrator_levin_reprepare (sbilv);
}

/**
 * ncm_sbessel_integrator_levin_get_turning_knot_margin:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Returns: the value of #NcmSBesselIntegratorLevin:turning-knot-margin.
 */
gdouble
ncm_sbessel_integrator_levin_get_turning_knot_margin (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->turning_knot_margin;
}

/**
 * ncm_sbessel_integrator_levin_set_tau_constraint_guard:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @guard: factor above the forcing's scale that rejects a tau-constraint solve
 *
 * Sets #NcmSBesselIntegratorLevin:tau-constraint-guard. Read at every solve; no
 * cached state depends on it.
 *
 */
void
ncm_sbessel_integrator_levin_set_tau_constraint_guard (NcmSBesselIntegratorLevin *sbilv, gdouble guard)
{
  g_assert_cmpfloat (guard, >=, 1.0);
  sbilv->tau_constraint_guard = guard;
}

/**
 * ncm_sbessel_integrator_levin_get_tau_constraint_guard:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Returns: the value of #NcmSBesselIntegratorLevin:tau-constraint-guard.
 */
gdouble
ncm_sbessel_integrator_levin_get_tau_constraint_guard (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->tau_constraint_guard;
}

/**
 * ncm_sbessel_integrator_levin_set_tau_constraint_order_fraction:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @fraction: largest tau-constraint working order as a fraction of $N_\min$
 *
 * Sets #NcmSBesselIntegratorLevin:tau-constraint-order-fraction. Takes effect on
 * the next solve; cached operators are not touched.
 */
void
ncm_sbessel_integrator_levin_set_tau_constraint_order_fraction (NcmSBesselIntegratorLevin *sbilv, gdouble fraction)
{
  g_assert_cmpfloat (fraction, >, 0.0);
  g_assert_cmpfloat (fraction, <=, 1.0);
  sbilv->tau_constraint_order_fraction = fraction;
}

/**
 * ncm_sbessel_integrator_levin_get_tau_constraint_order_fraction:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Returns: the value of #NcmSBesselIntegratorLevin:tau-constraint-order-fraction.
 */
gdouble
ncm_sbessel_integrator_levin_get_tau_constraint_order_fraction (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->tau_constraint_order_fraction;
}

/**
 * ncm_sbessel_integrator_levin_set_dead_edge_cells:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @enable: whether to use cached cells for dead-junction end pieces
 *
 * Sets #NcmSBesselIntegratorLevin:dead-edge-cells.
 */
void
ncm_sbessel_integrator_levin_set_dead_edge_cells (NcmSBesselIntegratorLevin *sbilv, gboolean enable)
{
  sbilv->dead_edge_cells = enable;
}

/**
 * ncm_sbessel_integrator_levin_get_dead_edge_cells:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Returns: the value of #NcmSBesselIntegratorLevin:dead-edge-cells.
 */
gboolean
ncm_sbessel_integrator_levin_get_dead_edge_cells (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->dead_edge_cells;
}

/**
 * ncm_sbessel_integrator_levin_get_n_constraint_fallbacks:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Returns: how many panel operators the tau-constraint guard has sent back to
 * Dirichlet data since construction.
 */
guint
ncm_sbessel_integrator_levin_get_n_constraint_fallbacks (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->n_constraint_fallbacks;
}

/**
 * ncm_sbessel_integrator_levin_get_n_panel_solves:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Number of panel solves performed by this integrator, cumulative like the
 * fallback counter.
 *
 * Returns: the panel solve count.
 */
guint
ncm_sbessel_integrator_levin_get_n_panel_solves (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->n_panel_solves;
}

/**
 * ncm_sbessel_integrator_levin_get_n_tau_solves:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Number of panel solves that used, and kept, the tau constraint. Together with
 * ncm_sbessel_integrator_levin_get_n_panel_solves() this measures how often the
 * constraint rule applies to a given forcing.
 *
 * Returns: the tau-constraint solve count.
 */
guint
ncm_sbessel_integrator_levin_get_n_tau_solves (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->n_tau_solves;
}

/**
 * ncm_sbessel_integrator_levin_get_n_locked_eligible_solves:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Number of Dirichlet panel solves on a panel the constraint rule admits, whose
 * forcing would have passed the order check. A fallback is not recorded on the panel
 * operator, so this counts a panel whose resting constraint came from a different
 * block range rather than one lost to an earlier k.
 *
 * Returns: the locked-eligible solve count.
 */
guint
ncm_sbessel_integrator_levin_get_n_locked_eligible_solves (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->n_locked_eligible_solves;
}

/*
 * Applies the per-panel constraint rule to an operator that has just been given the
 * bounds [a, b] for multipoles up to ell_max. No-op when the rule is disabled.
 */

/* Oscillation count of the panel above the block's highest turning point. */
static gdouble
_ncm_sbessel_integrator_levin_osc (gdouble a, gdouble b, guint ell_max)
{
  const gdouble nu_max = sqrt (ell_max * (ell_max + 1.0));
  const gdouble span   = b - GSL_MAX (a, nu_max);

  return (span > 0.0) ? 2.0 * span / M_PI : 0.0;
}

static void
_ncm_sbessel_integrator_levin_apply_constraint (NcmSBesselIntegratorLevin *sbilv, NcmSBesselOdeOperator *op, gdouble a, gdouble b, guint ell_max)
{
  if (sbilv->tau_constraint_min_osc <= 0.0)
    return;

  if (_ncm_sbessel_integrator_levin_osc (a, b, ell_max) > sbilv->tau_constraint_min_osc)
  {
    ncm_sbessel_ode_operator_set_constraint (op, NCM_SBESSEL_ODE_CONSTRAINT_TAU);

    return;
  }

  /*
   * A panel the tau rule turns down goes to the pinned constraint, with the pins at
   * the peak of the homogeneous spectrum. Above the turning point the Chebyshev
   * coefficients of $x j_\ell$ and $x y_\ell$ behave as $|J_n(\Delta)|$, which peaks
   * at $n = \Delta - 0.81\Delta^{1/3}$, so pinning the adjacent pair there sets the
   * two conditions where the homogeneous content is largest and the multiple the
   * truncation can admit, $a_n(u_p)/a_n(v)$, is smallest. Two adjacent indices cover
   * both parities, which is what makes the choice insensitive to where the panel sits
   * in phase.
   *
   * $\Delta = (b - \max(a,\nu))/2$ is the half-width of the oscillatory part of the
   * panel, taken over the same span the oscillation count uses and not over the full
   * $b-a$. Below the turning point the coefficients have already collapsed, so a
   * $\Delta$ from the full width puts the pins where there is no homogeneous content
   * left to pin and the two conditions go nearly degenerate: at $\ell = 500$ that
   * costs O(1) relative error on every window of the certified table. A panel too
   * short for the peak to leave room keeps Dirichlet data.
   */
  {
    const gdouble nu_max = sqrt (ell_max * (ell_max + 1.0));
    const gdouble delta  = 0.5 * (b - GSL_MAX (a, nu_max));
    const glong peak     = (glong) round (delta - 0.81 * cbrt (delta));

    if (peak >= 4)
    {
      ncm_sbessel_ode_operator_set_pinned_constraint (op, peak, peak + 1);

      return;
    }
  }

  ncm_sbessel_ode_operator_set_constraint (op, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
}

/*
 * The tau constraint is valid while the working order sits on the plateau of the
 * homogeneous spectrum, |a_n(v)| = O(1) for n below N_min ~ 2 span / pi. Past it
 * a_n(v) collapses super-exponentially and the truncation admits a homogeneous
 * multiple of size a_n(u_p) / a_n(v). The working order is set by the forcing, so
 * this check needs the right-hand side and runs at solve time. An unconverged
 * forcing fit (order at the cap) is the usual way to get here.
 */
static gboolean
_ncm_sbessel_integrator_levin_tau_constraint_order_ok (NcmSBesselIntegratorLevin *sbilv, NcmSBesselOdeOperator *op, gdouble a, gdouble b, guint ell_max)
{
  const glong order   = ncm_sbessel_ode_operator_get_tau_constraint_order (op, sbilv->rhs->len);
  const gdouble n_min = _ncm_sbessel_integrator_levin_osc (a, b, ell_max);

  return order < sbilv->tau_constraint_order_fraction * n_min;
}

/**
 * ncm_sbessel_integrator_levin_set_tau_constraint_min_osc:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @min_osc: oscillation count above which a panel uses the tau constraint; 0 disables
 *
 * Sets #NcmSBesselIntegratorLevin:tau-constraint-min-osc. Cached panel operators are
 * reconfigured on the next solve so the rule reaches them.
 *
 */
void
ncm_sbessel_integrator_levin_set_tau_constraint_min_osc (NcmSBesselIntegratorLevin *sbilv, gdouble min_osc)
{
  g_assert_cmpfloat (min_osc, >=, 0.0);

  if (sbilv->tau_constraint_min_osc == min_osc)
    return;

  sbilv->tau_constraint_min_osc = min_osc;

  _ncm_sbessel_integrator_levin_reprepare (sbilv);
}

/**
 * ncm_sbessel_integrator_levin_get_tau_constraint_min_osc:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Returns: the value of #NcmSBesselIntegratorLevin:tau-constraint-min-osc.
 */
gdouble
ncm_sbessel_integrator_levin_get_tau_constraint_min_osc (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->tau_constraint_min_osc;
}

/**
 * ncm_sbessel_integrator_levin_set_tau_constraint:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @tau_constraint: whether panel solves use the tau constraint
 *
 * Selects the constraint used by every panel solve. The panel operators carry the
 * constraint they were built with, so the cached ones are discarded here and rebuilt on
 * the next solve; setting the flag on the solver alone would leave them untouched.
 *
 */
void
ncm_sbessel_integrator_levin_set_tau_constraint (NcmSBesselIntegratorLevin *sbilv, gboolean tau_constraint)
{
  ncm_sbessel_ode_solver_set_default_constraint (sbilv->ode_solver,
                                                 tau_constraint ? NCM_SBESSEL_ODE_CONSTRAINT_TAU :
                                                 NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);

  _ncm_sbessel_integrator_levin_reprepare (sbilv);
}

/**
 * ncm_sbessel_integrator_levin_get_tau_constraint:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Returns: whether panel solves use the tau constraint.
 */
gboolean
ncm_sbessel_integrator_levin_get_tau_constraint (NcmSBesselIntegratorLevin *sbilv)
{
  return ncm_sbessel_ode_solver_get_default_constraint (sbilv->ode_solver) == NCM_SBESSEL_ODE_CONSTRAINT_TAU;
}

/**
 * ncm_sbessel_integrator_levin_peek_ode_solver:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gives access to the #NcmSBesselOdeSolver the integrator builds its panel operators
 * from. Settings applied here reach operators created or reconfigured afterwards.
 *
 * Returns: (transfer none): the integrator's solver
 */
NcmSBesselOdeSolver *
ncm_sbessel_integrator_levin_peek_ode_solver (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->ode_solver;
}

void
ncm_sbessel_integrator_levin_set_reltol (NcmSBesselIntegratorLevin *sbilv, gdouble reltol)
{
  g_assert_cmpfloat (reltol, >, 0.0);

  if (sbilv->reltol == reltol)
    return;

  sbilv->reltol = reltol;
  ncm_sbessel_ode_solver_set_tolerance (sbilv->ode_solver, reltol);

  if (sbilv->operators == NULL)
  {
    _ncm_sbessel_integrator_levin_prepare_knots_operators (sbilv, sbilv->alloc_ell_min, sbilv->alloc_ell_max);
  }
  else
  {
    /* Each operator holds its own copy of the tolerance, so all of them have to be
     * moved to the new one. Their intervals and multipole range do not change. */
    const guint ell_min = sbilv->alloc_ell_min;
    const guint ell_max = sbilv->alloc_ell_max;
    guint i;

    g_hash_table_remove_all (sbilv->edge_operators);

    for (i = 0; i < sbilv->operators->len; i++)
    {
      NcmSBesselOdeOperator *op = g_ptr_array_index (sbilv->operators, i);
      const gdouble x_a         = g_array_index (sbilv->knots, gdouble, i);
      const gdouble x_b         = g_array_index (sbilv->knots, gdouble, i + 1);

      ncm_sbessel_ode_solver_reconfigure_operator (sbilv->ode_solver, op, x_a, x_b, ell_min, ell_max);
      _ncm_sbessel_integrator_levin_apply_constraint (sbilv, op, x_a, x_b, ell_max);
    }

    ncm_sbessel_ode_solver_reconfigure_operator (sbilv->ode_solver, sbilv->ode_operator_temp_a, 0.0, 1.0, ell_min, ell_max);
    ncm_sbessel_ode_solver_reconfigure_operator (sbilv->ode_solver, sbilv->ode_operator_temp_b, 0.0, 1.0, ell_min, ell_max);
    sbilv->ode_operator_temp_a_valid = FALSE;
    sbilv->ode_operator_temp_b_valid = FALSE;
  }
}

/**
 * ncm_sbessel_integrator_levin_get_reltol:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the ODE solve tolerance.
 *
 * Returns: the ODE solve relative tolerance
 */
gdouble
ncm_sbessel_integrator_levin_get_reltol (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->reltol;
}

/**
 * ncm_sbessel_integrator_levin_set_cheb_min_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @cheb_min_order: minimum Chebyshev order
 *
 * Sets the minimum order of Chebyshev decomposition for RHS computation.
 */
void
ncm_sbessel_integrator_levin_set_cheb_min_order (NcmSBesselIntegratorLevin *sbilv, guint cheb_min_order)
{
  sbilv->cheb_min_order = cheb_min_order;
}

/**
 * ncm_sbessel_integrator_levin_get_cheb_min_order:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the minimum Chebyshev order for RHS computation.
 *
 * Returns: the minimum Chebyshev order
 */
guint
ncm_sbessel_integrator_levin_get_cheb_min_order (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->cheb_min_order;
}

/**
 * ncm_sbessel_integrator_levin_set_cheb_reltol:
 * @sbilv: a #NcmSBesselIntegratorLevin
 * @cheb_reltol: Chebyshev decomposition relative tolerance
 *
 * Sets the tolerance of the integrand's Chebyshev fit. The fit is redone per
 * panel, so this applies from the next integration and rebuilds nothing,
 * unlike ncm_sbessel_integrator_levin_set_reltol(). The looser of the two
 * bounds the result.
 */
void
ncm_sbessel_integrator_levin_set_cheb_reltol (NcmSBesselIntegratorLevin *sbilv, gdouble cheb_reltol)
{
  g_assert_cmpfloat (cheb_reltol, >, 0.0);

  sbilv->cheb_reltol = cheb_reltol;
}

/**
 * ncm_sbessel_integrator_levin_get_cheb_reltol:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the tolerance of the integrand's Chebyshev fit.
 *
 * Returns: the Chebyshev fit relative tolerance
 */
gdouble
ncm_sbessel_integrator_levin_get_cheb_reltol (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->cheb_reltol;
}

/**
 * ncm_sbessel_integrator_levin_get_x_knots_min:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the minimum knot value.
 *
 * Returns: the minimum knot value
 */
gdouble
ncm_sbessel_integrator_levin_get_x_knots_min (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->x_knots_min;
}

/**
 * ncm_sbessel_integrator_levin_get_x_knots_max:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the maximum knot value.
 *
 * Returns: the maximum knot value
 */
gdouble
ncm_sbessel_integrator_levin_get_x_knots_max (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->x_knots_max;
}

/**
 * ncm_sbessel_integrator_levin_get_n_knots:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Number of knots in the base log-spaced grid, the value
 * #NcmSBesselIntegratorLevin:n-knots was built with. The grid a block is actually solved
 * on carries one more when
 * #NcmSBesselIntegratorLevin:turning-knot-margin puts a knot above that block's turning
 * point.
 *
 * Returns: the number of base knots
 */
guint
ncm_sbessel_integrator_levin_get_n_knots (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->n_knots;
}

/**
 * ncm_sbessel_integrator_levin_get_ell_cache_max:
 * @sbilv: a #NcmSBesselIntegratorLevin
 *
 * Gets the maximum ell value for the precomputed spherical Bessel cache.
 *
 * Returns: the maximum ell value for cache
 */
guint
ncm_sbessel_integrator_levin_get_ell_cache_max (NcmSBesselIntegratorLevin *sbilv)
{
  return sbilv->ell_cache_max;
}

