/***************************************************************************
 *            test_nc_xcor_kquad_common.c
 *
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

/*
 * What the kernel-space outer integral tests share: the cosmology, distances, power
 * spectrum and Levin integrator they are built on, the top-hat kernels they integrate,
 * and an independent reference for the integral itself.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_integration.h>

#include "test_nc_xcor_kquad_common.h"

/*
 * Flat XCDM with w = -1, a BBKS spectrum capped well below its default, distances to
 * z = 6 and one Levin integrator covering multipoles 0 to 8, all prepared.
 */
void
test_nc_xcor_kquad_env_init (TestNcXcorKQuadEnv *env)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_H0,      70.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_C, 0.255);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.7);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_B, 0.045);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W, -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "Omegak", 0.0, NULL);

  env->cosmo = cosmo;
  env->dist  = nc_distance_new (TEST_NC_XCOR_KQUAD_ZMAX);
  env->ps    = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                      NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  env->sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 8));

  ncm_powspec_set_kmax (env->ps, TEST_NC_XCOR_KQUAD_KMAX);

  nc_distance_prepare (env->dist, cosmo);
  ncm_powspec_prepare (env->ps, NCM_MODEL (cosmo));
}

void
test_nc_xcor_kquad_env_clear (TestNcXcorKQuadEnv *env)
{
  ncm_sbessel_integrator_free (env->sbi);
  ncm_powspec_free (env->ps);
  nc_distance_free (env->dist);
  nc_hicosmo_free (env->cosmo);
}

/*
 * A top-hat kernel on [chi_lower, chi_upper] in the given tier, prepared. Its adaptive
 * apparatus is capped rather than run to convergence: what these tests check is the
 * outer integral of whatever closure the kernel builds, not the closure's accuracy.
 */
NcXcorKernel *
test_nc_xcor_kquad_tophat (TestNcXcorKQuadEnv *env, gdouble chi_lower, gdouble chi_upper, gint l_limber)
{
  NcXcorKernel *xclk = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_new_full (env->dist, env->ps,
                                                                                chi_lower, chi_upper,
                                                                                env->sbi));

  nc_xcor_kernel_set_max_border_expansions (xclk, 1);
  nc_xcor_kernel_set_max_iter (xclk, 4);
  nc_xcor_kernel_set_reltol (xclk, 1.0e-3);
  nc_xcor_kernel_set_peak_epsilon (xclk, 1.0e-4);
  nc_xcor_kernel_set_panel_order_cap (xclk, 12);
  nc_xcor_kernel_set_lmax (xclk, 16);
  nc_xcor_kernel_set_l_limber (xclk, l_limber);

  nc_xcor_kernel_prepare (xclk, env->cosmo);

  return xclk;
}

/*
 * A top-hat kernel on [chi_lower, chi_upper] in the given tier with its adaptive
 * apparatus at the defaults, converged rather than capped, for the checks that are
 * about accuracy. Positive tolerances replace the defaults; sbi is its integrator.
 */
NcXcorKernel *
test_nc_xcor_kquad_tophat_converged (TestNcXcorKQuadEnv *env, NcmSBesselIntegrator *sbi, gdouble chi_lower, gdouble chi_upper, gint l_limber, gdouble reltol, gdouble peak_epsilon)
{
  NcXcorKernel *xclk = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_new_full (env->dist, env->ps,
                                                                                chi_lower, chi_upper, sbi));

  if (reltol > 0.0)
    nc_xcor_kernel_set_reltol (xclk, reltol);

  if (peak_epsilon > 0.0)
    nc_xcor_kernel_set_peak_epsilon (xclk, peak_epsilon);

  nc_xcor_kernel_set_l_limber (xclk, l_limber);
  nc_xcor_kernel_prepare (xclk, env->cosmo);

  return xclk;
}

/* Gauss-Legendre orders of the reference: one to learn the block's peak, then a ladder */
#define TEST_REF_SCALE_ORDER 64
static const guint test_ref_orders[] = {8, 16, 32, 64, 128, 264};

#define TEST_REF_N_ORDERS G_N_ELEMENTS (test_ref_orders)

/*
 * The reference's rules, fixed at construction, and its work space, kept and sized to
 * the block on each evaluation; xclki1 and xclki2 are the pair being evaluated.
 */
struct _TestNcXcorKQuadRef
{
  gsl_integration_glfixed_table *scale_rule;
  gsl_integration_glfixed_table *rules[TEST_REF_N_ORDERS];
  GArray *edges;
  GArray *cells;
  GArray *W1;
  GArray *W2;
  GArray *cur;
  GArray *prev;
  GArray *sum;
  GArray *comp;
  NcXcorKernelIntegrand *xclki1;
  NcXcorKernelIntegrand *xclki2;
};

TestNcXcorKQuadRef *
test_nc_xcor_kquad_ref_new (void)
{
  TestNcXcorKQuadRef *ref = g_new0 (TestNcXcorKQuadRef, 1);
  guint o;

  ref->scale_rule = gsl_integration_glfixed_table_alloc (TEST_REF_SCALE_ORDER);

  for (o = 0; o < TEST_REF_N_ORDERS; o++)
    ref->rules[o] = gsl_integration_glfixed_table_alloc (test_ref_orders[o]);

  ref->edges = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ref->cells = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ref->W1    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ref->W2    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ref->cur   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ref->prev  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ref->sum   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ref->comp  = g_array_new (FALSE, FALSE, sizeof (gdouble));

  return ref;
}

void
test_nc_xcor_kquad_ref_free (TestNcXcorKQuadRef *ref)
{
  guint o;

  for (o = 0; o < TEST_REF_N_ORDERS; o++)
    gsl_integration_glfixed_table_free (ref->rules[o]);

  gsl_integration_glfixed_table_free (ref->scale_rule);
  g_array_unref (ref->edges);
  g_array_unref (ref->cells);
  g_array_unref (ref->W1);
  g_array_unref (ref->W2);
  g_array_unref (ref->cur);
  g_array_unref (ref->prev);
  g_array_unref (ref->sum);
  g_array_unref (ref->comp);
  g_free (ref);
}

#define TEST_REF_AT(a, i) (g_array_index ((a), gdouble, (i)))

static void _test_ref_cells (TestNcXcorKQuadRef *ref);
static gdouble _test_ref_peak (TestNcXcorKQuadRef *ref);
static gdouble _test_ref_cell_converged (TestNcXcorKQuadRef *ref, gdouble a, gdouble b, gdouble atol);
static void _test_ref_add (TestNcXcorKQuadRef *ref);

/*
 * The outer integral 2 / (pi RH^3) int dk k^2 W1 W2 of every component into @cl, over
 * the intersection of the two integrands' ranges, @xclki2 %NULL for an auto spectrum.
 * Cells are the union of both integrands' breakpoints. On each, Gauss-Legendre is
 * raised until it moves by less than 1e-14 of the block's peak: an absolute test,
 * because a cell whose own value nearly cancels can never settle a relative one. Cells
 * are summed with Neumaier's compensation, since a far-separated cross cancels across
 * them. Returns the largest last move of any cell, as a fraction of the peak.
 */
gdouble
test_nc_xcor_kquad_ref_eval (TestNcXcorKQuadRef *ref, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, gdouble RH, NcmVector *cl)
{
  const guint len            = nc_xcor_kernel_integrand_get_len (xclki1);
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (RH));
  gdouble peak, worst_move = 0.0;
  guint i, c;

  g_assert_cmpuint (ncm_vector_len (cl), ==, len);

  ref->xclki1 = xclki1;
  ref->xclki2 = xclki2;

  g_array_set_size (ref->W1, len);
  g_array_set_size (ref->W2, len);
  g_array_set_size (ref->cur, len);
  g_array_set_size (ref->prev, len);
  g_array_set_size (ref->sum, len);
  g_array_set_size (ref->comp, len);

  _test_ref_cells (ref);
  peak = _test_ref_peak (ref);

  for (c = 0; c < len; c++)
  {
    TEST_REF_AT (ref->sum, c)  = 0.0;
    TEST_REF_AT (ref->comp, c) = 0.0;
  }

  if (peak > 0.0)
  {
    for (i = 0; i + 1 < ref->cells->len; i++)
    {
      const gdouble move = _test_ref_cell_converged (ref, TEST_REF_AT (ref->cells, i), TEST_REF_AT (ref->cells, i + 1), 1.0e-14 * peak);

      worst_move = GSL_MAX (worst_move, move / peak);
      _test_ref_add (ref);
    }
  }

  for (c = 0; c < len; c++)
    ncm_vector_set (cl, c, const_factor * (TEST_REF_AT (ref->sum, c) + TEST_REF_AT (ref->comp, c)));

  return worst_move;
}

static void _test_ref_breakpoints (NcXcorKernelIntegrand *xclki, GArray *edges);
static gint _test_ref_cmp_double (gconstpointer a, gconstpointer b);

/*
 * The cell edges into ref->cells: the ends of the integrands' common range and,
 * strictly inside it, the sorted union of their breakpoints. Empty when the ranges do
 * not meet.
 */
static void
_test_ref_cells (TestNcXcorKQuadRef *ref)
{
  gdouble k_min, k_max;
  guint i;

  g_array_set_size (ref->edges, 0);
  g_array_set_size (ref->cells, 0);

  nc_xcor_kernel_integrand_get_range (ref->xclki1, &k_min, &k_max);
  _test_ref_breakpoints (ref->xclki1, ref->edges);

  if (ref->xclki2 != NULL)
  {
    gdouble a2, b2;

    nc_xcor_kernel_integrand_get_range (ref->xclki2, &a2, &b2);
    k_min = GSL_MAX (k_min, a2);
    k_max = GSL_MIN (k_max, b2);
    _test_ref_breakpoints (ref->xclki2, ref->edges);
  }

  if (k_min < k_max)
  {
    g_array_sort (ref->edges, _test_ref_cmp_double);
    g_array_append_val (ref->cells, k_min);

    for (i = 0; i < ref->edges->len; i++)
    {
      const gdouble e = TEST_REF_AT (ref->edges, i);

      if ((e > TEST_REF_AT (ref->cells, ref->cells->len - 1)) && (e < k_max))
        g_array_append_val (ref->cells, e);
    }

    g_array_append_val (ref->cells, k_max);
  }
}

/*
 * Appends to edges the abscissas in k between which every component of xclki is a
 * single polynomial: each component's range ends, and its panel edges, which are in
 * k / s_i, times its scale s_i; or the knots of a spline closure.
 */
static void
_test_ref_breakpoints (NcXcorKernelIntegrand *xclki, GArray *edges)
{
  const guint len      = nc_xcor_kernel_integrand_get_len (xclki);
  const guint n_panels = nc_xcor_kernel_integrand_get_n_panels (xclki);
  NcmVector *knots     = nc_xcor_kernel_integrand_peek_knots (xclki);
  guint c, i;

  for (c = 0; c < len; c++)
  {
    gdouble a, b;

    nc_xcor_kernel_integrand_get_range_comp (xclki, c, &a, &b);
    g_array_append_val (edges, a);
    g_array_append_val (edges, b);
  }

  for (c = 0; c < len; c++)
  {
    const gdouble s = nc_xcor_kernel_integrand_get_scale (xclki, c);

    for (i = 0; i < n_panels; i++)
    {
      NcmMatrix *coeffs;
      gdouble a, b;

      nc_xcor_kernel_integrand_peek_panel (xclki, i, &coeffs, &a, &b);
      a *= s;
      b *= s;
      g_array_append_val (edges, a);
      g_array_append_val (edges, b);
    }
  }

  if (knots != NULL)
  {
    for (i = 0; i < ncm_vector_len (knots); i++)
    {
      const gdouble k = ncm_vector_get (knots, i);

      g_array_append_val (edges, k);
    }
  }
}

static gint
_test_ref_cmp_double (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x > y) - (x < y);
}

static void _test_ref_cell (TestNcXcorKQuadRef *ref, gsl_integration_glfixed_table *rule, gdouble a, gdouble b);

/* The block's peak, max over components of |sum over cells|, by the scale rule alone. */
static gdouble
_test_ref_peak (TestNcXcorKQuadRef *ref)
{
  const guint len = ref->sum->len;
  gdouble peak    = 0.0;
  guint i, c;

  for (c = 0; c < len; c++)
    TEST_REF_AT (ref->sum, c) = 0.0;

  for (i = 0; i + 1 < ref->cells->len; i++)
  {
    _test_ref_cell (ref, ref->scale_rule, TEST_REF_AT (ref->cells, i), TEST_REF_AT (ref->cells, i + 1));

    for (c = 0; c < len; c++)
      TEST_REF_AT (ref->sum, c) += TEST_REF_AT (ref->cur, c);
  }

  for (c = 0; c < len; c++)
    peak = GSL_MAX (peak, fabs (TEST_REF_AT (ref->sum, c)));

  return peak;
}

/*
 * One cell into ref->cur, the Gauss-Legendre order raised over test_ref_orders until two
 * consecutive orders differ by less than atol in every component. Returns that last
 * difference, the order's own when the ladder ran out.
 */
static gdouble
_test_ref_cell_converged (TestNcXcorKQuadRef *ref, gdouble a, gdouble b, gdouble atol)
{
  const guint len = ref->cur->len;
  gdouble move    = GSL_POSINF;
  guint o, c;

  for (o = 0; o < TEST_REF_N_ORDERS; o++)
  {
    _test_ref_cell (ref, ref->rules[o], a, b);

    if (o > 0)
    {
      move = 0.0;

      for (c = 0; c < len; c++)
        move = GSL_MAX (move, fabs (TEST_REF_AT (ref->cur, c) - TEST_REF_AT (ref->prev, c)));

      if (move < atol)
        break;
    }

    memcpy (ref->prev->data, ref->cur->data, sizeof (gdouble) * len);
  }

  return move;
}

/* The sum over one cell of w k^2 W1 W2 by a fixed Gauss-Legendre rule, into ref->cur. */
static void
_test_ref_cell (TestNcXcorKQuadRef *ref, gsl_integration_glfixed_table *rule, gdouble a, gdouble b)
{
  const guint len = ref->cur->len;
  gdouble *W1     = (gdouble *) ref->W1->data;
  gdouble *W2     = (ref->xclki2 != NULL) ? (gdouble *) ref->W2->data : W1;
  guint i, c;

  for (c = 0; c < len; c++)
    TEST_REF_AT (ref->cur, c) = 0.0;

  for (i = 0; i < rule->n; i++)
  {
    gdouble k, w;

    gsl_integration_glfixed_point (a, b, i, &k, &w, rule);
    nc_xcor_kernel_integrand_eval (ref->xclki1, k, W1);

    if (ref->xclki2 != NULL)
      nc_xcor_kernel_integrand_eval (ref->xclki2, k, W2);

    for (c = 0; c < len; c++)
      TEST_REF_AT (ref->cur, c) += w * k * k * W1[c] * W2[c];
  }
}

/* Adds ref->cur to ref->sum with Neumaier's compensation, kept in ref->comp. */
static void
_test_ref_add (TestNcXcorKQuadRef *ref)
{
  const guint len = ref->cur->len;
  guint c;

  for (c = 0; c < len; c++)
  {
    const gdouble s    = TEST_REF_AT (ref->sum, c);
    const gdouble x    = TEST_REF_AT (ref->cur, c);
    const gdouble cand = s + x;

    TEST_REF_AT (ref->comp, c) += (fabs (s) >= fabs (x)) ? (s - cand) + x : (x - cand) + s;
    TEST_REF_AT (ref->sum, c)   = cand;
  }
}

