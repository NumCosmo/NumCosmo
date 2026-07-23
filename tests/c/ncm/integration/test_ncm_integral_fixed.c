/***************************************************************************
 *            test_ncm_integral_fixed.c
 *
 *  Fri May 23 19:08:43 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2026 <caiolimadeoliveira@pm.me>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_math.h>

typedef struct _TestNcmIntegralFixedConfig
{
  gulong n_nodes;
  gulong rule_n;
  gdouble xl;
  gdouble xu;
} TestNcmIntegralFixedConfig;

static void test_ncm_integral_fixed_get_nodes (gconstpointer pdata);
static void test_ncm_integral_fixed_integ_vec_mult (gconstpointer pdata);
static void test_ncm_integral_fixed_vec_mult_constant (gconstpointer pdata);
static void test_ncm_integral_fixed_calibrate (void);

gint
main (gint argc, gchar *argv[])
{
  static const TestNcmIntegralFixedConfig configs[] = {
    { 5, 3,  0.0, 1.0}, /* odd rule_n, simple interval */
    { 5, 4,  0.0, 1.0}, /* even rule_n */
    {10, 5,  1.0, 3.0}, /* odd rule_n, shifted interval */
    {10, 4,  1.0, 3.0}, /* even rule_n, shifted interval */
    { 3, 7, -1.0, 2.0}, /* odd rule_n >= 4 for exact x^5 */
    { 3, 8, -1.0, 2.0}, /* even rule_n >= 4 */
  };
  const gint n_configs = (gint) (sizeof (configs) / sizeof (configs[0]));
  gint i;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  for (i = 0; i < n_configs; i++)
  {
    gchar *path;

    path = g_strdup_printf ("/ncm/integral_fixed/n%lu_r%lu/get_nodes", configs[i].n_nodes, configs[i].rule_n);
    g_test_add_data_func (path, &configs[i], test_ncm_integral_fixed_get_nodes);
    g_free (path);

    path = g_strdup_printf ("/ncm/integral_fixed/n%lu_r%lu/integ_vec_mult", configs[i].n_nodes, configs[i].rule_n);
    g_test_add_data_func (path, &configs[i], test_ncm_integral_fixed_integ_vec_mult);
    g_free (path);

    path = g_strdup_printf ("/ncm/integral_fixed/n%lu_r%lu/vec_mult_constant", configs[i].n_nodes, configs[i].rule_n);
    g_test_add_data_func (path, &configs[i], test_ncm_integral_fixed_vec_mult_constant);
    g_free (path);
  }

  g_test_add_func ("/ncm/integral_fixed/calibrate", test_ncm_integral_fixed_calibrate);

  g_test_run ();

  return 0;
}

static gdouble
_f_one (gdouble x, gpointer userdata)
{
  NCM_UNUSED (userdata);
  NCM_UNUSED (x);

  return 1.0;
}

static gdouble
_f_x (gdouble x, gpointer userdata)
{
  NCM_UNUSED (userdata);

  return x;
}

static gdouble
_f_x2 (gdouble x, gpointer userdata)
{
  NCM_UNUSED (userdata);

  return x * x;
}

static gdouble
_f_x3 (gdouble x, gpointer userdata)
{
  NCM_UNUSED (userdata);

  return x * x * x;
}

/* Checks that get_nodes returns the same abscissae used by calc_nodes.
 * Strategy: calc_nodes with F(x)=1 gives int_nodes[k] = w_k.
 *           calc_nodes with F(x)=x gives int_nodes[k] = w_k * x_k.
 *           get_nodes returns x_k.
 *           So: (w_k*x_k) / w_k should equal x_k from get_nodes.           */
static void
test_ncm_integral_fixed_get_nodes (gconstpointer pdata)
{
  const TestNcmIntegralFixedConfig *cfg = (const TestNcmIntegralFixedConfig *) pdata;
  const gulong n_nodes                  = cfg->n_nodes;
  const gulong rule_n                   = cfg->rule_n;
  const gdouble xl                      = cfg->xl;
  const gdouble xu                      = cfg->xu;
  const gulong total_n                  = (n_nodes - 1) * rule_n;
  NcmIntegralFixed *intf_w              = ncm_integral_fixed_new (n_nodes, rule_n, xl, xu);
  NcmIntegralFixed *intf_wx             = ncm_integral_fixed_new (n_nodes, rule_n, xl, xu);
  NcmVector *nodes                      = ncm_vector_new (total_n);
  gsl_function F_one                    = {_f_one, NULL};
  gsl_function F_x                      = {_f_x, NULL};
  gulong k;

  ncm_integral_fixed_calc_nodes (intf_w, &F_one);
  ncm_integral_fixed_calc_nodes (intf_wx, &F_x);
  ncm_integral_fixed_get_nodes (intf_wx, nodes);

  for (k = 0; k < total_n; k++)
  {
    const gdouble w_k    = intf_w->int_nodes[k];
    const gdouble wx_k   = intf_wx->int_nodes[k];
    const gdouble x_k    = ncm_vector_get (nodes, k);
    const gdouble x_from = (fabs (w_k) > 0.0) ? wx_k / w_k : 0.0;

    if (fabs (w_k) > 0.0)
      ncm_assert_cmpdouble_e (x_from, ==, x_k, 1.0e-14, 0.0);

    /* All nodes must lie within [xl, xu] */
    g_assert_cmpfloat (x_k, >=, xl);
    g_assert_cmpfloat (x_k, <=, xu);
  }

  /* Analytic check: sum of weights * delta_x/2 == xu - xl */
  {
    const gdouble delta_x = (xu - xl) / (gdouble) (n_nodes - 1);
    gdouble sum_w         = 0.0;
    gulong i;

    for (i = 0; i < total_n; i++)
      sum_w += intf_w->int_nodes[i];

    ncm_assert_cmpdouble_e (sum_w * delta_x * 0.5, ==, xu - xl, 1.0e-13, 0.0);
  }

  ncm_integral_fixed_free (intf_w);
  ncm_integral_fixed_free (intf_wx);
  ncm_vector_free (nodes);
}

/* Checks that integ_vec_mult agrees with integ_mult and with the analytic value.
 * Uses F(x) = x^2 (stored in intf) and G(x) = x^3 (evaluated at nodes).
 * The integral ∫_xl^xu x^2 * x^3 dx = ∫ x^5 dx = (xu^6 - xl^6) / 6.
 * The GL rule is exact for polynomials of degree < 2*rule_n - 1, so for
 * rule_n >= 4 (degree 5 fits) this must agree to ~1e-13.                   */
static void
test_ncm_integral_fixed_integ_vec_mult (gconstpointer pdata)
{
  const TestNcmIntegralFixedConfig *cfg = (const TestNcmIntegralFixedConfig *) pdata;
  const gulong n_nodes                  = cfg->n_nodes;
  const gulong rule_n                   = cfg->rule_n;
  const gdouble xl                      = cfg->xl;
  const gdouble xu                      = cfg->xu;
  const gulong total_n                  = (n_nodes - 1) * rule_n;
  NcmIntegralFixed *intf                = ncm_integral_fixed_new (n_nodes, rule_n, xl, xu);
  NcmVector *nodes                      = ncm_vector_new (total_n);
  NcmVector *g_at_nodes                 = ncm_vector_new (total_n);
  gsl_function F_x2                     = {_f_x2, NULL};
  gsl_function G_x3                     = {_f_x3, NULL};
  gdouble res_mult, res_vec_mult;
  gulong k;

  ncm_integral_fixed_calc_nodes (intf, &F_x2);
  ncm_integral_fixed_get_nodes (intf, nodes);

  for (k = 0; k < total_n; k++)
  {
    const gdouble x_k = ncm_vector_get (nodes, k);

    ncm_vector_set (g_at_nodes, k, x_k * x_k * x_k);
  }

  res_mult     = ncm_integral_fixed_integ_mult (intf, &G_x3);
  res_vec_mult = ncm_integral_fixed_integ_vec_mult (intf, g_at_nodes);

  /* The two APIs must agree to floating-point precision */
  ncm_assert_cmpdouble (res_vec_mult, ==, res_mult);

  /* For rule_n >= 4 the GL rule integrates x^5 exactly */
  if (rule_n >= 4)
  {
    const gdouble analytic = (gsl_pow_6 (xu) - gsl_pow_6 (xl)) / 6.0;

    ncm_assert_cmpdouble_e (res_vec_mult, ==, analytic, 1.0e-12, 0.0);
  }

  ncm_integral_fixed_free (intf);
  ncm_vector_free (nodes);
  ncm_vector_free (g_at_nodes);
}

typedef struct _GaussArg
{
  gdouble mu;
  gdouble sigma;
} GaussArg;

/* A narrow Gaussian weight F(x): the "spiky P(z)" stand-in. */
static gdouble
_f_gauss (gdouble x, gpointer userdata)
{
  const GaussArg *a = (const GaussArg *) userdata;
  const gdouble u   = (x - a->mu) / a->sigma;

  return exp (-0.5 * u * u) / (sqrt (2.0 * M_PI) * a->sigma);
}

/* A smooth, gently tilted probe G(x): the "shape factor" stand-in. */
static gdouble
_g_tilt (gdouble x, gpointer userdata)
{
  NCM_UNUSED (userdata);

  return 1.0 + 0.3 * x;
}

/* High-resolution reference (matches the one calibrate uses internally) and the
 * pass test for a single config, so the test can independently re-check both
 * acceptance criteria. */
static gboolean
_calib_config_passes (gsl_function *F, gsl_function *G, gdouble xl, gdouble xu,
                      gulong n_nodes, gulong rule_n, gdouble reltol,
                      gdouble I_ref, gdouble exact_F)
{
  NcmIntegralFixed *t = ncm_integral_fixed_new (n_nodes, rule_n, xl, xu);
  gdouble I_fg, mass;

  ncm_integral_fixed_calc_nodes (t, F);
  I_fg = ncm_integral_fixed_integ_mult (t, G);
  mass = ncm_integral_fixed_nodes_eval (t);
  ncm_integral_fixed_free (t);

  return (fabs (I_fg - I_ref) / fabs (I_ref) < reltol) &&
         (fabs (mass - exact_F) / fabs (exact_F) < reltol);
}

/* Calibrates a fixed GL rule to integrate F*G of a narrow Gaussian weight times a
 * tilted probe. Checks: (i) the selected config integrates F*G to the target
 * tolerance against an independent ultra-high-res reference; (ii) the missed-mass
 * guard holds; (iii) the total node count stays under the cap; (iv) bisection
 * minimality - the next-smaller n_nodes for the chosen rule fails the criteria. */
static void
test_ncm_integral_fixed_calibrate (void)
{
  const gdouble xl     = -6.0;
  const gdouble xu     = 6.0;
  GaussArg garg        = { 0.4, 0.25 };
  gsl_function F       = {_f_gauss, &garg};
  gsl_function G       = {_g_tilt, NULL};
  const gdouble reltol = 1.0e-7;
  const gulong max_tot = 4000;
  /* Analytic mass of F over [xl, xu]. */
  const gdouble exact_F = 0.5 * (erf ((xu - garg.mu) / (M_SQRT2 * garg.sigma)) -
                                 erf ((xl - garg.mu) / (M_SQRT2 * garg.sigma)));
  guint n_nodes = 0, rule_n = 0;
  NcmIntegralFixed *intf;
  gdouble I_sel, I_truth, mass_sel;

  intf = ncm_integral_fixed_calibrate (&F, &G, xl, xu, reltol, exact_F, max_tot, &n_nodes, &rule_n);

  g_assert_nonnull (intf);
  g_assert_cmpuint (n_nodes, >=, 2);
  g_assert_cmpuint (rule_n, >=, 1);
  g_assert_cmpuint ((n_nodes - 1) * rule_n, <=, max_tot);

  /* Integral of F*G at the selected config. */
  I_sel = ncm_integral_fixed_integ_mult (intf, &G);

  /* Independent ultra-high-resolution reference. */
  {
    NcmIntegralFixed *truth = ncm_integral_fixed_new (1000, 7, xl, xu);

    ncm_integral_fixed_calc_nodes (truth, &F);
    I_truth = ncm_integral_fixed_integ_mult (truth, &G);
    ncm_integral_fixed_free (truth);
  }

  ncm_assert_cmpdouble_e (I_sel, ==, I_truth, 10.0 * reltol, 0.0);

  /* Missed-mass guard: INT F at the selected config matches the analytic mass. */
  mass_sel = ncm_integral_fixed_nodes_eval (intf);
  ncm_assert_cmpdouble_e (mass_sel, ==, exact_F, 10.0 * reltol, 0.0);

  /* Bisection minimality: the selected config passes both criteria, while one
   * fewer node (same rule) does not - so n_nodes is the exact threshold. */
  g_assert_true (_calib_config_passes (&F, &G, xl, xu, n_nodes, rule_n, reltol, I_truth, exact_F));

  if (n_nodes > 2)
    g_assert_false (_calib_config_passes (&F, &G, xl, xu, n_nodes - 1, rule_n, reltol, I_truth, exact_F));

  ncm_integral_fixed_free (intf);
}

/* Trivial sanity: F(x)=1, G(x)=c => integral = c*(xu-xl) */
static void
test_ncm_integral_fixed_vec_mult_constant (gconstpointer pdata)
{
  const TestNcmIntegralFixedConfig *cfg = (const TestNcmIntegralFixedConfig *) pdata;
  const gulong n_nodes                  = cfg->n_nodes;
  const gulong rule_n                   = cfg->rule_n;
  const gdouble xl                      = cfg->xl;
  const gdouble xu                      = cfg->xu;
  const gulong total_n                  = (n_nodes - 1) * rule_n;
  const gdouble c                       = 3.7;
  NcmIntegralFixed *intf                = ncm_integral_fixed_new (n_nodes, rule_n, xl, xu);
  NcmVector *g_at_nodes                 = ncm_vector_new (total_n);
  gsl_function F_one                    = {_f_one, NULL};
  gulong k;

  ncm_integral_fixed_calc_nodes (intf, &F_one);

  for (k = 0; k < total_n; k++)
    ncm_vector_set (g_at_nodes, k, c);

  ncm_assert_cmpdouble_e (ncm_integral_fixed_integ_vec_mult (intf, g_at_nodes), ==, c * (xu - xl), 1.0e-14, 0.0);

  ncm_integral_fixed_free (intf);
  ncm_vector_free (g_at_nodes);
}

