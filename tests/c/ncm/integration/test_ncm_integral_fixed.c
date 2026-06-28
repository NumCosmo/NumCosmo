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

