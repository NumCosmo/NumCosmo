/***************************************************************************
 *            test_ncm_sbessel_integrator_levin.c
 *
 *  Wed October 14 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_sbessel_integrator_levin.c
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
 * NcmSBesselIntegratorLevin checks that need no outside reference: entry points against
 * each other, closed forms, the panel records, the constraint counters and the aborts.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include "test_ncm_sbessel_arb_cases.h"

void test_ncm_sbessel_integrator_levin_arb_constant (void);
void test_ncm_sbessel_integrator_levin_arb_linear (void);
void test_ncm_sbessel_integrator_levin_arb_quadratic (void);
void test_ncm_sbessel_integrator_levin_arb_rational (void);
void test_ncm_sbessel_integrator_levin_arb_moving_edge (void);
void test_ncm_sbessel_integrator_levin_arb_accepted_edge (void);
void test_ncm_sbessel_integrator_levin_arb_rejected_edge (void);
void test_ncm_sbessel_integrator_levin_arb_set_reltol (void);
void test_ncm_sbessel_integrator_levin_arb_repeated_set_reltol (void);
void test_ncm_sbessel_integrator_levin_arb_ell_zero (void);
void test_ncm_sbessel_integrator_levin_arb_k_scaling (void);
void test_ncm_sbessel_integrator_levin_arb_narrow (void);
void test_ncm_sbessel_integrator_levin_arb_high_ell (void);
void test_ncm_sbessel_integrator_levin_arb_oscillatory (void);
void test_ncm_sbessel_integrator_levin_arb_evanescent (void);
void test_ncm_sbessel_integrator_levin_truth (gconstpointer data);
void test_ncm_sbessel_integrator_levin_ell_range_consistency (void);
void test_ncm_sbessel_integrator_levin_zero_bessel_batch (void);
void test_ncm_sbessel_integrator_levin_ell_range_caching (void);
void test_ncm_sbessel_integrator_levin_block_size (gconstpointer data);
void test_ncm_sbessel_integrator_levin_loose_tolerances_agree (void);
void test_ncm_sbessel_integrator_levin_accuracy_monotonic (void);
void test_ncm_sbessel_integrator_levin_record_off_by_default (void);
void test_ncm_sbessel_integrator_levin_record_sum (void);
void test_ncm_sbessel_integrator_levin_record_tile (void);
void test_ncm_sbessel_integrator_levin_record_cleared (void);
void test_ncm_sbessel_integrator_levin_record_disable (void);
void test_ncm_sbessel_integrator_levin_deriv_constant (void);
void test_ncm_sbessel_integrator_levin_deriv_batched (void);
void test_ncm_sbessel_integrator_levin_deriv_zero (void);
void test_ncm_sbessel_integrator_levin_deriv_recurrence (void);
void test_ncm_sbessel_integrator_levin_tau_panel_qualifies (void);
void test_ncm_sbessel_integrator_levin_tau_smooth (void);
void test_ncm_sbessel_integrator_levin_tau_falls_back (gconstpointer data);
void test_ncm_sbessel_integrator_levin_tau_rule_off (void);
void test_ncm_sbessel_integrator_levin_turning_default (void);
void test_ncm_sbessel_integrator_levin_turning_round_trip (void);
void test_ncm_sbessel_integrator_levin_turning_removes_fallbacks (gconstpointer data);
void test_ncm_sbessel_integrator_levin_turning_same_integral (void);
void test_ncm_sbessel_integrator_levin_turning_low_ell (void);
void test_ncm_sbessel_integrator_levin_turning_not_locked (gconstpointer data);
void test_ncm_sbessel_integrator_levin_turning_rule_off (void);
void test_ncm_sbessel_integrator_levin_shared_second (void);
void test_ncm_sbessel_integrator_levin_shared_identical (void);
void test_ncm_sbessel_integrator_levin_shared_turning_row (void);
void test_ncm_sbessel_integrator_levin_traps (void);
void test_ncm_sbessel_integrator_levin_invalid_reltol_zero (void);
void test_ncm_sbessel_integrator_levin_invalid_reltol_negative (void);
void test_ncm_sbessel_integrator_levin_invalid_cheb_reltol_zero (void);
void test_ncm_sbessel_integrator_levin_invalid_cheb_reltol_negative (void);
void test_ncm_sbessel_integrator_levin_invalid_deriv (void);

static const gchar *test_truth_tables[] = {"gauss", "rational"};
static const guint test_block_sizes[]   = {1, 4, 8, 16};
static const gdouble test_tau_widths[]  = {30.0, 100.0};
static const guint test_turning_ells[]  = {500, 1000};

gint
main (gint argc, gchar *argv[])
{
  guint i;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);

  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/constant", test_ncm_sbessel_integrator_levin_arb_constant);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/linear", test_ncm_sbessel_integrator_levin_arb_linear);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/quadratic", test_ncm_sbessel_integrator_levin_arb_quadratic);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/rational", test_ncm_sbessel_integrator_levin_arb_rational);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/moving_edge", test_ncm_sbessel_integrator_levin_arb_moving_edge);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/accepted_edge", test_ncm_sbessel_integrator_levin_arb_accepted_edge);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/rejected_edge", test_ncm_sbessel_integrator_levin_arb_rejected_edge);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/set_reltol", test_ncm_sbessel_integrator_levin_arb_set_reltol);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/repeated_set_reltol", test_ncm_sbessel_integrator_levin_arb_repeated_set_reltol);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/ell_zero", test_ncm_sbessel_integrator_levin_arb_ell_zero);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/k_scaling", test_ncm_sbessel_integrator_levin_arb_k_scaling);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/narrow", test_ncm_sbessel_integrator_levin_arb_narrow);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/high_ell", test_ncm_sbessel_integrator_levin_arb_high_ell);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/oscillatory", test_ncm_sbessel_integrator_levin_arb_oscillatory);
  g_test_add_func ("/ncm/sbessel_integrator_levin/arb/evanescent", test_ncm_sbessel_integrator_levin_arb_evanescent);

  for (i = 0; i < G_N_ELEMENTS (test_truth_tables); i++)
  {
    gchar *path = g_strdup_printf ("/ncm/sbessel_integrator_levin/truth/%s", test_truth_tables[i]);

    g_test_add_data_func (path, test_truth_tables[i], test_ncm_sbessel_integrator_levin_truth);
    g_free (path);
  }

  g_test_add_func ("/ncm/sbessel_integrator_levin/ell_range_consistency", test_ncm_sbessel_integrator_levin_ell_range_consistency);
  g_test_add_func ("/ncm/sbessel_integrator_levin/zero_bessel_batch", test_ncm_sbessel_integrator_levin_zero_bessel_batch);
  g_test_add_func ("/ncm/sbessel_integrator_levin/ell_range_caching", test_ncm_sbessel_integrator_levin_ell_range_caching);

  for (i = 0; i < G_N_ELEMENTS (test_block_sizes); i++)
  {
    gchar *path = g_strdup_printf ("/ncm/sbessel_integrator_levin/block_size/%u", test_block_sizes[i]);

    g_test_add_data_func (path, &test_block_sizes[i], test_ncm_sbessel_integrator_levin_block_size);
    g_free (path);
  }

  g_test_add_func ("/ncm/sbessel_integrator_levin/resolution_floor/loose_tolerances_agree", test_ncm_sbessel_integrator_levin_loose_tolerances_agree);
  g_test_add_func ("/ncm/sbessel_integrator_levin/resolution_floor/accuracy_monotonic", test_ncm_sbessel_integrator_levin_accuracy_monotonic);

  g_test_add_func ("/ncm/sbessel_integrator_levin/record/off_by_default", test_ncm_sbessel_integrator_levin_record_off_by_default);
  g_test_add_func ("/ncm/sbessel_integrator_levin/record/sum", test_ncm_sbessel_integrator_levin_record_sum);
  g_test_add_func ("/ncm/sbessel_integrator_levin/record/tile", test_ncm_sbessel_integrator_levin_record_tile);
  g_test_add_func ("/ncm/sbessel_integrator_levin/record/cleared", test_ncm_sbessel_integrator_levin_record_cleared);
  g_test_add_func ("/ncm/sbessel_integrator_levin/record/disable", test_ncm_sbessel_integrator_levin_record_disable);

  g_test_add_func ("/ncm/sbessel_integrator_levin/deriv/constant", test_ncm_sbessel_integrator_levin_deriv_constant);
  g_test_add_func ("/ncm/sbessel_integrator_levin/deriv/batched", test_ncm_sbessel_integrator_levin_deriv_batched);
  g_test_add_func ("/ncm/sbessel_integrator_levin/deriv/zero", test_ncm_sbessel_integrator_levin_deriv_zero);
  g_test_add_func ("/ncm/sbessel_integrator_levin/deriv/recurrence", test_ncm_sbessel_integrator_levin_deriv_recurrence);

  g_test_add_func ("/ncm/sbessel_integrator_levin/tau/panel_qualifies", test_ncm_sbessel_integrator_levin_tau_panel_qualifies);
  g_test_add_func ("/ncm/sbessel_integrator_levin/tau/smooth", test_ncm_sbessel_integrator_levin_tau_smooth);

  for (i = 0; i < G_N_ELEMENTS (test_tau_widths); i++)
  {
    gchar *path = g_strdup_printf ("/ncm/sbessel_integrator_levin/tau/falls_back/%g", test_tau_widths[i]);

    g_test_add_data_func (path, &test_tau_widths[i], test_ncm_sbessel_integrator_levin_tau_falls_back);
    g_free (path);
  }

  g_test_add_func ("/ncm/sbessel_integrator_levin/tau/rule_off", test_ncm_sbessel_integrator_levin_tau_rule_off);

  g_test_add_func ("/ncm/sbessel_integrator_levin/turning/default", test_ncm_sbessel_integrator_levin_turning_default);
  g_test_add_func ("/ncm/sbessel_integrator_levin/turning/round_trip", test_ncm_sbessel_integrator_levin_turning_round_trip);

  for (i = 0; i < G_N_ELEMENTS (test_turning_ells); i++)
  {
    gchar *path = g_strdup_printf ("/ncm/sbessel_integrator_levin/turning/removes_fallbacks/%u", test_turning_ells[i]);

    g_test_add_data_func (path, &test_turning_ells[i], test_ncm_sbessel_integrator_levin_turning_removes_fallbacks);
    g_free (path);
  }

  g_test_add_func ("/ncm/sbessel_integrator_levin/turning/same_integral", test_ncm_sbessel_integrator_levin_turning_same_integral);
  g_test_add_func ("/ncm/sbessel_integrator_levin/turning/low_ell", test_ncm_sbessel_integrator_levin_turning_low_ell);

  for (i = 0; i < G_N_ELEMENTS (test_turning_ells); i++)
  {
    gchar *path = g_strdup_printf ("/ncm/sbessel_integrator_levin/turning/not_locked/%u", test_turning_ells[i]);

    g_test_add_data_func (path, &test_turning_ells[i], test_ncm_sbessel_integrator_levin_turning_not_locked);
    g_free (path);
  }

  g_test_add_func ("/ncm/sbessel_integrator_levin/turning/rule_off", test_ncm_sbessel_integrator_levin_turning_rule_off);

  g_test_add_func ("/ncm/sbessel_integrator_levin/shared/second", test_ncm_sbessel_integrator_levin_shared_second);
  g_test_add_func ("/ncm/sbessel_integrator_levin/shared/identical", test_ncm_sbessel_integrator_levin_shared_identical);
  g_test_add_func ("/ncm/sbessel_integrator_levin/shared/turning_row", test_ncm_sbessel_integrator_levin_shared_turning_row);

  g_test_add_func ("/ncm/sbessel_integrator_levin/traps", test_ncm_sbessel_integrator_levin_traps);
  g_test_add_func ("/ncm/sbessel_integrator_levin/invalid/reltol_zero/subprocess", test_ncm_sbessel_integrator_levin_invalid_reltol_zero);
  g_test_add_func ("/ncm/sbessel_integrator_levin/invalid/reltol_negative/subprocess", test_ncm_sbessel_integrator_levin_invalid_reltol_negative);
  g_test_add_func ("/ncm/sbessel_integrator_levin/invalid/cheb_reltol_zero/subprocess", test_ncm_sbessel_integrator_levin_invalid_cheb_reltol_zero);
  g_test_add_func ("/ncm/sbessel_integrator_levin/invalid/cheb_reltol_negative/subprocess", test_ncm_sbessel_integrator_levin_invalid_cheb_reltol_negative);
  g_test_add_func ("/ncm/sbessel_integrator_levin/invalid/deriv/subprocess", test_ncm_sbessel_integrator_levin_invalid_deriv);

  g_test_run ();
}

static gdouble
_test_absmax (NcmVector *v)
{
  gdouble absmin, absmax;

  ncm_vector_get_absminmax (v, &absmin, &absmax);

  return absmax;
}

/* |actual - desired| <= atol + rtol |desired|, the numpy assert_allclose criterion */
static void
_test_assert_close (gdouble actual, gdouble desired, gdouble rtol, gdouble atol)
{
  g_assert_true (gsl_finite (actual));
  g_assert_cmpfloat (fabs (actual - desired), <=, atol + rtol * fabs (desired));
}

static gdouble
_test_exp_decay (gpointer user_data, gdouble x, gdouble k)
{
  const gdouble *rate = (const gdouble *) user_data;

  return exp (-(*rate) * x);
}

static gdouble
_test_counted_one (gpointer user_data, gdouble x, gdouble k)
{
  guint *calls = (guint *) user_data;

  (*calls)++;

  return 1.0;
}

static gdouble
_test_one (gpointer user_data, gdouble x, gdouble k)
{
  return 1.0;
}

typedef struct _TestGauss
{
  gdouble center;
  gdouble sigma;
  gdouble amplitude;
} TestGauss;

static gdouble
_test_gauss (gpointer user_data, gdouble x, gdouble k)
{
  const TestGauss *g = (const TestGauss *) user_data;
  const gdouble z    = (x - g->center) / g->sigma;

  return g->amplitude * exp (-0.5 * z * z);
}

/* Certified integrals of make_sbessel_arb_cases.py levin, see test_ncm_sbessel_arb_cases.h */
static const TestArbCell test_arb_cells[] = {
  /* BEGIN make_sbessel_arb_cases.py levin */
  {"constant", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 0, 1.0, 0.6021586306762569},
  {"constant", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 1, 1.0, 0.7958237222715151},
  {"constant", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 5, 1.0, 0.48260284581290397},
  {"constant", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 10, 1.0, 0.34015366131471403},
  {"constant", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 20, 1.0, 0.07200515204667218},
  {"linear", TEST_ARB_POWER, 1.0, 20.0, 1.0, 0.0, 0, FALSE, 0, 1.0, 0.13222024405474772},
  {"linear", TEST_ARB_POWER, 1.0, 20.0, 1.0, 0.0, 0, FALSE, 1, 1.0, 0.5306843647565257},
  {"linear", TEST_ARB_POWER, 1.0, 20.0, 1.0, 0.0, 0, FALSE, 5, 1.0, 1.9425860132872468},
  {"linear", TEST_ARB_POWER, 1.0, 20.0, 1.0, 0.0, 0, FALSE, 10, 1.0, 3.1777343624297174},
  {"linear", TEST_ARB_POWER, 1.0, 20.0, 1.0, 0.0, 0, FALSE, 20, 1.0, 1.330664477477039},
  {"quadratic", TEST_ARB_POWER, 1.0, 20.0, 2.0, 0.0, 0, FALSE, 0, 1.0, -7.549864664479969},
  {"quadratic", TEST_ARB_POWER, 1.0, 20.0, 2.0, 0.0, 0, FALSE, 1, 1.0, -17.152993541635162},
  {"quadratic", TEST_ARB_POWER, 1.0, 20.0, 2.0, 0.0, 0, FALSE, 5, 1.0, -3.7098496929237466},
  {"quadratic", TEST_ARB_POWER, 1.0, 20.0, 2.0, 0.0, 0, FALSE, 10, 1.0, 25.794284875856388},
  {"rational", TEST_ARB_RATIONAL, 1.0, 20.0, 10.5, 3.0, 0, FALSE, 0, 1.0, -12.391811110003745},
  {"rational", TEST_ARB_RATIONAL, 1.0, 20.0, 10.5, 3.0, 0, FALSE, 1, 1.0, 3.2008147246960106},
  {"rational", TEST_ARB_RATIONAL, 1.0, 20.0, 10.5, 3.0, 0, FALSE, 5, 1.0, -15.439707699282732},
  {"rational", TEST_ARB_RATIONAL, 1.0, 20.0, 10.5, 3.0, 0, FALSE, 10, 1.0, 29.801551483932744},
  {"rational", TEST_ARB_RATIONAL, 1.0, 20.0, 10.5, 3.0, 0, FALSE, 20, 1.0, 0.08489097486799253},
  {"rational", TEST_ARB_RATIONAL, 1.0, 20.0, 10.5, 3.0, 0, FALSE, 200, 20.0, 0.4148178074571125},
  {"moving_edge_0", TEST_ARB_EXP, 5.0, 50.0, 33.333333333333336, 0.0, 0, FALSE, 0, 1.0, 0.009119978431775066},
  {"moving_edge_0", TEST_ARB_EXP, 5.0, 50.0, 33.333333333333336, 0.0, 0, FALSE, 1, 1.0, -0.16417347888847092},
  {"moving_edge_0", TEST_ARB_EXP, 5.0, 50.0, 33.333333333333336, 0.0, 0, FALSE, 2, 1.0, -0.10432400526491802},
  {"moving_edge_0", TEST_ARB_EXP, 5.0, 50.0, 33.333333333333336, 0.0, 0, FALSE, 3, 1.0, 0.08752243602152396},
  {"moving_edge_0", TEST_ARB_EXP, 5.0, 50.0, 33.333333333333336, 0.0, 0, FALSE, 4, 1.0, 0.2555904588029378},
  {"moving_edge_0", TEST_ARB_EXP, 5.0, 50.0, 33.333333333333336, 0.0, 0, FALSE, 5, 1.0, 0.34648365298352607},
  {"moving_edge_0", TEST_ARB_EXP, 5.0, 50.0, 33.333333333333336, 0.0, 0, FALSE, 6, 1.0, 0.3706808798675939},
  {"moving_edge_0", TEST_ARB_EXP, 5.0, 50.0, 33.333333333333336, 0.0, 0, FALSE, 7, 1.0, 0.35429189497075153},
  {"moving_edge_1", TEST_ARB_EXP, 6.7, 43.0, 33.333333333333336, 0.0, 0, FALSE, 0, 1.0, 0.11132862017108842},
  {"moving_edge_1", TEST_ARB_EXP, 6.7, 43.0, 33.333333333333336, 0.0, 0, FALSE, 1, 1.0, 0.0514075231243722},
  {"moving_edge_1", TEST_ARB_EXP, 6.7, 43.0, 33.333333333333336, 0.0, 0, FALSE, 2, 1.0, -0.09750506779987338},
  {"moving_edge_1", TEST_ARB_EXP, 6.7, 43.0, 33.333333333333336, 0.0, 0, FALSE, 3, 1.0, -0.12946822430690377},
  {"moving_edge_1", TEST_ARB_EXP, 6.7, 43.0, 33.333333333333336, 0.0, 0, FALSE, 4, 1.0, -0.016869271809440943},
  {"moving_edge_1", TEST_ARB_EXP, 6.7, 43.0, 33.333333333333336, 0.0, 0, FALSE, 5, 1.0, 0.13936825342229892},
  {"moving_edge_1", TEST_ARB_EXP, 6.7, 43.0, 33.333333333333336, 0.0, 0, FALSE, 6, 1.0, 0.24850488506975554},
  {"moving_edge_1", TEST_ARB_EXP, 6.7, 43.0, 33.333333333333336, 0.0, 0, FALSE, 7, 1.0, 0.29197853219169384},
  {"moving_edge_2", TEST_ARB_EXP, 9.5, 32.0, 33.333333333333336, 0.0, 0, FALSE, 0, 1.0, -0.0879525949168579},
  {"moving_edge_2", TEST_ARB_EXP, 9.5, 32.0, 33.333333333333336, 0.0, 0, FALSE, 1, 1.0, -0.009908381689475356},
  {"moving_edge_2", TEST_ARB_EXP, 9.5, 32.0, 33.333333333333336, 0.0, 0, FALSE, 2, 1.0, 0.08859703860058944},
  {"moving_edge_2", TEST_ARB_EXP, 9.5, 32.0, 33.333333333333336, 0.0, 0, FALSE, 3, 1.0, 0.05261904530206776},
  {"moving_edge_2", TEST_ARB_EXP, 9.5, 32.0, 33.333333333333336, 0.0, 0, FALSE, 4, 1.0, -0.06015245320021174},
  {"moving_edge_2", TEST_ARB_EXP, 9.5, 32.0, 33.333333333333336, 0.0, 0, FALSE, 5, 1.0, -0.10923788735438252},
  {"moving_edge_2", TEST_ARB_EXP, 9.5, 32.0, 33.333333333333336, 0.0, 0, FALSE, 6, 1.0, -0.052066984908389374},
  {"moving_edge_2", TEST_ARB_EXP, 9.5, 32.0, 33.333333333333336, 0.0, 0, FALSE, 7, 1.0, 0.05042494149045293},
  {"moving_edge_3", TEST_ARB_EXP, 10.5, 96.0, 33.333333333333336, 0.0, 0, FALSE, 0, 1.0, -0.0395063336597935},
  {"moving_edge_3", TEST_ARB_EXP, 10.5, 96.0, 33.333333333333336, 0.0, 0, FALSE, 1, 1.0, -0.06053206247432944},
  {"moving_edge_3", TEST_ARB_EXP, 10.5, 96.0, 33.333333333333336, 0.0, 0, FALSE, 2, 1.0, 0.023646311006344813},
  {"moving_edge_3", TEST_ARB_EXP, 10.5, 96.0, 33.333333333333336, 0.0, 0, FALSE, 3, 1.0, 0.07428521787473621},
  {"moving_edge_3", TEST_ARB_EXP, 10.5, 96.0, 33.333333333333336, 0.0, 0, FALSE, 4, 1.0, 0.0238668133674281},
  {"moving_edge_3", TEST_ARB_EXP, 10.5, 96.0, 33.333333333333336, 0.0, 0, FALSE, 5, 1.0, -0.06022051335629506},
  {"moving_edge_3", TEST_ARB_EXP, 10.5, 96.0, 33.333333333333336, 0.0, 0, FALSE, 6, 1.0, -0.08985303866903385},
  {"moving_edge_3", TEST_ARB_EXP, 10.5, 96.0, 33.333333333333336, 0.0, 0, FALSE, 7, 1.0, -0.04413403157505907},
  {"accepted_edge", TEST_ARB_TOPHAT, 1.0, 2.0, 0.0, 0.0, 0, FALSE, 0, 1.0, 0.6593299064355118},
  {"accepted_edge", TEST_ARB_TOPHAT, 1.0, 2.0, 0.0, 0.0, 0, FALSE, 1, 1.0, 0.38682227139505565},
  {"accepted_edge", TEST_ARB_TOPHAT, 1.0, 2.0, 0.0, 0.0, 0, FALSE, 2, 1.0, 0.12832130915740367},
  {"accepted_edge", TEST_ARB_TOPHAT, 1.0, 2.0, 0.0, 0.0, 0, FALSE, 3, 1.0, 0.030526685853749243},
  {"accepted_edge", TEST_ARB_TOPHAT, 1.0, 2.0, 0.0, 0.0, 0, FALSE, 4, 1.0, 0.0057388279129687095},
  {"accepted_edge", TEST_ARB_TOPHAT, 1.0, 2.0, 0.0, 0.0, 0, FALSE, 5, 1.0, 0.0008982701648965717},
  {"accepted_edge", TEST_ARB_TOPHAT, 1.0, 2.0, 0.0, 0.0, 0, FALSE, 6, 1.0, 0.00012090747281368287},
  {"accepted_edge", TEST_ARB_TOPHAT, 1.0, 2.0, 0.0, 0.0, 0, FALSE, 7, 1.0, 1.4304072407906955e-05},
  {"rejected_edge", TEST_ARB_COS, 1.0, 2.0, 500.0, 0.0, 0, FALSE, 0, 1.0, 0.0015370690742047299},
  {"set_reltol", TEST_ARB_EXP, 100.0, 10000.0, 1000.0, 0.0, 0, FALSE, 0, 1.0, 0.007750496513501052},
  {"set_reltol", TEST_ARB_EXP, 100.0, 10000.0, 1000.0, 0.0, 0, FALSE, 1, 1.0, -0.004589534901865097},
  {"set_reltol", TEST_ARB_EXP, 100.0, 10000.0, 1000.0, 0.0, 0, FALSE, 2, 1.0, -0.007890476440616185},
  {"set_reltol", TEST_ARB_EXP, 100.0, 10000.0, 1000.0, 0.0, 0, FALSE, 3, 1.0, 0.004197348115723649},
  {"set_reltol", TEST_ARB_EXP, 100.0, 10000.0, 1000.0, 0.0, 0, FALSE, 4, 1.0, 0.008189812011263732},
  {"set_reltol", TEST_ARB_EXP, 100.0, 10000.0, 1000.0, 0.0, 0, FALSE, 5, 1.0, -0.0034639808864253},
  {"set_reltol", TEST_ARB_EXP, 100.0, 10000.0, 1000.0, 0.0, 0, FALSE, 6, 1.0, -0.008579950301882813},
  {"set_reltol", TEST_ARB_EXP, 100.0, 10000.0, 1000.0, 0.0, 0, FALSE, 7, 1.0, 0.0023527834949275653},
  {"repeated_set_reltol", TEST_ARB_EXP, 12.0, 47.0, 100.0, 0.0, 0, FALSE, 0, 1.0, 0.07102420363141924},
  {"repeated_set_reltol", TEST_ARB_EXP, 12.0, 47.0, 100.0, 0.0, 0, FALSE, 1, 1.0, -0.04201163150805379},
  {"repeated_set_reltol", TEST_ARB_EXP, 12.0, 47.0, 100.0, 0.0, 0, FALSE, 2, 1.0, -0.08221534646397724},
  {"repeated_set_reltol", TEST_ARB_EXP, 12.0, 47.0, 100.0, 0.0, 0, FALSE, 3, 1.0, 0.013426384614975292},
  {"repeated_set_reltol", TEST_ARB_EXP, 12.0, 47.0, 100.0, 0.0, 0, FALSE, 4, 1.0, 0.09328165963961846},
  {"repeated_set_reltol", TEST_ARB_EXP, 12.0, 47.0, 100.0, 0.0, 0, FALSE, 5, 1.0, 0.048487996142564084},
  {"ell_zero", TEST_ARB_TOPHAT, 1.0, 10.0, 0.0, 0.0, 0, FALSE, 0, 1.0, 0.712264523851691},
  {"k_scaling", TEST_ARB_TOPHAT, 1.0, 10.0, 0.0, 0.0, 0, FALSE, 5, 0.1, 0.0001557694944990179},
  {"k_scaling", TEST_ARB_TOPHAT, 1.0, 10.0, 0.0, 0.0, 0, FALSE, 5, 1.0, 0.6484700039890495},
  {"k_scaling", TEST_ARB_TOPHAT, 1.0, 10.0, 0.0, 0.0, 0, FALSE, 5, 10.0, -0.01113331884326862},
  {"k_scaling", TEST_ARB_TOPHAT, 1.0, 10.0, 0.0, 0.0, 0, FALSE, 5, 100.0, -4.653748260250706e-05},
  {"narrow", TEST_ARB_POWER, 5.0, 5.001, 1.0, 0.0, 0, FALSE, 5, 1.0, 0.0005342563413119415},
  {"high_ell", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 100, 100.0, 0.0008907936235087944},
  {"high_ell", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 200, 200.0, 0.0003121418486470063},
  {"high_ell", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 300, 300.0, 0.00016794746443265403},
  {"oscillatory", TEST_ARB_EXP, 1.0, 100.0, 100.0, 0.0, 0, FALSE, 10, 500.0, 3.2810031358293054e-06},
  {"evanescent", TEST_ARB_GAUSS, 1.0, 8.0, 4.0, 1.2, 2, FALSE, 30, 0.5, 8.723201535126188e-26},
  /* END make_sbessel_arb_cases.py levin */
};

static const TestArbCell *
_test_arb_case (const gchar *name, guint *n)
{
  return test_arb_cells_find (test_arb_cells, G_N_ELEMENTS (test_arb_cells), name, n);
}

/* The integrand of one case, counting its calls and the range of x it was evaluated at */
typedef struct _TestArbEval
{
  const TestArbCell *cell;
  guint calls;
  gdouble x_min;
  gdouble x_max;
} TestArbEval;

static TestArbEval
_test_arb_eval_new (const TestArbCell *cell)
{
  TestArbEval ev = {cell, 0, GSL_POSINF, GSL_NEGINF};

  return ev;
}

static gdouble
_test_arb_integrand (gpointer user_data, gdouble x, gdouble k)
{
  TestArbEval *ev = (TestArbEval *) user_data;

  ev->calls++;
  ev->x_min = GSL_MIN (ev->x_min, x);
  ev->x_max = GSL_MAX (ev->x_max, x);

  return test_arb_cell_f (ev->cell, x);
}

/*
 * Each row of case @name with a fresh single-multipole integrator at the defaults,
 * held to |I - I_arb| <= atol + rtol |I_arb|. The default reltol is 1e-13 and
 * cheb-reltol 1e-8, so rtol 1e-8 is the tolerance the integrator is asked for.
 */
static void
_test_arb_check_singles (const gchar *name, gdouble rtol, gdouble atol)
{
  guint n, i;
  const TestArbCell *c = _test_arb_case (name, &n);

  for (i = 0; i < n; i++)
  {
    NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (c[i].ell, c[i].ell));
    TestArbEval ev            = _test_arb_eval_new (&c[i]);
    NcmVector *res            = ncm_vector_new (1);

    ncm_sbessel_integrator_integrate_deriv (sbi, &_test_arb_integrand, c[i].a, c[i].b, c[i].k, c[i].deriv, res, &ev);
    _test_assert_close (ncm_vector_get (res, 0), c[i].value, rtol, atol);

    ncm_vector_free (res);
    ncm_sbessel_integrator_free (sbi);
  }
}

void
test_ncm_sbessel_integrator_levin_arb_constant (void)
{
  _test_arb_check_singles ("constant", 1.0e-8, 1.0e-20);
}

void
test_ncm_sbessel_integrator_levin_arb_linear (void)
{
  _test_arb_check_singles ("linear", 1.0e-8, 1.0e-20);
}

void
test_ncm_sbessel_integrator_levin_arb_quadratic (void)
{
  _test_arb_check_singles ("quadratic", 1.0e-8, 1.0e-20);
}

void
test_ncm_sbessel_integrator_levin_arb_rational (void)
{
  _test_arb_check_singles ("rational", 1.0e-8, 1.0e-20);
}

/*
 * Moving edge cells are reused without sampling outside [a, b]: four ranges of
 * exp (-x / s) on one integrator with 9 knots on [1, 1e4], ell 0..7.
 */
void
test_ncm_sbessel_integrator_levin_arb_moving_edge (void)
{
  const gchar *cases[]             = {"moving_edge_0", "moving_edge_1", "moving_edge_2", "moving_edge_3"};
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new_full (0, 7, 1.0, 1.0e4, 9, 500, 1.0e-9, 2, 1.0e-10);
  NcmVector *res                   = ncm_vector_new (8);
  guint i, j, n;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    const TestArbCell *c = _test_arb_case (cases[i], &n);
    TestArbEval ev       = _test_arb_eval_new (c);

    g_assert_cmpuint (n, ==, 8);
    ncm_sbessel_integrator_integrate (NCM_SBESSEL_INTEGRATOR (sbilv), &_test_arb_integrand, c->a, c->b, 1.0, res, &ev);

    g_assert_cmpfloat (ev.x_min, >=, nextafter (c->a, GSL_NEGINF));
    g_assert_cmpfloat (ev.x_max, <=, nextafter (c->b, GSL_POSINF));

    for (j = 0; j < n; j++)
      _test_assert_close (ncm_vector_get (res, c[j].ell), c[j].value, 1.0e-8, 1.0e-13);
  }

  ncm_vector_free (res);
  ncm_sbessel_integrator_levin_free (sbilv);
}

/*
 * A smooth edge is fitted once: the first knot is the left end of [1, 2], leaving one
 * right edge, so the kernel is called for one fit of 9 nodes plus the 5 probes that
 * decide whether the edge is a dead-junction cell.
 */
void
test_ncm_sbessel_integrator_levin_arb_accepted_edge (void)
{
  guint n, j;
  const TestArbCell *c             = _test_arb_case ("accepted_edge", &n);
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new_full (0, 7, 1.0, 1.0e4, 9, 500, 1.0e-10, 2, 1.0e-11);
  NcmVector *res                   = ncm_vector_new (8);
  TestArbEval ev                   = _test_arb_eval_new (c);

  g_assert_cmpuint (n, ==, 8);
  ncm_sbessel_integrator_integrate (NCM_SBESSEL_INTEGRATOR (sbilv), &_test_arb_integrand, c->a, c->b, 1.0, res, &ev);

  g_assert_cmpuint (ev.calls, ==, 9 + 5);

  for (j = 0; j < n; j++)
    _test_assert_close (ncm_vector_get (res, c[j].ell), c[j].value, 1.0e-9, 1.0e-14);

  ncm_vector_free (res);
  ncm_sbessel_integrator_levin_free (sbilv);
}

/* A rejected extrapolation reuses the fit of 1025 nodes: only the 5 dead-junction probes are added */
void
test_ncm_sbessel_integrator_levin_arb_rejected_edge (void)
{
  guint n;
  const TestArbCell *c             = _test_arb_case ("rejected_edge", &n);
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new_full (0, 7, 1.0, 1.0e4, 9, 500, 1.0e-10, 2, 1.0e-11);
  NcmVector *res                   = ncm_vector_new (8);
  TestArbEval ev                   = _test_arb_eval_new (c);

  g_assert_cmpuint (n, ==, 1);
  g_assert_cmpuint (c->ell, ==, 0);
  ncm_sbessel_integrator_integrate (NCM_SBESSEL_INTEGRATOR (sbilv), &_test_arb_integrand, c->a, c->b, 1.0, res, &ev);

  g_assert_cmpuint (ev.calls, ==, 1025 + 5);
  _test_assert_close (ncm_vector_get (res, 0), c->value, 1.0e-8, 1.0e-13);

  ncm_vector_free (res);
  ncm_sbessel_integrator_levin_free (sbilv);
}

/* The largest |I - I_arb| over the rows of one block, divided by the largest |I_arb| */
static gdouble
_test_arb_block_error (const TestArbCell *c, guint n, NcmVector *res)
{
  gdouble err = 0.0, peak = 0.0;
  guint j;

  for (j = 0; j < n; j++)
  {
    err  = GSL_MAX (err, fabs (ncm_vector_get (res, c[j].ell) - c[j].value));
    peak = GSL_MAX (peak, fabs (c[j].value));
  }

  return err / peak;
}

/*
 * A reltol set after construction reaches the operators, which are built with the
 * multipole range and must be replaced. On exp (-x / 1000) over [100, 1e4], ell 0..7,
 * reltol 1e-4 leaves the result at about 1e-6 of the peak; after set_reltol (1e-11) it
 * is about 1e-12. cheb-reltol is 1e-12 throughout, so reltol is the setting that limits
 * the result.
 */
void
test_ncm_sbessel_integrator_levin_arb_set_reltol (void)
{
  guint n;
  const TestArbCell *c             = _test_arb_case ("set_reltol", &n);
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new_full (0, 7, 1.0e-4, 1.0e6, 21, 1200, 1.0e-4, 2, 1.0e-12);
  NcmVector *res                   = ncm_vector_new (8);
  TestArbEval ev                   = _test_arb_eval_new (c);
  gdouble loose, tight;

  g_assert_cmpuint (n, ==, 8);
  ncm_sbessel_integrator_integrate (NCM_SBESSEL_INTEGRATOR (sbilv), &_test_arb_integrand, c->a, c->b, 1.0, res, &ev);
  loose = _test_arb_block_error (c, n, res);

  ncm_sbessel_integrator_levin_set_reltol (sbilv, 1.0e-11);
  g_assert_cmpfloat (ncm_sbessel_integrator_levin_get_reltol (sbilv), ==, 1.0e-11);

  ncm_sbessel_integrator_integrate (NCM_SBESSEL_INTEGRATOR (sbilv), &_test_arb_integrand, c->a, c->b, 1.0, res, &ev);
  tight = _test_arb_block_error (c, n, res);

  g_assert_cmpfloat (tight, <, loose);
  g_assert_cmpfloat (tight, <, 1.0e-9);

  ncm_vector_free (res);
  ncm_sbessel_integrator_levin_free (sbilv);
}

/* Rebuilding the operators eight times in a row keeps the result right */
void
test_ncm_sbessel_integrator_levin_arb_repeated_set_reltol (void)
{
  guint n, i, j;
  const TestArbCell *c             = _test_arb_case ("repeated_set_reltol", &n);
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new (0, 5);
  NcmVector *res                   = ncm_vector_new (6);
  TestArbEval ev                   = _test_arb_eval_new (c);

  g_assert_cmpuint (n, ==, 6);

  for (i = 0; i < 8; i++)
  {
    ncm_sbessel_integrator_levin_set_reltol (sbilv, 1.0e-9 * (1.0 + 0.1 * i));
    ncm_sbessel_integrator_integrate (NCM_SBESSEL_INTEGRATOR (sbilv), &_test_arb_integrand, c->a, c->b, 1.0, res, &ev);

    for (j = 0; j < n; j++)
      _test_assert_close (ncm_vector_get (res, c[j].ell), c[j].value, 1.0e-7, 1.0e-14);
  }

  ncm_vector_free (res);
  ncm_sbessel_integrator_levin_free (sbilv);
}

void
test_ncm_sbessel_integrator_levin_arb_ell_zero (void)
{
  _test_arb_check_singles ("ell_zero", 1.0e-8, 1.0e-20);
}

/* ell 5 at k = 0.1, 1, 10 and 100 */
void
test_ncm_sbessel_integrator_levin_arb_k_scaling (void)
{
  _test_arb_check_singles ("k_scaling", 1.0e-8, 1.0e-20);
}

/* x on [5, 5.001] */
void
test_ncm_sbessel_integrator_levin_arb_narrow (void)
{
  _test_arb_check_singles ("narrow", 1.0e-8, 1.0e-15);
}

/* 1 on [1, 20] at ell = k = 100, 200 and 300 */
void
test_ncm_sbessel_integrator_levin_arb_high_ell (void)
{
  _test_arb_check_singles ("high_ell", 1.0e-8, 1.0e-15);
}

/* exp (-x / 100) on [1, 100] at ell 10, k 500 */
void
test_ncm_sbessel_integrator_levin_arb_oscillatory (void)
{
  _test_arb_check_singles ("oscillatory", 1.0e-8, 1.0e-15);
}

/* The order 2 integral at ell 30, k 0.5, where it is 8.7e-26: relative, not absolute, precision */
void
test_ncm_sbessel_integrator_levin_arb_evanescent (void)
{
  _test_arb_check_singles ("evanescent", 1.0e-9, 0.0);
}

/*
 * The 501 x 101 truth tables <name>_jl_500.bin of the data-file release (from
 * tests/tools/make_sbessel_jl_bin.py), table[l][i] = int_a^b f(x) j_l(kvals[i] x) dx
 * with f the Gaussian or rational shape of ncm_sbessel_integrator_integrate_gaussian()
 * and _rational(), computed in blocks of 8 multipoles. Every cell is held to
 * |I - I_table| <= 1e-8 times the largest |table| of its multipole, the default
 * cheb-reltol: relative error is not usable near the zeros of j_l.
 */
void
test_ncm_sbessel_integrator_levin_truth (gconstpointer data)
{
  const gchar *name         = (const gchar *) data;
  gchar *basename           = g_strdup_printf ("%s_jl_500.bin", name);
  gchar *filename           = nc_data_download_release_file (basename);
  NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmVarDict *vd            = ncm_serialize_var_dict_from_variant_file (ser, filename, TRUE);
  const gboolean rational   = (strcmp (name, "rational") == 0);
  const guint block         = 8;
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, block - 1));
  gdouble center, std, lb, ub;
  GArray *kvals, *lvals;
  NcmMatrix *table;
  NcmVector *peak;
  guint ell_min, n_ell, n_k, ell0;

  g_assert_true (ncm_var_dict_get_double (vd, "center", &center));
  g_assert_true (ncm_var_dict_get_double (vd, "std", &std));
  g_assert_true (ncm_var_dict_get_double (vd, "lower-bound", &lb));
  g_assert_true (ncm_var_dict_get_double (vd, "upper-bound", &ub));
  g_assert_true (ncm_var_dict_get_double_array (vd, "kvals", &kvals));
  g_assert_true (ncm_var_dict_get_int_array (vd, "lvals", &lvals));
  g_assert_true (ncm_var_dict_get_object (vd, "table", ser, (GObject **) &table));

  ell_min = g_array_index (lvals, gint, 0);
  n_ell   = lvals->len;
  n_k     = kvals->len;
  g_assert_cmpuint (ncm_matrix_nrows (table), ==, n_ell);
  g_assert_cmpuint (ncm_matrix_ncols (table), ==, n_k);

  peak = ncm_vector_new (n_ell);

  for (ell0 = 0; ell0 < n_ell; ell0++)
  {
    NcmVector *row = ncm_matrix_get_row (table, ell0);

    ncm_vector_set (peak, ell0, _test_absmax (row));
    g_assert_cmpfloat (ncm_vector_get (peak, ell0), >, 0.0);
    ncm_vector_free (row);
  }

  for (ell0 = 0; ell0 < n_ell; ell0 += block)
  {
    const guint ell1 = GSL_MIN (ell0 + block, n_ell) - 1;
    NcmVector *res   = ncm_vector_new (ell1 - ell0 + 1);
    guint i, j;

    ncm_sbessel_integrator_set_ell_range (sbi, ell_min + ell0, ell_min + ell1);

    for (i = 0; i < n_k; i++)
    {
      const gdouble k = g_array_index (kvals, gdouble, i);

      if (rational)
        ncm_sbessel_integrator_integrate_rational (sbi, center, std, lb, ub, k, res);
      else
        ncm_sbessel_integrator_integrate_gaussian (sbi, center, std, lb, ub, k, res);

      for (j = ell0; j <= ell1; j++)
      {
        g_assert_true (gsl_finite (ncm_vector_get (res, j - ell0)));
        g_assert_cmpfloat (fabs (ncm_vector_get (res, j - ell0) - ncm_matrix_get (table, j, i)), <=, 1.0e-8 * ncm_vector_get (peak, j));
      }
    }

    ncm_vector_free (res);
  }

  g_array_unref (kvals);
  g_array_unref (lvals);
  ncm_vector_free (peak);
  ncm_matrix_free (table);
  ncm_var_dict_unref (vd);
  ncm_serialize_free (ser);
  ncm_sbessel_integrator_free (sbi);
  g_free (filename);
  g_free (basename);
}

/* integrate over a block agrees with integrate_ell for each multipole */
void
test_ncm_sbessel_integrator_levin_ell_range_consistency (void)
{
  const guint ell_min       = 5;
  const guint ell_max       = 15;
  const gdouble rate        = 0.1;
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (ell_min, ell_max));
  NcmVector *res            = ncm_vector_new (ell_max - ell_min + 1);
  guint ell;

  ncm_sbessel_integrator_integrate (sbi, &_test_exp_decay, 1.0, 20.0, 1.0, res, (gpointer) & rate);

  for (ell = ell_min; ell <= ell_max; ell++)
  {
    const gdouble single = ncm_sbessel_integrator_integrate_ell (sbi, &_test_exp_decay, 1.0, 20.0, 1.0, ell, (gpointer) & rate);

    _test_assert_close (ncm_vector_get (res, ell - ell_min), single, 1.0e-10, 1.0e-14);
  }

  ncm_vector_free (res);
  ncm_sbessel_integrator_free (sbi);
}

/* A panel whose Bessel batch is identically zero (ell 500 on [1, 2]) never calls the kernel */
void
test_ncm_sbessel_integrator_levin_zero_bessel_batch (void)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (500, 500));
  NcmVector *res            = ncm_vector_new (1);
  guint calls               = 0;

  ncm_sbessel_integrator_integrate (sbi, &_test_counted_one, 1.0, 2.0, 1.0, res, &calls);

  g_assert_cmpuint (calls, ==, 0);
  g_assert_cmpfloat (ncm_vector_get (res, 0), ==, 0.0);

  ncm_vector_free (res);
  ncm_sbessel_integrator_free (sbi);
}

/* A new multipole range changes the result and agrees with integrate_ell on it */
void
test_ncm_sbessel_integrator_levin_ell_range_caching (void)
{
  const gdouble rate        = 0.1;
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 5));
  NcmVector *res1           = ncm_vector_new (6);
  NcmVector *res2           = ncm_vector_new (6);
  gboolean differ           = FALSE;
  guint i;

  ncm_sbessel_integrator_integrate (sbi, &_test_exp_decay, 1.0, 10.0, 1.0, res1, (gpointer) & rate);
  ncm_sbessel_integrator_set_ell_range (sbi, 10, 15);
  ncm_sbessel_integrator_integrate (sbi, &_test_exp_decay, 1.0, 10.0, 1.0, res2, (gpointer) & rate);

  /* numpy allclose defaults: rtol 1e-5, atol 1e-8 */
  for (i = 0; i < 6; i++)
    differ = differ || (fabs (ncm_vector_get (res1, i) - ncm_vector_get (res2, i)) > 1.0e-8 + 1.0e-5 * fabs (ncm_vector_get (res2, i)));

  g_assert_true (differ);

  for (i = 0; i < 6; i++)
  {
    const gdouble single = ncm_sbessel_integrator_integrate_ell (sbi, &_test_exp_decay, 1.0, 10.0, 1.0, 10 + i, (gpointer) & rate);

    _test_assert_close (ncm_vector_get (res2, i), single, 1.0e-10, 1.0e-14);
  }

  ncm_vector_free (res1);
  ncm_vector_free (res2);
  ncm_sbessel_integrator_free (sbi);
}

/* ell 10 to 25 in blocks of the given size (each divides 16) agree with one multipole at a time */
void
test_ncm_sbessel_integrator_levin_block_size (gconstpointer data)
{
  const guint block_size    = *(const guint *) data;
  const guint ell_min       = 10;
  const guint ell_max       = 25;
  const gdouble rate        = 0.2;
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (ell_min, ell_min + block_size - 1));
  NcmSBesselIntegrator *ref = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (ell_min, ell_min));
  NcmVector *res            = ncm_vector_new (block_size);
  guint ell0;

  for (ell0 = ell_min; ell0 <= ell_max; ell0 += block_size)
  {
    const guint ell1 = ell0 + block_size - 1;
    guint ell;

    ncm_sbessel_integrator_set_ell_range (sbi, ell0, ell1);
    ncm_sbessel_integrator_integrate (sbi, &_test_exp_decay, 1.0, 10.0, 1.0, res, (gpointer) & rate);

    for (ell = ell0; ell <= ell1; ell++)
    {
      const gdouble single = ncm_sbessel_integrator_integrate_ell (ref, &_test_exp_decay, 1.0, 10.0, 1.0, ell, (gpointer) & rate);

      _test_assert_close (ncm_vector_get (res, ell - ell0), single, 1.0e-10, 1.0e-14);
    }
  }

  ncm_vector_free (res);
  ncm_sbessel_integrator_free (sbi);
  ncm_sbessel_integrator_free (ref);
}

/*
 * The Gaussian (center 5, sigma 1) on [0.1, 10] at ell 2, k 500, with reltol and
 * cheb-reltol both set to @reltol and every other setting at its default. Each panel
 * spans many oscillations, where the adaptive QR's decay test can stop at too small
 * an order.
 */
static gdouble
_test_oscillatory_gaussian (gdouble reltol)
{
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new_full (2, 2,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_X_KNOTS_MIN,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_X_KNOTS_MAX,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_N_KNOTS,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_ELL_CACHE_MAX,
                                                                            reltol,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_MIN_ORDER,
                                                                            reltol);
  const gdouble val = ncm_sbessel_integrator_integrate_gaussian_ell (NCM_SBESSEL_INTEGRATOR (sbilv), 5.0, 1.0, 0.1, 10.0, 500.0, 2);

  ncm_sbessel_integrator_levin_free (sbilv);

  return val;
}

/* Every request from 1e-4 to 1e-10 is within 1 % of the 1e-12 result */
void
test_ncm_sbessel_integrator_levin_loose_tolerances_agree (void)
{
  const gdouble reltol_a[] = {1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10};
  const gdouble converged  = _test_oscillatory_gaussian (1.0e-12);
  guint i;

  for (i = 0; i < G_N_ELEMENTS (reltol_a); i++)
    _test_assert_close (_test_oscillatory_gaussian (reltol_a[i]), converged, 1.0e-2, 0.0);
}

/*
 * The error against the 1e-12 result grows by at most 3x from 1e-4 to 1e-6 and from 1e-6
 * to 1e-8. Below about 1e-7 the refinement lands on discrete levels that depend on the
 * BLAS kernel, so 1e-11 gets an absolute bound instead.
 */
void
test_ncm_sbessel_integrator_levin_accuracy_monotonic (void)
{
  const gdouble converged = _test_oscillatory_gaussian (1.0e-12);
  const gdouble err_4     = fabs (_test_oscillatory_gaussian (1.0e-4) / converged - 1.0);
  const gdouble err_6     = fabs (_test_oscillatory_gaussian (1.0e-6) / converged - 1.0);
  const gdouble err_8     = fabs (_test_oscillatory_gaussian (1.0e-8) / converged - 1.0);
  const gdouble err_11    = fabs (_test_oscillatory_gaussian (1.0e-11) / converged - 1.0);

  g_assert_cmpfloat (err_6, <=, 3.0 * err_4);
  g_assert_cmpfloat (err_8, <=, 3.0 * err_6);
  g_assert_cmpfloat (err_11, <, 1.0e-5);
}

/* A recording integrator at ell, reltol and cheb-reltol 1e-12, the rest at the defaults */
static NcmSBesselIntegratorLevin *
_test_recording_integrator (guint ell)
{
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new_full (ell, ell,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_X_KNOTS_MIN,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_X_KNOTS_MAX,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_N_KNOTS,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_ELL_CACHE_MAX,
                                                                            1.0e-12,
                                                                            NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_CHEB_MIN_ORDER,
                                                                            1.0e-12);

  ncm_sbessel_integrator_levin_set_record_panels (sbilv, TRUE);

  return sbilv;
}

void
test_ncm_sbessel_integrator_levin_record_off_by_default (void)
{
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new (2, 2);

  g_assert_false (ncm_sbessel_integrator_levin_get_record_panels (sbilv));
  g_assert_cmpuint (ncm_sbessel_integrator_levin_get_n_panel_records (sbilv), ==, 0);

  ncm_sbessel_integrator_levin_free (sbilv);
}

/* The recorded panel contributions add up to the result */
void
test_ncm_sbessel_integrator_levin_record_sum (void)
{
  NcmSBesselIntegratorLevin *sbilv = _test_recording_integrator (2);
  const gdouble val                = ncm_sbessel_integrator_integrate_gaussian_ell (NCM_SBESSEL_INTEGRATOR (sbilv), 5.0, 1.0, 0.1, 10.0, 50.0, 2);
  const guint n                    = ncm_sbessel_integrator_levin_get_n_panel_records (sbilv);
  gdouble total                    = 0.0;
  guint i;

  g_assert_cmpuint (n, >, 1);

  for (i = 0; i < n; i++)
    total += ncm_sbessel_integrator_levin_get_panel_contrib (sbilv, i);

  _test_assert_close (total, val, 1.0e-12, 0.0);

  ncm_sbessel_integrator_levin_free (sbilv);
}

/* The recorded panels abut, cover [k a, k b] in order and carry the multipole */
void
test_ncm_sbessel_integrator_levin_record_tile (void)
{
  const gdouble k                  = 50.0;
  const gdouble a                  = 0.1;
  const gdouble b                  = 10.0;
  NcmSBesselIntegratorLevin *sbilv = _test_recording_integrator (2);
  guint n, i;

  ncm_sbessel_integrator_integrate_gaussian_ell (NCM_SBESSEL_INTEGRATOR (sbilv), 5.0, 1.0, a, b, k, 2);
  n = ncm_sbessel_integrator_levin_get_n_panel_records (sbilv);
  g_assert_cmpuint (n, >, 0);

  _test_assert_close (ncm_sbessel_integrator_levin_get_panel_a (sbilv, 0), k * a, 1.0e-12, 0.0);
  _test_assert_close (ncm_sbessel_integrator_levin_get_panel_b (sbilv, n - 1), k * b, 1.0e-12, 0.0);

  for (i = 1; i < n; i++)
    _test_assert_close (ncm_sbessel_integrator_levin_get_panel_a (sbilv, i), ncm_sbessel_integrator_levin_get_panel_b (sbilv, i - 1), 1.0e-12, 0.0);

  for (i = 0; i < n; i++)
    g_assert_cmpint (ncm_sbessel_integrator_levin_get_panel_ell (sbilv, i), ==, 2);

  ncm_sbessel_integrator_levin_free (sbilv);
}

/* The records describe the most recent call only */
void
test_ncm_sbessel_integrator_levin_record_cleared (void)
{
  NcmSBesselIntegratorLevin *sbilv = _test_recording_integrator (2);
  guint first;

  ncm_sbessel_integrator_integrate_gaussian_ell (NCM_SBESSEL_INTEGRATOR (sbilv), 5.0, 1.0, 0.1, 10.0, 50.0, 2);
  first = ncm_sbessel_integrator_levin_get_n_panel_records (sbilv);
  ncm_sbessel_integrator_integrate_gaussian_ell (NCM_SBESSEL_INTEGRATOR (sbilv), 5.0, 1.0, 0.1, 10.0, 50.0, 2);

  g_assert_cmpuint (ncm_sbessel_integrator_levin_get_n_panel_records (sbilv), ==, first);

  ncm_sbessel_integrator_levin_free (sbilv);
}

void
test_ncm_sbessel_integrator_levin_record_disable (void)
{
  NcmSBesselIntegratorLevin *sbilv = _test_recording_integrator (2);

  ncm_sbessel_integrator_integrate_gaussian_ell (NCM_SBESSEL_INTEGRATOR (sbilv), 5.0, 1.0, 0.1, 10.0, 50.0, 2);
  g_assert_cmpuint (ncm_sbessel_integrator_levin_get_n_panel_records (sbilv), >, 0);

  ncm_sbessel_integrator_levin_set_record_panels (sbilv, FALSE);
  g_assert_cmpuint (ncm_sbessel_integrator_levin_get_n_panel_records (sbilv), ==, 0);

  ncm_sbessel_integrator_levin_free (sbilv);
}

/* j_l'(x), as scipy's spherical_jn(derivative=True) computes it */
static gdouble
_test_jl_deriv (guint ell, gdouble x)
{
  if (ell == 0)
    return -ncm_sf_sbessel (1, x);

  return ncm_sf_sbessel (ell - 1, x) - (ell + 1.0) / x * ncm_sf_sbessel (ell, x);
}

/*
 * For K = 1 the derivative integrals reduce to endpoint values:
 * int_a^b j_l'(k x) dx = [j_l(k x)]_a^b / k and int_a^b j_l''(k x) dx = [j_l'(k x)]_a^b / k.
 */
void
test_ncm_sbessel_integrator_levin_deriv_constant (void)
{
  const guint ell_a[] = {0, 1, 2, 5, 20, 60};
  const gdouble k_a[] = {0.5, 3.0, 40.0};
  const gdouble a     = 0.3;
  const gdouble b     = 12.0;
  NcmVector *res      = ncm_vector_new (1);
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (ell_a); i++)
  {
    NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (ell_a[i], ell_a[i]));

    for (j = 0; j < G_N_ELEMENTS (k_a); j++)
    {
      const gdouble k  = k_a[j];
      const gdouble d1 = (ncm_sf_sbessel (ell_a[i], k * b) - ncm_sf_sbessel (ell_a[i], k * a)) / k;
      const gdouble d2 = (_test_jl_deriv (ell_a[i], k * b) - _test_jl_deriv (ell_a[i], k * a)) / k;

      ncm_sbessel_integrator_integrate_deriv (sbi, &_test_one, a, b, k, 1, res, NULL);
      _test_assert_close (ncm_vector_get (res, 0), d1, 1.0e-11, 1.0e-16);

      ncm_sbessel_integrator_integrate_deriv (sbi, &_test_one, a, b, k, 2, res, NULL);
      _test_assert_close (ncm_vector_get (res, 0), d2, 1.0e-11, 1.0e-16);
    }

    ncm_sbessel_integrator_free (sbi);
  }

  ncm_vector_free (res);
  ncm_mpsf_sbessel_free_cache ();
}

/* An ell 2..33 block reproduces single-multipole results for both derivative orders */
void
test_ncm_sbessel_integrator_levin_deriv_batched (void)
{
  const guint ell_min       = 2;
  const guint ell_max       = 33;
  const guint ell_a[]       = {2, 9, 33};
  TestGauss g               = {4.0, 1.5, 1.0};
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (ell_min, ell_max));
  NcmVector *res_batch      = ncm_vector_new (ell_max - ell_min + 1);
  NcmVector *res_single     = ncm_vector_new (1);
  guint deriv, i;

  for (deriv = 1; deriv <= 2; deriv++)
  {
    ncm_sbessel_integrator_integrate_deriv (sbi, &_test_gauss, 0.5, 9.0, 7.0, deriv, res_batch, &g);

    for (i = 0; i < G_N_ELEMENTS (ell_a); i++)
    {
      NcmSBesselIntegrator *single = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (ell_a[i], ell_a[i]));

      ncm_sbessel_integrator_integrate_deriv (single, &_test_gauss, 0.5, 9.0, 7.0, deriv, res_single, &g);
      _test_assert_close (ncm_vector_get (res_batch, ell_a[i] - ell_min), ncm_vector_get (res_single, 0), 1.0e-10, 1.0e-16);

      ncm_sbessel_integrator_free (single);
    }
  }

  ncm_vector_free (res_batch);
  ncm_vector_free (res_single);
  ncm_sbessel_integrator_free (sbi);
}

/* Order zero dispatches to ncm_sbessel_integrator_integrate(), bit for bit */
void
test_ncm_sbessel_integrator_levin_deriv_zero (void)
{
  TestGauss g               = {4.0, 1.5, 1.0};
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 8));
  NcmVector *res_a          = ncm_vector_new (9);
  NcmVector *res_b          = ncm_vector_new (9);
  guint i;

  ncm_sbessel_integrator_integrate_deriv (sbi, &_test_gauss, 0.5, 9.0, 3.0, 0, res_a, &g);
  ncm_sbessel_integrator_integrate (sbi, &_test_gauss, 0.5, 9.0, 3.0, res_b, &g);

  for (i = 0; i < 9; i++)
    g_assert_cmpfloat (ncm_vector_get (res_a, i), ==, ncm_vector_get (res_b, i));

  ncm_vector_free (res_a);
  ncm_vector_free (res_b);
  ncm_sbessel_integrator_free (sbi);
}

static gdouble
_test_gauss_over_y2 (gpointer user_data, gdouble x, gdouble k)
{
  return _test_gauss (user_data, x, k) / gsl_pow_2 (k * x);
}

static gdouble
_test_gauss_over_y (gpointer user_data, gdouble x, gdouble k)
{
  return _test_gauss (user_data, x, k) / (k * x);
}

/*
 * j_l''(y) = (l (l - 1) / y^2 - 1) j_l(y) + (2 / y) j_{l+1}(y), so the order 2 integral
 * equals a combination of three order 0 integrals, a separate code path.
 */
void
test_ncm_sbessel_integrator_levin_deriv_recurrence (void)
{
  const guint ell           = 12;
  const gdouble a           = 0.8;
  const gdouble b           = 10.0;
  const gdouble k           = 5.0;
  TestGauss g               = {4.0, 1.5, 1.0};
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (ell, ell + 1));
  NcmVector *res            = ncm_vector_new (2);
  gdouble i_f, i_f_y2, i_f_y1;

  ncm_sbessel_integrator_integrate_deriv (sbi, &_test_gauss, a, b, k, 2, res, &g);

  i_f    = ncm_sbessel_integrator_integrate_ell (sbi, &_test_gauss, a, b, k, ell, &g);
  i_f_y2 = ncm_sbessel_integrator_integrate_ell (sbi, &_test_gauss_over_y2, a, b, k, ell, &g);
  i_f_y1 = ncm_sbessel_integrator_integrate_ell (sbi, &_test_gauss_over_y, a, b, k, ell + 1, &g);

  _test_assert_close (ncm_vector_get (res, 0), ell * (ell - 1.0) * i_f_y2 - i_f + 2.0 * i_f_y1, 1.0e-9, 0.0);

  ncm_vector_free (res);
  ncm_sbessel_integrator_free (sbi);
}

/*
 * The tau constraint is valid while the working order n satisfies n_F << n < N_min ~
 * 2 span / pi. These checks use one panel, x in [3162.28, 3672.79] at k = 1, ell
 * 200..207, and a Gaussian forcing of amplitude 1e-16 centred on it with width
 * span / width_fraction.
 */
#define TEST_TAU_A 3162.28
#define TEST_TAU_B 3672.79
#define TEST_TAU_ELL_MIN 200
#define TEST_TAU_ELL_MAX 207
#define TEST_TAU_N_ELL (TEST_TAU_ELL_MAX - TEST_TAU_ELL_MIN + 1)

static TestGauss
_test_tau_forcing (gdouble width_fraction)
{
  TestGauss g = {0.5 * (TEST_TAU_A + TEST_TAU_B), (TEST_TAU_B - TEST_TAU_A) / width_fraction, 1.0e-16};

  return g;
}

/* Integrates the tau test panel at the given tau-constraint-min-osc; returns the fallback count */
static guint
_test_tau_integrate (TestGauss *g, gdouble min_osc, NcmVector *res)
{
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new (TEST_TAU_ELL_MIN, TEST_TAU_ELL_MAX);
  guint fallbacks;

  ncm_sbessel_integrator_levin_set_tau_constraint_min_osc (sbilv, min_osc);
  ncm_sbessel_integrator_integrate (NCM_SBESSEL_INTEGRATOR (sbilv), &_test_gauss, TEST_TAU_A, TEST_TAU_B, 1.0, res, g);
  fallbacks = ncm_sbessel_integrator_levin_get_n_constraint_fallbacks (sbilv);

  ncm_sbessel_integrator_levin_free (sbilv);

  return fallbacks;
}

/* The test panel holds more oscillations than the default tau-constraint-min-osc */
void
test_ncm_sbessel_integrator_levin_tau_panel_qualifies (void)
{
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new (TEST_TAU_ELL_MIN, TEST_TAU_ELL_MAX);

  g_assert_cmpfloat (2.0 / M_PI * (TEST_TAU_B - TEST_TAU_A), >, ncm_sbessel_integrator_levin_get_tau_constraint_min_osc (sbilv));

  ncm_sbessel_integrator_levin_free (sbilv);
}

/* A wide forcing (width fraction 10) is solved with tau, no fallback, and agrees with Dirichlet */
void
test_ncm_sbessel_integrator_levin_tau_smooth (void)
{
  TestGauss g          = _test_tau_forcing (10.0);
  NcmVector *dirichlet = ncm_vector_new (TEST_TAU_N_ELL);
  NcmVector *tau       = ncm_vector_new (TEST_TAU_N_ELL);
  gboolean differ      = FALSE;
  gdouble peak;
  guint i;

  _test_tau_integrate (&g, 0.0, dirichlet);
  g_assert_cmpuint (_test_tau_integrate (&g, 200.0, tau), ==, 0);

  peak = _test_absmax (dirichlet);

  /* Identical bits would mean the tau path was never taken */
  for (i = 0; i < TEST_TAU_N_ELL; i++)
  {
    differ = differ || (ncm_vector_get (tau, i) != ncm_vector_get (dirichlet, i));
    _test_assert_close (ncm_vector_get (tau, i), ncm_vector_get (dirichlet, i), 0.0, 1.0e-8 * peak);
  }

  g_assert_true (differ);

  ncm_vector_free (dirichlet);
  ncm_vector_free (tau);
}

/*
 * A forcing whose order exceeds 0.75 N_min is redone with Dirichlet. The results are far
 * below int |F j_l| <= 1e-16 width sqrt(2 pi) / x0 (|j_l(x)| <= 1/x), so they are compared
 * on that scale; the fallback repeats the Dirichlet solve, so they agree to rounding.
 */
void
test_ncm_sbessel_integrator_levin_tau_falls_back (gconstpointer data)
{
  const gdouble width_fraction = *(const gdouble *) data;
  TestGauss g                  = _test_tau_forcing (width_fraction);
  const gdouble scale          = 1.0e-16 * g.sigma * sqrt (2.0 * M_PI) / g.center;
  NcmVector *dirichlet         = ncm_vector_new (TEST_TAU_N_ELL);
  NcmVector *tau               = ncm_vector_new (TEST_TAU_N_ELL);
  guint i;

  _test_tau_integrate (&g, 0.0, dirichlet);
  g_assert_cmpuint (_test_tau_integrate (&g, 200.0, tau), >, 0);

  for (i = 0; i < TEST_TAU_N_ELL; i++)
    _test_assert_close (ncm_vector_get (tau, i), ncm_vector_get (dirichlet, i), 0.0, 1.0e-14 * scale);

  ncm_vector_free (dirichlet);
  ncm_vector_free (tau);
}

/* With tau-constraint-min-osc 0 every panel is Dirichlet and no fallback is counted */
void
test_ncm_sbessel_integrator_levin_tau_rule_off (void)
{
  TestGauss g    = _test_tau_forcing (100.0);
  NcmVector *res = ncm_vector_new (TEST_TAU_N_ELL);

  g_assert_cmpuint (_test_tau_integrate (&g, 0.0, res), ==, 0);

  ncm_vector_free (res);
}

/*
 * The knot at turning-knot-margin times nu: the tau guard cannot certify a panel that
 * straddles the turning point, since its bound has min |x^2 - nu^2| in the denominator.
 * These checks integrate a Gaussian (center 3e4, sigma 6e3) on x in [1e2, 1e5], k = 1,
 * ell ell..ell + 7.
 */
typedef struct _TestTurningRun
{
  guint n_tau_solves;
  guint n_fallbacks;
  guint n_locked_eligible;
} TestTurningRun;

static TestTurningRun
_test_turning_run (guint ell, gdouble margin, gdouble min_osc, NcmVector *res)
{
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new (ell, ell + 7);
  TestGauss g                      = {3.0e4, 6.0e3, 1.0};
  TestTurningRun run;

  ncm_sbessel_integrator_levin_set_turning_knot_margin (sbilv, margin);
  ncm_sbessel_integrator_levin_set_tau_constraint_min_osc (sbilv, min_osc);

  ncm_sbessel_integrator_integrate (NCM_SBESSEL_INTEGRATOR (sbilv), &_test_gauss, 1.0e2, 1.0e5, 1.0, res, &g);

  run.n_tau_solves      = ncm_sbessel_integrator_levin_get_n_tau_solves (sbilv);
  run.n_fallbacks       = ncm_sbessel_integrator_levin_get_n_constraint_fallbacks (sbilv);
  run.n_locked_eligible = ncm_sbessel_integrator_levin_get_n_locked_eligible_solves (sbilv);

  ncm_sbessel_integrator_levin_free (sbilv);

  return run;
}

void
test_ncm_sbessel_integrator_levin_turning_default (void)
{
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new (2, 9);

  g_assert_cmpfloat (ncm_sbessel_integrator_levin_get_turning_knot_margin (sbilv), ==, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TURNING_KNOT_MARGIN);
  g_assert_cmpfloat (ncm_sbessel_integrator_levin_get_turning_knot_margin (sbilv), >, 1.0);

  ncm_sbessel_integrator_levin_free (sbilv);
}

/* The margin reads back, zero (off) included */
void
test_ncm_sbessel_integrator_levin_turning_round_trip (void)
{
  const gdouble margin_a[]         = {1.2, 0.0, 1.05};
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new (2, 9);
  guint i;

  for (i = 0; i < G_N_ELEMENTS (margin_a); i++)
  {
    ncm_sbessel_integrator_levin_set_turning_knot_margin (sbilv, margin_a[i]);
    g_assert_cmpfloat (ncm_sbessel_integrator_levin_get_turning_knot_margin (sbilv), ==, margin_a[i]);
  }

  ncm_sbessel_integrator_levin_free (sbilv);
}

/* Without the knot the straddling panel is refused; with it nothing is, and one more panel uses tau */
void
test_ncm_sbessel_integrator_levin_turning_removes_fallbacks (gconstpointer data)
{
  const guint ell        = *(const guint *) data;
  NcmVector *res         = ncm_vector_new (8);
  const TestTurningRun o = _test_turning_run (ell, 0.0, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_MIN_OSC, res);
  const TestTurningRun n = _test_turning_run (ell, 1.05, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_MIN_OSC, res);

  g_assert_cmpuint (o.n_fallbacks, >, 0);
  g_assert_cmpuint (n.n_fallbacks, ==, 0);
  g_assert_cmpuint (n.n_tau_solves, >, o.n_tau_solves);

  ncm_vector_free (res);
}

/* The knot changes which constraint runs, not the integral */
void
test_ncm_sbessel_integrator_levin_turning_same_integral (void)
{
  const guint ell_a[] = {200, 1000};
  NcmVector *off      = ncm_vector_new (8);
  NcmVector *on       = ncm_vector_new (8);
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (ell_a); i++)
  {
    gdouble peak;

    _test_turning_run (ell_a[i], 0.0, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_MIN_OSC, off);
    _test_turning_run (ell_a[i], 1.05, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_MIN_OSC, on);
    peak = _test_absmax (off);

    for (j = 0; j < 8; j++)
      _test_assert_close (ncm_vector_get (on, j), ncm_vector_get (off, j), 0.0, 1.0e-8 * peak);
  }

  ncm_vector_free (off);
  ncm_vector_free (on);
}

/* Below tau-constraint-min-osc the straddling panel is Dirichlet either way: no knot, same bits */
void
test_ncm_sbessel_integrator_levin_turning_low_ell (void)
{
  NcmVector *off = ncm_vector_new (8);
  NcmVector *on  = ncm_vector_new (8);
  guint j;

  _test_turning_run (20, 0.0, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_MIN_OSC, off);
  _test_turning_run (20, 1.05, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_MIN_OSC, on);

  for (j = 0; j < 8; j++)
    g_assert_cmpfloat (ncm_vector_get (on, j), ==, ncm_vector_get (off, j));

  ncm_vector_free (off);
  ncm_vector_free (on);
}

/* A guard fallback is counted once, as a fallback; with one block nothing is locked-eligible */
void
test_ncm_sbessel_integrator_levin_turning_not_locked (gconstpointer data)
{
  const guint ell          = *(const guint *) data;
  NcmVector *res           = ncm_vector_new (8);
  const TestTurningRun run = _test_turning_run (ell, 0.0, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TAU_CONSTRAINT_MIN_OSC, res);

  g_assert_cmpuint (run.n_fallbacks, >, 0);
  g_assert_cmpuint (run.n_locked_eligible, ==, 0);

  ncm_vector_free (res);
}

/* The knot is gated on the same oscillation count as the rule: rule off, no tau solve */
void
test_ncm_sbessel_integrator_levin_turning_rule_off (void)
{
  NcmVector *res           = ncm_vector_new (8);
  const TestTurningRun run = _test_turning_run (1000, NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TURNING_KNOT_MARGIN, 0.0, res);

  g_assert_cmpuint (run.n_tau_solves, ==, 0);

  ncm_vector_free (res);
}

/*
 * Integrators over one knot grid share their j_l rows, read-only. These checks integrate
 * a Gaussian (center 3e3, sigma 6e2) on x in [1e2, 1e4], k = 1, over 8 multipoles.
 */
static void
_test_shared_run (NcmSBesselIntegratorLevin *sbilv, NcmVector *res)
{
  TestGauss g = {3.0e3, 6.0e2, 1.0};

  ncm_sbessel_integrator_integrate (NCM_SBESSEL_INTEGRATOR (sbilv), &_test_gauss, 1.0e2, 1.0e4, 1.0, res, &g);
}

static void
_test_assert_equal_vectors (NcmVector *a, NcmVector *b)
{
  guint i;

  for (i = 0; i < ncm_vector_len (a); i++)
    g_assert_cmpfloat (ncm_vector_get (a, i), ==, ncm_vector_get (b, i));
}

/* A second integrator on the same grid, another block, leaves the first one's results unchanged */
void
test_ncm_sbessel_integrator_levin_shared_second (void)
{
  NcmSBesselIntegratorLevin *first  = ncm_sbessel_integrator_levin_new (200, 207);
  NcmSBesselIntegratorLevin *second = ncm_sbessel_integrator_levin_new (600, 607);
  NcmVector *before                 = ncm_vector_new (8);
  NcmVector *other                  = ncm_vector_new (8);
  NcmVector *after                  = ncm_vector_new (8);

  _test_shared_run (first, before);
  _test_shared_run (second, other);
  _test_shared_run (first, after);
  _test_assert_equal_vectors (after, before);

  ncm_vector_free (before);
  ncm_vector_free (other);
  ncm_vector_free (after);
  ncm_sbessel_integrator_levin_free (first);
  ncm_sbessel_integrator_levin_free (second);
}

void
test_ncm_sbessel_integrator_levin_shared_identical (void)
{
  NcmSBesselIntegratorLevin *a = ncm_sbessel_integrator_levin_new (200, 207);
  NcmSBesselIntegratorLevin *b = ncm_sbessel_integrator_levin_new (200, 207);
  NcmVector *res_a             = ncm_vector_new (8);
  NcmVector *res_b             = ncm_vector_new (8);

  _test_shared_run (a, res_a);
  _test_shared_run (b, res_b);
  _test_assert_equal_vectors (res_a, res_b);

  ncm_vector_free (res_a);
  ncm_vector_free (res_b);
  ncm_sbessel_integrator_levin_free (a);
  ncm_sbessel_integrator_levin_free (b);
}

/* Blocks insert their turning knot at different x, so that row is per instance */
void
test_ncm_sbessel_integrator_levin_shared_turning_row (void)
{
  NcmSBesselIntegratorLevin *low  = ncm_sbessel_integrator_levin_new (400, 407);
  NcmSBesselIntegratorLevin *high = ncm_sbessel_integrator_levin_new (900, 907);
  NcmVector *low_first            = ncm_vector_new (8);
  NcmVector *high_res             = ncm_vector_new (8);
  NcmVector *low_again            = ncm_vector_new (8);

  _test_shared_run (low, low_first);
  _test_shared_run (high, high_res);
  _test_shared_run (low, low_again);
  _test_assert_equal_vectors (low_again, low_first);

  ncm_vector_free (low_first);
  ncm_vector_free (high_res);
  ncm_vector_free (low_again);
  ncm_sbessel_integrator_levin_free (low);
  ncm_sbessel_integrator_levin_free (high);
}

void
test_ncm_sbessel_integrator_levin_traps (void)
{
  g_test_trap_subprocess ("/ncm/sbessel_integrator_levin/invalid/reltol_zero/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*assertion failed*");

  g_test_trap_subprocess ("/ncm/sbessel_integrator_levin/invalid/reltol_negative/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*assertion failed*");

  g_test_trap_subprocess ("/ncm/sbessel_integrator_levin/invalid/cheb_reltol_zero/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*assertion failed*");

  g_test_trap_subprocess ("/ncm/sbessel_integrator_levin/invalid/cheb_reltol_negative/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*assertion failed*");

  g_test_trap_subprocess ("/ncm/sbessel_integrator_levin/invalid/deriv/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not supported*");
}

void
test_ncm_sbessel_integrator_levin_invalid_reltol_zero (void)
{
  ncm_sbessel_integrator_levin_set_reltol (ncm_sbessel_integrator_levin_new (0, 3), 0.0);
}

void
test_ncm_sbessel_integrator_levin_invalid_reltol_negative (void)
{
  ncm_sbessel_integrator_levin_set_reltol (ncm_sbessel_integrator_levin_new (0, 3), -1.0e-8);
}

void
test_ncm_sbessel_integrator_levin_invalid_cheb_reltol_zero (void)
{
  ncm_sbessel_integrator_levin_set_cheb_reltol (ncm_sbessel_integrator_levin_new (0, 3), 0.0);
}

void
test_ncm_sbessel_integrator_levin_invalid_cheb_reltol_negative (void)
{
  ncm_sbessel_integrator_levin_set_cheb_reltol (ncm_sbessel_integrator_levin_new (0, 3), -1.0e-8);
}

/* Orders above 2 are not implemented anywhere */
void
test_ncm_sbessel_integrator_levin_invalid_deriv (void)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 3));
  NcmVector *res            = ncm_vector_new (4);

  ncm_sbessel_integrator_integrate_deriv (sbi, &_test_one, 1.0, 2.0, 1.0, 3, res, NULL);
}

