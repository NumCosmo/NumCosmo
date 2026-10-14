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

static const guint test_block_sizes[]  = {1, 4, 8, 16};
static const gdouble test_tau_widths[] = {30.0, 100.0};
static const guint test_turning_ells[] = {500, 1000};

gint
main (gint argc, gchar *argv[])
{
  guint i;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);

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

/* |actual - desired| <= atol + rtol |desired|, the numpy assert_allclose criterion */
static gdouble
_test_absmax (NcmVector *v)
{
  gdouble absmin, absmax;

  ncm_vector_get_absminmax (v, &absmin, &absmax);

  return absmax;
}

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

