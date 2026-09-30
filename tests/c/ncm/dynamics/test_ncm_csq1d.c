/***************************************************************************
 *            test_ncm_csq1d.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_csq1d.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_trig.h>
#include <glib.h>
#include <glib-object.h>

#include "test_ncm_csq1d_bessel.h"

static void
test_ncm_csq1d_properties (void)
{
  /* Construction defaults, then every setter against its getter and property. */
  TestCSQ1DBessel *b = test_csq1d_bessel_new (2.0, 1.0, TRUE);
  NcmCSQ1D *csq1d    = NCM_CSQ1D (b);
  gdouble d;
  gboolean save;
  NcmCSQ1DInitialStateType ict;

  g_assert_cmpfloat (ncm_csq1d_get_reltol (csq1d), ==, NCM_DEFAULT_PRECISION);
  g_assert_cmpfloat (ncm_csq1d_get_abstol (csq1d), ==, 0.0);
  g_assert_cmpfloat (ncm_csq1d_get_ti (csq1d), ==, 0.0);
  g_assert_cmpfloat (ncm_csq1d_get_tf (csq1d), ==, 1.0);
  g_assert_cmpfloat (ncm_csq1d_get_adiab_threshold (csq1d), ==, 1.0);
  g_assert_cmpfloat (ncm_csq1d_get_prop_threshold (csq1d), ==, 0.1);
  g_assert_true (ncm_csq1d_get_save_evol (csq1d));
  g_assert_cmpint (ncm_csq1d_get_initial_condition_type (csq1d), ==, NCM_CSQ1D_INITIAL_CONDITION_TYPE_AD_HOC);
  g_assert_cmpfloat (ncm_csq1d_get_vacuum_reltol (csq1d), ==, 1.0e-5);
  g_assert_cmpfloat (ncm_csq1d_get_vacuum_max_time (csq1d), ==, 1.0);

  ncm_csq1d_set_reltol (csq1d, 1.0e-6);
  ncm_csq1d_set_abstol (csq1d, 1.0e-7);
  ncm_csq1d_set_ti (csq1d, -100.0);
  ncm_csq1d_set_tf (csq1d, -1.0e-3);
  ncm_csq1d_set_adiab_threshold (csq1d, 1.0e-3);
  ncm_csq1d_set_prop_threshold (csq1d, 2.0e-3);
  ncm_csq1d_set_save_evol (csq1d, FALSE);
  ncm_csq1d_set_initial_condition_type (csq1d, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC2);
  ncm_csq1d_set_vacuum_reltol (csq1d, 1.0e-8);
  ncm_csq1d_set_vacuum_max_time (csq1d, -10.0);

  g_assert_cmpfloat (ncm_csq1d_get_reltol (csq1d), ==, 1.0e-6);
  g_assert_cmpfloat (ncm_csq1d_get_abstol (csq1d), ==, 1.0e-7);
  g_assert_cmpfloat (ncm_csq1d_get_ti (csq1d), ==, -100.0);
  g_assert_cmpfloat (ncm_csq1d_get_tf (csq1d), ==, -1.0e-3);
  g_assert_cmpfloat (ncm_csq1d_get_adiab_threshold (csq1d), ==, 1.0e-3);
  g_assert_cmpfloat (ncm_csq1d_get_prop_threshold (csq1d), ==, 2.0e-3);
  g_assert_false (ncm_csq1d_get_save_evol (csq1d));
  g_assert_cmpint (ncm_csq1d_get_initial_condition_type (csq1d), ==, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC2);
  g_assert_cmpfloat (ncm_csq1d_get_vacuum_reltol (csq1d), ==, 1.0e-8);
  g_assert_cmpfloat (ncm_csq1d_get_vacuum_max_time (csq1d), ==, -10.0);

  g_object_get (csq1d, "reltol", &d, NULL);
  g_assert_cmpfloat (d, ==, 1.0e-6);
  g_object_get (csq1d, "abstol", &d, NULL);
  g_assert_cmpfloat (d, ==, 1.0e-7);
  g_object_get (csq1d, "ti", &d, NULL);
  g_assert_cmpfloat (d, ==, -100.0);
  g_object_get (csq1d, "tf", &d, NULL);
  g_assert_cmpfloat (d, ==, -1.0e-3);
  g_object_get (csq1d, "adiab-threshold", &d, NULL);
  g_assert_cmpfloat (d, ==, 1.0e-3);
  g_object_get (csq1d, "prop-threshold", &d, NULL);
  g_assert_cmpfloat (d, ==, 2.0e-3);
  g_object_get (csq1d, "save-evol", &save, NULL);
  g_assert_false (save);
  g_object_get (csq1d, "vacuum-type", &ict, NULL);
  g_assert_cmpint (ict, ==, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC2);
  g_object_get (csq1d, "vacuum-reltol", &d, NULL);
  g_assert_cmpfloat (d, ==, 1.0e-8);
  g_object_get (csq1d, "vacuum-max-time", &d, NULL);
  g_assert_cmpfloat (d, ==, -10.0);

  g_object_set (csq1d, "vacuum-type", NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC4, "save-evol", TRUE, NULL);
  g_assert_cmpint (ncm_csq1d_get_initial_condition_type (csq1d), ==, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC4);
  g_assert_true (ncm_csq1d_get_save_evol (csq1d));

  ncm_csq1d_free (csq1d);
}

/* The Bessel system prepared at t in [-1e4, -1e-3] with the vacuum set before -10. */
static TestCSQ1DBessel *
_test_ncm_csq1d_bessel_prepared_new (NcmCSQ1DInitialStateType ict)
{
  TestCSQ1DBessel *b = test_csq1d_bessel_new (2.0, 1.0, TRUE);
  NcmCSQ1D *csq1d    = NCM_CSQ1D (b);

  ncm_csq1d_set_ti (csq1d, -1.0e4);
  ncm_csq1d_set_tf (csq1d, -1.0e-3);
  ncm_csq1d_set_reltol (csq1d, 1.0e-10);
  ncm_csq1d_set_abstol (csq1d, 0.0);
  ncm_csq1d_set_save_evol (csq1d, TRUE);
  ncm_csq1d_set_initial_condition_type (csq1d, ict);
  ncm_csq1d_set_vacuum_max_time (csq1d, -10.0);
  ncm_csq1d_set_vacuum_reltol (csq1d, 1.0e-8);

  return b;
}

static void
test_ncm_csq1d_prepare_aborts_subprocess (void)
{
  const gchar *which = g_getenv ("TEST_NCM_CSQ1D_ABORT");
  TestCSQ1DBessel *b;

  if (g_strcmp0 (which, "ad_hoc") == 0)
  {
    b = _test_ncm_csq1d_bessel_prepared_new (NCM_CSQ1D_INITIAL_CONDITION_TYPE_AD_HOC);
    ncm_csq1d_prepare (NCM_CSQ1D (b), NULL);
  }
  else if (g_strcmp0 (which, "nonadiab2") == 0)
  {
    b = _test_ncm_csq1d_bessel_prepared_new (NCM_CSQ1D_INITIAL_CONDITION_TYPE_NONADIABATIC2);
    ncm_csq1d_prepare (NCM_CSQ1D (b), NULL);
  }
  else if (g_strcmp0 (which, "ad_hoc_before") == 0)
  {
    NcmCSQ1DState *s0 = ncm_csq1d_state_new ();
    NcmCSQ1DState *s1 = ncm_csq1d_state_new ();

    b = _test_ncm_csq1d_bessel_prepared_new (NCM_CSQ1D_INITIAL_CONDITION_TYPE_AD_HOC);
    ncm_csq1d_compute_adiab (NCM_CSQ1D (b), NULL, -200.0, s0, NULL, NULL);
    ncm_csq1d_set_init_cond (NCM_CSQ1D (b), NULL, NCM_CSQ1D_EVOL_STATE_ADIABATIC, s0);
    ncm_csq1d_prepare (NCM_CSQ1D (b), NULL);
    ncm_csq1d_eval_at (NCM_CSQ1D (b), NULL, -300.0, s1);
  }
  else if (g_strcmp0 (which, "init_adiab") == 0)
  {
    /* At t = -2 the adiabatic alpha is -1.74, beyond the default threshold 1. */
    b = _test_ncm_csq1d_bessel_prepared_new (NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC4);
    ncm_csq1d_set_init_cond_adiab (NCM_CSQ1D (b), NULL, -2.0);
  }
}

static void
_test_ncm_csq1d_trap_abort (const gchar *which, const gchar *message)
{
  g_setenv ("TEST_NCM_CSQ1D_ABORT", which, TRUE);
  g_test_trap_subprocess ("/ncm/csq1d/prepare/aborts/subprocess", 0, G_TEST_SUBPROCESS_DEFAULT);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr (message);
  g_unsetenv ("TEST_NCM_CSQ1D_ABORT");
}

static void
test_ncm_csq1d_prepare_aborts (void)
{
  /* The documented refusals of ncm_csq1d_prepare() and ncm_csq1d_set_init_cond_adiab(). */
  _test_ncm_csq1d_trap_abort ("ad_hoc", "*initial conditions must be set*");
  _test_ncm_csq1d_trap_abort ("nonadiab2", "*not implemented*");
  _test_ncm_csq1d_trap_abort ("init_adiab", "*is not a valid adiabatic time*");
  _test_ncm_csq1d_trap_abort ("ad_hoc_before", "*is before the initial conditions*");
}

/* The phase of the exact positive-frequency mode, -arg H^(1)_a (-k t), for t < 0. */
static gdouble
_test_theta_exact (const gdouble a, const gdouble k, const gdouble t)
{
  const gdouble x = -k * t;

  return -atan2 (gsl_sf_bessel_Ynu (a, x), gsl_sf_bessel_Jnu (a, x));
}

static gdouble
_test_theta (NcmCSQ1D *csq1d, const gdouble t)
{
  return ncm_csq1d_eval_int_nu (csq1d, NULL, t) + ncm_csq1d_eval_delta_theta_at (csq1d, t);
}

/* Largest error of theta (t2) - theta (t1) against the exact difference, modulo 2 pi. */
static gdouble
_test_phase_pair_error (NcmCSQ1D *csq1d, const gdouble a, const gdouble k, const gdouble *t1, guint n1, const gdouble *t2, guint n2)
{
  gdouble err = 0.0;
  guint i, j;

  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < n2; j++)
    {
      const gdouble d = (_test_theta (csq1d, t2[j]) - _test_theta (csq1d, t1[i]))
                        - (_test_theta_exact (a, k, t2[j]) - _test_theta_exact (a, k, t1[i]));

      err = GSL_MAX (err, fabs (remainder (d, 2.0 * M_PI)));
    }
  }

  return err;
}

static void
test_ncm_csq1d_phase_hankel (void)
{
  /*
   * theta = int nu + delta theta, both with origin at the start of the numerical
   * evolution t_ode, against the exact Hankel phase: phase differences between times
   * after t_ode, on both sides of it, and before it, where the solution is the adiabatic
   * vacuum. Both integrals vanish at t_ode, int nu is k (t - t_ode) and theta is
   * continuous there.
   */
  const gdouble as[] = {2.0, 0.5};
  const gdouble k    = 1.0;
  const gdouble ti   = -1.0e4;
  guint l;

  for (l = 0; l < G_N_ELEMENTS (as); l++)
  {
    const gdouble a    = as[l];
    TestCSQ1DBessel *b = test_csq1d_bessel_new (a, k, TRUE);
    NcmCSQ1D *csq1d    = NCM_CSQ1D (b);
    gdouble t_ode, t_after[8], t_before[8], t_end[2], t_mid[1], e_nu = 0.0;
    guint i;

    ncm_csq1d_set_ti (csq1d, ti);
    ncm_csq1d_set_tf (csq1d, -1.0e-1);
    ncm_csq1d_set_reltol (csq1d, 1.0e-10);
    ncm_csq1d_set_abstol (csq1d, 0.0);
    ncm_csq1d_set_save_evol (csq1d, TRUE);
    ncm_csq1d_set_initial_condition_type (csq1d, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC4);
    ncm_csq1d_set_vacuum_max_time (csq1d, -10.0);
    ncm_csq1d_set_vacuum_reltol (csq1d, 1.0e-8);
    ncm_csq1d_prepare (csq1d, NULL);
    ncm_csq1d_prepare_phase_splines (csq1d, NULL);

    /* The vacuum time ncm_csq1d_prepare() finds, where the evolution starts. */
    g_assert_true (ncm_csq1d_find_adiab_time_limit (csq1d, NULL, ti, -10.0, 1.0e-8, &t_ode));

    {
      NcmCSQ1DState *s_evol  = ncm_csq1d_state_new ();
      NcmCSQ1DState *s_adiab = ncm_csq1d_state_new ();
      const gdouble t_after0 = t_ode * (1.0 - 1.0e-9);

      /* The evolution starts on the vacuum; measured 3.6e-15. */
      ncm_csq1d_eval_at (csq1d, NULL, t_after0, s_evol);
      ncm_csq1d_compute_adiab_frame (csq1d, NULL, NCM_CSQ1D_FRAME_ORIG, t_after0, s_adiab, NULL, NULL);
      g_assert_cmpfloat (ncm_csq1d_state_compute_distance (s_evol, s_adiab), <, 4.0e-14);

      ncm_csq1d_state_free (s_evol);
      ncm_csq1d_state_free (s_adiab);
    }

    for (i = 0; i < 8; i++)
    {
      t_after[i]  = -9.0 * pow (0.3 / 9.0, i / 7.0);
      t_before[i] = ti * 0.9 * pow (1.1 * t_ode / (0.9 * ti), i / 7.0);
    }

    t_end[0] = -5.0;
    t_end[1] = -0.2;
    t_mid[0] = 0.95 * t_ode;

    g_assert_cmpfloat (ncm_csq1d_eval_int_nu (csq1d, NULL, t_ode), ==, 0.0);
    g_assert_cmpfloat (ncm_csq1d_eval_delta_theta_at (csq1d, t_ode), ==, 0.0);

    for (i = 0; i < 8; i++)
    {
      e_nu = GSL_MAX (e_nu, fabs (ncm_csq1d_eval_int_nu (csq1d, NULL, t_before[i]) - k * (t_before[i] - t_ode)));
      e_nu = GSL_MAX (e_nu, fabs (ncm_csq1d_eval_int_nu (csq1d, NULL, t_after[i]) - k * (t_after[i] - t_ode)));
    }

    /* Measured, worst of a = 2 and 1/2, in radians: int nu 8.4e-10; pairs after 2.3e-8,
     * straddling 4.9e-8, before 3.6e-9. */
    g_assert_cmpfloat (e_nu, <, 1.0e-8);
    g_assert_cmpfloat (_test_phase_pair_error (csq1d, a, k, t_after, 7, &t_end[1], 1), <, 3.0e-7);
    g_assert_cmpfloat (_test_phase_pair_error (csq1d, a, k, t_before, 8, t_end, 2), <, 5.0e-7);
    g_assert_cmpfloat (_test_phase_pair_error (csq1d, a, k, t_before, 8, t_mid, 1), <, 4.0e-8);

    {
      const gdouble dt = 1.0e-9 * fabs (t_ode);

      /* Measured 1.7e-11. */
      g_assert_cmpfloat (fabs (_test_theta (csq1d, t_ode + dt) - _test_theta (csq1d, t_ode - dt) - 2.0 * k * dt), <, 2.0e-10);
    }

    ncm_csq1d_free (csq1d);
  }
}

/*
 * The exact state of the Bessel system at t < 0, from the positive-frequency mode
 * phi = C u^(-a) H^(1)_a (k u), u = -t, normalized by the Wronskian (|C|^2 = pi / 4),
 * with P_phi = C k u^(1 + a) H^(1)_(a + 1) (k u): J11 = 2 |phi|^2, J22 = 2 |P_phi|^2,
 * J12 = 2 Re (phi P_phi^*), and alpha = -asinh J12, gamma = ln (J22 / cosh alpha).
 */
static void
_test_bessel_exact_state (const gdouble a, const gdouble k, const gdouble t, NcmCSQ1DState *state)
{
  const gdouble u     = -t;
  const gdouble x     = k * u;
  const gdouble Ja    = gsl_sf_bessel_Jnu (a, x);
  const gdouble Ya    = gsl_sf_bessel_Ynu (a, x);
  const gdouble Ja1   = gsl_sf_bessel_Jnu (a + 1.0, x);
  const gdouble Ya1   = gsl_sf_bessel_Ynu (a + 1.0, x);
  const gdouble J12   = 0.5 * M_PI * x * (Ja * Ja1 + Ya * Ya1);
  const gdouble J22   = 0.5 * M_PI * k * k * pow (u, 2.0 + 2.0 * a) * (Ja1 * Ja1 + Ya1 * Ya1);
  const gdouble alpha = -asinh (J12);

  ncm_csq1d_state_set_ag (state, NCM_CSQ1D_FRAME_ORIG, t, alpha, log (J22) - gsl_sf_lncosh (alpha));
}

/*
 * Relative error of the complex structure of @s against @r: J11 and J22 relative to
 * themselves, J12 relative to sqrt (J11 J22) = sqrt (1 + J12^2). In strongly squeezed
 * states a small relative error is a large hyperbolic distance, so the distance is not
 * the accuracy measure there.
 */
static gdouble
_test_J_relerr (NcmCSQ1DState *s, NcmCSQ1DState *r)
{
  gdouble s11, s12, s22, r11, r12, r22;

  ncm_csq1d_state_get_J (s, &s11, &s12, &s22);
  ncm_csq1d_state_get_J (r, &r11, &r12, &r22);

  return GSL_MAX (GSL_MAX (fabs (s11 / r11 - 1.0), fabs (s22 / r22 - 1.0)), fabs (s12 - r12) / sqrt (r11 * r22));
}

/* Largest error of the evolved complex structure over the evolution times. */
static gdouble
_test_evolution_error (NcmCSQ1D *csq1d, const gdouble a, const gdouble k)
{
  NcmCSQ1DState *s_evol  = ncm_csq1d_state_new ();
  NcmCSQ1DState *s_exact = ncm_csq1d_state_new ();
  GArray *t_a            = ncm_csq1d_get_time_array (csq1d, NULL);
  gdouble err            = 0.0;
  guint i;

  for (i = 0; i < t_a->len; i++)
  {
    const gdouble t = g_array_index (t_a, gdouble, i);

    ncm_csq1d_eval_at (csq1d, NULL, t, s_evol);
    _test_bessel_exact_state (a, k, t, s_exact);
    err = GSL_MAX (err, _test_J_relerr (s_evol, s_exact));
  }

  g_array_unref (t_a);
  ncm_csq1d_state_free (s_evol);
  ncm_csq1d_state_free (s_exact);

  return err;
}

static void
test_ncm_csq1d_evolution_hankel (void)
{
  /*
   * From the fourth and second order adiabatic vacua the evolved state follows the exact
   * Bessel state down to t = -1e-3, deep in the non-adiabatic regime, within the accuracy
   * the vacuum was set at.
   */
  TestCSQ1DBessel *b4 = test_csq1d_bessel_new (2.0, 1.0, TRUE);
  TestCSQ1DBessel *b2 = test_csq1d_bessel_new (2.0, 1.0, TRUE);
  NcmCSQ1D *c4        = NCM_CSQ1D (b4);
  NcmCSQ1D *c2        = NCM_CSQ1D (b2);

  ncm_csq1d_set_ti (c4, -1.0e4);
  ncm_csq1d_set_tf (c4, -1.0e-3);
  ncm_csq1d_set_reltol (c4, 1.0e-10);
  ncm_csq1d_set_abstol (c4, 0.0);
  ncm_csq1d_set_save_evol (c4, TRUE);
  ncm_csq1d_set_initial_condition_type (c4, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC4);
  ncm_csq1d_set_vacuum_max_time (c4, -10.0);
  ncm_csq1d_set_vacuum_reltol (c4, 1.0e-8);
  ncm_csq1d_prepare (c4, NULL);

  ncm_csq1d_set_ti (c2, -1.0e5);
  ncm_csq1d_set_tf (c2, -1.0e-3);
  ncm_csq1d_set_reltol (c2, 1.0e-10);
  ncm_csq1d_set_abstol (c2, 0.0);
  ncm_csq1d_set_save_evol (c2, TRUE);
  ncm_csq1d_set_initial_condition_type (c2, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC2);
  ncm_csq1d_set_vacuum_max_time (c2, -1.0);
  ncm_csq1d_set_vacuum_reltol (c2, 1.0e-3);
  ncm_csq1d_prepare (c2, NULL);

  /* Measured 2.4e-8 and 2.5e-8. */
  g_assert_cmpfloat (_test_evolution_error (c4, 2.0, 1.0), <, 3.0e-7);
  g_assert_cmpfloat (_test_evolution_error (c2, 2.0, 1.0), <, 3.0e-7);

  ncm_csq1d_free (c4);
  ncm_csq1d_free (c2);
}

/* Error of the adiabatic vacuum at t against the exact state, and its estimate. */
static gdouble
_test_adiab_error (NcmCSQ1D *csq1d, const gdouble a, const gdouble k, const gdouble t, gdouble *estimate)
{
  NcmCSQ1DState *s_ad = ncm_csq1d_state_new ();
  NcmCSQ1DState *s_ex = ncm_csq1d_state_new ();
  gdouble ar, gr, err;

  ncm_csq1d_compute_adiab_frame (csq1d, NULL, NCM_CSQ1D_FRAME_ORIG, t, s_ad, &ar, &gr);
  _test_bessel_exact_state (a, k, t, s_ex);
  err = _test_J_relerr (s_ad, s_ex);

  if (estimate != NULL)
    estimate[0] = GSL_MAX (ar, gr);

  ncm_csq1d_state_free (s_ad);
  ncm_csq1d_state_free (s_ex);

  return err;
}

static void
test_ncm_csq1d_adiab_hankel (void)
{
  /*
   * Against the exact Bessel state the fourth order vacuum has an error falling as
   * |t|^-5, below its estimate, and at the time ncm_csq1d_find_adiab_time_limit()
   * places it the error is below the tolerance asked for. The second order has an
   * error falling as |t|^-3.
   */
  const gdouble a     = 2.0;
  const gdouble k     = 1.0;
  TestCSQ1DBessel *b4 = test_csq1d_bessel_new (a, k, TRUE);
  TestCSQ1DBessel *b2 = test_csq1d_bessel_new (a, k, TRUE);
  NcmCSQ1D *c4        = NCM_CSQ1D (b4);
  NcmCSQ1D *c2        = NCM_CSQ1D (b2);
  const gdouble ts[]  = {-3000.0, -300.0, -100.0, -30.0};

  /* Below 1e-12 the comparison reaches the accuracy of the GSL Bessel functions used as
   * reference, 4e-14 at x ~ 7000. */
  const gdouble precs[] = {1.0e-12, 1.0e-10, 1.0e-6};
  guint i;

  ncm_csq1d_set_initial_condition_type (c4, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC4);
  ncm_csq1d_set_initial_condition_type (c2, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC2);

  for (i = 0; i < G_N_ELEMENTS (ts); i++)
  {
    gdouble est;
    const gdouble err = _test_adiab_error (c4, a, k, ts[i], &est);

    g_assert_cmpfloat (err, <, est);
  }

  /* Measured slopes: 420 = 3.33^5 against 412 and 27 = 3^3 exactly. */
  ncm_assert_cmpdouble_e (_test_adiab_error (c4, a, k, -30.0, NULL) / _test_adiab_error (c4, a, k, -100.0, NULL), ==, pow (100.0 / 30.0, 5.0), 0.1, 0.0);
  ncm_assert_cmpdouble_e (_test_adiab_error (c2, a, k, -1000.0, NULL) / _test_adiab_error (c2, a, k, -3000.0, NULL), ==, 27.0, 0.01, 0.0);

  for (i = 0; i < G_N_ELEMENTS (precs); i++)
  {
    gdouble t_v;

    g_assert_true (ncm_csq1d_find_adiab_time_limit (c4, NULL, -1.0e4, -10.0, precs[i], &t_v));
    g_assert_cmpfloat (_test_adiab_error (c4, a, k, t_v, NULL), <=, precs[i]);
  }

  ncm_csq1d_free (c4);
  ncm_csq1d_free (c2);
}

static void
test_ncm_csq1d_adiab_finders (void)
{
  /*
   * find_adiab_time_limit returns t1 when both ends are adiabatic and FALSE when neither
   * is. For the Bessel system |F1| = 2.5 / |t| is smallest at the lower end, so
   * find_adiab_max returns it, with the lower border there and the upper one where
   * |F1 - F1_min| = epsilon; when F1 changes by less than epsilon the borders are the
   * ends.
   */
  TestCSQ1DBessel *b = test_csq1d_bessel_new (2.0, 1.0, TRUE);
  NcmCSQ1D *csq1d    = NCM_CSQ1D (b);
  gdouble t_v, t_min, F1_min, t_Bl, t_Bu;

  ncm_csq1d_set_initial_condition_type (csq1d, NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC4);

  g_assert_true (ncm_csq1d_find_adiab_time_limit (csq1d, NULL, -1.0e4, -100.0, 1.0e-2, &t_v));
  g_assert_cmpfloat (t_v, ==, -100.0);
  g_assert_false (ncm_csq1d_find_adiab_time_limit (csq1d, NULL, -100.0, -10.0, 1.0e-16, &t_v));

  t_min = ncm_csq1d_find_adiab_max (csq1d, NULL, -1.0e3, -10.0, 0.1, &F1_min, &t_Bl, &t_Bu);
  ncm_assert_cmpdouble_e (t_min, ==, -1.0e3, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (F1_min, ==, 2.5 / t_min, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (t_Bl, ==, -1.0e3, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (t_Bu, ==, 2.5 / (-0.1 + 2.5 / t_min), 1.0e-7, 0.0);

  ncm_csq1d_find_adiab_max (csq1d, NULL, -1.0e3, -500.0, 0.1, &F1_min, &t_Bl, &t_Bu);
  ncm_assert_cmpdouble_e (t_Bl, ==, -1.0e3, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (t_Bu, ==, -500.0, 1.0e-14, 0.0);

  ncm_csq1d_free (csq1d);
}

static void
test_ncm_csq1d_evolution_ad_hoc (void)
{
  /*
   * Ad hoc conditions equal to the adiabatic vacuum, at the vacuum time, give the same
   * solution as the adiabatic vacuum; eval_at reads the evolution from the time of the
   * conditions on.
   */
  TestCSQ1DBessel *ba = _test_ncm_csq1d_bessel_prepared_new (NCM_CSQ1D_INITIAL_CONDITION_TYPE_ADIABATIC4);
  TestCSQ1DBessel *bh = _test_ncm_csq1d_bessel_prepared_new (NCM_CSQ1D_INITIAL_CONDITION_TYPE_AD_HOC);
  NcmCSQ1D *ca        = NCM_CSQ1D (ba);
  NcmCSQ1D *ch        = NCM_CSQ1D (bh);
  NcmCSQ1DState *s0   = ncm_csq1d_state_new ();
  NcmCSQ1DState *sa   = ncm_csq1d_state_new ();
  NcmCSQ1DState *sh   = ncm_csq1d_state_new ();
  const gdouble ts[]  = {-100.0, -10.0, -1.0, -0.2, -1.0e-2};
  gdouble t_v;
  guint i;

  ncm_csq1d_prepare (ca, NULL);
  g_assert_true (ncm_csq1d_find_adiab_time_limit (ca, NULL, -1.0e4, -10.0, 1.0e-8, &t_v));

  ncm_csq1d_compute_adiab (ca, NULL, t_v, s0, NULL, NULL);
  ncm_csq1d_set_init_cond (ch, NULL, NCM_CSQ1D_EVOL_STATE_ADIABATIC, s0);
  ncm_csq1d_prepare (ch, NULL);

  for (i = 0; i < G_N_ELEMENTS (ts); i++)
  {
    ncm_csq1d_eval_at (ca, NULL, ts[i], sa);
    ncm_csq1d_eval_at (ch, NULL, ts[i], sh);
    g_assert_cmpfloat (_test_J_relerr (sh, sa), <, 1.0e-14);
  }

  ncm_csq1d_state_free (s0);
  ncm_csq1d_state_free (sa);
  ncm_csq1d_state_free (sh);
  ncm_csq1d_free (ca);
  ncm_csq1d_free (ch);
}

static void
test_ncm_csq1d_defaults (void)
{
  /*
   * A subclass implementing only xi, nu and F1 gets nu^2 = nu nu, m = e^xi / nu and
   * F2 = F1' / (2 nu) from NcmCSQ1D; the Bessel system has them in closed form.
   */
  const gdouble a          = 2.0;
  const gdouble k          = 1.3;
  TestCSQ1DBesselMin *bmin = test_csq1d_bessel_min_new (a, k, TRUE);
  NcmCSQ1D *csq1d          = NCM_CSQ1D (bmin);
  gdouble max_F2_err       = 0.0;
  guint i;

  for (i = 0; i < 50; i++)
  {
    const gdouble t = -pow (10.0, -1.0 + 5.0 * i / 49.0);

    ncm_assert_cmpdouble_e (ncm_csq1d_eval_nu2 (csq1d, NULL, t), ==, k * k, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_csq1d_eval_m (csq1d, NULL, t), ==, test_csq1d_bessel_m (a, k, -1.0, t), 1.0e-13, 0.0);

    max_F2_err = GSL_MAX (max_F2_err, fabs (ncm_csq1d_eval_F2 (csq1d, NULL, t) / test_csq1d_bessel_F2 (a, k, -1.0, t) - 1.0));
  }

  /* Measured 3.0e-13 over t in [-1e4, -0.1]. */
  g_assert_cmpfloat (max_F2_err, <, 3.0e-12);

  ncm_csq1d_free (csq1d);
}

/* Points (alpha, gamma) covering the branches of the distance: small, large alpha, large
 * gamma differences. */
static const gdouble _test_ag[][2] =
{
  {0.0, 0.0}, {0.3, -0.2}, {-0.7, 0.4}, {1.5, 0.1}, {-2.5, -1.3}, {0.05, 3.0}, {4.0, -4.0},
};

/* The hyperbolic distance from the upper half-plane coordinates. */
static gdouble
_test_dist_half_plane (NcmCSQ1DState *s0, NcmCSQ1DState *s1)
{
  gdouble x0, lny0, x1, lny1;

  ncm_csq1d_state_get_poincare_half_plane (s0, &x0, &lny0);
  ncm_csq1d_state_get_poincare_half_plane (s1, &x1, &lny1);

  {
    const gdouble y0 = exp (lny0);
    const gdouble y1 = exp (lny1);

    return acosh (1.0 + (gsl_pow_2 (x0 - x1) + gsl_pow_2 (y0 - y1)) / (2.0 * y0 * y1));
  }
}

/* The hyperbolic distance from the Poincare disc coordinates. */
static gdouble
_test_dist_disc (NcmCSQ1DState *s0, NcmCSQ1DState *s1)
{
  gdouble x0, y0, x1, y1;

  ncm_csq1d_state_get_poincare_disc (s0, &x0, &y0);
  ncm_csq1d_state_get_poincare_disc (s1, &x1, &y1);

  return acosh (1.0 + 2.0 * (gsl_pow_2 (x0 - x1) + gsl_pow_2 (y0 - y1)) / ((1.0 - x0 * x0 - y0 * y0) * (1.0 - x1 * x1 - y1 * y1)));
}

static void
test_ncm_csq1d_state_maps (void)
{
  /*
   * The (chi, U+-) maps invert each other, with chi = sinh alpha and
   * U+- = ln cosh alpha +- gamma. The complex structure J has unit determinant,
   * e^U- = J11 and e^U+ = J22. The vector (phi, P_phi) rebuilds J as
   * J11 = 2 |phi|^2, J22 = 2 |P_phi|^2, J12 = 2 Re (phi P_phi^*), with Wronskian
   * phi P_phi^* - phi^* P_phi = i, in the phase where phi is real and positive.
   */
  NcmCSQ1DState *s = ncm_csq1d_state_new ();
  NcmCSQ1DState *r = ncm_csq1d_state_new ();
  guint i;

  for (i = 0; i < G_N_ELEMENTS (_test_ag); i++)
  {
    const gdouble alpha = _test_ag[i][0];
    const gdouble gamma = _test_ag[i][1];
    gdouble chi, Up, Um, a, g, J11, J12, J22, phi[2], Pphi[2];

    ncm_csq1d_state_set_ag (s, NCM_CSQ1D_FRAME_ORIG, 1.0, alpha, gamma);

    ncm_csq1d_state_get_up (s, &chi, &Up);
    ncm_assert_cmpdouble_e (chi, ==, sinh (alpha), 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (Up, ==, log (cosh (alpha)) + gamma, 1.0e-14, 1.0e-15);
    ncm_csq1d_state_set_up (r, NCM_CSQ1D_FRAME_ORIG, 1.0, chi, Up);
    ncm_csq1d_state_get_ag (r, &a, &g);
    ncm_assert_cmpdouble_e (a, ==, alpha, 1.0e-14, 1.0e-15);
    ncm_assert_cmpdouble_e (g, ==, gamma, 1.0e-14, 1.0e-14);

    ncm_csq1d_state_get_um (s, &chi, &Um);
    ncm_assert_cmpdouble_e (Um, ==, log (cosh (alpha)) - gamma, 1.0e-14, 1.0e-15);
    ncm_csq1d_state_set_um (r, NCM_CSQ1D_FRAME_ORIG, 1.0, chi, Um);
    ncm_csq1d_state_get_ag (r, &a, &g);
    ncm_assert_cmpdouble_e (a, ==, alpha, 1.0e-14, 1.0e-15);
    ncm_assert_cmpdouble_e (g, ==, gamma, 1.0e-14, 1.0e-14);

    ncm_csq1d_state_get_J (s, &J11, &J12, &J22);
    ncm_assert_cmpdouble_e (J11 * J22 - J12 * J12, ==, 1.0, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (J11, ==, exp (Um), 1.0e-13, 0.0);
    ncm_assert_cmpdouble_e (J22, ==, exp (Up), 1.0e-13, 0.0);

    ncm_csq1d_state_get_phi_Pphi (s, phi, Pphi);
    ncm_assert_cmpdouble_e (2.0 * (phi[0] * phi[0] + phi[1] * phi[1]), ==, J11, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (2.0 * (Pphi[0] * Pphi[0] + Pphi[1] * Pphi[1]), ==, J22, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (2.0 * (phi[0] * Pphi[0] + phi[1] * Pphi[1]), ==, J12, 1.0e-14, 1.0e-15);
    ncm_assert_cmpdouble_e (2.0 * (phi[1] * Pphi[0] - phi[0] * Pphi[1]), ==, 1.0, 1.0e-14, 0.0);

    /* phi is real and positive. */
    g_assert_cmpfloat (phi[1], ==, 0.0);
    g_assert_cmpfloat (phi[0], >, 0.0);

    /* Half-plane x = chi e^-U+, ln y = -U+; disc x = chi / (1 + (e^U+ + e^U-) / 2),
     * y = -(e^U+ - e^U-) / 2 / (1 + (e^U+ + e^U-) / 2). */
    {
      gdouble x, lny, y;

      ncm_csq1d_state_get_poincare_half_plane (s, &x, &lny);
      ncm_assert_cmpdouble_e (x, ==, chi * exp (-Up), 1.0e-14, 1.0e-15);
      ncm_assert_cmpdouble_e (lny, ==, -Up, 1.0e-14, 1.0e-15);

      ncm_csq1d_state_get_poincare_disc (s, &x, &y);
      ncm_assert_cmpdouble_e (x, ==, chi / (1.0 + 0.5 * (exp (Up) + exp (Um))), 1.0e-14, 1.0e-15);
      ncm_assert_cmpdouble_e (y, ==, -0.5 * (exp (Up) - exp (Um)) / (1.0 + 0.5 * (exp (Up) + exp (Um))), 1.0e-14, 1.0e-15);
    }
  }

  ncm_csq1d_state_free (s);
  ncm_csq1d_state_free (r);
}

static void
test_ncm_csq1d_state_distance (void)
{
  /* The distance agrees with those of the half-plane and disc coordinates, is
   * symmetric, and vanishes between equal points. */
  NcmCSQ1DState *s0 = ncm_csq1d_state_new ();
  NcmCSQ1DState *s1 = ncm_csq1d_state_new ();
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (_test_ag); i++)
  {
    ncm_csq1d_state_set_ag (s0, NCM_CSQ1D_FRAME_ORIG, 1.0, _test_ag[i][0], _test_ag[i][1]);
    g_assert_cmpfloat (ncm_csq1d_state_compute_distance (s0, s0), ==, 0.0);

    for (j = 0; j < G_N_ELEMENTS (_test_ag); j++)
    {
      gdouble d;

      if (i == j)
        continue;

      ncm_csq1d_state_set_ag (s1, NCM_CSQ1D_FRAME_ORIG, 1.0, _test_ag[j][0], _test_ag[j][1]);
      d = ncm_csq1d_state_compute_distance (s0, s1);

      ncm_assert_cmpdouble_e (ncm_csq1d_state_compute_distance (s1, s0), ==, d, 1.0e-14, 0.0);
      /* Measured 6.7e-16 and 6.9e-15. */
      ncm_assert_cmpdouble_e (_test_dist_half_plane (s0, s1), ==, d, 1.0e-14, 0.0);
      ncm_assert_cmpdouble_e (_test_dist_disc (s0, s1), ==, d, 1.0e-13, 0.0);
    }
  }

  ncm_csq1d_state_free (s0);
  ncm_csq1d_state_free (s1);
}

static void
test_ncm_csq1d_state_distance_exact (void)
{
  /*
   * Along alpha at fixed gamma the distance is |delta alpha|; along gamma at fixed
   * alpha it is 2 asinh (cosh (alpha) sinh (delta gamma / 2)). Both down to tiny
   * separations and at large alpha, where cancellations lose them.
   */
  const gdouble alphas[] = {0.0, 0.5, 1.5, 5.0, 50.0, 300.0};
  const gdouble deltas[] = {1.0e-12, 1.0e-8, 1.0e-4, 0.1, 3.0};
  NcmCSQ1DState *s0      = ncm_csq1d_state_new ();
  NcmCSQ1DState *s1      = ncm_csq1d_state_new ();
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (alphas); i++)
  {
    const gdouble alpha = alphas[i];
    const gdouble gamma = -0.4;

    ncm_csq1d_state_set_ag (s0, NCM_CSQ1D_FRAME_ORIG, 1.0, alpha, gamma);
    g_assert_cmpfloat (ncm_csq1d_state_compute_distance (s0, s0), ==, 0.0);

    for (j = 0; j < G_N_ELEMENTS (deltas); j++)
    {
      const gdouble d = deltas[j];

      ncm_csq1d_state_set_ag (s1, NCM_CSQ1D_FRAME_ORIG, 1.0, alpha + d, gamma);
      ncm_assert_cmpdouble_e (ncm_csq1d_state_compute_distance (s0, s1), ==, fabs ((alpha + d) - alpha), 1.0e-12, 0.0);

      ncm_csq1d_state_set_ag (s1, NCM_CSQ1D_FRAME_ORIG, 1.0, alpha, gamma + d);
      ncm_assert_cmpdouble_e (ncm_csq1d_state_compute_distance (s0, s1), ==,
                              2.0 * asinh (cosh (alpha) * sinh (0.5 * ((gamma + d) - gamma))), 1.0e-12, 0.0);
    }
  }

  ncm_csq1d_state_free (s0);
  ncm_csq1d_state_free (s1);
}

static void
test_ncm_csq1d_state_circle (void)
{
  /* The circle of radius r around a point lies at distance r from it, for any angle;
   * the angle 0 moves alpha by r at fixed gamma. */
  const gdouble radii[] = {1.0e-4, 1.0e-3, 0.3, 2.0, 30.0, 1.0e4};
  NcmCSQ1DState *s      = ncm_csq1d_state_new ();
  NcmCSQ1DState *c      = ncm_csq1d_state_new ();
  gdouble max_err       = 0.0;
  guint i, j, l;

  for (i = 0; i < G_N_ELEMENTS (_test_ag); i++)
  {
    ncm_csq1d_state_set_ag (s, NCM_CSQ1D_FRAME_ORIG, 1.0, _test_ag[i][0], _test_ag[i][1]);

    for (j = 0; j < G_N_ELEMENTS (radii); j++)
    {
      gdouble a, g;

      ncm_csq1d_state_get_circle (s, radii[j], 0.0, c);
      ncm_csq1d_state_get_ag (c, &a, &g);
      ncm_assert_cmpdouble_e (a, ==, _test_ag[i][0] + radii[j], 1.0e-13, 1.0e-15);
      ncm_assert_cmpdouble_e (g, ==, _test_ag[i][1], 1.0e-14, 1.0e-15);

      for (l = 0; l < 25; l++)
      {
        ncm_csq1d_state_get_circle (s, radii[j], -4.0 * M_PI + 8.0 * M_PI * l / 24.0 + 0.1, c);
        max_err = GSL_MAX (max_err, fabs (ncm_csq1d_state_compute_distance (s, c) / radii[j] - 1.0));
      }
    }
  }

  /* Measured 5.5e-11, at the smallest radius, over r in [1e-4, 1e4]. */
  g_assert_cmpfloat (max_err, <, 5.0e-10);

  ncm_csq1d_state_free (s);
  ncm_csq1d_state_free (c);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/csq1d/defaults", &test_ncm_csq1d_defaults);
  g_test_add_func ("/ncm/csq1d/properties", &test_ncm_csq1d_properties);
  g_test_add_func ("/ncm/csq1d/phase/hankel", &test_ncm_csq1d_phase_hankel);
  g_test_add_func ("/ncm/csq1d/evolution/hankel", &test_ncm_csq1d_evolution_hankel);
  g_test_add_func ("/ncm/csq1d/evolution/ad_hoc", &test_ncm_csq1d_evolution_ad_hoc);
  g_test_add_func ("/ncm/csq1d/adiab/hankel", &test_ncm_csq1d_adiab_hankel);
  g_test_add_func ("/ncm/csq1d/adiab/finders", &test_ncm_csq1d_adiab_finders);
  g_test_add_func ("/ncm/csq1d/prepare/aborts", &test_ncm_csq1d_prepare_aborts);
  g_test_add_func ("/ncm/csq1d/prepare/aborts/subprocess", &test_ncm_csq1d_prepare_aborts_subprocess);
  g_test_add_func ("/ncm/csq1d/state/maps", &test_ncm_csq1d_state_maps);
  g_test_add_func ("/ncm/csq1d/state/distance", &test_ncm_csq1d_state_distance);
  g_test_add_func ("/ncm/csq1d/state/distance/exact", &test_ncm_csq1d_state_distance_exact);
  g_test_add_func ("/ncm/csq1d/state/circle", &test_ncm_csq1d_state_circle);

  g_test_run ();
}

