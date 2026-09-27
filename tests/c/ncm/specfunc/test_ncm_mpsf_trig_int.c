/***************************************************************************
 *            test_ncm_mpsf_trig_int.c
 *
 *  Fri Nov 11 17:31:11 2022
 *  Copyright  2022  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2022 <vitenti@uel.br>
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

#include <gsl/gsl_sf_expint.h>

typedef struct _TestNcmMPSFSTrigInt
{
  guint ntests;
} TestNcmMPSFSTrigInt;

void test_ncm_mpsf_trig_int_new (TestNcmMPSFSTrigInt *test, gconstpointer pdata);
void test_ncm_mpsf_trig_int_free (TestNcmMPSFSTrigInt *test, gconstpointer pdata);

void test_ncm_mpsf_trig_int_sin_cmp_gsl (TestNcmMPSFSTrigInt *test, gconstpointer pdata);

void test_ncm_mpsf_trig_int_traps (TestNcmMPSFSTrigInt *test, gconstpointer pdata);
void test_ncm_mpsf_trig_int_invalid_st (TestNcmMPSFSTrigInt *test, gconstpointer pdata);
void test_ncm_mpsf_trig_int_odd (void);
void test_ncm_mpsf_trig_int_threads (void);

#define NTOT 100
#define XMAX 25.0
#define L 80

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/sf/trig_int/sin/cmp/gsl", TestNcmMPSFSTrigInt, NULL,
              &test_ncm_mpsf_trig_int_new,
              &test_ncm_mpsf_trig_int_sin_cmp_gsl,
              &test_ncm_mpsf_trig_int_free);

  g_test_add ("/ncm/sf/trig_int/traps", TestNcmMPSFSTrigInt, NULL,
              &test_ncm_mpsf_trig_int_new,
              &test_ncm_mpsf_trig_int_traps,
              &test_ncm_mpsf_trig_int_free);

  g_test_add ("/ncm/sf/trig_int/invalid/st/subprocess", TestNcmMPSFSTrigInt, NULL,
              &test_ncm_mpsf_trig_int_new,
              &test_ncm_mpsf_trig_int_invalid_st,
              &test_ncm_mpsf_trig_int_free);

  g_test_add_func ("/ncm/sf/trig_int/odd", &test_ncm_mpsf_trig_int_odd);
  g_test_add_func ("/ncm/sf/trig_int/threads", &test_ncm_mpsf_trig_int_threads);

  g_test_run ();

  ncm_mpsf_sin_int_free_cache ();
}

void
test_ncm_mpsf_trig_int_new (TestNcmMPSFSTrigInt *test, gconstpointer pdata)
{
}

void
test_ncm_mpsf_trig_int_free (TestNcmMPSFSTrigInt *test, gconstpointer pdata)
{
}

void
test_ncm_mpsf_trig_int_sin_cmp_gsl (TestNcmMPSFSTrigInt *test, gconstpointer pdata)
{
  guint i;

  for (i = 0; i < NTOT; i++)
  {
    const gdouble x      = 1.0 * pow (10.0, 0.0 + XMAX / (NTOT - 1.0) * i);
    const gdouble ncm_Si = ncm_sf_sin_int (x);
    const gdouble gsl_Si = gsl_sf_Si (x);

    ncm_assert_cmpdouble_e (ncm_Si, ==, gsl_Si, 1.0e-7, 0.0);
  }
}

void
test_ncm_mpsf_trig_int_traps (TestNcmMPSFSTrigInt *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/sf/trig_int/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_mpsf_trig_int_invalid_st (TestNcmMPSFSTrigInt *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

/* Si(-x) = -Si(x) exactly, on both branches: x = 1e3 and 5e4 take the asymptotic
 * series at 53 bits, the others the Taylor series. */
void
test_ncm_mpsf_trig_int_odd (void)
{
  const gdouble xs[] = {0.25, 3.0, 17.5, 1.0e3, 5.0e4};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (xs); i++)
  {
    g_assert_cmpfloat (ncm_sf_sin_int (-xs[i]), ==, -ncm_sf_sin_int (xs[i]));
    ncm_assert_cmpdouble_e (ncm_sf_sin_int (-xs[i]), ==, -gsl_sf_Si (xs[i]), 1.0e-7, 0.0);
  }
}

#define TRIG_NTHREADS 4
#define TRIG_NX 64

typedef struct _TrigWorker
{
  gdouble res[TRIG_NX];
} TrigWorker;

static gdouble
_trig_x (guint i)
{
  return 0.5 + 60.0 * i / (TRIG_NX - 1.0);
}

static gpointer
_trig_worker (gpointer data)
{
  TrigWorker *w = (TrigWorker *) data;
  guint i;

  for (i = 0; i < TRIG_NX; i++)
  {
    MPFR_DECL_INIT (res, 128);

    mpq_t q;

    mpq_init (q);
    mpq_set_d (q, _trig_x (i));
    ncm_mpsf_sin_int_mpfr (q, res, MPFR_RNDN);
    w->res[i] = mpfr_get_d (res, MPFR_RNDN);
    mpq_clear (q);
  }

  return NULL;
}

/* Concurrent evaluations give the serial results bit for bit. */
void
test_ncm_mpsf_trig_int_threads (void)
{
  TrigWorker serial;
  TrigWorker workers[TRIG_NTHREADS];
  GThread *threads[TRIG_NTHREADS];
  guint t, i;

  _trig_worker (&serial);

  for (t = 0; t < TRIG_NTHREADS; t++)
    threads[t] = g_thread_new ("trig-int", &_trig_worker, &workers[t]);

  for (t = 0; t < TRIG_NTHREADS; t++)
    g_thread_join (threads[t]);

  for (t = 0; t < TRIG_NTHREADS; t++)
    for (i = 0; i < TRIG_NX; i++)
      g_assert_cmpfloat (workers[t].res[i], ==, serial.res[i]);
}

