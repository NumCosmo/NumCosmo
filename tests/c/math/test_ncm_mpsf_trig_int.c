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

  g_test_run ();
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

