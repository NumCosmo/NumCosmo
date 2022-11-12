/***************************************************************************
 *            test_ncm_mpsf_0F1.c
 *
 *  Fri Nov 11 12:11:41 2022
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

#include <gsl/gsl_sf_hyperg.h>

typedef struct _TestNcmMPSFS0F1
{
  guint ntests;
} TestNcmMPSFS0F1;

void test_ncm_mpsf_0F1_new (TestNcmMPSFS0F1 *test, gconstpointer pdata);
void test_ncm_mpsf_0F1_free (TestNcmMPSFS0F1 *test, gconstpointer pdata);

void test_ncm_mpsf_0F1_cmp_gsl (TestNcmMPSFS0F1 *test, gconstpointer pdata);

void test_ncm_mpsf_0F1_traps (TestNcmMPSFS0F1 *test, gconstpointer pdata);
void test_ncm_mpsf_0F1_invalid_st (TestNcmMPSFS0F1 *test, gconstpointer pdata);

#define NTOT 100
#define XMAX 2.0
#define L 80

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_set_nonfatal_assertions ();
  
  g_test_add ("/ncm/sf/0F1/cmp/gsl", TestNcmMPSFS0F1, NULL,
              &test_ncm_mpsf_0F1_new,
              &test_ncm_mpsf_0F1_cmp_gsl,
              &test_ncm_mpsf_0F1_free);
  
  g_test_add ("/ncm/sf/0F1/traps", TestNcmMPSFS0F1, NULL,
              &test_ncm_mpsf_0F1_new,
              &test_ncm_mpsf_0F1_traps,
              &test_ncm_mpsf_0F1_free);
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/sf/0F1/invalid/st/subprocess", TestNcmMPSFS0F1, NULL,
              &test_ncm_mpsf_0F1_new,
              &test_ncm_mpsf_0F1_invalid_st,
              &test_ncm_mpsf_0F1_free);
#endif
  g_test_run ();
}

void
test_ncm_mpsf_0F1_new (TestNcmMPSFS0F1 *test, gconstpointer pdata)
{
}

void
test_ncm_mpsf_0F1_free (TestNcmMPSFS0F1 *test, gconstpointer pdata)
{
  ncm_mpsf_0F1_free_cache ();
}

void
test_ncm_mpsf_0F1_cmp_gsl (TestNcmMPSFS0F1 *test, gconstpointer pdata)
{
  guint i, j;
  
  for (j = 0; j <= L; j++)
  {
    const gdouble nu = g_test_rand_double_range (0.0, 1.0 * L);
    for (i = 0; i < NTOT; i++)
    {
      const gdouble x       = j * 1.0 * pow (10.0, 0.0 + XMAX / (NTOT - 1.0) * i);
      const gdouble ncm_0F1 = ncm_sf_0F1 (nu, x);
      const gdouble gsl_0F1 = gsl_sf_hyperg_0F1 (nu, x);
      
      ncm_assert_cmpdouble_e (ncm_0F1, ==, gsl_0F1, 1.0e-7, 0.0);
    }
  }
}

void
test_ncm_mpsf_0F1_traps (TestNcmMPSFS0F1 *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/sf/0F1/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_mpsf_0F1_invalid_st (TestNcmMPSFS0F1 *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

