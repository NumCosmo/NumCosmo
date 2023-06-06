/***************************************************************************
 *            test_ncm_cfg.c
 *
 *  Mon Jun 05 12:04:44 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_cfg.c
 *
 * Copyright (C) 2023 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

typedef struct _TesNcmCfg
{
  guint place_holder;
} TesNcmCfg;

void test_ncm_cfg_new (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_free (TesNcmCfg *test, gconstpointer pdata);

void test_ncm_cfg_misc (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_logfile_set_logstream (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_logfile_on_off (TesNcmCfg *test, gconstpointer pdata);

void test_ncm_cfg_traps (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_invalid (TesNcmCfg *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/cfg/misc", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_misc,
              &test_ncm_cfg_free);

  g_test_add ("/ncm/cfg/logfile/set_logstream", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_logfile_set_logstream,
              &test_ncm_cfg_free);

  g_test_add ("/ncm/cfg/logfile/on_off", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_logfile_on_off,
              &test_ncm_cfg_free);

  g_test_add ("/ncm/cfg/traps", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_traps,
              &test_ncm_cfg_free);

#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/cfg/logfile/subprocess", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_invalid,
              &test_ncm_cfg_free);
#endif
  g_test_run ();
}

void
test_ncm_cfg_new (TesNcmCfg *test, gconstpointer pdata)
{
  test->place_holder = 0;
}

void
test_ncm_cfg_free (TesNcmCfg *test, gconstpointer pdata)
{
  /* NcmDiff *diff = test->diff; */
}

void
test_ncm_cfg_misc (TesNcmCfg *test, gconstpointer pdata)
{
  /* Testing MPI */
  {
    guint nslaves = ncm_cfg_mpi_nslaves ();

    g_assert_cmpint (nslaves, >=, 0);
  }

  /* Test get full path */
  {
    gchar *full_path = ncm_cfg_get_fullpath ("test_full_path_%d.txt", 1);

    g_assert_cmpstr (full_path, >=, ".numcosmo/test_full_path_1.txt");

    g_free (full_path);
  }

  /* Test string to comment */
  {
    gchar *comment = ncm_cfg_string_to_comment ("test string to comment, this is a very long comment that "
                                                "should be truncated to 80 characters, but it is not, so it "
                                                "will be truncated by the user.");

    g_assert_cmpstr (comment, >=, "###############################################################################");
    g_free (comment);
  }

  /* Setting n threads */
  {
    ncm_cfg_set_openmp_nthreads (1);
    ncm_cfg_set_openblas_nthreads (1);
    ncm_cfg_set_blis_nthreads (1);
    ncm_cfg_set_mkl_nthreads (1);
  }
}

void
test_ncm_cfg_logfile_set_logstream (TesNcmCfg *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    ncm_cfg_set_logstream (stderr);
    ncm_message ("This message should be printed in stderr %d", 1);

    return;
  }

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_stderr ("This message should be printed in stderr 1");
}

void
test_ncm_cfg_logfile_on_off (TesNcmCfg *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    ncm_cfg_set_logstream (stderr);

    ncm_cfg_logfile (FALSE);
    ncm_message ("This message should not be printed in stderr %d", 1);
    ncm_cfg_logfile (TRUE);
    ncm_message ("This message should be printed in stderr %d", 1);

    return;
  }

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_stderr ("This message should be printed in stderr 1");
}

void
test_ncm_cfg_traps (TesNcmCfg *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/cfg/logfile/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_cfg_invalid (TesNcmCfg *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

