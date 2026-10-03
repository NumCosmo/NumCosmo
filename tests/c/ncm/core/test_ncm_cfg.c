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
#include <fftw3.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <glib/gstdio.h>

typedef struct _TesNcmCfg
{
  guint place_holder;
} TesNcmCfg;

static gchar *test_ncm_cfg_executable;

static void
_test_ncm_cfg_spawn (const gchar *test_path, gchar **envp, gboolean should_pass, const gchar *expected_output)
{
  gchar *argv[] = {test_ncm_cfg_executable, "-p", (gchar *) test_path, NULL};
  gchar *output = NULL;
  gchar *errors = NULL;
  gint status;
  GError *error = NULL;

  g_assert_true (g_spawn_sync (NULL, argv, envp, 0, NULL, NULL, &output, &errors, &status, &error));
  g_assert_no_error (error);

  if (should_pass)
  {
    g_assert_true (g_spawn_check_wait_status (status, &error));
    g_assert_no_error (error);
    g_assert_null (strstr (output, "GLib-DEBUG: setenv()/putenv()"));
    g_assert_null (strstr (output, "GLib-DEBUG: unsetenv()"));
    g_assert_null (strstr (errors, "GLib-DEBUG: setenv()/putenv()"));
    g_assert_null (strstr (errors, "GLib-DEBUG: unsetenv()"));
  }
  else
  {
    g_assert_false (g_spawn_check_wait_status (status, &error));
    g_clear_error (&error);
  }

  if (expected_output != NULL)
    g_assert_true (strstr (output, expected_output) != NULL || strstr (errors, expected_output) != NULL);

  g_free (output);
  g_free (errors);
}

void test_ncm_cfg_new (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_free (TesNcmCfg *test, gconstpointer pdata);

void test_ncm_cfg_misc (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_fftw_planner (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_logfile_set_logstream (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_logfile_on_off (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_logfile_str_on_off (TesNcmCfg *test, gconstpointer pdata);

void test_ncm_cfg_traps (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_invalid (TesNcmCfg *test, gconstpointer pdata);
void test_ncm_cfg_string_ww (void);
void test_ncm_cfg_command_line (void);
void test_ncm_cfg_enum (void);
void test_ncm_cfg_keyfile (void);
void test_ncm_cfg_paths (void);
void test_ncm_cfg_data_filename (void);
void test_ncm_cfg_array_variant (void);
void test_ncm_cfg_version (void);

static void
_test_ncm_cfg_error_logger (const gchar *message)
{
  printf ("CAUGHT: %s\n", message);
  fflush (stdout);
}

/* The error logger receives criticals also while ncm_cfg_logfile() is off */
static void
test_ncm_cfg_error_log_handler (void)
{
  g_test_trap_subprocess ("/ncm/cfg/error_log_handler/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stdout ("*CAUGHT: a critical with the log off*");
}

static void
test_ncm_cfg_error_log_handler_subprocess (void)
{
  ncm_cfg_set_error_log_handler (&_test_ncm_cfg_error_logger);
  ncm_cfg_logfile (FALSE);
  g_log ("NUMCOSMO", G_LOG_LEVEL_CRITICAL, "a critical with the log off");
}

/* Every member is printed, also values that are negative or after a gap */
static void
test_ncm_cfg_enum_print_all_gaps (void)
{
  g_test_trap_subprocess ("/ncm/cfg/enum_print_all_gaps/subprocess", 0, 0);
  g_test_trap_assert_passed ();
  g_test_trap_assert_stdout ("*TEST_GAP_NEG*TEST_GAP_FIVE*");
}

static void
test_ncm_cfg_enum_print_all_gaps_subprocess (void)
{
  static const GEnumValue values[] = {
    {-1, "TEST_GAP_NEG", "neg"},
    {5, "TEST_GAP_FIVE", "five"},
    {0, NULL, NULL}
  };
  GType type = g_enum_register_static ("TestNcmCfgGapEnum", values);

  ncm_cfg_enum_print_all (type, "Gap enum");
}

/* A key is new once per process and per planner flag */
static void
test_ncm_cfg_fftw_plan_begin_end (void)
{
  GError *error = NULL;
  gboolean first;

  ncm_cfg_set_fftw_default_flag_str ("estimate", 10.0, &error);
  g_assert_no_error (error);

  first = ncm_cfg_fftw_plan_begin ("test_ncm_cfg_plan_%d", 7);
  ncm_cfg_fftw_plan_end (first);
  g_assert_true (first);

  first = ncm_cfg_fftw_plan_begin ("test_ncm_cfg_plan_%d", 7);
  ncm_cfg_fftw_plan_end (first);
  g_assert_false (first);

  first = ncm_cfg_fftw_plan_begin ("test_ncm_cfg_plan_%d", 8);
  ncm_cfg_fftw_plan_end (first);
  g_assert_true (first);

  ncm_cfg_set_fftw_default_flag_str ("measure", 10.0, &error);
  g_assert_no_error (error);

  first = ncm_cfg_fftw_plan_begin ("test_ncm_cfg_plan_%d", 7);
  ncm_cfg_fftw_plan_end (first);
  g_assert_true (first);

  ncm_cfg_set_fftw_default_flag_str ("estimate", 10.0, &error);
  g_assert_no_error (error);
}

/* Without NCM_FFTW_PLANNER_TIMELIMIT the limit in effect, and reported, is 10 s */
static void
test_ncm_cfg_fftw_timelimit_default (void)
{
  gchar **envp = g_get_environ ();

  envp = g_environ_unsetenv (envp, "NCM_FFTW_PLANNER");
  envp = g_environ_unsetenv (envp, "NCM_FFTW_PLANNER_TIMELIMIT");
  _test_ncm_cfg_spawn ("/ncm/cfg/fftw_timelimit_default/subprocess", envp, TRUE, NULL);
  g_strfreev (envp);
}

static void
test_ncm_cfg_fftw_timelimit_default_subprocess (void)
{
  g_assert_cmpfloat (ncm_cfg_get_fftw_timelimit (), ==, 10.0);
}

/* OMP_NUM_THREADS above OMP_THREAD_LIMIT makes ncm_cfg_init() warn; GTest makes the warning fatal */
static void
test_ncm_cfg_omp_thread_limit (void)
{
#ifdef _OPENMP
  gchar **envp = g_get_environ ();

  envp = g_environ_setenv (envp, "OMP_NUM_THREADS", "2", TRUE);
  envp = g_environ_setenv (envp, "OMP_THREAD_LIMIT", "1", TRUE);
  _test_ncm_cfg_spawn ("/ncm/cfg/omp_thread_limit/subprocess", envp, FALSE,
                       "OMP_NUM_THREADS (2) exceeds OMP_THREAD_LIMIT (1)");
  g_strfreev (envp);
#else
  g_test_skip ("built without OpenMP");
#endif /* _OPENMP */
}

static void
test_ncm_cfg_omp_thread_limit_subprocess (void)
{
}

typedef enum
{
  TEST_NCM_CFG_PATH_DEFAULT,
  TEST_NCM_CFG_PATH_XDG,
  TEST_NCM_CFG_PATH_OVERRIDE,
} TestNcmCfgPath;

static void
_test_ncm_cfg_path (TestNcmCfgPath mode)
{
  if (g_getenv ("NCM_TEST_CFG_CHILD") != NULL)
  {
    gchar *expected;
    gchar *path;
    gchar *expected_path;
    GError *error = NULL;

    ncm_cfg_init ();

    if (mode == TEST_NCM_CFG_PATH_OVERRIDE)
      expected = g_strdup (g_getenv ("NUMCOSMO_HOME"));
    else if (mode == TEST_NCM_CFG_PATH_XDG)
      expected = g_build_filename (g_getenv ("XDG_DATA_HOME"), "numcosmo", NULL);
    else
      expected = g_build_filename (g_getenv ("HOME"), ".numcosmo", NULL);

    g_assert_cmpstr (ncm_cfg_get_fullpath_base (), ==, expected);
    g_assert_true (g_file_test (expected, G_FILE_TEST_IS_DIR));

    path          = ncm_cfg_get_fullpath ("test_%d.txt", 3);
    expected_path = g_build_filename (expected, "test_3.txt", NULL);
    g_assert_cmpstr (path, ==, expected_path);
    g_assert_false (ncm_cfg_exists ("test_%d.txt", 3));
    g_assert_true (g_file_set_contents (path, "test", -1, &error));
    g_assert_no_error (error);
    g_assert_true (ncm_cfg_exists ("test_%d.txt", 3));
    g_assert_cmpint (g_remove (path), ==, 0);

    g_free (expected_path);
    g_free (path);
    g_free (expected);
  }
  else
  {
    gchar *tmp_dir = g_dir_make_tmp ("test_ncm_cfg_path_XXXXXX", NULL);
    gchar *xdg_dir;
    gchar *override_dir;
    gchar *chosen_dir;
    gchar **envp = g_get_environ ();
    const gchar *test_path;

    g_assert_nonnull (tmp_dir);
    xdg_dir      = g_build_filename (tmp_dir, "xdg", NULL);
    override_dir = g_build_filename (tmp_dir, "custom", "data", NULL);
    envp         = g_environ_setenv (envp, "HOME", tmp_dir, TRUE);
    envp         = g_environ_unsetenv (envp, "XDG_DATA_HOME");
    envp         = g_environ_unsetenv (envp, "NUMCOSMO_HOME");
    envp         = g_environ_setenv (envp, "NCM_TEST_CFG_CHILD", "1", TRUE);

    if (mode != TEST_NCM_CFG_PATH_DEFAULT)
      envp = g_environ_setenv (envp, "XDG_DATA_HOME", xdg_dir, TRUE);

    if (mode == TEST_NCM_CFG_PATH_OVERRIDE)
      envp = g_environ_setenv (envp, "NUMCOSMO_HOME", override_dir, TRUE);

    test_path = (mode == TEST_NCM_CFG_PATH_DEFAULT) ? "/ncm/cfg/path/default" :
                (mode == TEST_NCM_CFG_PATH_XDG) ? "/ncm/cfg/path/xdg" : "/ncm/cfg/path/override";
    _test_ncm_cfg_spawn (test_path, envp, TRUE, NULL);

    if (mode == TEST_NCM_CFG_PATH_OVERRIDE)
    {
      chosen_dir = override_dir;
      g_assert_cmpint (g_rmdir (chosen_dir), ==, 0);
      chosen_dir = g_build_filename (tmp_dir, "custom", NULL);
      g_assert_cmpint (g_rmdir (chosen_dir), ==, 0);
      g_free (chosen_dir);
    }
    else
    {
      chosen_dir = g_build_filename (tmp_dir, mode == TEST_NCM_CFG_PATH_XDG ? "xdg/numcosmo" : ".numcosmo", NULL);
      g_assert_cmpint (g_rmdir (chosen_dir), ==, 0);
      g_free (chosen_dir);

      if (mode == TEST_NCM_CFG_PATH_XDG)
        g_assert_cmpint (g_rmdir (xdg_dir), ==, 0);
    }

    g_assert_cmpint (g_rmdir (tmp_dir), ==, 0);

    g_strfreev (envp);
    g_free (tmp_dir);
    g_free (xdg_dir);
    g_free (override_dir);
  }
}

static void
test_ncm_cfg_path_default (void)
{
  _test_ncm_cfg_path (TEST_NCM_CFG_PATH_DEFAULT);
}

static void
test_ncm_cfg_path_xdg (void)
{
  _test_ncm_cfg_path (TEST_NCM_CFG_PATH_XDG);
}

static void
test_ncm_cfg_path_override (void)
{
  _test_ncm_cfg_path (TEST_NCM_CFG_PATH_OVERRIDE);
}

static gpointer
_test_ncm_cfg_fftw_plan_worker (gpointer data)
{
  const guint id = GPOINTER_TO_UINT (data);
  gdouble *in    = fftw_alloc_real (64);
  gdouble *out   = fftw_alloc_real (64);
  guint i;

  for (i = 0; i < 50; i++)
  {
    const guint n        = 8 + ((id + i) % 7) * 8;
    const gboolean first = ncm_cfg_fftw_plan_begin ("test_ncm_cfg_threads_redft00_%u", n);
    fftw_plan plan       = fftw_plan_r2r_1d (n, in, out, FFTW_REDFT00, ncm_cfg_get_fftw_default_flag ());

    ncm_cfg_fftw_plan_end (first);
    ncm_cfg_fftw_plan_destroy (plan);
  }

  fftw_free (in);
  fftw_free (out);

  return NULL;
}

/* Planning and destroying from several threads goes through the planning lock */
static void
test_ncm_cfg_fftw_plan_destroy (void)
{
  GThread *threads[4];
  GError *error = NULL;
  guint i;

  ncm_cfg_fftw_plan_destroy (NULL);

  ncm_cfg_set_fftw_default_flag_str ("estimate", 10.0, &error);
  g_assert_no_error (error);

  for (i = 0; i < G_N_ELEMENTS (threads); i++)
    threads[i] = g_thread_new ("fftw-plan", &_test_ncm_cfg_fftw_plan_worker, GUINT_TO_POINTER (i));

  for (i = 0; i < G_N_ELEMENTS (threads); i++)
    g_thread_join (threads[i]);
}

gint
main (gint argc, gchar *argv[])
{
  test_ncm_cfg_executable = g_find_program_in_path (argv[0]);
  g_assert_nonnull (test_ncm_cfg_executable);
  g_test_init (&argc, &argv, NULL);

  if (g_getenv ("NCM_TEST_CFG_CHILD") == NULL)
  {
    ncm_cfg_init_full_ptr (&argc, &argv);
    ncm_cfg_enable_gsl_err_handler ();
  }

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/cfg/misc", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_misc,
              &test_ncm_cfg_free);

  g_test_add ("/ncm/cfg/fftw_planner", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_fftw_planner,
              &test_ncm_cfg_free);

  g_test_add ("/ncm/cfg/logfile/set_logstream", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_logfile_set_logstream,
              &test_ncm_cfg_free);

  g_test_add ("/ncm/cfg/logfile/on_off", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_logfile_on_off,
              &test_ncm_cfg_free);

  g_test_add ("/ncm/cfg/logfile_str/on_off", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_logfile_str_on_off,
              &test_ncm_cfg_free);

  g_test_add_func ("/ncm/cfg/string_ww", &test_ncm_cfg_string_ww);
  g_test_add_func ("/ncm/cfg/command_line", &test_ncm_cfg_command_line);
  g_test_add_func ("/ncm/cfg/enum", &test_ncm_cfg_enum);
  g_test_add_func ("/ncm/cfg/fftw_plan_begin_end", &test_ncm_cfg_fftw_plan_begin_end);
  g_test_add_func ("/ncm/cfg/fftw_plan_destroy", &test_ncm_cfg_fftw_plan_destroy);
  g_test_add_func ("/ncm/cfg/fftw_timelimit_default", &test_ncm_cfg_fftw_timelimit_default);
  g_test_add_func ("/ncm/cfg/fftw_timelimit_default/subprocess", &test_ncm_cfg_fftw_timelimit_default_subprocess);
  g_test_add_func ("/ncm/cfg/omp_thread_limit", &test_ncm_cfg_omp_thread_limit);
  g_test_add_func ("/ncm/cfg/omp_thread_limit/subprocess", &test_ncm_cfg_omp_thread_limit_subprocess);
  g_test_add_func ("/ncm/cfg/path/default", &test_ncm_cfg_path_default);
  g_test_add_func ("/ncm/cfg/path/xdg", &test_ncm_cfg_path_xdg);
  g_test_add_func ("/ncm/cfg/path/override", &test_ncm_cfg_path_override);
  g_test_add_func ("/ncm/cfg/enum_print_all_gaps", &test_ncm_cfg_enum_print_all_gaps);
  g_test_add_func ("/ncm/cfg/enum_print_all_gaps/subprocess", &test_ncm_cfg_enum_print_all_gaps_subprocess);
  g_test_add_func ("/ncm/cfg/error_log_handler", &test_ncm_cfg_error_log_handler);
  g_test_add_func ("/ncm/cfg/error_log_handler/subprocess", &test_ncm_cfg_error_log_handler_subprocess);
  g_test_add_func ("/ncm/cfg/keyfile", &test_ncm_cfg_keyfile);
  g_test_add_func ("/ncm/cfg/paths", &test_ncm_cfg_paths);
  g_test_add_func ("/ncm/cfg/data_filename", &test_ncm_cfg_data_filename);
  g_test_add_func ("/ncm/cfg/array_variant", &test_ncm_cfg_array_variant);
  g_test_add_func ("/ncm/cfg/version", &test_ncm_cfg_version);

  g_test_add ("/ncm/cfg/traps", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_traps,
              &test_ncm_cfg_free);

  g_test_add ("/ncm/cfg/logfile/subprocess", TesNcmCfg, NULL,
              &test_ncm_cfg_new,
              &test_ncm_cfg_invalid,
              &test_ncm_cfg_free);

  return g_test_run ();
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

    gchar *expected = g_build_filename (ncm_cfg_get_fullpath_base (), "test_full_path_1.txt", NULL);

    g_assert_cmpstr (full_path, ==, expected);

    g_free (expected);
    g_free (full_path);
  }

  /* Test string to comment */
  {
    gchar *comment = ncm_cfg_string_to_comment ("test string to comment, this is a very long comment that "
                                                "should be truncated to 80 characters, but it is not, so it "
                                                "will be truncated by the user.");

    g_assert_true (g_str_has_prefix (comment, "###############################################################################\n\n  test string"));
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
test_ncm_cfg_fftw_planner (TesNcmCfg *test, gconstpointer pdata)
{
  const gchar *flags[] = {"estimate", "measure", "patient", "exhaustive"};
  GError *error        = NULL;
  guint i;

  /* String round-trip for every planner flag. */
  for (i = 0; i < G_N_ELEMENTS (flags); i++)
  {
    ncm_cfg_set_fftw_default_flag_str (flags[i], 60.0, &error);
    g_assert_no_error (error);
    g_assert_cmpstr (ncm_cfg_get_fftw_default_flag_str (), ==, flags[i]);
    g_assert_cmpfloat (ncm_cfg_get_fftw_timelimit (), ==, 60.0);
  }

  /* Fallback path used by the -Dfftw-planner build option: with NCM_FFTW_PLANNER
   * unset, the fallback string must take effect. Guards against silent breakage of
   * the compile-time default plumbing (NUMCOSMO_FFTW_PLAN). */
  g_unsetenv ("NCM_FFTW_PLANNER");
  g_unsetenv ("NCM_FFTW_PLANNER_TIMELIMIT");
  ncm_cfg_set_fftw_default_from_env_str ("estimate", 30.0, &error);
  g_assert_no_error (error);
  g_assert_cmpstr (ncm_cfg_get_fftw_default_flag_str (), ==, "estimate");
  g_assert_cmpfloat (ncm_cfg_get_fftw_timelimit (), ==, 30.0);

  /* Environment override wins over the fallback. */
  g_setenv ("NCM_FFTW_PLANNER", "measure", TRUE);
  g_setenv ("NCM_FFTW_PLANNER_TIMELIMIT", "12.5", TRUE);
  ncm_cfg_set_fftw_default_from_env_str ("estimate", 30.0, &error);
  g_assert_no_error (error);
  g_assert_cmpstr (ncm_cfg_get_fftw_default_flag_str (), ==, "measure");
  g_assert_cmpfloat (ncm_cfg_get_fftw_timelimit (), ==, 12.5);
  g_unsetenv ("NCM_FFTW_PLANNER");
  g_unsetenv ("NCM_FFTW_PLANNER_TIMELIMIT");

  /* Invalid planner string must raise an error and leave the flag unchanged. */
  ncm_cfg_set_fftw_default_flag_str ("nonsense", 60.0, &error);
  g_assert_error (error, NCM_CFG_ERROR, NCM_CFG_ERROR_INVALID_FFTW_FLAG_STRING);
  g_clear_error (&error);

  /* Restore the build default for the remaining tests. */
  ncm_cfg_set_fftw_default_flag_str ("estimate", 60.0, NULL);
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
  g_test_trap_assert_stderr ("*This message should be printed in stderr 1*");
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
  g_test_trap_assert_stderr ("*This message should be printed in stderr 1*");
}

void
test_ncm_cfg_logfile_str_on_off (TesNcmCfg *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    ncm_cfg_set_logstream (stderr);

    ncm_cfg_logfile (FALSE);
    ncm_message_str ("This message should not be printed in stderr 1");
    ncm_cfg_logfile (TRUE);
    ncm_message_str ("This message should be printed in stderr 1");

    return;
  }

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_stderr ("*This message should be printed in stderr 1*");
}

void
test_ncm_cfg_traps (TesNcmCfg *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/cfg/logfile/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_cfg_invalid (TesNcmCfg *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

void
test_ncm_cfg_string_ww (void)
{
  gchar *ww = ncm_string_ww ("aaa bbb ccc ddd eee", "> ", "  ", 12);

  /* The first line holds as many words as fit in ncols minus the first prefix */
  g_assert_cmpstr (ww, ==, "> aaa bbb\n   ccc ddd\n   eee\n");
  g_free (ww);
}

void
test_ncm_cfg_command_line (void)
{
  gchar *argv[] = {"numcosmo", "run", "two words", "-x"};
  gchar *cmd    = ncm_cfg_command_line (argv, G_N_ELEMENTS (argv));

  g_assert_cmpstr (cmd, ==, "numcosmo run 'two words' -x");
  g_free (cmd);
}

void
test_ncm_cfg_enum (void)
{
  const GType t = NCM_TYPE_CFG_ERROR;

  g_assert_cmpint (ncm_cfg_get_enum_by_id_name_nick (t, "1")->value, ==, NCM_CFG_ERROR_INVALID_FFTW_FLAG_STRING);
  g_assert_cmpint (ncm_cfg_get_enum_by_id_name_nick (t, "NCM_CFG_ERROR_INVALID_FFTW_TIMELIMIT")->value, ==, NCM_CFG_ERROR_INVALID_FFTW_TIMELIMIT);
  g_assert_cmpint (ncm_cfg_get_enum_by_id_name_nick (t, "flag")->value, ==, NCM_CFG_ERROR_INVALID_FFTW_FLAG);
  g_assert_null (ncm_cfg_get_enum_by_id_name_nick (t, "not-a-nick"));
  g_assert_null (ncm_cfg_get_enum_by_id_name_nick (t, "7"));

  g_assert_cmpstr (ncm_cfg_enum_get_value (t, 2)->value_nick, ==, "timelimit");
  g_assert_null (ncm_cfg_enum_get_value (t, 7));
}

void
test_ncm_cfg_keyfile (void)
{
  gboolean flag          = TRUE;
  gint number            = 42;
  gdouble x              = 2.5;
  gchar *name            = g_strdup ("value");
  gchar **list           = g_strsplit ("a,b", ",", -1);
  GOptionEntry entries[] = {
    {"flag", 0, 0, G_OPTION_ARG_NONE, &flag, "A flag", NULL},
    {"number", 0, 0, G_OPTION_ARG_INT, &number, "A number", NULL},
    {"x", 0, 0, G_OPTION_ARG_DOUBLE, &x, "A double", NULL},
    {"name", 0, 0, G_OPTION_ARG_STRING, &name, "A string", NULL},
    {"list", 0, 0, G_OPTION_ARG_STRING_ARRAY, &list, "A list", NULL},
    { NULL, 0, 0, 0, NULL, NULL, NULL },
  };
  GKeyFile *kfile = g_key_file_new ();
  gchar *argv[16];
  gchar *argv_owned[16];
  gint argc     = 1;
  GError *error = NULL;

  ncm_cfg_entries_to_keyfile (kfile, "group", entries);
  g_assert_cmpint (g_key_file_get_integer (kfile, "group", "number", NULL), ==, 42);
  g_assert_true (g_key_file_get_boolean (kfile, "group", "flag", NULL));

  /* Keyfile back to arguments: parsing them restores the values */
  argv[0] = g_strdup ("prog");
  ncm_cfg_keyfile_to_arg (kfile, "group", entries, argv, &argc);
  g_assert_cmpint (argc, ==, 1 + 1 + 2 + 2 + 2 + 4);

  /* g_option_context_parse() reorders argv, keep the strings to free them */
  memcpy (argv_owned, argv, argc * sizeof (gchar *));

  flag   = FALSE;
  number = 0;
  x      = 0.0;
  g_clear_pointer (&name, g_free);
  g_clear_pointer (&list, g_strfreev);

  {
    GOptionContext *ctx = g_option_context_new (NULL);
    gchar **argv_p      = argv;

    g_option_context_add_main_entries (ctx, entries, NULL);
    g_assert_true (g_option_context_parse (ctx, &argc, &argv_p, &error));
    g_assert_no_error (error);
    g_option_context_free (ctx);
  }

  g_assert_true (flag);
  g_assert_cmpint (number, ==, 42);
  g_assert_cmpfloat (x, ==, 2.5);
  g_assert_cmpstr (name, ==, "value");
  g_assert_cmpuint (g_strv_length (list), ==, 2);
  g_assert_cmpstr (list[1], ==, "b");

  {
    gint i;

    for (i = 0; i < 12; i++)
      g_free (argv_owned[i]);
  }

  g_free (name);
  g_strfreev (list);
  g_key_file_unref (kfile);
}

void
test_ncm_cfg_paths (void)
{
  const gchar *base = ncm_cfg_get_fullpath_base ();
  gchar *path       = ncm_cfg_get_fullpath ("sub_%d.txt", 3);
  gchar *expected   = g_build_filename (base, "sub_3.txt", NULL);

  g_assert_true (g_file_test (base, G_FILE_TEST_IS_DIR));
  g_assert_cmpstr (path, ==, expected);
  g_assert_false (ncm_cfg_exists ("test_ncm_cfg_file_that_does_not_exist_%d", 17));

  g_free (path);
  g_free (expected);
}

void
test_ncm_cfg_data_filename (void)
{
  gchar *data_dir = ncm_cfg_get_data_directory ();
  gchar *found    = ncm_cfg_get_data_filename ("BBN_spline2d.obj", TRUE);
  gchar *expected = g_build_filename (data_dir, "BBN_spline2d.obj", NULL);

  g_assert_true (g_file_test (data_dir, G_FILE_TEST_IS_DIR));
  g_assert_cmpstr (found, ==, expected);
  g_assert_null (ncm_cfg_get_data_filename ("test_ncm_cfg_no_such_data_file", FALSE));

  g_free (data_dir);
  g_free (found);
  g_free (expected);
}

void
test_ncm_cfg_array_variant (void)
{
  GArray *a            = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *b            = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const gdouble vals[] = {1.0, -2.0, 3.5};
  GVariant *var;

  g_array_append_vals (a, vals, 3);
  var = ncm_cfg_array_to_variant (a, G_VARIANT_TYPE_DOUBLE);
  g_assert_true (g_variant_is_of_type (var, G_VARIANT_TYPE ("ad")));

  ncm_cfg_array_set_variant (b, var);
  g_assert_cmpuint (b->len, ==, 3);
  g_assert_cmpfloat (g_array_index (b, gdouble, 2), ==, 3.5);

  g_variant_unref (var);
  g_array_unref (a);
  g_array_unref (b);
}

void
test_ncm_cfg_version (void)
{
  guint major, minor, micro;
  const guint v = ncm_cfg_get_version (&major, &minor, &micro);
  gchar *vstr   = ncm_cfg_get_version_string ();
  gchar *vexp   = g_strdup_printf ("%u.%u.%u", major, minor, micro);

  g_assert_cmpuint (v, ==, 10000 * major + 100 * minor + micro);
  g_assert_cmpstr (vstr, ==, vexp);
  g_assert_true (ncm_cfg_version_check (major, minor, micro));
  g_assert_false (ncm_cfg_version_check (major, minor, micro + 1));
  g_assert_true (ncm_cfg_version_check (0, 0, 0));
  g_assert_nonnull (ncm_cfg_get_commit_hash ());

  g_free (vstr);
  g_free (vexp);
}

