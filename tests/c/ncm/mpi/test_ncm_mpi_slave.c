/***************************************************************************
 *            test_ncm_mpi_slave.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mpi_slave.c
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
 * The worker side of the NcmMPIJob protocol on its own. The test does not call
 * ncm_cfg_init(), so no rank enters the worker loop by itself: rank 1 calls
 * ncm_mpi_slave_serve_job() and rank 0 speaks the protocol with raw MPI calls. After
 * each job the worker sends rank 0 what ncm_mpi_slave_serve_job() returned.
 *
 * The same program checks the worker's refusals. Each one aborts the whole MPI job, so
 * `--abort=<case>` runs a single one under mpiexec, and `--spawn-aborts`, run without
 * mpiexec, launches each case through the mpiexec in NCM_TEST_MPIEXEC and checks the
 * exit status and the message.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>
#include <mpi.h>
#include <string.h>

#include "ncm/mpi/ncm_mpi_slave.h"
#include "test_ncm_mpi_job_shape.h"

#define TEST_WORKER (1)
#define TEST_TAG_REPORT (1000)

static gboolean _test_killed = FALSE;

static void
_test_send_cmd (gint cmd)
{
  MPI_Send (&cmd, 1, MPI_INT, TEST_WORKER, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD);
}

static void
_test_send_job (NcmMPIJob *mpi_job)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  GVariant *job_ser = ncm_serialize_to_variant (ser, G_OBJECT (mpi_job));

  _test_send_cmd (NCM_MPI_CTRL_SLAVE_INIT);
  MPI_Send (g_variant_get_data (job_ser), g_variant_get_size (job_ser), MPI_BYTE, TEST_WORKER, NCM_MPI_CTRL_TAG_JOB, MPI_COMM_WORLD);

  g_variant_unref (job_ser);
  ncm_serialize_free (ser);
}

static void
_test_send_work (const gdouble *input, gint input_len)
{
  _test_send_cmd (NCM_MPI_CTRL_SLAVE_WORK);
  MPI_Send (input, input_len, MPI_DOUBLE, TEST_WORKER, NCM_MPI_CTRL_TAG_WORK_INPUT, MPI_COMM_WORLD);
}

static void
_test_recv_return (gdouble *ret, gint return_len)
{
  MPI_Status status;
  gint count = 0;

  MPI_Recv (ret, return_len, MPI_DOUBLE, TEST_WORKER, NCM_MPI_CTRL_TAG_WORK_RETURN, MPI_COMM_WORLD, &status);
  MPI_Get_count (&status, MPI_DOUBLE, &count);
  g_assert_cmpint (count, ==, return_len);
}

static gint
_test_recv_report (void)
{
  gint report = -1;

  MPI_Recv (&report, 1, MPI_INT, TEST_WORKER, TEST_TAG_REPORT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  return report;
}

static NcmMPIJob *
_test_shape_new (guint input_len, guint return_len)
{
  return g_object_new (TEST_TYPE_MPI_JOB_SHAPE, "input-len", input_len, "return-len", return_len, NULL);
}

/* The return of the shape job at input i with entries i + j / 2. */
static void
_test_shape_input (gdouble *input, guint input_len, guint i)
{
  guint j;

  for (j = 0; j < input_len; j++)
    input[j] = i + 0.5 * j;
}

static void
_test_shape_assert_return (const gdouble *ret, guint input_len, guint return_len, guint i)
{
  gdouble s = 0.0;
  guint j, k;

  for (j = 0; j < input_len; j++)
    s += (j + 1.0) * (i + 0.5 * j);

  for (k = 0; k < return_len; k++)
    g_assert_cmpfloat (ret[k], ==, (k + 1.0) * s);
}

/* Sends n inputs one at a time, each followed by its return. */
static void
_test_run_shape (guint input_len, guint return_len, guint n)
{
  NcmMPIJob *mpi_job = _test_shape_new (input_len, return_len);
  gdouble *input     = g_new (gdouble, input_len);
  gdouble *ret       = g_new (gdouble, return_len);
  guint i;

  _test_send_job (mpi_job);

  for (i = 0; i < n; i++)
  {
    _test_shape_input (input, input_len, i);
    _test_send_work (input, input_len);
    _test_recv_return (ret, return_len);
    _test_shape_assert_return (ret, input_len, return_len, i);
  }

  _test_send_cmd (NCM_MPI_CTRL_SLAVE_FREE);
  g_assert_cmpint (_test_recv_report (), ==, TRUE);

  g_free (input);
  g_free (ret);
  ncm_mpi_job_free (mpi_job);
}

static void
test_ncm_mpi_slave_work (void)
{
  _test_run_shape (3, 5, 20);
}

static void
test_ncm_mpi_slave_reinit (void)
{
  /* After FREE the worker takes a new job, of another type or shape. */
  NcmRNG *rng        = ncm_rng_seeded_new (NULL, 20260930);
  NcmMPIJobTest *mjt = ncm_mpi_job_test_new ();
  NcmVector *vec     = NULL;
  guint i;

  _test_run_shape (3, 5, 4);

  ncm_mpi_job_test_set_rand_vector (mjt, 10, rng);
  g_object_get (mjt, "vector", &vec, NULL);
  _test_send_job (NCM_MPI_JOB (mjt));

  for (i = 0; i < 10; i++)
  {
    const gdouble input = 9 - i;
    gdouble ret         = 0.0;

    _test_send_work (&input, 1);
    _test_recv_return (&ret, 1);
    g_assert_cmpfloat (ret, ==, ncm_vector_get (vec, 9 - i));
  }

  _test_send_cmd (NCM_MPI_CTRL_SLAVE_FREE);
  g_assert_cmpint (_test_recv_report (), ==, TRUE);

  _test_run_shape (7, 2, 4);

  ncm_vector_free (vec);
  ncm_mpi_job_test_free (mjt);
  ncm_rng_free (rng);
}

static void
test_ncm_mpi_slave_pending (void)
{
  /*
   * Returns far above the eager limit stay pending until rank 0 posts their receives,
   * which it does only after FREE: the worker must send them all before it releases
   * the job and reports.
   */
  const guint input_len  = 2;
  const guint return_len = 200000;
  const guint n          = 3;
  NcmMPIJob *mpi_job     = _test_shape_new (input_len, return_len);
  gdouble *ret           = g_new (gdouble, return_len);
  guint i;

  _test_send_job (mpi_job);

  for (i = 0; i < n; i++)
  {
    gdouble input[2];

    _test_shape_input (input, input_len, i);
    _test_send_work (input, input_len);
  }

  _test_send_cmd (NCM_MPI_CTRL_SLAVE_FREE);

  /* No report while the returns wait for their receives. */
  for (i = 0; i < 20; i++)
  {
    gint flag = 0;

    MPI_Iprobe (TEST_WORKER, TEST_TAG_REPORT, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
    g_assert_false (flag);
    g_usleep (10000);
  }

  for (i = 0; i < n; i++)
  {
    _test_recv_return (ret, return_len);
    _test_shape_assert_return (ret, input_len, return_len, i);
  }

  g_assert_cmpint (_test_recv_report (), ==, TRUE);

  g_free (ret);
  ncm_mpi_job_free (mpi_job);
}

static void
test_ncm_mpi_slave_kill (void)
{
  /* KILL ends the job, and the worker reports that it stops. */
  NcmMPIJob *mpi_job  = _test_shape_new (1, 1);
  const gdouble input = 1.0;
  gdouble ret         = 0.0;

  _test_send_job (mpi_job);
  _test_send_work (&input, 1);
  _test_recv_return (&ret, 1);

  _test_send_cmd (NCM_MPI_CTRL_SLAVE_KILL);
  _test_killed = TRUE;
  g_assert_cmpint (_test_recv_report (), ==, FALSE);

  ncm_mpi_job_free (mpi_job);
}

/* The refusals: each case aborts the worker with the message it is checked against. */
typedef struct _TestAbortCase
{
  const gchar *name;
  const gchar *message;
} TestAbortCase;

static const TestAbortCase abort_cases[] =
{
  {"init_twice",       "already initialized"},
  {"work_before_init", "received work"},
  {"unknown_command",  "unknown MPI message"},
};

static void
_test_abort_case_run (const gchar *name)
{
  NcmMPIJob *mpi_job = _test_shape_new (1, 1);
  gint i;

  if (g_strcmp0 (name, "init_twice") == 0)
  {
    _test_send_job (mpi_job);
    _test_send_job (mpi_job);
  }
  else if (g_strcmp0 (name, "work_before_init") == 0)
  {
    _test_send_cmd (NCM_MPI_CTRL_SLAVE_WORK);
  }
  else if (g_strcmp0 (name, "unknown_command") == 0)
  {
    _test_send_cmd (NCM_MPI_CTRL_SLAVE_LEN + 7);
  }
  else
  {
    g_error ("unknown abort case `%s'.", name);
  }

  /* The worker should abort the MPI job; if it has not after ten seconds, stop it and
   * exit normally, which the driver reports as a failure. */
  for (i = 0; i < 1000; i++)
  {
    gint flag = 0;

    MPI_Iprobe (TEST_WORKER, TEST_TAG_REPORT, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);

    if (flag)
      break;

    g_usleep (10000);
  }

  ncm_mpi_job_free (mpi_job);
}

static void
test_ncm_mpi_slave_abort (gconstpointer pdata)
{
  const TestAbortCase *ac = pdata;
  const gchar *mpiexec    = g_getenv ("NCM_TEST_MPIEXEC");
  const gchar *self       = g_getenv ("NCM_TEST_SELF");
  gchar *case_arg         = g_strdup_printf ("--abort=%s", ac->name);
  gchar *argv[]           = {(gchar *) mpiexec, (gchar *) "-n", (gchar *) "2", (gchar *) self, case_arg, NULL};
  gchar *err              = NULL;
  gint wait_status        = 0;
  GError *error           = NULL;

  g_assert_nonnull (mpiexec);
  g_assert_nonnull (self);

  g_assert_true (g_spawn_sync (NULL, argv, NULL, G_SPAWN_STDOUT_TO_DEV_NULL, NULL, NULL, NULL, &err, &wait_status, &error));
  g_assert_no_error (error);

#if GLIB_CHECK_VERSION (2, 70, 0)
  g_assert_false (g_spawn_check_wait_status (wait_status, NULL));
#else
  g_assert_false (g_spawn_check_exit_status (wait_status, NULL));
#endif /* GLIB_CHECK_VERSION(2,70,0) */
  g_assert_nonnull (strstr (err, ac->message));

  g_free (err);
  g_free (case_arg);
}

static gint
_test_spawn_aborts (gint argc, gchar *argv[])
{
  guint i;

  g_test_init (&argc, &argv, NULL);

  for (i = 0; i < G_N_ELEMENTS (abort_cases); i++)
  {
    gchar *path = g_strdup_printf ("/ncm/mpi/slave/abort/%s", abort_cases[i].name);

    g_test_add_data_func (path, &abort_cases[i], &test_ncm_mpi_slave_abort);
    g_free (path);
  }

  return g_test_run ();
}

gint
main (gint argc, gchar *argv[])
{
  const gchar *abort_case = NULL;
  gint rank               = 0;
  gint i;

  for (i = 1; i < argc; i++)
  {
    if (g_strcmp0 (argv[i], "--spawn-aborts") == 0)
    {
      argv[i] = argv[argc - 1];
      argc--;

      return _test_spawn_aborts (argc, argv);
    }

    if (g_str_has_prefix (argv[i], "--abort="))
      abort_case = argv[i] + strlen ("--abort=");
  }

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  ncm_cfg_register_objects ();
  ncm_cfg_register_obj (TEST_TYPE_MPI_JOB_SHAPE);

  if (rank == TEST_WORKER)
  {
    gint more = TRUE;

    while (more)
    {
      more = ncm_mpi_slave_serve_job ();
      MPI_Send (&more, 1, MPI_INT, 0, TEST_TAG_REPORT, MPI_COMM_WORLD);
    }
  }
  else if (rank == 0)
  {
    if (abort_case != NULL)
    {
      _test_abort_case_run (abort_case);
    }
    else
    {
      g_test_init (&argc, &argv, NULL);

      g_test_add_func ("/ncm/mpi/slave/work", &test_ncm_mpi_slave_work);
      g_test_add_func ("/ncm/mpi/slave/reinit", &test_ncm_mpi_slave_reinit);
      g_test_add_func ("/ncm/mpi/slave/pending", &test_ncm_mpi_slave_pending);
      g_test_add_func ("/ncm/mpi/slave/kill", &test_ncm_mpi_slave_kill);

      g_test_run ();
    }

    /* A filtered run may not reach the kill test; the worker still has to stop. */
    if (!_test_killed)
    {
      _test_send_cmd (NCM_MPI_CTRL_SLAVE_KILL);
      _test_recv_report ();
    }
  }

  MPI_Finalize ();

  return 0;
}

