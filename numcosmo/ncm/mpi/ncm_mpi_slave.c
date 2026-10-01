/***************************************************************************
 *            ncm_mpi_slave.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mpi_slave.c
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
 * The worker side of the NcmMPIJob protocol, see #NcmMPIJobCtrlTag for the messages. A
 * worker serves one job at a time: NCM_MPI_CTRL_SLAVE_INIT brings the serialized job,
 * each NCM_MPI_CTRL_SLAVE_WORK one input whose return is sent back without waiting, and
 * NCM_MPI_CTRL_SLAVE_FREE or NCM_MPI_CTRL_SLAVE_KILL ends the job once every return
 * has been sent. After FREE the worker waits for the next job; after KILL it stops.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/mpi/ncm_mpi_slave.h"
#include "ncm/mpi/ncm_mpi_job.h"
#include "ncm/core/ncm_serialize.h"

#ifdef HAVE_MPI
#include <mpi.h>

extern NcmMPIJobCtrl _mpi_ctrl;

typedef struct _NcmMPISlaveReturn
{
  gpointer obj;
  gpointer buf;
} NcmMPISlaveReturn;

typedef struct _NcmMPISlave
{
  NcmSerialize *ser;
  NcmMPIJob *mpi_job;
  gpointer input;
  gpointer input_buf;
  gint input_len;
  gint input_size;
  gint return_len;
  gint return_size;
  MPI_Datatype input_dtype;
  MPI_Datatype return_dtype;
  GArray *ret_requests; /* MPI_Request of the returns still being sent */
  GArray *rets;         /* NcmMPISlaveReturn, in the order of ret_requests */
} NcmMPISlave;

static void _ncm_mpi_slave_init_job (NcmMPISlave *slave);
static void _ncm_mpi_slave_work (NcmMPISlave *slave);
static void _ncm_mpi_slave_release_sent (NcmMPISlave *slave, gboolean wait);
static void _ncm_mpi_slave_release_return (NcmMPISlave *slave, guint i);
static void _ncm_mpi_slave_clear (NcmMPISlave *slave);

#endif /* HAVE_MPI */

/*
 * ncm_mpi_slave_run:
 *
 * Serves jobs until NCM_MPI_CTRL_SLAVE_KILL.
 */
void
ncm_mpi_slave_run (void)
{
#ifdef HAVE_MPI
  NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Starting slave!\n", _mpi_ctrl.size, _mpi_ctrl.rank);

  while (ncm_mpi_slave_serve_job ())
    ;

  NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Dying slave!\n", _mpi_ctrl.size, _mpi_ctrl.rank);
#else
  g_error ("ncm_mpi_slave_run: MPI unsupported.");

#endif /* HAVE_MPI */
}

/*
 * ncm_mpi_slave_serve_job:
 *
 * Serves one job, from the commands before its NCM_MPI_CTRL_SLAVE_INIT to its
 * NCM_MPI_CTRL_SLAVE_FREE or NCM_MPI_CTRL_SLAVE_KILL, and releases it.
 *
 * Returns: FALSE after NCM_MPI_CTRL_SLAVE_KILL.
 */
gboolean
ncm_mpi_slave_serve_job (void)
{
#ifdef HAVE_MPI
  NcmMPISlave slave = {
    ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP),
    NULL, NULL, NULL,
    0, 0, 0, 0,
    MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
    g_array_new (FALSE, FALSE, sizeof (MPI_Request)),
    g_array_new (FALSE, TRUE, sizeof (NcmMPISlaveReturn)),
  };
  gboolean end  = FALSE;
  gboolean kill = FALSE;

  while (!end)
  {
    MPI_Status status;
    gint cmd = 0;

    NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Waiting for command...\n", _mpi_ctrl.size, _mpi_ctrl.rank);
    MPI_Recv (&cmd, 1, MPI_INT, NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD, &status);
    NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Received %d\n", _mpi_ctrl.size, _mpi_ctrl.rank, cmd);

    switch (cmd)
    {
      case NCM_MPI_CTRL_SLAVE_INIT:
        _ncm_mpi_slave_init_job (&slave);
        break;
      case NCM_MPI_CTRL_SLAVE_WORK:
        _ncm_mpi_slave_work (&slave);
        break;
      case NCM_MPI_CTRL_SLAVE_FREE:
        end = TRUE;
        break;
      case NCM_MPI_CTRL_SLAVE_KILL:
        end  = TRUE;
        kill = TRUE;
        break;
      /* LCOV_EXCL_START */
      default:
        g_error ("ncm_mpi_slave_serve_job: unknown MPI message `%d' to slave %d", cmd, _mpi_ctrl.rank);
        break;
        /* LCOV_EXCL_STOP */
    }

    _ncm_mpi_slave_release_sent (&slave, FALSE);
  }

  _ncm_mpi_slave_clear (&slave);

  return !kill;

#else
  g_error ("ncm_mpi_slave_serve_job: MPI unsupported.");

  return FALSE;

#endif /* HAVE_MPI */
}

#ifdef HAVE_MPI

static void
_ncm_mpi_slave_init_job (NcmMPISlave *slave)
{
  MPI_Status status;
  GVariant *job_ser = NULL;
  gint job_size     = 0;
  gint job_recv     = 0;
  gchar *job        = NULL;

  if (slave->mpi_job != NULL)
    g_error ("ncm_mpi_slave_serve_job: slave %d already initialized.", _mpi_ctrl.rank);

  MPI_Probe (NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_JOB, MPI_COMM_WORLD, &status);
  MPI_Get_count (&status, MPI_BYTE, &job_size);

  job = g_new (gchar, job_size);
  MPI_Recv (job, job_size, MPI_BYTE, NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_JOB, MPI_COMM_WORLD, &status);
  MPI_Get_count (&status, MPI_BYTE, &job_recv);

  NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Slave object received size %d.\n", _mpi_ctrl.size, _mpi_ctrl.rank, job_recv);

  g_assert_cmpint (job_recv, ==, job_size);

  job_ser        = g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE), job, job_size, TRUE, g_free, job);
  slave->mpi_job = NCM_MPI_JOB (ncm_serialize_from_variant (slave->ser, job_ser));

  g_variant_unref (job_ser);

  ncm_mpi_job_work_init (slave->mpi_job);

  slave->input_dtype  = ncm_mpi_job_input_datatype (slave->mpi_job, &slave->input_len, &slave->input_size);
  slave->return_dtype = ncm_mpi_job_return_datatype (slave->mpi_job, &slave->return_len, &slave->return_size);
  slave->input        = ncm_mpi_job_create_input (slave->mpi_job);
  slave->input_buf    = ncm_mpi_job_get_input_buffer (slave->mpi_job, slave->input);
}

static void
_ncm_mpi_slave_work (NcmMPISlave *slave)
{
  NcmMPISlaveReturn r;
  MPI_Request request;
  MPI_Status status;
  gint input_recv = 0;

  if (slave->mpi_job == NULL)
    g_error ("ncm_mpi_slave_serve_job: uninitialized slave `%d' received work.", _mpi_ctrl.rank);

  MPI_Recv (slave->input_buf, slave->input_len, slave->input_dtype, NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_WORK_INPUT, MPI_COMM_WORLD, &status);
  MPI_Get_count (&status, slave->input_dtype, &input_recv);

  g_assert_cmpint (input_recv, ==, slave->input_len);

  ncm_mpi_job_unpack_input (slave->mpi_job, slave->input_buf, slave->input);

  r.obj = ncm_mpi_job_create_return (slave->mpi_job);
  ncm_mpi_job_run (slave->mpi_job, slave->input, r.obj);
  r.buf = ncm_mpi_job_pack_return (slave->mpi_job, r.obj);

  MPI_Isend (r.buf, slave->return_len, slave->return_dtype, NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_WORK_RETURN, MPI_COMM_WORLD, &request);

  g_array_append_val (slave->ret_requests, request);
  g_array_append_val (slave->rets, r);
}

/* Releases the returns whose sends completed, or waits for all of them. */
static void
_ncm_mpi_slave_release_sent (NcmMPISlave *slave, gboolean wait)
{
  gint i;

  if (wait && (slave->ret_requests->len > 0))
    MPI_Waitall (slave->ret_requests->len, (MPI_Request *) slave->ret_requests->data, MPI_STATUSES_IGNORE);

  for (i = (gint) slave->ret_requests->len - 1; i >= 0; i--)
  {
    gint done = wait;

    if (!wait)
      MPI_Test (&g_array_index (slave->ret_requests, MPI_Request, i), &done, MPI_STATUS_IGNORE);

    if (done)
      _ncm_mpi_slave_release_return (slave, i);
  }
}

static void
_ncm_mpi_slave_release_return (NcmMPISlave *slave, guint i)
{
  NcmMPISlaveReturn r = g_array_index (slave->rets, NcmMPISlaveReturn, i);

  ncm_mpi_job_destroy_return_buffer (slave->mpi_job, r.obj, r.buf);
  ncm_mpi_job_destroy_return (slave->mpi_job, r.obj);

  g_array_remove_index_fast (slave->ret_requests, i);
  g_array_remove_index_fast (slave->rets, i);
}

static void
_ncm_mpi_slave_clear (NcmMPISlave *slave)
{
  _ncm_mpi_slave_release_sent (slave, TRUE);

  if (slave->input_buf != NULL)
    ncm_mpi_job_destroy_input_buffer (slave->mpi_job, slave->input, slave->input_buf);

  if (slave->input != NULL)
    ncm_mpi_job_destroy_input (slave->mpi_job, slave->input);

  ncm_mpi_job_clear (&slave->mpi_job);

  g_array_unref (slave->ret_requests);
  g_array_unref (slave->rets);
  ncm_serialize_free (slave->ser);
}

#endif /* HAVE_MPI */

