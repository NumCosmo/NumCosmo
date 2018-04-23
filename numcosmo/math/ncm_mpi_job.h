/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_mpi_job.h
 *
 *  Sat April 21 18:33:41 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mpi_job.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_MPI_JOB_H_
#define _NCM_MPI_JOB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_serialize.h>

G_BEGIN_DECLS

#define NCM_TYPE_MPI_JOB             (ncm_mpi_job_get_type ())
#define NCM_MPI_JOB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MPI_JOB, NcmMPIJob))
#define NCM_MPI_JOB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MPI_JOB, NcmMPIJobClass))
#define NCM_IS_MPI_JOB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MPI_JOB))
#define NCM_IS_MPI_JOB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MPI_JOB))
#define NCM_MPI_JOB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MPI_JOB, NcmMPIJobClass))

typedef struct _NcmMPIJobClass NcmMPIJobClass;
typedef struct _NcmMPIJob NcmMPIJob;
typedef struct _NcmMPIJobPrivate NcmMPIJobPrivate;

struct _NcmMPIJobClass
{
	/*< private >*/
	GObjectClass parent_class;
	void *(*work_init) (NcmMPIJob *mpi_job);
	void *(*work_clear) (NcmMPIJob *mpi_job);
	GObject *(*run) (NcmMPIJob *mpi_job, GObject *input);
	NcmVector *(*run_vector) (NcmMPIJob *mpi_job, NcmVector *input);
};

struct _NcmMPIJob
{
	/*< private >*/
	GObject parent_instance;
	NcmMPIJobPrivate *priv;
};

/**************************************************************************************/

typedef struct _NcmMPIJobCtrl 
{
	gint rank;
	gint size;
	gint nslaves;
	gint working_slaves;
} NcmMPIJobCtrl;

#define NCM_MPI_CTRL_MASTER_ID (0)

/**
 * NcmMPIJobCtrlMsg:
 * @NCM_MPI_CTRL_SLAVE_INIT: FIXME
 * @NCM_MPI_CTRL_SLAVE_FREE: FIXME
 * @NCM_MPI_CTRL_SLAVE_WORK: FIXME
 * @NCM_MPI_CTRL_SLAVE_WORK_VECTOR: FIXME
 * 
 * FIXME
 * 
 */ 
typedef enum _NcmMPIJobCtrlMsg 
{
	NCM_MPI_CTRL_SLAVE_INIT = 0,
	NCM_MPI_CTRL_SLAVE_FREE,
	NCM_MPI_CTRL_SLAVE_WORK,
	NCM_MPI_CTRL_SLAVE_WORK_VECTOR,
	/* < private > */
	NCM_MPI_CTRL_SLAVE_LEN, /*< skip >*/
} NcmMPIJobCtrlMsg;

/**
 * NcmMPIJobCtrlTag:
 * @NCM_MPI_CTRL_TAG_CMD: FIXME
 * @NCM_MPI_CTRL_TAG_JOB: FIXME
 * @NCM_MPI_CTRL_TAG_WORK_INPUT: FIXME
 * @NCM_MPI_CTRL_TAG_WORK_RETURN: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcmMPIJobCtrlTag
{
	NCM_MPI_CTRL_TAG_CMD = 0,
	NCM_MPI_CTRL_TAG_JOB,
	NCM_MPI_CTRL_TAG_WORK_INPUT,
	NCM_MPI_CTRL_TAG_WORK_RETURN,
	/* < private > */
	NCM_MPI_CTRL_TAG_LEN, /*< skip >*/
} NcmMPIJobCtrlTag;

/**************************************************************************************/

GType ncm_mpi_job_get_type (void) G_GNUC_CONST;

NcmMPIJob *ncm_mpi_job_ref (NcmMPIJob *mpi_job);

void ncm_mpi_job_free (NcmMPIJob *mpi_job);
void ncm_mpi_job_clear (NcmMPIJob **mpi_job);

void ncm_mpi_job_set_vec_sizes (NcmMPIJob *mpi_job, const guint ret_size, const guint in_size);
guint ncm_mpi_job_get_input_vec_size (NcmMPIJob *mpi_job);
guint ncm_mpi_job_get_return_vec_size (NcmMPIJob *mpi_job);

void ncm_mpi_job_work_init (NcmMPIJob *mpi_job);
void ncm_mpi_job_work_clear (NcmMPIJob *mpi_job);

GObject *ncm_mpi_job_run (NcmMPIJob *mpi_job, GObject *input);
NcmVector *ncm_mpi_job_run_vector (NcmMPIJob *mpi_job, NcmVector *input);

void ncm_mpi_job_init_all_slaves (NcmMPIJob *mpi_job, NcmSerialize *ser);
GPtrArray *ncm_mpi_job_run_vectors (NcmMPIJob *mpi_job, GPtrArray *input_array);

G_END_DECLS

#endif /* _NCM_MPI_JOB_H_ */
