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
#ifdef NUMCOSMO_HAVE_MPI
#  ifndef NUMCOSMO_GIR_SCAN
#    include <mpi.h>
typedef MPI_Datatype NcmMPIDatatype;
#  else
typedef gpointer NcmMPIDatatype;
#  endif /* NUMCOSMO_GIR_SCAN */
#else
typedef gpointer NcmMPIDatatype;
#endif /* NUMCOSMO_HAVE_MPI */

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
	void (*work_init) (NcmMPIJob *mpi_job);
	void (*work_clear) (NcmMPIJob *mpi_job);
	NcmMPIDatatype (*input_datatype) (NcmMPIJob *mpi_job, gint *len, gint *size);
	NcmMPIDatatype (*return_datatype) (NcmMPIJob *mpi_job, gint *len, gint *size);
	gpointer (*create_input) (NcmMPIJob *mpi_job);
	gpointer (*create_return) (NcmMPIJob *mpi_job);
	void (*destroy_input) (NcmMPIJob *mpi_job, gpointer input);
	void (*destroy_return) (NcmMPIJob *mpi_job, gpointer ret);
	gpointer (*get_input_buffer) (NcmMPIJob *mpi_job, gpointer input);
	gpointer (*get_return_buffer) (NcmMPIJob *mpi_job, gpointer ret);
	void (*destroy_input_buffer) (NcmMPIJob *mpi_job, gpointer input, gpointer buf);
	void (*destroy_return_buffer) (NcmMPIJob *mpi_job, gpointer ret, gpointer buf);
	gpointer (*pack_input) (NcmMPIJob *mpi_job, gpointer input);
	gpointer (*pack_return) (NcmMPIJob *mpi_job, gpointer ret);
	void (*unpack_input) (NcmMPIJob *mpi_job, gpointer buf, gpointer input);
	void (*unpack_return) (NcmMPIJob *mpi_job, gpointer buf, gpointer ret);
	void (*run) (NcmMPIJob *mpi_job, gpointer input, gpointer ret);
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
	gint initialized;
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
 * @NCM_MPI_CTRL_SLAVE_KILL: FIXME
 * @NCM_MPI_CTRL_SLAVE_WORK: FIXME
 * 
 * FIXME
 * 
 */ 
typedef enum _NcmMPIJobCtrlMsg 
{
	NCM_MPI_CTRL_SLAVE_INIT = 0,
	NCM_MPI_CTRL_SLAVE_FREE,
	NCM_MPI_CTRL_SLAVE_KILL,
	NCM_MPI_CTRL_SLAVE_WORK,
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

void ncm_mpi_job_work_init (NcmMPIJob *mpi_job);
void ncm_mpi_job_work_clear (NcmMPIJob *mpi_job);

NcmMPIDatatype ncm_mpi_job_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size);
NcmMPIDatatype ncm_mpi_job_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size);

gpointer ncm_mpi_job_create_input (NcmMPIJob *mpi_job);
gpointer ncm_mpi_job_create_return (NcmMPIJob *mpi_job);

void ncm_mpi_job_destroy_input (NcmMPIJob *mpi_job, gpointer input);
void ncm_mpi_job_destroy_return (NcmMPIJob *mpi_job, gpointer ret);

gpointer ncm_mpi_job_get_input_buffer (NcmMPIJob *mpi_job, gpointer input);
gpointer ncm_mpi_job_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret);

void ncm_mpi_job_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf);
void ncm_mpi_job_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf);

gpointer ncm_mpi_job_pack_input (NcmMPIJob *mpi_job, gpointer input);
gpointer ncm_mpi_job_pack_return (NcmMPIJob *mpi_job, gpointer ret);

void ncm_mpi_job_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input);
void ncm_mpi_job_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret);

void ncm_mpi_job_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret);

void ncm_mpi_job_init_all_slaves (NcmMPIJob *mpi_job, NcmSerialize *ser);
void ncm_mpi_job_run_array (NcmMPIJob *mpi_job, GPtrArray *input_array, GPtrArray *ret_array);

void ncm_mpi_job_free_all_slaves (NcmMPIJob *mpi_job);

/*#define NCM_MPI_DEBUG 1*/
#ifdef NCM_MPI_DEBUG
#define NCM_MPI_JOB_DEBUG_PRINT(fmt, ...)  printf (fmt, ##__VA_ARGS__)
#else
#define NCM_MPI_JOB_DEBUG_PRINT(fmt, ...)
#endif

G_END_DECLS

#endif /* _NCM_MPI_JOB_H_ */
