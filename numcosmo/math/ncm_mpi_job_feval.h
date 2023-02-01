/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_mpi_job_feval.h
 *
 *  Fri February 21 11:16:10 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_feval.h
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_MPI_JOB_FEVAL_H_
#define _NCM_MPI_JOB_FEVAL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mpi_job.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_MPI_JOB_FEVAL             (ncm_mpi_job_feval_get_type ())
#define NCM_MPI_JOB_FEVAL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MPI_JOB_FEVAL, NcmMPIJobFEval))
#define NCM_MPI_JOB_FEVAL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MPI_JOB_FEVAL, NcmMPIJobFEvalClass))
#define NCM_IS_MPI_JOB_FEVAL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MPI_JOB_FEVAL))
#define NCM_IS_MPI_JOB_FEVAL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MPI_JOB_FEVAL))
#define NCM_MPI_JOB_FEVAL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MPI_JOB_FEVAL, NcmMPIJobFEvalClass))

typedef struct _NcmMPIJobFEvalClass NcmMPIJobFEvalClass;
typedef struct _NcmMPIJobFEval NcmMPIJobFEval;
typedef struct _NcmMPIJobFEvalPrivate NcmMPIJobFEvalPrivate;

struct _NcmMPIJobFEvalClass
{
  /*< private >*/
  NcmMPIJobClass parent_class;
};

struct _NcmMPIJobFEval
{
  /*< private >*/
  NcmMPIJob parent_instance;
  NcmMPIJobFEvalPrivate *priv;
};

GType ncm_mpi_job_feval_get_type (void) G_GNUC_CONST;

NcmMPIJobFEval *ncm_mpi_job_feval_new (NcmFit *fit, NcmObjArray *func_oa);
NcmMPIJobFEval *ncm_mpi_job_feval_ref (NcmMPIJobFEval *mjfeval);

void ncm_mpi_job_feval_free (NcmMPIJobFEval *mjfeval);
void ncm_mpi_job_feval_clear (NcmMPIJobFEval **mjfeval);

G_END_DECLS

#endif /* _NCM_MPI_JOB_FEVAL_H_ */

