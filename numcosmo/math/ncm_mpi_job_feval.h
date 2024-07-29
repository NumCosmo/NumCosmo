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

#define NCM_TYPE_MPI_JOB_FEVAL (ncm_mpi_job_feval_get_type ())

G_DECLARE_FINAL_TYPE (NcmMPIJobFEval, ncm_mpi_job_feval, NCM, MPI_JOB_FEVAL, NcmMPIJob)

NcmMPIJobFEval *ncm_mpi_job_feval_new (NcmFit * fit, NcmObjArray * func_oa);
NcmMPIJobFEval *ncm_mpi_job_feval_ref (NcmMPIJobFEval *mjfeval);

void ncm_mpi_job_feval_free (NcmMPIJobFEval *mjfeval);
void ncm_mpi_job_feval_clear (NcmMPIJobFEval **mjfeval);

G_END_DECLS

#endif /* _NCM_MPI_JOB_FEVAL_H_ */

