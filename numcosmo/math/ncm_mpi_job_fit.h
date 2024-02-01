/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_mpi_job_fit.h
 *
 *  Mon April 23 17:04:24 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_fit.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_MPI_JOB_FIT_H_
#define _NCM_MPI_JOB_FIT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mpi_job.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_MPI_JOB_FIT (ncm_mpi_job_fit_get_type ())

G_DECLARE_FINAL_TYPE (NcmMPIJobFit, ncm_mpi_job_fit, NCM, MPI_JOB_FIT, NcmMPIJob)

NcmMPIJobFit *ncm_mpi_job_fit_new (NcmFit * fit, NcmObjArray * func_oa);
NcmMPIJobFit *ncm_mpi_job_fit_ref (NcmMPIJobFit *mjfit);

void ncm_mpi_job_fit_free (NcmMPIJobFit *mjfit);
void ncm_mpi_job_fit_clear (NcmMPIJobFit **mjfit);

G_END_DECLS

#endif /* _NCM_MPI_JOB_FIT_H_ */

