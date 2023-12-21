/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_mpi_job_mcmc.h
 *
 *  Fri April 27 16:32:34 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_mcmc.h
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

#ifndef _NCM_MPI_JOB_MCMC_H_
#define _NCM_MPI_JOB_MCMC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mpi_job.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_MPI_JOB_MCMC (ncm_mpi_job_mcmc_get_type ())

G_DECLARE_FINAL_TYPE (NcmMPIJobMCMC, ncm_mpi_job_mcmc, NCM, MPI_JOB_MCMC, NcmMPIJob)

NcmMPIJobMCMC *ncm_mpi_job_mcmc_new (NcmFit * fit, NcmObjArray * func_oa);
NcmMPIJobMCMC *ncm_mpi_job_mcmc_ref (NcmMPIJobMCMC *mjmcmc);

void ncm_mpi_job_mcmc_free (NcmMPIJobMCMC *mjmcmc);
void ncm_mpi_job_mcmc_clear (NcmMPIJobMCMC **mjmcmc);

G_END_DECLS

#endif /* _NCM_MPI_JOB_MCMC_H_ */

