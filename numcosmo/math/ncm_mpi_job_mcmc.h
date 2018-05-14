/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_mpi_job_mcmc.h
 *
 *  Fri April 27 16:32:34 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_mcmc.h
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

#ifndef _NCM_MPI_JOB_MCMC_H_
#define _NCM_MPI_JOB_MCMC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mpi_job.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_MPI_JOB_MCMC             (ncm_mpi_job_mcmc_get_type ())
#define NCM_MPI_JOB_MCMC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MPI_JOB_MCMC, NcmMPIJobMCMC))
#define NCM_MPI_JOB_MCMC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MPI_JOB_MCMC, NcmMPIJobMCMCClass))
#define NCM_IS_MPI_JOB_MCMC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MPI_JOB_MCMC))
#define NCM_IS_MPI_JOB_MCMC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MPI_JOB_MCMC))
#define NCM_MPI_JOB_MCMC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MPI_JOB_MCMC, NcmMPIJobMCMCClass))

typedef struct _NcmMPIJobMCMCClass NcmMPIJobMCMCClass;
typedef struct _NcmMPIJobMCMC NcmMPIJobMCMC;
typedef struct _NcmMPIJobMCMCPrivate NcmMPIJobMCMCPrivate;

struct _NcmMPIJobMCMCClass
{
	/*< private >*/
	NcmMPIJobClass parent_class;
};

struct _NcmMPIJobMCMC
{
	/*< private >*/
	NcmMPIJob parent_instance;
	NcmMPIJobMCMCPrivate *priv;
};

GType ncm_mpi_job_mcmc_get_type (void) G_GNUC_CONST;

NcmMPIJobMCMC *ncm_mpi_job_mcmc_new (NcmFit *fit, NcmObjArray *func_oa);
NcmMPIJobMCMC *ncm_mpi_job_mcmc_ref (NcmMPIJobMCMC *mjmcmc);

void ncm_mpi_job_mcmc_free (NcmMPIJobMCMC *mjmcmc);
void ncm_mpi_job_mcmc_clear (NcmMPIJobMCMC **mjmcmc);

G_END_DECLS

#endif /* _NCM_MPI_JOB_MCMC_H_ */

