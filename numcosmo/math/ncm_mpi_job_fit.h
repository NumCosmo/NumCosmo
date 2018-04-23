/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_mpi_job_fit.h
 *
 *  Mon April 23 17:04:24 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_fit.h
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

#ifndef _NCM_MPI_JOB_FIT_H_
#define _NCM_MPI_JOB_FIT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mpi_job.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_MPI_JOB_FIT             (ncm_mpi_job_fit_get_type ())
#define NCM_MPI_JOB_FIT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MPI_JOB_FIT, NcmMPIJobFit))
#define NCM_MPI_JOB_FIT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MPI_JOB_FIT, NcmMPIJobFitClass))
#define NCM_IS_MPI_JOB_FIT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MPI_JOB_FIT))
#define NCM_IS_MPI_JOB_FIT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MPI_JOB_FIT))
#define NCM_MPI_JOB_FIT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MPI_JOB_FIT, NcmMPIJobFitClass))

typedef struct _NcmMPIJobFitClass NcmMPIJobFitClass;
typedef struct _NcmMPIJobFit NcmMPIJobFit;
typedef struct _NcmMPIJobFitPrivate NcmMPIJobFitPrivate;

struct _NcmMPIJobFitClass
{
	/*< private >*/
	NcmMPIJobClass parent_class;
};

struct _NcmMPIJobFit
{
	/*< private >*/
	NcmMPIJob parent_instance;
	NcmMPIJobFitPrivate *priv;
};

/**
 * NcmMPIJobFitType:
 * @NCM_MPI_JOB_FIT_TYPE_M2LNL_VAL: computes $-2\ln(L)$.
 * 
 * #NcmMPIJobFit job types.
 * 
 */ 
typedef enum _NcmMPIJobFitType
{
	NCM_MPI_JOB_FIT_TYPE_M2LNL_VAL = 0,
	/* < private > */
	NCM_MPI_JOB_FIT_TYPE_LEN, /*< skip >*/
} NcmMPIJobFitType;

GType ncm_mpi_job_fit_get_type (void) G_GNUC_CONST;

NcmMPIJobFit *ncm_mpi_job_fit_new (NcmFit *fit, NcmMPIJobFitType job_type);
NcmMPIJobFit *ncm_mpi_job_fit_ref (NcmMPIJobFit *mjfit);

void ncm_mpi_job_fit_free (NcmMPIJobFit *mjfit);
void ncm_mpi_job_fit_clear (NcmMPIJobFit **mjfit);

G_END_DECLS

#endif /* _NCM_MPI_JOB_FIT_H_ */

