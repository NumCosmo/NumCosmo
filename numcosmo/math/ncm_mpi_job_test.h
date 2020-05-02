/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_mpi_job_test.h
 *
 *  Sun April 22 14:47:43 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_test.h
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

#ifndef _NCM_MPI_JOB_TEST_H_
#define _NCM_MPI_JOB_TEST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mpi_job.h>
#include <numcosmo/math/ncm_rng.h>

G_BEGIN_DECLS

#define NCM_TYPE_MPI_JOB_TEST             (ncm_mpi_job_test_get_type ())
#define NCM_MPI_JOB_TEST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MPI_JOB_TEST, NcmMPIJobTest))
#define NCM_MPI_JOB_TEST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MPI_JOB_TEST, NcmMPIJobTestClass))
#define NCM_IS_MPI_JOB_TEST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MPI_JOB_TEST))
#define NCM_IS_MPI_JOB_TEST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MPI_JOB_TEST))
#define NCM_MPI_JOB_TEST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MPI_JOB_TEST, NcmMPIJobTestClass))

typedef struct _NcmMPIJobTestClass NcmMPIJobTestClass;
typedef struct _NcmMPIJobTest NcmMPIJobTest;
typedef struct _NcmMPIJobTestPrivate NcmMPIJobTestPrivate;

struct _NcmMPIJobTestClass
{
  /*< private >*/
  NcmMPIJobClass parent_class;
};

struct _NcmMPIJobTest
{
  /*< private >*/
  NcmMPIJob parent_instance;
  NcmMPIJobTestPrivate *priv;
};

GType ncm_mpi_job_test_get_type (void) G_GNUC_CONST;

NcmMPIJobTest *ncm_mpi_job_test_new (void);
NcmMPIJobTest *ncm_mpi_job_test_ref (NcmMPIJobTest *mjt);

void ncm_mpi_job_test_free (NcmMPIJobTest *mjt);
void ncm_mpi_job_test_clear (NcmMPIJobTest **mjt);

void ncm_mpi_job_test_set_rand_vector (NcmMPIJobTest *mjt, const guint len, NcmRNG *rng);

G_END_DECLS

#endif /* _NCM_MPI_JOB_TEST_H_ */

