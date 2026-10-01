/***************************************************************************
 *            test_ncm_mpi_job_shape.h
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mpi_job_shape.h
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


#ifndef _TEST_NCM_MPI_JOB_SHAPE_H_
#define _TEST_NCM_MPI_JOB_SHAPE_H_

#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

/*
 * TestMPIJobShape: a job whose input and return messages have different lengths and
 * travel in the default pooled buffers of NcmMPIJob, which no library job uses. The
 * return is ret[k] = (k + 1) sum_j (j + 1) input[j]. With float-input the input travels
 * as MPI_FLOAT, so the return message is longer in bytes than return-len elements of
 * the input datatype.
 */
#define TEST_TYPE_MPI_JOB_SHAPE (test_mpi_job_shape_get_type ())

G_DECLARE_FINAL_TYPE (TestMPIJobShape, test_mpi_job_shape, TEST, MPI_JOB_SHAPE, NcmMPIJob)

G_END_DECLS

#endif /* _TEST_NCM_MPI_JOB_SHAPE_H_ */

