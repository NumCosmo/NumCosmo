/***************************************************************************
 *            ncm_mpi_slave.h
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mpi_slave.h
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

/*
 * The worker side of the NcmMPIJob protocol. Internal: ncm_cfg_init() runs it on every
 * rank but the master, and tests may call it directly.
 */

#ifndef _NCM_MPI_SLAVE_H_
#define _NCM_MPI_SLAVE_H_

#include <glib.h>

G_BEGIN_DECLS

void ncm_mpi_slave_run (void);

gboolean ncm_mpi_slave_serve_job (void);

G_END_DECLS

#endif /* _NCM_MPI_SLAVE_H_ */

