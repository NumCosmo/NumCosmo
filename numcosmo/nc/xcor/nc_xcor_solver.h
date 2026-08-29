/***************************************************************************
 *            nc_xcor_solver.h
 *
 *  Thu August 06 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

#ifndef _NC_XCOR_SOLVER_H_
#define _NC_XCOR_SOLVER_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/xcor/nc_xcor.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_SOLVER (nc_xcor_solver_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorSolver, nc_xcor_solver, NC, XCOR_SOLVER, GObject)

NcXcorSolver *nc_xcor_solver_new (void);
NcXcorSolver *nc_xcor_solver_ref (NcXcorSolver *solver);
void nc_xcor_solver_free (NcXcorSolver *solver);
void nc_xcor_solver_clear (NcXcorSolver **solver);

guint nc_xcor_solver_register_kernel (NcXcorSolver *solver, NcXcorKernel *xclk);
void nc_xcor_solver_request_cl (NcXcorSolver *solver, guint kernel_id_1, guint kernel_id_2, guint lmin, guint lmax);

guint nc_xcor_solver_get_n_kernels (NcXcorSolver *solver);
guint nc_xcor_solver_get_n_requests (NcXcorSolver *solver);
NcXcorKernel *nc_xcor_solver_peek_kernel (NcXcorSolver *solver, guint kernel_id);
void nc_xcor_solver_get_request (NcXcorSolver *solver, guint request_index, guint *kernel_id_1, guint *kernel_id_2, guint *lmin, guint *lmax);

void nc_xcor_solver_clear_requests (NcXcorSolver *solver);

void nc_xcor_solver_plan_blocks (NcXcorSolver *solver, guint default_block_size);
guint nc_xcor_solver_get_n_blocks (NcXcorSolver *solver);
void nc_xcor_solver_get_block (NcXcorSolver *solver, guint block_index, guint *lmin, guint *lmax);

void nc_xcor_solver_solve (NcXcorSolver *solver, NcXcor *xc, NcHICosmo *cosmo);
NcmVector *nc_xcor_solver_get_result (NcXcorSolver *solver, guint request_index);

void nc_xcor_solver_set_integrator (NcXcorSolver *solver, NcmSBesselIntegrator *sbi);
NcmSBesselIntegrator *nc_xcor_solver_peek_block_integrator (NcXcorSolver *solver, guint block_index);

G_END_DECLS

#endif /* _NC_XCOR_SOLVER_H_ */

