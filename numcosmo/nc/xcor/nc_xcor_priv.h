/***************************************************************************
 *            nc_xcor_priv.h
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

#ifndef _NC_XCOR_PRIV_H_
#define _NC_XCOR_PRIV_H_

#include <numcosmo/nc/xcor/nc_xcor.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel.h>

G_BEGIN_DECLS

/*
 * Not public API: shared only between nc_xcor.c and nc_xcor_solver.c, so
 * NcXcorSolver can reuse nc_xcor.c's KERNEL_CUBATURE outer-integral
 * machinery with integrand(s) it built and cached itself, instead of
 * letting nc_xcor_compute() rebuild them once per pair.
 */
void _nc_xcor_kernel_integrate_block_cubature (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp);

/*
 * Assembles one pair's Cl block from a joint integrand covering several
 * kernels on a shared knot set, by exact 5-node Gauss-Legendre over the knot
 * panels. Component (kernel_id, il) lives at kernel_id * nell + il, with nell
 * taken from @vp. Lets NcXcorSolver build one joint integrand per ell block and
 * read every requested pair out of it, instead of one build per pair.
 */
void _nc_xcor_kernel_fixed_assemble (NcXcor *xc, NcXcorKernelIntegrand *xclki, guint kernel_id_1, guint kernel_id_2, NcmVector *vp);

/* Fails loudly when NcXcor:reltol asks for more than the kernel's closure carries. */
void _nc_xcor_check_kernel_tolerance (NcXcor *xc, NcXcorKernel *xclk);

G_END_DECLS

#endif /* _NC_XCOR_PRIV_H_ */

