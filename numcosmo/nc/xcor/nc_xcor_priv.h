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
 * Same for %NC_XCOR_METHOD_KERNEL_FIXED: exact 5-node Gauss-Legendre over the
 * common refinement of the two integrands' own knot sets. Takes the same
 * arguments as the cubature version above and is interchangeable with it, so
 * NcXcorSolver drives both methods through one cached-integrand path.
 */
void _nc_xcor_kernel_integrate_block_fixed (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp);

/* Fails loudly when NcXcor:reltol asks for more than the kernel's closure carries. */
void _nc_xcor_check_kernel_tolerance (NcXcor *xc, NcXcorKernel *xclk);

/*
 * TRUE when both kernels run in the Limber tier over a block starting at @lmin
 * and their redshift supports do not overlap, in which case their Cl vanishes
 * identically and the kernel-space methods must not integrate it.
 *
 * A Limber kernel is supported only where xi = (l + 1/2) / k lies inside its
 * own radial range, so two disjoint bins have disjoint support in k and their
 * product is zero -- the same statement as the Limber-z tier's overlap test.
 * The non-Limber tier is the opposite case: there the two kernels couple only
 * through the outer k integral and disjoint bins do correlate, which is why
 * this must be asked per tier and not once per method.
 */
gboolean _nc_xcor_kernels_limber_disjoint (NcXcorKernel *xclk1, NcXcorKernel *xclk2, gboolean isauto, guint lmin);

G_END_DECLS

#endif /* _NC_XCOR_PRIV_H_ */

