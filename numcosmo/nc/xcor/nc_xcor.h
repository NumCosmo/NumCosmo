/***************************************************************************
 *            nc_xcor.h
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_XCOR_H_
#define _NC_XCOR_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_c.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/nc/background/nc_distance.h>
#include <numcosmo/nc/background/nc_hicosmo.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel.h>
#include <numcosmo/ncm/powspec/ncm_powspec.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR (nc_xcor_get_type ())

G_DECLARE_FINAL_TYPE (NcXcor, nc_xcor, NC, XCOR, GObject)

/**
 * NcXcorMethod:
 * @NC_XCOR_METHOD_LIMBER_Z_GSL: Use GSL numerical integration
 * @NC_XCOR_METHOD_LIMBER_Z_CUBATURE: Use cubature numerical integration
 * @NC_XCOR_METHOD_KERNEL_GSL: Use GSL numerical integration over kernel variables
 * @NC_XCOR_METHOD_KERNEL_CUBATURE: Use cubature numerical integration over kernel variables
 * @NC_XCOR_METHOD_KERNEL_EXACT: Integrate the kernel closures exactly, on the
 *   union of their own knots
 * @NC_XCOR_METHOD_KERNEL_GSL_BLOCK: As %NC_XCOR_METHOD_KERNEL_GSL, but on the
 *   block closure the other kernel-space methods use
 *
 * Methods to compute integrals.
 *
 * %NC_XCOR_METHOD_KERNEL_GSL and %NC_XCOR_METHOD_KERNEL_GSL_BLOCK run the same
 * quadrature -- QUADPACK's qagp, broken on the closures' own knots -- over two
 * different integrands, and they are separate methods because the integrands
 * differ. The first fits a closure to one multipole at a time; the second
 * shares one closure across a whole #NcXcor:ell-batch-size block, fitted to an
 * $L^2$ norm over the block, which is what %NC_XCOR_METHOD_KERNEL_CUBATURE and
 * %NC_XCOR_METHOD_KERNEL_EXACT integrate. On a far-separated pair of top-hat
 * bins the two closures give $-1.05 \times 10^{-8}$ against $3.51 \times
 * 10^{-8}$ at the library's default tolerances, and agree from $10^{-6}$ down.
 * Use the block form to compare quadratures on one integrand, and the
 * per-multipole form as the independent check that does not share a closure
 * with anything.
 *
 * %NC_XCOR_METHOD_KERNEL_EXACT needs no tolerance and cannot fail to converge.
 * The name is meant literally rather than as "fixed-order": on the closures it
 * is handed, the quadrature carries no error at all, and refining its panels
 * changes nothing beyond rounding.
 * It uses the same per-kernel closures as %NC_XCOR_METHOD_KERNEL_CUBATURE and
 * differs only in the outer quadrature: each kernel's $W(k)$ is a cubic spline,
 * so on the common refinement of a pair's two knot sets the outer integrand
 * $k^2 W_i W_j$ is a degree-8 polynomial on every panel, and a 5-node
 * Gauss-Legendre rule integrates it exactly. The adaptive alternatives target a
 * tolerance the integrand may not be able to support, and abort when they
 * cannot reach it.
 *
 * Measured over 28 pairs of 7 top-hat bins it is also slightly faster than
 * %NC_XCOR_METHOD_KERNEL_CUBATURE (1.05x at $\ell = 0$, 1.18x over
 * $\ell = 0\dots26$), since the exact rule replaces adaptive refinement on
 * splines that have already been built.
 *
 */
typedef enum _NcXcorMethod /*< prefix=NC_XCOR_METHOD >*/
{
  NC_XCOR_METHOD_LIMBER_Z_GSL = 0,
  NC_XCOR_METHOD_LIMBER_Z_CUBATURE,
  NC_XCOR_METHOD_KERNEL_GSL,
  NC_XCOR_METHOD_KERNEL_CUBATURE,
  NC_XCOR_METHOD_KERNEL_EXACT,
  NC_XCOR_METHOD_KERNEL_GSL_BLOCK,
} NcXcorMethod;

#define NC_XCOR_PRECISION (1.0e-6)

GType nc_xcor_kinetic_get_type (void) G_GNUC_CONST;

NcXcor *nc_xcor_new (NcDistance *dist, NcmPowspec *ps, NcXcorMethod meth);
NcXcor *nc_xcor_ref (NcXcor *xc);
void nc_xcor_free (NcXcor *xc);
void nc_xcor_clear (NcXcor **xc);

NcXcorMethod nc_xcor_get_meth (NcXcor *xc);

void nc_xcor_set_closure_type (NcXcor *xc, NcXcorKernelClosure closure_type);
NcXcorKernelClosure nc_xcor_get_closure_type (NcXcor *xc);

void nc_xcor_set_reltol (NcXcor *xc, const gdouble reltol);
gdouble nc_xcor_get_reltol (NcXcor *xc);

void nc_xcor_set_ell_batch_size (NcXcor *xc, const guint ell_batch_size);
guint nc_xcor_get_ell_batch_size (NcXcor *xc);

void nc_xcor_prepare (NcXcor *xc, NcHICosmo *cosmo);

void nc_xcor_compute (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp);
void nc_xcor_compute_full (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp, NcmVector *vp_err);

void nc_xcor_integrate_block (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcXcorMethod meth, NcmVector *vp, NcmVector *vp_err);

const gchar *nc_xcor_method_get_name (NcXcorMethod meth);
gboolean nc_xcor_method_has_error_estimate (NcXcorMethod meth);
gboolean nc_xcor_method_is_kernel_space (NcXcorMethod meth);

G_END_DECLS

#endif /* _NC_XCOR_H_ */

