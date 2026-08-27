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
 * @NC_XCOR_METHOD_KERNEL_FIXED: Use fixed-knot Gauss-Legendre over kernel variables
 *
 * Methods to compute integrals.
 *
 * %NC_XCOR_METHOD_KERNEL_FIXED needs no tolerance and cannot fail to converge.
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
typedef enum _NcXcorMethod
{
  NC_XCOR_METHOD_LIMBER_Z_GSL = 0,
  NC_XCOR_METHOD_LIMBER_Z_CUBATURE,
  NC_XCOR_METHOD_KERNEL_GSL,
  NC_XCOR_METHOD_KERNEL_CUBATURE,
  NC_XCOR_METHOD_KERNEL_FIXED,
} NcXcorMethod;

#define NC_XCOR_PRECISION (1.0e-6)

GType nc_xcor_kinetic_get_type (void) G_GNUC_CONST;

NcXcor *nc_xcor_new (NcDistance *dist, NcmPowspec *ps, NcXcorMethod meth);
NcXcor *nc_xcor_ref (NcXcor *xc);
void nc_xcor_free (NcXcor *xc);
void nc_xcor_clear (NcXcor **xc);

NcXcorMethod nc_xcor_get_meth (NcXcor *xc);

void nc_xcor_set_reltol (NcXcor *xc, const gdouble reltol);
gdouble nc_xcor_get_reltol (NcXcor *xc);

void nc_xcor_set_ell_batch_size (NcXcor *xc, const guint ell_batch_size);
guint nc_xcor_get_ell_batch_size (NcXcor *xc);

void nc_xcor_prepare (NcXcor *xc, NcHICosmo *cosmo);

void nc_xcor_compute (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp);

G_END_DECLS

#endif /* _NC_XCOR_H_ */

