/***************************************************************************
 *            nc_xcor.h
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/xcor/nc_xcor_limber_kernel.h>
#include <numcosmo/math/ncm_powspec.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR (nc_xcor_get_type ())

G_DECLARE_FINAL_TYPE (NcXcor, nc_xcor, NC, XCOR, GObject)

/**
 * NcXcorLimberMethod:
 * @NC_XCOR_LIMBER_METHOD_GSL: Use GSL numerical integration
 * @NC_XCOR_LIMBER_METHOD_CUBATURE: Use cubature numerical integration
 *
 * Methods to compute Limber-approximated integrals.
 *
 */
typedef enum _NcXcorLimberMethod
{
  NC_XCOR_LIMBER_METHOD_GSL = 0,
  NC_XCOR_LIMBER_METHOD_CUBATURE,
} NcXcorLimberMethod;

#define NC_XCOR_PRECISION (1.0e-5)

GType nc_xcor_kinetic_get_type (void) G_GNUC_CONST;

NcXcor *nc_xcor_new (NcDistance *dist, NcmPowspec *ps, NcXcorLimberMethod meth);
NcXcor *nc_xcor_ref (NcXcor *xc);
void nc_xcor_free (NcXcor *xc);
void nc_xcor_clear (NcXcor **xc);

void nc_xcor_prepare (NcXcor *xc, NcHICosmo *cosmo);

void nc_xcor_limber (NcXcor *xc, NcXcorLimberKernel *xclk1, NcXcorLimberKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp);

G_END_DECLS

#endif /* _NC_XCOR_H_ */

