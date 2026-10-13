/***************************************************************************
 *            test_nc_xcor_closure_check.h
 *
 *  Sun Oct 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#ifndef _TEST_NC_XCOR_CLOSURE_CHECK_H_
#define _TEST_NC_XCOR_CLOSURE_CHECK_H_

#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

/* The worst deviation of a closure from the direct values, over its bound, where it
 * occurred, and over how many points */
typedef struct _TestNcXcorClosureCheck
{
  gdouble worst;
  gint at_ell;
  gdouble at_k;
  guint n_points;
} TestNcXcorClosureCheck;

gdouble test_nc_xcor_closure_direct (NcXcorKernel *xclk, NcHICosmo *cosmo, NcmSBesselIntegrator *sbi, gint ell, gdouble k);
TestNcXcorClosureCheck test_nc_xcor_closure_check (NcXcorKernel *xclk, NcHICosmo *cosmo, NcXcorKernelIntegrand *closure, gint lmin);

G_END_DECLS

#endif /* _TEST_NC_XCOR_CLOSURE_CHECK_H_ */

