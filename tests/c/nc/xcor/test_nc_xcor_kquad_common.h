/***************************************************************************
 *            test_nc_xcor_kquad_common.h
 *
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

#ifndef _TEST_NC_XCOR_KQUAD_COMMON_H_
#define _TEST_NC_XCOR_KQUAD_COMMON_H_

#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

#define TEST_NC_XCOR_KQUAD_ZMAX 6.0
#define TEST_NC_XCOR_KQUAD_KMAX 0.5

/**
 * TestNcXcorKQuadEnv:
 * @cosmo: flat XCDM cosmology
 * @dist: distances to TEST_NC_XCOR_KQUAD_ZMAX
 * @ps: BBKS linear power spectrum capped at TEST_NC_XCOR_KQUAD_KMAX
 * @sbi: Levin integrator shared by the kernels, multipoles 0 to 8
 */
typedef struct _TestNcXcorKQuadEnv
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;
  NcmSBesselIntegrator *sbi;
} TestNcXcorKQuadEnv;

void test_nc_xcor_kquad_env_init (TestNcXcorKQuadEnv *env);
void test_nc_xcor_kquad_env_clear (TestNcXcorKQuadEnv *env);

NcXcorKernel *test_nc_xcor_kquad_tophat (TestNcXcorKQuadEnv *env, gdouble chi_lower, gdouble chi_upper, gint l_limber);
NcXcorKernel *test_nc_xcor_kquad_tophat_converged (TestNcXcorKQuadEnv *env, NcmSBesselIntegrator *sbi, gdouble chi_lower, gdouble chi_upper, gint l_limber, gdouble reltol, gdouble peak_epsilon);

typedef struct _TestNcXcorKQuadRef TestNcXcorKQuadRef;

TestNcXcorKQuadRef *test_nc_xcor_kquad_ref_new (void);
void test_nc_xcor_kquad_ref_free (TestNcXcorKQuadRef *ref);
gdouble test_nc_xcor_kquad_ref_eval (TestNcXcorKQuadRef *ref, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, gdouble RH, NcmVector *cl);

G_END_DECLS

#endif /* _TEST_NC_XCOR_KQUAD_COMMON_H_ */

