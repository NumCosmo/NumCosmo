/***************************************************************************
 *            test_ncm_fit_esmcmc_parity.h
 *
 *  Tue September 15 19:10:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#ifndef _TEST_NCM_FIT_ESMCMC_PARITY_H_
#define _TEST_NCM_FIT_ESMCMC_PARITY_H_

#include <glib.h>
#include <numcosmo/numcosmo.h>

/*
 * APES with the production configuration (VKDE, Gaussian kernel, out-of-sample CV with the
 * kernel chosen by the same objective, centre shrinkage, uniform weights, 12 points per
 * dimension) on a fixed-seed
 * 3-d MVND target: one initial ensemble plus one iteration, so that the arms can be
 * compared before the chain amplifies rounding differences. Shared by the plain and the
 * MPI test binaries so that every parallel mode is measured against the same serial run.
 */
NcmMSetCatalog *test_ncm_fit_esmcmc_parity_apes_catalog (gboolean use_threads, gboolean use_mpi);

/*
 * Compares the two catalogs row by row with ncm_assert_cmpdouble_e (); reltol = 0 and
 * abstol = 0 demand bit identity.
 */
void test_ncm_fit_esmcmc_parity_compare (NcmMSetCatalog *ref, NcmMSetCatalog *other, const gdouble reltol, const gdouble abstol);

#endif /* _TEST_NCM_FIT_ESMCMC_PARITY_H_ */

