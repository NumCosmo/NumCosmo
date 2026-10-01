/***************************************************************************
 *            test_ncm_csq1d_bessel.h
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_csq1d_bessel.h
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

#ifndef _TEST_NCM_CSQ1D_BESSEL_H_
#define _TEST_NCM_CSQ1D_BESSEL_H_

#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

/*
 * The Bessel system of the CSQ1D tests: m = (s t)^(1 + 2 a) and nu = k, with s = -1 for
 * t < 0 (the adiabatic side, t -> -infinity) and s = 1 for t > 0. Its equation of motion
 * is phi'' + (1 + 2 a) phi' / t + k^2 phi = 0, solved by (s t)^(-a) times Bessel
 * functions of order a in k t.
 *
 * TestCSQ1DBessel implements every method it has a closed form for; the functions
 * below give those forms to the tests. TestCSQ1DBesselMin implements only xi, nu and
 * F1, leaving the rest to the NcmCSQ1D defaults.
 */
#define TEST_TYPE_CSQ1D_BESSEL (test_csq1d_bessel_get_type ())
G_DECLARE_FINAL_TYPE (TestCSQ1DBessel, test_csq1d_bessel, TEST, CSQ1D_BESSEL, NcmCSQ1D)

#define TEST_TYPE_CSQ1D_BESSEL_MIN (test_csq1d_bessel_min_get_type ())
G_DECLARE_FINAL_TYPE (TestCSQ1DBesselMin, test_csq1d_bessel_min, TEST, CSQ1D_BESSEL_MIN, NcmCSQ1D)

TestCSQ1DBessel *test_csq1d_bessel_new (gdouble a, gdouble k, gboolean adiab);
TestCSQ1DBesselMin *test_csq1d_bessel_min_new (gdouble a, gdouble k, gboolean adiab);

gdouble test_csq1d_bessel_m (gdouble a, gdouble k, gdouble s, gdouble t);
gdouble test_csq1d_bessel_xi (gdouble a, gdouble k, gdouble s, gdouble t);
gdouble test_csq1d_bessel_F1 (gdouble a, gdouble k, gdouble s, gdouble t);
gdouble test_csq1d_bessel_F2 (gdouble a, gdouble k, gdouble s, gdouble t);

G_END_DECLS

#endif /* _TEST_NCM_CSQ1D_BESSEL_H_ */

