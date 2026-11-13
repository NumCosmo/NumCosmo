/***************************************************************************
 *            test_nc_galaxy_shape_marginal_arb_cases.h
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_shape_marginal_arb_cases.h
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

/*
 * Values of ln P (eps_obs | g) certified by the Arb reference
 * tests/tools/nc_galaxy_shape_marginal_arb.c, checked by the NcGalaxyShapeFactorSpectral
 * tests. The table is printed by tests/tools/make_shape_marginal_arb_cases.py.
 */

#ifndef _TEST_NC_GALAXY_SHAPE_MARGINAL_ARB_CASES_H
#define _TEST_NC_GALAXY_SHAPE_MARGINAL_ARB_CASES_H

#include <glib.h>
#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

typedef enum _TestShapeArbPop
{
  TEST_SHAPE_ARB_POP_GAUSS,
  TEST_SHAPE_ARB_POP_BETA,
} TestShapeArbPop;

/*
 * One certified value: the shear convention, the population (Gauss: sigma = p1; Beta:
 * alpha = p1, beta = p2), std_noise, g, eps_obs and ln P.
 */
typedef struct _TestShapeArbCell
{
  NcGalaxyWLObsEllipConv conv;
  TestShapeArbPop pop;
  gdouble p1;
  gdouble p2;
  gdouble std_noise;
  gdouble g1;
  gdouble g2;
  gdouble e1;
  gdouble e2;
  gdouble lnP;
} TestShapeArbCell;

static const TestShapeArbCell test_shape_arb_cells[] = {
  /* BEGIN make_shape_marginal_arb_cases.py */
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_GAUSS, 0.3, 0.0, 0.1, 0.3, 0.0, 0.0, 0.5, -1.402431117754455},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_GAUSS, 0.3, 0.0, 0.01, 0.6, 0.0, 0.5, 0.0, -2.536716442503229},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_GAUSS, 0.05, 0.0, 0.005, 0.1, 0.0, 0.05, 0.0, -0.31255850401115876},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_GAUSS, 0.05, 0.0, 0.03, 0.6, 0.0, -0.9, 0.0, -182.8521316909928},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_GAUSS, 0.2, 0.0, 0.05, 0.2, -0.4, 0.3, 0.25, -3.978049030486719},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_GAUSS, 0.3, 0.0, 1.0, 3.0, 0.0, 2.0, 0.0, -2.872477807137722},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_GAUSS, 0.3, 0.0, 0.1, 1.5, 0.0, 1.2, 0.0, -4.056650475069435},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.0, 1.0, 0.03, 0.6, 0.0, 0.0, 0.5, -4.003703139471936},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 4.0, 1.0, 0.1, 0.6, 0.0, 0.9, 0.0, 0.7997452547544199},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.0, 1.0, 0.01, 0.3, 0.0, 0.0, 1.2, -206.2067527885707},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 4.0, 1.0, 0.3, 1.5, 0.0, -0.5, 0.0, -4.309400768872423},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.57, 1.65, 0.1, 0.3, 0.0, 0.5, 0.0, 0.5334663856834809},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.92, 1.98, 0.03, 0.6, 0.0, 0.0, 0.5, -2.5176649373938598},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.57, 1.65, 0.05, -0.1, 0.25, 0.2, -0.3, -1.532207192020498},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_GAUSS, 0.05, 0.0, 0.1, 0.005111126056639749, 0.0, 0.0, 0.9, -29.85714990728496},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.57, 1.65, 0.1, 0.01, 0.0, 0.0, 0.02, 0.2155439914200594},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 4.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -2.1641169160898808},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_GAUSS, 0.05, 0.0, 0.03, 0.3, -0.219359538352843, 0.99, 0.0, -79.94280295160708},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_GAUSS, 0.05, 0.0, 0.005, 0.3227416423252818, 0.46019288007369547, 0.99, 0.0, -189.91837555163738},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_GAUSS, 0.4, 0.0, 0.01, -0.8659144781362805, 0.24534081713806422, -0.9340398200927157, -0.3281304839254761, -5.125302056071098},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_GAUSS, 0.2, 0.0, 0.03, -0.368509388057976, 0.4734985014898527, -0.8999998321433729, 0.0005496743587696258, -8.44179742358457},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.5, 3.0, 0.1, 0.0020353798786394857, -0.009790670495406814, -0.5855653805208657, -0.7982563404918585, -5.007649495155437},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 4.0, 1.0, 0.03, -0.5386369764538615, 0.26432973271397636, 0.800491458947086, 0.4113555933164953, -2.7305199739899457},
  /* END make_shape_marginal_arb_cases.py */
};

G_END_DECLS

#endif /* _TEST_NC_GALAXY_SHAPE_MARGINAL_ARB_CASES_H */

