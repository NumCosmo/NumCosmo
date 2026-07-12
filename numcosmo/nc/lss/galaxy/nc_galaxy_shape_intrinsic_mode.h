/***************************************************************************
 *            nc_galaxy_shape_intrinsic_mode.h
 *
 *  Fri Jul 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_intrinsic_mode.h
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _NC_GALAXY_SHAPE_INTRINSIC_MODE_H_
#define _NC_GALAXY_SHAPE_INTRINSIC_MODE_H_

#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_pop.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#ifndef NUMCOSMO_GIR_SCAN

/*
 * NcGalaxyShapeIntrinsicMode:
 *
 * The joint mode of the intrinsic-ellipticity integrand,
 *   P_pop(chi_I) N_2(eps_obs - f_g(chi_I); std_noise^2),
 * over chi_I in the unit disc, plus the Hessian of its log at that point
 * (Cartesian chi_I_1, chi_I_2), for building a local Gaussian/Laplace
 * description of the peak. Not a GObject: a plain, stateless result struct
 * used by #NcGalaxyShapeFactorLaplace (the analytic approximation itself).
 */
typedef struct _NcGalaxyShapeIntrinsicMode
{
  gdouble rho;
  gdouble theta;
  gdouble hxx;
  gdouble hxy;
  gdouble hyy;
  gdouble ln_peak;
} NcGalaxyShapeIntrinsicMode;

void nc_galaxy_shape_intrinsic_mode_find (complex double (*apply_shear)(complex double, complex double),
                                          complex double (*apply_shear_inv) (complex double, complex double),
                                          NcGalaxyShapePop * pop, NcGalaxyShapePopData * pop_data,
                                          complex double g, complex double eps_obs, gdouble std_noise,
                                          NcGalaxyShapeIntrinsicMode * mode);

void nc_galaxy_shape_intrinsic_mode_find_trace_det (NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                                    complex double g, complex double eps_obs, gdouble std_noise,
                                                    NcGalaxyShapeIntrinsicMode *mode);

void nc_galaxy_shape_intrinsic_mode_find_trace (NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                                complex double g, complex double eps_obs, gdouble std_noise,
                                                NcGalaxyShapeIntrinsicMode *mode);

gdouble nc_galaxy_shape_intrinsic_mode_laplace (const NcGalaxyShapeIntrinsicMode *mode);

#endif /* NUMCOSMO_GIR_SCAN */

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_INTRINSIC_MODE_H_ */

