/***************************************************************************
 *            nc_galaxy_shape_factor_reach_private.h
 *
 *  Mon Nov 16 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_reach_private.h
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

/*
 * The reachable range of the shear of a galaxy over the box of a fit, shared by the shape
 * factors that restrict their per-galaxy work to it (#NcGalaxyShapeFactorMomentsTilt and
 * #NcGalaxyShapeFactorSpectral). Not public API.
 */

#ifndef _NC_GALAXY_SHAPE_FACTOR_REACH_PRIVATE_H_
#define _NC_GALAXY_SHAPE_FACTOR_REACH_PRIVATE_H_

#ifndef NUMCOSMO_GIR_SCAN

#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/nc/background/nc_hicosmo.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>
#include <numcosmo/nc/lss/halo/nc_halo_density_profile.h>
#include <numcosmo/nc/lss/wl/nc_wl_surface_mass_density.h>
#include <numcosmo/ncm/algebra/ncm_quaternion.h>

G_BEGIN_DECLS

/*
 * The bound of the reduced shear a galaxy can reach over the box of a fit. prepare() takes
 * the box from the models of the fit: the halo center over the bounds of ra and dec (or
 * their values when fixed), the mass up to the upper bound of log10MDelta, and, with
 * sample_c, a free cDelta over its bounds. The bound depends on the box only, never on the
 * current point of the fit, so a value computed once holds for the whole chain.
 *
 * holds: the box covers every free parameter of the fit. Only log10MDelta, ra, dec, the
 * parameters of the shape population and, with sample_c, cDelta may be free; any other free
 * parameter (the cosmology, the redshift of the halo, the other parameters of the profile)
 * moves the shear in ways the bound does not follow. models_ready: the models the bound
 * needs are in the fit. gen: changes whenever the box or the fixed models change, so a
 * galaxy recomputes its bound only then (NcGalaxyShapeFactorReachMemo).
 *
 * The profile copy dp_box is shared by the galaxies, which may be bounded in parallel, and
 * a profile updates state lazily on use: it is used under lock.
 */
typedef struct _NcGalaxyShapeFactorReach
{
  gboolean sample_c;
  gboolean models_ready;
  gboolean holds;
  guint gen;
  guint64 key;
  NcHICosmo *cosmo;
  NcWLSurfaceMassDensity *smd;
  NcHaloDensityProfile *dp_box;
  gdouble z_cl;
  gdouble prefactor; /* projected radius per sin of the angle, at z_cl */
  NcmTriVec centre;  /* unit vector of the center of the box of the halo center */
  gdouble delta;     /* largest projected distance from it to a corner of the box */
  gboolean has_c;
  guint c_idx;
  gdouble c_lo;
  gdouble c_hi;
  GMutex lock;
} NcGalaxyShapeFactorReach;

/* The bound of one galaxy and what it was computed for */
typedef struct _NcGalaxyShapeFactorReachMemo
{
  guint gen;
  gdouble z_max;
  gdouble ra;
  gdouble dec;
  gdouble m;
  gdouble c1;
  gdouble c2;
  gdouble value;
} NcGalaxyShapeFactorReachMemo;

/* Concentrations at which a free cDelta is sampled, and the factor that covers the shear between them */
#define NC_GALAXY_SHAPE_FACTOR_REACH_NC (9)
#define NC_GALAXY_SHAPE_FACTOR_REACH_C_MARGIN (1.05)

void _nc_galaxy_shape_factor_reach_init (NcGalaxyShapeFactorReach *reach, const gboolean sample_c);
void _nc_galaxy_shape_factor_reach_clear (NcGalaxyShapeFactorReach *reach);
void _nc_galaxy_shape_factor_reach_prepare (NcGalaxyShapeFactorReach *reach, NcmMSet *mset);
gdouble _nc_galaxy_shape_factor_reach_eval (NcGalaxyShapeFactorReach *reach, NcGalaxyShapeFactorReachMemo *memo, NcGalaxyShapeFactorData *data, const gdouble z_max);

/* The line of the shears of a galaxy for another profile (nc_galaxy_shape_factor.c), such as
 * dp_box: the calibration of its redshift nodes at the upper bound of the mass */
void nc_galaxy_shape_factor_data_prep_optzs (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, NcHaloDensityProfile *dp, const gdouble r_s, const gdouble rho_s, NcWLSurfaceMassDensityOptzs *optzs);
NcGalaxyShapeFactorLine *nc_galaxy_shape_factor_line_new_full (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const NcWLSurfaceMassDensityOptzs *optzs, const gdouble z_lo, const gdouble z_hi);

G_END_DECLS

#endif /* NUMCOSMO_GIR_SCAN */

#endif /* _NC_GALAXY_SHAPE_FACTOR_REACH_PRIVATE_H_ */

