/***************************************************************************
 *            nc_galaxy_redshift_binning.h
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_binning.h
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

#ifndef _NC_GALAXY_REDSHIFT_BINNING_H_
#define _NC_GALAXY_REDSHIFT_BINNING_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/spline/ncm_spline.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_pop.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_pop_lsst_srd.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_obs_sel.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_BINNING (nc_galaxy_redshift_binning_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyRedshiftBinning, nc_galaxy_redshift_binning, NC, GALAXY_REDSHIFT_BINNING, GObject)

NcGalaxyRedshiftBinning *nc_galaxy_redshift_binning_new (void);
NcGalaxyRedshiftBinning *nc_galaxy_redshift_binning_ref (NcGalaxyRedshiftBinning *gsdrb);

void nc_galaxy_redshift_binning_free (NcGalaxyRedshiftBinning *gsdrb);
void nc_galaxy_redshift_binning_clear (NcGalaxyRedshiftBinning **gsdrb);

void nc_galaxy_redshift_binning_set_reltol (NcGalaxyRedshiftBinning *gsdrb, const gdouble reltol);
gdouble nc_galaxy_redshift_binning_get_reltol (NcGalaxyRedshiftBinning *gsdrb);

void nc_galaxy_redshift_binning_set_zp_support_max (NcGalaxyRedshiftBinning *gsdrb, const gdouble zp_support_max);
gdouble nc_galaxy_redshift_binning_get_zp_support_max (NcGalaxyRedshiftBinning *gsdrb);

NcmSpline *nc_galaxy_redshift_binning_compute_dndz (NcGalaxyRedshiftBinning *gsdrb, NcGalaxyRedshiftPop *population, NcGalaxyRedshiftObsSel *observable_population, const gdouble zp_min, const gdouble zp_max);
NcmSpline *nc_galaxy_redshift_binning_compute_dndz_on_nodes (NcGalaxyRedshiftBinning *gsdrb, NcGalaxyRedshiftPop *population, NcGalaxyRedshiftObsSel *observable_population, const gdouble zp_min, const gdouble zp_max, NcmVector *z_nodes);

void nc_galaxy_redshift_binning_prepare (NcGalaxyRedshiftBinning *gsdrb, NcGalaxyRedshiftPop *population, NcGalaxyRedshiftObsSel *observable_population);

gdouble nc_galaxy_redshift_binning_eval_pzp (NcGalaxyRedshiftBinning *gsdrb, const gdouble zp);
NcmVector *nc_galaxy_redshift_binning_compute_equal_area_photoz_bins (NcGalaxyRedshiftBinning *gsdrb, const guint n_bins, const gdouble zp_max);

NcmVector *nc_galaxy_redshift_binning_lsst_srd_edges (NcGalaxyRedshiftPopLSSTSRDType type, NcGalaxyRedshiftPop **population, NcGalaxyRedshiftObsSel **observable_population);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_BINNING_H_ */

