/***************************************************************************
 *            nc_hcosmo_bi_dist.h
 *
 *  Sat Oct 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hcosmo_bi_dist.h
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

#ifndef _NC_HCOSMO_BI_DIST_H_
#define _NC_HCOSMO_BI_DIST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/background/nc_hcosmo_bi.h>

G_BEGIN_DECLS

#define NC_TYPE_HCOSMO_BI_DIST (nc_hcosmo_bi_dist_get_type ())

G_DECLARE_FINAL_TYPE (NcHCosmoBIDist, nc_hcosmo_bi_dist, NC, HCOSMO_BI_DIST, GObject)

NcHCosmoBIDist *nc_hcosmo_bi_dist_new (const gdouble zmax);
NcHCosmoBIDist *nc_hcosmo_bi_dist_ref (NcHCosmoBIDist *dist);

void nc_hcosmo_bi_dist_free (NcHCosmoBIDist *dist);
void nc_hcosmo_bi_dist_clear (NcHCosmoBIDist **dist);

void nc_hcosmo_bi_dist_set_zmax (NcHCosmoBIDist *dist, const gdouble zmax);
gdouble nc_hcosmo_bi_dist_get_zmax (NcHCosmoBIDist *dist);
void nc_hcosmo_bi_dist_set_reltol (NcHCosmoBIDist *dist, const gdouble reltol);
gdouble nc_hcosmo_bi_dist_get_reltol (NcHCosmoBIDist *dist);

void nc_hcosmo_bi_dist_prepare (NcHCosmoBIDist *dist, NcHCosmoBI *cosmo);
void nc_hcosmo_bi_dist_prepare_if_needed (NcHCosmoBIDist *dist, NcHCosmoBI *cosmo);

gdouble nc_hcosmo_bi_dist_get_t_lb_max (NcHCosmoBIDist *dist);
guint nc_hcosmo_bi_dist_get_n_nodes (NcHCosmoBIDist *dist);

void nc_hcosmo_bi_dist_sky_to_axes (NcHCosmoBIDist *dist, const gdouble px, const gdouble py, const gdouble pz, gdouble *n1, gdouble *n2, gdouble *n3);
gdouble nc_hcosmo_bi_dist_t_lb_emission (NcHCosmoBIDist *dist, const gdouble z, const gdouble n1, const gdouble n2, const gdouble n3);
gdouble nc_hcosmo_bi_dist_luminosity (NcHCosmoBIDist *dist, const gdouble z, const gdouble n1, const gdouble n2, const gdouble n3);
gdouble nc_hcosmo_bi_dist_dmodulus_sky (NcHCosmoBIDist *dist, const gdouble z, const gdouble px, const gdouble py, const gdouble pz);

#define NC_HCOSMO_BI_DIST_DEFAULT_RELTOL (1.0e-11)

G_END_DECLS

#endif /* _NC_HCOSMO_BI_DIST_H_ */

