/***************************************************************************
 *            nc_hicosmo_de_wpoly_library.h
 *
 *  Mon Sep 28 16:00:00 2026
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

#ifndef _NC_HICOSMO_DE_WPOLY_LIBRARY_H_
#define _NC_HICOSMO_DE_WPOLY_LIBRARY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/background/nc_hicosmo_de_wpoly.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_WPOLY_LIBRARY (nc_hicosmo_de_wpoly_library_get_type ())

G_DECLARE_FINAL_TYPE (NcHICosmoDEWPolyLibrary, nc_hicosmo_de_wpoly_library, NC, HICOSMO_DE_WPOLY_LIBRARY, GObject)

NcHICosmoDEWPolyLibrary *nc_hicosmo_de_wpoly_library_new (gulong seed, guint degree, const gdouble z_f, guint n_per_class);
NcHICosmoDEWPolyLibrary *nc_hicosmo_de_wpoly_library_ref (NcHICosmoDEWPolyLibrary *lib);
void nc_hicosmo_de_wpoly_library_free (NcHICosmoDEWPolyLibrary *lib);
void nc_hicosmo_de_wpoly_library_clear (NcHICosmoDEWPolyLibrary **lib);

void nc_hicosmo_de_wpoly_library_prepare (NcHICosmoDEWPolyLibrary *lib);

guint nc_hicosmo_de_wpoly_library_get_max_class (NcHICosmoDEWPolyLibrary *lib);
guint nc_hicosmo_de_wpoly_library_len (NcHICosmoDEWPolyLibrary *lib, guint m);
gulong nc_hicosmo_de_wpoly_library_get_n_draws (NcHICosmoDEWPolyLibrary *lib);

gulong nc_hicosmo_de_wpoly_library_get_draw_index (NcHICosmoDEWPolyLibrary *lib, guint m, guint i);
gdouble nc_hicosmo_de_wpoly_library_get_amplitude (NcHICosmoDEWPolyLibrary *lib, guint m, guint i);
GArray *nc_hicosmo_de_wpoly_library_peek_turning_points (NcHICosmoDEWPolyLibrary *lib, guint m, guint i);
NcmVector *nc_hicosmo_de_wpoly_library_peek_coeffs (NcHICosmoDEWPolyLibrary *lib, guint m, guint i);
NcHICosmoDEWPoly *nc_hicosmo_de_wpoly_library_get_fiducial (NcHICosmoDEWPolyLibrary *lib, guint m, guint i);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_WPOLY_LIBRARY_H_ */

