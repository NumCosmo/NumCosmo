/***************************************************************************
 *            nc_xcor_ssc_sij.h
 *
 *  Thu August 13 2026
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

#ifndef _NC_XCOR_SSC_SIJ_H_
#define _NC_XCOR_SSC_SIJ_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/background/nc_distance.h>
#include <numcosmo/nc/background/nc_hicosmo.h>
#include <numcosmo/nc/xcor/nc_xcor.h>
#include <numcosmo/ncm/powspec/ncm_powspec.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_SSC_SIJ (nc_xcor_ssc_sij_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorSSCSij, nc_xcor_ssc_sij, NC, XCOR_SSC_SIJ, GObject)

NcXcorSSCSij *nc_xcor_ssc_sij_new (NcDistance * dist, NcmPowspec * ps, NcmVector * z_edges);
NcXcorSSCSij *nc_xcor_ssc_sij_ref (NcXcorSSCSij *ssc_sij);
void nc_xcor_ssc_sij_free (NcXcorSSCSij *ssc_sij);
void nc_xcor_ssc_sij_clear (NcXcorSSCSij **ssc_sij);

guint nc_xcor_ssc_sij_get_nbins (NcXcorSSCSij *ssc_sij);
guint nc_xcor_ssc_sij_get_lmax (NcXcorSSCSij *ssc_sij);
gdouble nc_xcor_ssc_sij_get_fsky (NcXcorSSCSij *ssc_sij);

void nc_xcor_ssc_sij_set_mask_cl (NcXcorSSCSij *ssc_sij, NcmVector *mask_cl);
NcmVector *nc_xcor_ssc_sij_peek_mask_cl (NcXcorSSCSij *ssc_sij);

void nc_xcor_ssc_sij_set_area (NcXcorSSCSij *ssc_sij, gdouble area);
gdouble nc_xcor_ssc_sij_get_area (NcXcorSSCSij *ssc_sij);

void nc_xcor_ssc_sij_set_block_size (NcXcorSSCSij *ssc_sij, guint block_size);
guint nc_xcor_ssc_sij_get_block_size (NcXcorSSCSij *ssc_sij);

void nc_xcor_ssc_sij_set_reltol (NcXcorSSCSij *ssc_sij, gdouble reltol);
gdouble nc_xcor_ssc_sij_get_reltol (NcXcorSSCSij *ssc_sij);

void nc_xcor_ssc_sij_set_scaled_abstol (NcXcorSSCSij *ssc_sij, gdouble scaled_abstol);
gdouble nc_xcor_ssc_sij_get_scaled_abstol (NcXcorSSCSij *ssc_sij);

void nc_xcor_ssc_sij_set_method (NcXcorSSCSij *ssc_sij, NcXcorMethod method);
NcXcorMethod nc_xcor_ssc_sij_get_method (NcXcorSSCSij *ssc_sij);

void nc_xcor_ssc_sij_prepare (NcXcorSSCSij *ssc_sij, NcHICosmo *cosmo);
void nc_xcor_ssc_sij_prepare_if_needed (NcXcorSSCSij *ssc_sij, NcHICosmo *cosmo);

NcmMatrix *nc_xcor_ssc_sij_peek_matrix (NcXcorSSCSij *ssc_sij);
NcmMatrix *nc_xcor_ssc_sij_eval (NcXcorSSCSij *ssc_sij, NcHICosmo *cosmo);

NcmVector *nc_xcor_ssc_sij_mask_cl_fullsky (void);

G_END_DECLS

#endif /* _NC_XCOR_SSC_SIJ_H_ */

