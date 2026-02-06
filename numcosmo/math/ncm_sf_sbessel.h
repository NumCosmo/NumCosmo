/***************************************************************************
 *            ncm_sf_sbessel.h
 *
 *  Wed Mar 10 17:15:46 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NCM_SF_SBESSEL_H
#define _NCM_SF_SBESSEL_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

#define NCM_TYPE_SF_SBESSEL_ARRAY (ncm_sf_sbessel_array_get_type ())

G_DECLARE_FINAL_TYPE (NcmSFSBesselArray, ncm_sf_sbessel_array, NCM, SF_SBESSEL_ARRAY, GObject)

NcmSFSBesselArray *ncm_sf_sbessel_array_new (void);
NcmSFSBesselArray *ncm_sf_sbessel_array_new_full (guint lmax, gdouble threshold);
NcmSFSBesselArray *ncm_sf_sbessel_array_ref (NcmSFSBesselArray *sba);
void ncm_sf_sbessel_array_free (NcmSFSBesselArray *sba);
void ncm_sf_sbessel_array_clear (NcmSFSBesselArray **sba);

void ncm_sf_sbessel_array_eval (NcmSFSBesselArray *sba, guint ell, gdouble x, gdouble *jl_x);
GArray *ncm_sf_sbessel_array_eval1 (NcmSFSBesselArray *sba, guint ell, gdouble x);
guint ncm_sf_sbessel_array_eval_ell_cutoff (NcmSFSBesselArray *sba, gdouble x);
guint ncm_sf_sbessel_array_get_lmax (NcmSFSBesselArray *sba);
gdouble ncm_sf_sbessel_array_get_threshold (NcmSFSBesselArray *sba);

/* Stand alone functions */

gdouble ncm_sf_sbessel (gulong l, gdouble x);

void ncm_sf_sbessel_taylor (gulong l, gdouble x, gdouble *djl);
NcmSpline *ncm_sf_sbessel_spline (gulong l, gdouble xi, gdouble xf, gdouble reltol);

G_END_DECLS

#endif /* _NCM_SF_SBESSEL_H */

