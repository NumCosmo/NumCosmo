/***************************************************************************
 *            ncm_spline_hermite.h
 *
 *  Mon Oct 12 2026
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

#ifndef _NCM_SPLINE_HERMITE_H_
#define _NCM_SPLINE_HERMITE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/spline/ncm_spline.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_HERMITE (ncm_spline_hermite_get_type ())

G_DECLARE_FINAL_TYPE (NcmSplineHermite, ncm_spline_hermite, NCM, SPLINE_HERMITE, NcmSpline)

NcmSplineHermite *ncm_spline_hermite_new (void);
NcmSplineHermite *ncm_spline_hermite_new_full (NcmVector *xv, NcmVector *yv, NcmVector *d1v, NcmVector *d2v, gboolean init);

void ncm_spline_hermite_set_derivs (NcmSplineHermite *sh, NcmVector *d1v, NcmVector *d2v, gboolean init);
NcmVector *ncm_spline_hermite_peek_d1v (NcmSplineHermite *sh);
NcmVector *ncm_spline_hermite_peek_d2v (NcmSplineHermite *sh);

G_END_DECLS

#endif /* _NCM_SPLINE_HERMITE_H_ */

