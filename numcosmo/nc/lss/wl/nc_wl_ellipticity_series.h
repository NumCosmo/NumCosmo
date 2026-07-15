/***************************************************************************
 *            nc_wl_ellipticity_series.h
 *
 *  Wed Jul 15 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_wl_ellipticity_series.h
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
#ifndef _NC_WL_ELLIPTICITY_SERIES_H_
#define _NC_WL_ELLIPTICITY_SERIES_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_laurent_series.h>

G_BEGIN_DECLS

/*
 * NcWLEllipticitySeriesTrace / NcWLEllipticitySeriesTraceDet: the truncated
 * g-Taylor series counterparts of nc_wl_ellipticity.h's own closed-form
 * kernels, one object per ellipticity convention (matching that file's own
 * _trace/_trace_det split): chi_I(chi_L,g) via apply_shear_inv, and the
 * (linear, not log) Jacobian via det_jac -- see nc_wl_ellipticity_series.c's
 * own module doc comment for the derivation and dev-notes reference.
 *
 * Both are ordinary long-lived GObjects: constructed once with a fixed
 * truncation order (CONSTRUCT_ONLY "order" property), then eval()'d
 * repeatedly at whatever rho a caller needs, with get_chi()/get_jac()
 * reading back the most recent result (transfer-none, valid until the next
 * eval() call). All scratch this needs (both the chi/jac series themselves
 * and every intermediate the derivation uses) is owned internally and
 * allocated once at construction time.
 */

#define NC_TYPE_WL_ELLIPTICITY_SERIES_TRACE (nc_wl_ellipticity_series_trace_get_type ())

G_DECLARE_FINAL_TYPE (NcWLEllipticitySeriesTrace, nc_wl_ellipticity_series_trace, NC, WL_ELLIPTICITY_SERIES_TRACE, GObject)

NcWLEllipticitySeriesTrace *nc_wl_ellipticity_series_trace_new (guint order);
NcWLEllipticitySeriesTrace *nc_wl_ellipticity_series_trace_ref (NcWLEllipticitySeriesTrace *ser);
void nc_wl_ellipticity_series_trace_free (NcWLEllipticitySeriesTrace *ser);
void nc_wl_ellipticity_series_trace_clear (NcWLEllipticitySeriesTrace **ser);

guint nc_wl_ellipticity_series_trace_get_order (NcWLEllipticitySeriesTrace *ser);
void nc_wl_ellipticity_series_trace_eval (NcWLEllipticitySeriesTrace *ser, gdouble rho);
NcmLaurentSeriesTPS *nc_wl_ellipticity_series_trace_get_chi (NcWLEllipticitySeriesTrace *ser);
NcmLaurentSeriesTPS *nc_wl_ellipticity_series_trace_get_jac (NcWLEllipticitySeriesTrace *ser);
NcmLaurentSeriesTPS *nc_wl_ellipticity_series_trace_get_abs_sq (NcWLEllipticitySeriesTrace *ser);

#define NC_TYPE_WL_ELLIPTICITY_SERIES_TRACE_DET (nc_wl_ellipticity_series_trace_det_get_type ())

G_DECLARE_FINAL_TYPE (NcWLEllipticitySeriesTraceDet, nc_wl_ellipticity_series_trace_det, NC, WL_ELLIPTICITY_SERIES_TRACE_DET, GObject)

NcWLEllipticitySeriesTraceDet *nc_wl_ellipticity_series_trace_det_new (guint order);
NcWLEllipticitySeriesTraceDet *nc_wl_ellipticity_series_trace_det_ref (NcWLEllipticitySeriesTraceDet *ser);
void nc_wl_ellipticity_series_trace_det_free (NcWLEllipticitySeriesTraceDet *ser);
void nc_wl_ellipticity_series_trace_det_clear (NcWLEllipticitySeriesTraceDet **ser);

guint nc_wl_ellipticity_series_trace_det_get_order (NcWLEllipticitySeriesTraceDet *ser);
void nc_wl_ellipticity_series_trace_det_eval (NcWLEllipticitySeriesTraceDet *ser, gdouble rho);
NcmLaurentSeriesTPS *nc_wl_ellipticity_series_trace_det_get_chi (NcWLEllipticitySeriesTraceDet *ser);
NcmLaurentSeriesTPS *nc_wl_ellipticity_series_trace_det_get_jac (NcWLEllipticitySeriesTraceDet *ser);
NcmLaurentSeriesTPS *nc_wl_ellipticity_series_trace_det_get_abs_sq (NcWLEllipticitySeriesTraceDet *ser);

G_END_DECLS

#endif /* _NC_WL_ELLIPTICITY_SERIES_H_ */

