/***************************************************************************
 *            ncm_sf_sbessel_int.h
 *
 *  Wed Mar 10 17:16:29 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_SF_SBESSEL_INT_H
#define _NCM_SF_SBESSEL_INT_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_sf_sbessel.h>

G_BEGIN_DECLS

typedef struct _NcmSFSphericalBesselIntegRecur NcmSFSphericalBesselIntegRecur;

/**
 * NcmSFSphericalBesselIntegRecur:
 *
 * FIXME
 */
struct _NcmSFSphericalBesselIntegRecur
{
  /*< private >*/
  NcmSFSBesselRecur *jlrec;
  gdouble *int_jl_xn[4];
  gdouble *int_jlp1_xn[4];
  gboolean prepared;
};

typedef struct _NcmSFSphericalBesselIntSpline NcmSFSphericalBesselIntSpline;

/**
 * NcmSFSphericalBesselIntSpline:
 *
 * FIXME
 */
struct _NcmSFSphericalBesselIntSpline
{
  /*< private >*/
  NcmGrid *x_grid;
  gdouble *x_data;
  gdouble *dd_jl;
  gdouble *dd_jlp1;
  gdouble *dd_int_jl_x0;
  gdouble *dd_int_jl_x1;
  gdouble *dd_int_jl_x2;
  gdouble *dd_int_jl_x3;
  NcmSpline *jl;
  NcmSpline *jlp1;
  NcmSpline *int_jl_xn[4];
  NcmSFSphericalBesselIntegRecur *xnjlrec;
  gboolean prepared;
};

gdouble ncm_sf_sbessel_jl_xj_integral (gint l, gint j, gdouble x);

NcmSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_new (NcmSFSBesselRecur *jlrec, NcmGrid *x_grid);
NcmSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_new_from_section (NcmGridSection *x_sec);
NcmSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_load (gchar *filename);
NcmSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_read (FILE *f);
NcmSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_cached_new (glong l, NcmGridSection *x_sec);

void ncm_sf_sbessel_jl_xj_integral_recur_free (NcmSFSphericalBesselIntegRecur *xnjlrec, gboolean free_grid);
void ncm_sf_sbessel_jl_xj_integral_recur_set (NcmSFSphericalBesselIntegRecur *xnjlrec, glong l);
glong ncm_sf_sbessel_jl_xj_integral_recur_next (NcmSFSphericalBesselIntegRecur *xnjlrec);
glong ncm_sf_sbessel_jl_xj_integral_recur_previous (NcmSFSphericalBesselIntegRecur *xnjlrec);
glong ncm_sf_sbessel_jl_xj_integral_recur_goto (NcmSFSphericalBesselIntegRecur *xnjlrec, glong l);

void ncm_sf_sbessel_jl_xj_integral_recur_save (NcmSFSphericalBesselIntegRecur *xnjlrec, gchar *filename);
void ncm_sf_sbessel_jl_xj_integral_recur_write (NcmSFSphericalBesselIntegRecur *xnjlrec, FILE *f);

NcmSFSphericalBesselIntSpline *ncm_sf_sbessel_jl_xj_integrate_spline_new (NcmSFSphericalBesselIntegRecur *xnjlrec, gboolean init);
NcmSFSphericalBesselIntSpline *ncm_sf_sbessel_jl_xj_integrate_spline_cached_new (glong l, NcmGridSection *x_sec, gboolean init);

void ncm_sf_sbessel_jl_xj_integrate_spline_reset (NcmSFSphericalBesselIntSpline *int_jlspline);
void ncm_sf_sbessel_jl_xj_integrate_spline_set (NcmSFSphericalBesselIntSpline *int_jlspline, glong l);
void ncm_sf_sbessel_jl_xj_integrate_spline_next (NcmSFSphericalBesselIntSpline *int_jlspline);
void ncm_sf_sbessel_jl_xj_integrate_spline_previous (NcmSFSphericalBesselIntSpline *int_jlspline);
void ncm_sf_sbessel_jl_xj_integrate_spline_goto (NcmSFSphericalBesselIntSpline *int_jlspline, glong l);
void ncm_sf_sbessel_jl_xj_integral_a_b (NcmSFSphericalBesselIntSpline *int_jlspline, gdouble x0, gdouble x1, gdouble w, gdouble *xnjl_rules, gdouble *xndjl_rules, gdouble *xnd2jl_rules);
gdouble ncm_sf_sbessel_jl_xj_integrate_spline_eval (NcmSFSphericalBesselIntSpline *int_jlspline, gint d, gdouble x);
gdouble ncm_sf_sbessel_jl_xj_integral_spline (NcmSFSphericalBesselIntSpline *int_jlspline, NcmSpline *s0, NcmSpline *s1, NcmSpline *s2, gdouble w);

G_END_DECLS

#endif /* _NCM_SF_SBESSEL_INT_H */
