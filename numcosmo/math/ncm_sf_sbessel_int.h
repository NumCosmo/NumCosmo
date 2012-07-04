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
#include <string.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <glib.h>

G_BEGIN_DECLS

typedef struct _NcSFSphericalBesselIntegRecur NcSFSphericalBesselIntegRecur;

/**
 * NcSFSphericalBesselIntegRecur:
 *
 * FIXME
 */
struct _NcSFSphericalBesselIntegRecur
{
  /*< private >*/
  NcSFSBesselRecur *jlrec;
  gdouble *int_jl_xn[4];
  gdouble *int_jlp1_xn[4];
  gboolean prepared;
};

typedef struct _NcSFSphericalBesselIntSpline NcSFSphericalBesselIntSpline;

/**
 * NcSFSphericalBesselIntSpline:
 *
 * FIXME
 */
struct _NcSFSphericalBesselIntSpline
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
  NcSFSphericalBesselIntegRecur *xnjlrec;
  gboolean prepared;
};

gdouble ncm_sf_sbessel_jl_xj_integral (gint l, gint j, gdouble x);

NcSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_new (NcSFSBesselRecur *jlrec, NcmGrid *x_grid);
NcSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_new_from_section (NcmGridSection *x_sec);
NcSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_load (gchar *filename);
NcSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_read (FILE *f);
NcSFSphericalBesselIntegRecur *ncm_sf_sbessel_jl_xj_integral_recur_cached_new (glong l, NcmGridSection *x_sec);

void ncm_sf_sbessel_jl_xj_integral_recur_free (NcSFSphericalBesselIntegRecur *xnjlrec, gboolean free_grid);
void ncm_sf_sbessel_jl_xj_integral_recur_set (NcSFSphericalBesselIntegRecur *xnjlrec, glong l);
glong ncm_sf_sbessel_jl_xj_integral_recur_next (NcSFSphericalBesselIntegRecur *xnjlrec);
glong ncm_sf_sbessel_jl_xj_integral_recur_previous (NcSFSphericalBesselIntegRecur *xnjlrec);
glong ncm_sf_sbessel_jl_xj_integral_recur_goto (NcSFSphericalBesselIntegRecur *xnjlrec, glong l);

void ncm_sf_sbessel_jl_xj_integral_recur_save (NcSFSphericalBesselIntegRecur *xnjlrec, gchar *filename);
void ncm_sf_sbessel_jl_xj_integral_recur_write (NcSFSphericalBesselIntegRecur *xnjlrec, FILE *f);

NcSFSphericalBesselIntSpline *ncm_sf_sbessel_jl_xj_integrate_spline_new (NcSFSphericalBesselIntegRecur *xnjlrec, gboolean init);
NcSFSphericalBesselIntSpline *ncm_sf_sbessel_jl_xj_integrate_spline_cached_new (glong l, NcmGridSection *x_sec, gboolean init);

void ncm_sf_sbessel_jl_xj_integrate_spline_reset (NcSFSphericalBesselIntSpline *int_jlspline);
void ncm_sf_sbessel_jl_xj_integrate_spline_set (NcSFSphericalBesselIntSpline *int_jlspline, glong l);
void ncm_sf_sbessel_jl_xj_integrate_spline_next (NcSFSphericalBesselIntSpline *int_jlspline);
void ncm_sf_sbessel_jl_xj_integrate_spline_previous (NcSFSphericalBesselIntSpline *int_jlspline);
void ncm_sf_sbessel_jl_xj_integrate_spline_goto (NcSFSphericalBesselIntSpline *int_jlspline, glong l);
void ncm_sf_sbessel_jl_xj_integral_a_b (NcSFSphericalBesselIntSpline *int_jlspline, gdouble x0, gdouble x1, gdouble w, gdouble *xnjl_rules, gdouble *xndjl_rules, gdouble *xnd2jl_rules);
gdouble ncm_sf_sbessel_jl_xj_integrate_spline_eval (NcSFSphericalBesselIntSpline *int_jlspline, gint d, gdouble x);
gdouble ncm_sf_sbessel_jl_xj_integral_spline (NcSFSphericalBesselIntSpline *int_jlspline, NcmSpline *s0, NcmSpline *s1, NcmSpline *s2, gdouble w);

G_END_DECLS

#endif /* _NCM_SF_SBESSEL_INT_H */
