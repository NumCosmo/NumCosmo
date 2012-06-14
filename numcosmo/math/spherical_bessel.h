/***************************************************************************
 *            spherical-bessel.h
 *
 *  Wed Mar 10 17:15:46 2010
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
 
/**
 * @file
 * @brief FIXME
 *
 * FIXME
 */

#ifndef _NC_SF_SPHERICAL_BESSEL_H
#define _NC_SF_SPHERICAL_BESSEL_H
#include <string.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <glib.h>

G_BEGIN_DECLS

typedef struct _NcSFSpherBesselRecur NcSFSpherBesselRecur;

/**
 * NcSFSpherBesselRecur:
 * 
 * FIXME
 */
struct _NcSFSpherBesselRecur
{
  /*< private >*/
  gint32 l;
  NcmGrid *x_grid;
  gdouble *jl;
  gdouble *jlp1;
  gboolean prepared;
};

NcSFSpherBesselRecur *ncm_sf_spherical_bessel_recur_new (NcmGrid *x_grid);
NcSFSpherBesselRecur *ncm_sf_spherical_bessel_recur_read (FILE *f);

void ncm_sf_spherical_bessel_recur_free (NcSFSpherBesselRecur *jlrec, gboolean free_grid);
void ncm_sf_spherical_bessel_recur_set (NcSFSpherBesselRecur *jlrec, glong l);
void ncm_sf_spherical_bessel_recur_next (NcSFSpherBesselRecur *jlrec);
void ncm_sf_spherical_bessel_recur_previous (NcSFSpherBesselRecur *jlrec);
void ncm_sf_spherical_bessel_recur_goto (NcSFSpherBesselRecur *jlrec, glong l);
void ncm_sf_spherical_bessel_taylor_coeff_jl_jlp1 (NcSFSpherBesselRecur *jlrec, guint n, gdouble *djl, gdouble *djlp1);
void ncm_sf_spherical_bessel_recur_write (NcSFSpherBesselRecur *jlrec, FILE *f);

gdouble ncm_sf_spherical_bessel (gulong l, gdouble x);
void ncm_sf_spherical_bessel_taylor (gulong l, gdouble x, gdouble *djl);
void ncm_sf_spherical_bessel_deriv (gulong l, gdouble x, gdouble jl, gdouble jlp1, gdouble *djl);

NcmSpline *ncm_sf_spherical_bessel_spline (gulong l, gdouble xi, gdouble xf, gdouble reltol);

G_END_DECLS

#endif /* _NC_SF_SPHERICAL_BESSEL_H */
