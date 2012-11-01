/***************************************************************************
 *            ncm_sf_sbessel.h
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

#ifndef _NCM_SF_SBESSEL_H
#define _NCM_SF_SBESSEL_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/grid_one.h>
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

typedef struct _NcSFSBesselRecur NcSFSBesselRecur;

/**
 * NcSFSBesselRecur:
 *
 * FIXME
 */
struct _NcSFSBesselRecur
{
  /*< private >*/
  gint32 l;
  NcmGrid *x_grid;
  gdouble *jl;
  gdouble *jlp1;
  gboolean prepared;
};

NcSFSBesselRecur *ncm_sf_sbessel_recur_new (NcmGrid *x_grid);
NcSFSBesselRecur *ncm_sf_sbessel_recur_read (FILE *f);

void ncm_sf_sbessel_recur_free (NcSFSBesselRecur *jlrec, gboolean free_grid);
void ncm_sf_sbessel_recur_set (NcSFSBesselRecur *jlrec, glong l);
void ncm_sf_sbessel_recur_next (NcSFSBesselRecur *jlrec);
void ncm_sf_sbessel_recur_previous (NcSFSBesselRecur *jlrec);
void ncm_sf_sbessel_recur_goto (NcSFSBesselRecur *jlrec, glong l);
void ncm_sf_sbessel_taylor_coeff_jl_jlp1 (NcSFSBesselRecur *jlrec, guint n, gdouble *djl, gdouble *djlp1);
void ncm_sf_sbessel_recur_write (NcSFSBesselRecur *jlrec, FILE *f);

gdouble ncm_sf_sbessel (gulong l, gdouble x);
void ncm_sf_sbessel_taylor (gulong l, gdouble x, gdouble *djl);
void ncm_sf_sbessel_deriv (gulong l, gdouble x, gdouble jl, gdouble jlp1, gdouble *djl);

NcmSpline *ncm_sf_sbessel_spline (gulong l, gdouble xi, gdouble xf, gdouble reltol);

G_END_DECLS

#endif /* _NCM_SF_SBESSEL_H */
