/***************************************************************************
 *            hydrodyn_adiabatic.h
 *
 *  Mon Jun  7 14:54:02 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

#ifndef _NC_HYDRODYN_ADIABATIC_H
#define _NC_HYDRODYN_ADIABATIC_H

#include <glib.h>
#include <glib-object.h>

#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

typedef struct _NcPertHydrodyn
{
  NcHICosmo *model;
  gpointer cvode;
  gboolean initialized;
  N_Vector y;
  N_Vector yQ;
  gdouble reltol;
  gdouble abstol;
  long double k;
  long double alpha0;
  long double alphai;
  long double alphaf;
  NcmSpline *pw_spline;
} NcPertHydrodyn;

void nc_hydrodyn_adiabatic_stub (void);

G_END_DECLS

#endif /* _NC_HYDRODYN_ADIABATIC_H */
