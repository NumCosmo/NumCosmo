/***************************************************************************
 *            poly.h
 *
 *  Wed Aug 13 21:31:16 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NC_POLY_H
#define _NC_POLY_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_vector.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

gsl_vector *ncm_poly_new (gint degree);
gdouble ncm_poly_eval (gsl_vector *poly, gdouble x);
gdouble ncm_poly_eval_diff (gsl_vector *poly, guint n, gdouble x);
gdouble ncm_poly_eval_int (gsl_vector *poly, gdouble x);
gdouble ncm_poly_eval_int_P_over_x (gsl_vector *poly, gdouble x);

G_END_DECLS

#endif /* _NC_POLY_H */
