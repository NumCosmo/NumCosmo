/***************************************************************************
 *            dividedifference.h
 *
 *  Wed Mar 17 16:20:46 2010
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

#ifndef _NC_DIVIDEDIFFERENCE_H
#define _NC_DIVIDEDIFFERENCE_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

void ncm_interp_dd_init (const gdouble *vx, gdouble *dd, const gint np, const gint nf);
gdouble ncm_interp_dd_eval (const gdouble *vx, const gdouble *dd, const gdouble x, const gint np, const gint nf);

void ncm_interp_dd_init_2_4 (const gdouble *vx, gdouble *dd);
gdouble ncm_interp_dd_eval_2_4 (const gdouble *vx, const gdouble *dd, const gdouble x);

G_END_DECLS

#endif /* _NC_DIVIDEDIFFERENCE_H */
