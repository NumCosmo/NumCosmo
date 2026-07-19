/***************************************************************************
 *            ncm_poly_roots.h
 *
 *  Mon Jul 6 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_poly_roots.h
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _NCM_POLY_ROOTS_H_
#define _NCM_POLY_ROOTS_H_

#include <glib.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

gint ncm_poly_roots_real_quadratic (const gdouble a[3], gdouble roots[2]);

gint ncm_poly_roots_real_cubic (const gdouble a[4], gdouble roots[3]);
gint ncm_poly_roots_real_quartic (const gdouble a[5], gdouble roots[4]);
gint ncm_poly_roots_real_quartic_or_lower (const gdouble a[5], gdouble roots[4]);

G_END_DECLS

#endif /* _NCM_POLY_ROOTS_H_ */

