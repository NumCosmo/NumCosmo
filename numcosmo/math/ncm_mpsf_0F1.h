/***************************************************************************
 *            ncm_mpsf_0F1.h
 *
 *  Thu October 13 22:17:46 2011
 *  Copyright  2011  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_MPSF_0F1_H
#define _NCM_MPSF_0F1_H

#include <glib.h>
#include <glib-object.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gmp.h>
#include <mpfr.h>
#endif

G_BEGIN_DECLS

void ncm_mpsf_0F1_q (mpq_t b, mpq_t x, mpfr_ptr res, mp_rnd_t rnd);
void ncm_mpsf_0F1_d (gdouble b, gdouble x, mpfr_ptr res, mp_rnd_t rnd);
gdouble ncm_sf_0F1 (gdouble b, gdouble x);

G_END_DECLS

#endif /* _NCM_MPSF_0F1_H */
