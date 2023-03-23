/***************************************************************************
 *            ncm_mpsf_trig_int.h
 *
 *  Tue Feb  2 22:16:05 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NCM_MPSF_TRIG_INT_H
#define _NCM_MPSF_TRIG_INT_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gmp.h>
#include <mpfr.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

void ncm_mpsf_sin_int_mpfr (mpq_t q, mpfr_ptr res, mp_rnd_t rnd);
gdouble ncm_sf_sin_int (gdouble x);

G_END_DECLS

#endif /* _NCM_MPSF_TRIG_INT_H */

