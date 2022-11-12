/***************************************************************************
 *            ncm_mpsf_sbessel.h
 *
 *  Mon Nov 23 10:05:14 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_MPSF_SBESSEL_H
#define _NCM_MPSF_SBESSEL_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

typedef struct _NcmMpsfSBesselRecur NcmMpsfSBesselRecur;

/**
 * NcmMpsfSBesselRecur:
 *
 * FIXME
 */
struct _NcmMpsfSBesselRecur
{
  /*< private >*/
  guint32 prec;
  gint32 l;
  mpq_t q;
  mpfr_t x;
  mpfr_t jl[2];
  mpfr_t temp;
};

void ncm_mpsf_sbessel (gulong l, mpq_t q, mpfr_ptr res, mp_rnd_t rnd);
void ncm_mpsf_sbessel_d (gulong l, gdouble x, mpfr_ptr res, mp_rnd_t rnd);
void ncm_mpsf_sbessel_free_cache (void);

G_END_DECLS

#endif /* _NCM_MPSF_SBESSEL_H */

