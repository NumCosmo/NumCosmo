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

#include <stdio.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <gmp.h>
#include <mpfr.h>
#endif

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

NcmMpsfSBesselRecur *ncm_mpsf_sbessel_recur_new (gulong prec);
NcmMpsfSBesselRecur *ncm_mpsf_sbessel_recur_read (FILE *f);
void ncm_mpsf_sbessel_recur_free (NcmMpsfSBesselRecur *jlrec);
void ncm_mpsf_sbessel_recur_set_q (NcmMpsfSBesselRecur *jlrec, glong l, mpq_ptr q, mp_rnd_t rnd);
void ncm_mpsf_sbessel_recur_set_d (NcmMpsfSBesselRecur *jlrec, glong l, gdouble x, mp_rnd_t rnd);
void ncm_mpsf_sbessel_recur_next (NcmMpsfSBesselRecur *jlrec, mp_rnd_t rnd);
void ncm_mpsf_sbessel_recur_previous (NcmMpsfSBesselRecur *jlrec, mp_rnd_t rnd);
void ncm_mpsf_sbessel_recur_goto (NcmMpsfSBesselRecur *jlrec, glong l, mp_rnd_t rnd);
void ncm_mpsf_sbessel_recur_write (NcmMpsfSBesselRecur *jlrec, FILE *f);

G_END_DECLS

#endif /* _NCM_MPSF_SBESSEL_H */
