/***************************************************************************
 *            mp-spherical-bessel.h
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

#ifndef _NC_SF_MP_SPHERICAL_BESSEL_H
#define _NC_SF_MP_SPHERICAL_BESSEL_H
#include <string.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <glib.h>

G_BEGIN_DECLS

typedef struct _NcmSfMpSpherBesselRecur NcmSfMpSpherBesselRecur;

/**
 * NcmSfMpSpherBesselRecur:
 * 
 * FIXME
 */
struct _NcmSfMpSpherBesselRecur
{
  /*< private >*/
  guint32 prec;
  gint32 l;
  mpq_t q;
  mpfr_t x;
  mpfr_t jl[2];
  mpfr_t temp;
};

void ncm_sf_mp_spherical_bessel (gulong l, mpq_t q, mpfr_ptr res, mp_rnd_t rnd);
void ncm_sf_mp_spherical_bessel_d (gulong l, gdouble x, mpfr_ptr res, mp_rnd_t rnd);

NcmSfMpSpherBesselRecur *ncm_sf_mp_spherical_bessel_recur_new (gulong prec);
NcmSfMpSpherBesselRecur *ncm_sf_mp_spherical_bessel_recur_read (FILE *f);
void ncm_sf_mp_spherical_bessel_recur_free (NcmSfMpSpherBesselRecur *jlrec);
void ncm_sf_mp_spherical_bessel_recur_set_q (NcmSfMpSpherBesselRecur *jlrec, glong l, mpq_ptr q, mp_rnd_t rnd);
void ncm_sf_mp_spherical_bessel_recur_set_d (NcmSfMpSpherBesselRecur *jlrec, glong l, gdouble x, mp_rnd_t rnd);
void ncm_sf_mp_spherical_bessel_recur_next (NcmSfMpSpherBesselRecur *jlrec, mp_rnd_t rnd);
void ncm_sf_mp_spherical_bessel_recur_previous (NcmSfMpSpherBesselRecur *jlrec, mp_rnd_t rnd);
void ncm_sf_mp_spherical_bessel_recur_goto (NcmSfMpSpherBesselRecur *jlrec, glong l, mp_rnd_t rnd);
void ncm_sf_mp_spherical_bessel_recur_write (NcmSfMpSpherBesselRecur *jlrec, FILE *f);

G_END_DECLS

#endif /* _NC_SF_MP_SPHERICAL_BESSEL_H */
