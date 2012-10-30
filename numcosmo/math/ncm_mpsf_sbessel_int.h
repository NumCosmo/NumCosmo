/***************************************************************************
 *            ncm_mpsf_sbessel_int.h
 *
 *  Thu Feb 11 23:24:17 2010
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

#ifndef _NCM_MPSF_SBESSEL_INT_H
#define _NCM_MPSF_SBESSEL_INT_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_mpsf_sbessel.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/grid_one.h>

G_BEGIN_DECLS

typedef struct _NcmMpsfSBesselIntegRecur NcmMpsfSBesselIntegRecur;

/**
 * NcmMpsfSBesselIntegRecur:
 *
 * FIXME
 */
struct _NcmMpsfSBesselIntegRecur
{
  /*< private >*/
  NcmMpsfSBesselRecur *jlrec;
  guint32 prec;
  mpq_ptr q;
  mpfr_ptr x;
  mpfr_t int_jl_xn[4];
  mpfr_t int_jlp1_xn[4];
  mpfr_t x_pow_n[4];
  mpfr_t temp1;
  mpfr_t temp2;
};

typedef struct _NcmMpsfSBesselIntSpline NcmMpsfSBesselIntSpline;

/**
 * NcmMpsfSBesselIntSpline:
 *
 * FIXME
 */
struct _NcmMpsfSBesselIntSpline
{
  /*< private >*/
  NcmSpline *s;
  NcmMpsfSBesselIntegRecur **xnjlrec;
  guint32 *map_ij2r;
  guint32 row;
  guint32 rules_count;
  guint32 prec;
  NcmGrid *x;
  NcmGrid *k;
  glong l;
  GHashTable *int_hash;
  mpfr_t rules[4];
  mpfr_t crules[4];
  gboolean prepared;
};

#define NCM_MPSF_SBESSEL_INT_MAP(int_jlspline,xi,kj) ((int_jlspline)->map_ij2r[((int_jlspline)->row * kj + xi)])

void ncm_mpsf_sbessel_jl_xj_integral (gint l, gint j, gdouble x, mpfr_t res, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integral_q (gint l, gint j, mpq_t q, mpfr_t res, mp_rnd_t rnd);
gdouble ncm_mpsf_sbessel_integrate (NcmMpsfSBesselIntSpline *int_jlspline, NcmSpline *s, gint l, guint ki, guint xi, gint diff);

NcmMpsfSBesselIntegRecur *ncm_mpsf_sbessel_jl_xj_integral_recur_new (gulong prec, NcmMpsfSBesselRecur *jlrec);
NcmMpsfSBesselIntegRecur *ncm_mpsf_sbessel_jl_xj_integral_recur_read (FILE *f);
void ncm_mpsf_sbessel_jl_xj_integral_recur_free (NcmMpsfSBesselIntegRecur *xnjlrec);
void ncm_mpsf_sbessel_jl_xj_integral_recur_next (NcmMpsfSBesselIntegRecur *xnjlrec, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integral_recur_previous (NcmMpsfSBesselIntegRecur *xnjlrec, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integral_recur_goto (NcmMpsfSBesselIntegRecur *xnjlrec, glong l, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integral_recur_calc_djl (NcmMpsfSBesselIntegRecur *xnjlrec, glong n, mpfr_ptr rule, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integral_recur_calc_d2jl (NcmMpsfSBesselIntegRecur *xnjlrec, glong n, mpfr_ptr rule, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integral_recur_write (NcmMpsfSBesselIntegRecur *xnjlrec, FILE *f);
void ncm_mpsf_sbessel_jl_xj_integral_a_b (NcmMpsfSBesselIntSpline *int_jlspline, guint ki, guint xi, mp_rnd_t rnd);

NcmMpsfSBesselIntSpline *ncm_mpsf_sbessel_jl_xj_integrate_spline_new (gulong prec, NcmGrid *x, NcmGrid *k);
NcmMpsfSBesselIntSpline *ncm_mpsf_sbessel_jl_xj_integrate_spline_new_from_sections (gulong prec, NcmGridSection *x_secs, NcmGridSection *k_secs);
NcmMpsfSBesselIntSpline *ncm_mpsf_sbessel_jl_xj_integrate_spline_cached_new (gulong prec, glong l, NcmGridSection *x_secs, NcmGridSection *k_secs, mp_rnd_t rnd);
NcmMpsfSBesselIntSpline *ncm_mpsf_sbessel_jl_xj_integrate_spline_load (gchar *filename);

void ncm_mpsf_sbessel_jl_xj_integrate_spline_prepare (NcmMpsfSBesselIntSpline *int_jlspline, glong l, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integrate_spline_prepare_new (NcmMpsfSBesselIntSpline *int_jlspline, glong l, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integrate_spline_save (NcmMpsfSBesselIntSpline *int_jlspline, gchar *filename);
void ncm_mpsf_sbessel_jl_xj_integrate_spline_next (NcmMpsfSBesselIntSpline *int_jlspline, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integrate_spline_previous (NcmMpsfSBesselIntSpline *int_jlspline, mp_rnd_t rnd);
void ncm_mpsf_sbessel_jl_xj_integrate_spline_goto (NcmMpsfSBesselIntSpline *int_jlspline, gulong l, mp_rnd_t rnd);

G_END_DECLS

#endif /* _NCM_MPSF_SBESSEL_INT_H */
