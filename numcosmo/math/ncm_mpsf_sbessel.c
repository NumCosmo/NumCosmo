/***************************************************************************
 *            ncm_mpsf_sbessel.c
 *
 *  Mon Nov 23 10:04:27 2009
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

/**
 * SECTION:ncm_mpsf_sbessel
 * @title: Multiple Precision Spherical Bessel j_l
 * @short_description: Spherical bessel implementation with support for multiple precision calculation
 *
 * FIXME
 */


#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>


#include <string.h>
#include <stdio.h>
#include <glib.h>
#include <gmp.h>
#include <mpfr.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_bessel.h>

/**
 * ncm_mpsf_sbessel_recur_new: (skip)
 * @prec: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMpsfSBesselRecur *
ncm_mpsf_sbessel_recur_new (gulong prec)
{
  NcmMpsfSBesselRecur *jlrec = g_slice_new (NcmMpsfSBesselRecur);
  jlrec->prec = prec;

  mpq_init (jlrec->q);
  mpfr_inits2 (prec, jlrec->x, jlrec->jl[0], jlrec->jl[1], jlrec->temp, NULL);

  return jlrec;
}

/**
 * ncm_mpsf_sbessel_recur_set_q: (skip)
 * @jlrec: a #NcmMpsfSBesselRecur
 * @l: FIXME
 * @q: FIXME
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_recur_set_q (NcmMpsfSBesselRecur *jlrec, glong l, mpq_ptr q, mp_rnd_t rnd)
{
  jlrec->l = l;

  mpq_set (jlrec->q, q);
  mpfr_set_q (jlrec->x, q, rnd);

  ncm_mpsf_sbessel (l + 0, q, jlrec->jl[0], rnd);
  ncm_mpsf_sbessel (l + 1, q, jlrec->jl[1], rnd);
}

/**
 * ncm_mpsf_sbessel_recur_set_d: (skip)
 * @jlrec: a #NcmMpsfSBesselRecur
 * @l: FIXME
 * @x: FIXME
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_recur_set_d (NcmMpsfSBesselRecur *jlrec, glong l, gdouble x, mp_rnd_t rnd)
{
  jlrec->l = l;

  ncm_rational_coarce_double (x, jlrec->q);
  mpfr_set_q (jlrec->x, jlrec->q, rnd);

  ncm_mpsf_sbessel (l + 0, jlrec->q, jlrec->jl[0], rnd);
  ncm_mpsf_sbessel (l + 1, jlrec->q, jlrec->jl[1], rnd);
}

/**
 * ncm_mpsf_sbessel_recur_free:
 * @jlrec: a #NcmMpsfSBesselRecur
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_recur_free (NcmMpsfSBesselRecur *jlrec)
{
  mpfr_clears (jlrec->x, jlrec->jl[0], jlrec->jl[1], jlrec->temp, NULL);
  mpq_clear (jlrec->q);

  g_slice_free (NcmMpsfSBesselRecur, jlrec);
}

/**
 * ncm_mpsf_sbessel_recur_next: (skip)
 * @jlrec: a #NcmMpsfSBesselRecur
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_recur_next (NcmMpsfSBesselRecur *jlrec, mp_rnd_t rnd)
{
  if (mpfr_sgn (jlrec->x) != 0)
  {
    mpfr_mul_ui (jlrec->temp, jlrec->jl[1], 2 * jlrec->l + 3, rnd);
    mpfr_div (jlrec->temp, jlrec->temp, jlrec->x, rnd);
    mpfr_sub (jlrec->temp, jlrec->temp, jlrec->jl[0], rnd);

    mpfr_swap (jlrec->jl[0], jlrec->jl[1]);
    mpfr_set (jlrec->jl[1], jlrec->temp, rnd);
  }
  jlrec->l++;
}

/**
 * ncm_mpsf_sbessel_recur_previous: (skip)
 * @jlrec: a #NcmMpsfSBesselRecur
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_recur_previous (NcmMpsfSBesselRecur *jlrec, mp_rnd_t rnd)
{
  if (mpfr_sgn (jlrec->x) != 0)
  {
    mpfr_mul_ui (jlrec->temp, jlrec->jl[0], 2 * jlrec->l + 1, rnd);
    mpfr_div (jlrec->temp, jlrec->temp, jlrec->x, rnd);
    mpfr_sub (jlrec->temp, jlrec->temp, jlrec->jl[1], rnd);

    mpfr_swap (jlrec->jl[0], jlrec->jl[1]);
    mpfr_set (jlrec->jl[0], jlrec->temp, rnd);
  }

  jlrec->l--;
}

/**
 * ncm_mpsf_sbessel_recur_goto: (skip)
 * @jlrec: a #NcmMpsfSBesselRecur
 * @l: FIXME
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_recur_goto (NcmMpsfSBesselRecur *jlrec, glong l, mp_rnd_t rnd)
{
  glong sign = GSL_SIGN (l - jlrec->l);
  glong sub = labs(l - jlrec->l);
  glong i;
  if (sub == 0)
    return;
  if (sign == 1)
    for (i = 0; i < sub; i++)
      ncm_mpsf_sbessel_recur_next (jlrec, rnd);
  else
    for (i = 0; i < sub; i++)
      ncm_mpsf_sbessel_recur_previous (jlrec, rnd);
}

/**
 * ncm_mpsf_sbessel_recur_write:
 * @jlrec: a #NcmMpsfSBesselRecur
 * @f: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_recur_write (NcmMpsfSBesselRecur *jlrec, FILE *f)
{
  NCM_WRITE_UINT32(f, jlrec->prec);
  NCM_WRITE_INT32(f, jlrec->l);
  ncm_mpq_out_raw (f, jlrec->q);
  ncm_mpfr_out_raw (f, jlrec->x);
  ncm_mpfr_out_raw (f, jlrec->jl[0]);
  ncm_mpfr_out_raw (f, jlrec->jl[1]);
}

/**
 * ncm_mpsf_sbessel_recur_read: (skip)
 * @f: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmMpsfSBesselRecur *
ncm_mpsf_sbessel_recur_read (FILE *f)
{
  guint32 prec;
  NcmMpsfSBesselRecur *jlrec;
  NCM_READ_UINT32(f, prec);
  jlrec = ncm_mpsf_sbessel_recur_new (prec);
  NCM_READ_INT32(f, jlrec->l);
  ncm_mpq_inp_raw (jlrec->q, f);
  ncm_mpfr_inp_raw (jlrec->x, f);
  ncm_mpfr_inp_raw (jlrec->jl[0], f);
  ncm_mpfr_inp_raw (jlrec->jl[1], f);

  return jlrec;
}

typedef struct __binsplit_spherical_bessel
{
  gulong l;
  glong sincos;
  mpq_t mq2_2;
  mpfr_t sin;
  mpfr_t cos;
} _binsplit_spherical_bessel;

static gpointer
_besselj_bs_alloc (void)
{
  _binsplit_spherical_bessel *bs_data = g_slice_new (_binsplit_spherical_bessel);
  mpq_init (bs_data->mq2_2);
  mpfr_init (bs_data->sin);
  mpfr_init (bs_data->cos);
  return ncm_binsplit_alloc ((gpointer)bs_data);
}

static void
_besselj_bs_free (gpointer p)
{
  NcmBinSplit *bs = (NcmBinSplit *)p;
  _binsplit_spherical_bessel *bs_data = (_binsplit_spherical_bessel *)bs->userdata;
  mpq_clear (bs_data->mq2_2);
  mpfr_clear (bs_data->sin);
  mpfr_clear (bs_data->cos);
  g_slice_free (_binsplit_spherical_bessel, bs_data);
  /* Leak we dont have a free function for binsplit FIXME:LEAK */
}

NcmBinSplit **
_ncm_mpsf_sbessel_get_bs ()
{
  static GStaticMutex create_lock = G_STATIC_MUTEX_INIT;
  static NcmMemoryPool *mp = NULL;

  g_static_mutex_lock (&create_lock);
  if (mp == NULL)
    mp = ncm_memory_pool_new (_besselj_bs_alloc, _besselj_bs_free);
  g_static_mutex_unlock (&create_lock);

  return ncm_memory_pool_get (mp);
}

#define NC_BINSPLIT_EVAL_NAME binsplit_spherical_bessel_taylor
#define _mq2_2 (((_binsplit_spherical_bessel *)data)->mq2_2)
#define _l (((_binsplit_spherical_bessel *)data)->l)

NCM_BINSPLIT_DECL(binsplit_spherical_bessel_taylor_p,v,u,n,data)
{
  if (n == 0)
    mpz_set (v, u);
  else
    mpz_mul (v, u, mpq_numref (_mq2_2));
}
#define _BINSPLIT_FUNC_P binsplit_spherical_bessel_taylor_p

NCM_BINSPLIT_DECL(binsplit_spherical_bessel_taylor_q,v,u,n,data)
{
  if (n == 0)
    mpz_set (v, u);
  else
  {
    mpz_mul_ui (v, u, n * (2L * (n + _l) + 1L));
    mpz_mul (v, v, mpq_denref(_mq2_2));
  }
}
#define _BINSPLIT_FUNC_Q binsplit_spherical_bessel_taylor_q

#define _BINSPLIT_FUNC_B NCM_BINSPLIT_DENC_NULL
#define _BINSPLIT_FUNC_A NCM_BINSPLIT_DENC_NULL

#include "binsplit_eval.c"
#undef _mq2_2

/* Assymptotic expansion 4F1 */

#define NC_BINSPLIT_EVAL_NAME binsplit_spherical_bessel_assympt
#define _mq2_2 (((_binsplit_spherical_bessel *)data)->mq2_2)
#define _l (((_binsplit_spherical_bessel *)data)->l)
#define _sincos (((_binsplit_spherical_bessel *)data)->sincos)

NCM_BINSPLIT_DECL(binsplit_spherical_bessel_assympt_p,v,u,n,data)
{
  if (n == 0)
    mpz_set (v, u);
  else
  {
    mpz_mul_ui (v, u, (_sincos - _l + 2L * n - 2L) * (_sincos - _l + 2L* n - 1L));
    mpz_mul_ui (v, v, (_sincos + 1L + _l + 2L * n - 2L) * (_sincos + 1L + _l + 2L* n - 1L));
    mpz_mul (v, v, mpq_numref (_mq2_2));
  }
}
#define _BINSPLIT_FUNC_P binsplit_spherical_bessel_assympt_p

NCM_BINSPLIT_DECL(binsplit_spherical_bessel_assympt_q,v,u,n,data)
{
  if (n == 0)
    mpz_set (v, u);
  else
  {
    mpz_mul_ui (v, u, (_sincos + 1L + 2L * n - 2L) * (_sincos + 1L + 2L* n - 1L));
    mpz_mul (v, v, mpq_denref(_mq2_2));
  }
}
#define _BINSPLIT_FUNC_Q binsplit_spherical_bessel_assympt_q

#define _BINSPLIT_FUNC_B NCM_BINSPLIT_DENC_NULL
#define _BINSPLIT_FUNC_A NCM_BINSPLIT_DENC_NULL

#include "binsplit_eval.c"
#undef _mq2_2

static void
_taylor_mpfr (gulong l, mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  NcmBinSplit **bs_ptr = _ncm_mpsf_sbessel_get_bs ();
  NcmBinSplit *bs = *bs_ptr;
  _binsplit_spherical_bessel *data = (_binsplit_spherical_bessel *) bs->userdata;
  gulong n;

  data->l = l;
  mpq_mul (data->mq2_2, q, q);
  mpq_neg (data->mq2_2, data->mq2_2);
  mpq_div_2exp (data->mq2_2, data->mq2_2, 1);

  //mpfr_printf ("# Taylor %ld %Qd | %Qd\n", l, q, data->mq2_2);

  ncm_binsplit_eval_prec (bs, binsplit_spherical_bessel_taylor, 10, mpfr_get_prec (res));

  //mpfr_printf ("# Taylor %ld %Qd | %Zd %Zd\n", l, q, bs->T, bs->Q);

  mpfr_set_q (res, q, rnd);
  mpfr_pow_ui (res, res, l, rnd);
  mpfr_mul_z (res, res, bs->T, rnd);
  mpfr_div_z (res, res, bs->Q, rnd);

  for (n = 1; n <= l; n++)
    mpfr_div_ui (res, res, 2L * n + 1, rnd);

  ncm_memory_pool_return (bs_ptr);
  return;
}

static void
_assympt_mpfr (gulong l, mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  NcmBinSplit **bs_ptr = _ncm_mpsf_sbessel_get_bs ();
  NcmBinSplit *bs = *bs_ptr;
  _binsplit_spherical_bessel *data = (_binsplit_spherical_bessel *) bs->userdata;
  gulong prec = mpfr_get_prec (res);
#define sin_x data->sin
#define cos_x data->cos
  mpfr_set_prec (sin_x, prec);
  mpfr_set_prec (cos_x, prec);

  mpfr_set_q (res, q, rnd);
  mpfr_sin_cos (sin_x, cos_x, res, rnd);

  switch (l % 4)
  {
    case 0:
      break;
    case 1:
      mpfr_swap (sin_x, cos_x);
      mpfr_neg (sin_x, sin_x, rnd);
      break;
    case 2:
      mpfr_neg (sin_x, sin_x, rnd);
      mpfr_neg (cos_x, cos_x, rnd);
      break;
    case 3:
      mpfr_swap (sin_x, cos_x);
      mpfr_neg (cos_x, cos_x, rnd);
      break;
  }

  if (l > 0)
  {
    mpfr_mul_ui (cos_x, cos_x, l * (l + 1), rnd);
    mpfr_div (cos_x, cos_x, res, rnd);
    mpfr_div (cos_x, cos_x, res, rnd);
    mpfr_div_2ui (cos_x, cos_x, 1, rnd);
  }

  mpfr_div (sin_x, sin_x, res, rnd);

  data->l = l;
  mpq_inv (data->mq2_2, q);
  mpq_mul (data->mq2_2, data->mq2_2, data->mq2_2);
  mpq_neg (data->mq2_2, data->mq2_2);
  mpq_div_2exp (data->mq2_2, data->mq2_2, 2);

  data->sincos = 0;
  binsplit_spherical_bessel_assympt (bs, 0, (l + 1) / 2 + (l + 1) % 2);
  mpfr_mul_z (sin_x, sin_x, bs->T, rnd);
  mpfr_div_z (sin_x, sin_x, bs->Q, rnd);

  data->sincos = 1;
  if (l > 0)
  {
    binsplit_spherical_bessel_assympt (bs, 0, l / 2 + l % 2);
    mpfr_mul_z (cos_x, cos_x, bs->T, rnd);
    mpfr_div_z (cos_x, cos_x, bs->Q, rnd);
    mpfr_add (res, sin_x, cos_x, rnd);
  }
  else
    mpfr_set (res, sin_x, rnd);

  ncm_memory_pool_return (bs_ptr);
  return;
}

/**
 * ncm_mpsf_sbessel: (skip)
 * @l: FIXME
 * @q: FIXME
 * @res: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel (gulong l, mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  gdouble x = mpq_get_d (q);
  if (x == 0)
  {
    if (l == 0)
      mpfr_set_ui (res, 1, rnd);
    else
      mpfr_set_ui (res, 0, rnd);
    return;
  }

  if (fabs(x) < l)
    _taylor_mpfr (l, q, res, rnd);
  else
    _assympt_mpfr (l, q, res, rnd);
}

/**
 * ncm_mpsf_sbessel_d: (skip)
 * @l: FIXME
 * @x: FIXME
 * @res: FIXME
 * @rnd: FIXME
 *
 * FIXME
 */
void
ncm_mpsf_sbessel_d (gulong l, gdouble x, mpfr_ptr res, mp_rnd_t rnd)
{
  mpq_t q;
  mpq_init (q);
  ncm_rational_coarce_double (x, q);
  ncm_mpsf_sbessel (l, q, res, rnd);
  mpq_clear (q);
}
