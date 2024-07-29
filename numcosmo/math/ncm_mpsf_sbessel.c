/***************************************************************************
 *            ncm_mpsf_sbessel.c
 *
 *  Mon Nov 23 10:04:27 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
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

/**
 * SECTION:ncm_mpsf_sbessel
 * @title: NcmMpsfSBessel
 * @short_description: Multiple precision spherical bessel implementation.
 *
 * Implementation of multiple precision spherical Bessel functions using the GNU MPFR
 * library. This module utilizes binary splitting to compute the functions, employing
 * both the Taylor series and asymptotic expansion methods. It ensures high precision,
 * making it suitable for accurate computations in scenarios involving spherical Bessel
 * functions.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mpsf_sbessel.h"
#include "math/ncm_util.h"
#include "math/ncm_cfg.h"
#include "math/ncm_binsplit.h"
#include "math/ncm_memory_pool.h"

typedef struct __binsplit_spherical_bessel
{
  gulong l;
  glong sincos;
  mpq_t mq2_2;
  mpfr_t sin;
  mpfr_t cos;
} _binsplit_spherical_bessel;

static gpointer
_besselj_bs_alloc (gpointer userdata)
{
  _binsplit_spherical_bessel *bs_data = g_slice_new (_binsplit_spherical_bessel);

  NCM_UNUSED (userdata);

  mpq_init (bs_data->mq2_2);
  mpfr_init (bs_data->sin);
  mpfr_init (bs_data->cos);

  return ncm_binsplit_alloc ((gpointer) bs_data);
}

static void
_besselj_bs_free (gpointer p)
{
  NcmBinSplit *bs                     = (NcmBinSplit *) p;
  _binsplit_spherical_bessel *bs_data = (_binsplit_spherical_bessel *) bs->userdata;

  mpq_clear (bs_data->mq2_2);
  mpfr_clear (bs_data->sin);
  mpfr_clear (bs_data->cos);
  g_slice_free (_binsplit_spherical_bessel, bs_data);
  /* Leak we dont have a free function for binsplit FIXME:LEAK */
}

G_LOCK_DEFINE_STATIC (__create_lock);

static NcmMemoryPool *__mp = NULL;

NcmBinSplit **
_ncm_mpsf_sbessel_get_bs (void)
{
  G_LOCK (__create_lock);

  if (__mp == NULL)
    __mp = ncm_memory_pool_new (_besselj_bs_alloc, NULL, _besselj_bs_free);

  G_UNLOCK (__create_lock);

  return ncm_memory_pool_get (__mp);
}

#define NC_BINSPLIT_EVAL_NAME binsplit_spherical_bessel_taylor
#define _mq2_2 (((_binsplit_spherical_bessel *) data)->mq2_2)
#define _l (((_binsplit_spherical_bessel *) data)->l)

NCM_BINSPLIT_DECL (binsplit_spherical_bessel_taylor_p, v, u, n, data)
{
  if (n == 0)
    mpz_set (v, u);
  else
    mpz_mul (v, u, mpq_numref (_mq2_2));
}
#define _BINSPLIT_FUNC_P binsplit_spherical_bessel_taylor_p

NCM_BINSPLIT_DECL (binsplit_spherical_bessel_taylor_q, v, u, n, data)
{
  if (n == 0)
  {
    mpz_set (v, u);
  }
  else
  {
    mpz_mul_ui (v, u, n * (2L * (n + _l) + 1L));
    mpz_mul (v, v, mpq_denref (_mq2_2));
  }
}
#define _BINSPLIT_FUNC_Q binsplit_spherical_bessel_taylor_q

#define _BINSPLIT_FUNC_B NCM_BINSPLIT_DENC_NULL
#define _BINSPLIT_FUNC_A NCM_BINSPLIT_DENC_NULL

#include "ncm_binsplit_eval.c"
#undef _mq2_2

/* Assymptotic expansion 4F1 */

#define NC_BINSPLIT_EVAL_NAME binsplit_spherical_bessel_assympt
#define _mq2_2 (((_binsplit_spherical_bessel *) data)->mq2_2)
#define _l (((_binsplit_spherical_bessel *) data)->l)
#define _sincos (((_binsplit_spherical_bessel *) data)->sincos)

NCM_BINSPLIT_DECL (binsplit_spherical_bessel_assympt_p, v, u, n, data)
{
  if (n == 0)
  {
    mpz_set (v, u);
  }
  else
  {
    mpz_mul_ui (v, u, (_sincos - _l + 2L * n - 2L) * (_sincos - _l + 2L * n - 1L));
    mpz_mul_ui (v, v, (_sincos + 1L + _l + 2L * n - 2L) * (_sincos + 1L + _l + 2L * n - 1L));
    mpz_mul (v, v, mpq_numref (_mq2_2));
  }
}
#define _BINSPLIT_FUNC_P binsplit_spherical_bessel_assympt_p

NCM_BINSPLIT_DECL (binsplit_spherical_bessel_assympt_q, v, u, n, data)
{
  if (n == 0)
  {
    mpz_set (v, u);
  }
  else
  {
    mpz_mul_ui (v, u, (_sincos + 1L + 2L * n - 2L) * (_sincos + 1L + 2L * n - 1L));
    mpz_mul (v, v, mpq_denref (_mq2_2));
  }
}
#define _BINSPLIT_FUNC_Q binsplit_spherical_bessel_assympt_q

#define _BINSPLIT_FUNC_B NCM_BINSPLIT_DENC_NULL
#define _BINSPLIT_FUNC_A NCM_BINSPLIT_DENC_NULL

#include "ncm_binsplit_eval.c"
#undef _mq2_2

static void
_taylor_mpfr (gulong l, mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  NcmBinSplit **bs_ptr             = _ncm_mpsf_sbessel_get_bs ();
  NcmBinSplit *bs                  = *bs_ptr;
  _binsplit_spherical_bessel *data = (_binsplit_spherical_bessel *) bs->userdata;
  gulong n;

  data->l = l;
  mpq_mul (data->mq2_2, q, q);
  mpq_neg (data->mq2_2, data->mq2_2);
  mpq_div_2exp (data->mq2_2, data->mq2_2, 1);

  /*mpfr_printf ("# Taylor %ld %Qd | %Qd\n", l, q, data->mq2_2); */

  ncm_binsplit_eval_prec (bs, binsplit_spherical_bessel_taylor, 10, mpfr_get_prec (res));

  /*mpfr_printf ("# Taylor %ld %Qd | %Zd %Zd\n", l, q, bs->T, bs->Q); */

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
  NcmBinSplit **bs_ptr             = _ncm_mpsf_sbessel_get_bs ();
  NcmBinSplit *bs                  = *bs_ptr;
  _binsplit_spherical_bessel *data = (_binsplit_spherical_bessel *) bs->userdata;
  gulong prec                      = mpfr_get_prec (res);

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
  {
    mpfr_set (res, sin_x, rnd);
  }

  ncm_memory_pool_return (bs_ptr);

  return;
}

/**
 * ncm_mpsf_sbessel: (skip)
 * @l: $\ell$ Spherical Bessel $j_\ell$ parameters as a rational number $\ell = q_\ell$
 * @q: argument as a rational number $x = q_x$
 * @res: mpfr variable containing the result $j_\ell(x)$
 * @rnd: mpfr rounding mode
 *
 * Computes the Spherical Bessel function $j_\ell(x)$.
 *
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

  if (fabs (x) < l)
    _taylor_mpfr (l, q, res, rnd);
  else
    _assympt_mpfr (l, q, res, rnd);
}

/**
 * ncm_mpsf_sbessel_d: (skip)
 * @l: $\ell$ Spherical Bessel $j_\ell$ parameters
 * @x: function argument $x$
 * @res: mpfr variable containing the result $j_\ell(x)$
 * @rnd: mpfr rounding mode
 *
 * Computes the Spherical Bessel function $j_\ell(x)$.
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

/**
 * ncm_mpsf_sbessel_free_cache:
 *
 * Frees all buffers created to compute
 * ncm_mpsf_sbessel functions.
 *
 */
void
ncm_mpsf_sbessel_free_cache (void)
{
  G_LOCK (__create_lock);

  if (__mp != NULL)
  {
    ncm_memory_pool_free (__mp, TRUE);
    __mp = NULL;
  }

  G_UNLOCK (__create_lock);
}

