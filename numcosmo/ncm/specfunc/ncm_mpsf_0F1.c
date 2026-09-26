/***************************************************************************
 *            ncm_mpsf_0F1.c
 *
 *  Thu October 13 22:16:28 2011
 *  Copyright  2011  Sandro Dias Pinto Vitenti
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
 * NcmMpsf0F1:
 *
 * Arbitrary-precision confluent hypergeometric limit function ${}_0F_1$.
 *
 * Computes ${}_0F_1(b; x) = \sum_{n=0}^\infty x^n / [(b)_n\, n!]$ from its Taylor series
 * by binary splitting, see #NcmBinSplit, summing until the next term is below the
 * precision of the result. The splitting buffers come from a thread-safe pool, released
 * by ncm_mpsf_0F1_free_cache().
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/specfunc/ncm_mpsf_0F1.h"
#include "ncm/specfunc/ncm_binsplit.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_util.h"

typedef struct __binsplit_0F1
{
  mpq_t b;
  mpz_t xn_bd;
  mpz_t xd;
  mpz_t tmp;
} _binsplit_0F1;

static gpointer
_besselj_bs_alloc (gpointer userdata)
{
  _binsplit_0F1 *bs_data = g_slice_new (_binsplit_0F1);

  NCM_UNUSED (userdata);

  mpq_init (bs_data->b);
  mpz_init (bs_data->xn_bd);
  mpz_init (bs_data->xd);
  mpz_init (bs_data->tmp);

  return ncm_binsplit_alloc ((gpointer) bs_data);
}

static void
_besselj_bs_free (gpointer p)
{
  NcmBinSplit *bs        = (NcmBinSplit *) p;
  _binsplit_0F1 *bs_data = (_binsplit_0F1 *) bs->userdata;

  mpq_clear (bs_data->b);
  mpz_clear (bs_data->xn_bd);
  mpz_clear (bs_data->xd);
  mpz_clear (bs_data->tmp);
  g_slice_free (_binsplit_0F1, bs_data);
  ncm_binsplit_free (bs);
}

static NcmMemoryPool *__mp = NULL;

G_LOCK_DEFINE_STATIC (__create_lock);

/**
 * _ncm_mpsf_0F1_get_bs: (skip)
 *
 * Takes a #NcmBinSplit for the ${}_0F_1$ series from the pool; it must be given back
 * with ncm_memory_pool_return().
 *
 * Returns: a pointer to the #NcmBinSplit pointer.
 */
NcmBinSplit **
_ncm_mpsf_0F1_get_bs (void)
{
  G_LOCK (__create_lock);

  if (__mp == NULL)
    __mp = ncm_memory_pool_new (_besselj_bs_alloc, NULL, _besselj_bs_free);

  G_UNLOCK (__create_lock);

  return ncm_memory_pool_get (__mp);
}

#define NC_BINSPLIT_EVAL_NAME binsplit_0F1_taylor
#define _xn_bd (((_binsplit_0F1 *) data)->xn_bd)
#define _xd (((_binsplit_0F1 *) data)->xd)
#define _b (((_binsplit_0F1 *) data)->b)
#define _tmp (((_binsplit_0F1 *) data)->tmp)

NCM_BINSPLIT_DECL (binsplit_0F1_taylor_p, v, u, n, data)
{
  if (n == 0)
    mpz_set (v, u);
  else
    mpz_mul (v, u, _xn_bd);
}
#define _BINSPLIT_FUNC_P binsplit_0F1_taylor_p

NCM_BINSPLIT_DECL (binsplit_0F1_taylor_q, v, u, n, data)
{
  if (n == 0)
  {
    mpz_set (v, u);
  }
  else
  {
    mpz_set (_tmp, mpq_numref (_b));
    mpz_addmul_ui (_tmp, mpq_denref (_b), n - 1L);
    mpz_mul (v, u, _tmp);
    mpz_mul_ui (v, v, n);
    mpz_mul (v, v, _xd);
  }
}
#define _BINSPLIT_FUNC_Q binsplit_0F1_taylor_q

#define _BINSPLIT_FUNC_B NCM_BINSPLIT_DENC_NULL
#define _BINSPLIT_FUNC_A NCM_BINSPLIT_DENC_NULL

#include "ncm/specfunc/ncm_binsplit_eval.c"
#undef _xn_bd
#undef _xd
#undef _b
#undef _tmp

static void
_taylor_0F1 (mpq_t b, mpq_t x, mpfr_ptr res, mp_rnd_t rnd)
{
  NcmBinSplit **bs_ptr = _ncm_mpsf_0F1_get_bs ();
  NcmBinSplit *bs      = *bs_ptr;
  _binsplit_0F1 *data  = (_binsplit_0F1 *) bs->userdata;

  mpq_set (data->b, b);
  mpz_set (data->xd, mpq_denref (x));
  mpz_mul (data->xn_bd, mpq_numref (x), mpq_denref (b));

  ncm_binsplit_eval_prec (bs, binsplit_0F1_taylor, 10, mpfr_get_prec (res));

  mpfr_set_z (res, bs->T, rnd);
  mpfr_div_z (res, res, bs->Q, rnd);

  ncm_memory_pool_return (bs_ptr);

  return;
}

/**
 * ncm_mpsf_0F1_q: (skip)
 * @b: the parameter $b$
 * @q: the argument $x$
 * @res: the output, at its own precision
 * @rnd: the rounding mode
 *
 * Computes ${}_0F_1(b; x)$ into @res.
 */
void
ncm_mpsf_0F1_q (mpq_t b, mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  _taylor_0F1 (b, q, res, rnd);
}

/**
 * ncm_mpsf_0F1_d: (skip)
 * @b: the parameter $b$
 * @x: the argument $x$
 * @res: the output, at its own precision
 * @rnd: the rounding mode
 *
 * Same as ncm_mpsf_0F1_q() with @b and @x converted to rationals that agree with them to
 * $10^{-15}$, see ncm_rational_coarse_double().
 */
void
ncm_mpsf_0F1_d (gdouble b, gdouble x, mpfr_ptr res, mp_rnd_t rnd)
{
  mpq_t xq, bq;

  mpq_init (xq);
  mpq_init (bq);
  ncm_rational_coarse_double (b, bq);
  ncm_rational_coarse_double (x, xq);
  ncm_mpsf_0F1_q (bq, xq, res, rnd);
  mpq_clear (xq);
  mpq_clear (bq);
}

/**
 * ncm_sf_0F1: (skip)
 * @b: the parameter $b$
 * @x: the argument $x$
 *
 * Same as ncm_mpsf_0F1_d() at 53 bits, rounded to nearest.
 *
 * Returns: ${}_0F_1(b; x)$.
 */
gdouble
ncm_sf_0F1 (gdouble b, gdouble x)
{
  MPFR_DECL_INIT (res, 53);

  gdouble res_d;

  ncm_mpsf_0F1_d (b, x, res, GMP_RNDN);
  res_d = mpfr_get_d (res, GMP_RNDN);

  return res_d;
}

/**
 * ncm_mpsf_0F1_free_cache:
 *
 * Frees the pool of splitting buffers of the ${}_0F_1$ functions.
 */
void
ncm_mpsf_0F1_free_cache (void)
{
  G_LOCK (__create_lock);

  if (__mp != NULL)
  {
    ncm_memory_pool_free (__mp, TRUE);
    __mp = NULL;
  }

  G_UNLOCK (__create_lock);
}

