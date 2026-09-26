/***************************************************************************
 *            ncm_mpsf_trig_int.c
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

/**
 * NcmMpsfTrigInt:
 *
 * Arbitrary-precision sine integral $\mathrm{Si}(x) = \int_0^x \sin t / t\,\mathrm{d}t$.
 *
 * Sums, by binary splitting (see #NcmBinSplit), either the Taylor series or, for large
 * $x > 0$, the asymptotic series of $\mathrm{Si}(x) = \pi/2 - f(x)\cos x - g(x)\sin x$.
 * The asymptotic series diverges; it is used only when its smallest term, near the
 * index $n$ with $2n \approx x$, is below the precision of the result. For $x \ge
 * 2^{p}$, with $p$ the precision in bits, the result is $\pi/2$, and negative arguments
 * use $\mathrm{Si}(-x) = -\mathrm{Si}(x)$. The splitting buffers come from thread-safe
 * pools, released by ncm_mpsf_sin_int_free_cache().
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/specfunc/ncm_mpsf_trig_int.h"
#include "ncm/specfunc/ncm_binsplit.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"
#include "ncm/core/ncm_memory_pool.h"

#define NC_BINSPLIT_EVAL_NAME binsplit_sin_integral_taylor
#define _mx2 (((mpq_ptr) data))

NCM_BINSPLIT_DECL (binsplit_sin_integral_taylor_p, v, u, n, data)
{
  if (n == 0)
    mpz_set (v, u);
  else
    mpz_mul (v, u, mpq_numref (_mx2));
}
#define _BINSPLIT_FUNC_P binsplit_sin_integral_taylor_p

NCM_BINSPLIT_DECL (binsplit_sin_integral_taylor_q, v, u, n, data)
{
  if (n == 0)
  {
    mpz_set (v, u);
  }
  else
  {
    mpz_mul_ui (v, u, n * (2L * n + 1L));
    mpz_mul (v, v, mpq_denref (_mx2));
  }
}
#define _BINSPLIT_FUNC_Q binsplit_sin_integral_taylor_q

NCM_BINSPLIT_DECL (binsplit_sin_integral_taylor_b, v, u, n, data)
{
  NCM_UNUSED (data);
  mpz_mul_ui (v, u, (2L * n + 1L));
}
#define _BINSPLIT_FUNC_B binsplit_sin_integral_taylor_b
#define _HAS_FUNC_B

#define _BINSPLIT_FUNC_A NCM_BINSPLIT_DENC_NULL

#include "ncm/specfunc/ncm_binsplit_eval.c"
#undef _mx2

#define NC_BINSPLIT_EVAL_NAME binsplit_sin_integral_assym
/* -x^2/2, despite the name */
#define _m2x2 (((mpq_ptr *) data)[0])
#define _sincos (GPOINTER_TO_INT (((gpointer *) data)[1]))

NCM_BINSPLIT_DECL (binsplit_sin_integral_assym_p, v, u, n, data)
{
  if (n == 0)
  {
    mpz_set (v, u);
  }
  else
  {
    mpz_mul (v, u, mpq_denref (_m2x2));
    mpz_mul_ui (v, v, n * (2L * n + _sincos));
  }
}
#define _BINSPLIT_FUNC_P binsplit_sin_integral_assym_p

NCM_BINSPLIT_DECL (binsplit_sin_integral_assym_q, v, u, n, data)
{
  if (n == 0)
    mpz_set (v, u);
  else
    mpz_mul (v, u, mpq_numref (_m2x2));
}
#define _BINSPLIT_FUNC_Q binsplit_sin_integral_assym_q

#define _BINSPLIT_FUNC_B NCM_BINSPLIT_DENC_NULL
#define _BINSPLIT_FUNC_A NCM_BINSPLIT_DENC_NULL

#include "ncm/specfunc/ncm_binsplit_eval.c"
#undef _m2x2

/* Pools of splitting buffers: the Taylor series takes -x^2/2 as user data, the
 * asymptotic one an array {-x^2/2, sincos} */

static gpointer
_taylor_bs_alloc (gpointer userdata)
{
  mpq_ptr mq2 = g_slice_new (__mpq_struct);

  NCM_UNUSED (userdata);
  mpq_init (mq2);

  return ncm_binsplit_alloc (mq2);
}

static void
_taylor_bs_free (gpointer p)
{
  NcmBinSplit *bs = (NcmBinSplit *) p;
  mpq_ptr mq2     = (mpq_ptr) bs->userdata;

  mpq_clear (mq2);
  g_slice_free (__mpq_struct, mq2);
  ncm_binsplit_free (bs);
}

static gpointer
_assym_bs_alloc (gpointer userdata)
{
  gpointer *data = g_new0 (gpointer, 2);

  NCM_UNUSED (userdata);

  return ncm_binsplit_alloc (data);
}

static void
_assym_bs_free (gpointer p)
{
  NcmBinSplit *bs = (NcmBinSplit *) p;

  g_free (bs->userdata);
  ncm_binsplit_free (bs);
}

G_LOCK_DEFINE_STATIC (__create_lock);

static NcmMemoryPool *__mp_taylor = NULL;
static NcmMemoryPool *__mp_assym  = NULL;

static NcmBinSplit **
_ncm_mpsf_trig_int_get_bs (gboolean taylor)
{
  NcmMemoryPool *mp;

  G_LOCK (__create_lock);

  if (__mp_taylor == NULL)
  {
    __mp_taylor = ncm_memory_pool_new (_taylor_bs_alloc, NULL, _taylor_bs_free);
    __mp_assym  = ncm_memory_pool_new (_assym_bs_alloc, NULL, _assym_bs_free);
  }

  mp = taylor ? __mp_taylor : __mp_assym;

  G_UNLOCK (__create_lock);

  return ncm_memory_pool_get (mp);
}

static void
_taylor_mpfr (mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  NcmBinSplit **bs_ptr = _ncm_mpsf_trig_int_get_bs (TRUE);
  NcmBinSplit *bs      = *bs_ptr;
  mpq_ptr mq2          = (mpq_ptr) bs->userdata;

  mpq_mul (mq2, q, q);
  mpq_neg (mq2, mq2);
  mpq_div_2exp (mq2, mq2, 1);

  ncm_binsplit_eval_prec (bs, binsplit_sin_integral_taylor, 10, mpfr_get_prec (res));
  ncm_binsplit_get (bs, res);
  mpfr_mul_q (res, res, q, rnd);

  ncm_memory_pool_return (bs_ptr);
}

static void
_assym_mpfr (mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  glong prec = mpfr_get_prec (res);
  glong mprec;

  mpfr_set_q (res, q, rnd);
  mprec = res->_mpfr_exp;

  if ((prec - mprec) > 0)
  {
    NcmBinSplit **bs_ptr = _ncm_mpsf_trig_int_get_bs (FALSE);
    NcmBinSplit *bs      = *bs_ptr;
    gpointer *data       = (gpointer *) bs->userdata;
    gulong nf;
    mpq_t mq2_2;

    mpq_init (mq2_2);
    mpq_mul (mq2_2, q, q);
    mpq_neg (mq2_2, mq2_2);
    mpq_div_2exp (mq2_2, mq2_2, 1);

    MPFR_DECL_INIT (sin_x, prec);
    MPFR_DECL_INIT (cos_x, prec);

    mpfr_sin_cos (sin_x, cos_x, res, rnd);
    mpfr_div_q (sin_x, sin_x, q, rnd);
    mpfr_div_q (sin_x, sin_x, q, rnd);
    mpfr_div_q (cos_x, cos_x, q, rnd);

    nf = ceil (fabs (prec * M_LN2 / (log (fabs (mpq_get_d (mq2_2))))));

    if (nf == 0)
      nf = 4;

    data[0] = mq2_2;
    data[1] = GINT_TO_POINTER (-1);
    ncm_binsplit_eval_prec (bs, binsplit_sin_integral_assym, nf, prec - mprec);
    mpfr_mul_z (cos_x, cos_x, bs->T, rnd);
    mpfr_div_z (cos_x, cos_x, bs->Q, rnd);

    data[1] = GINT_TO_POINTER (1);
    ncm_binsplit_eval_prec (bs, binsplit_sin_integral_assym, nf, prec - mprec);
    mpfr_mul_z (sin_x, sin_x, bs->T, rnd);
    mpfr_div_z (sin_x, sin_x, bs->Q, rnd);

    mpfr_const_pi (res, rnd);
    mpfr_div_2ui (res, res, 1, rnd);
    mpfr_sub (res, res, sin_x, rnd);
    mpfr_sub (res, res, cos_x, rnd);

    data[0] = NULL;
    mpq_clear (mq2_2);
    ncm_memory_pool_return (bs_ptr);
  }
  else
  {
    mpfr_const_pi (res, rnd);
    mpfr_div_2ui (res, res, 1, rnd);
  }
}

/**
 * ncm_mpsf_sin_int_mpfr: (skip)
 * @q: the argument $x$
 * @res: the output, at its own precision
 * @rnd: the rounding mode
 *
 * Computes $\mathrm{Si}(x)$ into @res, see #NcmMpsfTrigInt.
 */
void
ncm_mpsf_sin_int_mpfr (mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  if (mpq_sgn (q) < 0)
  {
    /* Si is odd; negation swaps the directed rounding modes */
    const mp_rnd_t rnd_m = (rnd == MPFR_RNDU) ? MPFR_RNDD : ((rnd == MPFR_RNDD) ? MPFR_RNDU : rnd);
    mpq_t mq;

    mpq_init (mq);
    mpq_neg (mq, q);
    ncm_mpsf_sin_int_mpfr (mq, res, rnd_m);
    mpfr_neg (res, res, rnd);
    mpq_clear (mq);

    return;
  }

  {
    const gdouble x     = mpq_get_d (q);
    const gdouble dnmax = 1.0 / 4.0 * (-5.0 + sqrt (1.0 + 4.0 * x * x));
    const gdouble nmax  = (dnmax > 0) ? ceil (dnmax) : 0;
    const gulong prec   = mpfr_get_prec (res);


    if (nmax == 0)
    {
      _taylor_mpfr (q, res, rnd);
    }
    else
    {
      const gdouble maxsize = ((2.0 * nmax + 0.0) * log (x) - lgamma (1.0 + 2.0 * nmax)) / M_LN2;

      if (maxsize > prec)
        _assym_mpfr (q, res, rnd);
      else
        _taylor_mpfr (q, res, rnd);
    }
  }
}

/**
 * ncm_mpsf_sin_int_free_cache:
 *
 * Frees the pools of splitting buffers of the sine integral.
 */
void
ncm_mpsf_sin_int_free_cache (void)
{
  G_LOCK (__create_lock);

  if (__mp_taylor != NULL)
  {
    ncm_memory_pool_free (__mp_taylor, TRUE);
    ncm_memory_pool_free (__mp_assym, TRUE);
    __mp_taylor = NULL;
    __mp_assym  = NULL;
  }

  G_UNLOCK (__create_lock);
}

/**
 * ncm_sf_sin_int:
 * @x: the argument $x$
 *
 * Same as ncm_mpsf_sin_int_mpfr() at 53 bits, rounded to nearest, with @x converted to a
 * rational that agrees with it to $10^{-15}$, see ncm_rational_coarse_double().
 *
 * Returns: $\mathrm{Si}(x)$.
 */
gdouble
ncm_sf_sin_int (gdouble x)
{
  MPFR_DECL_INIT (res, 53);

  mpq_t xq;
  gdouble res_d;

  mpq_init (xq);

  ncm_rational_coarse_double (x, xq);
  ncm_mpsf_sin_int_mpfr (xq, res, GMP_RNDN);
  mpq_clear (xq);

  res_d = mpfr_get_d (res, GMP_RNDN);

  return res_d;
}

