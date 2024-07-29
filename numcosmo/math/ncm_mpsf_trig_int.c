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
 * SECTION:ncm_mpsf_trig_int
 * @title: NcmMpsfTrigInt
 * @short_description: Multiple precision sin integral implementation.
 *
 * Implementation of multiple precision sine integral using the GNU MPFR library.
 * This implementation employs binary splitting to compute the integral utilizing
 * the Taylor series and asymptotic expansion methods.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mpsf_trig_int.h"
#include "math/ncm_binsplit.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

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

#include "ncm_binsplit_eval.c"
#undef _mx2

#define NC_BINSPLIT_EVAL_NAME binsplit_sin_integral_assym
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

#include "ncm_binsplit_eval.c"
#undef _m2x2

static void
_taylor_mpfr (mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  static NcmBinSplit *bs = NULL;
  static mpq_ptr mq2     = NULL;

  if (mq2 == NULL)
  {
    mq2 = g_slice_new (__mpq_struct);
    mpq_init (mq2);
  }

  mpq_mul (mq2, q, q);
  mpq_neg (mq2, mq2);
  mpq_div_2exp (mq2, mq2, 1);

  if (bs == NULL)
    bs = ncm_binsplit_alloc ((gpointer) mq2);

  ncm_binsplit_eval_prec (bs, binsplit_sin_integral_taylor, 10, mpfr_get_prec (res));
  ncm_binsplit_get (bs, res);
  mpfr_mul_q (res, res, q, rnd);
}

static void
_assym_mpfr (mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  static NcmBinSplit *bs = NULL;
  glong prec             = mpfr_get_prec (res);
  glong mprec;

  mpfr_set_q (res, q, rnd);
  mprec = res->_mpfr_exp;

  if ((prec - mprec) > 0)
  {
    gulong nf;
    static gpointer *data = NULL;
    mpq_t mq2_2;

    mpq_init (mq2_2);
    mpq_mul (mq2_2, q, q);
    mpq_neg (mq2_2, mq2_2);
    mpq_div_2exp (mq2_2, mq2_2, 1);
    MPFR_DECL_INIT (sin_x, prec);
    MPFR_DECL_INIT (cos_x, prec);

    if (data == NULL)
      data = g_slice_alloc (2 * sizeof (gpointer));

    if (bs == NULL)
      bs = ncm_binsplit_alloc (data);

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

    mpq_clear (mq2_2);
  }
  else
  {
    mpfr_const_pi (res, rnd);
    mpfr_div_2ui (res, res, 1, rnd);
  }
}

/**
 * ncm_mpsf_sin_int_mpfr: (skip)
 * @q: argument as a rational number $x = q$
 * @res: mpfr variable containing the result $\mathrm{Si}(x)$
 * @rnd: mpfr rounding mode
 *
 * Computes the sine integral
 * $$\mathrm{Si}(x) = \int_0^x \frac{\sin x'}{x'}\mathrm{d}x.$$
 *
 */
void
ncm_mpsf_sin_int_mpfr (mpq_t q, mpfr_ptr res, mp_rnd_t rnd)
{
  const gdouble x     = mpq_get_d (q);
  const gdouble dnmax = 1.0 / 4.0 * (-5.0 + sqrt (1.0 + 4.0 * x * x));
  const gdouble nmax  = (dnmax > 0) ? ceil (dnmax) : 0;
  const gulong prec   = mpfr_get_prec (res);

  /*printf ("# dnmax %g nmax %g\n", dnmax, nmax); */

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

/**
 * ncm_sf_sin_int:
 * @x: value of the argument $x$
 *
 * Computes the sine integral
 * $$\mathrm{Si}(x) = \int_0^x \frac{\sin x'}{x'}\mathrm{d}x.$$
 *
 * Returns: the value of $\mathrm{Si}(x)$.
 */
gdouble
ncm_sf_sin_int (gdouble x)
{
  MPFR_DECL_INIT (res, 53);

  mpq_t xq;
  gdouble res_d;

  mpq_init (xq);

  ncm_rational_coarce_double (x, xq);
  ncm_mpsf_sin_int_mpfr (xq, res, GMP_RNDN);
  mpq_clear (xq);

  res_d = mpfr_get_d (res, GMP_RNDN);

  return res_d;
}

