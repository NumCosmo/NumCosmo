/***************************************************************************
 *            binsplit.h
 *
 *  Tue Jan 19 18:28:16 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com>
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
 * SECTION:binsplit
 * @title: NcmBinSplit
 * @short_description: Binnary splitting algorithms used to evaluate sums fast and with arbitrary precision.
 * @stability: Stable
 * @include: numcosmo/math/binsplit.h
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/binsplit.h"
#include "math/ncm_cfg.h"

static gboolean one_init = FALSE;
mpz_t NCM_BINSPLIT_ONE;

/**
 * ncm_binsplit_alloc: (skip)
 * @userdata: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmBinSplit *
ncm_binsplit_alloc (gpointer userdata)
{
  NcmBinSplit *bs = g_slice_new (NcmBinSplit);
  bs->userdata = userdata;

  mpz_init (bs->P);
  mpz_init (bs->Q);
  mpz_init (bs->B);
  mpz_init (bs->T);
  mpz_init (bs->temp1);
  mpz_init (bs->temp2);
  mpz_init (bs->temp3);
  mpz_init (bs->temp4);
  bs->bs[0] = bs->bs[1] = NULL;

  mpz_set_ui (bs->B, 1);

  {
    G_LOCK_DEFINE_STATIC (create_lock);
    G_LOCK (create_lock);
    if (!one_init)
    {
      mpz_init (NCM_BINSPLIT_ONE);
      mpz_set_ui (NCM_BINSPLIT_ONE, 1L);
      one_init = TRUE;
    }
    G_UNLOCK (create_lock);
  }

  return bs;
}

/**
 * ncm_binsplit_test_next: (skip)
 * @bs: a #NcmBinSplit
 * @bs_eval: a #NcmBinSplitEval
 * @nt: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
glong
ncm_binsplit_test_next (NcmBinSplit *bs, NcmBinSplitEval bs_eval, gulong nt)
{
  size_t s1, s2;
  if (bs->bs[0] == NULL)
    bs->bs[0] = ncm_binsplit_alloc (bs->userdata);

  bs_eval (bs->bs[0], bs->n2, bs->n2 + nt);

  mpz_mul (bs->temp1, bs->bs[0]->B, bs->bs[0]->Q);
  mpz_mul (bs->temp1, bs->temp1, bs->T);

  mpz_mul (bs->temp2, bs->B, bs->P);
  mpz_mul (bs->temp2, bs->temp2, bs->bs[0]->T);

  s1 = mpz_sizeinbase(bs->temp1, 2);
  s2 = mpz_sizeinbase(bs->temp2, 2);

//  mpz_mul (bs->temp3, bs->Q, bs->B);
//  mpfr_printf ("# %Zd/%Zd\n", bs->T, bs->temp3);
//  printf ("# testing [%lu %lu) %zu %zu %zd | %.15e\n", bs->n2, bs->n2 + nt, s1, s2, s1-s2, ncm_binsplit_get_d (bs, GMP_RNDN));

  return (s1 - s2);
}

/**
 * ncm_binsplit_join:
 * @bs: a #NcmBinSplit
 * @bs_l: a #NcmBinSplit
 * @bs_r: a #NcmBinSplit
 *
 * FIXME
 *
*/
void
ncm_binsplit_join (NcmBinSplit *bs, NcmBinSplit *bs_l, NcmBinSplit *bs_r)
{
  mpz_ptr temp5 = bs->P; /* To be used as temporary variable, but only after
                          * finishing all operations involving P from both sides
                          */
  g_assert ((bs != bs_l) || (bs != bs_r));
  g_assert (bs_l->n2 == bs_r->n1);

  mpz_mul (bs->temp1, bs_l->P, bs_r->P);
  mpz_mul (bs->temp2, bs_l->Q, bs_r->Q);
  mpz_mul (bs->temp3, bs_l->B, bs_r->B);

  mpz_mul (bs->temp4, bs_l->B, bs_l->P);
  mpz_mul (bs->temp4, bs->temp4, bs_r->T);

  mpz_mul (temp5, bs_r->B, bs_r->Q);
  mpz_mul (temp5, temp5, bs_l->T);

  mpz_add (bs->temp4, bs->temp4, temp5);

  mpz_swap (bs->P, bs->temp1);
  mpz_swap (bs->Q, bs->temp2);
  mpz_swap (bs->B, bs->temp3);
  mpz_swap (bs->T, bs->temp4);

  bs->n1 = bs_l->n1;
  bs->n2 = bs_r->n2;

  return;
}

/**
 * ncm_binsplit_eval_join: (skip)
 * @bs: a #NcmBinSplit
 * @bs_eval: a #NcmBinSplitEval
 * @nt: FIXME
 *
 * FIXME
 *
*/
void
ncm_binsplit_eval_join (NcmBinSplit *bs, NcmBinSplitEval bs_eval, gulong nt)
{
  gulong n3 = bs->n2 + nt;
  gulong nl = (bs->n2 + n3) / 2;

  if (bs->bs[0] == NULL)
    bs->bs[0] = ncm_binsplit_alloc (bs->userdata);
  if (bs->bs[1] == NULL)
    bs->bs[1] = ncm_binsplit_alloc (bs->userdata);

  if (nt == 1)
  {
    bs_eval (bs->bs[1], nl, n3);
    ncm_binsplit_join (bs, bs, bs->bs[1]);
  }
  else
  {
    bs_eval (bs->bs[0], bs->n2, nl);
    bs_eval (bs->bs[1], nl, n3);

    ncm_binsplit_join (bs->bs[0], bs->bs[0], bs->bs[1]);
    ncm_binsplit_join (bs, bs, bs->bs[0]);
  }
  return;
}

/**
 * ncm_binsplit_eval_prec: (skip)
 * @bs: a #NcmBinSplit
 * @bs_eval: a #NcmBinSplitEval
 * @step: FIXME
 * @prec: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gulong
ncm_binsplit_eval_prec (NcmBinSplit *bs, NcmBinSplitEval bs_eval, gulong step, glong prec)
{
  gulong tstep = 0;
  bs_eval (bs, 0, step);
  tstep += step;

  while (ncm_binsplit_test_next (bs, bs_eval, 1) < prec)
  {
    ncm_binsplit_eval_join (bs, bs_eval, step);
    tstep += step;
  }

  return tstep;
}

/**
 * ncm_binsplit_get: (skip)
 * @bs: a #NcmBinSplit
 * @res: FIXME
 *
 * FIXME
 *
*/
void
ncm_binsplit_get (NcmBinSplit *bs, mpfr_t res)
{
  mpfr_set_z (res, bs->T, GMP_RNDD);
  mpfr_div_z (res, res, bs->B, GMP_RNDD);
  mpfr_div_z (res, res, bs->Q, GMP_RNDD);
  return;
}

/**
 * ncm_binsplit_get_q: (skip)
 * @bs: a #NcmBinSplit
 * @q: FIXME
 *
 * FIXME
 *
*/
void
ncm_binsplit_get_q (NcmBinSplit *bs, mpq_t q)
{
  mpz_set (mpq_numref(q), bs->T);
  mpz_set (mpq_denref(q), bs->B);
  mpz_mul (mpq_denref(q), mpq_denref(q), bs->Q);
  mpq_canonicalize (q);

  return;
}

/**
 * ncm_binsplit_get_d: (skip)
 * @bs: a #NcmBinSplit
 * @rnd: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_binsplit_get_d (NcmBinSplit *bs, mp_rnd_t rnd)
{
  MPFR_DECL_INIT (res, 64);
  mpfr_set_z (res, bs->T, rnd);
  mpfr_div_z (res, res, bs->B, rnd);
  mpfr_div_z (res, res, bs->Q, rnd);
  return mpfr_get_d (res, rnd);
}
