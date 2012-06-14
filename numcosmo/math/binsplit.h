/***************************************************************************
 *            binsplit.h
 *
 *  Tue Jan 19 18:28:16 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
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

#ifndef _NC_BINSPLIT_H
#define _NC_BINSPLIT_H
#include <string.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <glib.h>

G_BEGIN_DECLS

typedef struct _NcmBinSplit NcmBinSplit;

/**
 * NcmBinSplit:
 * 
 * FIXME
 */
struct _NcmBinSplit
{
  /*< private >*/
  gpointer userdata;
  gulong n1;
  gulong n2;
  mpz_t P;
  mpz_t Q;
  mpz_t B;
  mpz_t T;
  mpz_t temp1;
  mpz_t temp2;
  mpz_t temp3;
  mpz_t temp4;
  NcmBinSplit *bs[2];
};

extern mpz_t NCM_BINSPLIT_ONE;

typedef void (*NcmBinSplitEval) (NcmBinSplit *bs, gulong n1, gulong n2);

NcmBinSplit *ncm_binsplit_alloc (gpointer userdata);
glong ncm_binsplit_test_next (NcmBinSplit *bs, NcmBinSplitEval bs_eval, gulong nt);
void ncm_binsplit_join (NcmBinSplit *bs, NcmBinSplit *bs_l, NcmBinSplit *bs_r);
void ncm_binsplit_eval_join (NcmBinSplit *bs, NcmBinSplitEval bs_eval, gulong nt);
gulong ncm_binsplit_eval_prec (NcmBinSplit *bs, NcmBinSplitEval bs_eval, gulong step, glong prec);
void ncm_binsplit_get (NcmBinSplit *bs, mpfr_t res);
void ncm_binsplit_get_q (NcmBinSplit *bs, mpq_t q);
gdouble ncm_binsplit_get_d (NcmBinSplit *bs, mp_rnd_t rnd);

#define NCM_BINSPLIT_DECL(name,v,u,n,data) static inline void name (mpz_t v, mpz_t u, gulong n, gpointer data)
#define NCM_BINSPLIT_DENC_NULL(a,b,c,d) 

G_END_DECLS

#endif /* _NC_BINSPLIT_H */

