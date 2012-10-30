/***************************************************************************
 *            ncm_mpsf_sbessel_int.c
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

/**
 * SECTION:ncm_mpsf_sbessel_int
 * @title: Multiple Precision Spherical Bessel j_l integral
 * @short_description: Spherical bessel integrals implementation with support for multiple precision calculation
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mpsf_sbessel_int.h"
#include "math/mpq_tree.h"
#include "math/util.h"
#include "math/ncm_cfg.h"
#include "math/binsplit.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_mpsf_trig_int.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_new: (skip)
 * @prec: FIXME
 * @jlrec: a #NcmMpsfSBesselRecur
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmMpsfSBesselIntegRecur *
ncm_mpsf_sbessel_jl_xj_integral_recur_new (gulong prec, NcmMpsfSBesselRecur *jlrec)
{
  NcmMpsfSBesselIntegRecur *xnjlrec = g_slice_new (NcmMpsfSBesselIntegRecur);

  if (jlrec == NULL)
    xnjlrec->jlrec = ncm_mpsf_sbessel_recur_new (prec);
  else
    xnjlrec->jlrec = jlrec;

  g_assert (prec == xnjlrec->jlrec->prec);

  xnjlrec->prec = prec;

  mpfr_inits2 (prec,
    xnjlrec->int_jl_xn[0], xnjlrec->int_jl_xn[1], xnjlrec->int_jl_xn[2], xnjlrec->int_jl_xn[3],
    xnjlrec->int_jlp1_xn[0], xnjlrec->int_jlp1_xn[1], xnjlrec->int_jlp1_xn[2], xnjlrec->int_jlp1_xn[3],
    xnjlrec->x_pow_n[0], xnjlrec->x_pow_n[1], xnjlrec->x_pow_n[2], xnjlrec->x_pow_n[3],
    xnjlrec->temp1, xnjlrec->temp2,
    NULL);

  xnjlrec->x = xnjlrec->jlrec->x;
  xnjlrec->q = xnjlrec->jlrec->q;

  return xnjlrec;
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_set_q: (skip)
 * @xnjlrec: a NcmMpsfSBesselIntegRecur
 * @l: FIXME
 * @q: FIXME
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_jl_xj_integral_recur_set_q (NcmMpsfSBesselIntegRecur *xnjlrec, glong l, mpq_ptr q, mp_rnd_t rnd)
{
  gint i;

  ncm_mpsf_sbessel_recur_set_q (xnjlrec->jlrec, l, q, rnd);

  mpfr_set_ui (xnjlrec->x_pow_n[0], 1, rnd);
  mpfr_mul (xnjlrec->x_pow_n[1], xnjlrec->x_pow_n[0], xnjlrec->x, rnd);
  mpfr_mul (xnjlrec->x_pow_n[2], xnjlrec->x_pow_n[1], xnjlrec->x, rnd);
  mpfr_mul (xnjlrec->x_pow_n[3], xnjlrec->x_pow_n[2], xnjlrec->x, rnd);

  for (i = 0; i < 4; i++)
  {
    ncm_mpsf_sbessel_jl_xj_integral_q (l + 0, i, xnjlrec->q, xnjlrec->int_jl_xn[i], rnd);
    ncm_mpsf_sbessel_jl_xj_integral_q (l + 1, i, xnjlrec->q, xnjlrec->int_jlp1_xn[i], rnd);
  }
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_set_d: (skip)
 * @xnjlrec: a NcmMpsfSBesselIntegRecur
 * @l: FIXME
 * @x: FIXME
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_jl_xj_integral_recur_set_d (NcmMpsfSBesselIntegRecur *xnjlrec, glong l, gdouble x, mp_rnd_t rnd)
{
  gint i;
  ncm_mpsf_sbessel_recur_set_d (xnjlrec->jlrec, l, x, rnd);

  mpfr_set_ui (xnjlrec->x_pow_n[0], 1, rnd);
  mpfr_mul (xnjlrec->x_pow_n[1], xnjlrec->x_pow_n[0], xnjlrec->x, rnd);
  mpfr_mul (xnjlrec->x_pow_n[2], xnjlrec->x_pow_n[1], xnjlrec->x, rnd);
  mpfr_mul (xnjlrec->x_pow_n[3], xnjlrec->x_pow_n[2], xnjlrec->x, rnd);

  for (i = 0; i < 4; i++)
  {
    ncm_mpsf_sbessel_jl_xj_integral_q (l + 0, i, xnjlrec->q, xnjlrec->int_jl_xn[i], rnd);
    ncm_mpsf_sbessel_jl_xj_integral_q (l + 1, i, xnjlrec->q, xnjlrec->int_jlp1_xn[i], rnd);
  }
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_free: (skip)
 * @xnjlrec: a NcmMpsfSBesselIntegRecur
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_jl_xj_integral_recur_free (NcmMpsfSBesselIntegRecur *xnjlrec)
{
  ncm_mpsf_sbessel_recur_free (xnjlrec->jlrec);
  mpfr_clears (
    xnjlrec->int_jl_xn[0], xnjlrec->int_jl_xn[1], xnjlrec->int_jl_xn[2], xnjlrec->int_jl_xn[3],
    xnjlrec->int_jlp1_xn[0], xnjlrec->int_jlp1_xn[1], xnjlrec->int_jlp1_xn[2], xnjlrec->int_jlp1_xn[3],
    xnjlrec->temp1, xnjlrec->temp2,
    NULL);
  g_slice_free (NcmMpsfSBesselIntegRecur, xnjlrec);
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_next: (skip)
 * @xnjlrec: a NcmMpsfSBesselIntegRecur
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_jl_xj_integral_recur_next (NcmMpsfSBesselIntegRecur *xnjlrec, mp_rnd_t rnd)
{
  gint i;
  const glong l = xnjlrec->jlrec->l;

  ncm_mpsf_sbessel_recur_next (xnjlrec->jlrec, rnd);

  for (i = 0; i < 4; i++)
  {
    mpfr_mul (xnjlrec->temp1, xnjlrec->x_pow_n[i], xnjlrec->jlrec->jl[0], rnd);

    mpfr_mul_ui (xnjlrec->temp1, xnjlrec->temp1, 2 * l + 3, rnd);
    mpfr_mul_ui (xnjlrec->temp2, xnjlrec->int_jl_xn[i], i + l + 1, rnd);

    mpfr_swap (xnjlrec->int_jl_xn[i], xnjlrec->int_jlp1_xn[i]);
    mpfr_sub (xnjlrec->int_jlp1_xn[i], xnjlrec->temp1, xnjlrec->temp2, rnd);
    mpfr_div_si (xnjlrec->int_jlp1_xn[i], xnjlrec->int_jlp1_xn[i], i - l - 2, rnd);
  }
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_previous: (skip)
 * @xnjlrec: a NcmMpsfSBesselIntegRecur
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_jl_xj_integral_recur_previous (NcmMpsfSBesselIntegRecur *xnjlrec, mp_rnd_t rnd)
{
  gint i;
  const glong l = xnjlrec->jlrec->l;

  ncm_mpsf_sbessel_recur_previous (xnjlrec->jlrec, rnd);

  for (i = 0; i < 4; i++)
  {
    mpfr_mul (xnjlrec->temp1, xnjlrec->x_pow_n[i], xnjlrec->jlrec->jl[1], rnd);

    mpfr_mul_ui (xnjlrec->temp1, xnjlrec->temp1, 2 * l + 1, rnd);
    mpfr_mul_si (xnjlrec->temp2, xnjlrec->int_jlp1_xn[i], i - l - 1, rnd);

    mpfr_swap (xnjlrec->int_jl_xn[i], xnjlrec->int_jlp1_xn[i]);
    mpfr_sub (xnjlrec->int_jl_xn[i], xnjlrec->temp1, xnjlrec->temp2, rnd);
    mpfr_div_ui (xnjlrec->int_jl_xn[i], xnjlrec->int_jl_xn[i], i + l, rnd);
  }
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_goto: (skip)
 * @xnjlrec: a NcmMpsfSBesselIntegRecur
 * @l: FIXME
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_jl_xj_integral_recur_goto (NcmMpsfSBesselIntegRecur *xnjlrec, glong l, mp_rnd_t rnd)
{
  glong sign = GSL_SIGN (l - xnjlrec->jlrec->l);
  glong sub = labs(l - xnjlrec->jlrec->l);
  glong i;
  if (sub == 0)
    return;
  if (sign == 1)
    for (i = 0; i < sub; i++)
      ncm_mpsf_sbessel_jl_xj_integral_recur_next (xnjlrec, rnd);
  else
    for (i = 0; i < sub; i++)
      ncm_mpsf_sbessel_jl_xj_integral_recur_previous (xnjlrec, rnd);
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_calc_djl: (skip)
 * @xnjlrec: a NcmMpsfSBesselIntegRecur
 * @n: FIXME
 * @rule: FIXME
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_jl_xj_integral_recur_calc_djl (NcmMpsfSBesselIntegRecur *xnjlrec, glong n, mpfr_ptr rule, mp_rnd_t rnd)
{
  g_assert (n >= 0 && n < 4);

  mpfr_mul (rule, xnjlrec->x_pow_n[n], xnjlrec->jlrec->jl[0], rnd);
  if (n > 0)
  {
    mpfr_mul_ui (xnjlrec->temp1, xnjlrec->int_jl_xn[n-1], n, rnd);
    mpfr_sub (rule, rule, xnjlrec->temp1, rnd);
  }
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_calc_d2jl: (skip)
 * @xnjlrec: a NcmMpsfSBesselIntegRecur
 * @n: FIXME
 * @rule: FIXME
 * @rnd: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_jl_xj_integral_recur_calc_d2jl (NcmMpsfSBesselIntegRecur *xnjlrec, glong n, mpfr_ptr rule, mp_rnd_t rnd)
{
  const glong l = xnjlrec->jlrec->l;
  g_assert (n >= 0 && n < 4);

  if (l == n)
    mpfr_set_ui (rule, 0, rnd);
  else if (n == 0)
  {
    if (mpfr_sgn (xnjlrec->x) == 0)
      mpfr_set_ui (rule, 0, rnd);
    else
    {
      mpfr_div (rule, xnjlrec->jlrec->jl[0], xnjlrec->x, rnd);
      mpfr_mul_si (rule, rule, l, rnd);
    }
  }
  else if (n > 0)
  {
    mpfr_mul (rule, xnjlrec->x_pow_n[n-1], xnjlrec->jlrec->jl[0], rnd);
    mpfr_mul_si (rule, rule, l - n, rnd);
  }

  mpfr_mul (xnjlrec->temp1, xnjlrec->x_pow_n[n], xnjlrec->jlrec->jl[1], rnd);
  mpfr_sub (rule, rule, xnjlrec->temp1, rnd);

  if (n > 1)
  {
    mpfr_mul_ui (xnjlrec->temp1, xnjlrec->int_jl_xn[n-2], n * (n - 1), rnd);
    mpfr_add (rule, rule, xnjlrec->temp1, rnd);
  }
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_write: (skip)
 * @xnjlrec: a NcmMpsfSBesselIntegRecur
 * @f: FIXME
 *
 * FIXME
 *
*/
void
ncm_mpsf_sbessel_jl_xj_integral_recur_write (NcmMpsfSBesselIntegRecur *xnjlrec, FILE *f)
{
  gint i;
  ncm_mpsf_sbessel_recur_write (xnjlrec->jlrec, f);
  for (i = 0; i < 4; i++)
  {
    ncm_mpfr_out_raw (f, xnjlrec->int_jl_xn[i]);
    ncm_mpfr_out_raw (f, xnjlrec->int_jlp1_xn[i]);
    ncm_mpfr_out_raw (f, xnjlrec->x_pow_n[i]);
  }
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_recur_read: (skip)
 * @f: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmMpsfSBesselIntegRecur *
ncm_mpsf_sbessel_jl_xj_integral_recur_read (FILE *f)
{
  gint i;
  NcmMpsfSBesselRecur *jlrec = ncm_mpsf_sbessel_recur_read (f);
  NcmMpsfSBesselIntegRecur *xnjlrec =
    ncm_mpsf_sbessel_jl_xj_integral_recur_new (jlrec->prec, jlrec);

  xnjlrec->prec = xnjlrec->prec;
  xnjlrec->x = xnjlrec->x;
  xnjlrec->q = xnjlrec->q;

  for (i = 0; i < 4; i++)
  {
    ncm_mpfr_inp_raw (xnjlrec->int_jl_xn[i], f);
    ncm_mpfr_inp_raw (xnjlrec->int_jlp1_xn[i], f);
    ncm_mpfr_inp_raw (xnjlrec->x_pow_n[i], f);
  }

  return xnjlrec;
}

/*********************************************************************************************************
 *
 *
 *********************************************************************************************************/

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_new: (skip)
 * @prec: FIXME
 * @x: a #NcmGrid
 * @k: a #NcmGrid
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmMpsfSBesselIntSpline *
ncm_mpsf_sbessel_jl_xj_integrate_spline_new (gulong prec, NcmGrid *x, NcmGrid *k)
{
  NcmMpsfSBesselIntSpline *int_jlspline = g_slice_new (NcmMpsfSBesselIntSpline);
  gulong map_size;
  int_jlspline->x = x;
  int_jlspline->k = k;
  map_size = int_jlspline->x->nnodes * int_jlspline->k->nnodes;

  mpfr_inits2 (prec,
    int_jlspline->rules[0], int_jlspline->rules[1], int_jlspline->rules[2], int_jlspline->rules[3],
    int_jlspline->crules[0], int_jlspline->crules[1], int_jlspline->crules[2], int_jlspline->crules[3],
    NULL);

  int_jlspline->prec = prec;

  int_jlspline->l = 0;

  int_jlspline->int_hash = ncm_mpq_hash_new ();

  int_jlspline->xnjlrec = g_slice_alloc0 (map_size * sizeof (NcmMpsfSBesselIntegRecur *));
  int_jlspline->map_ij2r = g_slice_alloc0 (map_size * sizeof (guint));
  int_jlspline->row = int_jlspline->x->nnodes;

  int_jlspline->prepared = FALSE;
  return int_jlspline;
}

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_new_from_sections: (skip)
 * @prec: FIXME
 * @x_secs: a #NcmGridSection
 * @k_secs: a #NcmGridSection
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmMpsfSBesselIntSpline *
ncm_mpsf_sbessel_jl_xj_integrate_spline_new_from_sections (gulong prec, NcmGridSection *x_secs, NcmGridSection *k_secs)
{
  NcmGrid *x = ncm_grid_new_from_sections (x_secs);
  NcmGrid *k = ncm_grid_new_from_sections (k_secs);
  return ncm_mpsf_sbessel_jl_xj_integrate_spline_new (prec, x, k);
}

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_cached_new: (skip)
 * @prec: FIXME
 * @l: FIXME
 * @x_secs: a #NcmGridSection
 * @k_secs: a #NcmGridSection
 * @rnd: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmMpsfSBesselIntSpline *
ncm_mpsf_sbessel_jl_xj_integrate_spline_cached_new (gulong prec, glong l, NcmGridSection *x_secs, NcmGridSection *k_secs, mp_rnd_t rnd)
{
  NcmMpsfSBesselIntSpline *int_jlspline;
  gchar *name_x = ncm_grid_get_name (x_secs);
  gchar *name_k = ncm_grid_get_name (k_secs);
  gchar *rname = g_strconcat (name_x, name_k, NULL);
  gchar *rname_hash = g_compute_checksum_for_string (G_CHECKSUM_MD5, rname, strlen (rname));
  gchar *filename = g_strdup_printf ("xnjl_rule_%ld_%lu_%s.dat", l, prec, rname_hash);
  printf ("# looking for cache (%s) name (%s)\n", filename, rname);

  if (ncm_cfg_exists (filename))
    int_jlspline = ncm_mpsf_sbessel_jl_xj_integrate_spline_load (filename);
  else
  {
    int_jlspline = ncm_mpsf_sbessel_jl_xj_integrate_spline_new_from_sections (prec, x_secs,  k_secs);
    ncm_mpsf_sbessel_jl_xj_integrate_spline_prepare (int_jlspline, l, rnd);
    ncm_mpsf_sbessel_jl_xj_integrate_spline_save (int_jlspline, filename);
  }

  g_free (name_x);
  g_free (name_k);
  g_free (rname);
  g_free (rname_hash);
  g_free (filename);
  return int_jlspline;
}

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_prepare_new: (skip)
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @l: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integrate_spline_prepare_new (NcmMpsfSBesselIntSpline *int_jlspline, glong l, mp_rnd_t rnd)
{
  gdouble xi, xf, dx;
  guint np, i;
  xi = mpq_get_d (int_jlspline->x->nodes[0]) * mpq_get_d (int_jlspline->k->nodes[0]);
  xf = mpq_get_d (int_jlspline->x->nodes[int_jlspline->x->nnodes-1]) * mpq_get_d (int_jlspline->k->nodes[int_jlspline->k->nnodes-1]);
  dx = (xf-xi);
  np = ceil (dx / (2.0 * M_PI) * 10.0);
  for (i = 0; i < np; i++)
  {
    gdouble x = xi + dx / (np - 1.0) * i;
    NcmMpsfSBesselIntegRecur *xnjlrec = ncm_mpsf_sbessel_jl_xj_integral_recur_new (int_jlspline->prec, NULL);
    ncm_mpsf_sbessel_jl_xj_integral_recur_set_d (xnjlrec, l, x, rnd);
    printf ("%u %u %.15g\n", i, np, x);
  }
  g_assert_not_reached ();
}

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_prepare: (skip)
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @l: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integrate_spline_prepare (NcmMpsfSBesselIntSpline *int_jlspline, glong l, mp_rnd_t rnd)
{
  gulong i, j;
  GTimer *bench = g_timer_new ();
  gdouble elapsed[5] = {0, 0, 0, 0, 0};
  mpq_t q;
  GHashTable *rhash = g_hash_table_new (NULL, NULL);
  mpq_init (q);
  int_jlspline->l = l;
  int_jlspline->rules_count = 0;

  for (j = 0; j < int_jlspline->k->nnodes; j++)
  {
    for (i = 0; i < int_jlspline->x->nnodes; i++)
    {
      NcmMpsfSBesselIntegRecur *xnjlrec;
      mpq_mul (q, int_jlspline->x->nodes[i], int_jlspline->k->nodes[j]);
      //mpfr_printf ("# %lu %lu | %Qd\n", i, j, q); fflush (stdout);
      if ((xnjlrec = g_hash_table_lookup (int_jlspline->int_hash, q)) == NULL)
      {
        xnjlrec = ncm_mpsf_sbessel_jl_xj_integral_recur_new (int_jlspline->prec, NULL);
        ncm_mpsf_sbessel_jl_xj_integral_recur_set_q (xnjlrec, l, q, rnd);
        mpq_ptr qq = (mpq_ptr)g_slice_new (mpq_t);
        mpq_init (qq);
        mpq_set (qq, q);
        g_hash_table_insert (int_jlspline->int_hash, qq, xnjlrec);

        int_jlspline->xnjlrec[int_jlspline->rules_count] = xnjlrec;

        g_hash_table_insert (rhash, xnjlrec, GINT_TO_POINTER(int_jlspline->rules_count));

        NCM_MPSF_SBESSEL_INT_MAP (int_jlspline,i,j) = int_jlspline->rules_count;
        int_jlspline->rules_count++;
      }
      else
        NCM_MPSF_SBESSEL_INT_MAP(int_jlspline,i,j) =
          GPOINTER_TO_INT(g_hash_table_lookup (rhash, xnjlrec));
    }

    elapsed[0] = elapsed[1];
    elapsed[1] = elapsed[2];
    elapsed[2] = elapsed[3];
    elapsed[3] = elapsed[4];
    elapsed[4] = g_timer_elapsed (bench, NULL);

    printf ("# %lu %lu | total %fs | estimated left %fs %fs\n", i, j, elapsed[4], elapsed[4] * (int_jlspline->k->nnodes - j - 1) / (j + 1.0), (elapsed[4]-elapsed[3]) * (int_jlspline->k->nnodes - j - 1) + 0.0 * (elapsed[4] - 2.0 * elapsed[3] + elapsed[2]) * gsl_pow_2 (int_jlspline->k->nnodes - j - 1) / 2.0 ); fflush (stdout);
  }

  g_hash_table_unref (rhash);
  int_jlspline->prepared = TRUE;
  g_timer_destroy (bench);
  mpq_clear (q);
}

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_save:
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @filename: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integrate_spline_save (NcmMpsfSBesselIntSpline *int_jlspline, gchar *filename)
{
  gulong i, j;
  g_assert (int_jlspline->prepared == TRUE);
  FILE *f = ncm_cfg_fopen (filename, "w");
  gchar magic[5] = "NcSBR";

  if (fwrite (magic, 1, 5, f) != 5)
    g_error ("ncm_mpsf_sbessel_jl_xj_integrate_spline_save: io error");

  NCM_WRITE_UINT32(f, int_jlspline->prec);
  ncm_grid_write (int_jlspline->x, f);
  ncm_grid_write (int_jlspline->k, f);

  for (j = 0; j < int_jlspline->k->nnodes; j++)
    for (i = 0; i < int_jlspline->x->nnodes; i++)
      NCM_WRITE_UINT32(f, NCM_MPSF_SBESSEL_INT_MAP(int_jlspline,i,j));

  NCM_WRITE_UINT32(f, int_jlspline->rules_count);

  for (i = 0; i < int_jlspline->rules_count; i++)
    ncm_mpsf_sbessel_jl_xj_integral_recur_write (int_jlspline->xnjlrec[i], f);

  fclose (f);
}

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_load: (skip)
 * @filename: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmMpsfSBesselIntSpline *
ncm_mpsf_sbessel_jl_xj_integrate_spline_load (gchar *filename)
{
  gulong i, j;
  FILE *f = ncm_cfg_fopen (filename, "r");
  gchar magic[5];
  NcmGrid *x_grid, *k_grid;
  NcmMpsfSBesselIntSpline *int_jlspline;
  guint32 prec;

  if (fread (magic, 1, 5, f) != 5)
    g_error ("ncm_mpsf_sbessel_jl_xj_integrate_spline_load: io error");
  NCM_READ_UINT32(f, prec);
  x_grid = ncm_grid_read (f);
  k_grid = ncm_grid_read (f);
  int_jlspline = ncm_mpsf_sbessel_jl_xj_integrate_spline_new (prec, x_grid, k_grid);

  for (j = 0; j < int_jlspline->k->nnodes; j++)
    for (i = 0; i < int_jlspline->x->nnodes; i++)
      NCM_READ_UINT32(f, NCM_MPSF_SBESSEL_INT_MAP(int_jlspline,i,j));

  NCM_READ_UINT32(f, int_jlspline->rules_count);

  for (i = 0; i < int_jlspline->rules_count; i++)
    int_jlspline->xnjlrec[i] = ncm_mpsf_sbessel_jl_xj_integral_recur_read (f);

  fclose (f);
  return int_jlspline;
}

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_next: (skip)
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integrate_spline_next (NcmMpsfSBesselIntSpline *int_jlspline, mp_rnd_t rnd)
{
  guint i;
  for (i = 0; i < int_jlspline->rules_count; i++)
    ncm_mpsf_sbessel_jl_xj_integral_recur_next (int_jlspline->xnjlrec[i], rnd);
}

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_previous: (skip)
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integrate_spline_previous (NcmMpsfSBesselIntSpline *int_jlspline, mp_rnd_t rnd)
{
  guint i;
  for (i = 0; i < int_jlspline->rules_count; i++)
    ncm_mpsf_sbessel_jl_xj_integral_recur_previous (int_jlspline->xnjlrec[i], rnd);
}

/**
 * ncm_mpsf_sbessel_jl_xj_integrate_spline_goto: (skip)
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @l: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integrate_spline_goto (NcmMpsfSBesselIntSpline *int_jlspline, gulong l, mp_rnd_t rnd)
{
  guint i;
  for (i = 0; i < int_jlspline->rules_count; i++)
  {
//    if (i % 1000 == 0 ){ printf ("# Rule going %u -> %lu %u\n", int_jlspline->xnjlrec[i]->jlrec->l, l, i); fflush (stdout); }
    ncm_mpsf_sbessel_jl_xj_integral_recur_goto (int_jlspline->xnjlrec[i], l, rnd);
  }
}

/*************************************************************************************
 * Code to calculate the integral of  x^n j_l , taylor series version x <= l
 *************************************************************************************/

#define NC_BINSPLIT_EVAL_NAME binsplit_int_xn_jl
#define _l (((binsplit_int_xn_jl_data *)data)->l)
#define _j (((binsplit_int_xn_jl_data *)data)->j)
#define _x (((binsplit_int_xn_jl_data *)data)->q)

typedef struct _binsplit_int_xn_jl_data
{
  glong l;
  glong j;
  mpq_ptr q;
} binsplit_int_xn_jl_data;

NCM_BINSPLIT_DECL(xnjl_int_p,v,u,n,data)
{
  if (n == 0)
    mpz_set (v, u);
  else
  {
    mpz_mul (v, u, mpq_numref (_x));
    mpz_mul_2exp (v, v, 1);
  }
}
#define _BINSPLIT_FUNC_P xnjl_int_p

NCM_BINSPLIT_DECL(xnjl_int_q,v,u,n,data)
{
  if (n == 0)
    mpz_set (v, u);
  else
  {
    mpz_mul_ui (v, u, (2L * _l + 2L * n + 1L) * n);
    mpz_mul (v, v, mpq_denref(_x));
  }
}
#define _BINSPLIT_FUNC_Q xnjl_int_q

NCM_BINSPLIT_DECL(xnjl_int_b,v,u,n,data)
{
  mpz_mul_ui (v, u, 1L + _j + _l + 2L * n);
}
#define _BINSPLIT_FUNC_B xnjl_int_b
#define _HAS_FUNC_B

#define _BINSPLIT_FUNC_A NCM_BINSPLIT_DENC_NULL

#include "binsplit_eval.c"
#undef _l
#undef _j
#undef _xn
#undef _xd
#undef _x2exp

inline static gdouble
_taylor_get_nmax (gint l, gint j, gdouble mx2_4)
{
  gdouble rs[3] = {0.0, 0.0, 0.0};
  gint nr = gsl_poly_solve_cubic (
    (8.0 + j + 3.0 * l) / 2.0,
    (21.0 + 5.0 * j + (15.0 + 2.0 * (j + l)) * l - 4.0 * fabs(mx2_4)) / 4.0,
    (9.0 + 3.0 * j + (9.0 + 2.0 * (j + l)) * l - 2.0 * (1.0 + j + l) * fabs(mx2_4)) / 4.0,
    &rs[0], &rs[1], &rs[2]);
  return ceil(rs[nr - 1]);
}

static void
_xn_jl_taylor (gint l, gint j, mpq_t q, mpfr_t res, mp_rnd_t rnd)
{
  static NcmBinSplit *bs = NULL;
  const gdouble x = mpq_get_d (q);
  const gdouble mx2_4 = -x * x / 4.0;
  const gdouble nmax = _taylor_get_nmax (l, j, mx2_4);
  gulong prec = mpfr_get_prec (res);
  static binsplit_int_xn_jl_data *data = NULL;
  mpfr_t sqrt_pi, gamma_l_3_2, x_pow_1_l_j;
  mpq_t mq2_4;

  g_assert (x >= 0);
//  g_assert (x <= l);

  if (data == NULL)
    data = g_slice_new (binsplit_int_xn_jl_data);
  if (bs == NULL)
    bs = ncm_binsplit_alloc ((gpointer)data);

  mpfr_inits2 (prec, sqrt_pi, gamma_l_3_2, x_pow_1_l_j, NULL);
  mpfr_const_pi (sqrt_pi, rnd);
  mpfr_sqrt (sqrt_pi, sqrt_pi, rnd);
  mpq_init (mq2_4);

  mpq_mul (mq2_4, q, q);
  mpq_div_2exp (mq2_4, mq2_4, 2);
  mpq_neg (mq2_4, mq2_4);

  data->l = l;
  data->j = j;
  data->q = mq2_4;

  ncm_binsplit_eval_prec (bs, binsplit_int_xn_jl, (nmax < 10) ? 10 : nmax, prec);
  ncm_binsplit_get (bs, res);

  mpfr_mul (res, res, sqrt_pi, rnd);
  mpfr_set_ui (gamma_l_3_2, 2L * l + 3L, rnd);
  mpfr_div_2ui (gamma_l_3_2, gamma_l_3_2, 1, rnd);
  mpfr_gamma (gamma_l_3_2, gamma_l_3_2, rnd);
  mpfr_div (res, res, gamma_l_3_2, rnd);

  mpfr_set_q (x_pow_1_l_j, q, rnd);
  mpfr_pow_ui (x_pow_1_l_j, x_pow_1_l_j, 1 + l + j, rnd);
  mpfr_mul (res, res, x_pow_1_l_j, rnd);

  mpfr_div_2ui (res, res, l + 1, rnd);

  mpq_clear (mq2_4);
  mpfr_clears (sqrt_pi, gamma_l_3_2, x_pow_1_l_j, NULL);
}


/*************************************************************************************
 * Code to calculate the integral of  x^n j_l , finite series version x > l
 *************************************************************************************/

static void
_xn_jl_finite_x (gint l, gint j, mpq_t q, mpz_t sum_sin, mpz_t sum_cos, mpz_t cor)
{
  static mpz_t term, qn2;
  static gboolean init_mpz = TRUE;
  mpz_ptr sum = cor;
  gulong n, signs[2] = {0, 0};
  gulong first_term, second_term;
  if (init_mpz)
  {
    mpz_inits (term, sum, qn2, NULL);
    init_mpz = FALSE;
  }

  first_term = ((1-l) < 0) ? ((4 - ((l - 1) % 4)) % 4) : ((1 - l) % 4);
  second_term = (first_term + 1) % 4;

  signs[first_term  % 2] = (first_term  / 2) == 0 ? 1 : 0; /* There is an overall -1 sign          */
  signs[second_term % 2] = (second_term / 2) == 0 ? 1 : 0; /* Thats why the signs are changed here */

  mpz_set_ui (term, 1L);
  mpz_set_ui (sum, 1L);
  if (first_term % 2 == 0)
  {
    mpz_set_si (sum_sin, signs[first_term  % 2]++ == 0 ? 1 : -1);
    mpz_set_si (sum_cos, 0);
  }
  else
  {
    mpz_set_si (sum_sin, 0);
    mpz_set_si (sum_cos, signs[first_term  % 2]++ == 0 ? 1 : -1);
  }

  mpz_mul (qn2, mpq_numref (q), mpq_numref (q));

  for (n = 0; n < j - 1; n++)
  {
    glong wsum = abs (l - (n + 1) - 1) % 2;
    mpz_ptr lsum = (wsum == 0) ? sum_sin : sum_cos;
    mpz_mul_ui (term, term, (l - n) * (1L + l + n));
    mpz_divexact_ui (term, term, 2L * (1L + n));
    mpz_mul (term, term, mpq_denref (q));

    mpz_mul_ui (sum, sum, (j - 1 - n));
    mpz_mul (sum, sum, mpq_denref (q));

    mpz_add (sum, sum, term);

    mpz_mul (lsum, lsum, qn2);
    ((signs[wsum]++) % 2 == 0) ? mpz_add (lsum, lsum, sum) : mpz_sub (lsum, lsum, sum);
  }

  if (abs (l-j-2) % 2 == 0) /* Last index n = j - 1 */
    mpz_mul (sum_cos, sum_cos, mpq_numref(q));
  else
    mpz_mul (sum_sin, sum_sin, mpq_numref(q));
}

static void
_xn_jl_finite_inverse_x (gint l, gint j, mpq_t q, mpz_t sum_sin, mpz_t sum_cos, mpz_t sum_sinint, mpz_t cor, gulong prec, gulong *lk)
{
  static mpz_t term, termc, sumc, sum, term_const, four_xn2, qd_jm1;
  static gboolean init_mpz = TRUE;
  gboolean converged[2] = {FALSE, FALSE};
  glong signs[2] = {1, 1};
  gulong n, poch_jp1_j = 1L;
  gulong last_k[2] = {0, 0};
  gulong first_term, second_term;
  const gboolean have_sum = ((l-j) % 2 == 0);
  glong wsum;
  mpz_ptr suml;

  if (init_mpz)
  {
    mpz_inits (term, termc, sumc, sum, term_const, four_xn2, qd_jm1, NULL);
    init_mpz = FALSE;
  }

  second_term = (4 - ((l - j) % 4)) % 4;                    /* Analysing k = 2 to get the second signs */
  first_term = (second_term + 1) % 4;                       /* Analysing k = 1 to get the first signs  */
  signs[first_term  % 2] = (first_term  / 2) == 0 ? 1 : -1; /* FIXME (doc this) */
  signs[second_term % 2] = (second_term / 2) == 0 ? 1 : -1; /* sin() */

  mpz_set (four_xn2, mpq_numref (q));
  mpz_mul (four_xn2, four_xn2, four_xn2);
  mpz_mul_ui (four_xn2, four_xn2, 4L);

  mpz_bin_uiui (term, l+j, l-j);
  for (n = 0L; n < j; n++)
    poch_jp1_j *= (j + 1L + n);
  mpz_mul_ui (term, term, poch_jp1_j);

  if (j-1 > 0)
  {
    mpz_set (qd_jm1, mpq_denref(q));
    mpz_pow_ui (qd_jm1, qd_jm1, j-1);
  }

  if (have_sum)
  {
    gulong lmj_2 = (l-j) / 2;
    gulong lpj_2 = (l+j) / 2;
    mpz_bin_uiui (term_const, l + j, lpj_2);

    for (n = 0; n < j; n++)
      mpz_mul_ui (term_const, term_const, lmj_2 + 1 + n);
    mpz_set (sum_sinint, term_const);
    if (lmj_2 % 2 == 1)
      mpz_neg (term_const, term_const);

    mpz_mul (term_const, term_const, mpq_denref (q));
    if (j-1 > 0)
      mpz_mul (term_const, term_const, qd_jm1);
  }
  else
  {
    mpz_set (termc, term);
    mpz_mul_si (termc, termc, (j - l) * (1L + j + l));
    mpz_divexact_ui (termc, termc, (1L + j));
    if (j-1 > 0)
      mpz_divexact (cor, cor, qd_jm1);

    if (((l-j-1) / 2 ) % 2 == 1)
    {
      mpz_neg (termc, termc);
      mpz_neg (cor, cor);
    }

    mpz_set (sumc, termc);
    mpz_mul_2exp (cor, cor, j + 1);

    mpz_add (cor, cor, sumc);
  }

  mpz_mul (term, term, mpq_denref (q));
  if (j-1 > 0)
    mpz_mul (term, term, qd_jm1);
  mpz_mul_2exp (term, term, l - j);

  mpz_set (sum, term);

  mpz_mul_2exp (sum_cos, sum_cos, l - 1);
  mpz_mul_2exp (sum_sin, sum_sin, l - 1);

  wsum = (first_term % 2);
  suml = (wsum == 0) ? sum_sin : sum_cos;
  mpz_mul_2exp (suml, suml, 1);
  mpz_mul (suml, suml, mpq_numref(q));

  (signs[wsum] == 1) ? mpz_sub (suml, suml, sum) : mpz_add (suml, suml, sum);
  if (have_sum)
    (signs[wsum] == -1) ? mpz_sub (suml, suml, term_const) : mpz_add (suml, suml, term_const);
  signs[wsum] *= -1;
  last_k[wsum] = 1;

  for (n = 0; n < l-j-1; n++)
  {
    gulong k = n + 1L;
    wsum = (abs(2 + j - (k + 1) - l) % 2);
    suml = (wsum == 0) ? sum_sin : sum_cos;

    if (have_sum)
    {
      mpz_mul (term_const, term_const, mpq_denref (q));
      mpz_mul_ui (term_const, term_const, 2L * k);
    }
    else
    {
      mpz_mul_si (termc, termc, (j - l + k) * (1L + j + l + k));
      mpz_divexact_ui (termc, termc, (1L + j + k) * (1L + k));
      mpz_mul_ui (termc, termc, k);

      mpz_mul_si (sumc, sumc, (j - l + k) * (1L + j + l + k));
      mpz_divexact_ui (sumc, sumc, (1L + j + k) * (1L + k));
      mpz_mul_ui (sumc, sumc, (k + 1L));
      mpz_add (sumc, sumc, termc);

      mpz_mul_ui (cor, cor, 2L * (k + 1L));
      mpz_add (cor, cor, sumc);
    }

    if (!(converged[0] && converged[1]))
    {
      mpz_mul_si (term, term, (j - l + n) * (1L + j + l + n));
      mpz_divexact_ui (term, term, (1L + j + n) * (1L + n));
      mpz_mul_ui (term, term, k);
      mpz_mul (term, term, mpq_denref (q));

      mpz_mul_ui (sum, sum, 2L * k);
      mpz_mul (sum, sum, mpq_denref (q));

      mpz_add (sum, sum, term);
    }
    else if (have_sum)
      break;

    if (!(converged[0] && converged[1]))
    {
      mpz_mul (suml, suml, four_xn2);
      if ((glong)(mpz_sizeinbase(suml, 2) - mpz_sizeinbase(sum, 2)) > (glong)prec)
        converged[wsum] = TRUE;
      (signs[wsum] == 1) ? mpz_sub (suml, suml, sum) : mpz_add (suml, suml, sum);
      if (have_sum)
        (signs[wsum] == -1) ? mpz_sub (suml, suml, term_const) : mpz_add (suml, suml, term_const);
      signs[wsum] *= -1;
      last_k[wsum] = k + 1;
    }
  }

  if (last_k[0] && last_k[1] && (abs (last_k[0] - last_k[1]) != 1))
    g_error ("Series for sin and cos dont converge together %lu %lu\n", last_k[0], last_k[1]);

  if (last_k[1] > last_k[0])
  {
    mpz_mul_2exp (sum_sin, sum_sin, 1);
    mpz_mul (sum_sin, sum_sin, mpq_numref (q));
    *lk = last_k[1];
  }
  else
  {
    mpz_mul_2exp (sum_cos, sum_cos, 1);
    mpz_mul (sum_cos, sum_cos, mpq_numref (q));
    *lk = last_k[0];
  }
}


static void
_xn_jl_finite (gint l, gint j, mpq_t q, mpfr_t res, mp_rnd_t rnd)
{
  static mpz_t sum_sin, sum_cos, sum_sinint, cor, qd_jm1;
  static gboolean init_mpz = TRUE;
  gulong prec = mpfr_get_prec (res);
  gulong lk = 0;
  mpfr_t cos_x, sin_x, sinint_x, xx, pre_cos, pre_sin, xn_pow, corfact;
  mpfr_inits2 (prec, cos_x, sin_x, sinint_x, xx, pre_cos, pre_sin, xn_pow, corfact, NULL);

  if (init_mpz)
  {
    mpz_inits (sum_sin, sum_cos, sum_sinint, cor, qd_jm1, NULL);
    init_mpz = FALSE;
  }
  else
  {
    mpz_set_ui (sum_sin, 0);
    mpz_set_ui (sum_cos, 0);
    mpz_set_ui (sum_sinint, 0);
    mpz_set_ui (cor, 0);
  }

  mpfr_set_q (xx, q, rnd);
  mpfr_sin_cos (sin_x, cos_x, xx, rnd);

  if (j > 0)
    _xn_jl_finite_x (l, j, q, sum_sin, sum_cos, cor);

  if (l > j)
    _xn_jl_finite_inverse_x (l, j, q, sum_sin, sum_cos, sum_sinint, cor, prec, &lk);
  else
  {
    if (l == j)
    {
      gulong lmj_2 = (l-j) / 2;
      gulong lpj_2 = (l+j) / 2;
      gulong n;
      mpz_bin_uiui (sum_sinint, l + j, lpj_2);
      for (n = 0; n < j; n++)
        mpz_mul_ui (sum_sinint, sum_sinint, lmj_2 + 1 + n);
    }
    if (j-1 > 0)
    {
      mpz_pow_ui (qd_jm1, mpq_denref(q), j-1);
      mpz_divexact (cor, cor, qd_jm1);
    }
    if ((abs((l-j-1) / 2) ) % 2 == 1)
      mpz_neg (cor, cor);

  }

  mpfr_set_z (pre_sin, sum_sin, rnd);
  mpfr_set_z (pre_cos, sum_cos, rnd);

  mpfr_mul (pre_sin, pre_sin, sin_x, rnd);
  mpfr_mul (pre_cos, pre_cos, cos_x, rnd);
  mpfr_add (res, pre_cos, pre_sin, rnd);

  if (l > j)
  {
    mpfr_div_2ui (res, res, l - 1 + lk, rnd);
    mpfr_set_z (xn_pow, mpq_numref (q), rnd);
    mpfr_pow_ui (xn_pow, xn_pow, lk, rnd);
    mpfr_div (res, res, xn_pow, rnd);
  }

  if (j > 1)
  {
    mpfr_set_z (xn_pow, mpq_denref (q), rnd);
    mpfr_pow_ui (xn_pow, xn_pow, j-1, rnd);
    mpfr_div (res, res, xn_pow, rnd);
  }

  if ((l+j) % 2 == 1)
  {
    if (l > j)
    {
      mpfr_fac_ui (corfact, l-j, rnd);
      mpfr_ui_div (corfact, 1, corfact, rnd);
      mpfr_div_2ui (corfact, corfact, l, rnd);
    }
    else
      mpfr_set_ui (corfact, 1, rnd);

    mpfr_mul_z (corfact, corfact, cor, rnd);
    mpfr_sub (res, res, corfact, rnd);
  }
  else if (l >= j)
  {
    ncm_mpsf_sin_int_mpfr (q, sinint_x, rnd);
    mpfr_mul_z (sinint_x, sinint_x, sum_sinint, rnd);
    mpfr_div_2ui (sinint_x, sinint_x, l, rnd);
    mpfr_add (res, res, sinint_x, rnd);
  }

  mpfr_clears (cos_x, sin_x, sinint_x, xx, pre_cos, pre_sin, xn_pow, corfact, NULL);
}


/*************************************************************************************
 * Code to calculate the integral of  x^n j_l, interface function
 *************************************************************************************/

/**
 * ncm_mpsf_sbessel_jl_xj_integral_q: (skip)
 * @l: FIXME
 * @j: FIXME
 * @q: FIXME
 * @res: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integral_q (gint l, gint j, mpq_t q, mpfr_t res, mp_rnd_t rnd)
{
  gboolean taylor;
  gint sign = mpz_sgn (mpq_numref (q));

  g_assert (sign >= 0);

  if (sign == 0)
  {
    mpfr_set_ui (res, 0, rnd);
    return;
  }

  taylor = (mpq_cmp_ui (q, l, 1) <= 0);

  if (taylor)
    _xn_jl_taylor (l, j, q, res, rnd);
  else
    _xn_jl_finite (l, j, q, res, rnd);

  if (FALSE)
    mpfr_printf ("# l %d j %d q %Qd x %.15e res % .36Re | used assym %s\n", l, j, q, mpq_get_d (q), res, !taylor ? "yes" : "no");
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral: (skip)
 * @l: FIXME
 * @j: FIXME
 * @x: FIXME
 * @res: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integral (gint l, gint j, gdouble x, mpfr_t res, mp_rnd_t rnd)
{
  mpq_t q;
  mpq_init (q);
  ncm_rational_coarce_double (x, q);
  ncm_mpsf_sbessel_jl_xj_integral_q (l, j, q, res, rnd);
  mpq_clear (q);
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_a_b: (skip)
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @ki: FIXME
 * @xi: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integral_a_b (NcmMpsfSBesselIntSpline *int_jlspline, guint ki, guint xi, mp_rnd_t rnd)
{
  gulong prec = int_jlspline->prec;
  MPFR_DECL_INIT (w_pow_ip1, prec);
  MPFR_DECL_INIT (ww, prec);
  gint i;
  guint indexi = NCM_MPSF_SBESSEL_INT_MAP (int_jlspline, xi, ki);
  guint indexip1 = NCM_MPSF_SBESSEL_INT_MAP (int_jlspline, xi + 1, ki);

  mpfr_set_q (ww, int_jlspline->k->nodes[ki], rnd);
  mpfr_set (w_pow_ip1, ww, rnd);

  for (i = 0; i < 4; i++)
  {
    mpfr_sub (int_jlspline->rules[i], int_jlspline->xnjlrec[indexip1]->int_jl_xn[i], int_jlspline->xnjlrec[indexi]->int_jl_xn[i], rnd);
    mpfr_div (int_jlspline->rules[i], int_jlspline->rules[i], w_pow_ip1, rnd);

    mpfr_mul (w_pow_ip1, w_pow_ip1, ww, rnd);
  }
}

/**
 * ncm_mpsf_sbessel_djl_xj_integral_a_b: (skip)
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @ki: FIXME
 * @xi: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_djl_xj_integral_a_b (NcmMpsfSBesselIntSpline *int_jlspline, guint ki, guint xi, mp_rnd_t rnd)
{
  gulong prec = int_jlspline->prec;
  MPFR_DECL_INIT (w_pow_ip1, prec);
  MPFR_DECL_INIT (ww, prec);
  gint i;
  guint indexi = NCM_MPSF_SBESSEL_INT_MAP (int_jlspline, xi, ki);
  guint indexip1 = NCM_MPSF_SBESSEL_INT_MAP (int_jlspline, xi + 1, ki);

  mpfr_set_q (ww, int_jlspline->k->nodes[ki], rnd);
  mpfr_set (w_pow_ip1, ww, rnd);

  for (i = 0; i < 4; i++)
  {
    MPFR_DECL_INIT (temp, prec);

    ncm_mpsf_sbessel_jl_xj_integral_recur_calc_djl (int_jlspline->xnjlrec[indexip1], i, int_jlspline->rules[i], rnd);
    ncm_mpsf_sbessel_jl_xj_integral_recur_calc_djl (int_jlspline->xnjlrec[indexi], i, temp, rnd);

    mpfr_sub (int_jlspline->rules[i], int_jlspline->rules[i], temp, rnd);
    mpfr_div (int_jlspline->rules[i], int_jlspline->rules[i], w_pow_ip1, rnd);

    mpfr_mul (w_pow_ip1, w_pow_ip1, ww, rnd);
  }
}

/**
 * ncm_mpsf_sbessel_d2jl_xj_integral_a_b: (skip)
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @ki: FIXME
 * @xi: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_d2jl_xj_integral_a_b (NcmMpsfSBesselIntSpline *int_jlspline, guint ki, guint xi, mp_rnd_t rnd)
{
  gulong prec = int_jlspline->prec;
  MPFR_DECL_INIT (w_pow_ip1, prec);
  MPFR_DECL_INIT (ww, prec);
  gint i;
  guint indexi = NCM_MPSF_SBESSEL_INT_MAP (int_jlspline, xi, ki);
  guint indexip1 = NCM_MPSF_SBESSEL_INT_MAP (int_jlspline, xi + 1, ki);

  mpfr_set_q (ww, int_jlspline->k->nodes[ki], rnd);
  mpfr_set (w_pow_ip1, ww, rnd);

  for (i = 0; i < 4; i++)
  {
    MPFR_DECL_INIT (temp, prec);
    ncm_mpsf_sbessel_jl_xj_integral_recur_calc_d2jl (int_jlspline->xnjlrec[indexip1], i, int_jlspline->rules[i], rnd);
    ncm_mpsf_sbessel_jl_xj_integral_recur_calc_d2jl (int_jlspline->xnjlrec[indexi], i, temp, rnd);

    mpfr_sub (int_jlspline->rules[i], int_jlspline->rules[i], temp, rnd);
    mpfr_div (int_jlspline->rules[i], int_jlspline->rules[i], w_pow_ip1, rnd);

    mpfr_mul (w_pow_ip1, w_pow_ip1, ww, rnd);
  }
}

/**
 * ncm_mpsf_sbessel_jl_xj_integral_a_b_center: (skip)
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @ki: FIXME
 * @xi: FIXME
 * @rnd: FIXME
 *
 * FIXME
*/
void
ncm_mpsf_sbessel_jl_xj_integral_a_b_center (NcmMpsfSBesselIntSpline *int_jlspline, guint ki, guint xi, mp_rnd_t rnd)
{
  MPFR_DECL_INIT (d, int_jlspline->prec);
  mpfr_t *rules = int_jlspline->rules;
  mpfr_t *crules = int_jlspline->crules;

  mpfr_set_q (d, int_jlspline->x->nodes[xi], rnd);


  mpfr_set (crules[0], rules[0], rnd);

  mpfr_mul (crules[1], d, rules[0], rnd);
  mpfr_sub (crules[1], rules[1], crules[1], rnd);

  mpfr_add (crules[2], crules[1], rules[1], rnd);
  mpfr_mul (crules[2], crules[2], d, rnd);
  mpfr_sub (crules[2], rules[2], crules[2], rnd);

  mpfr_fms (crules[3], rules[1], d, crules[2], rnd);
  mpfr_sub (crules[3], crules[3], rules[2], rnd);
  mpfr_sub (crules[3], crules[3], rules[2], rnd);
  mpfr_mul (crules[3], crules[3], d, rnd);
  mpfr_add (crules[3], rules[3], crules[3], rnd);
}

typedef struct
{
  double * c;
  double * g;
  double * diag;
  double * offdiag;
} cspline_state_t;

static inline void
coeff_calc (const double c_array[], double dy, double dx, size_t index, double *b, double *c, double *d)
{
  const double c_i = c_array[index];
  const double c_ip1 = c_array[index + 1];
  *b = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
  *c = c_i;
  *d = (c_ip1 - c_i) / (3.0 * dx);
}

/**
 * ncm_mpsf_sbessel_integrate:
 * @int_jlspline: a #NcmMpsfSBesselIntSpline
 * @s: a #NcmSpline
 * @l: FIXME
 * @ki: FIXME
 * @xi: FIXME
 * @diff: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_mpsf_sbessel_integrate (NcmMpsfSBesselIntSpline *int_jlspline, NcmSpline *s, gint l, guint ki, guint xi, gint diff)
{
  gdouble res = 0.0;
  gdouble a, b;
  gdouble y_lo, y_hi;
  gdouble dx, dy;

  y_lo = ncm_vector_get (s->yv, xi);
  y_hi = ncm_vector_get (s->yv, xi + 1);
  a = ncm_vector_get (s->xv, xi);
  b = ncm_vector_get (s->xv, xi + 1);
  dx = b - a;
  dy = y_hi - y_lo;

  /*printf ("# Preparing the pure integrals %d [%d %.15g]\n", diff, l, w);fflush(stdout);*/
  switch (diff)
  {
    case 0:
      ncm_mpsf_sbessel_jl_xj_integral_a_b (int_jlspline, ki, xi, GMP_RNDN);
      break;
    case 1:
      ncm_mpsf_sbessel_djl_xj_integral_a_b (int_jlspline, ki, xi, GMP_RNDN);
      break;
    case 2:
      ncm_mpsf_sbessel_d2jl_xj_integral_a_b (int_jlspline, ki, xi, GMP_RNDN);
      break;
    default:
      g_error ("Higher derivatives not implemented.");
  }

//printf ("#d%d [%.5f %.5f %d] MP %.15e %.15e %.15e %.15e\n", diff, mpq_get_d (int_jlspline->x->nodes[xi]),
//    mpq_get_d (int_jlspline->k->nodes[ki]), l, mpfr_get_d (int_jlspline->rules[0], GMP_RNDN), mpfr_get_d (int_jlspline->rules[1], GMP_RNDN), mpfr_get_d (int_jlspline->rules[2], GMP_RNDN), mpfr_get_d (int_jlspline->rules[3], GMP_RNDN));

  ncm_mpsf_sbessel_jl_xj_integral_a_b_center (int_jlspline, ki, xi, GMP_RNDN);

  {
		NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
		const cspline_state_t *state = (const cspline_state_t *) (sg->interp)->state;
    double b_i, c_i, d_i;

    coeff_calc(state->c, dy, dx, xi,  &b_i, &c_i, &d_i);

    //printf ("%d %d %d (%.15g %.15g) % .15g % .15g % .15g % .15g\n", l, ki, xi, a, b, y_lo, b_i, c_i, d_i);

    res += y_lo * mpfr_get_d (int_jlspline->crules[0], GMP_RNDD);
    res += b_i  * mpfr_get_d (int_jlspline->crules[1], GMP_RNDD);
    res += c_i  * mpfr_get_d (int_jlspline->crules[2], GMP_RNDD);
    res += d_i  * mpfr_get_d (int_jlspline->crules[3], GMP_RNDD);
  }

  return res;
}
