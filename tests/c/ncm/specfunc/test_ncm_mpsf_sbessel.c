/***************************************************************************
 *            test_ncm_mpsf_sbessel.c
 *
 *  Sat September 26 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mpsf_sbessel.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_sf_bessel.h>
#include <mpfr.h>

#define REF_PREC 1024
#define TEST_PREC 256

/* j_l(x) at REF_PREC bits by upward recurrence from j_0 = sin x / x and
 * j_1 = sin x / x^2 - cos x / x, independent of both branches of ncm_mpsf_sbessel. */
static void
_ref_sbessel (gulong l, mpq_t q, mpfr_t res)
{
  mpfr_t x, jm, j0, jp, s, c;
  gulong n;

  mpfr_inits2 (REF_PREC, x, jm, j0, jp, s, c, (mpfr_ptr) NULL);
  mpfr_set_q (x, q, MPFR_RNDN);
  mpfr_sin_cos (s, c, x, MPFR_RNDN);

  mpfr_div (jm, s, x, MPFR_RNDN);
  mpfr_div (j0, jm, x, MPFR_RNDN);
  mpfr_div (jp, c, x, MPFR_RNDN);
  mpfr_sub (j0, j0, jp, MPFR_RNDN);

  if (l == 0)
  {
    mpfr_set (res, jm, MPFR_RNDN);
  }
  else
  {
    for (n = 1; n < l; n++)
    {
      /* j_{n+1} = (2n + 1) / x j_n - j_{n-1} */
      mpfr_mul_ui (jp, j0, 2 * n + 1, MPFR_RNDN);
      mpfr_div (jp, jp, x, MPFR_RNDN);
      mpfr_sub (jp, jp, jm, MPFR_RNDN);
      mpfr_swap (jm, j0);
      mpfr_swap (j0, jp);
    }

    mpfr_set (res, j0, MPFR_RNDN);
  }

  mpfr_clears (x, jm, j0, jp, s, c, (mpfr_ptr) NULL);
}

/* Both branches, Taylor for |x| < l and the closed form for |x| >= l, agree with the
 * reference to the precision of the result (measured within one bit of 2^-256). */
static void
test_ncm_mpsf_sbessel_branches (void)
{
  const struct
  {
    gulong l;
    glong num;
    gulong den;
  } cases[] = {
    {0, 3, 10}, {0, 29, 4}, {1, 5, 2}, {1, -29, 4}, {10, 19, 2}, {10, 21, 2}, {10, 3, 1}, {25, 99, 4}, {25, 401, 4},
  };

  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    mpq_t q;
    mpfr_t got, ref, diff;

    mpq_init (q);
    mpq_set_si (q, cases[i].num, cases[i].den);
    mpfr_init2 (got, TEST_PREC);
    mpfr_init2 (ref, REF_PREC);
    mpfr_init2 (diff, REF_PREC);

    ncm_mpsf_sbessel (cases[i].l, q, got, MPFR_RNDN);
    _ref_sbessel (cases[i].l, q, ref);

    mpfr_sub (diff, got, ref, MPFR_RNDN);
    mpfr_div (diff, diff, ref, MPFR_RNDN);
    g_assert_cmpfloat (fabs (mpfr_get_d (diff, MPFR_RNDN)), <, ldexp (1.0, -TEST_PREC + 4));

    mpfr_clears (got, ref, diff, (mpfr_ptr) NULL);
    mpq_clear (q);
  }

  ncm_mpsf_sbessel_free_cache ();
}

/* The double-argument form agrees with GSL (measured 1.1e-15 at worst). */
static void
test_ncm_mpsf_sbessel_gsl (void)
{
  const gdouble xs[] = {0.1, 1.0, 3.7, 12.5, 40.0};
  const gulong ls[]  = {0, 1, 5, 20};
  guint i, j;

  for (j = 0; j < G_N_ELEMENTS (ls); j++)
  {
    for (i = 0; i < G_N_ELEMENTS (xs); i++)
    {
      MPFR_DECL_INIT (res, 53);

      ncm_mpsf_sbessel_d (ls[j], xs[i], res, MPFR_RNDN);
      ncm_assert_cmpdouble_e (mpfr_get_d (res, MPFR_RNDN), ==, gsl_sf_bessel_jl (ls[j], xs[i]), 1.0e-14, 1.0e-300);
    }
  }

  ncm_mpsf_sbessel_free_cache ();
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/mpsf/sbessel/branches", &test_ncm_mpsf_sbessel_branches);
  g_test_add_func ("/ncm/mpsf/sbessel/gsl", &test_ncm_mpsf_sbessel_gsl);

  g_test_run ();
}

