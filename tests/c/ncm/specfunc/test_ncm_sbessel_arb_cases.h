/***************************************************************************
 *            test_ncm_sbessel_arb_cases.h
 *
 *  Wed October 14 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_sbessel_arb_cases.h
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

/*
 * The rows of the certified tables printed by tests/tools/make_sbessel_arb_cases.py:
 * I = int_a^b f(x) j_ell^(d)(k x) dx from nc_xcor_kernel_analytic_arb, relative radius
 * 1e-25, rounded to double. f is written here exactly as the tool writes it: TOPHAT 1;
 * POWER x^n, n = p1; COS cos (p1 x); EXP exp (-x / p1); RATIONAL
 * x^2 / (1 + ((x - p1) / p2)^2)^3; GAUSS exp (-((x - p1) / p2)^2 / 2). The rows of one
 * case are contiguous, multipoles in increasing order. Each test holds its own table.
 */

#ifndef _TEST_NCM_SBESSEL_ARB_CASES_H_
#define _TEST_NCM_SBESSEL_ARB_CASES_H_

#include <numcosmo/numcosmo.h>
#include <string.h>

G_BEGIN_DECLS

typedef enum _TestArbShape
{
  TEST_ARB_TOPHAT,
  TEST_ARB_POWER,
  TEST_ARB_COS,
  TEST_ARB_EXP,
  TEST_ARB_RATIONAL,
  TEST_ARB_GAUSS,
} TestArbShape;

typedef struct _TestArbCell
{
  const gchar *name;
  TestArbShape shape;
  gdouble a;
  gdouble b;
  gdouble p1;
  gdouble p2;
  guint deriv;
  guint ell;
  gdouble k;
  gdouble value;
} TestArbCell;

/* The first row of case @name among @n_cells rows and, in @n, how many rows it has */
static inline const TestArbCell *
test_arb_cells_find (const TestArbCell *cells, guint n_cells, const gchar *name, guint *n)
{
  const TestArbCell *first = NULL;
  guint i;

  *n = 0;

  for (i = 0; i < n_cells; i++)
  {
    if (strcmp (cells[i].name, name) == 0)
    {
      if (first == NULL)
        first = &cells[i];

      (*n)++;
    }
  }

  g_assert_nonnull (first);

  return first;
}

/* f (x) of the row's shape */
static inline gdouble
test_arb_cell_f (const TestArbCell *c, gdouble x)
{
  switch (c->shape)
  {
    case TEST_ARB_TOPHAT:
      return 1.0;

    case TEST_ARB_POWER:
      return gsl_pow_int (x, (gint) c->p1);

    case TEST_ARB_COS:
      return cos (c->p1 * x);

    case TEST_ARB_EXP:
      return exp (-x / c->p1);

    case TEST_ARB_RATIONAL:
      return x * x / gsl_pow_3 (1.0 + gsl_pow_2 ((x - c->p1) / c->p2));

    case TEST_ARB_GAUSS:
      return exp (-0.5 * gsl_pow_2 ((x - c->p1) / c->p2));

    default:
      g_assert_not_reached ();

      return 0.0;
  }
}

G_END_DECLS

#endif /* _TEST_NCM_SBESSEL_ARB_CASES_H_ */

