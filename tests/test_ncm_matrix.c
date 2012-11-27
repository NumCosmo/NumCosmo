/***************************************************************************
 *            test_ncm_matrix.c
 *
 *  Sat April 21 14:30:26 2012
 *  Copyright  2012  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

void test_ncm_matrix_new (void);
void test_ncm_matrix_new_gsl (void);
void test_ncm_matrix_new_array (void);
void test_ncm_matrix_new_data_slice (void);
void test_ncm_matrix_new_data_malloc (void);
void test_ncm_matrix_new_data_static (void);
void test_ncm_matrix_new_data_static_tda (void);
void test_ncm_matrix_operations (void);
void test_ncm_matrix_free (void);
void test_ncm_matrix_submatrix (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/numcosmo/ncm_matrix/new", &test_ncm_matrix_new);
  g_test_add_func ("/numcosmo/ncm_matrix/new_gsl", &test_ncm_matrix_new_gsl);
  g_test_add_func ("/numcosmo/ncm_matrix/new_array", &test_ncm_matrix_new_array);
  g_test_add_func ("/numcosmo/ncm_matrix/new_data_slice", &test_ncm_matrix_new_data_slice);
  g_test_add_func ("/numcosmo/ncm_matrix/new_data_malloc", &test_ncm_matrix_new_data_malloc);
  g_test_add_func ("/numcosmo/ncm_matrix/new_data_static", &test_ncm_matrix_new_data_static);
  g_test_add_func ("/numcosmo/ncm_matrix/new_data_static_tda", &test_ncm_matrix_new_data_static_tda);
  g_test_add_func ("/numcosmo/ncm_matrix/operations", &test_ncm_matrix_operations);
  g_test_add_func ("/numcosmo/ncm_matrix/submatrix", &test_ncm_matrix_submatrix);
  g_test_add_func ("/numcosmo/ncm_matrix/free", &test_ncm_matrix_free);

  g_test_run ();
}

static NcmMatrix *m;

#define _NCM_MATRIX_TEST_NROW 10
#define _NCM_MATRIX_TEST_NCOL 7

void
test_ncm_matrix_new_sanity (NcmMatrix *cm)
{
  guint i;

  g_assert (NCM_IS_MATRIX (cm));
  g_assert (g_object_is_floating (cm));
  ncm_matrix_ref (cm);
  g_assert (!g_object_is_floating (cm));

  for (i = 0; i < 10 * _NCM_MATRIX_TEST_NROW; i++)
  {
	const guint nr = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NROW - 1);
	const guint nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 1);
	const gdouble d = g_test_rand_double ();
	ncm_matrix_set (cm, nr, nc, d);
	g_assert (ncm_matrix_get (cm, nr, nc) == d);
  }
}

void
test_ncm_matrix_new (void)
{
  m = ncm_matrix_new (_NCM_MATRIX_TEST_NROW, _NCM_MATRIX_TEST_NCOL);
  test_ncm_matrix_new_sanity (m);
}

void
test_ncm_matrix_operations (void)
{
  NcmMatrix *cm = ncm_matrix_dup (m);
  guint i;
  test_ncm_matrix_new_sanity (cm);

  for (i = 0; i < 10 * _NCM_MATRIX_TEST_NROW; i++)
  {
	const guint nr = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NROW - 1);
	const guint nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 1);
	const gdouble d = g_test_rand_double ();
	ncm_matrix_set (m, nr, nc, d);
	g_assert (ncm_matrix_get (m, nr, nc) == d);
  }

  for (i = 0; i < 10 * _NCM_MATRIX_TEST_NROW; i++)
  {
	guint nr, nc;
	gdouble *d;

	nr = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NROW - 1);
	nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 1);
	d = ncm_matrix_ptr (m, nr, nc);
	(*d) *= g_test_rand_double ();

	g_assert (ncm_matrix_get (m, nr, nc) == *d);
  }

  {
	NcmMatrix *subm = ncm_matrix_get_submatrix (m, 3, 0, _NCM_MATRIX_TEST_NROW - 3, _NCM_MATRIX_TEST_NCOL);
	NcmMatrix *sm_cp = ncm_matrix_dup (subm);
	ncm_matrix_transpose (subm);
	for (i = 0; i < 10 * _NCM_MATRIX_TEST_NROW; i++)
	{
	  guint nr, nc;
	  gdouble d;

	  nr = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NROW - 4);
	  nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 1);
	  d = ncm_matrix_get (sm_cp, nr, nc);

	  g_assert (ncm_matrix_get (subm, nc, nr) == d);
	}
	ncm_matrix_free (sm_cp);
	ncm_matrix_free (subm);
  }

  for (i = 0; i < 10 * _NCM_MATRIX_TEST_NROW; i++)
  {
	guint nr = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NROW - 1);
	guint nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 1);

	ncm_matrix_set_identity (m);
	if (nr == nc)
	  g_assert (ncm_matrix_get (m, nr, nc) == 1.0);
	else
	  g_assert (ncm_matrix_get (m, nr, nc) == 0.0);
  }

  {
	ncm_matrix_set_zero (m);
	for (i = 0; i < 10 * _NCM_MATRIX_TEST_NROW; i++)
	{
	  guint nr = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NROW - 1);
	  guint nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 1);

	  g_assert (ncm_matrix_get (m, nr, nc) == 0.0);
	}
  }

  for (i = 0; i < 10 * _NCM_MATRIX_TEST_NROW; i++)
  {
	guint nr, nc;
	gdouble d, d1 = g_test_rand_double ();

	nr = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NROW - 1);
	nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 1);
	d = ncm_matrix_get (m, nr, nc);
	ncm_matrix_scale (m, d1);
	g_assert (ncm_matrix_get (m, nr, nc) == d * d1);
  }

  {
	NcmVector *v = ncm_vector_new (_NCM_MATRIX_TEST_NROW);
	for (i = 0; i < _NCM_MATRIX_TEST_NROW; i++)
	{
	  ncm_vector_set (v, i, 8.2 + i);
	}
	guint nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 1);
	ncm_matrix_set_col (m, nc, v);
	for (i = 0; i < _NCM_MATRIX_TEST_NROW; i++)
	  g_assert (ncm_matrix_get (m, i, nc) == ncm_vector_get (v, i));
  }

  ncm_matrix_memcpy (m, cm);
  for (i = 0; i < _NCM_MATRIX_TEST_NROW; i++)
  {
	guint j;
	for (j = 0; j < _NCM_MATRIX_TEST_NCOL; j++)
	  g_assert (ncm_matrix_get (m, i, j) == ncm_matrix_get (cm, i, j));
  }

}

void
test_ncm_matrix_new_gsl (void)
{
  gsl_matrix *gm = gsl_matrix_alloc (_NCM_MATRIX_TEST_NROW, _NCM_MATRIX_TEST_NCOL);
  NcmMatrix *mm = ncm_matrix_new_gsl (gm);
  test_ncm_matrix_new_sanity (mm);
  g_assert (NCM_MATRIX_NROWS (mm) == gm->size1 && NCM_MATRIX_NCOLS (mm) == gm->size2);

  ncm_matrix_free (mm);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	ncm_matrix_free (mm);
	exit (0);
  }
  g_test_trap_assert_failed ();

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	g_assert (mm->gm == NULL);
	exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_matrix_new_array (void)
{
  GArray *ga = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL);
  NcmMatrix *mm;
  g_array_set_size (ga, _NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL);
  mm = ncm_matrix_new_array (ga, _NCM_MATRIX_TEST_NCOL);
  test_ncm_matrix_new_sanity (mm);
  g_array_unref (ga);

  g_assert (NCM_MATRIX_NROWS (mm) == ga->len / _NCM_MATRIX_TEST_NCOL);

  g_assert (ga == ncm_matrix_get_array (mm));
  g_array_unref (ga);

  ncm_matrix_free (mm);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	ncm_matrix_free (mm);
	exit (0);
  }
  g_test_trap_assert_failed ();

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	g_array_unref (mm->a);
	exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_matrix_new_data_slice (void)
{
  gdouble *d = g_slice_alloc (_NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL * sizeof (gdouble));
  NcmMatrix *mm;
  mm = ncm_matrix_new_data_slice (d, _NCM_MATRIX_TEST_NROW, _NCM_MATRIX_TEST_NCOL);
  test_ncm_matrix_new_sanity (mm);

  g_assert ((NCM_MATRIX_NROWS (mm) * NCM_MATRIX_NCOLS (mm)) == (_NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL));

  ncm_matrix_free (mm);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	ncm_matrix_free (mm);
	exit (0);
  }
  g_test_trap_assert_failed ();

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	g_assert (NCM_MATRIX_DATA (mm) == NULL);
	exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_matrix_new_data_malloc (void)
{
  gdouble *d = g_malloc ( _NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL * sizeof (gdouble));
  NcmMatrix *mm;
  mm = ncm_matrix_new_data_malloc (d, _NCM_MATRIX_TEST_NROW, _NCM_MATRIX_TEST_NCOL);
  test_ncm_matrix_new_sanity (mm);

  g_assert ((NCM_MATRIX_NROWS (mm) * NCM_MATRIX_NCOLS (mm)) == (_NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL));

  ncm_matrix_free (mm);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	ncm_matrix_free (mm);
	exit (0);
  }
  g_test_trap_assert_failed ();

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	g_assert (NCM_MATRIX_DATA (mm) == NULL);
	exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_matrix_new_data_static (void)
{
  gdouble d[_NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL];
  NcmMatrix *mm;
  mm = ncm_matrix_new_data_static (d, _NCM_MATRIX_TEST_NROW, _NCM_MATRIX_TEST_NCOL);
  test_ncm_matrix_new_sanity (mm);

  g_assert ((NCM_MATRIX_NROWS (mm) * NCM_MATRIX_NCOLS (mm)) == (_NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL));

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	ncm_matrix_free (mm);
	exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_matrix_new_data_static_tda (void)
{
  gdouble d[_NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL * 2];
  NcmMatrix *mm;
  guint i;
  mm = ncm_matrix_new_data_static_tda (d, _NCM_MATRIX_TEST_NROW, _NCM_MATRIX_TEST_NCOL, _NCM_MATRIX_TEST_NCOL * 2);
  test_ncm_matrix_new_sanity (mm);

  g_assert ((NCM_MATRIX_NROWS (mm) * NCM_MATRIX_NCOLS (mm)) == (_NCM_MATRIX_TEST_NROW * _NCM_MATRIX_TEST_NCOL));

  for (i = 0; i < 10 * _NCM_MATRIX_TEST_NROW; i++)
  {
	const guint nr = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NROW - 1);
	const guint nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 1);
	const gdouble val = g_test_rand_double ();
	ncm_matrix_set (mm, nr, nc, val);
	g_assert (ncm_matrix_get (mm, nr, nc) == val);
  }

  {
	NcmVector *v = ncm_matrix_get_col (mm, _NCM_MATRIX_TEST_NCOL - 1);
	g_assert (ncm_vector_len (v) == NCM_MATRIX_NROWS (mm));
	for (i = 0; i < _NCM_MATRIX_TEST_NROW; i++)
	  g_assert (ncm_vector_get (v, i) == ncm_matrix_get (mm, i, _NCM_MATRIX_TEST_NCOL - 1));
	ncm_vector_free (v);
  }

  {
	NcmVector *v = ncm_matrix_get_row (mm, _NCM_MATRIX_TEST_NROW - 1);
	g_assert (ncm_vector_len (v) == NCM_MATRIX_NCOLS (mm));
	for (i = 0; i < _NCM_MATRIX_TEST_NCOL; i++)
	  g_assert (ncm_vector_get (v, i) == ncm_matrix_get (mm, _NCM_MATRIX_TEST_NROW - 1, i));
	ncm_vector_free (v);
  }

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	ncm_matrix_free (mm);
	exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_matrix_submatrix (void)
{
  NcmMatrix *sm = ncm_matrix_get_submatrix (m, 5, 3, _NCM_MATRIX_TEST_NROW - 5, _NCM_MATRIX_TEST_NCOL - 3);
  guint ntests = 20 * _NCM_MATRIX_TEST_NROW;

  g_assert (NCM_MATRIX_NROWS (sm) == (_NCM_MATRIX_TEST_NROW - 5) && NCM_MATRIX_NCOLS (sm) == (_NCM_MATRIX_TEST_NCOL - 3));

  while (ntests--)
  {
	guint nr = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NROW - 5);
	guint nc = g_test_rand_int_range (0, _NCM_MATRIX_TEST_NCOL - 3);
	ncm_matrix_set (sm, nr, nc, g_test_rand_double ());
	g_assert_cmpfloat (ncm_matrix_get (sm, nr, nc), ==, ncm_matrix_get (m, nr + 5, nc + 3));
  }

  g_assert (NCM_IS_MATRIX (sm));
  g_assert (g_object_is_floating (sm));
  ncm_matrix_ref (sm);
  g_assert (!g_object_is_floating (sm));

  ncm_matrix_free (m);
  g_assert (G_IS_OBJECT (m));
  ncm_matrix_ref (m);

  ncm_matrix_free (sm);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	ncm_matrix_free (sm);
	exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_matrix_free (void)
{
  ncm_matrix_free (m);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
	ncm_matrix_free (m);
	exit (0);
  }
  g_test_trap_assert_failed ();
}
