/***************************************************************************
 *            ncm_lapack.c
 *
 *  Sun March 18 22:33:15 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_lapack
 * @title: Lapack Helper C Functions
 * @short_description: Encapsulated fortran lapack functions
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_lapack.h"
#include "math/ncm_util.h"

#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#ifdef HAVE_MKL_LAPACKE_H
#  include <mkl_lapacke.h>
#elif defined HAVE_CLAPACK_H
#  include <clapack.h>
#endif

#ifdef HAVE_MKL_LAPACK_H
#  include <mkl_lapack.h>
#elif HAVE_LAPACK
void dptsv_ (gint *N, gint *NRHS, gdouble *d, gdouble *e, gdouble *b, gint *ldb, gint *info);
void dpotrf_ (const char *uplo, const gint *n, double *a, const gint *lda, gint *info);
#endif /* HAVE_LAPACK */

/**
 * ncm_lapack_dptsv:
 * @d: FIXME
 * @e: FIXME
 * @b: FIXME
 * @x: FIXME
 * @size: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gint
ncm_lapack_dptsv (gdouble *d, gdouble *e, gdouble *b, gdouble *x, guint size)
{
#ifdef HAVE_LAPACKE
  lapack_int info = LAPACKE_dptsv (LAPACK_ROW_MAJOR, size, 1, d, e, b, 1);
  if (x != b)
    memcpy (x, b, sizeof (gdouble) * size);
  return info;
#elif defined HAVE_LAPACK
	gint N = size;
	gint NRHS = 1;
	gint LDB = N;
	gint info;
	dptsv_ (&N, &NRHS, d, e, b, &LDB, &info);
  if (x != b)
    memcpy (x, b, sizeof (gdouble) * size);
	return info;
#else
  gsl_vector_view b_vec       = gsl_vector_view_array (b, size);
  gsl_vector_view diag_vec    = gsl_vector_view_array (d, size);
  gsl_vector_view offdiag_vec = gsl_vector_view_array (e, size - 1);
  gsl_vector_view x_vec       = gsl_vector_view_array (x, size);
  gint status = gsl_linalg_solve_symm_tridiag (&diag_vec.vector,
                                               &offdiag_vec.vector,
                                               &b_vec.vector,
                                               &x_vec.vector);
  NCM_TEST_GSL_RESULT ("ncm_lapack_dptsv[gsl_linalg_solve_symm_tridiag]", status);
  return status; /* THAT'S NOT OK FIXME */
#endif /* HAVE_LAPACK */
}

/**
 * ncm_lapack_dpotrf:
 * @uplo: FIXME
 * @size: FIXME
 * @a: FIXME
 * @lda: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gint 
ncm_lapack_dpotrf (gchar uplo, guint size, gdouble *a, guint lda)
{
#ifdef HAVE_LAPACKE
  lapack_int info = LAPACKE_dpotrf (LAPACK_ROW_MAJOR, uplo, size, a, lda);
  return info;
#elif defined HAVE_CLAPACK
  gint ret = clapack_dpotrf (CblasRowMajor, 
                             uplo == 'U' ? CblasUpper : CblasLower, 
                             size, a, size);
  if (ret < 0)
    g_error ("ncm_lapack_dpotrf: invalid parameter %d", -ret);
  else if (ret > 0)
    g_error ("ncm_lapack_dpotrf: the leading minor of order %d is not positive definite", ret);
  return ret;
#elif defined HAVE_LAPACK
  gint info = 0;
  gint n = size;
  gint LDA = lda;
  gchar UPLO = uplo == 'L' ? 'U' : 'L';
  dpotrf_ (&UPLO, &n, a, &LDA, &info);  
  return info;
#else /* Fall back to gsl cholesky */
  gint ret;
  gsl_matrix_view mv = gsl_matrix_view_array_with_tda (a, size, size, lda);
  ret = gsl_linalg_cholesky_decomp (&mv.matrix);
  NCM_TEST_GSL_RESULT("gsl_linalg_cholesky_decomp", ret);
  return ret; /* THAT'S NOT OK FIXME */
#endif
}
