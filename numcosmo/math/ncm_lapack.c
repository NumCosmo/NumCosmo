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
 * @title: NcmLapack
 * @short_description: Encapsulated LAPACK functions. 
 *
 * This object is dedicated to encapsulate functions from <ulink url="http://www.netlib.org/lapack/">LAPACK</ulink> choosing the most suitable backend.
 * 
 * Priority order: (1) LAPACKE, (2) CLAPACK, (3) LAPACK and (4) GSL.
 * 
 * The description of each function follows its respective LAPACK documentation.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_lapack.h"
#include "math/ncm_matrix.h"
#include "math/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <string.h>
#include <gsl/gsl_vector.h>

#ifdef HAVE_BLAS
#  ifdef HAVE_MKL_CBLAS_H
#    include <mkl_cblas.h>
#  elif defined (HAVE_CBLAS_H)
#    include <cblas.h>
#  else
#    include <gsl/gsl_cblas.h>
#    include <gsl/gsl_linalg.h>
#  endif
#endif

#ifdef HAVE_MKL_LAPACKE_H
#  include <mkl_lapacke.h>
#elif defined HAVE_LAPACKE
#  include <lapacke.h>
#elif defined HAVE_CLAPACK_H
#  include <clapack.h>
#endif

#ifdef HAVE_MKL_LAPACK_H
#  include <mkl_lapack.h>
#elif HAVE_LAPACKE
#elif HAVE_LAPACK
void dptsv_ (gint *N, gint *NRHS, gdouble *d, gdouble *e, gdouble *b, gint *ldb, gint *info);
void dpotrf_ (const char *uplo, const gint *n, double *a, const gint *lda, gint *info);
void dpotri_ (const char *uplo, const gint *n, double *a, const gint *lda, gint *info);
void dpotrs_ (const char *uplo, const gint *n, const gint *nrhs, double *a, const gint *lda, double *b, const gint *ldb, gint *info);
void dposv_ (const char *uplo, const gint *n, const gint *nrhs, double *a, const gint *lda, double *b, const gint *ldb, gint *info);
void dggglm_ (const gint *N, const gint *M, const gint *P, gdouble *X, const gint *LDA, gdouble *L, const gint *LDB, 
              gdouble *d, gdouble *p, gdouble *y, gdouble *work, const gint *lwork, gint *info);
#endif /* HAVE_LAPACK */
#endif /* NUMCOSMO_GIR_SCAN */

#define _NCM_LAPACK_CONV_UPLO(uplo) (uplo == 'L' ? 'U' : 'L')
#define _NCM_LAPACK_CONV_TRANS(trans) (trans == 'N' ? 'T' : 'N')

/**
 * ncm_lapack_dptsv:
 * @d: array of doubles with dimension @size
 * @e: array of doubles with dimension @size -1
 * @b: array of doubles with dimension @size
 * @x: array of doubles with dimension @size
 * @size: The order of the matrix $A$ (>= 0)
 *
 * This function computes the solution to a real system of linear equations
 * $A*X = B$ (B = @b), where $A$ is an N-by-N (N = @size) symmetric positive definite tridiagonal
 * matrix, and $X$ and $B$ are N-by-NRHS (NRHS = 1) matrices.
 *
 * $A$ is factored as $A = L*D*L^T$, and the factored form of $A$ is then
 * used to solve the system of equations.
 *
 * Returns: i = 0:  successful exit
 * 
 *        < 0:  -i, the i-th argument had an illegal value
 * 
 *        > 0:   i, the leading minor of order i is not
 *               positive definite, and the solution has not been
 *               computed.  The factorization has not been completed
 *               unless i = N.
 */
gint
ncm_lapack_dptsv (gdouble *d, gdouble *e, gdouble *b, gdouble *x, guint size)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
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
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @size: The order of the matrix @a, @size >= 0
 * @a: array of doubles with dimension (@size, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1,@size)
 *
 * This function computes the Cholesky factorization of a real symmetric
 * positive definite matrix @a.
 * 
 * The factorization has the form
 * $A = U^T * U$, if @uplo = 'U', or
 * $A = L  * L^T$, if @uplo = 'L',
 * where A = @a, $U$ is an upper triangular matrix and $L$ is lower triangular.
 *
 * Returns: i = 0:  successful exit
 * 
 *          < 0:  -i, the i-th argument had an illegal value
 * 
 *          > 0:   i, the leading minor of order i is not
 *                positive definite, and the factorization could not be
 *                completed.
 */
gint 
ncm_lapack_dpotrf (gchar uplo, guint size, gdouble *a, guint lda)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dpotrf (LAPACK_ROW_MAJOR, uplo, size, a, lda);
  return info;
#elif defined (HAVE_CLAPACK) && defined (NUMCOSMO_PREFER_LAPACKE)
  gint ret = clapack_dpotrf (CblasRowMajor, 
                             uplo == 'U' ? CblasUpper : CblasLower, 
                             size, a, lda);
  return ret;
#elif defined HAVE_LAPACK
  gint info = 0;
  gint n = size;
  gint LDA = lda;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);
  dpotrf_ (&uplo, &n, a, &LDA, &info);  
  return info;
#else /* Fall back to gsl cholesky */
  gint ret;
  gsl_matrix_view mv = gsl_matrix_view_array_with_tda (a, size, size, lda);
  ret = gsl_linalg_cholesky_decomp (&mv.matrix);
  NCM_TEST_GSL_RESULT ("gsl_linalg_cholesky_decomp", ret);
  return ret; /* THAT'S NOT OK FIXME */
#endif
}

/**
 * ncm_lapack_dpotri:
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @size: The order of the matrix @a, @size >= 0
 * @a: array of doubles with dimension (@size, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1,@size)
 *
 * This function computes the inverse of a real symmetric positive
 * definite matrix @a = A using the Cholesky factorization 
 * $A = U^T*U$ or $A = L*L^T$ computed by ncm_lapack_dpotrf().
 *
 * Returns: i = 0:  successful exit
 * 
 *          < 0:  -i, the i-th argument had an illegal value
 *
 *          > 0: the (i,i) element of the factor U
 *            or L is zero, and the inverse could not be computed.
 */
gint 
ncm_lapack_dpotri (gchar uplo, guint size, gdouble *a, guint lda)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dpotri (LAPACK_ROW_MAJOR, uplo, size, a, lda);
  return info;
#elif defined (HAVE_CLAPACK) && defined (NUMCOSMO_PREFER_LAPACKE)
  gint ret = clapack_dpotri (CblasRowMajor, 
                             uplo == 'U' ? CblasUpper : CblasLower, 
                             size, a, lda);
  return ret;
#elif defined HAVE_LAPACK
  gint info = 0;
  gint n = size;
  gint LDA = lda;
  
  uplo = _NCM_LAPACK_CONV_UPLO (uplo);
  dpotri_ (&uplo, &n, a, &LDA, &info);  
  return info;
#else /* Fall back to gsl cholesky */
  gint ret;
  gsl_matrix_view mv = gsl_matrix_view_array_with_tda (a, size, size, lda);
  ret = gsl_linalg_cholesky_invert (&mv.matrix);
  NCM_TEST_GSL_RESULT ("gsl_linalg_cholesky_decomp", ret);
  return ret; /* THAT'S NOT OK FIXME */
#endif
}

/**
 * ncm_lapack_dpotrs:
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @size: The order of the matrix @a, @size >= 0
 * @nrhs: Number of right-hand-side vectors to solve
 * @a: array of doubles with dimension (@size, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1, @size)
 * @b: array of doubles with dimension (@size, @ldb)
 * @ldb: The leading dimension of the array @b, @ldb >= max (1, @size)
 *
 * This function computes the solution of $A X = B$ for a real symmetric positive
 * definite matrix @a = A using the Cholesky factorization $A = U^T*U$ or $A = L*L^T$
 * already performed by ncm_lapack_dpotrf().
 * On entry @b contain the vectors $B$ and on exit @b contain the solutions if the return
 * is 0.
 *
 * Returns: i = 0:  successful exit
 * 
 *          < 0:  -i, the i-th argument had an illegal value
 *
 *          > 0: the (i,i) element of the factor U
 *            or L is zero, and the inverse could not be computed.
 */
gint 
ncm_lapack_dpotrs (gchar uplo, guint size, guint nrhs, gdouble *a, guint lda, gdouble *b, guint ldb)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dpotrs (LAPACK_ROW_MAJOR, uplo, size, nrhs, a, lda, b, ldb);
  return info;
#elif defined (HAVE_CLAPACK) && defined (NUMCOSMO_PREFER_LAPACKE)
	gint ret = clapack_dpotrs (CblasRowMajor, 
	                          uplo == 'U' ? CblasUpper : CblasLower, 
	                          size, a, lda, b, ldb);
  return ret;
#elif defined HAVE_LAPACK
  gint info = 0;
  gint n    = size;
  gint NRHS = nrhs;
  gint LDA  = lda;
  gint LDB  = ldb;
  
  uplo = _NCM_LAPACK_CONV_UPLO (uplo);
  dpotrs_ (&uplo, &n, &NRHS, a, &LDA, b, &LDB, &info);  
  return info;
#else /* No fall back */
	g_error ("ncm_lapack_dpotrs: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dposv:
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @size: The order of the matrix @a, @size >= 0
 * @nrhs: Number of right-hand-side vectors to solve
 * @a: array of doubles with dimension (@size, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1, @size)
 * @b: array of doubles with dimension (@size, @ldb)
 * @ldb: The leading dimension of the array @b, @ldb >= max (1, @size)
 *
 * This function computes the solution of $A X = B$ for a real symmetric positive
 * definite matrix @a = A using the Cholesky factorization $A = U^T*U$ or $A = L*L^T$.
 * On entry @b contain the vectors $B$ and on exit @b contain the solutions if the return
 * is 0.
 *
 * Returns: i = 0:  successful exit
 * 
 *          < 0:  -i, the i-th argument had an illegal value
 *
 *          > 0: the (i,i) element of the factor U
 *            or L is zero, and the inverse could not be computed.
 */
gint 
ncm_lapack_dposv (gchar uplo, guint size, guint nrhs, gdouble *a, guint lda, gdouble *b, guint ldb)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dposv (LAPACK_ROW_MAJOR, uplo, size, nrhs, a, lda, b, ldb);
  return info;
#elif defined (HAVE_CLAPACK) && defined (NUMCOSMO_PREFER_LAPACKE)
	gint ret = clapack_dposv (CblasRowMajor, 
	                          uplo == 'U' ? CblasUpper : CblasLower, 
	                          size, a, lda, b, ldb);
  return ret;
#elif defined HAVE_LAPACK
  gint info = 0;
  gint n    = size;
  gint NRHS = nrhs;
  gint LDA  = lda;
  gint LDB  = ldb;
  
  uplo = _NCM_LAPACK_CONV_UPLO (uplo);
  dposv_ (&uplo, &n, &NRHS, a, &LDA, b, &LDB, &info);  
  return info;
#else /* No fall back */
	g_error ("ncm_lapack_dposv: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dggglm_alloc:
 * @L: a #NcmMatrix
 * @X: a #NcmMatrix
 * @p: a #NcmVector
 * @d: a #NcmVector
 * @y: a #NcmVector
 * 
 * Calculates and allocs memory to solve the system
 * determined by the parameters.
 * 
 * This function is expect the matrix @X and @L to be row-major.
 * 
 * Returns: (transfer full) (array) (element-type double): the newly allocated workspace
 */
GArray *
ncm_lapack_dggglm_alloc (NcmMatrix *L, NcmMatrix *X, NcmVector *p, NcmVector *d, NcmVector *y)
{
#ifdef HAVE_LAPACK
  gint N   = ncm_matrix_nrows (L);
  gint M   = ncm_matrix_ncols (X);
  gint P   = ncm_matrix_ncols (L);
  gint LDA = N;
  gint LDB = N;
  gdouble work;
  gint lwork = -1;
  gint info = 0;
  
  g_assert_cmpint (N, ==, ncm_matrix_nrows (X));
  g_assert_cmpint (N, ==, ncm_vector_len (d));
  g_assert_cmpint (M, ==, ncm_vector_len (p));
  g_assert_cmpint (P, ==, ncm_vector_len (y));
  
  dggglm_ (&N, &M, &P, 
           ncm_matrix_data (X), 
           &LDA, 
           ncm_matrix_data (L), 
           &LDB, 
           ncm_vector_data (d), 
           ncm_vector_data (p),
           ncm_vector_data (y),
           &work,
           &lwork,
           &info);

  if (info != 0)
  {
    g_error ("ncm_lapack_dggglm_alloc: cannot estimate size for dggglm.");
  }

  {
    GArray *a = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), work);
    g_array_set_size (a, work);
    return a;
  }
#else
  g_error ("ncm_lapack_dggglm_alloc: lapack support is necessary.");
  return NULL;
#endif
}

/**
 * ncm_lapack_dggglm_run:
 * @ws: (in) (array) (element-type double): a workspace
 * @L: a #NcmMatrix
 * @X: a #NcmMatrix
 * @p: a #NcmVector
 * @d: a #NcmVector
 * @y: a #NcmVector
 * 
 * Runs the dggglm function using the workspace @ws.
 * 
 * This function is expect the matrix @X and @L to be row-major.
 * 
 */
gint
ncm_lapack_dggglm_run (GArray *ws, NcmMatrix *L, NcmMatrix *X, NcmVector *p, NcmVector *d, NcmVector *y)
{
#ifdef HAVE_LAPACK
  gint N   = ncm_matrix_nrows (L);
  gint M   = ncm_matrix_ncols (X);
  gint P   = ncm_matrix_ncols (L);
  gint LDA = N;
  gint LDB = N;
  gdouble *work = &g_array_index (ws, gdouble, 0);
  gint lwork    = ws->len;
  gint info     = 0;
  
  g_assert_cmpint (N, ==, ncm_matrix_nrows (X));
  g_assert_cmpint (N, ==, ncm_vector_len (d));
  g_assert_cmpint (M, ==, ncm_vector_len (p));
  g_assert_cmpint (P, ==, ncm_vector_len (y));

  dggglm_ (&N, &M, &P, 
           ncm_matrix_data (X), 
           &LDA,
           ncm_matrix_data (L), 
           &LDB, 
           ncm_vector_data (d), 
           ncm_vector_data (p),
           ncm_vector_data (y),
           work,
           &lwork,
           &info);
  return info;
#else
  g_error ("ncm_lapack_dggglm_alloc: lapack support is necessary.");
  return -1;
#endif
}

void dtrsv_ (char *UL, char *T, char *D, gint *N, double *A, gint *LDA,
             double *X, gint *INCX);

/**
 * ncm_lapack_dtrsv:
 * @uplo: FIXME
 * @trans: FIXME
 * @diag: FIXME
 * @A: FIXME
 * @v: FIXME
 * 
 * Runs the dtrsv function.
 * 
 */
void
ncm_lapack_dtrsv (gchar uplo, gchar trans, gchar diag, NcmMatrix *A, NcmVector *v)
{
  gint N      = ncm_vector_len (v);
  gint stride = ncm_vector_stride (v);
  trans       = _NCM_LAPACK_CONV_TRANS (trans);
  uplo        = _NCM_LAPACK_CONV_UPLO (uplo);
  
  dtrsv_ (&uplo, &trans, &diag, &N, ncm_matrix_data (A), &N, ncm_vector_data (v), &stride);
}
