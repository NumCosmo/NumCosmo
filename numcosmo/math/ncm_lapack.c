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

/*#define NUMCOSMO_PREFER_LAPACKE 1*/

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
#include "math/ncm_flapack.h"
#endif /* HAVE_LAPACK */
#endif /* NUMCOSMO_GIR_SCAN */

#define _NCM_LAPACK_CONV_UPLO(uplo) (uplo == 'L' ? 'U' : 'L')
#define _NCM_LAPACK_CONV_TRANS(trans) (trans == 'N' ? 'T' : 'N')

G_DEFINE_BOXED_TYPE (NcmLapackWS, ncm_lapack_ws, ncm_lapack_ws_dup, ncm_lapack_ws_free);

/**
 * ncm_lapack_ws_new:
 * 
 * Creates a new Lapack workspace object.
 * 
 * Returns:(transfer full): a newly created #NcmLapackWS.
 */ 
NcmLapackWS *
ncm_lapack_ws_new (void)
{
  NcmLapackWS *ws = g_new (NcmLapackWS, 1);

  ws->work  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ws->iwork = g_array_new (FALSE, FALSE, sizeof (gint));
  
  return ws;
}

/**
 * ncm_lapack_ws_dup:
 * @ws: a #NcmLapackWS
 * 
 * Duplicates a Lapack workspace object.
 * 
 * Returns: (transfer full): a copy of @ws.
 */ 
NcmLapackWS *
ncm_lapack_ws_dup (NcmLapackWS *ws)
{
  NcmLapackWS *ws_dup = ncm_lapack_ws_new ();

  g_array_set_size (ws_dup->work,  ws->work->len);
  g_array_set_size (ws_dup->iwork, ws->iwork->len);

  return ws_dup;
}

/**
 * ncm_lapack_ws_free:
 * @ws: a #NcmLapackWS
 * 
 * Frees a Lapack workspace object.
 * 
 */ 
void
ncm_lapack_ws_free (NcmLapackWS *ws)
{
  g_array_unref (ws->work);
  g_array_unref (ws->iwork);
  g_free (ws);
}

/**
 * ncm_lapack_ws_clear:
 * @ws: a #NcmLapackWS
 * 
 * Clears a Lapack workspace object.
 * 
 */ 
void
ncm_lapack_ws_clear (NcmLapackWS **ws)
{
  g_assert (ws != NULL);
  if (*ws != NULL)
  {
    g_array_unref (ws[0]->work);
    g_array_unref (ws[0]->iwork);
    g_free (ws[0]);
    ws[0] = NULL;
  }
}

/**
 * ncm_lapack_dptsv:
 * @d: array of doubles with dimension @n
 * @e: array of doubles with dimension @n -1
 * @b: array of doubles with dimension @n
 * @x: array of doubles with dimension @n
 * @n: The order of the matrix $A$ (>= 0)
 *
 * This function computes the solution to a real system of linear equations
 * $A*X = B$ (B = @b), where $A$ is an N-by-N (N = @n) symmetric positive definite tridiagonal
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
ncm_lapack_dptsv (gdouble *d, gdouble *e, gdouble *b, gdouble *x, gint n)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dptsv (LAPACK_ROW_MAJOR, n, 1, d, e, b, 1);
  if (x != b)
    memcpy (x, b, sizeof (gdouble) * n);
  return info;
#elif defined HAVE_LAPACK
	gint NRHS = 1;
	gint LDB  = n;
	gint info;
	
	dptsv_ (&n, &NRHS, d, e, b, &LDB, &info);
  if (x != b)
    memcpy (x, b, sizeof (gdouble) * n);
	
	return info;
#else
  gsl_vector_view b_vec       = gsl_vector_view_array (b, n);
  gsl_vector_view diag_vec    = gsl_vector_view_array (d, n);
  gsl_vector_view offdiag_vec = gsl_vector_view_array (e, n - 1);
  gsl_vector_view x_vec       = gsl_vector_view_array (x, n);
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
 * @n: The order of the matrix @a, @n >= 0
 * @a: array of doubles with dimension (@n, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1,@n)
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
ncm_lapack_dpotrf (gchar uplo, gint n, gdouble *a, gint lda)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dpotrf (LAPACK_ROW_MAJOR, uplo, n, a, lda);
  return info;
#elif defined (HAVE_CLAPACK) && defined (NUMCOSMO_PREFER_LAPACKE)
  gint ret = clapack_dpotrf (CblasRowMajor, 
                             uplo == 'U' ? CblasUpper : CblasLower, 
                             n, a, lda);
  return ret;
#elif defined HAVE_LAPACK
  gint info = 0;
  uplo      = _NCM_LAPACK_CONV_UPLO (uplo);
	
  dpotrf_ (&uplo, &n, a, &lda, &info);  

	return info;
#else /* Fall back to gsl cholesky */
  gint ret;
  gsl_matrix_view mv = gsl_matrix_view_array_with_tda (a, n, n, lda);
  ret = gsl_linalg_cholesky_decomp (&mv.matrix);
  NCM_TEST_GSL_RESULT ("gsl_linalg_cholesky_decomp", ret);
  return ret; /* THAT'S NOT OK FIXME */
#endif
}

/**
 * ncm_lapack_dpotri:
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @n: The order of the matrix @a, @n >= 0
 * @a: array of doubles with dimension (@n, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1,@n)
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
ncm_lapack_dpotri (gchar uplo, gint n, gdouble *a, gint lda)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dpotri (LAPACK_ROW_MAJOR, uplo, n, a, lda);
  return info;
#elif defined (HAVE_CLAPACK) && defined (NUMCOSMO_PREFER_LAPACKE)
  gint ret = clapack_dpotri (CblasRowMajor, 
                             uplo == 'U' ? CblasUpper : CblasLower, 
                             n, a, lda);
  return ret;
#elif defined HAVE_LAPACK
  gint info = 0;
  uplo      = _NCM_LAPACK_CONV_UPLO (uplo);
	
  dpotri_ (&uplo, &n, a, &lda, &info);  
	
  return info;
#else /* Fall back to gsl cholesky */
  gint ret;
  gsl_matrix_view mv = gsl_matrix_view_array_with_tda (a, n, n, lda);
  ret = gsl_linalg_cholesky_invert (&mv.matrix);
  NCM_TEST_GSL_RESULT ("gsl_linalg_cholesky_decomp", ret);
  return ret; /* THAT'S NOT OK FIXME */
#endif
}

/**
 * ncm_lapack_dpotrs:
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @n: The order of the matrix @a, @n >= 0
 * @nrhs: Number of right-hand-side vectors to solve
 * @a: array of doubles with dimension (@n, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1, @n)
 * @b: array of doubles with dimension (@n, @ldb)
 * @ldb: The leading dimension of the array @b, @ldb >= max (1, @n)
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
ncm_lapack_dpotrs (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *b, gint ldb)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dpotrs (LAPACK_ROW_MAJOR, uplo, n, nrhs, a, lda, b, ldb);
  return info;
#elif defined (HAVE_CLAPACK) && defined (NUMCOSMO_PREFER_LAPACKE)
	gint ret = clapack_dpotrs (CblasRowMajor, 
	                          uplo == 'U' ? CblasUpper : CblasLower, 
	                          n, a, lda, b, ldb);
  return ret;
#elif defined HAVE_LAPACK
  gint info = 0;
  uplo      = _NCM_LAPACK_CONV_UPLO (uplo);
	
  dpotrs_ (&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);  

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dpotrs: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dposv:
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @n: The order of the matrix @a, @n >= 0
 * @nrhs: Number of right-hand-side vectors to solve
 * @a: array of doubles with dimension (@n, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1, @n)
 * @b: array of doubles with dimension (@n, @ldb)
 * @ldb: The leading dimension of the array @b, @ldb >= max (1, @n)
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
ncm_lapack_dposv (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *b, gint ldb)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dposv (LAPACK_ROW_MAJOR, uplo, n, nrhs, a, lda, b, ldb);
  return info;
#elif defined (HAVE_CLAPACK) && defined (NUMCOSMO_PREFER_LAPACKE)
	gint ret = clapack_dposv (CblasRowMajor, 
	                          uplo == 'U' ? CblasUpper : CblasLower, 
	                          n, a, lda, b, ldb);
  return ret;
#elif defined HAVE_LAPACK
  gint info = 0;
  uplo      = _NCM_LAPACK_CONV_UPLO (uplo);
	
  dposv_ (&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);  

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dposv: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dsytrf:
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @n: The order of the matrix @a, @n >= 0
 * @a: array of doubles with dimension (@n, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1, @n)
 * @ipiv: Information about decomposition swaps and blocks
 * @ws: a #NcmLapackWS
 *
 * This function computes the factorization of a real symmetric
 * matrix @a, using the Bunch-Kaufman diagonal pivoting method.
 * 
 * Returns: i = 0:  successful exit
 * 
 *          < 0:  -i, the i-th argument had an illegal value
 *
 *          > 0: the (i,i) element of the factor U
 *            or L is zero, and the inverse could not be computed.
 */
gint 
ncm_lapack_dsytrf (gchar uplo, gint n, gdouble *a, gint lda, gint *ipiv, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  gdouble lwork_size;
  lapack_int info;

  info = LAPACKE_dsytrf_work (LAPACK_ROW_MAJOR, uplo, n, a, lda, ipiv, &lwork_size, -1);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);
  info = LAPACKE_dsytrf_work (LAPACK_ROW_MAJOR, uplo, n, a, lda, ipiv, ws->work->data, ws->work->len);  

  return info;
#elif defined HAVE_LAPACK
  gdouble lwork_size;
  gint lwork = -1;
  gint info  = 0;

  uplo      = _NCM_LAPACK_CONV_UPLO (uplo);
	
  dsytrf_ (&uplo, &n, a, &lda, ipiv, &lwork_size, &lwork, &info);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);

  lwork = ws->work->len;
  dsytrf_ (&uplo, &n, a, &lda, ipiv, &g_array_index (ws->work, gdouble, 0), &lwork, &info);

	return info;
#else /* No fall back. */
	g_error ("ncm_lapack_dsytrf: no lapack support!");
	return -1;
#endif
}

/**
 * ncm_lapack_dsytrs:
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @n: The order of the matrix @a, @n >= 0
 * @nrhs: Number of right-hand-side vectors to solve
 * @a: array of doubles with dimension (@n, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1, @n)
 * @ipiv: Information about decomposition swaps and blocks
 * @b: array of doubles with dimension (@n, @ldb)
 * @ldb: The leading dimension of the array @b, @ldb >= max (1, @n)
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
ncm_lapack_dsytrs (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gint *ipiv, gdouble *b, gint ldb)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dsytrs (LAPACK_ROW_MAJOR, uplo, n, nrhs, a, lda, ipiv, b, ldb);
  return info;
#elif defined HAVE_LAPACK
  gint info = 0;
  uplo      = _NCM_LAPACK_CONV_UPLO (uplo);
	
  dsytrs_ (&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);  

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dsytri:
 * @uplo: 'U' upper triangle of @a is stored; 'L' lower triangle of @a is stored
 * @n: The order of the matrix @a, @n >= 0
 * @a: array of doubles with dimension (@n, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1, @n)
 * @ipiv: Information about decomposition swaps and blocks
 * @ws: a #NcmLapackWS
 *
 * This function compute the inverse of a real symmetric indefinite matrix @a using
 * the factorization @a =	U*D*U**T or @a =	L*D*L**T computed by ncm_lapack_dsytrf().
 *
 * Returns: i = 0:  successful exit
 * 
 *          < 0:  -i, the i-th argument had an illegal value
 *
 *          > 0: the (i,i) element of the factor U
 *            or L is zero, and the inverse could not be computed.
 */
gint 
ncm_lapack_dsytri (gchar uplo, gint n, gdouble *a, gint lda, gint *ipiv, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info;

  if (ws->work->len < n)
    g_array_set_size (ws->work, n);
  
  info = LAPACKE_dsytri (LAPACK_ROW_MAJOR, uplo, n, nrhs, a, lda, ipiv, &g_array_index (ws->work, gdouble, 0));
  return info;
#elif defined HAVE_LAPACK
  gint info = 0;
  uplo      = _NCM_LAPACK_CONV_UPLO (uplo);
	
  if (ws->work->len < n)
    g_array_set_size (ws->work, n);

  dsytri_ (&uplo, &n, a, &lda, ipiv, &g_array_index (ws->work, gdouble, 0), &info);  

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dsysvxx:
 * @fact: FACT is CHARACTER*1
 * Specifies whether or not the factored form of the matrix A is
 * supplied on entry, and if not, whether the matrix A should be
 * equilibrated before it is factored.
 * - = 'F':  On entry, AF and IPIV contain the factored form of A.
 *   If EQUED is not 'N', the matrix A has been
 *   equilibrated with scaling factors given by S.
 *   A, AF, and IPIV are not modified.
 * - = 'N':  The matrix A will be copied to AF and factored.
 * - = 'E':  The matrix A will be equilibrated if necessary, then
 *   copied to AF and factored.
 * @uplo: UPLO is CHARACTER*1
 * - = 'U':  Upper triangle of A is stored;
 * - = 'L':  Lower triangle of A is stored.
 * @n: N is INTEGER
 * The number of linear equations, i.e., the order of the
 * matrix A.  N >= 0.
 * @nrhs: NRHS is INTEGER
 * The number of right hand sides, i.e., the number of columns
 * of the matrices B and X.  NRHS >= 0.
 * @a: A is DOUBLE PRECISION array, dimension (LDA,N)
 * The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
 * upper triangular part of A contains the upper triangular
 * part of the matrix A, and the strictly lower triangular
 * part of A is not referenced.  If UPLO = 'L', the leading
 * N-by-N lower triangular part of A contains the lower
 * triangular part of the matrix A, and the strictly upper
 * triangular part of A is not referenced.
 * On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by
 * diag(S)*A*diag(S).
 * @lda: LDA is INTEGER
 * The leading dimension of the array A.  LDA >= max(1,N).
 * @af: AF is DOUBLE PRECISION array, dimension (LDAF,N)
 * - If FACT = 'F', then AF is an input argument and on entry
 *   contains the block diagonal matrix D and the multipliers
 *   used to obtain the factor U or L from the factorization A =
 *   U*D*U**T or A = L*D*L**T as computed by DSYTRF.
 * - If FACT = 'N', then AF is an output argument and on exit
 *   returns the block diagonal matrix D and the multipliers
 *   used to obtain the factor U or L from the factorization A =
 *   U*D*U**T or A = L*D*L**T.
 * @ldaf: LDA is INTEGER
 * The leading dimension of the array A.  LDA >= max(1,N).
 * @ipiv: IPIV is INTEGER array, dimension (N)
 * If FACT = 'F', then IPIV is an input argument and on entry
 * contains details of the interchanges and the block
 * structure of D, as determined by DSYTRF.  If IPIV(k) > 0,
 * then rows and columns k and IPIV(k) were interchanged and
 * D(k,k) is a 1-by-1 diagonal block.  If UPLO = 'U' and
 * IPIV(k) = IPIV(k-1) < 0, then rows and columns k-1 and
 * -IPIV(k) were interchanged and D(k-1:k,k-1:k) is a 2-by-2
 * diagonal block.  If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0,
 * then rows and columns k+1 and -IPIV(k) were interchanged
 * and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
 * If FACT = 'N', then IPIV is an output argument and on exit
 * contains details of the interchanges and the block
 * structure of D, as determined by DSYTRF.
 * @equed: EQUED is CHARACTER*1
 * Specifies the form of equilibration that was done.
 * - = 'N':  No equilibration (always true if FACT = 'N').
 * - = 'Y':  Both row and column equilibration, i.e., A has been
 *   replaced by diag(S) * A * diag(S).
 *   EQUED is an input argument if FACT = 'F'; otherwise, it is an
 *   output argument.
 * @s: S is DOUBLE PRECISION array, dimension (N)
 * The scale factors for A.  If EQUED = 'Y', A is multiplied on
 * the left and right by diag(S).  S is an input argument if FACT =
 * 'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED
 * = 'Y', each element of S must be positive.  If S is output, each
 * element of S is a power of the radix. If S is input, each element
 * of S should be a power of the radix to ensure a reliable solution
 * and error estimates. Scaling by powers of the radix does not cause
 * rounding errors unless the result underflows or overflows.
 * Rounding errors during scaling lead to refining with a matrix that
 * is not equivalent to the input matrix, producing error estimates
 * that may not be reliable.
 * @b: B is DOUBLE PRECISION array, dimension (LDB,NRHS)
 * On entry, the N-by-NRHS right hand side matrix B.
 * On exit,
 * - if EQUED = 'N', B is not modified;
 * - if EQUED = 'Y', B is overwritten by diag(S)*B;
 * @ldb: LDB is INTEGER
 * The leading dimension of the array B.  LDB >= max(1,N).
 * @x: X is DOUBLE PRECISION array, dimension (LDX,NRHS)
 * If INFO = 0, the N-by-NRHS solution matrix X to the original
 * system of equations.  Note that A and B are modified on exit if
 * EQUED .ne. 'N', and the solution to the equilibrated system is
 * inv(diag(S))*X.
 * @ldx: LDX is INTEGER
 * The leading dimension of the array X.  LDX >= max(1,N).
 * @rcond: RCOND is DOUBLE PRECISION
 * Reciprocal scaled condition number.  This is an estimate of the
 * reciprocal Skeel condition number of the matrix A after
 * equilibration (if done).  If this is less than the machine
 * precision (in particular, if it is zero), the matrix is singular
 * to working precision.  Note that the error may still be small even
 * if this number is very small and the matrix appears ill-
 * conditioned.
 * @rpvgrw: RPVGRW is DOUBLE PRECISION
 * Reciprocal pivot growth.  On exit, this contains the reciprocal
 * pivot growth factor norm(A)/norm(U). The "max absolute element"
 * norm is used.  If this is much less than 1, then the stability of
 * the LU factorization of the (equilibrated) matrix A could be poor.
 * This also means that the solution X, estimated condition numbers,
 * and error bounds could be unreliable. If factorization fails with
 * 0<INFO<=N, then this contains the reciprocal pivot growth factor
 * for the leading INFO columns of A.
 * @berr: BERR is DOUBLE PRECISION array, dimension (NRHS)
 * Componentwise relative backward error.  This is the
 * componentwise relative backward error of each solution vector X(j)
 * (i.e., the smallest relative change in any element of A or B that
 * makes X(j) an exact solution).
 * @n_err_bnds: N_ERR_BNDS is INTEGER
 * Number of error bounds to return for each right hand side
 * and each type (normwise or componentwise).  See ERR_BNDS_NORM and
 * ERR_BNDS_COMP below.
 * @err_bnds_norm: ERR_BNDS_NORM is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)
 * For each right-hand side, this array contains information about
 * various error bounds and condition numbers corresponding to the
 * normwise relative error, which is defined as follows:
 * Normwise relative error in the ith solution vector:
 *         max_j (abs(XTRUE(j,i) - X(j,i)))
 *        ------------------------------
 *              max_j abs(X(j,i))
 * The array is indexed by the type of error information as described
 * below. There currently are up to three pieces of information
 * returned.
 * The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
 * right-hand side.
 * The second index in ERR_BNDS_NORM(:,err) contains the following
 * three fields:
 * - err = 1 "Trust/don't trust" boolean. Trust the answer if the
 *   reciprocal condition number is less than the threshold
 *   sqrt(n) * dlamch('Epsilon').
 * - err = 2 "Guaranteed" error bound: The estimated forward error,
 *   almost certainly within a factor of 10 of the true error
 *   so long as the next entry is greater than the threshold
 *   sqrt(n) * dlamch('Epsilon'). This error bound should only
 *   be trusted if the previous boolean is true.
 * - err = 3  Reciprocal condition number: Estimated normwise
 *   reciprocal condition number.  Compared with the threshold
 *   sqrt(n) * dlamch('Epsilon') to determine if the error
 *   estimate is "guaranteed". These reciprocal condition
 *   numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
 *   appropriately scaled matrix Z.
 *   Let Z = S*A, where S scales each row by a power of the
 *   radix so all absolute row sums of Z are approximately 1.
 * See Lapack Working Note 165 for further details and extra
 * cautions.
 * @err_bnds_comp: ERR_BNDS_COMP is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)
 * For each right-hand side, this array contains information about
 * various error bounds and condition numbers corresponding to the
 * componentwise relative error, which is defined as follows:
 * Componentwise relative error in the ith solution vector:
 *                abs(XTRUE(j,i) - X(j,i))
 *         max_j ----------------------
 *                     abs(X(j,i))
 * The array is indexed by the right-hand side i (on which the
 * componentwise relative error depends), and the type of error
 * information as described below. There currently are up to three
 * pieces of information returned for each right-hand side. If
 * componentwise accuracy is not requested (PARAMS(3) = 0.0), then
 * ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most
 * the first (:,N_ERR_BNDS) entries are returned.
 * The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
 * right-hand side.
 * The second index in ERR_BNDS_COMP(:,err) contains the following
 * three fields:
 * - err = 1 "Trust/don't trust" boolean. Trust the answer if the
 *   reciprocal condition number is less than the threshold
 *   sqrt(n) * dlamch('Epsilon').
 * - err = 2 "Guaranteed" error bound: The estimated forward error,
 *   almost certainly within a factor of 10 of the true error
 *   so long as the next entry is greater than the threshold
 *   sqrt(n) * dlamch('Epsilon'). This error bound should only
 *   be trusted if the previous boolean is true.
 * - err = 3  Reciprocal condition number: Estimated componentwise
 *   reciprocal condition number.  Compared with the threshold
 *   sqrt(n) * dlamch('Epsilon') to determine if the error
 *   estimate is "guaranteed". These reciprocal condition
 *   numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
 *   appropriately scaled matrix Z.
 *   Let Z = S*(A*diag(x)), where x is the solution for the
 *   current right-hand side and S scales each row of
 *   A*diag(x) by a power of the radix so all absolute row
 *   sums of Z are approximately 1.
 * See Lapack Working Note 165 for further details and extra
 * cautions.
 * @nparams: NPARAMS is INTEGER
 * Specifies the number of parameters set in PARAMS.  If .LE. 0, the
 * PARAMS array is never referenced and default values are used.
 * @params: PARAMS is DOUBLE PRECISION array, dimension (NPARAMS)
 * Specifies algorithm parameters.  If an entry is .LT. 0.0, then
 * that entry will be filled with default value used for that
 * parameter.  Only positions up to NPARAMS are accessed; defaults
 * are used for higher-numbered parameters.
 * - PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
 *   refinement or not.
 *   Default: 1.0D+0
 *   - = 0.0 : No refinement is performed, and no error bounds are
 *     computed.
 *   - = 1.0 : Use the extra-precise refinement algorithm.
 *     (other values are reserved for future use)
 * - PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
 *   computations allowed for refinement.
 *   Default: 10
 *   Aggressive: Set to 100 to permit convergence using approximate
 *   factorizations or factorizations other than LU. If
 *   the factorization uses a technique other than
 *   Gaussian elimination, the guarantees in
 *   err_bnds_norm and err_bnds_comp may no longer be
 *   trustworthy.
 * - PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
 *   will attempt to find a solution with small componentwise
 *   relative error in the double-precision algorithm.  Positive
 *   is true, 0.0 is false.
 *   Default: 1.0 (attempt componentwise convergence)
 * @work: WORK is DOUBLE PRECISION array, dimension (4*N)
 * @iwork: IWORK is INTEGER array, dimension (N)
 * 
 * # Purpose #
 * 
 * DSYSVXX uses the diagonal pivoting factorization to compute the
 * solution to a double precision system of linear equations A * X = B, where A
 * is an N-by-N symmetric matrix and X and B are N-by-NRHS matrices.
 *
 * If requested, both normwise and maximum componentwise error bounds
 * are returned. DSYSVXX will return a solution with a tiny
 * guaranteed error (O(eps) where eps is the working machine
 * precision) unless the matrix is very ill-conditioned, in which
 * case a warning is returned. Relevant condition numbers also are
 * calculated and returned.
 *
 * DSYSVXX accepts user-provided factorizations and equilibration
 * factors; see the definitions of the FACT and EQUED options.
 * Solving with refinement and using a factorization from a previous
 * DSYSVXX call will also produce a solution with either O(eps)
 * errors or warnings, but we cannot make that claim for general
 * user-provided factorizations and equilibration factors if they
 * differ from what DSYSVXX would itself produce.
 * 
 * # Description #
 * 
 * The following steps are performed:
 * 
 * 1. If FACT = 'E', double precision scaling factors are computed to equilibrate
 *    the system:
 *    - diag(S)*A*diag(S)     *inv(diag(S))*X = diag(S)*B
 *    Whether or not the system will be equilibrated depends on the
 *    scaling of the matrix A, but if equilibration is used, A is
 *    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.
 * 2. If FACT = 'N' or 'E', the LU decomposition is used to factor
 *    the matrix A (after equilibration if FACT = 'E') as
 *    - A = U * D * U**T,  if UPLO = 'U', or
 *    - A = L * D * L**T,  if UPLO = 'L',
 *    where U (or L) is a product of permutation and unit upper (lower)
 *    triangular matrices, and D is symmetric and block diagonal with
 *    1-by-1 and 2-by-2 diagonal blocks.
 * 3. If some D(i,i)=0, so that D is exactly singular, then the
 *    routine returns with INFO = i. Otherwise, the factored form of A
 *    is used to estimate the condition number of the matrix A (see
 *    argument RCOND).  If the reciprocal of the condition number is
 *    less than machine precision, the routine still goes on to solve
 *    for X and compute error bounds as described below.
 * 4. The system of equations is solved for X using the factored form
 *    of A.
 * 5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero),
 *    the routine will use iterative refinement to try to get a small
 *    error and error bounds.  Refinement calculates the residual to at
 *    least twice the working precision.
 * 6. If equilibration was used, the matrix X is premultiplied by
 *    diag(R) so that it solves the original system before
 *    equilibration.
 * 
 * Some optional parameters are bundled in the PARAMS array.  These
 * settings determine how refinement is performed, but often the
 * defaults are acceptable.  If the defaults are acceptable, users
 * can pass NPARAMS = 0 which prevents the source code from accessing
 * the PARAMS argument.
 * 
 * Returns: INFO is INTEGER
 * - = 0:  Successful exit. The solution to every right-hand side is
 *   guaranteed.
 * - < 0:  If INFO = -i, the i-th argument had an illegal value
 * - > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
 *   has been completed, but the factor U is exactly singular, so
 *   the solution and error bounds could not be computed. RCOND = 0
 *   is returned.
 * - = N+J: The solution corresponding to the Jth right-hand side is
 *   not guaranteed. The solutions corresponding to other right-
 *   hand sides K with K > J may not be guaranteed as well, but
 *   only the first such right-hand side is reported. If a small
 *   componentwise error is not requested (PARAMS(3) = 0.0) then
 *   the Jth right-hand side is the first with a normwise error
 *   bound that is not guaranteed (the smallest J such
 *   that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
 *   the Jth right-hand side is the first with either a normwise or
 *   componentwise error bound that is not guaranteed (the smallest
 *   J such that either ERR_BNDS_NORM(J,1) = 0.0 or
 *   ERR_BNDS_COMP(J,1) = 0.0). See the definition of
 *   ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
 *   about all of the right-hand sides check ERR_BNDS_NORM or
 *   ERR_BNDS_COMP.
 */
gint 
ncm_lapack_dsysvxx (gchar fact, gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *af, gint ldaf, gint *ipiv, gchar *equed, gdouble *s, gdouble *b, gint ldb, gdouble *x, gint ldx, gdouble *rcond, gdouble *rpvgrw, gdouble *berr, const gint n_err_bnds, gdouble *err_bnds_norm, gdouble *err_bnds_comp, const gint nparams, gdouble *params, gdouble *work, gint *iwork)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dsysvxx_work (LAPACK_ROW_MAJOR, fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork);
  return info;
#elif defined (HAVE_LAPACK) && defined (HAVE_DSYSVXX_)
  gint info = 0;
  uplo      = _NCM_LAPACK_CONV_UPLO (uplo);
	
	dsysvxx_ (&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, s, b, &ldb, x, &ldx, rcond, rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dsyevr:
 * @jobz: FIXME
 * @range: FIXME
 * @uplo: FIXME
 * @n: FIXME
 * @a: FIXME
 * @lda: FIXME
 * @vl: FIXME
 * @vu: FIXME
 * @il: FIXME
 * @iu: FIXME
 * @abstol: FIXME
 * @m: FIXME
 * @w: FIXME
 * @z: FIXME
 * @ldz: FIXME
 * @isuppz: FIXME
 * @ws: a #NcmLapackWS
 * 
 * FIXME
 * 
 * Returns: FIXME
 */ 
gint 
ncm_lapack_dsyevr (gchar jobz, gchar range, gchar uplo, gint n, gdouble *a, gint lda, gdouble vl, gdouble vu, gint il, gint iu, gdouble abstol, gint *m, gdouble *w, gdouble *z, gint ldz, gint *isuppz, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  gint liwork_size;
  gdouble lwork_size;
  
  info = LAPACKE_dsyevr_work (LAPACK_ROW_MAJOR, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, &lwork_size, -1, &liwork_size, -1);
  if (ws->work->len < lwork_size)
    g_array_set_size (ws->work, lwork_size);
  if (ws->iwork->len < liwork_size)
    g_array_set_size (ws->iwork, liwork_size);
  info = LAPACKE_dsyevr_work (LAPACK_ROW_MAJOR, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, &g_array_index (ws->work, gdouble, 0), ws->work->len, &g_array_index (ws->iwork, gint, 0), ws->iwork->len);

  return info;
#elif defined (HAVE_LAPACK) && defined (HAVE_DSYEVR_)
  gint lwork  = -1;
  gint liwork = -1;
  gint info   = 0;
  gint liwork_size;
  gdouble lwork_size;

  uplo        = _NCM_LAPACK_CONV_UPLO (uplo);
	
	dsyevr_ (&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &lwork_size, &lwork, &liwork_size, &liwork, &info);
  if (ws->work->len < lwork_size)
    g_array_set_size (ws->work, lwork_size);
  if (ws->iwork->len < liwork_size)
    g_array_set_size (ws->iwork, liwork_size);
  lwork  = lwork_size;
  liwork = liwork_size;
  dsyevr_ (&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &g_array_index (ws->work, gdouble, 0), &lwork, &g_array_index (ws->iwork, gint, 0), &liwork, &info);
  
	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsyevr: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dsyevd:
 * @jobz: FIXME
 * @uplo: FIXME
 * @n: FIXME
 * @a: FIXME
 * @lda: FIXME
 * @w: FIXME
 * @ws: a #NcmLapackWS
 * 
 * FIXME
 * 
 * Returns: FIXME
 */ 
gint 
ncm_lapack_dsyevd (gchar jobz, gchar uplo, gint n, gdouble *a, gint lda, gdouble *w, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  gint liwork_size;
  gdouble lwork_size;
  
  info = LAPACKE_dsyevd_work (LAPACK_ROW_MAJOR, jobz, uplo, n, a, lda, w, &lwork_size, -1, &liwork_size, -1);
  if (ws->work->len < lwork_size)
    g_array_set_size (ws->work, lwork_size);
  if (ws->iwork->len < liwork_size)
    g_array_set_size (ws->iwork, liwork_size);
  info = LAPACKE_dsyevr_work (LAPACK_ROW_MAJOR, jobz, uplo, n, a, lda, w, &g_array_index (ws->work, gdouble, 0), ws->work->len, &g_array_index (ws->iwork, gint, 0), ws->iwork->len);

  return info;
#elif defined (HAVE_LAPACK) && defined (HAVE_DSYEVD_)
  gint lwork  = -1;
  gint liwork = -1;
  gint info   = 0;
  gint liwork_size;
  gdouble lwork_size;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);
	
	dsyevd_ (&jobz, &uplo, &n, a, &lda, w, &lwork_size, &lwork, &liwork_size, &liwork, &info);
  if (ws->work->len < lwork_size)
    g_array_set_size (ws->work, lwork_size);
  if (ws->iwork->len < liwork_size)
    g_array_set_size (ws->iwork, liwork_size);
  lwork  = lwork_size;
  liwork = liwork_size;
  dsyevd_ (&jobz, &uplo, &n, a, &lda, w, &g_array_index (ws->work, gdouble, 0), &lwork, &g_array_index (ws->iwork, gint, 0), &liwork, &info);
  
	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsyevr: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dsysvx:
 * @fact: FACT is CHARACTER*1
 * Specifies whether or not the factored form of the matrix A is
 * supplied on entry, and if not, whether the matrix A should be
 * equilibrated before it is factored.
 * - = 'F':  On entry, AF and IPIV contain the factored form of A.
 *   If EQUED is not 'N', the matrix A has been
 *   equilibrated with scaling factors given by S.
 *   A, AF, and IPIV are not modified.
 * - = 'N':  The matrix A will be copied to AF and factored.
 * @uplo: UPLO is CHARACTER*1
 * - = 'U':  Upper triangle of A is stored;
 * - = 'L':  Lower triangle of A is stored.
 * @n: N is INTEGER
 * The number of linear equations, i.e., the order of the
 * matrix A.  N >= 0.
 * @nrhs: NRHS is INTEGER
 * The number of right hand sides, i.e., the number of columns
 * of the matrices B and X.  NRHS >= 0.
 * @a: A is DOUBLE PRECISION array, dimension (LDA,N)
 * The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
 * upper triangular part of A contains the upper triangular
 * part of the matrix A, and the strictly lower triangular
 * part of A is not referenced.  If UPLO = 'L', the leading
 * N-by-N lower triangular part of A contains the lower
 * triangular part of the matrix A, and the strictly upper
 * triangular part of A is not referenced.
 * @lda: LDA is INTEGER
 * The leading dimension of the array A.  LDA >= max(1,N).
 * @af: AF is DOUBLE PRECISION array, dimension (LDAF,N)
 * - If FACT = 'F', then AF is an input argument and on entry
 *   contains the block diagonal matrix D and the multipliers
 *   used to obtain the factor U or L from the factorization A =
 *   U*D*U**T or A = L*D*L**T as computed by DSYTRF.
 * - If FACT = 'N', then AF is an output argument and on exit
 *   returns the block diagonal matrix D and the multipliers
 *   used to obtain the factor U or L from the factorization A =
 *   U*D*U**T or A = L*D*L**T.
 * @ldaf: LDA is INTEGER
 * The leading dimension of the array A.  LDA >= max(1,N).
 * @ipiv: IPIV is INTEGER array, dimension (N)
 * If FACT = 'F', then IPIV is an input argument and on entry
 * contains details of the interchanges and the block
 * structure of D, as determined by DSYTRF.  If IPIV(k) > 0,
 * then rows and columns k and IPIV(k) were interchanged and
 * D(k,k) is a 1-by-1 diagonal block.  If UPLO = 'U' and
 * IPIV(k) = IPIV(k-1) < 0, then rows and columns k-1 and
 * -IPIV(k) were interchanged and D(k-1:k,k-1:k) is a 2-by-2
 * diagonal block.  If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0,
 * then rows and columns k+1 and -IPIV(k) were interchanged
 * and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
 * If FACT = 'N', then IPIV is an output argument and on exit
 * contains details of the interchanges and the block
 * structure of D, as determined by DSYTRF.
 * @b: B is DOUBLE PRECISION array, dimension (LDB,NRHS)
 * On entry, the N-by-NRHS right hand side matrix B.
 * On exit,
 * - if EQUED = 'N', B is not modified;
 * - if EQUED = 'Y', B is overwritten by diag(S)*B;
 * @ldb: LDB is INTEGER
 * The leading dimension of the array B.  LDB >= max(1,N).
 * @x: X is DOUBLE PRECISION array, dimension (LDX,NRHS)
 * If INFO = 0, the N-by-NRHS solution matrix X to the original
 * system of equations.  Note that A and B are modified on exit if
 * EQUED .ne. 'N', and the solution to the equilibrated system is
 * inv(diag(S))*X.
 * @ldx: LDX is INTEGER
 * The leading dimension of the array X.  LDX >= max(1,N).
 * @rcond: RCOND is DOUBLE PRECISION
 * Reciprocal scaled condition number.  This is an estimate of the
 * reciprocal Skeel condition number of the matrix A after
 * equilibration (if done).  If this is less than the machine
 * precision (in particular, if it is zero), the matrix is singular
 * to working precision.  Note that the error may still be small even
 * if this number is very small and the matrix appears ill-
 * conditioned.
 * @ferr: FERR is DOUBLE PRECISION array, dimension (NRHS)
 * The estimated forward error bound for each solution vector
 * X(j) (the j-th column of the solution matrix X).
 * If XTRUE is the true solution corresponding to X(j), FERR(j)
 * is an estimated upper bound for the magnitude of the largest
 * element in (X(j) - XTRUE) divided by the magnitude of the
 * largest element in X(j).  The estimate is as reliable as
 * the estimate for RCOND, and is almost always a slight
 * overestimate of the true error..
 * @berr: BERR is DOUBLE PRECISION array, dimension (NRHS)
 * Componentwise relative backward error.  This is the
 * componentwise relative backward error of each solution vector X(j)
 * (i.e., the smallest relative change in any element of A or B that
 * makes X(j) an exact solution).
 * @work: WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
 * On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 * @lwork: LWORK is INTEGER
 * The length of WORK.  LWORK >= max(1,3*N), and for best
 * performance, when FACT = 'N', LWORK >= max(1,3*N,N*NB), where
 * NB is the optimal blocksize for DSYTRF.
 * If LWORK = -1, then a workspace query is assumed; the routine
 * only calculates the optimal size of the WORK array, returns
 * this value as the first entry of the WORK array, and no error
 * message related to LWORK is issued by XERBLA.
 * @iwork: IWORK is INTEGER array, dimension (N)
 * 
 * # Purpose #
 * 
 * DSYSVX uses the diagonal pivoting factorization to compute the
 * solution to a real system of linear equations A * X = B,
 * where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
 * matrices.
 * 
 * Error bounds on the solution and a condition estimate are also
 * provided.
 * 
 * # Description #
 * 
 * The following steps are performed:
 * 
 * 1. If FACT = 'N', the diagonal pivoting method is used to factor A.
 *    The form of the factorization is
 * 		- A = U * D * U**T,  if UPLO = 'U', or
 *    - A = L * D * L**T,  if UPLO = 'L',
 *    where U (or L) is a product of permutation and unit upper (lower)
 *    triangular matrices, and D is symmetric and block diagonal with
 *    1-by-1 and 2-by-2 diagonal blocks.
 * 2. If some D(i,i)=0, so that D is exactly singular, then the routine
 *    returns with INFO = i. Otherwise, the factored form of A is used
 *    to estimate the condition number of the matrix A.  If the
 *    reciprocal of the condition number is less than machine precision,
 *    INFO = N+1 is returned as a warning, but the routine still goes on
 *    to solve for X and compute error bounds as described below.
 * 3. The system of equations is solved for X using the factored form
 *    of A.
 * 4. Iterative refinement is applied to improve the computed solution
 *    matrix and calculate error bounds and backward error estimates
 *    for it.
 * 
 * Returns: INFO is INTEGER
 * - = 0: successful exit
 * - < 0: if INFO = -i, the i-th argument had an illegal value
 * - > 0: if INFO = i, and i is
 * - <= N:  D(i,i) is exactly zero.  The factorization
 *   has been completed but the factor D is exactly
 *   singular, so the solution and error bounds could
 *   not be computed. RCOND = 0 is returned.
 * - = N+1: D is nonsingular, but RCOND is less than machine
 *   precision, meaning that the matrix is singular
 *   to working precision.  Nevertheless, the
 *   solution and error bounds are computed because
 *   there are a number of situations where the
 *   computed solution can be more accurate than the
 *   value of RCOND would suggest.
 */
gint 
ncm_lapack_dsysvx (gchar fact, gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *af, gint ldaf, gint *ipiv, gdouble *b, gint ldb, gdouble *x, gint ldx, gdouble *rcond, gdouble *ferr, gdouble *berr, gdouble *work, gint lwork, gint *iwork)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dsysvx_work (LAPACK_ROW_MAJOR, fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, lwork, iwork);
  return info;
#elif defined (HAVE_LAPACK) && defined (HAVE_DSYSVX_)
  gint info = 0;
  uplo      = _NCM_LAPACK_CONV_UPLO (uplo);
	
	dsysvx_ (&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, &lwork, iwork, &info);

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dgeev:
 * @jobvl: @n left eigenvectors of @a are not computed, 'V' left eigenvectors of @a are computed
 * @jobvr: @n right eigenvectors of @a are not computed, 'V' right eigenvectors of @a are computed
 * @n: The order of the matrix @a, @n >= 0
 * @a: array of doubles with dimension (@n, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1, @n)
 * @wr: contain the real part of the computed eigenvalues
 * @wi: contain the imaginary part of the computed eigenvalues
 * @vl: if @jobvl = 'V', the left eigenvectors $u(j)$ are stored one after another in the rows of @vl, in the same order as their eigenvalues
 * @ldvl: the leading dimension of the array @vl
 * @vr: if @jobvr = 'V', the left eigenvectors $v(j)$ are stored one after another in the rows of @vr, in the same order as their eigenvalues
 * @ldvr: the leading dimension of the array @vr
 * @work: work area, must have @lwork allocated doubles
 * @lwork: work area size
 *
 * This function computes the eigensystem for a real matrix @a = A.
 * 
 * Calling this function with lwork == -1 computed the ideal @lwork in @work[0].
 *
 * Returns: i = 0:  successful exit
 * 
 *          < 0:  -i, the i-th argument had an illegal value
 *
 *          > 0: the (i,i) element of the factor U
 *            or L is zero, and the inverse could not be computed.
 */
gint 
ncm_lapack_dgeev (gchar jobvl, gchar jobvr, gint n, gdouble *a, gint lda, gdouble *wr, gdouble *wi, gdouble *vl, gint ldvl, gdouble *vr, gint ldvr, gdouble *work, gint lwork)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dgeev_work (LAPACK_ROW_MAJOR, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork);
  return info;
#elif defined HAVE_LAPACK
  gint info = 0;
	
	/* swap L <=> R : col-major <=> row-major */
  dgeev_ (&jobvr, &jobvl, &n, a, &lda, wr, wi, vr, &ldvr, vl, &ldvl, work, &lwork, &info);  

	return info;
#else /* No fall back. */
	g_error ("ncm_lapack_dgeev: no lapack support!");
	return -1;
#endif
}

/**
 * ncm_lapack_dgeevx:
 * @balanc: FIXME
 * @jobvl: @n left eigenvectors of @a are not computed, 'V' left eigenvectors of @a are computed
 * @jobvr: @n right eigenvectors of @a are not computed, 'V' right eigenvectors of @a are computed
 * @sense: FIXME
 * @n: The order of the matrix @a, @n >= 0
 * @a: array of doubles with dimension (@n, @lda)
 * @lda: The leading dimension of the array @a, @lda >= max (1, @n)
 * @wr: contain the real part of the computed eigenvalues
 * @wi: contain the imaginary part of the computed eigenvalues
 * @vl: if @jobvl = 'V', the left eigenvectors $u(j)$ are stored one after another in the rows of @vl, in the same order as their eigenvalues
 * @ldvl: the leading dimension of the array @vl
 * @vr: if @jobvr = 'V', the left eigenvectors $v(j)$ are stored one after another in the rows of @vr, in the same order as their eigenvalues
 * @ldvr: the leading dimension of the array @vr
 * @ilo: FIXME
 * @ihi: FIXME
 * @scale: FIXME
 * @abnrm: FIXME
 * @rconde: FIXME
 * @rcondv: FIXME
 * @work: work area, must have @lwork allocated doubles
 * @lwork: work area size
 * @iwork: FIXME
 *
 * This function computes the eigensystem for a real matrix @a = A.
 * 
 * Calling this function with lwork == -1 computed the ideal @lwork in @work[0].
 *
 * Returns: i = 0:  successful exit
 * 
 *          < 0:  -i, the i-th argument had an illegal value
 *
 *          > 0: the (i,i) element of the factor U
 *            or L is zero, and the inverse could not be computed.
 */
gint 
ncm_lapack_dgeevx (gchar balanc, gchar jobvl, gchar jobvr, gchar sense, gint n, gdouble *a, gint lda, gdouble *wr, gdouble *wi, gdouble *vl, gint ldvl, gdouble *vr, gint ldvr, gint *ilo, gint *ihi, gdouble *scale, gdouble *abnrm, gdouble *rconde, gdouble *rcondv, gdouble *work, gint lwork, gint *iwork)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  lapack_int info = LAPACKE_dgeevx_work (LAPACK_ROW_MAJOR, balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork);
  return info;
#elif defined HAVE_LAPACK
  gint info = 0;
	
	/* swap L <=> R : col-major <=> row-major */
	dgeevx_ (&balanc, &jobvr, &jobvl, &sense, &n, a, &lda, wr, wi, vr, &ldvr, vl, &ldvl, ilo, ihi, scale, abnrm, rconde, rcondv, work, &lwork, iwork, &info);

	return info;
#else /* No fall back. */
	g_error ("ncm_lapack_dgeev: no lapack support!");
	return -1;
#endif
}

/**
 * ncm_lapack_dgeqrf:
 * @m: M is INTEGER
 * The number of rows of the matrix A.  M >= 0.
 * @n: N is INTEGER
 * The number of columns of the matrix A.  N >= 0.
 * @a: A is DOUBLE PRECISION array, dimension (LDA,N)
 * On entry, the M-by-N matrix A.
 * On exit, the elements on and above the diagonal of the array
 * contain the min(M,N)-by-N upper trapezoidal matrix R (R is
 * upper triangular if m >= n); the elements below the diagonal,
 * with the array TAU, represent the orthogonal matrix Q as a
 * product of min(m,n) elementary reflectors (see Further
 * Details).
 * @lda: LDA is INTEGER
 * The leading dimension of the array A.  LDA >= max(1,M).
 * @tau: TAU is DOUBLE PRECISION array, dimension (min(M,N))
 * The scalar factors of the elementary reflectors (see Further
 * Details).
 * @ws: a #NcmLapackWS
 * 
 * DGEQRF computes a QR factorization of a real M-by-N matrix A:
 * A = Q * R.
 * 
 * Returns: = 0:  successful exit
 * < 0:  if INFO = -i, the i-th argument had an illegal value
 */ 
gint 
ncm_lapack_dgeqrf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  gdouble lwork_size;
  lapack_int info;

  info = LAPACKE_dgeqrf_work (LAPACK_ROW_MAJOR, m, n, a, lda, tau, &lwork_size, -1);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);
  info = LAPACKE_dgeqrf_work (LAPACK_ROW_MAJOR, m, n, a, lda, tau, &g_array_index (ws->work, gdouble, 0), ws->work->len);

  return info;
#elif defined (HAVE_LAPACK) && defined (HAVE_DGELQF_) /* To account for row-major => col-major QR => LQ */
  gdouble lwork_size;
  gint lwork = -1;
  gint info  = 0;
  
	dgelqf_ (&m, &n, a, &lda, tau, &lwork_size, &lwork, &info);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);
  lwork = ws->work->len;
	dgelqf_ (&m, &n, a, &lda, tau, &g_array_index (ws->work, gdouble, 0), &lwork, &info);

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dgerqf:
 * @m: M is INTEGER
 * The number of rows of the matrix A.  M >= 0.
 * @n: N is INTEGER
 * The number of columns of the matrix A.  N >= 0.
 * @a: A is DOUBLE PRECISION array, dimension (LDA,N)
 * On entry, the M-by-N matrix A.
 * On exit, the elements on and above the diagonal of the array
 * contain the min(M,N)-by-N upper trapezoidal matrix R (R is
 * upper triangular if m >= n); the elements below the diagonal,
 * with the array TAU, represent the orthogonal matrix Q as a
 * product of min(m,n) elementary reflectors (see Further
 * Details).
 * @lda: LDA is INTEGER
 * The leading dimension of the array A.  LDA >= max(1,M).
 * @tau: TAU is DOUBLE PRECISION array, dimension (min(M,N))
 * The scalar factors of the elementary reflectors (see Further
 * Details).
 * @ws: a #NcmLapackWS
 * 
 * DGERQF computes a RQ factorization of a real M-by-N matrix A:
 * A = R * Q.
 * 
 * Returns: = 0:  successful exit
 * < 0:  if INFO = -i, the i-th argument had an illegal value
 */ 
gint 
ncm_lapack_dgerqf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  gdouble lwork_size;
  lapack_int info;

  info = LAPACKE_dgerqf_work (LAPACK_ROW_MAJOR, m, n, a, lda, tau, &lwork_size, -1);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);
  info = LAPACKE_dgerqf_work (LAPACK_ROW_MAJOR, m, n, a, lda, tau, &g_array_index (ws->work, gdouble, 0), ws->work->len);

  return info;
#elif defined (HAVE_LAPACK) && defined (HAVE_DGEQLF_) /* To account for row-major => col-major RQ => QL */
  gdouble lwork_size;
  gint lwork = -1;
  gint info  = 0;
  
	dgeqlf_ (&m, &n, a, &lda, tau, &lwork_size, &lwork, &info);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);
  lwork = ws->work->len;
	dgeqlf_ (&m, &n, a, &lda, tau, &g_array_index (ws->work, gdouble, 0), &lwork, &info);

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dgeqlf:
 * @m: M is INTEGER
 * The number of rows of the matrix A.  M >= 0.
 * @n: N is INTEGER
 * The number of columns of the matrix A.  N >= 0.
 * @a: A is DOUBLE PRECISION array, dimension (LDA,N)
 * On entry, the M-by-N matrix A.
 * On exit, the elements on and above the diagonal of the array
 * contain the min(M,N)-by-N upper trapezoidal matrix R (R is
 * upper triangular if m >= n); the elements below the diagonal,
 * with the array TAU, represent the orthogonal matrix Q as a
 * product of min(m,n) elementary reflectors (see Further
 * Details).
 * @lda: LDA is INTEGER
 * The leading dimension of the array A.  LDA >= max(1,M).
 * @tau: TAU is DOUBLE PRECISION array, dimension (min(M,N))
 * The scalar factors of the elementary reflectors (see Further
 * Details).
 * @ws: a #NcmLapackWS
 * 
 * DGEQLF computes a QL factorization of a real M-by-N matrix A:
 * A = Q * L.
 * 
 * Returns: = 0:  successful exit
 * < 0:  if INFO = -i, the i-th argument had an illegal value
 */ 
gint 
ncm_lapack_dgeqlf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  gdouble lwork_size;
  lapack_int info;

  info = LAPACKE_dgeqlf_work (LAPACK_ROW_MAJOR, m, n, a, lda, tau, &lwork_size, -1);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);
  info = LAPACKE_dgeqlf_work (LAPACK_ROW_MAJOR, m, n, a, lda, tau, &g_array_index (ws->work, gdouble, 0), ws->work->len);

  return info;
#elif defined (HAVE_LAPACK) && defined (HAVE_DGERQF_) /* To account for row-major => col-major QL => RQ */
  gdouble lwork_size;
  gint lwork = -1;
  gint info  = 0;
  
	dgerqf_ (&m, &n, a, &lda, tau, &lwork_size, &lwork, &info);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);
  lwork = ws->work->len;
	dgerqf_ (&m, &n, a, &lda, tau, &g_array_index (ws->work, gdouble, 0), &lwork, &info);

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");
#endif
}

/**
 * ncm_lapack_dgelqf:
 * @m: M is INTEGER
 * The number of rows of the matrix A.  M >= 0.
 * @n: N is INTEGER
 * The number of columns of the matrix A.  N >= 0.
 * @a: A is DOUBLE PRECISION array, dimension (LDA,N)
 * On entry, the M-by-N matrix A.
 * On exit, the elements on and above the diagonal of the array
 * contain the min(M,N)-by-N upper trapezoidal matrix R (R is
 * upper triangular if m >= n); the elements below the diagonal,
 * with the array TAU, represent the orthogonal matrix Q as a
 * product of min(m,n) elementary reflectors (see Further
 * Details).
 * @lda: LDA is INTEGER
 * The leading dimension of the array A.  LDA >= max(1,M).
 * @tau: TAU is DOUBLE PRECISION array, dimension (min(M,N))
 * The scalar factors of the elementary reflectors (see Further
 * Details).
 * @ws: a #NcmLapackWS
 * 
 * DGELQF computes a LQ factorization of a real M-by-N matrix A:
 * A = L * Q.
 * 
 * Returns: = 0:  successful exit
 * < 0:  if INFO = -i, the i-th argument had an illegal value
 */ 
gint 
ncm_lapack_dgelqf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACKE) && defined (NUMCOSMO_PREFER_LAPACKE)
  gdouble lwork_size;
  lapack_int info;

  info = LAPACKE_dgelqf_work (LAPACK_ROW_MAJOR, m, n, a, lda, tau, &lwork_size, -1);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);
  info = LAPACKE_dgelqf_work (LAPACK_ROW_MAJOR, m, n, a, lda, tau, &g_array_index (ws->work, gdouble, 0), ws->work->len);

  return info;
#elif defined (HAVE_LAPACK) && defined (HAVE_DGEQRF_) /* To account for row-major => col-major LQ => QR */
  gdouble lwork_size;
  gint lwork = -1;
  gint info  = 0;
  
	dgeqrf_ (&m, &n, a, &lda, tau, &lwork_size, &lwork, &info);
  if (lwork_size > ws->work->len)
    g_array_set_size (ws->work, lwork_size);
  lwork = ws->work->len;
	dgeqrf_ (&m, &n, a, &lda, tau, &g_array_index (ws->work, gdouble, 0), &lwork, &info);

	return info;
#else /* No fall back */
	g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");
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

void dtrsv_ (char *UL, char *T, char *D, gint *N, gdouble *A, gint *LDA,
             gdouble *X, gint *INCX);

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
