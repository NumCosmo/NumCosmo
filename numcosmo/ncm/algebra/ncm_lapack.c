/***************************************************************************
 *            ncm_lapack.c
 *
 *  Sun March 18 22:33:15 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * NcmLapack:
 *
 * Wrappers of LAPACK routines for the row-major data of #NcmMatrix.
 *
 * LAPACK stores matrices in column-major order, so it reads a row-major array as the
 * transpose. Each wrapper states how it accounts for that: for symmetric matrices the
 * triangle @uplo is converted to the other one, some routines are replaced by their
 * transposed counterpart (QR by LQ, for instance), and the rest read the array in
 * column-major order as passed. The arguments otherwise follow the
 * [LAPACK documentation](https://www.netlib.org/lapack/explore-html/). Routines that need a
 * workspace take a #NcmLapackWS and size it with a workspace query. Without LAPACK, only
 * dptsv, dpotrf and dpotri fall back to GSL; the others abort.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_lapack.h"
#include "ncm/algebra/ncm_matrix.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#ifdef HAVE_LAPACK
#include "ncm/algebra/ncm_flapack.h"
#endif /* HAVE_LAPACK */
#endif /* NUMCOSMO_GIR_SCAN */

#define _NCM_LAPACK_CONV_UPLO(uplo) (uplo == 'L' ? 'U' : 'L')
#define _NCM_LAPACK_CONV_TRANS(trans) (trans == 'N' ? 'T' : 'N')

G_DEFINE_BOXED_TYPE (NcmLapackWS, ncm_lapack_ws, ncm_lapack_ws_dup, ncm_lapack_ws_free)

/**
 * ncm_lapack_ws_new:
 *
 * Creates an empty workspace, grown as needed by the wrappers that take one.
 *
 * Returns: (transfer full): a new #NcmLapackWS.
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
 * Frees @ws.
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
 * If *@ws is not %NULL, frees it and sets *@ws to %NULL.
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
 * @d: the diagonal
 * @e: the off-diagonal, of length @n - 1
 * @b: the right-hand side
 * @x: the solution
 * @n: order of the matrix
 *
 * Solves $A x = b$ for the symmetric positive definite tridiagonal $A$ with diagonal @d and
 * off-diagonal @e, with LAPACK DPTSV, which factors $A = L D L^\intercal$ and overwrites @d, @e
 * and @b. Without LAPACK, uses gsl_linalg_solve_symm_tridiag().
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dptsv (gdouble *d, gdouble *e, gdouble *b, gdouble *x, gint n)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DPTSV_)
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
  gint status                 = gsl_linalg_solve_symm_tridiag (&diag_vec.vector,
                                                               &offdiag_vec.vector,
                                                               &b_vec.vector,
                                                               &x_vec.vector);

  NCM_TEST_GSL_RESULT ("ncm_lapack_dptsv[gsl_linalg_solve_symm_tridiag]", status);

  return status;

#endif /* HAVE_LAPACK */
}

/**
 * ncm_lapack_dpotrf:
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @a: the matrix
 * @lda: its leading dimension
 *
 * Replaces the @uplo triangle of the symmetric positive definite @a by its Cholesky factor,
 * with LAPACK DPOTRF. Without LAPACK, uses gsl_linalg_cholesky_decomp(), which aborts
 * instead of returning an error.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dpotrf (gchar uplo, gint n, gdouble *a, gint lda)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DPOTRF_)
  gint info = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dpotrf_ (&uplo, &n, a, &lda, &info);

  return info;

#else /* Fall back to gsl cholesky */
  gint ret;
  gsl_matrix_view mv = gsl_matrix_view_array_with_tda (a, n, n, lda);

  ret = gsl_linalg_cholesky_decomp (&mv.matrix);
  NCM_TEST_GSL_RESULT ("gsl_linalg_cholesky_decomp", ret);

  return ret;

#endif
}

/**
 * ncm_lapack_dpotri:
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @a: the Cholesky factor
 * @lda: its leading dimension
 *
 * Replaces the Cholesky factor from ncm_lapack_dpotrf() by the @uplo triangle of the
 * inverse, with LAPACK DPOTRI. Without LAPACK, uses gsl_linalg_cholesky_invert(), which
 * aborts instead of returning an error.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dpotri (gchar uplo, gint n, gdouble *a, gint lda)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DPOTRI_)
  gint info = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dpotri_ (&uplo, &n, a, &lda, &info);

  return info;

#else /* Fall back to gsl cholesky */
  gint ret;
  gsl_matrix_view mv = gsl_matrix_view_array_with_tda (a, n, n, lda);

  ret = gsl_linalg_cholesky_invert (&mv.matrix);
  NCM_TEST_GSL_RESULT ("gsl_linalg_cholesky_decomp", ret);

  return ret;

#endif
}

/**
 * ncm_lapack_dpotrs:
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @nrhs: number of right-hand sides
 * @a: the Cholesky factor
 * @lda: its leading dimension
 * @b: the right-hand sides, one per row
 * @ldb: their leading dimension
 *
 * Solves $A x = b$ for each right-hand side with the Cholesky factor from
 * ncm_lapack_dpotrf(), with LAPACK DPOTRS; @b is overwritten by the solutions. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dpotrs (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *b, gint ldb)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DPOTRS_)
  gint info = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dpotrs_ (&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dpotrs: lapack not present, no fallback implemented.");

#endif
}

/**
 * ncm_lapack_dposv:
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @nrhs: number of right-hand sides
 * @a: the matrix
 * @lda: its leading dimension
 * @b: the right-hand sides, one per row
 * @ldb: their leading dimension
 *
 * Solves $A x = b$ for the symmetric positive definite @a, with LAPACK DPOSV: @a is
 * overwritten by its Cholesky factor and @b by the solutions. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dposv (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *b, gint ldb)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DPOSV_)
  gint info = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dposv_ (&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dposv: lapack not present, no fallback implemented.");

#endif
}

/**
 * ncm_lapack_dgesv:
 * @n: order of the matrix
 * @nrhs: number of right-hand sides
 * @a: the matrix, in column-major order
 * @lda: its leading dimension
 * @ipiv: the pivot indices
 * @b: the right-hand sides, one per row
 * @ldb: their leading dimension
 *
 * Solves $A x = b$ for a general @a by LU factorization with partial pivoting, with LAPACK
 * DGESV. @a is read in column-major order, so a row-major #NcmMatrix gives the solution for
 * its transpose. @a is overwritten by the factorization and @b by the solutions. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dgesv (gint n, gint nrhs, gdouble *a, gint lda, gint *ipiv, gdouble *b, gint ldb)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGESV_)
  gint info = 0;

  dgesv_ (&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dgesv: lapack not present, no fallback implemented.");

#endif
}

/**
 * ncm_lapack_dsytrf:
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @a: the matrix
 * @lda: its leading dimension
 * @ipiv: the pivoting information
 * @ws: a #NcmLapackWS
 *
 * Replaces the @uplo triangle of the symmetric @a by its Bunch-Kaufman factorization, with
 * LAPACK DSYTRF. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dsytrf (gchar uplo, gint n, gdouble *a, gint lda, gint *ipiv, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DSYTRF_)
  gdouble lwork_size;
  gint lwork = -1;
  gint info  = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

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
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @nrhs: number of right-hand sides
 * @a: the factorization
 * @lda: its leading dimension
 * @ipiv: the pivoting information
 * @b: the right-hand sides, one per row
 * @ldb: their leading dimension
 *
 * Solves $A x = b$ with the factorization from ncm_lapack_dsytrf(), with LAPACK DSYTRS; @b
 * is overwritten by the solutions. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dsytrs (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gint *ipiv, gdouble *b, gint ldb)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DSYTRS_)
  gint info = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dsytrs_ (&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");

#endif
}

/**
 * ncm_lapack_dsytri:
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @a: the factorization
 * @lda: its leading dimension
 * @ipiv: the pivoting information
 * @ws: a #NcmLapackWS
 *
 * Replaces the factorization from ncm_lapack_dsytrf() by the @uplo triangle of the inverse,
 * with LAPACK DSYTRI. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dsytri (gchar uplo, gint n, gdouble *a, gint lda, gint *ipiv, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DSYTRI_)
  gint info = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  g_assert_cmpint (n, >=, 0);

  if (ws->work->len < (guint) n)
    g_array_set_size (ws->work, n);

  dsytri_ (&uplo, &n, a, &lda, ipiv, &g_array_index (ws->work, gdouble, 0), &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");

#endif
}

/**
 * ncm_lapack_dsysvxx:
 * @fact: 'F', 'N' or 'E', see LAPACK DSYSVXX
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @nrhs: number of right-hand sides
 * @a: the matrix
 * @lda: its leading dimension
 * @af: the factorization
 * @ldaf: its leading dimension
 * @ipiv: the pivoting information
 * @equed: the equilibration done
 * @s: the scale factors
 * @b: the right-hand sides, one per row
 * @ldb: their leading dimension
 * @x: the solutions
 * @ldx: their leading dimension
 * @rcond: the reciprocal condition number
 * @rpvgrw: the reciprocal pivot growth factor
 * @berr: the backward errors
 * @n_err_bnds: number of error bounds
 * @err_bnds_norm: the normwise error bounds
 * @err_bnds_comp: the componentwise error bounds
 * @nparams: number of parameters
 * @params: the algorithm parameters
 * @work: the workspace
 * @iwork: the integer workspace
 *
 * Solves $A x = b$ for the symmetric @a with iterative refinement and error bounds, with
 * LAPACK DSYSVXX. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dsysvxx (gchar fact, gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *af, gint ldaf, gint *ipiv, gchar *equed, gdouble *s, gdouble *b, gint ldb, gdouble *x, gint ldx, gdouble *rcond, gdouble *rpvgrw, gdouble *berr, const gint n_err_bnds, gdouble *err_bnds_norm, gdouble *err_bnds_comp, const gint nparams, gdouble *params, gdouble *work, gint *iwork)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DSYSVXX_)
  gint info = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dsysvxx_ (&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, s, b, &ldb, x, &ldx, rcond, rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dsytrs: lapack not present, no fallback implemented.");

#endif
}

/**
 * ncm_lapack_dsyevr:
 * @jobz: 'N' for eigenvalues only, 'V' for eigenvectors too
 * @range: 'A', 'V' or 'I': all eigenvalues, those in (@vl, @vu], or those of index @il to @iu
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @a: the matrix
 * @lda: its leading dimension
 * @vl: lower bound of the eigenvalues for @range 'V'
 * @vu: upper bound of the eigenvalues for @range 'V'
 * @il: index of the smallest eigenvalue for @range 'I'
 * @iu: index of the largest eigenvalue for @range 'I'
 * @abstol: absolute tolerance of the eigenvalues
 * @m: number of eigenvalues found
 * @w: the eigenvalues, in increasing order
 * @z: the eigenvectors, one per row
 * @ldz: their leading dimension
 * @isuppz: the support of the eigenvectors
 * @ws: a #NcmLapackWS
 *
 * Computes selected eigenvalues and, optionally, eigenvectors of the symmetric @a by the
 * relatively robust representations algorithm, with LAPACK DSYEVR; @a is overwritten. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dsyevr (gchar jobz, gchar range, gchar uplo, gint n, gdouble *a, gint lda, gdouble vl, gdouble vu, gint il, gint iu, gdouble abstol, gint *m, gdouble *w, gdouble *z, gint ldz, gint *isuppz, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DSYEVR_)
  gint lwork  = -1;
  gint liwork = -1;
  gint info   = 0;
  gint liwork_size;
  gdouble lwork_size;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dsyevr_ (&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &lwork_size, &lwork, &liwork_size, &liwork, &info);

  if (ws->work->len < lwork_size)
    g_array_set_size (ws->work, lwork_size);

  g_assert_cmpint (liwork_size, >=, 0);

  if (ws->iwork->len < (guint) liwork_size)
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
 * @jobz: 'N' for eigenvalues only, 'V' for eigenvectors too
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @a: the matrix
 * @lda: its leading dimension
 * @w: the eigenvalues, in increasing order
 * @ws: a #NcmLapackWS
 *
 * Computes all eigenvalues and, optionally, eigenvectors of the symmetric @a by divide and
 * conquer, with LAPACK DSYEVD; with @jobz 'V' the eigenvectors overwrite @a, one per row. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dsyevd (gchar jobz, gchar uplo, gint n, gdouble *a, gint lda, gdouble *w, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DSYEVD_)
  gint lwork  = -1;
  gint liwork = -1;
  gint info   = 0;
  gint liwork_size;
  gdouble lwork_size;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dsyevd_ (&jobz, &uplo, &n, a, &lda, w, &lwork_size, &lwork, &liwork_size, &liwork, &info);

  if (ws->work->len < lwork_size)
    g_array_set_size (ws->work, lwork_size);

  g_assert_cmpint (liwork_size, >=, 0);

  if (ws->iwork->len < (guint) liwork_size)
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
 * ncm_lapack_dsysv:
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @nrhs: number of right-hand sides
 * @a: the matrix
 * @lda: its leading dimension
 * @ipiv: the pivoting information
 * @b: the right-hand sides, one per row
 * @ldb: their leading dimension
 * @work: the workspace
 * @lwork: its length, or -1 for a workspace query
 *
 * Solves $A x = b$ for the symmetric @a, with LAPACK DSYSV. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dsysv (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gint *ipiv, gdouble *b, gint ldb, gdouble *work, gint lwork)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DSYSV_)
  gint info = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dsysv_ (&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dsysv: lapack not present, no fallback implemented.");

#endif
}

/**
 * ncm_lapack_dsysvx:
 * @fact: 'F' or 'N', whether @af and @ipiv hold the factorization
 * @uplo: 'U' or 'L', the triangle of @a stored, in the row-major sense
 * @n: order of the matrix
 * @nrhs: number of right-hand sides
 * @a: the matrix
 * @lda: its leading dimension
 * @af: the factorization
 * @ldaf: its leading dimension
 * @ipiv: the pivoting information
 * @b: the right-hand sides, one per row
 * @ldb: their leading dimension
 * @x: the solutions
 * @ldx: their leading dimension
 * @rcond: the reciprocal condition number
 * @ferr: the forward error bounds
 * @berr: the backward errors
 * @work: the workspace
 * @lwork: its length, or -1 for a workspace query
 * @iwork: the integer workspace
 *
 * Solves $A x = b$ for the symmetric @a, with error bounds and a condition estimate, with
 * LAPACK DSYSVX. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dsysvx (gchar fact, gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *af, gint ldaf, gint *ipiv, gdouble *b, gint ldb, gdouble *x, gint ldx, gdouble *rcond, gdouble *ferr, gdouble *berr, gdouble *work, gint lwork, gint *iwork)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DSYSVX_)
  gint info = 0;

  uplo = _NCM_LAPACK_CONV_UPLO (uplo);

  dsysvx_ (&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, &lwork, iwork, &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dsysvx: lapack not present, no fallback implemented.");

#endif
}

/**
 * ncm_lapack_dgeev:
 * @jobvl: 'N' or 'V', whether to compute the left eigenvectors
 * @jobvr: 'N' or 'V', whether to compute the right eigenvectors
 * @n: order of the matrix
 * @a: the matrix
 * @lda: its leading dimension
 * @wr: real parts of the eigenvalues
 * @wi: imaginary parts of the eigenvalues
 * @vl: the left eigenvectors, one per row
 * @ldvl: their leading dimension
 * @vr: the right eigenvectors, one per row
 * @ldvr: their leading dimension
 * @work: the workspace
 * @lwork: its length, or -1 for a workspace query
 *
 * Computes the eigenvalues and, optionally, the eigenvectors of the general @a, with LAPACK
 * DGEEV. LAPACK reads the transpose, whose left and right eigenvectors are the right and
 * left ones of @a, so the two are exchanged in the call. @a is overwritten. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dgeev (gchar jobvl, gchar jobvr, gint n, gdouble *a, gint lda, gdouble *wr, gdouble *wi, gdouble *vl, gint ldvl, gdouble *vr, gint ldvr, gdouble *work, gint lwork)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGEEV_)
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
 * @balanc: 'N', 'P', 'S' or 'B', the balancing
 * @jobvl: 'N' or 'V', whether to compute the left eigenvectors
 * @jobvr: 'N' or 'V', whether to compute the right eigenvectors
 * @sense: 'N', 'E', 'V' or 'B', the condition numbers computed
 * @n: order of the matrix
 * @a: the matrix
 * @lda: its leading dimension
 * @wr: real parts of the eigenvalues
 * @wi: imaginary parts of the eigenvalues
 * @vl: the left eigenvectors, one per row
 * @ldvl: their leading dimension
 * @vr: the right eigenvectors, one per row
 * @ldvr: their leading dimension
 * @ilo: first index of the balanced block
 * @ihi: last index of the balanced block
 * @scale: the balancing details
 * @abnrm: the one-norm of the balanced matrix
 * @rconde: reciprocal condition numbers of the eigenvalues
 * @rcondv: reciprocal condition numbers of the eigenvectors
 * @work: the workspace
 * @lwork: its length, or -1 for a workspace query
 * @iwork: the integer workspace
 *
 * Same as ncm_lapack_dgeev() with balancing and condition numbers, with LAPACK DGEEVX. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dgeevx (gchar balanc, gchar jobvl, gchar jobvr, gchar sense, gint n, gdouble *a, gint lda, gdouble *wr, gdouble *wi, gdouble *vl, gint ldvl, gdouble *vr, gint ldvr, gint *ilo, gint *ihi, gdouble *scale, gdouble *abnrm, gdouble *rconde, gdouble *rcondv, gdouble *work, gint lwork, gint *iwork)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGEEVX_)
  gint info = 0;

  /* swap L <=> R : col-major <=> row-major */
  dgeevx_ (&balanc, &jobvr, &jobvl, &sense, &n, a, &lda, wr, wi, vr, &ldvr, vl, &ldvl, ilo, ihi, scale, abnrm, rconde, rcondv, work, &lwork, iwork, &info);

  return info;

#else /* No fall back. */
  g_error ("ncm_lapack_dgeevx: no lapack support!");

  return -1;

#endif
}

/**
 * ncm_lapack_dgeqrf:
 * @m: number of rows seen by LAPACK, the columns of the row-major @a
 * @n: number of columns seen by LAPACK, the rows of the row-major @a
 * @a: the matrix
 * @lda: its leading dimension
 * @tau: the scalar factors of the elementary reflectors
 * @ws: a #NcmLapackWS
 *
 * Computes the QR factorization $A = Q R$ of the row-major @a by calling LAPACK DGELQF on
 * the transpose that LAPACK reads; @a is overwritten by the factors, as described by LAPACK
 * DGELQF. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dgeqrf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGELQF_) /* To account for row-major => col-major QR => LQ */
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
 * @m: number of rows seen by LAPACK, the columns of the row-major @a
 * @n: number of columns seen by LAPACK, the rows of the row-major @a
 * @a: the matrix
 * @lda: its leading dimension
 * @tau: the scalar factors of the elementary reflectors
 * @ws: a #NcmLapackWS
 *
 * Computes the RQ factorization $A = R Q$ of the row-major @a by calling LAPACK DGEQLF on
 * the transpose that LAPACK reads; @a is overwritten by the factors, as described by LAPACK
 * DGEQLF. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dgerqf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGEQLF_) /* To account for row-major => col-major RQ => QL */
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
 * @m: number of rows seen by LAPACK, the columns of the row-major @a
 * @n: number of columns seen by LAPACK, the rows of the row-major @a
 * @a: the matrix
 * @lda: its leading dimension
 * @tau: the scalar factors of the elementary reflectors
 * @ws: a #NcmLapackWS
 *
 * Computes the QL factorization $A = Q L$ of the row-major @a by calling LAPACK DGERQF on
 * the transpose that LAPACK reads; @a is overwritten by the factors, as described by LAPACK
 * DGERQF. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dgeqlf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGERQF_) /* To account for row-major => col-major QL => RQ */
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
 * @m: number of rows seen by LAPACK, the columns of the row-major @a
 * @n: number of columns seen by LAPACK, the rows of the row-major @a
 * @a: the matrix
 * @lda: its leading dimension
 * @tau: the scalar factors of the elementary reflectors
 * @ws: a #NcmLapackWS
 *
 * Computes the LQ factorization $A = L Q$ of the row-major @a by calling LAPACK DGEQRF on
 * the transpose that LAPACK reads; @a is overwritten by the factors, as described by LAPACK
 * DGEQRF. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dgelqf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGEQRF_) /* To account for row-major => col-major LQ => QR */
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
 * @L: a row-major #NcmMatrix
 * @X: a row-major #NcmMatrix
 * @p: a #NcmVector
 * @d: a #NcmVector
 * @y: a #NcmVector
 *
 * Allocates the workspace of ncm_lapack_dggglm_run() for these arguments, with a LAPACK
 * DGGGLM workspace query. Aborts without LAPACK.
 *
 * Returns: (transfer full) (array) (element-type double): the workspace.
 */
GArray *
ncm_lapack_dggglm_alloc (NcmMatrix *L, NcmMatrix *X, NcmVector *p, NcmVector *d, NcmVector *y)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGEQRF_)
  gint N   = ncm_matrix_nrows (L);
  gint M   = ncm_matrix_ncols (X);
  gint P   = ncm_matrix_ncols (L);
  gint LDA = N;
  gint LDB = N;
  gdouble work;
  gint lwork = -1;
  gint info  = 0;

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
    g_error ("ncm_lapack_dggglm_alloc: cannot estimate size for dggglm.");

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
 * @ws: (in) (array) (element-type double): the workspace from ncm_lapack_dggglm_alloc()
 * @L: a row-major #NcmMatrix
 * @X: a row-major #NcmMatrix
 * @p: a #NcmVector
 * @d: a #NcmVector
 * @y: a #NcmVector
 *
 * Solves the general Gauss-Markov linear model problem of LAPACK DGGGLM for the given
 * matrices and vectors. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dggglm_run (GArray *ws, NcmMatrix *L, NcmMatrix *X, NcmVector *p, NcmVector *d, NcmVector *y)
{
#ifdef HAVE_LAPACK
  gint N        = ncm_matrix_nrows (L);
  gint M        = ncm_matrix_ncols (X);
  gint P        = ncm_matrix_ncols (L);
  gint LDA      = N;
  gint LDB      = N;
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

/**
 * ncm_lapack_dgels:
 * @trans: 'N' or 'T', whether to use the transpose of @a, in the row-major sense
 * @m: number of rows of the row-major @a
 * @n: number of columns of the row-major @a
 * @nrhs: number of right-hand sides
 * @a: the matrix
 * @lda: its leading dimension
 * @b: the right-hand sides
 * @ldb: their leading dimension
 * @work: the workspace
 * @lwork: its length, or -1 for a workspace query
 *
 * Solves the least-squares or minimum-norm problem for the full-rank @a, with LAPACK
 * DGELS; @trans is converted and @m and @n are exchanged for the column-major call. @a is
 * overwritten by its QR or LQ factorization and @b by the solutions. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dgels (gchar trans, const gint m, const gint n, const gint nrhs, gdouble *a, const gint lda, gdouble *b, const gint ldb, double *work, const gint lwork)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGELS_)
  gint info = 0;

  trans = _NCM_LAPACK_CONV_TRANS (trans);

  dgels_ (&trans, &n, &m, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dgels: lapack not present, no fallback implemented.");

#endif
}

/**
 * ncm_lapack_dgelsd:
 * @m: number of rows seen by LAPACK
 * @n: number of columns seen by LAPACK
 * @nrhs: number of right-hand sides
 * @a: the matrix
 * @lda: its leading dimension
 * @b: the right-hand sides
 * @ldb: their leading dimension
 * @s: the singular values, in decreasing order
 * @rcond: the threshold below which singular values are treated as zero
 * @rank: the effective rank
 * @ws: a #NcmLapackWS
 *
 * Solves the least-squares problem by the singular value decomposition, with LAPACK DGELSD,
 * passing the arguments unchanged, so @a is read in column-major order. Aborts without LAPACK.
 *
 * Returns: the LAPACK `info`: zero on success, $-i$ if argument $i$ is invalid, positive for the
 * failure described in the LAPACK documentation.
 */
gint
ncm_lapack_dgelsd (const gint m, const gint n, const gint nrhs, gdouble *a, const gint lda, gdouble *b, const gint ldb, gdouble *s, gdouble *rcond, gint *rank, NcmLapackWS *ws)
{
#if defined (HAVE_LAPACK) && defined (HAVE_DGELS_)
  gint lwork = -1;
  gint info  = 0;
  gint liwork_size;
  gdouble lwork_size;

  dgelsd_ (&m, &n, &nrhs, a, &lda, b, &ldb, s, rcond, rank, &lwork_size, &lwork, &liwork_size, &info);

  if (ws->work->len < lwork_size)
    g_array_set_size (ws->work, lwork_size);

  g_assert_cmpint (liwork_size, >=, 0);

  if (ws->iwork->len < (guint) liwork_size)
    g_array_set_size (ws->iwork, liwork_size);

  lwork = lwork_size;
  dgelsd_ (&m, &n, &nrhs, a, &lda, b, &ldb, s, rcond, rank,
           &g_array_index (ws->work, gdouble, 0),
           &lwork,
           &g_array_index (ws->iwork, gint, 0),
           &info);

  return info;

#else /* No fall back */
  g_error ("ncm_lapack_dgelsd: lapack not present, no fallback implemented.");
#endif
}

