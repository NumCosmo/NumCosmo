/***************************************************************************
 *            ncm_lapack.h
 *
 *  Sun March 18 22:33:15 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_LAPACK_H_
#define _NCM_LAPACK_H_

#include <glib.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/build_cfg.h>

#define _NCM_LAPACK_CONV_UPLO(uplo) (uplo == 'L' ? 'U' : 'L')
#define _NCM_LAPACK_CONV_TRANS(trans) (trans == 'N' ? 'T' : 'N')

G_BEGIN_DECLS

typedef struct _NcmLapackWS NcmLapackWS;

struct _NcmLapackWS
{
  GArray *work;
  GArray *iwork;
};

GType ncm_lapack_ws_get_type (void) G_GNUC_CONST;

NcmLapackWS *ncm_lapack_ws_new (void);
NcmLapackWS *ncm_lapack_ws_dup (NcmLapackWS *ws);
void ncm_lapack_ws_free (NcmLapackWS *ws);
void ncm_lapack_ws_clear (NcmLapackWS **ws);

gint ncm_lapack_dptsv (gdouble *d, gdouble *e, gdouble *b, gdouble *x, gint n);
gint ncm_lapack_dpotrf (gchar uplo, gint n, gdouble *a, gint lda);
gint ncm_lapack_dpotri (gchar uplo, gint n, gdouble *a, gint lda);
gint ncm_lapack_dpotrs (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *b, gint ldb);
gint ncm_lapack_dposv (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *b, gint ldb);

gint ncm_lapack_dsytrf (gchar uplo, gint n, gdouble *a, gint lda, gint *ipiv, NcmLapackWS *ws);
gint ncm_lapack_dsytrs (gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gint *ipiv, gdouble *b, gint ldb);
gint ncm_lapack_dsytri (gchar uplo, gint n, gdouble *a, gint lda, gint *ipiv, NcmLapackWS *ws);
gint ncm_lapack_dsysvx (gchar fact, gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *af, gint ldaf, gint *ipiv, gdouble *b, gint ldb, gdouble *x, gint ldx, gdouble *rcond, gdouble *ferr, gdouble *berr, gdouble *work, gint lwork, gint *iwork);
gint ncm_lapack_dsysvxx (gchar fact, gchar uplo, gint n, gint nrhs, gdouble *a, gint lda, gdouble *af, gint ldaf, gint *ipiv, gchar *equed, gdouble *s, gdouble *b, gint ldb, gdouble *x, gint ldx, gdouble *rcond, gdouble *rpvgrw, gdouble *berr, const gint n_err_bnds, gdouble *err_bnds_norm, gdouble *err_bnds_comp, const gint nparams, gdouble *params, gdouble *work, gint *iwork);

gint ncm_lapack_dsyevr (gchar jobz, gchar range, gchar uplo, gint n, gdouble *a, gint lda, gdouble vl, gdouble vu, gint il, gint iu, gdouble abstol, gint *m, gdouble *w, gdouble *z, gint ldz, gint *isuppz, NcmLapackWS *ws);
gint ncm_lapack_dsyevd (gchar jobz, gchar uplo, gint n, gdouble *a, gint lda, gdouble *w, NcmLapackWS *ws);

gint ncm_lapack_dgeev (gchar jobvl, gchar jobvr, gint n, gdouble *a, gint lda, gdouble *wr, gdouble *wi, gdouble *vl, gint ldvl, gdouble *vr, gint ldvr, gdouble *work, gint lwork);
gint ncm_lapack_dgeevx (gchar balanc, gchar jobvl, gchar jobvr, gchar sense, gint n, gdouble *a, gint lda, gdouble *wr, gdouble *wi, gdouble *vl, gint ldvl, gdouble *vr, gint ldvr, gint *ilo, gint *ihi, gdouble *scale, gdouble *abnrm, gdouble *rconde, gdouble *rcondv, gdouble *work, gint lwork, gint *iwork);

gint ncm_lapack_dgeqrf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws);
gint ncm_lapack_dgerqf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws);
gint ncm_lapack_dgeqlf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws);
gint ncm_lapack_dgelqf (gint m, gint n, gdouble *a, gint lda, gdouble *tau, NcmLapackWS *ws);

GArray *ncm_lapack_dggglm_alloc (NcmMatrix *L, NcmMatrix *X, NcmVector *p, NcmVector *d, NcmVector *y);
gint ncm_lapack_dggglm_run (GArray *ws, NcmMatrix *L, NcmMatrix *X, NcmVector *p, NcmVector *d, NcmVector *y);

#define NCM_LAPACK_CHECK_INFO(func,info) G_STMT_START { if ((info) != 0) g_error ("# NcmLapack[%s] error %4d", func, (info)); } G_STMT_END 

G_END_DECLS

#endif /* _NCM_LAPACK_H_ */

#ifndef _NCM_LAPACK_INLINE_H_
#define _NCM_LAPACK_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_LAPACK_INLINE_H_ */
