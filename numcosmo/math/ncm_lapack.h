/***************************************************************************
 *            ncm_lapack.h
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

#ifndef _NCM_LAPACK_H_
#define _NCM_LAPACK_H_

#include <glib.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/build_cfg.h>

#define _NCM_LAPACK_CONV_UPLO(uplo) (uplo == 'L' ? 'U' : 'L')
#define _NCM_LAPACK_CONV_TRANS(trans) (trans == 'N' ? 'T' : 'N')

G_BEGIN_DECLS

gint ncm_lapack_dptsv (gdouble *d, gdouble *e, gdouble *b, gdouble *x, guint size);
gint ncm_lapack_dpotrf (gchar uplo, guint size, gdouble *a, guint lda);
gint ncm_lapack_dpotri (gchar uplo, guint size, gdouble *a, guint lda);
gint ncm_lapack_dpotrs (gchar uplo, guint size, guint nrhs, gdouble *a, guint lda, gdouble *b, guint ldb);
gint ncm_lapack_dposv (gchar uplo, guint size, guint nrhs, gdouble *a, guint lda, gdouble *b, guint ldb);

GArray *ncm_lapack_dggglm_alloc (NcmMatrix *L, NcmMatrix *X, NcmVector *p, NcmVector *d, NcmVector *y);
gint ncm_lapack_dggglm_run (GArray *ws, NcmMatrix *L, NcmMatrix *X, NcmVector *p, NcmVector *d, NcmVector *y);
void ncm_lapack_dtrsv (gchar uplo, gchar trans, gchar diag, NcmMatrix *A, NcmVector *v);

#define NCM_LAPACK_CHECK_INFO(func,info) if ((info) != 0) g_error ("Lapack[%s] error %d", func, (info))

G_END_DECLS

#endif /* _NCM_LAPACK_H_ */

#ifndef _NCM_LAPACK_INLINE_H_
#define _NCM_LAPACK_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_LAPACK_INLINE_H_ */
