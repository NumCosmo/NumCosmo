/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_nnls.h
 *
 *  Thu April 8 10:06:36 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_nnls.h
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_NNLS_H_
#define _NCM_NNLS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_NNLS             (ncm_nnls_get_type ())
#define NCM_NNLS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_NNLS, NcmNNLS))
#define NCM_NNLS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_NNLS, NcmNNLSClass))
#define NCM_IS_NNLS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_NNLS))
#define NCM_IS_NNLS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_NNLS))
#define NCM_NNLS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_NNLS, NcmNNLSClass))

typedef struct _NcmNNLSClass NcmNNLSClass;
typedef struct _NcmNNLS NcmNNLS;
typedef struct _NcmNNLSPrivate NcmNNLSPrivate;

struct _NcmNNLSClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmNNLS
{
  /*< private >*/
  GObject parent_instance;
  NcmNNLSPrivate *priv;
};

/**
 * NcmNNLSUMethod:
 * @NCM_NNLS_UMETHOD_NORMAL: Solve using normal equations and Cholesky decomposition
 * @NCM_NNLS_UMETHOD_NORMAL_LU: Solve using normal equations and LU decomposition
 * @NCM_NNLS_UMETHOD_QR: Solve using QR decomposition
 * @NCM_NNLS_UMETHOD_DGELSD: Solve using QR decomposition (Lapack's dgelsd)
 * @NCM_NNLS_UMETHOD_GSL: Solve using GSL's gsl_multifit_linear
 *
 * Method used to solve the intermediate unconstrained least-squares.
 */
typedef enum _NcmNNLSUMethod /*< enum,underscore_name=NCM_NNLS_UMETHOD >*/
{
  NCM_NNLS_UMETHOD_NORMAL,
  NCM_NNLS_UMETHOD_NORMAL_LU,
  NCM_NNLS_UMETHOD_QR,
  NCM_NNLS_UMETHOD_DGELSD,
  NCM_NNLS_UMETHOD_GSL,
  /* < private > */
  NCM_NNLS_UMETHOD_LEN, /*< skip >*/
} NcmNNLSUMethod;

GType ncm_nnls_get_type (void) G_GNUC_CONST;

NcmNNLS *ncm_nnls_new (guint nrows, guint ncols);
NcmNNLS *ncm_nnls_ref (NcmNNLS *nnls);

void ncm_nnls_free (NcmNNLS *nnls);
void ncm_nnls_clear (NcmNNLS **nnls);

void ncm_nnls_set_umethod (NcmNNLS *nnls, NcmNNLSUMethod umethod);
NcmNNLSUMethod ncm_nnls_get_umethod (NcmNNLS *nnls);

void ncm_nnls_set_reltol (NcmNNLS *nnls, const gdouble reltol);
gdouble ncm_nnls_get_reltol (NcmNNLS *nnls);

guint ncm_nnls_get_nrows (NcmNNLS *nnls);
guint ncm_nnls_get_ncols (NcmNNLS *nnls);

gdouble ncm_nnls_solve (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f);
gdouble ncm_nnls_solve_LH (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f);
gdouble ncm_nnls_solve_lowrankqp (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f);
gdouble ncm_nnls_solve_splx (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f);
gdouble ncm_nnls_solve_gsmo (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f);

NcmVector *ncm_nnls_get_residuals (NcmNNLS *nnls);

G_END_DECLS

#endif /* _NCM_NNLS_H_ */
