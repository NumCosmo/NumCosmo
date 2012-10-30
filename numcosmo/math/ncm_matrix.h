/***************************************************************************
 *            ncm_matrix.h
 *
 *  Thu January 05 20:18:45 2012
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

#ifndef _NCM_MATRIX_H_
#define _NCM_MATRIX_H_

#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_matrix.h>
#include <numcosmo/nc_macros.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_MATRIX             (ncm_matrix_get_type ())
#define NCM_MATRIX(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MATRIX, NcmMatrix))
#define NCM_MATRIX_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MATRIX, NcmMatrixClass))
#define NCM_IS_MATRIX(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MATRIX))
#define NCM_IS_MATRIX_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MATRIX))
#define NCM_MATRIX_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MATRIX, NcmMatrixClass))

typedef struct _NcmMatrixClass NcmMatrixClass;
typedef struct _NcmMatrix NcmMatrix;

/**
 * NcmMatrixInternal:
 * @NCM_MATRIX_SLICE: FIXME
 * @NCM_MATRIX_GSL_MATRIX: FIXME
 * @NCM_MATRIX_MALLOC: FIXME
 * @NCM_MATRIX_GARRAY: FIXME
 * @NCM_MATRIX_DERIVED: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcmMatrixInternal
{
  NCM_MATRIX_SLICE = 0,
  NCM_MATRIX_GSL_MATRIX,
  NCM_MATRIX_MALLOC,
  NCM_MATRIX_GARRAY,
  NCM_MATRIX_DERIVED,
} NcmMatrixInternal;

struct _NcmMatrixClass
{
  /*< private >*/
  GInitiallyUnownedClass parent_class;
};

/**
 * NcmMatrix:
 *
 * FIXME
 */
struct _NcmMatrix
{
  /*< private >*/
  GInitiallyUnowned parent_instance;
  gsl_matrix_view mv;
  gsl_matrix *gm;
  GArray *a;
  GObject *pobj;
  NcmMatrixInternal type;
};

GType ncm_matrix_get_type (void) G_GNUC_CONST;

#define NCM_MATRIX_GSL(cm) (&(cm)->mv.matrix)
#define NCM_MATRIX_COL_LEN(cm) ((cm)->mv.matrix.size1)
#define NCM_MATRIX_ROW_LEN(cm) ((cm)->mv.matrix.size2)
#define NCM_MATRIX_DATA(cm) ((cm)->mv.matrix.data)

#define NCM_MATRIX_NROWS(cm) ((cm)->mv.matrix.size1)
#define NCM_MATRIX_NCOLS(cm) ((cm)->mv.matrix.size2)

NcmMatrix *ncm_matrix_new (const gsize nrows, const gsize ncols);
NcmMatrix *ncm_matrix_new_gsl (gsl_matrix *gm);
NcmMatrix *ncm_matrix_new_array (GArray *a, const gsize ncols);
NcmMatrix *ncm_matrix_new_data_slice (gdouble *d, const gsize nrows, const gsize ncols);
NcmMatrix *ncm_matrix_new_data_malloc (gdouble *d, const gsize nrows, const gsize ncols);
NcmMatrix *ncm_matrix_new_data_static (gdouble *d, const gsize nrows, const gsize ncols);
NcmMatrix *ncm_matrix_new_data_static_tda (gdouble *d, const gsize nrows, const gsize ncols, const gsize tda);

NcmMatrix *ncm_matrix_get_submatrix (NcmMatrix *cm, const gsize k1, const gsize k2, const gsize nrows, const gsize ncols);
NcmVector *ncm_matrix_get_col (NcmMatrix *cm, const gsize col);
NcmVector *ncm_matrix_get_row (NcmMatrix *cm, const gsize row);

G_INLINE_FUNC const NcmMatrix *ncm_matrix_new_gsl_const (gsl_matrix *m);
G_INLINE_FUNC gdouble ncm_matrix_get (const NcmMatrix *cm, const guint i, const guint j);
G_INLINE_FUNC gdouble *ncm_matrix_ptr (NcmMatrix *cm, const guint i, const guint j);
G_INLINE_FUNC NcmMatrix *ncm_matrix_ref (NcmMatrix *cm);
G_INLINE_FUNC GArray *ncm_matrix_get_array (NcmMatrix *cm);
G_INLINE_FUNC void ncm_matrix_set (NcmMatrix *cm, const guint i, const guint j, const gdouble val);
G_INLINE_FUNC void ncm_matrix_transpose (NcmMatrix *cm);
G_INLINE_FUNC void ncm_matrix_set_identity (NcmMatrix *cm);
G_INLINE_FUNC void ncm_matrix_set_zero (NcmMatrix *cm);
G_INLINE_FUNC void ncm_matrix_scale (NcmMatrix *cm, const gdouble val);
G_INLINE_FUNC void ncm_matrix_memcpy (NcmMatrix *cm1, const NcmMatrix *cm2);
G_INLINE_FUNC void ncm_matrix_set_col (NcmMatrix *cm, const guint n, const NcmVector *cv);

NcmMatrix *ncm_matrix_copy (const NcmMatrix *cm);

void ncm_matrix_free (NcmMatrix *cm);

G_END_DECLS

#endif /* _NCM_MATRIX_H_ */

#ifndef _NCM_MATRIX_INLINE_H_
#define _NCM_MATRIX_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC const NcmMatrix *
ncm_matrix_new_gsl_const (gsl_matrix *m)
{
  return ncm_matrix_new_data_static_tda ((m)->data, (m)->size1, (m)->size2, (m)->tda);
}

G_INLINE_FUNC gdouble
ncm_matrix_get (const NcmMatrix *cm, const guint i, const guint j)
{
  return gsl_matrix_get (NCM_MATRIX_GSL (cm), i, j);
}

G_INLINE_FUNC gdouble *
ncm_matrix_ptr (NcmMatrix *cm, guint i, guint j)
{
  return gsl_matrix_ptr (NCM_MATRIX_GSL (cm), i, j);
}

G_INLINE_FUNC NcmMatrix *
ncm_matrix_ref (NcmMatrix *cm)
{
  return g_object_ref_sink (cm);
}

G_INLINE_FUNC void
ncm_matrix_set (NcmMatrix *cm, guint i, guint j, gdouble val)
{
  gsl_matrix_set (NCM_MATRIX_GSL (cm), i, j, val);
}

G_INLINE_FUNC void
ncm_matrix_transpose (NcmMatrix *cm)
{
  const gint ret = gsl_matrix_transpose (NCM_MATRIX_GSL (cm));
  NC_TEST_GSL_RESULT ("gsl_matrix_transpose", ret);
}

G_INLINE_FUNC void
ncm_matrix_set_identity (NcmMatrix *cm)
{
  gsl_matrix_set_identity (NCM_MATRIX_GSL (cm));
}

G_INLINE_FUNC void
ncm_matrix_set_zero (NcmMatrix *cm)
{
  gsl_matrix_set_zero (NCM_MATRIX_GSL (cm));
}

G_INLINE_FUNC void
ncm_matrix_scale (NcmMatrix *cm, gdouble val)
{
  gsl_matrix_scale (NCM_MATRIX_GSL (cm),val);
}

G_INLINE_FUNC void
ncm_matrix_memcpy (NcmMatrix *cm1, const NcmMatrix *cm2)
{
  gsl_matrix_memcpy (NCM_MATRIX_GSL(cm1), NCM_MATRIX_GSL(cm2));
}

G_INLINE_FUNC void
ncm_matrix_set_col (NcmMatrix *cm, const guint n, const NcmVector *cv)
{
  gsl_matrix_set_col (NCM_MATRIX_GSL (cm), n, ncm_vector_const_gsl (cv));
}

G_INLINE_FUNC GArray *
ncm_matrix_get_array (NcmMatrix *cm)
{
  g_assert (cm->a != NULL);
  return g_array_ref (cm->a);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_MATRIX_INLINE_H_ */
