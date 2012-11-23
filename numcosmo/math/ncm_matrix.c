/***************************************************************************
 *            ncm_matrix.c
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

/**
 * SECTION:ncm_matrix
 * @title: Matrix
 * @short_description: Allocation and access.
 *
 * This object defines the functions for allocating and accessing matrices.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_matrix.h"

G_DEFINE_TYPE (NcmMatrix, ncm_matrix, G_TYPE_INITIALLY_UNOWNED);

/**
 * ncm_matrix_new:
 * @nrows: number of rows.
 * @ncols: number of columns.
 *
 * This function allocates memory for a new #NcmMatrix of doubles
 * with @nrows rows and @ncols columns.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new (gsize nrows, gsize ncols)
{
  NcmMatrix *cm = g_object_new (NCM_TYPE_MATRIX, NULL);
  gdouble *d = g_slice_alloc (sizeof(gdouble) * nrows * ncols);
  cm->mv = gsl_matrix_view_array (d, nrows, ncols);
  cm->type = NCM_MATRIX_SLICE;

  return cm;
}

/**
 * ncm_matrix_new_sunk:
 * @nrows: number of rows.
 * @ncols: number of columns.
 *
 * This function allocates memory for a new #NcmMatrix of doubles
 * with @nrows rows and @ncols columns. Returns a sunk reference.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_sunk (gsize nrows, gsize ncols)
{
  return ncm_matrix_ref (ncm_matrix_new (nrows, ncols));
}

/**
 * ncm_matrix_new_gsl: (skip)
 * @gm: matrix from GNU Scientific Library (GSL) to be converted into a #NcmMatrix.
 *
 * This function saves @gm internally and frees it when it is no longer necessary.
 * The @gm matrix must not be freed.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_gsl (gsl_matrix *gm)
{
  NcmMatrix *cm = g_object_new (NCM_TYPE_MATRIX, NULL);
  cm->gm = gm;
  cm->mv = gsl_matrix_submatrix (gm, 0, 0, gm->size1, gm->size2);
  cm->type = NCM_MATRIX_GSL_MATRIX;
  return cm;
}

/**
 * ncm_matrix_new_gsl_static: (skip)
 * @gm: matrix from GNU Scientific Library (GSL) to be converted into a #NcmMatrix.
 *
 * This function saves @gm internally and does not frees it.
 * The @gm matrix must be valid during the life of the created #NcmMatrix.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_gsl_static (gsl_matrix *gm)
{
  NcmMatrix *cm = ncm_matrix_new_gsl (gm);
  cm->type = NCM_MATRIX_DERIVED;
  return cm;
}

/**
 * ncm_matrix_new_array:
 * @a: GArray of doubles to be converted into a #NcmMatrix.
 * @ncols: number of columns.
 *
 * The number of rows is defined dividing the lenght of @a by @ncols.
 * This function saves @a internally and frees it when it is no longer necessary.
 * The GArray @a must not be freed.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_array (GArray *a, gsize ncols)
{
  NcmMatrix *cm = g_object_new (NCM_TYPE_MATRIX, NULL);
  gsize nrows = a->len / ncols;
  gdouble *d = &g_array_index (a, gdouble, 0);

  g_assert (a->len % ncols == 0);

  cm->mv = gsl_matrix_view_array (d, nrows, ncols);
  cm->a = a;
  cm->type = NCM_MATRIX_GARRAY;
  g_array_ref (a);

  return cm;
}

/**
 * ncm_matrix_new_data_slice: (skip)
 * @d: pointer to the first double allocated.
 * @nrows: number of rows.
 * @ncols: number of columns.
 *
 * This function returns a #NcmMatrix of the array @d allocated using g_slice function.
 * It saves @d internally and frees it when it is no longer necessary.
 * The matrix has @nrows rows and @ncols columns.
 * The physical number of columns in memory is also given by @ncols.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_data_slice (gdouble *d, gsize nrows, gsize ncols)
{
  NcmMatrix *cm = g_object_new (NCM_TYPE_MATRIX, NULL);
  cm->mv = gsl_matrix_view_array (d, nrows, ncols);
  cm->a = NULL;
  cm->type = NCM_MATRIX_SLICE;

  return cm;
}

/**
 * ncm_matrix_new_data_malloc: (skip)
 * @d: pointer to the first double allocated.
 * @nrows: number of rows.
 * @ncols: number of columns.
 *
 * This function returns a #NcmMatrix of the array @d allocated using malloc.
 * It saves @d internally and frees it when it is no longer necessary.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_data_malloc (gdouble *d, gsize nrows, gsize ncols)
{
  NcmMatrix *cm = ncm_matrix_new_data_slice (d, nrows, ncols);
  cm->type = NCM_MATRIX_MALLOC;
  return cm;
}

/**
 * ncm_matrix_new_data_static: (skip)
 * @d: pointer to the first double allocated.
 * @nrows: number of rows.
 * @ncols: number of columns.
 *
 * This function returns a #NcmMatrix of the array @d.
 * The memory allocated is kept during all time life of the object and
 * must not be freed during this period.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_data_static (gdouble *d, gsize nrows, gsize ncols)
{
  NcmMatrix *cm = ncm_matrix_new_data_slice (d, nrows, ncols);
  cm->type = NCM_MATRIX_DERIVED;
  return cm;
}

/**
 * ncm_matrix_new_data_static_tda: (skip)
 * @d: pointer to the first double allocated.
 * @nrows: number of rows.
 * @ncols: number of columns.
 * @tda: physical number of columns which may differ from the corresponding dimension of the matrix.
 *
 * This function returns a #NcmMatrix of the array @d with a physical number of columns tda which may differ
 * from the corresponding dimension of the matrix. The matrix has @nrows rows and @ncols columns, and the physical
 * number of columns in memory is given by tda.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_data_static_tda (gdouble *d, gsize nrows, gsize ncols, gsize tda)
{
  NcmMatrix *cm = g_object_new (NCM_TYPE_MATRIX, NULL);
  cm->mv = gsl_matrix_view_array_with_tda (d, nrows, ncols, tda);
  cm->a = NULL;
  cm->type = NCM_MATRIX_DERIVED;
  return cm;
}

/**
 * ncm_matrix_get_submatrix:
 * @cm: a #NcmMatrix.
 * @k1: row index of the original matrix @cm.
 * @k2: column index of the original matrix @cm.
 * @nrows: number of rows of the submatrix.
 * @ncols: number of columns of the submatrix.
 *
 * This function returns a submatrix #NcmMatrix of the matrix @cm.
 * The upper-left element of the submatrix is the element (@k1,@k2) of the original matrix.
 * The submatrix has @nrows rows and @ncols columns.
 *
 * Returns: (transfer full): A #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_get_submatrix (NcmMatrix *cm, gsize k1, gsize k2, gsize nrows, gsize ncols)
{
  NcmMatrix *scm = g_object_new (NCM_TYPE_MATRIX, NULL);
  scm->mv = gsl_matrix_submatrix (NCM_MATRIX_GSL (cm), k1, k2, nrows, ncols);
  scm->pobj = G_OBJECT (cm);
  g_object_ref (cm);
  scm->type = NCM_MATRIX_DERIVED;

  return scm;
}

/**
 * ncm_matrix_get_col:
 * @cm: a #NcmMatrix.
 * @col: column index.
 *
 * This function returns the elements of the @col column of the matrix @cm
 * into a #NcmVector.
 *
 * Returns: (transfer full): A #NcmVector.
 */
NcmVector *
ncm_matrix_get_col (NcmMatrix *cm, const gsize col)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);
  cv->vv = gsl_matrix_column (NCM_MATRIX_GSL (cm), col);
  cv->a = NULL;
  cv->type = NCM_VECTOR_DERIVED;
  cv->pobj = G_OBJECT (cm);
  g_object_ref (cm);

  return cv;
}

/**
 * ncm_matrix_get_row:
 * @cm: a #NcmMatrix.
 * @row: row index.
 *
 * This function returns the elements of the @row row of the matrix @cm
 * into a #NcmVector.
 *
 * Returns: (transfer full): A #NcmVector.
 */
NcmVector *
ncm_matrix_get_row (NcmMatrix *cm, const gsize row)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);
  cv->vv = gsl_matrix_row (NCM_MATRIX_GSL (cm), row);
  cv->a = NULL;
  cv->type = NCM_VECTOR_DERIVED;
  cv->pobj = G_OBJECT (cm);
  g_object_ref (cm);

  return cv;
}

/**
 * ncm_matrix_free:
 * @cm: a #NcmMatrix.
 *
 * Atomically decrements the reference count of @cm by one. If the reference count drops to 0,
 * all memory allocated by @cm is released.
 *
 */
void
ncm_matrix_free (NcmMatrix *cm)
{
  if (g_object_is_floating (cm))
    g_object_ref_sink (cm);
  g_object_unref (cm);
}

/**
 * ncm_matrix_clear:
 * @cm: a #NcmMatrix.
 *
 * Atomically decrements the reference count of @cm by one. If the reference count drops to 0,
 * all memory allocated by @cm is released. The pointer is set to NULL.
 *
 */
void
ncm_matrix_clear (NcmMatrix **cm)
{
  if (*cm != NULL && g_object_is_floating (*cm))
    g_object_ref_sink (*cm);
  g_clear_object (cm);
}

static void
_ncm_matrix_dispose (GObject *object)
{
  NcmMatrix *cm = NCM_MATRIX (object);
  
  if (cm->a != NULL)
  {
    g_array_unref (cm->a);
    cm->a = NULL;
  }

  g_clear_object (&cm->pobj);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_matrix_parent_class)->dispose (object);
}

static void
_ncm_matrix_finalize (GObject *object)
{
  NcmMatrix *cm = NCM_MATRIX (object);
  switch (cm->type)
  {
	case NCM_MATRIX_SLICE:
	  g_slice_free1 (sizeof(gdouble) * NCM_MATRIX_NROWS (cm) * NCM_MATRIX_NCOLS (cm), NCM_MATRIX_DATA (cm));
	  NCM_MATRIX_DATA (cm) = NULL;
	  break;
	case NCM_MATRIX_GARRAY:
	  break;
	case NCM_MATRIX_MALLOC:
	  g_free (NCM_MATRIX_DATA (cm));
	  NCM_MATRIX_DATA (cm) = NULL;
	  break;
	case NCM_MATRIX_GSL_MATRIX:
	  gsl_matrix_free (cm->gm);
	  cm->gm = NULL;
	  break;
	case NCM_MATRIX_DERIVED:
	  NCM_MATRIX_DATA (cm) = NULL;
	  break;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_matrix_parent_class)->finalize (object);
}


/**
 * ncm_matrix_dup:
 * @cm: a constant #NcmMatrix
 *
 * Duplicates @cm setting the same values of the original propertities.
 *
 * Returns: (transfer full): A #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_dup (const NcmMatrix *cm)
{
  NcmMatrix *cm_cp = ncm_matrix_new (NCM_MATRIX_COL_LEN (cm), NCM_MATRIX_ROW_LEN (cm));
  ncm_matrix_memcpy (cm_cp, cm);
  return cm_cp;
}

/**
 * ncm_matrix_new_gsl_const: (skip)
 * @m: matrix from GNU Scientific Library (GSL).
 *
 * This function converts @m into a constant #NcmMatrix.
 *
 * Returns: A new constant #NcmMatrix.
 */
/**
 * ncm_matrix_get:
 * @cm: a constant #NcmMatrix.
 * @i: row index.
 * @j: column index.
 *
 *
 * Returns: The (@i,@j)-th element of the matrix @cm.
 */
/**
 * ncm_matrix_ptr:
 * @cm: a #NcmMatrix.
 * @i: row index.
 * @j: column index.
 *
 * Returns: A pointer to the (@i,@j)-th element of the matrix @cm.
 */
/**
 * ncm_matrix_ref:
 * @cm: a #NcmMatrix.
 *
 * Increase the reference count of @cm by one.
 *
 * Returns: (transfer full): @cm
 */
/**
 * ncm_matrix_set:
 * @cm: a #NcmMatrix.
 * @i: row index.
 * @j: column index.
 * @val: a double.
 *
 * This function sets the value of the (@i,@j)-th element of the matrix @cm to @val.
 *
 */
/**
 * ncm_matrix_transpose:
 * @cm: a #NcmMatrix.
 *
 * This function replaces the matrix @cm by its transpose by copying the elements of the matrix in-place.
 * The matrix must be square for this operation to be possible.
 *
 */
/**
 * ncm_matrix_set_identity:
 * @cm: a #NcmMatrix.
 *
 * This function sets the elements of the matrix @cm to the corresponding elements of the identity matrix,
 * i.e. a unit diagonal with all off-diagonal elements zero. This applies to both square and rectangular matrices.
 *
 */
/**
 * ncm_matrix_set_zero:
 * @cm: a #NcmMatrix.
 *
 * This function sets all the elements of the matrix @cm to zero.
 *
 */
/**
 * ncm_matrix_scale:
 * @cm: a #NcmMatrix.
 * @val: a double.
 *
 * This function multiplies the elements of the matrix @cm by the constant factor @val.
 * The result is stored in @cm.
 *
 */
/**
 * ncm_matrix_memcpy:
 * @cm1: a #NcmMatrix.
 * @cm2: a #NcmMatrix.
 *
 * This function copies the elements of the matrix @cm1 into the matrix @cm2.
 * The two matrices must have the same size.
 *
 */
/**
 * ncm_matrix_set_col:
 * @cm: a #NcmMatrix.
 * @n: column index.
 * @cv: a constant #NcmVector.
 *
 * This function copies the elements of the vector @cv into the @n-th column of the matrix @cm.
 * The length of the vector must be the same as the length of the column.
 *
 */
/**
 * ncm_matrix_get_array:
 * @cm: a #NcmMatrix.
 *
 * FIXME
 *
 * Returns: (transfer container) (element-type double): FIXME
 */

static void
ncm_matrix_init (NcmMatrix *m)
{
  memset (&m->mv, 0, sizeof (gsl_matrix_view));
  m->a = NULL;
  m->gm = NULL;
  m->pobj = NULL;
  m->type = 0;
}

static void
ncm_matrix_class_init (NcmMatrixClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  object_class->dispose = &_ncm_matrix_dispose;
  object_class->finalize = &_ncm_matrix_finalize;
}
