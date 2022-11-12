/***************************************************************************
 *            ncm_matrix.c
 *
 *  Thu January 05 20:18:45 2012
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
 * SECTION:ncm_matrix
 * @title: NcmMatrix
 * @short_description: Matrix object representing an array of doubles.
 *
 * This object defines the functions for allocating and accessing matrices.
 * Also includes several matrix operations.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_matrix.h"
#include "math/ncm_vector.h"
#include "math/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_linalg.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_VALS,
};

G_DEFINE_TYPE (NcmMatrix, ncm_matrix, G_TYPE_OBJECT);

static void
ncm_matrix_init (NcmMatrix *m)
{
  memset (&m->mv, 0, sizeof (gsl_matrix_view));
  m->pdata = NULL;
  m->pfree = NULL;
  m->type  = 0;
}

static void
_ncm_matrix_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMatrix *m = NCM_MATRIX (object);
  
  g_return_if_fail (NCM_IS_MATRIX (object));
  
  switch (prop_id)
  {
    case PROP_VALS:
    {
      GVariant *var = ncm_matrix_get_variant (m);
      
      g_value_take_variant (value, var);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_matrix_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMatrix *m = NCM_MATRIX (object);
  
  g_return_if_fail (NCM_IS_MATRIX (object));
  
  switch (prop_id)
  {
    case PROP_VALS:
    {
      GVariant *var = g_value_get_variant (value);
      
      ncm_matrix_set_from_variant (m, var);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_matrix_dispose (GObject *object)
{
  NcmMatrix *cm = NCM_MATRIX (object);
  
  if (cm->pdata != NULL)
  {
    g_assert (cm->pfree != NULL);
    cm->pfree (cm->pdata);
    cm->pdata = NULL;
    cm->pfree = NULL;
  }
  
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
      g_slice_free1 (sizeof (gdouble) * ncm_matrix_nrows (cm) * ncm_matrix_ncols (cm), ncm_matrix_data (cm));
      break;
    case NCM_MATRIX_GARRAY:
    case NCM_MATRIX_MALLOC:
    case NCM_MATRIX_GSL_MATRIX:
    case NCM_MATRIX_DERIVED:
      break;
    default:
      g_assert_not_reached ();
      break;
  }
  
  cm->mv.matrix.data = NULL;
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_matrix_parent_class)->finalize (object);
}

static void
ncm_matrix_class_init (NcmMatrixClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  
  object_class->set_property = &_ncm_matrix_set_property;
  object_class->get_property = &_ncm_matrix_get_property;
  object_class->dispose      = &_ncm_matrix_dispose;
  object_class->finalize     = &_ncm_matrix_finalize;

  /**   
   * NcmMatrix:values:
   *
   * GVariant representation of the matrix used to serialize the object.
   * 
   */
  g_object_class_install_property (object_class, PROP_VALS,
                                   g_param_spec_variant ("values", NULL, "values",
                                                         G_VARIANT_TYPE_ARRAY, NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_matrix_new:
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * This function allocates memory for a new #NcmMatrix of doubles
 * with @nrows rows and @ncols columns.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new (const guint nrows, const guint ncols)
{
  gdouble *d    = g_slice_alloc (sizeof (gdouble) * nrows * ncols);
  NcmMatrix *cm = ncm_matrix_new_full (d, nrows, ncols, ncols, NULL, NULL);
  
  cm->type = NCM_MATRIX_SLICE;
  
  return cm;
}

/**
 * ncm_matrix_new0:
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * This function allocates memory for a new #NcmMatrix of doubles
 * with @nrows rows and @ncols columns and sets all elements to zero.
 *
 * Returns: (transfer full): A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new0 (const guint nrows, const guint ncols)
{
  gdouble *d    = g_slice_alloc0 (sizeof (gdouble) * nrows * ncols);
  NcmMatrix *cm = ncm_matrix_new_full (d, nrows, ncols, ncols, NULL, NULL);
  
  cm->type = NCM_MATRIX_SLICE;
  
  return cm;
}

/**
 * ncm_matrix_new_full:
 * @d: pointer to the data
 * @nrows: number of rows
 * @ncols: number of columns
 * @tda: row trailing dimension
 * @pdata: (allow-none): descending data pointer
 * @pfree: (scope notified) (allow-none): free function to be called when destroying the matrix
 *
 * This function allocates memory for a new #NcmMatrix of doubles
 * with @nrows rows and @ncols columns.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_full (gdouble *d, guint nrows, guint ncols, guint tda, gpointer pdata, GDestroyNotify pfree)
{
  NcmMatrix *cm = g_object_new (NCM_TYPE_MATRIX, NULL);
  
  if (tda == ncols)
    cm->mv = gsl_matrix_view_array (d, nrows, ncols);
  else
    cm->mv = gsl_matrix_view_array_with_tda (d, nrows, ncols, tda);
  
  cm->type = NCM_MATRIX_DERIVED;
  g_assert ((pdata == NULL) || (pdata != NULL && pfree != NULL));
  cm->pdata = pdata;
  cm->pfree = pfree;
  
  return cm;
}

/**
 * ncm_matrix_new_gsl: (skip)
 * @gm: matrix from [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) to be converted into a #NcmMatrix
 *
 * This function saves @gm internally and frees it when it is no longer necessary.
 * The @gm matrix must not be freed.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_gsl (gsl_matrix *gm)
{
  NcmMatrix *cm = ncm_matrix_new_full (gm->data, gm->size1, gm->size2, gm->tda,
                                       gm, (GDestroyNotify) & gsl_matrix_free);
  
  cm->type = NCM_MATRIX_GSL_MATRIX;
  
  return cm;
}

/**
 * ncm_matrix_new_gsl_static: (skip)
 * @gm: matrix from [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) to be converted into a #NcmMatrix
 *
 * This function saves @gm internally and does not frees it.
 * The @gm matrix must be valid during the life of the created #NcmMatrix.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_gsl_static (gsl_matrix *gm)
{
  NcmMatrix *cm = ncm_matrix_new_full (gm->data, gm->size1, gm->size2, gm->tda,
                                       NULL, NULL);
  
  return cm;
}

/**
 * ncm_matrix_new_array:
 * @a: (array) (element-type double): GArray of doubles to be converted into a #NcmMatrix
 * @ncols: number of columns
 *
 * The number of rows is defined dividing the lenght of @a by @ncols.
 * This function saves @a internally and frees it when it is no longer necessary.
 * The GArray @a must not be freed.
 *
 * Returns: (transfer full): A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_array (GArray *a, guint ncols)
{
  g_assert_cmpuint (a->len % ncols, ==, 0);
  {
    NcmMatrix *cm = ncm_matrix_new_full (&g_array_index (a, gdouble, 0),
                                         a->len / ncols, ncols, ncols,
                                         g_array_ref (a),
                                         (GDestroyNotify) & g_array_unref);
    
    cm->type = NCM_MATRIX_GARRAY;
    
    return cm;
  }
}

/**
 * ncm_matrix_new_data_slice: (skip)
 * @d: pointer to the first double allocated
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * This function returns a #NcmMatrix of the array @d allocated using g_slice function.
 * It saves @d internally and frees it when it is no longer necessary.
 * The matrix has @nrows rows and @ncols columns.
 * The physical number of columns in memory is also given by @ncols.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_data_slice (gdouble *d, guint nrows, guint ncols)
{
  NcmMatrix *cm = ncm_matrix_new_full (d, nrows, ncols, ncols,
                                       NULL, NULL);
  
  cm->type = NCM_MATRIX_SLICE;
  
  return cm;
}

/**
 * ncm_matrix_new_data_malloc: (skip)
 * @d: pointer to the first double allocated
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * This function returns a #NcmMatrix of the array @d allocated using malloc.
 * It saves @d internally and frees it when it is no longer necessary.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_data_malloc (gdouble *d, guint nrows, guint ncols)
{
  NcmMatrix *cm = ncm_matrix_new_full (d, nrows, ncols, ncols,
                                       d, &g_free);
  
  cm->type = NCM_MATRIX_MALLOC;
  
  return cm;
}

/**
 * ncm_matrix_new_data_static: (skip)
 * @d: pointer to the first double allocated
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * This function returns a #NcmMatrix of the array @d.
 * The memory allocated is kept during all time life of the object and
 * must not be freed during this period.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_data_static (gdouble *d, guint nrows, guint ncols)
{
  NcmMatrix *cm = ncm_matrix_new_full (d, nrows, ncols, ncols,
                                       NULL, NULL);
  
  cm->type = NCM_MATRIX_DERIVED;
  
  return cm;
}

/**
 * ncm_matrix_new_data_static_tda: (skip)
 * @d: pointer to the first double allocated
 * @nrows: number of rows
 * @ncols: number of columns
 * @tda: physical number of columns which may differ from the corresponding dimension of the matrix
 *
 * This function returns a #NcmMatrix of the array @d with a physical number of columns tda which may differ
 * from the corresponding dimension of the matrix. The matrix has @nrows rows and @ncols columns, and the physical
 * number of columns in memory is given by tda.
 *
 * Returns: A new #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_new_data_static_tda (gdouble *d, guint nrows, guint ncols, guint tda)
{
  NcmMatrix *cm = ncm_matrix_new_full (d, nrows, ncols, tda,
                                       NULL, NULL);
  
  cm->type = NCM_MATRIX_DERIVED;
  
  return cm;
}

/**
 * ncm_matrix_new_variant:
 * @var: a variant of type "aad"
 *
 * Creates a new matrix using the values from @var.
 *
 * Returns: (transfer full): a #NcmMatrix with the values from @var.
 */
NcmMatrix *
ncm_matrix_new_variant (GVariant *var)
{
  NcmMatrix *cm = g_object_new (NCM_TYPE_MATRIX,
                                "values", var,
                                NULL);
  
  return cm;
}

/**
 * ncm_matrix_const_new_data:
 * @d: pointer to the first double allocated
 * @nrows: number of rows
 * @ncols: number of cols
 *
 * This function returns a constant #NcmMatrix of the array @d.
 * The memory allocated is kept during all time life of the object and
 * must not be freed during this period.
 *
 * Returns: A new constant #NcmMatrix.
 */
const NcmMatrix *
ncm_matrix_const_new_data (const gdouble *d, guint nrows, guint ncols)
{
  NcmMatrix *cm = g_object_new (NCM_TYPE_MATRIX, NULL);
  
  cm->mv   = gsl_matrix_view_array ((gdouble *) d, nrows, ncols);
  cm->type = NCM_MATRIX_DERIVED;
  
  return cm;
}

/**
 * ncm_matrix_const_new_variant:
 * @var: a variant of type "aad"
 *
 * Creates a new constant matrix using the same memory of @var.
 *
 * Returns: (transfer full): a #NcmMatrix with the values from @var.
 */
const NcmMatrix *
ncm_matrix_const_new_variant (GVariant *var)
{
  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE ("aad")));
  {
    GVariant *row      = g_variant_get_child_value (var, 0);
    guint nrows        = g_variant_n_children (var);
    guint ncols        = g_variant_n_children (row);
    gconstpointer data = g_variant_get_data (var);
    const NcmMatrix *m = ncm_matrix_const_new_data (data, nrows, ncols);
    
    NCM_MATRIX (m)->pdata = g_variant_ref_sink (var);
    NCM_MATRIX (m)->pfree = (GDestroyNotify) & g_variant_unref;
    
    return m;
  }
}

/**
 * ncm_matrix_ref:
 * @cm: a #NcmMatrix
 *
 * Increase the reference count of @cm by one.
 *
 * Returns: (transfer full): @cm.
 */
NcmMatrix *
ncm_matrix_ref (NcmMatrix *cm)
{
  return g_object_ref (cm);
}

/**
 * ncm_matrix_get_submatrix:
 * @cm: a #NcmMatrix
 * @k1: row index of the original matrix @cm
 * @k2: column index of the original matrix @cm
 * @nrows: number of rows of the submatrix
 * @ncols: number of columns of the submatrix
 *
 * This function returns a submatrix #NcmMatrix of the matrix @cm.
 * The upper-left element of the submatrix is the element (@k1,@k2) of the original matrix.
 * The submatrix has @nrows rows and @ncols columns.
 *
 * Returns: (transfer full): A #NcmMatrix.
 */
NcmMatrix *
ncm_matrix_get_submatrix (NcmMatrix *cm, guint k1, guint k2, guint nrows, guint ncols)
{
  NcmMatrix *scm = g_object_new (NCM_TYPE_MATRIX, NULL);
  
  g_assert_cmpuint (nrows + k1, <=, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncols + k2, <=, ncm_matrix_ncols (cm));
  
  scm->mv = gsl_matrix_submatrix (ncm_matrix_gsl (cm), k1, k2, nrows, ncols);
  
  scm->pdata = g_object_ref (cm);
  scm->pfree = g_object_unref;
  scm->type  = NCM_MATRIX_DERIVED;
  
  return scm;
}

/**
 * ncm_matrix_get_col:
 * @cm: a #NcmMatrix
 * @col: column index
 *
 * This function returns the elements of the @col column of the matrix @cm
 * into a #NcmVector.
 *
 * Returns: (transfer full): A #NcmVector.
 */
NcmVector *
ncm_matrix_get_col (NcmMatrix *cm, const guint col)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);
  
  cv->vv   = gsl_matrix_column (ncm_matrix_gsl (cm), col);
  cv->type = NCM_VECTOR_DERIVED;
  
  cv->pdata = g_object_ref (cm);
  cv->pfree = &g_object_unref;
  
  return cv;
}

/**
 * ncm_matrix_get_row:
 * @cm: a #NcmMatrix
 * @row: row index
 *
 * This function returns the elements of the @row row of the matrix @cm
 * into a #NcmVector.
 *
 * Returns: (transfer full): A #NcmVector.
 */
NcmVector *
ncm_matrix_get_row (NcmMatrix *cm, const guint row)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);
  
  cv->vv   = gsl_matrix_row (ncm_matrix_gsl (cm), row);
  cv->type = NCM_VECTOR_DERIVED;
  
  cv->pdata = g_object_ref (cm);
  cv->pfree = &g_object_unref;
  
  return cv;
}

/**
 * ncm_matrix_as_vector:
 * @cm: a #NcmMatrix
 *
 * Creates a vector containing the row-wise concatenation
 * of the matrix @cm. It requires a matrix with tda==ncols.
 *
 * Returns: (transfer full): A #NcmVector.
 */
NcmVector *
ncm_matrix_as_vector (NcmMatrix *cm)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  const guint len   = nrows * ncols;

  NcmVector *v = ncm_vector_new_full (ncm_matrix_data (cm), len, 1, ncm_matrix_ref (cm), (GDestroyNotify)&ncm_matrix_free);

  g_assert (ncm_matrix_tda (cm) == ncm_matrix_ncols (cm));

  return v;
}

/**
 * ncm_matrix_set_from_variant:
 * @cm: a #NcmMatrix
 * @var: a #GVariant of type "aad"
 *
 * This function sets the values of @cm using the variant @var.
 *
 */
void
ncm_matrix_set_from_variant (NcmMatrix *cm, GVariant *var)
{
  if (!g_variant_is_of_type (var, G_VARIANT_TYPE ("aad")))
  {
    g_error ("ncm_matrix_set_from_variant: Cannot convert `%s' variant to an array of arrays of doubles",
             g_variant_get_type_string (var));
  }
  else
  {
    GVariant *row = g_variant_get_child_value (var, 0);
    guint nrows   = g_variant_n_children (var);
    guint ncols   = g_variant_n_children (row);
    guint i;
    
    g_variant_unref (row);
    
    /* Sometimes we receive a NcmMatrix in the process of instantiation. */
    if ((ncm_matrix_nrows (cm) == 0) && (ncm_matrix_ncols (cm) == 0))
    {
      gdouble *d = g_slice_alloc (sizeof (gdouble) * nrows * ncols);
      
      cm->mv   = gsl_matrix_view_array (d, nrows, ncols);
      cm->type = NCM_MATRIX_SLICE;
    }
    else if ((nrows != ncm_matrix_nrows (cm)) || (ncols != ncm_matrix_ncols (cm)))
    {
      g_error ("ncm_matrix_set_from_variant: cannot set matrix values, variant contains (%u, %u) childs but matrix dimension is (%u, %u)", nrows, ncols, ncm_matrix_nrows (cm), ncm_matrix_ncols (cm));
    }
    
    for (i = 0; i < nrows; i++)
    {
      NcmVector *m_row = ncm_matrix_get_row (cm, i);
      
      row = g_variant_get_child_value (var, i);
      
      {
        gsize v_ncols           = 0;
        const gdouble *row_data = g_variant_get_fixed_array (row, &v_ncols, sizeof (gdouble));
        
        g_assert_cmpuint (v_ncols, ==, ncols);
        ncm_vector_set_data (m_row, row_data, ncols);
        
        g_variant_unref (row);
        ncm_vector_free (m_row);
      }
    }
  }
}

/**
 * ncm_matrix_get_variant:
 * @cm: a #NcmMatrix
 *
 * This function gets a variant of values taken from @cm.
 *
 * Returns: (transfer full): the newly created #GVariant.
 */
GVariant *
ncm_matrix_get_variant (NcmMatrix *cm)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  GVariant **rows   = g_new (GVariant *, nrows);
  GVariant *var;
  guint i;
  
  for (i = 0; i < nrows; i++)
  {
    rows[i] = g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), ncm_matrix_ptr (cm, i, 0), ncols, sizeof (gdouble));
  }
  
  var = g_variant_new_array (G_VARIANT_TYPE ("ad"), rows, nrows);
  g_free (rows);
  g_variant_ref_sink (var);
  
  return var;
}

/**
 * ncm_matrix_peek_variant:
 * @cm: a #NcmMatrix
 *
 * This function gets a variant of values taken from @cm using the same memory.
 * The matrix @cm should not be modified during the variant existance.
 *
 * Returns: (transfer full): the newly created #GVariant.
 */
GVariant *
ncm_matrix_peek_variant (NcmMatrix *cm)
{
  guint nrows = ncm_matrix_nrows (cm);
  guint ncols = ncm_matrix_ncols (cm);
  GVariant *var;
  GVariant **rows = g_new (GVariant *, nrows);
  guint row_size  = ncols * sizeof (gdouble);
  guint i         = 0;
  
  rows[i] = g_variant_new_from_data (G_VARIANT_TYPE ("ad"),
                                     ncm_matrix_ptr (cm, i, 0),
                                     row_size,
                                     TRUE,
                                     (GDestroyNotify) & ncm_matrix_free,
                                     ncm_matrix_ref (cm));
  
  for (i = 1; i < nrows; i++)
  {
    rows[i] = g_variant_new_from_data (G_VARIANT_TYPE ("ad"),
                                       ncm_matrix_ptr (cm, i, 0),
                                       row_size, TRUE, NULL, NULL);
  }
  
  var = g_variant_new_array (G_VARIANT_TYPE ("ad"), rows, nrows);
  
  g_free (rows);
  
  return g_variant_ref_sink (var);
}

/**
 * ncm_matrix_set_from_data:
 * @cm: a #NcmMatrix
 * @data: (array) (element-type double): Array of doubles
 *
 * This function sets the valuus of @cm using @data. Data
 * must have the same size as #NcmMatrix.
 *
 */
void
ncm_matrix_set_from_data (NcmMatrix *cm, gdouble *data)
{
  guint nrows = ncm_matrix_nrows (cm);
  guint ncols = ncm_matrix_ncols (cm);
  guint i, j;
  
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
    {
      ncm_matrix_set (cm, i, j, data[i * ncols + j]);
    }
  }
}

/**
 * ncm_matrix_set_from_array:
 * @cm: a #NcmMatrix
 * @a: (array) (element-type double): Array of doubles
 *
 * This function sets the valuus of @cm using @data. Data
 * must have the same size as #NcmMatrix.
 *
 */
void
ncm_matrix_set_from_array (NcmMatrix *cm, GArray *a)
{
  guint nrows = ncm_matrix_nrows (cm);
  guint ncols = ncm_matrix_ncols (cm);
  guint i, j;
  
  g_assert_cmpuint (nrows * ncols, ==, a->len);
  
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
    {
      ncm_matrix_set (cm, i, j, g_array_index (a, gdouble, i * ncols + j));
    }
  }
}

/**
 * ncm_matrix_free:
 * @cm: a #NcmMatrix
 *
 * Atomically decrements the reference count of @cm by one. If the reference count drops to 0,
 * all memory allocated by @cm is released.
 *
 */
void
ncm_matrix_free (NcmMatrix *cm)
{
  g_object_unref (cm);
}

/**
 * ncm_matrix_clear:
 * @cm: a #NcmMatrix
 *
 * Atomically decrements the reference count of @cm by one. If the reference count drops to 0,
 * all memory allocated by @cm is released. The pointer is set to NULL.
 *
 */
void
ncm_matrix_clear (NcmMatrix **cm)
{
  g_clear_object (cm);
}

/**
 * ncm_matrix_const_free:
 * @cm: a constant #NcmMatrix
 *
 * Atomically decrements the reference count of @cv by one. If the reference count drops to 0,
 * all memory allocated by @cv is released.
 *
 */
void
ncm_matrix_const_free (const NcmMatrix *cm)
{
  ncm_matrix_free (NCM_MATRIX (cm));
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
  NcmMatrix *cm_cp = ncm_matrix_new (ncm_matrix_col_len (cm), ncm_matrix_row_len (cm));
  
  ncm_matrix_memcpy (cm_cp, cm);
  
  return cm_cp;
}

/**
 * ncm_matrix_substitute:
 * @cm: a #NcmMatrix
 * @nm: (allow-none): a #NcmMatrix
 * @check_size: a boolean
 *
 * Substitute the matrix *@cm by @nm, first it unref *@cm if it is not NULL.
 * If @check_size is TRUE then check if the two matrix have the same size.
 *
 */
void
ncm_matrix_substitute (NcmMatrix **cm, NcmMatrix *nm, gboolean check_size)
{
  if (*cm == nm)
    return;
  
  if (*cm != NULL)
  {
    if ((nm != NULL) && check_size)
    {
      g_assert_cmpuint (ncm_matrix_nrows (*cm), ==, ncm_matrix_nrows (nm));
      g_assert_cmpuint (ncm_matrix_ncols (*cm), ==, ncm_matrix_ncols (nm));
    }
    
    ncm_matrix_clear (cm);
  }
  
  if (nm != NULL)
    *cm = ncm_matrix_ref (nm);
}

/**
 * ncm_matrix_add_mul:
 * @cm1: a #NcmMatrix
 * @a: a constant gdouble
 * @cm2: a #NcmMatrix
 *
 * This function performs the operation @cm1 = @a $\times$ @cm2 $+$ @cm1.
 *
 */
void
ncm_matrix_add_mul (NcmMatrix *cm1, const gdouble a, NcmMatrix *cm2)
{
  const gboolean no_pad_cm1 = ncm_matrix_gsl (cm1)->tda == ncm_matrix_ncols (cm1);
  const gboolean no_pad_cm2 = ncm_matrix_gsl (cm2)->tda == ncm_matrix_ncols (cm2);
  
  g_assert (ncm_matrix_ncols (cm1) == ncm_matrix_ncols (cm2));
  g_assert (ncm_matrix_nrows (cm1) == ncm_matrix_nrows (cm2));
  
  if (no_pad_cm1 && no_pad_cm2)
  {
    const gint N = ncm_matrix_ncols (cm2) * ncm_matrix_nrows (cm2);
    
    cblas_daxpy (N, a,
                 ncm_matrix_data (cm2), 1,
                 ncm_matrix_data (cm1), 1);
  }
  else
  {
    const gint N        = ncm_matrix_row_len (cm2);
    const guint cm2_tda = ncm_matrix_gsl (cm2)->tda;
    const guint cm1_tda = ncm_matrix_gsl (cm1)->tda;
    guint i;
    
    for (i = 0; i < ncm_matrix_nrows (cm2); i++)
    {
      cblas_daxpy (N, a,
                   &ncm_matrix_data (cm2)[cm2_tda * i], 1,
                   &ncm_matrix_data (cm1)[cm1_tda * i], 1);
    }
  }
}

/**
 * ncm_matrix_cmp:
 * @cm1: a constant #NcmMatrix
 * @cm2: a constant #NcmMatrix
 * @scale: a constant gdouble
 *
 * This function performes a comparison, component-wise, of the two matrices, given by $\left| \right. ($@cm1 $-$ @cm2 $)/($@scale $+$ @cm2$)\left. \right|$ and returns its maximum value.
 *
 * Returns: The maximum value of the operation $\left| \right. ($@cm1 $-$ @cm2 $)/($@scale $+$ @cm2$)\left. \right|$.
 */
gdouble
ncm_matrix_cmp (const NcmMatrix *cm1, const NcmMatrix *cm2, const gdouble scale)
{
  const guint nrows = ncm_matrix_nrows (cm1);
  const guint ncols = ncm_matrix_ncols (cm1);
  gdouble reltol = 0.0;
  gint i, j;
  
  g_assert_cmpuint (ncols, ==, ncm_matrix_ncols (cm2));
  g_assert_cmpuint (nrows, ==, ncm_matrix_nrows (cm2));
  
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
    {
      const gdouble cm1_ij    = ncm_matrix_get (cm1, i, j);
      const gdouble cm2_ij    = ncm_matrix_get (cm2, i, j);
      const gdouble reltol_ij = fabs ((cm1_ij - cm2_ij) / (scale + cm2_ij));
      
      reltol = MAX (reltol, reltol_ij);
    }
  }
  
  return reltol;
}

/**
 * ncm_matrix_cmp_diag:
 * @cm1: a #NcmMatrix
 * @cm2: a #NcmMatrix
 * @scale: a constant gdouble
 *
 * This function is similar to ncm_matrix_cmp(), but now only the diagonal elements are compared.
 *
 * Returns: The maximum value of the operation $\left| \right. ($@cm1 $-$ @cm2 $)/($@scale $+$ @cm2$)\left. \right|$ comparing only the diagonal elements.
 */
gdouble
ncm_matrix_cmp_diag (const NcmMatrix *cm1, const NcmMatrix *cm2, const gdouble scale)
{
  const guint nrows = ncm_matrix_nrows (cm1);
  const guint ncols = ncm_matrix_ncols (cm1);
  const gdouble len = MIN (nrows, ncols);
  gdouble reltol    = 0.0;
  gint i;
  
  g_assert_cmpuint (ncols, ==, ncm_matrix_ncols (cm2));
  g_assert_cmpuint (nrows, ==, ncm_matrix_nrows (cm2));
  
  for (i = 0; i < len; i++)
  {
    const gdouble cm1_ii    = ncm_matrix_get (cm1, i, i);
    const gdouble cm2_ii    = ncm_matrix_get (cm2, i, i);
    const gdouble reltol_ii = fabs ((cm1_ii - cm2_ii) / (scale + cm2_ii));
    
    reltol = MAX (reltol, reltol_ii);
  }
  
  return reltol;
}

/**
 * ncm_matrix_copy_triangle:
 * @cm: a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix
 *
 * If @UL == 'U' copy the upper triangle over the lower.
 * If @UL == 'L' copy the lower triangle over the upper.
 *
 */
void
ncm_matrix_copy_triangle (NcmMatrix *cm, gchar UL)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  guint i, j;
  
  if (nrows != ncols)
    g_error ("ncm_matrix_copy_triangle: only works on square a matrix [%ux%u]", nrows, ncols);
  
  if ((UL != 'U') && (UL != 'L'))
    g_error ("ncm_matrix_copy_triangle: expect U or L and received %c.", UL);
  
  if (UL == 'U')
  {
    for (i = 0; i < nrows; i++)
    {
      for (j = i + 1; j < ncols; j++)
      {
        ncm_matrix_set (cm, j, i, ncm_matrix_get (cm, i, j));
      }
    }
  }
  else if (UL == 'L')
  {
    for (i = 0; i < nrows; i++)
    {
      for (j = i + 1; j < ncols; j++)
      {
        ncm_matrix_set (cm, i, j, ncm_matrix_get (cm, j, i));
      }
    }
  }
  else
  {
    g_assert_not_reached ();
  }
}

/**
 * ncm_matrix_dsymm:
 * @cm: a #NcmMatrix
 * @UL: gchar indicating 'U'pper or 'L'ower matrix
 * @alpha: a constant gdouble
 * @A: a #NcmMatrix
 * @B: a #NcmMatrix
 * @beta: a constant gdouble
 *
 * This function performes the following operation:
 *
 * if @UL == 'U': $\left( \mathsf{cm} \leftarrow \alpha \mathbf{A} \mathbf{B} + \beta \, \mathsf{cm} \right)$;
 *
 * if @UL == 'L': $\left( \mathsf{cm} \leftarrow \alpha \mathbf{B} \mathbf{A} + \beta \, \mathsf{cm} \right)$.
 *
 * Where $\mathbf{A} = \mathbf{A}^\intercal$.
 *
 */
void
ncm_matrix_dsymm (NcmMatrix *cm, gchar UL, const gdouble alpha, NcmMatrix *A, NcmMatrix *B, const gdouble beta)
{
  g_assert (UL == 'U' || UL == 'L');
  g_assert_cmpuint (ncm_matrix_ncols (A), ==, ncm_matrix_ncols (B));
  g_assert_cmpuint (ncm_matrix_nrows (A), ==, ncm_matrix_nrows (B));
  g_assert_cmpuint (ncm_matrix_ncols (A), ==, ncm_matrix_ncols (cm));
  g_assert_cmpuint (ncm_matrix_nrows (A), ==, ncm_matrix_nrows (cm));
  
  cblas_dsymm (CblasRowMajor, CblasLeft, (UL == 'U') ? CblasUpper : CblasLower, ncm_matrix_nrows (cm), ncm_matrix_ncols (cm),
               alpha,
               ncm_matrix_data (A), ncm_matrix_tda (A),
               ncm_matrix_data (B), ncm_matrix_tda (B),
               beta,
               ncm_matrix_data (cm), ncm_matrix_tda (cm));
}

CBLAS_TRANSPOSE
_ncm_matrix_check_trans (const gchar *func_name, gchar T)
{
  switch (T)
  {
    case 'C':
    case 'T':
    
      return CblasTrans;
      
    case 'N':
    
      return CblasNoTrans;
      
      break;
    default:
      g_error ("%s: Unknown Trans type %c.", func_name, T);
      
      return 0;
  }
}

/**
 * ncm_matrix_dgemm:
 * @cm: a #NcmMatrix $C$
 * @TransA: char indicating 'T'ranspose or 'N'ot transposed matrix
 * @TransB: char indicating 'T'ranspose or 'N'ot transposed matrix
 * @alpha: $\alpha$
 * @A: a #NcmMatrix $A$
 * @B: a #NcmMatrix $B$
 * @beta: $\beta$
 *
 * Calculates $C = \alpha\mathrm{op}(A)\mathrm{op}(B) + \beta C$.
 *
 */
void
ncm_matrix_dgemm (NcmMatrix *cm, gchar TransA, gchar TransB, const gdouble alpha, NcmMatrix *A, NcmMatrix *B, const gdouble beta)
{
  CBLAS_TRANSPOSE cblas_TransA = _ncm_matrix_check_trans ("ncm_matrix_dgemm", TransA);
  CBLAS_TRANSPOSE cblas_TransB = _ncm_matrix_check_trans ("ncm_matrix_dgemm", TransB);
  const gsize opA_nrows        = (cblas_TransA == CblasNoTrans) ? ncm_matrix_nrows (A) : ncm_matrix_ncols (A);
  const gsize opA_ncols        = (cblas_TransA == CblasNoTrans) ? ncm_matrix_ncols (A) : ncm_matrix_nrows (A);
  const gsize opB_nrows        = (cblas_TransB == CblasNoTrans) ? ncm_matrix_nrows (B) : ncm_matrix_ncols (B);
  const gsize opB_ncols        = (cblas_TransB == CblasNoTrans) ? ncm_matrix_ncols (B) : ncm_matrix_nrows (B);
  
  g_assert_cmpuint (opA_ncols, ==, opB_nrows);
  
  g_assert_cmpuint (opA_nrows, ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (opB_ncols, ==, ncm_matrix_ncols (cm));
  
  cblas_dgemm (CblasRowMajor, cblas_TransA, cblas_TransB, opA_nrows, opB_ncols, opA_ncols,
               alpha,
               ncm_matrix_data (A), ncm_matrix_tda (A),
               ncm_matrix_data (B), ncm_matrix_tda (B),
               beta,
               ncm_matrix_data (cm), ncm_matrix_tda (cm));
}

/**
 * ncm_matrix_cholesky_decomp:
 * @cm: a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix
 *
 * Calculates in-place the Cholesky decomposition for a symmetric positive
 * definite matrix.
 *
 */
gint
ncm_matrix_cholesky_decomp (NcmMatrix *cm, gchar UL)
{
  gint ret = ncm_lapack_dpotrf (UL, ncm_matrix_nrows (cm), ncm_matrix_data (cm), ncm_matrix_tda (cm));
  
  return ret;
}

/**
 * ncm_matrix_cholesky_inverse:
 * @cm: a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix
 *
 * Calculates inplace the inverse of @cm that has been previously decomposed by
 * the Cholesky decomposition ncm_matrix_cholesky_decomp().
 *
 */
gint
ncm_matrix_cholesky_inverse (NcmMatrix *cm, gchar UL)
{
  gint ret = ncm_lapack_dpotri (UL, ncm_matrix_nrows (cm), ncm_matrix_data (cm), ncm_matrix_tda (cm));
  
  return ret;
}

/**
 * ncm_matrix_cholesky_lndet:
 * @cm: a #NcmMatrix
 *
 * Calculates determinant of a symmetric positive definite matrix,
 * that was previously decomposed using ncm_matrix_cholesky_decomp().
 *
 * Returns: the log determinant of @cm.
 */
gdouble
ncm_matrix_cholesky_lndet (NcmMatrix *cm)
{
  const gdouble lb = 1.0e-200;
  const gdouble ub = 1.0e+200;
  const guint n    = ncm_matrix_nrows (cm);
  gdouble detL     = 1.0;
  glong exponent   = 0;
  guint i;
  
  for (i = 0; i < n; i++)
  {
    const gdouble Lii   = fabs (ncm_matrix_get (cm, i, i));
    const gdouble ndetL = detL * Lii;
    
    if (G_UNLIKELY ((ndetL < lb) || (ndetL > ub)))
    {
      gint exponent_i = 0;
      
      detL      = frexp (ndetL, &exponent_i);
      exponent += exponent_i;
    }
    else
    {
      detL = ndetL;
    }
  }
  
  return 2.0 * (log (detL) + exponent * M_LN2);
}

/**
 * ncm_matrix_cholesky_solve:
 * @cm: a #NcmMatrix
 * @b: a #NcmVector
 * @UL: char indicating 'U'pper or 'L'ower matrix
 *
 * Calculates in-place the Cholesky decomposition for a symmetric positive
 * definite matrix and solve the system $A x = B$ where $A=$@cm and $B$=@b.
 *
 */
gint
ncm_matrix_cholesky_solve (NcmMatrix *cm, NcmVector *b, gchar UL)
{
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_vector_len (b));
  g_assert_cmpuint (ncm_vector_stride (b), ==, 1);
  
  return ncm_lapack_dposv (UL, ncm_matrix_nrows (cm), 1,
                           ncm_matrix_data (cm), ncm_matrix_tda (cm),
                           ncm_vector_data (b),  ncm_vector_len (b));
}

/**
 * ncm_matrix_cholesky_solve2:
 * @cm: a #NcmMatrix
 * @b: a #NcmVector
 * @UL: char indicating 'U'pper or 'L'ower matrix
 *
 * Using a previously computed Cholesky decomposition in @cm, through
 * ncm_matrix_cholesky_decomp(), solves the system $A x = B$ where
 * $A=$@cm and $B$=@b.
 *
 */
gint
ncm_matrix_cholesky_solve2 (NcmMatrix *cm, NcmVector *b, gchar UL)
{
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_vector_len (b));
  g_assert_cmpuint (ncm_vector_stride (b), ==, 1);
  
  return ncm_lapack_dpotrs (UL, ncm_matrix_nrows (cm), 1,
                            ncm_matrix_data (cm), ncm_matrix_tda (cm),
                            ncm_vector_data (b),  ncm_vector_len (b));
}

/**
 * ncm_matrix_nearPD:
 * @cm: a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix
 * @cholesky_decomp: if true substitue @cm for its Cholesky decomposition
 * @maxiter: maximum number of iterations
 *
 * Assuming that @cm is a symmetric matrix with data on @UL
 * side, computes the nearest positive definite matrix
 * in the Frobenius norm. See [Higham (2002)][XHigham2002].
 * The iterations stop when the Cholesky decomposition is valid.
 *
 * Returns: the return value of the last Cholesky decomposition.
 */
gint
ncm_matrix_nearPD (NcmMatrix *cm, gchar UL, gboolean cholesky_decomp, const guint maxiter)
{
  const guint n = ncm_matrix_ncols (cm);
  NcmVector *eva = ncm_vector_new (n);
  NcmVector *diag = ncm_vector_new (n);
  NcmMatrix *eve = ncm_matrix_new (n, n);
  NcmMatrix *D_S = ncm_matrix_new (n, n);
  NcmMatrix *R = ncm_matrix_new (n, n);
  GArray *isuppz = g_array_new (FALSE, FALSE, sizeof (gint));
  NcmLapackWS *ws = ncm_lapack_ws_new ();
  gint neva = 0;
  gint ret, i, iter;
  
  g_array_set_size (isuppz, 2 * n);
  
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_matrix_nrows (cm));
  
  ncm_matrix_set_zero (D_S);
  ncm_matrix_get_diag (cm, diag);
  
  iter = 0;
  
  while (TRUE)
  {
    gdouble min_pos_ev = GSL_POSINF;
    
    ncm_matrix_sub (cm, D_S);
    
    ncm_matrix_memcpy (R, cm);
    
    ret = ncm_lapack_dsyevr ('V', 'A', UL, n, ncm_matrix_data (cm), ncm_matrix_tda (cm),
                             0.0, 0.0,
                             0, 0, 0.0,
                             &neva, ncm_vector_data (eva),
                             ncm_matrix_data (eve), ncm_matrix_tda (eve),
                             &g_array_index (isuppz, gint, 0),
                             ws);
    g_assert_cmpint (ret, ==, 0);
    
    if (neva == 0)
      g_error ("ncm_matrix_nearPD: matrix is negative semi-definite.");
    
    for (i = 0; i < neva; i++)
    {
      if (ncm_vector_get (eva, i) > 0.0)
        min_pos_ev = MIN (ncm_vector_get (eva, i), min_pos_ev);
    }
    
    for (i = 0; i < neva; i++)
    {
      if (ncm_vector_get (eva, i) < 0.0)
        ncm_vector_set (eva, i, min_pos_ev * GSL_DBL_EPSILON);
      
      cblas_dscal (n, sqrt (ncm_vector_get (eva, i)), ncm_matrix_ptr (eve, i, 0), 1);
      /*ncm_matrix_mul_row (eve, i, sqrt (ncm_vector_get (eva, i))); */
    }
    
    /*ncm_vector_log_vals (eva, "EVA: ", "% 22.15g", TRUE);*/
    /*ncm_matrix_log_vals (eve, "EVE: ", "% 22.15g");*/
    
    cblas_dsyrk (CblasRowMajor, (UL == 'U') ? CblasUpper : CblasLower,
                 CblasTrans, n, neva,
                 1.0, ncm_matrix_data (eve), ncm_matrix_tda (eve),
                 0.0, ncm_matrix_data (cm), ncm_matrix_tda (cm));
    
    ncm_matrix_memcpy (D_S, cm);
    ncm_matrix_sub (D_S, R);
    
    ncm_matrix_set_diag (cm, diag);
    
    ncm_matrix_memcpy (R, cm);
    
    if ((ret = ncm_matrix_cholesky_decomp (R, UL)) == 0)
      break;
    
    if (iter > maxiter)
      break;
    
    iter++;
    /*printf ("CHOLESKY: %4d, ITER %6d\n", ncm_matrix_cholesky_decomp (R, UL), iter++);*/
    /*ncm_matrix_log_vals (cm, "NewCM: ", "% 22.15g");*/
  }
  
  if (cholesky_decomp)
    ncm_matrix_memcpy (cm, R);
  
  g_array_unref (isuppz);
  ncm_lapack_ws_free (ws);
  ncm_vector_free (eva);
  ncm_vector_free (diag);
  ncm_matrix_free (eve);
  ncm_matrix_free (R);
  ncm_matrix_free (D_S);
  
  return ret;
}

/**
 * ncm_matrix_sym_exp_cholesky:
 * @cm: $M$ a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix
 * @exp_cm_dec: on exit this matrix contain the upper triangular matrix $U$ where $\exp(M) = U^\intercal U$
 *
 * Assuming that @cm is a symmetric matrix with data on @UL
 * side, computes the matrix exponential of @cm and its cholesky
 * decomposition.
 */
void
ncm_matrix_sym_exp_cholesky (NcmMatrix *cm, gchar UL, NcmMatrix *exp_cm_dec)
{
  const guint n = ncm_matrix_ncols (cm);
  NcmVector *eva = ncm_vector_new (n);
  NcmMatrix *eve = exp_cm_dec;
  GArray *tau = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *isuppz = g_array_new (FALSE, FALSE, sizeof (gint));
  NcmLapackWS *ws = ncm_lapack_ws_new ();
  gint neva = 0;
  gint ret, i;
  
  g_array_set_size (isuppz, 2 * n);
  g_array_set_size (tau, n);
  
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (exp_cm_dec), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (exp_cm_dec), ==, ncm_matrix_nrows (exp_cm_dec));
  
  ret = ncm_lapack_dsyevr ('V', 'A', UL, n, ncm_matrix_data (cm), ncm_matrix_tda (cm),
                           0.0, 0.0,
                           0, 0, 0.0,
                           &neva, ncm_vector_data (eva),
                           ncm_matrix_data (eve), ncm_matrix_tda (eve),
                           &g_array_index (isuppz, gint, 0),
                           ws);
  g_assert_cmpint (ret, ==, 0);
  
  for (i = 0; i < n; i++)
    cblas_dscal (n, exp (0.5 * ncm_vector_get (eva, i)), ncm_matrix_ptr (eve, i, 0), 1);
  
  ret = ncm_lapack_dgeqrf (n, n, ncm_matrix_data (eve), ncm_matrix_tda (eve), &g_array_index (tau, gdouble, 0), ws);
  g_assert_cmpint (ret, ==, 0);
  
  g_array_unref (isuppz);
  g_array_unref (tau);
  ncm_lapack_ws_free (ws);
  ncm_vector_free (eva);
}

/**
 * ncm_matrix_sym_posdef_log:
 * @cm: $M$ a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix
 * @ln_cm: on exit this matrix contain the upper triangular matrix $U$ where $\exp(M) = U^\intercal U$
 *
 * Assuming that @cm is a symmetric matrix with data on @UL
 * side, computes the matrix logarithm of @cm.
 */
void
ncm_matrix_sym_posdef_log (NcmMatrix *cm, gchar UL, NcmMatrix *ln_cm)
{
  const guint n = ncm_matrix_ncols (cm);
  NcmVector *eva = ncm_vector_new (n);
  NcmMatrix *eve = ncm_matrix_new (n, n);
  NcmMatrix *temp = ncm_matrix_new (n, n);
  GArray *isuppz = g_array_new (FALSE, FALSE, sizeof (gint));
  NcmLapackWS *ws = ncm_lapack_ws_new ();
  gint neva = 0;
  gint ret, i;
  
  g_array_set_size (isuppz, 2 * n);
  
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (ln_cm), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (ln_cm), ==, ncm_matrix_nrows (ln_cm));
  
  ret = ncm_lapack_dsyevr ('V', 'A', UL, n, ncm_matrix_data (cm), ncm_matrix_tda (cm),
                           0.0, 0.0,
                           0, 0, 0.0,
                           &neva, ncm_vector_data (eva),
                           ncm_matrix_data (eve), ncm_matrix_tda (eve),
                           &g_array_index (isuppz, gint, 0),
                           ws);
  g_assert_cmpint (ret, ==, 0);
  
  ncm_matrix_memcpy (temp, eve);
  
  for (i = 0; i < n; i++)
  {
    const gdouble e_val = ncm_vector_get (eva, i);
    
    if (e_val <= 0.0)
      g_error ("ncm_matrix_sym_posdef_log: cannot compute the logarithm, matrix not positive definite [%d, % 22.15g].", i, e_val);
    
    cblas_dscal (n, log (e_val), ncm_matrix_ptr (eve, i, 0), 1);
  }
  
  cblas_dsyr2k (CblasRowMajor, (UL == 'U') ? CblasUpper : CblasLower,
                CblasTrans, n, n,
                0.5, ncm_matrix_data (temp), ncm_matrix_tda (temp),
                ncm_matrix_data (eve), ncm_matrix_tda (eve),
                0.0, ncm_matrix_data (ln_cm), ncm_matrix_tda (ln_cm));
  
  g_array_unref (isuppz);
  ncm_lapack_ws_free (ws);
  ncm_vector_free (eva);
  ncm_matrix_free (eve);
  ncm_matrix_free (temp);
}

/**
 * ncm_matrix_triang_to_sym:
 * @cm: $M$ a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix
 * @zero: whether it should first set to zero the other side of the matrix
 * @sym: a #NcmMatrix to store the result
 *
 * Assuming that @cm is a triangular square matrix with data on @UL
 * side, computes the symmetric matrix $M^\intercal \times M$ if
 * @cm is upper triangular or $M\times M^\intercal$ if it is
 * lower triangular.
 *
 * If @zero is TRUE it first sets to zero all elements above/below
 * the diagonal for UL == 'L'/'U'. It should be TRUE whenever @cm
 * has non-zero values at the other side.
 *
 */
void
ncm_matrix_triang_to_sym (NcmMatrix *cm, gchar UL, gboolean zero, NcmMatrix *sym)
{
  const guint n = ncm_matrix_ncols (cm);
  gint i;
  
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (sym), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (sym), ==, ncm_matrix_nrows (sym));
  
  if (zero)
  {
    if (UL == 'U')
    {
      for (i = 0; i < n; i++)
      {
        gint j;
        
        for (j = i + 1; j < n; j++)
        {
          ncm_matrix_set (cm, j, i, 0.0);
        }
      }
    }
    else if (UL == 'L')
    {
      for (i = 0; i < n; i++)
      {
        gint j;
        
        for (j = i + 1; j < n; j++)
        {
          ncm_matrix_set (cm, i, j, 0.0);
        }
      }
    }
    else
    {
      g_assert_not_reached ();
    }
  }
  
  ncm_matrix_memcpy (sym, cm);
  
  if (UL == 'U')
    cblas_dtrmm (CblasRowMajor,
                 CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n, n,
                 1.0, ncm_matrix_data (cm), ncm_matrix_tda (cm),
                 ncm_matrix_data (sym), ncm_matrix_tda (sym));
  else if (UL == 'L')
    cblas_dtrmm (CblasRowMajor,
                 CblasRight, CblasLower, CblasTrans, CblasNonUnit, n, n,
                 1.0, ncm_matrix_data (cm), ncm_matrix_tda (cm),
                 ncm_matrix_data (sym), ncm_matrix_tda (sym));
  else
    g_assert_not_reached ();
}

/**
 * ncm_matrix_square_to_sym:
 * @cm: $M$ a #NcmMatrix
 * @NT: char indicating 'N' or 'T'
 * @UL: char indicating 'U'pper or 'L'ower matrix
 * @sym: a #NcmMatrix to store the result
 *
 * Computes the symmetric matrix $M^\intercal \times M$ if
 * @NT == 'T' or $M\times M^\intercal$ if @NT == 'N'. The result
 * is stored in the upper/lower triangle if @UL='U'/'L'
 *
 *
 */
void
ncm_matrix_square_to_sym (NcmMatrix *cm, gchar NT, gchar UL, NcmMatrix *sym)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  const CBLAS_UPLO Uplo       = (UL == 'U') ? CblasUpper : CblasLower;
  CBLAS_TRANSPOSE Trans;
  gint n, k;

  if (NT == 'N')
  {
    Trans = CblasNoTrans;
    n     = nrows;
    k     = ncols;
  }
  else
  {
    Trans = CblasTrans;
    n     = ncols;
    k     = nrows;
  }

  g_assert_cmpuint (ncm_matrix_ncols (sym), ==, n);
  g_assert_cmpuint (ncm_matrix_ncols (sym), ==, ncm_matrix_nrows (sym));

  cblas_dsyrk (CblasRowMajor, Uplo, Trans, n, k,
               1.0, ncm_matrix_data (cm),  ncm_matrix_tda (cm),
               0.0, ncm_matrix_data (sym), ncm_matrix_tda (sym));
}

/**
 * ncm_matrix_update_vector:
 * @cm: $M$ a #NcmMatrix
 * @NT: char indicating 'N' or 'T'
 * @alpha: a double $\alpha$
 * @v: a #NcmVector to update
 * @beta: a double $\beta$
 * @u: a #NcmVector to store the result
 *
 * Computes the matrix - vector product $u = \alpha M v + \beta u$
 * if @NT == 'N' or $u = \alpha M^\intercal v + u$ if @NT == 'T'
 * and stores the result in @u.
 */
void
ncm_matrix_update_vector (NcmMatrix *cm, gchar NT, const gdouble alpha, NcmVector *v, const gdouble beta, NcmVector *u)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  CBLAS_TRANSPOSE Trans;

  if (NT == 'N')
  {
    Trans = CblasNoTrans;
    g_assert_cmpuint (nrows, ==, ncm_vector_len (u));
    g_assert_cmpuint (ncols, ==, ncm_vector_len (v));
  }
  else
  {
    Trans = CblasTrans;
    g_assert_cmpuint (nrows, ==, ncm_vector_len (v));
    g_assert_cmpuint (ncols, ==, ncm_vector_len (u));
  }

  cblas_dgemv (CblasRowMajor, Trans, nrows, ncols,
      alpha, ncm_matrix_data (cm), ncm_matrix_tda (cm),
      ncm_vector_data (v), ncm_vector_stride (v),
      beta, ncm_vector_data (u), ncm_vector_stride (u));
}

/**
 * ncm_matrix_sym_update_vector:
 * @cm: $M$ a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix
 * @alpha: a double $\alpha$
 * @v: a #NcmVector to update
 * @beta: a double $\beta$
 * @u: a #NcmVector to store the result
 *
 * Computes the matrix - vector product $u = \alpha M v + \beta u$
 * if @NT == 'N' or $u = M^\intercal v$ if @NT == 'T'
 * and stores the result in @u. This function assumes
 * that $M$ is symmetric and it´s stored in the Upper/Lower
 * triangle if @UL == 'U'/'L'.
 */
void
ncm_matrix_sym_update_vector (NcmMatrix *cm, gchar UL, const gdouble alpha, NcmVector *v, const gdouble beta, NcmVector *u)
{
  const guint nrows     = ncm_matrix_nrows (cm);
  const guint ncols     = ncm_matrix_ncols (cm);
  const CBLAS_UPLO Uplo = (UL == 'U') ? CblasUpper : CblasLower;

  g_assert_cmpuint (nrows, ==, ncols);

  cblas_dsymv (CblasRowMajor, Uplo, nrows,
               alpha, ncm_matrix_data (cm), ncm_matrix_tda (cm),
               ncm_vector_data (v), ncm_vector_stride (v),
               beta, ncm_vector_data (u), ncm_vector_stride (u));
}

/**
 * ncm_matrix_log_vals:
 * @cm: a #NcmMatrix
 * @prefix: the prefixed text
 * @format: double format
 *
 * Prints to the log the values of @cm.
 *
 */
void
ncm_matrix_log_vals (NcmMatrix *cm, gchar *prefix, gchar *format)
{
  guint i, j;
  
  for (i = 0; i < ncm_matrix_nrows (cm); i++)
  {
    g_message ("%s", prefix);
    
    for (j = 0; j < ncm_matrix_ncols (cm); j++)
    {
      g_message (" ");
      g_message (format, ncm_matrix_get (cm, i, j));
    }
    
    g_message ("\n");
  }
}

/**
 * ncm_matrix_fill_rand_cor:
 * @cm: a square #NcmMatrix
 * @cor_level: correlation level parameter
 * @rng: a #NcmRNG
 *
 * Overwrite @cm with a random correlation matrix, the
 * parameter @cor_level controls the correlation between
 * entries the lower @cor_level more correlated the entries
 * are.
 *
 */
void
ncm_matrix_fill_rand_cor (NcmMatrix *cm, const gdouble cor_level, NcmRNG *rng)
{
  const guint n   = ncm_matrix_nrows (cm);
  const guint nm1 = n - 1;
  
  g_assert_cmpfloat (cor_level, >, 0.0);
  g_assert_cmpuint (n, ==, ncm_matrix_ncols (cm));
  
  ncm_rng_lock (rng);
  {
    NcmMatrix *P = ncm_matrix_dup (cm);
    gint k;
    
    ncm_matrix_set_all (P, 0.0);
    ncm_matrix_set_identity (cm);
    
    for (k = 0; k < nm1; k++)
    {
      gint i;
      
      for (i = k + 1; i < n; i++)
      {
        gdouble p = (gsl_ran_beta (rng->r, cor_level, cor_level) - 0.5) * 2.0;
        gint l;
        
        ncm_matrix_set (P, k, i, p);
        
        for (l = k - 1; l >= 0; l--)
        {
          const gdouble Pli = ncm_matrix_get (P, l, i);
          const gdouble Plk = ncm_matrix_get (P, l, k);
          
          p = p * sqrt ((1.0 - gsl_pow_2 (Pli)) * (1.0 - gsl_pow_2 (Plk))) + Pli * Plk;
        }
        
        ncm_matrix_set (cm, k, i, p);
        ncm_matrix_set (cm, i, k, p);
      }
    }
    
    ncm_matrix_free (P);
  }
  
  ncm_rng_unlock (rng);
}

/**
 * ncm_matrix_fill_rand_cov:
 * @cm: a square #NcmMatrix
 * @sigma_min: mininum standard deviation
 * @sigma_max: maximum standard deviation
 * @cor_level: correlation level parameter
 * @rng: a #NcmRNG
 *
 * Overwrite @cm with a random covariance matrix, the
 * parameter @cor_level controls the correlation between
 * entries the lower @cor_level more correlated the entries
 * are.
 *
 */
void
ncm_matrix_fill_rand_cov (NcmMatrix *cm, const gdouble sigma_min, const gdouble sigma_max, const gdouble cor_level, NcmRNG *rng)
{
  const guint n = ncm_matrix_nrows (cm);
  gint k;
  
  g_assert_cmpfloat (sigma_min, >, 0.0);
  g_assert_cmpfloat (sigma_max, >, sigma_min);
  
  ncm_matrix_fill_rand_cor (cm, cor_level, rng);
  
  ncm_rng_lock (rng);
  
  for (k = 0; k < n; k++)
  {
    const gdouble sigma_k = ncm_rng_uniform_gen (rng, sigma_min, sigma_max);
    
    ncm_matrix_mul_col (cm, k, sigma_k);
    ncm_matrix_mul_row (cm, k, sigma_k);
  }
  
  ncm_rng_unlock (rng);
}

/**
 * ncm_matrix_fill_rand_cov2:
 * @cm: a square #NcmMatrix
 * @mu: mean #NcmVector
 * @reltol_min: mininum standard deviation
 * @reltol_max: maximum standard deviation
 * @cor_level: correlation level parameter
 * @rng: a #NcmRNG
 *
 * Overwrite @cm with a random covariance matrix, the
 * parameter @cor_level controls the correlation between
 * entries the lower @cor_level more correlated the entries
 * are.
 *
 */
void
ncm_matrix_fill_rand_cov2 (NcmMatrix *cm, NcmVector *mu, const gdouble reltol_min, const gdouble reltol_max, const gdouble cor_level, NcmRNG *rng)
{
  const guint n = ncm_matrix_nrows (cm);
  gint k;
  
  g_assert_cmpfloat (reltol_min, >, 0.0);
  g_assert_cmpfloat (reltol_max, >, reltol_min);
  g_assert_cmpuint (n, ==, ncm_vector_len (mu));
  
  ncm_matrix_fill_rand_cor (cm, cor_level, rng);
  
  ncm_rng_lock (rng);
  
  for (k = 0; k < n; k++)
  {
    const gdouble mu_k    = ncm_vector_get (mu, k);
    const gdouble r_k     = exp (ncm_rng_uniform_gen (rng, log (reltol_min), log (reltol_max)));
    const gdouble sigma_k = mu_k != 0.0 ? fabs (mu_k) * r_k : r_k;
    
    ncm_matrix_mul_col (cm, k, sigma_k);
    ncm_matrix_mul_row (cm, k, sigma_k);
  }
  
  ncm_rng_unlock (rng);
}

/**
 * ncm_matrix_cov2cor:
 * @cov: a square #NcmMatrix
 * @cor: the output matrix
 *
 * Convert a covariance matrix @cov to a correlation
 * matrix @cor. The matrices @cor and @cov can be the same
 * object.
 *
 */
void
ncm_matrix_cov2cor (const NcmMatrix *cov, NcmMatrix *cor)
{
  const guint n = ncm_matrix_nrows (cov);
  gint i;
  
  g_assert_cmpuint (ncm_matrix_ncols (cov), ==, n);
  
  if (cov != cor)
    ncm_matrix_memcpy (cor, cov);
  
  for (i = 0; i < n; i++)
  {
    NcmVector *row_i  = ncm_matrix_get_row (cor, i);
    NcmVector *col_i  = ncm_matrix_get_col (cor, i);
    const gdouble w_i = 1.0 / sqrt (fabs (ncm_matrix_get (cov, i, i)));
    
    ncm_vector_scale (row_i, w_i);
    ncm_vector_scale (col_i, w_i);

    ncm_vector_free (row_i);
    ncm_vector_free (col_i);
  }
}

/**
 * ncm_matrix_cov_dup_cor:
 * @cov: a square #NcmMatrix
 *
 * Convert a covariance matrix @cov to a newly allocated
 * correlation matrix.
 *
 * Returns: (transfer full): the newly allocated correlation matrix.
 */
NcmMatrix *
ncm_matrix_cov_dup_cor (const NcmMatrix *cov)
{
  NcmMatrix *cor = ncm_matrix_dup (cov);
  
  ncm_matrix_cov2cor (cov, cor);
  
  return cor;
}

/**
 * ncm_matrix_new_gsl_const: (skip)
 * @gm: matrix from [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
 *
 * This function converts @gm into a constant #NcmMatrix.
 *
 * Returns: A new constant #NcmMatrix.
 */

/**
 * ncm_matrix_get:
 * @cm: a constant #NcmMatrix
 * @i: row index
 * @j: column index
 *
 *
 * Returns: The (@i,@j)-th element of the matrix @cm.
 */

/**
 * ncm_matrix_get_colmajor:
 * @cm: a constant #NcmMatrix
 * @i: row index
 * @j: column index
 *
 * Gets the (@i,@j)-th component of @cm assuming
 * a [column-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
 *
 * All column-major methods should be used carefully, they are inconsistent with
 * most other methods and are used mainly to interface with Fortran sub-routines.
 *
 * Returns: The (@i,@j)-th element of the matrix @cm.
 */

/**
 * ncm_matrix_ptr:
 * @cm: a #NcmMatrix
 * @i: row index
 * @j: column index
 *
 * Returns: A pointer to the (@i,@j)-th element of the matrix @cm.
 */

/**
 * ncm_matrix_const_ptr:
 * @cm: a #NcmMatrix
 * @i: row index
 * @j: column index
 *
 * Returns: A constant pointer to the (@i,@j)-th element of the matrix @cm.
 */

/**
 * ncm_matrix_set:
 * @cm: a #NcmMatrix
 * @i: row index
 * @j: column index
 * @val: a double
 *
 * This function sets the value of the (@i,@j)-th element of the matrix @cm to @val.
 *
 */

/**
 * ncm_matrix_set_colmajor:
 * @cm: a #NcmMatrix
 * @i: row index
 * @j: column index
 * @val: a double
 *
 * This function sets the value of the (@i,@j)-th element of the matrix @cm to @val
 * considering it being in the [column-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
 *
 * All column-major methods should be used carefully, they are inconsistent with
 * most other methods and are used mainly to interface with Fortran sub-routines.
 *
 */

/**
 * ncm_matrix_addto:
 * @cm: a #NcmMatrix
 * @i: row index
 * @j: column index
 * @val: a double
 *
 * This function adds the value @val to the (@i,@j)-th element of the matrix @cm.
 *
 */

/**
 * ncm_matrix_transpose:
 * @cm: a #NcmMatrix
 *
 * This function replaces the matrix @cm by its transpose by copying the elements of the matrix in-place.
 * The matrix must be square for this operation to be possible.
 *
 */

/**
 * ncm_matrix_set_identity:
 * @cm: a #NcmMatrix
 *
 * This function sets the elements of the matrix @cm to the corresponding elements of the identity matrix,
 * i.e. a unit diagonal with all off-diagonal elements zero. This applies to both square and rectangular matrices.
 *
 */

/**
 * ncm_matrix_set_zero:
 * @cm: a #NcmMatrix
 *
 * This function sets all the elements of the matrix @cm to zero.
 *
 */

/**
 * ncm_matrix_set_all:
 * @cm: a #NcmMatrix
 * @val: a double
 *
 * This function sets all the elements of the matrix @cm to @val.
 *
 */

/**
 * ncm_matrix_add:
 * @cm1: a #NcmMatrix
 * @cm2: a #NcmMatrix
 *
 * This function adds the elements of the matrices @cm1 and @cm2.
 * The two matrices must have the same size.
 *
 */

/**
 * ncm_matrix_sub:
 * @cm1: a #NcmMatrix
 * @cm2: a #NcmMatrix
 *
 * This function subtracts the elements of the matrices @cm1 and @cm2.
 * The two matrices must have the same size.
 *
 */

/**
 * ncm_matrix_mul_elements:
 * @cm1: a #NcmMatrix
 * @cm2: a #NcmMatrix
 *
 * This function multiplies the elements of the matrices @cm1 and @cm2.
 * The two matrices must have the same size.
 *
 */

/**
 * ncm_matrix_div_elements:
 * @cm1: a #NcmMatrix
 * @cm2: a #NcmMatrix
 *
 * This function divides the elements of the matrices @cm1 and @cm2.
 * The two matrices must have the same size.
 *
 */

/**
 * ncm_matrix_scale:
 * @cm: a #NcmMatrix
 * @val: a double
 *
 * This function multiplies the elements of the matrix @cm by the constant factor @val.
 * The result is stored in @cm.
 *
 */

/**
 * ncm_matrix_add_constant:
 * @cm: a #NcmMatrix
 * @val: a double
 *
 * This function adds the the constant factor @val to the elements of the matrix @cm.
 * The result is stored in @cm.
 *
 */

/**
 * ncm_matrix_mul_row:
 * @cm: a #NcmMatrix
 * @row_i: row index
 * @val: a double
 *
 * This function multiplies row @row_i elements by @val.
 *
 */

/**
 * ncm_matrix_mul_col:
 * @cm: a #NcmMatrix
 * @col_i: column index
 * @val: a double
 *
 * This function multiplies column @col_i elements by @val.
 *
 */

/**
 * ncm_matrix_get_diag:
 * @cm: a #NcmMatrix
 * @diag: a #NcmVector
 *
 * This function copies de diagonal elements of the matrix @cm to the vector @diag.
 *
 */

/**
 * ncm_matrix_set_diag:
 * @cm: a #NcmMatrix
 * @diag: a #NcmVector
 *
 * This function copies de the elements of the vector @diag to the diagonal elements of the matrix @cm.
 *
 */

/**
 * ncm_matrix_memcpy:
 * @cm1: a #NcmMatrix
 * @cm2: a #NcmMatrix
 *
 * This function copies the elements of the matrix @cm2 into the matrix @cm1.
 * The two matrices must have the same size.
 *
 */

/**
 * ncm_matrix_memcpy_to_colmajor:
 * @cm1: a #NcmMatrix
 * @cm2: a #NcmMatrix
 *
 * This function copies the elements of the matrix @cm2 into the matrix @cm1.
 * The two matrices must have the same size. The elements are written in @cm1
 * in [column-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order)
 * order.
 *
 * All column-major methods should be used carefully, they are inconsistent with
 * most other methods and are used mainly to interface with Fortran sub-routines.
 */

/**
 * ncm_matrix_set_col:
 * @cm: a #NcmMatrix
 * @n: column index
 * @cv: a constant #NcmVector
 *
 * This function copies the elements of the vector @cv into the @n-th column of the matrix @cm.
 * The length of the vector must be the same as the length of the column.
 *
 */

/**
 * ncm_matrix_set_row:
 * @cm: a #NcmMatrix
 * @n: row index
 * @cv: a constant #NcmVector
 *
 * This function copies the elements of the vector @cv into the @n-th row of the matrix @cm.
 * The length of the vector must be the same as the length of the row.
 *
 */

/**
 * ncm_matrix_get_array:
 * @cm: a #NcmMatrix
 *
 * This function returns the array of @cv. It is only applied if the matrix @cm was created with ncm_matrix_new_array().
 *
 * Returns: (transfer container) (element-type double): A pointer to a double GArray.
 */

/**
 * ncm_matrix_fast_get:
 * @cm: a #NcmMatrix
 * @ij: element index of the #NcmMatrix base data
 *
 * This function returns the value of the @cm[@i,@j] element by direct access of its base data. Where @ij = i $\times$ @tda $+$ j.
 *
 * If the matrix was created with ncm_matrix_new() or ncm_matrix_new0() then @tda = @ncols.
 *
 * Returns: The (@i ,@j)-th element of @cm.
 */

/**
 * ncm_matrix_fast_set:
 * @cm: a #NcmMatrix
 * @ij: element index of the #NcmMatrix base data
 * @val: a double
 *
 * This function sets the value of the @cm[@i,@j] element to @val by direct access of its base data. Where @ij = i $\times$ @tda $+$ j.
 *
 * If the matrix was created with ncm_matrix_new() or ncm_matrix_new0() then @tda = @ncols.
 *
 */

/**
 * ncm_matrix_gsl: (skip)
 * @cm: a #NcmMatrix
 *
 * This function returns a pointer to the #gsl_matrix associated to the matrix @cm.
 *
 * Returns: A pointer to a #gsl_matrix.
 */

/**
 * ncm_matrix_const_gsl: (skip)
 * @cm: a #NcmMatrix
 *
 * This function returns a constant pointer to the #gsl_matrix associated to the matrix @cm.
 *
 * Returns: A constant pointer to a #gsl_matrix.
 */

/**
 * ncm_matrix_col_len:
 * @cm: a #NcmMatrix
 *
 * This function returns the number of elements in a column of the matrix @cm. The columns length.
 *
 * Returns: The columns length of @cm (a.k.a. @nrows).
 */

/**
 * ncm_matrix_row_len:
 * @cm: a #NcmMatrix
 *
 * This function returns the number of elements in a row of the matrix @cm. The rows length.
 *
 * Returns: The rows length of @cm (a.k.a. @ncols).
 */

/**
 * ncm_matrix_nrows:
 * @cm: a #NcmMatrix
 *
 * This function returns the number of elements in a row of the matrix @cm.
 *
 * Returns: The number of elements in a row of @cm.
 */

/**
 * ncm_matrix_ncols:
 * @cm: a #NcmMatrix
 *
 * This function returns the number of elements in a column of the matrix @cm.
 *
 * Returns: The number of elements in a column of @cm.
 */

/**
 * ncm_matrix_size:
 * @cm: a #NcmMatrix
 *
 * Calculates the total size of the matrix, @ncols $\times$ @nrows.
 *
 * Returns: Total size of the matrix.
 */

/**
 * ncm_matrix_tda:
 * @cm: a #NcmMatrix
 *
 * This functions returns the matrix @cm @tda value.
 *
 * Returns: The matrix tda.
 */

/**
 * ncm_matrix_data:
 * @cm: a #NcmMatrix
 *
 * This function returns a pointer to the matrix @cm base data.
 *
 * Returns: (transfer none): A pointer to @cm base data.
 */

/**
 * ncm_matrix_const_data:
 * @cm: a #NcmMatrix
 *
 * This function returns a constant pointer to the matrix @cm base data.
 *
 * Returns: (transfer none): A constant pointer to the matrix @cm base data.
 */

