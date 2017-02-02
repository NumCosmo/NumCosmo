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
 * @title: NcmMatrix
 * @short_description: Matrix object representing an array of doubles.
 *
 * This object defines the functions for allocating and accessing matrices.
 * Also includes serveral matrix operations.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_matrix.h"
#include "math/ncm_lapack.h"

#if (defined HAVE_CLAPACK_H) && (defined HAVE_CLAPACK_DPOTRF)
#include <clapack.h>
#else
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#endif

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
  m->type = 0;
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
      g_slice_free1 (sizeof(gdouble) * ncm_matrix_nrows (cm) * ncm_matrix_ncols (cm), ncm_matrix_data (cm));
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
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_matrix_set_property;
  object_class->get_property = &_ncm_matrix_get_property;
  object_class->dispose      = &_ncm_matrix_dispose;
  object_class->finalize     = &_ncm_matrix_finalize;

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
  gdouble *d = g_slice_alloc (sizeof(gdouble) * nrows * ncols);
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
  gdouble *d = g_slice_alloc0 (sizeof(gdouble) * nrows * ncols);
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
 * @gm: matrix from GNU Scientific Library (GSL) to be converted into a #NcmMatrix
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
                                       gm, (GDestroyNotify) &gsl_matrix_free);
  cm->type = NCM_MATRIX_GSL_MATRIX;
  return cm;
}

/**
 * ncm_matrix_new_gsl_static: (skip)
 * @gm: matrix from GNU Scientific Library (GSL) to be converted into a #NcmMatrix
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
                                         (GDestroyNotify) &g_array_unref);
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
 * Returns: (transfer full) :a #NcmMatrix with the values from @var.
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
  cm->mv = gsl_matrix_view_array ((gdouble *)d, nrows, ncols); 
  cm->type = NCM_MATRIX_DERIVED;
  return cm;
}

/**
 * ncm_matrix_const_new_variant:
 * @var: a variant of type "aad"
 * 
 * Creates a new constant matrix using the same memory of @var.
 * 
 * Returns: (transfer full) :a #NcmMatrix with the values from @var.
 */
const NcmMatrix *
ncm_matrix_const_new_variant (GVariant *var)
{
  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE ("aad")));
  {
    GVariant *row = g_variant_get_child_value (var, 0);
    guint nrows = g_variant_n_children (var);
    guint ncols = g_variant_n_children (row);
    gconstpointer data = g_variant_get_data (var);
    const NcmMatrix *m = ncm_matrix_const_new_data (data, nrows, ncols);

    NCM_MATRIX (m)->pdata = g_variant_ref_sink (var);
    NCM_MATRIX (m)->pfree = (GDestroyNotify) &g_variant_unref;

    return m;
  }
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
  scm->mv = gsl_matrix_submatrix (ncm_matrix_gsl (cm), k1, k2, nrows, ncols);

  scm->pdata = g_object_ref (cm);
  scm->pfree = g_object_unref;
  scm->type = NCM_MATRIX_DERIVED;

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
  cv->vv    = gsl_matrix_column (ncm_matrix_gsl (cm), col);
  cv->type  = NCM_VECTOR_DERIVED;

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
  cv->vv = gsl_matrix_row (ncm_matrix_gsl (cm), row);
  cv->type = NCM_VECTOR_DERIVED;

  cv->pdata = g_object_ref (cm);
  cv->pfree = &g_object_unref;

  return cv;
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
    guint nrows = g_variant_n_children (var);
    guint ncols = g_variant_n_children (row);
    guint i;
    g_variant_unref (row);

    /* Sometimes we receive a NcmMatrix in the process of instantiation. */
    if ((ncm_matrix_nrows (cm) == 0) && (ncm_matrix_ncols (cm) == 0))
    {
      gdouble *d = g_slice_alloc (sizeof (gdouble) * nrows * ncols);
      cm->mv = gsl_matrix_view_array (d, nrows, ncols);
      cm->type = NCM_MATRIX_SLICE;
    }
    else if (nrows != ncm_matrix_nrows (cm) || ncols != ncm_matrix_ncols (cm))
      g_error ("ncm_matrix_set_from_variant: cannot set matrix values, variant contains (%u, %u) childs but matrix dimension is (%u, %u)", nrows, ncols, ncm_matrix_nrows (cm), ncm_matrix_ncols (cm));

    for (i = 0; i < nrows; i++)
    {
      NcmVector *m_row = ncm_matrix_get_row (cm, i);
      row = g_variant_get_child_value (var, i);
      
      {
        gsize v_ncols = 0;
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
  const guint nrows     = ncm_matrix_nrows (cm);
  const guint ncols     = ncm_matrix_ncols (cm);
  GVariant **rows       = g_new (GVariant *, nrows);
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
  guint row_size = ncols * sizeof (gdouble);
  guint i = 0;

  rows[i] = g_variant_new_from_data (G_VARIANT_TYPE ("ad"),
                                     ncm_matrix_ptr (cm, i, 0),
                                     row_size,
                                     TRUE,
                                     (GDestroyNotify) &ncm_matrix_free,
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
    if (nm != NULL && check_size)
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
 * @cm: FIXME
 * @alpha: FIXME
 * @b: FIXME
 *
 * FIXME
 *
 */
void 
ncm_matrix_add_mul (NcmMatrix *cm, const gdouble alpha, NcmMatrix *b)
{
  const gboolean no_pad_cm = ncm_matrix_gsl (cm)->tda == ncm_matrix_ncols (cm);
  const gboolean no_pad_b = ncm_matrix_gsl (b)->tda == ncm_matrix_ncols (b);

  g_assert (ncm_matrix_ncols (cm) == ncm_matrix_ncols (b));
  g_assert (ncm_matrix_nrows (cm) == ncm_matrix_nrows (b));

  if (no_pad_cm && no_pad_b)
  {
    const gint N = ncm_matrix_ncols (b) * ncm_matrix_nrows (b);
    cblas_daxpy (N, alpha, 
                 ncm_matrix_data (b), 1,
                 ncm_matrix_data (cm), 1);
  }
  else
  {
    const gint N = ncm_matrix_row_len (b);
    const guint b_tda = ncm_matrix_gsl (b)->tda;
    const guint cm_tda = ncm_matrix_gsl (cm)->tda;
    guint i;

    for (i = 0; i < ncm_matrix_nrows (b); i++)
    {
      cblas_daxpy (N, alpha, 
                   &ncm_matrix_data (b)[b_tda * i], 1,
                   &ncm_matrix_data (cm)[cm_tda * i], 1);
    }
  }
}

/**
 * ncm_matrix_copy_triangle:
 * @cm: a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix 
 *
 * If @UL == 'U' copy the upper triangle over the lower.
 * If @UL == 'L' copy the lower triangle over the lower.
 *
 */
void 
ncm_matrix_copy_triangle (NcmMatrix *cm, gchar UL)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_nrows (cm);
  guint i, j;

  if (nrows != ncols)
    g_error ("ncm_matrix_copy_triangle: only works on square a matrix [%ux%u]", nrows, ncols);
  if (UL != 'U' && UL != 'L')
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
  else
  {
    for (i = 0; i < nrows; i++)
    {
      for (j = i + 1; j < ncols; j++)
      {
        ncm_matrix_set (cm, i, j, ncm_matrix_get (cm, j, i));
      }
    }
  }
}

/**
 * ncm_matrix_dsymm:
 * @cm: a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix 
 * @alpha: FIXME
 * @b: FIXME
 * @beta: FIXME
 * @c: FIXME
 *
 * FIXME
 *
 */
void 
ncm_matrix_dsymm (NcmMatrix *cm, gchar UL, const gdouble alpha, NcmMatrix *b, const gdouble beta, NcmMatrix *c)
{
  g_assert (UL == 'U' || UL == 'L');
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_matrix_ncols (b));
  g_assert_cmpuint (ncm_matrix_nrows (cm), ==, ncm_matrix_nrows (b));
  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_matrix_ncols (c));
  g_assert_cmpuint (ncm_matrix_nrows (cm), ==, ncm_matrix_nrows (c));

  cblas_dsymm (CblasRowMajor, CblasLeft, (UL == 'U') ? CblasUpper : CblasLower, ncm_matrix_nrows (cm), ncm_matrix_ncols (cm), 
               alpha, 
               ncm_matrix_data (cm), ncm_matrix_gsl (cm)->tda,
               ncm_matrix_data (b), ncm_matrix_gsl (b)->tda, 
               beta,
               ncm_matrix_data (c), ncm_matrix_gsl (c)->tda);
}

/**
 * ncm_matrix_cholesky_decomp:
 * @cm: a #NcmMatrix
 * @UL: char indicating 'U'pper or 'L'ower matrix
 *
 * Calculates inplace the Cholesky decomposition for a symmetric positive
 * definite matrix.
 * 
 */
gint
ncm_matrix_cholesky_decomp (NcmMatrix *cm, gchar UL)
{
  gint ret = ncm_lapack_dpotrf (UL, ncm_matrix_nrows (cm), ncm_matrix_data (cm), ncm_matrix_nrows (cm));
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
  gint ret = ncm_lapack_dpotri (UL, ncm_matrix_nrows (cm), ncm_matrix_data (cm), ncm_matrix_nrows (cm));
  return ret;
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
 * ncm_matrix_new_gsl_const: (skip)
 * @m: matrix from GNU Scientific Library (GSL)
 *
 * This function converts @m into a constant #NcmMatrix.
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
 * ncm_matrix_ref:
 * @cm: a #NcmMatrix
 *
 * Increase the reference count of @cm by one.
 *
 * Returns: (transfer full): @cm
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
 * ncm_matrix_scale:
 * @cm: a #NcmMatrix
 * @val: a double
 *
 * This function multiplies the elements of the matrix @cm by the constant factor @val.
 * The result is stored in @cm.
 *
 */
/**
 * ncm_matrix_memcpy:
 * @cm1: a #NcmMatrix
 * @cm2: a #NcmMatrix
 *
 * This function copies the elements of the matrix @cm1 into the matrix @cm2.
 * The two matrices must have the same size.
 *
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
 * ncm_matrix_get_array:
 * @cm: a #NcmMatrix
 *
 * FIXME
 *
 * Returns: (transfer container) (element-type double): FIXME
 */
/**
 * ncm_matrix_fast_get:
 * @cm: a #NcmMatrix
 * @ij: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_matrix_fast_set:
 * @cm: a #NcmMatrix
 * @ij: FIXME
 * @val: FIXME
 *
 * FIXME
 *
 */
/**
 * ncm_matrix_gsl: (skip)
 * @cm: a #NcmMatrix
 *
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_matrix_const_gsl: (skip)
 * @cm: a #NcmMatrix
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_matrix_col_len:
 * @cm: a #NcmMatrix
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_matrix_row_len:
 * @cm: a #NcmMatrix
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_matrix_nrows:
 * @cm: a #NcmMatrix
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_matrix_ncols:
 * @cm: a #NcmMatrix
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_matrix_size:
 * @cm: a #NcmMatrix
 *
 * Calculates the total size of the matrix, ncols * nrows.
 *
 * Returns: Total size of the matrix.
 */
/**
 * ncm_matrix_tda:
 * @cm: a #NcmMatrix
 *
 *
 * Returns: The matrix tda.
 */
/**
 * ncm_matrix_data:
 * @cm: a #NcmMatrix
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
/**
 * ncm_matrix_const_data:
 * @cm: a #NcmMatrix
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */

