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
 * NcmMatrix:
 *
 * Reference-counted row-major matrix of doubles, a view over a #gsl_matrix.
 *
 * Rows are separated by the trailing dimension, the tda, which is at least the number of
 * columns. The data can be allocated by the matrix or come from GSL, a #GArray, a #GVariant,
 * another #NcmMatrix (a submatrix) or a user array, with the ownership stated by each
 * constructor. Functions whose names contain `colmajor` read or write the data in
 * column-major order, for Fortran routines; the rest of the API is row-major. Wrappers of
 * BLAS and LAPACK take their character arguments ('N'/'T', 'U'/'L', 'L'/'R') in the row-major
 * sense.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_cblas.h" /* This must be included before any gsl header */
#include "ncm/algebra/ncm_matrix.h"
#include "ncm/algebra/ncm_vector.h"
#include "ncm/algebra/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_linalg.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_VALS,
};

G_DEFINE_TYPE (NcmMatrix, ncm_matrix, G_TYPE_OBJECT)

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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
                                                         G_VARIANT_TYPE ("aad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_matrix_new:
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * Allocates a matrix, not initialized.
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * Allocates a matrix with every element zero.
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @d: the data
 * @nrows: number of rows
 * @ncols: number of columns
 * @tda: distance between consecutive rows, in doubles
 * @pdata: (allow-none): data owned by the matrix
 * @pfree: (scope notified) (allow-none): function releasing @pdata
 *
 * Creates a matrix over @d; @pfree is called on @pdata when the matrix is finalized.
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @gm: a #gsl_matrix
 *
 * Creates a matrix over @gm, which it takes ownership of and frees with gsl_matrix_free().
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @gm: a #gsl_matrix
 *
 * Creates a matrix over @gm, which must outlive it.
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @a: (array) (element-type double): a #GArray of doubles
 * @ncols: number of columns
 *
 * Creates a matrix over the elements of @a, holding a reference to it, with @ncols columns
 * and the length of @a divided by @ncols rows. Aborts if @ncols does not divide the length.
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @d: the data, allocated with g_slice_alloc()
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * Creates a matrix over @d, which it takes ownership of.
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @d: the data, allocated with g_malloc() or malloc()
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * Creates a matrix over @d, which it takes ownership of and frees with g_free().
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @d: the data
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * Creates a matrix over @d, which must outlive it.
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @d: the data
 * @nrows: number of rows
 * @ncols: number of columns
 * @tda: distance between consecutive rows, in doubles
 *
 * Creates a matrix over @d, which must outlive it.
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @var: a #GVariant of type `aad`
 *
 * Creates a matrix with a copy of the elements of @var.
 *
 * Returns: (transfer full): a new #NcmMatrix.
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
 * @d: the data
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * Creates a constant matrix over @d, which must outlive it.
 *
 * Returns: (transfer full): a new constant #NcmMatrix.
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
 * @var: a #GVariant of type `aad`
 *
 * Creates a constant matrix over the data of @var, holding a reference to it.
 *
 * Returns: (transfer full): a new constant #NcmMatrix.
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
    NcmMatrix *m       = (NcmMatrix *) ncm_matrix_const_new_data (data, nrows, ncols);

    NCM_MATRIX (m)->pdata = g_variant_ref_sink (var);
    NCM_MATRIX (m)->pfree = (GDestroyNotify) & g_variant_unref;

    return m;
  }
}

/**
 * ncm_matrix_ref:
 * @cm: a #NcmMatrix
 *
 * Increases the reference count of @cm by one.
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
 * @k1: row of the first element
 * @k2: column of the first element
 * @nrows: number of rows
 * @ncols: number of columns
 *
 * Creates a view of the @nrows by @ncols block of @cm whose first element is (@k1, @k2),
 * holding a reference to @cm.
 *
 * Returns: (transfer full): the submatrix.
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
 * Creates a view of the column @col, holding a reference to @cm.
 *
 * Returns: (transfer full): the column as a #NcmVector.
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
 * Creates a view of the row @row, holding a reference to @cm.
 *
 * Returns: (transfer full): the row as a #NcmVector.
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
 * Creates a view of all the elements of @cm, row after row, holding a reference to @cm.
 * The tda of @cm must equal its number of columns.
 *
 * Returns: (transfer full): the #NcmVector.
 */
NcmVector *
ncm_matrix_as_vector (NcmMatrix *cm)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  const guint len   = nrows * ncols;

  NcmVector *v = ncm_vector_new_full (ncm_matrix_data (cm), len, 1, ncm_matrix_ref (cm), (GDestroyNotify) & ncm_matrix_free);

  g_assert (ncm_matrix_tda (cm) == ncm_matrix_ncols (cm));

  return v;
}

/**
 * ncm_matrix_set_from_variant:
 * @cm: a #NcmMatrix
 * @var: a #GVariant of type `aad`
 *
 * Sets the elements of @cm to those of @var. A matrix without elements is allocated with
 * the shape of @var. Aborts if @var has another type or shape.
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
 * Returns: (transfer full): a #GVariant of type `aad` with a copy of the elements.
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
 * Creates a #GVariant of type `aad` sharing the data of @cm, which must not change while
 * the variant exists.
 *
 * Returns: (transfer full): the #GVariant.
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
 * @data: (array) (element-type double): the values, row after row
 *
 * Sets the elements of @cm; @data must hold as many values as @cm has elements.
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
 * @a: (array) (element-type double): the values, row after row
 *
 * Sets the elements of @cm; @a must hold as many values as @cm has elements.
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
 * Decreases the reference count of @cm by one.
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
 * If *@cm is not %NULL, decreases its reference count by one and sets *@cm to %NULL.
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
 * Decreases the reference count of @cm by one.
 */
void
ncm_matrix_const_free (const NcmMatrix *cm)
{
  ncm_matrix_free (NCM_MATRIX ((NcmMatrix *) cm));
}

/**
 * ncm_matrix_dup:
 * @cm: a constant #NcmMatrix
 *
 * Returns: (transfer full): a newly allocated copy of @cm, with tda equal to its number of
 * columns.
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
 * @check_size: whether to require the same shape
 *
 * Replaces *@cm by a new reference to @nm, releasing the previous one. With @check_size,
 * aborts if both are set and their shapes differ.
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
 * @a: a double $a$
 * @cm2: a #NcmMatrix
 *
 * Sets $M_1 \to M_1 + a M_2$.
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
 * @scale: the scale $s$
 *
 * Returns: $\max_{ij} |(M_{1,ij} - M_{2,ij}) / (s + M_{2,ij})|$.
 */
gdouble
ncm_matrix_cmp (const NcmMatrix *cm1, const NcmMatrix *cm2, const gdouble scale)
{
  const guint nrows = ncm_matrix_nrows (cm1);
  const guint ncols = ncm_matrix_ncols (cm1);
  gdouble reltol    = 0.0;
  guint i, j;

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
 * @cm1: a constant #NcmMatrix
 * @cm2: a constant #NcmMatrix
 * @scale: the scale $s$
 *
 * Same as ncm_matrix_cmp() over the diagonal only.
 *
 * Returns: $\max_i |(M_{1,ii} - M_{2,ii}) / (s + M_{2,ii})|$.
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
 * ncm_matrix_zero_triangle:
 * @cm: a #NcmMatrix
 * @UL: 'U' or 'L', the triangle kept
 *
 * Sets to zero the strict lower triangle of the square @cm for @UL 'U', or its strict upper
 * triangle for 'L'. Factorizations leave the other triangle as it was, so a routine reading
 * the whole matrix needs it cleared. Aborts if @cm is not square or @UL is invalid.
 */
void
ncm_matrix_zero_triangle (NcmMatrix *cm, gchar UL)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  guint i, j;

  if (nrows != ncols)
    g_error ("ncm_matrix_zero_triangle: only works on a square matrix [%ux%u]", nrows, ncols);

  if ((UL != 'U') && (UL != 'L'))
    g_error ("ncm_matrix_zero_triangle: expect U or L and received %c.", UL);

  for (i = 0; i < nrows; i++)
  {
    for (j = i + 1; j < ncols; j++)
    {
      if (UL == 'U')
        ncm_matrix_set (cm, j, i, 0.0);
      else
        ncm_matrix_set (cm, i, j, 0.0);
    }
  }
}

/**
 * ncm_matrix_copy_triangle:
 * @cm: a #NcmMatrix
 * @UL: 'U' or 'L', the triangle copied
 *
 * Copies the upper triangle of the square @cm over the lower for @UL 'U', or the lower over
 * the upper for 'L'. Aborts if @cm is not square or @UL is invalid.
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
 * @cm: the result $C$
 * @UL: 'U' or 'L', the triangle used
 * @alpha: $\alpha$
 * @A: a symmetric #NcmMatrix $A$
 * @B: a #NcmMatrix $B$
 * @beta: $\beta$
 *
 * Sets $C \to \alpha A B + \beta C$, reading only the @UL triangle of $A$. The three
 * matrices are square and of the same size.
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

static CBLAS_UPLO
_ncm_matrix_check_uplo (const gchar *func_name, gchar UL)
{
  switch (UL)
  {
    case 'U':

      return CblasUpper;

    case 'L':

      return CblasLower;

    default:
      g_error ("%s: expect U or L and received %c.", func_name, UL);

      return 0;
  }
}

static CBLAS_SIDE
_ncm_matrix_check_side (const gchar *func_name, gchar Side)
{
  switch (Side)
  {
    case 'L':

      return CblasLeft;

    case 'R':

      return CblasRight;

    default:
      g_error ("%s: expect L or R and received %c.", func_name, Side);

      return 0;
  }
}

/**
 * ncm_matrix_dgemm:
 * @cm: the result $C$
 * @TransA: 'N' or 'T', whether to transpose $A$
 * @TransB: 'N' or 'T', whether to transpose $B$
 * @alpha: $\alpha$
 * @A: a #NcmMatrix $A$
 * @B: a #NcmMatrix $B$
 * @beta: $\beta$
 *
 * Sets $C \to \alpha\,\mathrm{op}(A)\,\mathrm{op}(B) + \beta C$.
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
 * ncm_matrix_dtrmm:
 * @cm: a #NcmMatrix $B$
 * @Side: 'L' or 'R', the side $A$ acts from
 * @UL: 'U' or 'L', whether $A$ is upper or lower triangular
 * @TransA: 'N' or 'T', whether to transpose $A$
 * @alpha: $\alpha$
 * @A: a triangular #NcmMatrix $A$
 *
 * Sets $B \to \alpha\,\mathrm{op}(A)\,B$ for @Side 'L', or $B \to \alpha B\,\mathrm{op}(A)$ for 'R'.
 * Only the @UL triangle of $A$ is read, its diagonal taken as stored.
 */
void
ncm_matrix_dtrmm (NcmMatrix *cm, gchar Side, gchar UL, gchar TransA, const gdouble alpha, NcmMatrix *A)
{
  const CBLAS_SIDE cblas_Side        = _ncm_matrix_check_side ("ncm_matrix_dtrmm", Side);
  const CBLAS_UPLO cblas_UL          = _ncm_matrix_check_uplo ("ncm_matrix_dtrmm", UL);
  const CBLAS_TRANSPOSE cblas_TransA = _ncm_matrix_check_trans ("ncm_matrix_dtrmm", TransA);
  const guint nrows                  = ncm_matrix_nrows (cm);
  const guint ncols                  = ncm_matrix_ncols (cm);

  g_assert_cmpuint (ncm_matrix_nrows (A), ==, ncm_matrix_ncols (A));
  g_assert_cmpuint (ncm_matrix_nrows (A), ==, (cblas_Side == CblasLeft) ? nrows : ncols);

  cblas_dtrmm (CblasRowMajor, cblas_Side, cblas_UL, cblas_TransA, CblasNonUnit, nrows, ncols,
               alpha,
               ncm_matrix_data (A), ncm_matrix_tda (A),
               ncm_matrix_data (cm), ncm_matrix_tda (cm));
}

/**
 * ncm_matrix_dtrsm:
 * @cm: a #NcmMatrix $B$
 * @Side: 'L' or 'R', the side $A$ acts from
 * @UL: 'U' or 'L', whether $A$ is upper or lower triangular
 * @TransA: 'N' or 'T', whether to transpose $A$
 * @alpha: $\alpha$
 * @A: a triangular #NcmMatrix $A$
 *
 * Sets $B \to \alpha\,\mathrm{op}(A)^{-1} B$ for @Side 'L', or $B \to \alpha B\,\mathrm{op}(A)^{-1}$
 * for 'R'. Only the @UL triangle of $A$ is read, its diagonal taken as stored.
 */
void
ncm_matrix_dtrsm (NcmMatrix *cm, gchar Side, gchar UL, gchar TransA, const gdouble alpha, NcmMatrix *A)
{
  const CBLAS_SIDE cblas_Side        = _ncm_matrix_check_side ("ncm_matrix_dtrsm", Side);
  const CBLAS_UPLO cblas_UL          = _ncm_matrix_check_uplo ("ncm_matrix_dtrsm", UL);
  const CBLAS_TRANSPOSE cblas_TransA = _ncm_matrix_check_trans ("ncm_matrix_dtrsm", TransA);
  const guint nrows                  = ncm_matrix_nrows (cm);
  const guint ncols                  = ncm_matrix_ncols (cm);

  g_assert_cmpuint (ncm_matrix_nrows (A), ==, ncm_matrix_ncols (A));
  g_assert_cmpuint (ncm_matrix_nrows (A), ==, (cblas_Side == CblasLeft) ? nrows : ncols);

  cblas_dtrsm (CblasRowMajor, cblas_Side, cblas_UL, cblas_TransA, CblasNonUnit, nrows, ncols,
               alpha,
               ncm_matrix_data (A), ncm_matrix_tda (A),
               ncm_matrix_data (cm), ncm_matrix_tda (cm));
}

/**
 * ncm_matrix_dtrmv:
 * @cm: a triangular #NcmMatrix $A$
 * @UL: 'U' or 'L', whether $A$ is upper or lower triangular
 * @Trans: 'N' or 'T', whether to transpose $A$
 * @v: a #NcmVector $v$
 *
 * Sets $v \to \mathrm{op}(A)\,v$. Only the @UL triangle of $A$ is read, its diagonal taken as
 * stored.
 */
void
ncm_matrix_dtrmv (NcmMatrix *cm, gchar UL, gchar Trans, NcmVector *v)
{
  const CBLAS_UPLO cblas_UL         = _ncm_matrix_check_uplo ("ncm_matrix_dtrmv", UL);
  const CBLAS_TRANSPOSE cblas_Trans = _ncm_matrix_check_trans ("ncm_matrix_dtrmv", Trans);
  const guint n                     = ncm_matrix_nrows (cm);

  g_assert_cmpuint (n, ==, ncm_matrix_ncols (cm));
  g_assert_cmpuint (n, ==, ncm_vector_len (v));

  cblas_dtrmv (CblasRowMajor, cblas_UL, cblas_Trans, CblasNonUnit, n,
               ncm_matrix_data (cm), ncm_matrix_tda (cm),
               ncm_vector_data (v), ncm_vector_stride (v));
}

/**
 * ncm_matrix_dtrsv:
 * @cm: a triangular #NcmMatrix $A$
 * @UL: 'U' or 'L', whether $A$ is upper or lower triangular
 * @Trans: 'N' or 'T', whether to transpose $A$
 * @v: a #NcmVector $v$
 *
 * Sets $v \to \mathrm{op}(A)^{-1} v$. Only the @UL triangle of $A$ is read, its diagonal taken
 * as stored.
 */
void
ncm_matrix_dtrsv (NcmMatrix *cm, gchar UL, gchar Trans, NcmVector *v)
{
  const CBLAS_UPLO cblas_UL         = _ncm_matrix_check_uplo ("ncm_matrix_dtrsv", UL);
  const CBLAS_TRANSPOSE cblas_Trans = _ncm_matrix_check_trans ("ncm_matrix_dtrsv", Trans);
  const guint n                     = ncm_matrix_nrows (cm);

  g_assert_cmpuint (n, ==, ncm_matrix_ncols (cm));
  g_assert_cmpuint (n, ==, ncm_vector_len (v));

  cblas_dtrsv (CblasRowMajor, cblas_UL, cblas_Trans, CblasNonUnit, n,
               ncm_matrix_data (cm), ncm_matrix_tda (cm),
               ncm_vector_data (v), ncm_vector_stride (v));
}

/**
 * ncm_matrix_dsyrk:
 * @cm: a square #NcmMatrix $C$
 * @UL: 'U' or 'L', the triangle of $C$ updated
 * @Trans: 'N' or 'T'
 * @alpha: $\alpha$
 * @A: a #NcmMatrix $A$
 * @beta: $\beta$
 *
 * Sets $C \to \alpha A A^\intercal + \beta C$ for @Trans 'N', or $C \to \alpha A^\intercal A + \beta C$
 * for 'T'. Only the @UL triangle of $C$ is written.
 */
void
ncm_matrix_dsyrk (NcmMatrix *cm, gchar UL, gchar Trans, const gdouble alpha, NcmMatrix *A, const gdouble beta)
{
  const CBLAS_UPLO cblas_UL         = _ncm_matrix_check_uplo ("ncm_matrix_dsyrk", UL);
  const CBLAS_TRANSPOSE cblas_Trans = _ncm_matrix_check_trans ("ncm_matrix_dsyrk", Trans);
  const guint n                     = (cblas_Trans == CblasNoTrans) ? ncm_matrix_nrows (A) : ncm_matrix_ncols (A);
  const guint k                     = (cblas_Trans == CblasNoTrans) ? ncm_matrix_ncols (A) : ncm_matrix_nrows (A);

  g_assert_cmpuint (ncm_matrix_nrows (cm), ==, ncm_matrix_ncols (cm));
  g_assert_cmpuint (ncm_matrix_nrows (cm), ==, n);

  cblas_dsyrk (CblasRowMajor, cblas_UL, cblas_Trans, n, k,
               alpha,
               ncm_matrix_data (A), ncm_matrix_tda (A),
               beta,
               ncm_matrix_data (cm), ncm_matrix_tda (cm));
}

/**
 * ncm_matrix_scale_rows:
 * @cm: a #NcmMatrix $M$
 * @s: a #NcmVector $s$, one entry per row
 *
 * Sets $M \to \mathrm{diag}(s)\,M$.
 */
void
ncm_matrix_scale_rows (NcmMatrix *cm, const NcmVector *s)
{
  const guint nrows = ncm_matrix_nrows (cm);
  guint i;

  g_assert_cmpuint (nrows, ==, ncm_vector_len (s));

  for (i = 0; i < nrows; i++)
    ncm_matrix_mul_row (cm, i, ncm_vector_get (s, i));
}

/**
 * ncm_matrix_scale_cols:
 * @cm: a #NcmMatrix $M$
 * @s: a #NcmVector $s$, one entry per column
 *
 * Sets $M \to M\,\mathrm{diag}(s)$.
 */
void
ncm_matrix_scale_cols (NcmMatrix *cm, const NcmVector *s)
{
  const guint ncols = ncm_matrix_ncols (cm);
  guint j;

  g_assert_cmpuint (ncols, ==, ncm_vector_len (s));

  for (j = 0; j < ncols; j++)
    ncm_matrix_mul_col (cm, j, ncm_vector_get (s, j));
}

/**
 * ncm_matrix_sub_row_vector:
 * @cm: a #NcmMatrix $M$
 * @v: a #NcmVector $v$, one entry per column
 *
 * Sets $M_{ij} \to M_{ij} - v_j$, subtracting @v from every row.
 */
void
ncm_matrix_sub_row_vector (NcmMatrix *cm, const NcmVector *v)
{
  const guint nrows  = ncm_matrix_nrows (cm);
  const guint ncols  = ncm_matrix_ncols (cm);
  const guint stride = ncm_vector_stride (v);
  const gdouble *vd  = ncm_vector_const_data (v);
  guint i, j;

  g_assert_cmpuint (ncols, ==, ncm_vector_len (v));

  for (i = 0; i < nrows; i++)
  {
    gdouble *row = ncm_matrix_ptr (cm, i, 0);

    for (j = 0; j < ncols; j++)
      row[j] -= vd[j * stride];
  }
}

/**
 * ncm_matrix_is_identity:
 * @cm: a square #NcmMatrix
 * @tol: absolute tolerance
 *
 * Aborts if @cm is not square.
 *
 * Returns: whether $\max_{ij} |M_{ij} - \delta_{ij}| <$ @tol.
 */
gboolean
ncm_matrix_is_identity (const NcmMatrix *cm, const gdouble tol)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  guint i, j;

  if (nrows != ncols)
    g_error ("ncm_matrix_is_identity: only works on a square matrix [%ux%u]", nrows, ncols);

  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
    {
      const gdouble dev = fabs (ncm_matrix_get (cm, i, j) - ((i == j) ? 1.0 : 0.0));

      if (!(dev < tol))
        return FALSE;
    }
  }

  return TRUE;
}

/**
 * ncm_matrix_cholesky_decomp:
 * @cm: a #NcmMatrix
 * @UL: 'U' or 'L', the triangle used
 *
 * Replaces the @UL triangle of the symmetric positive definite @cm by its Cholesky factor,
 * with LAPACK dpotrf.
 *
 * Returns: the LAPACK status: zero on success, positive if @cm is not positive definite.
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
 * @UL: 'U' or 'L', the triangle used
 *
 * Replaces the Cholesky factor in the @UL triangle of @cm, from
 * ncm_matrix_cholesky_decomp(), by that triangle of the inverse of the original matrix,
 * with LAPACK dpotri.
 *
 * Returns: the LAPACK status: zero on success.
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
 * @cm holds the Cholesky factor from ncm_matrix_cholesky_decomp(). The product of the
 * diagonal is rescaled as it is accumulated, so it cannot overflow.
 *
 * Returns: $\ln\det A = 2 \sum_i \ln |L_{ii}|$ of the original matrix $A$.
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
 * @b: a #NcmVector $b$ of stride one
 * @UL: 'U' or 'L', the triangle used
 *
 * Solves $A x = b$ for the symmetric positive definite $A$ in @cm, with LAPACK dposv:
 * @cm is replaced by the Cholesky factor in its @UL triangle and @b by $x$.
 *
 * Returns: the LAPACK status: zero on success.
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
 * @b: a #NcmVector $b$
 * @UL: 'U' or 'L', the triangle used
 *
 * Solves $A x = b$ with the Cholesky factor of $A$ already in @cm, from
 * ncm_matrix_cholesky_decomp(), with LAPACK dpotrs; @b is replaced by $x$.
 *
 * Returns: the LAPACK status: zero on success.
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
 * ncm_matrix_chol_chi2_cols:
 * @cm: a #NcmMatrix $X$ of size $d \times n_p$, one point per column
 * @theta: a #NcmVector $\theta$ of length $d$
 * @U: an upper triangular $d \times d$ #NcmMatrix
 * @work: a $d \times n_b$ #NcmMatrix
 * @chi2: a #NcmVector of length at least $n_p$
 *
 * Computes the squared Mahalanobis distance of each column of @cm from @theta under the
 * covariance $C = U^\intercal U$,
 * $$\chi^2_p = (x_p - \theta)^\intercal C^{-1} (x_p - \theta),$$
 * and stores it in the first $n_p$ entries of @chi2. The triangular solve
 * $y_p = (x_p - \theta) U^{-1}$ and the norm $|y_p|^2$ are done in one pass, on blocks of $n_b$
 * columns, with @work as the only scratch.
 *
 * Only the upper triangle of @U is read, its diagonal taken as stored. Since @work belongs to
 * the caller, several threads can share @cm and @U, each with its own @work. The four
 * arguments must be distinct objects that do not overlap.
 */
void
ncm_matrix_chol_chi2_cols (const NcmMatrix *cm, const NcmVector *theta, const NcmMatrix *U, NcmMatrix *work, NcmVector *chi2)
{
  const guint d               = ncm_matrix_nrows (cm);
  const guint np              = ncm_matrix_ncols (cm);
  const guint nb              = ncm_matrix_ncols (work);
  const guint tda_X           = ncm_matrix_tda (cm);
  const guint tda_U           = ncm_matrix_tda (U);
  const guint tda_W           = ncm_matrix_tda (work);
  const guint s_theta         = ncm_vector_stride (theta);
  const gdouble * restrict Xd = ncm_matrix_const_data (cm);
  const gdouble * restrict Ud = ncm_matrix_const_data (U);
  const gdouble * restrict td = ncm_vector_const_data (theta);
  gdouble * restrict Wd       = ncm_matrix_data (work);
  gdouble * restrict c2d      = ncm_vector_data (chi2);
  guint p0;

  g_assert_cmpuint (ncm_matrix_nrows (U), ==, d);
  g_assert_cmpuint (ncm_matrix_ncols (U), ==, d);
  g_assert_cmpuint (ncm_matrix_nrows (work), ==, d);
  g_assert_cmpuint (ncm_vector_len (theta), ==, d);
  g_assert_cmpuint (ncm_vector_len (chi2), >=, np);
  g_assert_cmpuint (ncm_vector_stride (chi2), ==, 1);
  g_assert_cmpuint (nb, >, 0);

  for (p0 = 0; p0 < np; p0 += nb)
  {
    const guint nbp       = MIN (nb, np - p0);
    gdouble * restrict c2 = &c2d[p0];
    guint p, j, k;

    for (j = 0; j < d; j++)
    {
      const gdouble theta_j        = td[j * s_theta];
      const gdouble * restrict X_j = &Xd[j * tda_X + p0];
      gdouble * restrict W_j       = &Wd[j * tda_W];

      for (p = 0; p < nbp; p++)
        W_j[p] = X_j[p] - theta_j;
    }

    for (p = 0; p < nbp; p++)
      c2[p] = 0.0;

    for (j = 0; j < d; j++)
    {
      const gdouble inv_U_jj       = 1.0 / Ud[j * tda_U + j];
      const gdouble * restrict U_j = &Ud[j * tda_U];
      gdouble * restrict W_j       = &Wd[j * tda_W];

      for (p = 0; p < nbp; p++)
      {
        const gdouble y_jp = W_j[p] * inv_U_jj;

        W_j[p] = y_jp;
        c2[p] += y_jp * y_jp;
      }

      for (k = j + 1; k < d; k++)
      {
        const gdouble U_jk     = U_j[k];
        gdouble * restrict W_k = &Wd[k * tda_W];

        for (p = 0; p < nbp; p++)
          W_k[p] -= U_jk * W_j[p];
      }
    }
  }
}

/**
 * ncm_matrix_nearPD:
 * @cm: a #NcmMatrix
 * @UL: 'U' or 'L', the triangle used
 * @cholesky_decomp: whether to leave the Cholesky factor in @cm
 * @maxiter: maximum number of iterations
 *
 * Replaces the symmetric @cm, stored in its @UL triangle, by the nearest positive definite
 * matrix in the Frobenius norm, [Higham (2002)](https://doi.org/10.1093/imanum/22.3.329),
 * iterating until its Cholesky decomposition succeeds or @maxiter is reached.
 *
 * Returns: the status of the last Cholesky decomposition, zero on success.
 */
gint
ncm_matrix_nearPD (NcmMatrix *cm, gchar UL, gboolean cholesky_decomp, const guint maxiter)
{
  const guint n   = ncm_matrix_ncols (cm);
  NcmVector *eva  = ncm_vector_new (n);
  NcmVector *diag = ncm_vector_new (n);
  NcmMatrix *eve  = ncm_matrix_new (n, n);
  NcmMatrix *D_S  = ncm_matrix_new (n, n);
  NcmMatrix *R    = ncm_matrix_new (n, n);
  GArray *isuppz  = g_array_new (FALSE, FALSE, sizeof (gint));
  NcmLapackWS *ws = ncm_lapack_ws_new ();
  gint neva       = 0;
  gint ret, i;
  guint iter;

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
 * ncm_matrix_cholesky_decomp_nearPD:
 * @cm: a symmetric #NcmMatrix
 * @decomp: a #NcmMatrix of the same size
 * @UL: 'U' or 'L', the triangle used
 * @maxiter: iterations of ncm_matrix_nearPD() allowed, zero for none
 * @repaired: (out) (nullable): whether the factor comes from ncm_matrix_nearPD()
 *
 * Stores in @decomp the Cholesky factor of @cm, reading the @UL triangle and leaving @cm
 * unchanged. If @cm is not positive definite to rounding, the factor of the nearest
 * positive definite matrix, from ncm_matrix_nearPD(), is stored instead.
 *
 * Returns: zero on success, otherwise the status of the last Cholesky decomposition.
 */
gint
ncm_matrix_cholesky_decomp_nearPD (const NcmMatrix *cm, NcmMatrix *decomp, gchar UL, const guint maxiter, gboolean *repaired)
{
  gint ret;

  ncm_matrix_memcpy (decomp, cm);
  ret = ncm_matrix_cholesky_decomp (decomp, UL);

  if (repaired != NULL)
    *repaired = FALSE;

  if ((ret != 0) && (maxiter > 0))
  {
    ncm_matrix_memcpy (decomp, cm);
    ret = ncm_matrix_nearPD (decomp, UL, TRUE, maxiter);

    if (repaired != NULL)
      *repaired = TRUE;
  }

  return ret;
}

/**
 * ncm_matrix_sym_exp_cholesky:
 * @cm: a symmetric #NcmMatrix $M$
 * @UL: 'U' or 'L', the triangle used
 * @exp_cm_dec: a #NcmMatrix of the same size
 *
 * Computes the matrix exponential of @cm from its eigendecomposition, and stores in
 * @exp_cm_dec the upper triangular $U$ with $\exp(M) = U^\intercal U$.
 */
void
ncm_matrix_sym_exp_cholesky (NcmMatrix *cm, gchar UL, NcmMatrix *exp_cm_dec)
{
  const guint n   = ncm_matrix_ncols (cm);
  NcmVector *eva  = ncm_vector_new (n);
  NcmMatrix *eve  = exp_cm_dec;
  GArray *tau     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *isuppz  = g_array_new (FALSE, FALSE, sizeof (gint));
  NcmLapackWS *ws = ncm_lapack_ws_new ();
  gint neva       = 0;
  gint ret;
  guint i;

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
 * @cm: a symmetric positive definite #NcmMatrix $M$
 * @UL: 'U' or 'L', the triangle used
 * @ln_cm: a #NcmMatrix of the same size
 *
 * Stores in @ln_cm the matrix logarithm $\ln M$, from the eigendecomposition of @cm.
 */
void
ncm_matrix_sym_posdef_log (NcmMatrix *cm, gchar UL, NcmMatrix *ln_cm)
{
  const guint n   = ncm_matrix_ncols (cm);
  NcmVector *eva  = ncm_vector_new (n);
  NcmMatrix *eve  = ncm_matrix_new (n, n);
  NcmMatrix *temp = ncm_matrix_new (n, n);
  GArray *isuppz  = g_array_new (FALSE, FALSE, sizeof (gint));
  NcmLapackWS *ws = ncm_lapack_ws_new ();
  gint neva       = 0;
  gint ret;
  guint i;

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
 * @cm: a triangular #NcmMatrix $M$
 * @UL: 'U' or 'L', whether $M$ is upper or lower triangular
 * @zero: whether to zero the other triangle first
 * @sym: a #NcmMatrix of the same size
 *
 * Stores in @sym the symmetric $M^\intercal M$ for an upper triangular @cm, or $M M^\intercal$ for a
 * lower one. Unless @zero is %TRUE, the other triangle of @cm must already be zero.
 */
void
ncm_matrix_triang_to_sym (NcmMatrix *cm, gchar UL, gboolean zero, NcmMatrix *sym)
{
  const guint n = ncm_matrix_ncols (cm);
  guint i;

  g_assert_cmpuint (ncm_matrix_ncols (cm), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (sym), ==, ncm_matrix_nrows (cm));
  g_assert_cmpuint (ncm_matrix_ncols (sym), ==, ncm_matrix_nrows (sym));

  if (zero)
  {
    if (UL == 'U')
    {
      for (i = 0; i < n; i++)
      {
        guint j;

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
        guint j;

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
 * @cm: a #NcmMatrix $M$
 * @NT: 'N' or 'T'
 * @UL: 'U' or 'L', the triangle of @sym written
 * @sym: a square #NcmMatrix
 *
 * Stores in the @UL triangle of @sym the symmetric $M M^\intercal$ for @NT 'N', or $M^\intercal M$ for
 * 'T'.
 */
void
ncm_matrix_square_to_sym (NcmMatrix *cm, gchar NT, gchar UL, NcmMatrix *sym)
{
  const guint nrows     = ncm_matrix_nrows (cm);
  const guint ncols     = ncm_matrix_ncols (cm);
  const CBLAS_UPLO Uplo = (UL == 'U') ? CblasUpper : CblasLower;
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
 * @cm: a #NcmMatrix $M$
 * @NT: 'N' or 'T', whether to transpose $M$
 * @alpha: $\alpha$
 * @v: a #NcmVector $v$
 * @beta: $\beta$
 * @u: a #NcmVector $u$
 *
 * Sets $u \to \alpha\,\mathrm{op}(M)\,v + \beta u$; any @NT other than 'N' transposes.
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
 * @cm: a symmetric #NcmMatrix $M$
 * @UL: 'U' or 'L', the triangle used
 * @alpha: $\alpha$
 * @v: a #NcmVector $v$
 * @beta: $\beta$
 * @u: a #NcmVector $u$
 *
 * Sets $u \to \alpha M v + \beta u$, reading only the @UL triangle of $M$.
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
 * @prefix: prefix of each row
 * @format: printf format of one element
 *
 * Logs the elements of @cm, one row per line.
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
 * @cor_level: the parameter $\beta > 0$
 * @rng: a #NcmRNG
 *
 * Replaces @cm by a random correlation matrix from the vine construction of
 * [Lewandowski, Kurowicka and Joe (2009)](https://doi.org/10.1016/j.jmva.2009.04.008): the
 * partial correlations are drawn from Beta($\beta$, $\beta$) mapped to $[-1, 1]$, so a smaller
 * @cor_level gives stronger correlations.
 */
void
ncm_matrix_fill_rand_cor (NcmMatrix *cm, const gdouble cor_level, NcmRNG *rng)
{
  const guint n   = ncm_matrix_nrows (cm);
  const guint nm1 = n - 1;

  g_assert_cmpfloat (cor_level, >, 0.0);
  g_assert_cmpuint (n, ==, ncm_matrix_ncols (cm));
  g_assert_cmpuint (n, >, 0);

  ncm_rng_lock (rng);
  {
    NcmMatrix *P = ncm_matrix_dup (cm);
    guint k;

    ncm_matrix_set_all (P, 0.0);
    ncm_matrix_set_identity (cm);

    for (k = 0; k < nm1; k++)
    {
      guint i;

      for (i = k + 1; i < n; i++)
      {
        gdouble p = (ncm_rng_beta_gen (rng, cor_level, cor_level) - 0.5) * 2.0;
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
 * @sigma_min: smallest standard deviation
 * @sigma_max: largest standard deviation
 * @cor_level: the parameter of ncm_matrix_fill_rand_cor()
 * @rng: a #NcmRNG
 *
 * Replaces @cm by a random covariance matrix: the correlations of
 * ncm_matrix_fill_rand_cor() and standard deviations drawn uniformly in
 * [@sigma_min, @sigma_max].
 */
void
ncm_matrix_fill_rand_cov (NcmMatrix *cm, const gdouble sigma_min, const gdouble sigma_max, const gdouble cor_level, NcmRNG *rng)
{
  const guint n = ncm_matrix_nrows (cm);
  guint k;

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
 * @mu: a #NcmVector $\mu$
 * @reltol_min: smallest relative error
 * @reltol_max: largest relative error
 * @cor_level: the parameter of ncm_matrix_fill_rand_cor()
 * @rng: a #NcmRNG
 *
 * Replaces @cm by a random covariance matrix: the correlations of
 * ncm_matrix_fill_rand_cor() and standard deviations $\sigma_k = |\mu_k| r_k$, or $r_k$ where
 * $\mu_k = 0$, with $r_k$ drawn log-uniformly in [@reltol_min, @reltol_max].
 */
void
ncm_matrix_fill_rand_cov2 (NcmMatrix *cm, NcmVector *mu, const gdouble reltol_min, const gdouble reltol_max, const gdouble cor_level, NcmRNG *rng)
{
  const guint n = ncm_matrix_nrows (cm);
  guint k;

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
 * @cor: a #NcmMatrix of the same size
 *
 * Stores in @cor the correlation matrix of the covariance @cov; they may be the same object.
 */
void
ncm_matrix_cov2cor (const NcmMatrix *cov, NcmMatrix *cor)
{
  const guint n = ncm_matrix_nrows (cov);
  guint i;

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
 * Returns: (transfer full): a new #NcmMatrix with the correlation matrix of the covariance @cov.
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
 * @gm: a #gsl_matrix
 *
 * Creates a constant matrix over @gm, which must outlive it.
 *
 * Returns: a new constant #NcmMatrix.
 */

/**
 * ncm_matrix_get:
 * @cm: a constant #NcmMatrix
 * @i: row index
 * @j: column index
 *
 * Returns: the element ($i$, $j$).
 */

/**
 * ncm_matrix_get_colmajor:
 * @cm: a constant #NcmMatrix
 * @i: row index
 * @j: column index
 *
 * Reads the data of @cm in column-major order.
 *
 * Returns: the element ($i$, $j$) in column-major order.
 */

/**
 * ncm_matrix_ptr:
 * @cm: a #NcmMatrix
 * @i: row index
 * @j: column index
 *
 * Returns: a pointer to the element ($i$, $j$).
 */

/**
 * ncm_matrix_const_ptr:
 * @cm: a constant #NcmMatrix
 * @i: row index
 * @j: column index
 *
 * Returns: a constant pointer to the element ($i$, $j$).
 */

/**
 * ncm_matrix_set:
 * @cm: a #NcmMatrix
 * @i: row index
 * @j: column index
 * @val: a double
 *
 * Sets the element ($i$, $j$) to @val.
 */

/**
 * ncm_matrix_set_colmajor:
 * @cm: a #NcmMatrix
 * @i: row index
 * @j: column index
 * @val: a double
 *
 * Sets the element ($i$, $j$), in column-major order, to @val.
 */

/**
 * ncm_matrix_addto:
 * @cm: a #NcmMatrix
 * @i: row index
 * @j: column index
 * @val: a double
 *
 * Adds @val to the element ($i$, $j$).
 */

/**
 * ncm_matrix_transpose:
 * @cm: a #NcmMatrix
 *
 * Transposes the square @cm in place.
 */

/**
 * ncm_matrix_transpose_memcpy:
 * @cm: a #NcmMatrix
 * @src: a #NcmMatrix
 *
 * Copies the transpose of @src into @cm, whose shape must be that of the transpose.
 */

/**
 * ncm_matrix_set_identity:
 * @cm: a #NcmMatrix
 *
 * Sets @cm, square or not, to one on the diagonal and zero elsewhere.
 */

/**
 * ncm_matrix_set_zero:
 * @cm: a #NcmMatrix
 *
 * Sets every element to zero.
 */

/**
 * ncm_matrix_set_all:
 * @cm: a #NcmMatrix
 * @val: a double
 *
 * Sets every element to @val.
 */

/**
 * ncm_matrix_add:
 * @cm1: a #NcmMatrix
 * @cm2: a constant #NcmMatrix
 *
 * Adds @cm2 to @cm1; they must have the same shape.
 */

/**
 * ncm_matrix_sub:
 * @cm1: a #NcmMatrix
 * @cm2: a constant #NcmMatrix
 *
 * Subtracts @cm2 from @cm1; they must have the same shape.
 */

/**
 * ncm_matrix_mul_elements:
 * @cm1: a #NcmMatrix
 * @cm2: a constant #NcmMatrix
 *
 * Multiplies @cm1 by @cm2, element by element; they must have the same shape.
 */

/**
 * ncm_matrix_div_elements:
 * @cm1: a #NcmMatrix
 * @cm2: a constant #NcmMatrix
 *
 * Divides @cm1 by @cm2, element by element; they must have the same shape.
 */

/**
 * ncm_matrix_scale:
 * @cm: a #NcmMatrix
 * @val: a double
 *
 * Multiplies every element by @val.
 */

/**
 * ncm_matrix_add_constant:
 * @cm: a #NcmMatrix
 * @val: a double
 *
 * Adds @val to every element.
 */

/**
 * ncm_matrix_mul_row:
 * @cm: a #NcmMatrix
 * @row_i: row index
 * @val: a double
 *
 * Multiplies the row @row_i by @val.
 */

/**
 * ncm_matrix_mul_col:
 * @cm: a #NcmMatrix
 * @col_i: column index
 * @val: a double
 *
 * Multiplies the column @col_i by @val.
 */

/**
 * ncm_matrix_get_diag:
 * @cm: a #NcmMatrix
 * @diag: a #NcmVector
 *
 * Copies the diagonal of @cm to the first entries of @diag.
 */

/**
 * ncm_matrix_set_diag:
 * @cm: a #NcmMatrix
 * @diag: a #NcmVector
 *
 * Sets the diagonal of @cm to the first entries of @diag.
 */

/**
 * ncm_matrix_memcpy:
 * @cm1: a #NcmMatrix
 * @cm2: a constant #NcmMatrix
 *
 * Copies @cm2 into @cm1; they must have the same shape.
 */

/**
 * ncm_matrix_memcpy_to_colmajor:
 * @cm1: a #NcmMatrix
 * @cm2: a constant #NcmMatrix
 *
 * Copies @cm2 into @cm1, writing the data of @cm1 in column-major order; they must have
 * the same shape.
 */

/**
 * ncm_matrix_set_col:
 * @cm: a #NcmMatrix
 * @n: column index
 * @cv: a constant #NcmVector
 *
 * Copies @cv into the column @n; its length must be the number of rows.
 */

/**
 * ncm_matrix_set_row:
 * @cm: a #NcmMatrix
 * @n: row index
 * @cv: a constant #NcmVector
 *
 * Copies @cv into the row @n; its length must be the number of columns.
 */

/**
 * ncm_matrix_get_array:
 * @cm: a #NcmMatrix
 *
 * @cm must have been created by ncm_matrix_new_array().
 *
 * Returns: (transfer full) (element-type double): a new reference to the #GArray of @cm.
 */

/**
 * ncm_matrix_dup_array:
 * @cm: a #NcmMatrix
 *
 * Returns: (transfer full) (element-type double): a new #GArray with a copy of the elements, row
 * after row.
 */

/**
 * ncm_matrix_fast_get:
 * @cm: a #NcmMatrix
 * @ij: index into the data
 *
 * Reads the data directly: the element ($i$, $j$) is at @ij $= i\,\mathrm{tda} + j$.
 *
 * Returns: the element at @ij.
 */

/**
 * ncm_matrix_fast_set:
 * @cm: a #NcmMatrix
 * @ij: index into the data
 * @val: a double
 *
 * Writes the data directly: the element ($i$, $j$) is at @ij $= i\,\mathrm{tda} + j$.
 */

/**
 * ncm_matrix_gsl: (skip)
 * @cm: a #NcmMatrix
 *
 * Returns: the #gsl_matrix of @cm.
 */

/**
 * ncm_matrix_const_gsl: (skip)
 * @cm: a constant #NcmMatrix
 *
 * Returns: the constant #gsl_matrix of @cm.
 */

/**
 * ncm_matrix_col_len:
 * @cm: a #NcmMatrix
 *
 * Same as ncm_matrix_nrows().
 *
 * Returns: the number of rows.
 */

/**
 * ncm_matrix_row_len:
 * @cm: a #NcmMatrix
 *
 * Same as ncm_matrix_ncols().
 *
 * Returns: the number of columns.
 */

/**
 * ncm_matrix_nrows:
 * @cm: a #NcmMatrix
 *
 * Returns: the number of rows.
 */

/**
 * ncm_matrix_ncols:
 * @cm: a #NcmMatrix
 *
 * Returns: the number of columns.
 */

/**
 * ncm_matrix_size:
 * @cm: a #NcmMatrix
 *
 * Returns: the number of elements, rows times columns.
 */

/**
 * ncm_matrix_tda:
 * @cm: a #NcmMatrix
 *
 * Returns: the distance between consecutive rows, in doubles.
 */

/**
 * ncm_matrix_data:
 * @cm: a #NcmMatrix
 *
 * Returns: (transfer none): a pointer to the first element.
 */

/**
 * ncm_matrix_const_data:
 * @cm: a constant #NcmMatrix
 *
 * Returns: (transfer none): a constant pointer to the first element.
 */

