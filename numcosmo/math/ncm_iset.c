/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_iset.c
 *
 *  Sun April 7 16:59:36 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_iset.c
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

/**
 * SECTION:ncm_iset
 * @title: NcmISet
 * @short_description: Index set object
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_iset.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmISetPrivate
{
  guint n;
  GQueue *iq;
  gboolean consistent;
  gboolean sorted;
  NcmVector *tmp;
  GArray *ptmp;
  GArray *atmp;
};

enum
{
  PROP_0,
  PROP_N,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmISet, ncm_iset, G_TYPE_OBJECT);

static gint _ncm_iset_cmp (gconstpointer a, gconstpointer b, gpointer user_data) {return GPOINTER_TO_INT (a) - GPOINTER_TO_INT (b); }

static void
ncm_iset_init (NcmISet *iset)
{
  NcmISetPrivate *const self = iset->priv = ncm_iset_get_instance_private (iset);
  self->n          = 0;
  self->iq         = g_queue_new ();
  self->consistent = TRUE;
  self->sorted     = FALSE;
  self->tmp        = NULL;
  self->ptmp       = NULL;
  self->atmp       = NULL;
}

static void _ncm_iset_set_max_size (NcmISet *iset, guint n);

static void
_ncm_iset_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmISet *iset = NCM_ISET (object);
  g_return_if_fail (NCM_IS_ISET (object));

  switch (prop_id)
  {
    case PROP_N:
      _ncm_iset_set_max_size (iset, g_value_get_uint (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_iset_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmISet *iset = NCM_ISET (object);
  g_return_if_fail (NCM_IS_ISET (object));

  switch (prop_id)
  {
    case PROP_N:
      g_value_set_uint (value, ncm_iset_get_max_size (iset));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_iset_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_iset_parent_class)->constructed (object);
  {
    /*NcmISet *iset = NCM_ISET (object);*/

  }
}

static void
_ncm_iset_dispose (GObject *object)
{
  NcmISet *iset = NCM_ISET (object);
  NcmISetPrivate * const self = iset->priv;
  
  ncm_vector_clear (&self->tmp);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_iset_parent_class)->dispose (object);
}

static void
_ncm_iset_finalize (GObject *object)
{
  NcmISet *iset = NCM_ISET (object);
  NcmISetPrivate * const self = iset->priv;

  g_array_unref (self->ptmp);
  g_array_unref (self->atmp);
  g_queue_free (self->iq);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_iset_parent_class)->finalize (object);
}

static void
ncm_iset_class_init (NcmISetClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_iset_set_property;
  object_class->get_property = &_ncm_iset_get_property;
  object_class->constructed  = &_ncm_iset_constructed;
  object_class->dispose      = &_ncm_iset_dispose;
  object_class->finalize     = &_ncm_iset_finalize;

  g_object_class_install_property (object_class,
                                   PROP_N,
                                   g_param_spec_uint ("max-index",
                                                      NULL,
                                                      "Maximum index",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void
_ncm_iset_set_max_size (NcmISet *iset, guint n)
{
  NcmISetPrivate * const self = iset->priv;
  self->n = n;
  ncm_vector_clear (&self->tmp);

  self->tmp  = ncm_vector_new (n);
  self->ptmp = g_array_new (FALSE, FALSE, sizeof (size_t));
  self->atmp = g_array_new (FALSE, FALSE, sizeof (gint));
}

/**
 * ncm_iset_new:
 * @n: maximum index
 * 
 * Creates a new #NcmISet object.
 * 
 * Returns: a new #NcmISet.
 */
NcmISet *
ncm_iset_new (guint n)
{
  NcmISet *iset = g_object_new (NCM_TYPE_ISET,
                                "max-index", n,
                                NULL);
  return iset;
}

/**
 * ncm_iset_ref:
 * @iset: a #NcmISet
 *
 * Increase the reference of @iset by one.
 *
 * Returns: (transfer full): @iset.
 */
NcmISet *
ncm_iset_ref (NcmISet *iset)
{
  return g_object_ref (iset);
}

/**
 * ncm_iset_free:
 * @iset: a #NcmISet
 *
 * Decrease the reference count of @iset by one.
 *
 */
void
ncm_iset_free (NcmISet *iset)
{
  g_object_unref (iset);
}

/**
 * ncm_iset_clear:
 * @iset: a #NcmISet
 *
 * Decrease the reference count of @iset by one, and sets the pointer *@iset to
 * NULL.
 *
 */
void
ncm_iset_clear (NcmISet **iset)
{
  g_clear_object (iset);
}

guint
ncm_iset_get_max_size (NcmISet *iset)
{
  NcmISetPrivate * const self = iset->priv;
  return self->n;
}

/**
 * ncm_iset_add_range:
 * @iset: a #NcmISet
 * @ii: initial index $i_i$
 * @fi: final index $i_f$
 *
 * Adds the interval $(i_i, i_f]$ to the set.
 * Note that $i_f$ is not included.
 *
 */
void
ncm_iset_add_range (NcmISet *iset, gint ii, gint fi)
{
  NcmISetPrivate *const self = iset->priv;
  gint i;

  for (i = ii; i < fi; i++)
  {
    g_queue_push_tail (self->iq, GINT_TO_POINTER (i));
  }
  self->sorted = FALSE;
}

/**
 * ncm_iset_add:
 * @iset: a #NcmISet
 * @i: index $i$
 *
 * Adds the index $i$ to the set.
 *
 */
void
ncm_iset_add (NcmISet *iset, gint i)
{
  NcmISetPrivate *const self = iset->priv;
  if (self->consistent)
    g_assert (g_queue_find (self->iq, GINT_TO_POINTER (i)) == NULL);

  g_queue_push_tail (self->iq, GINT_TO_POINTER (i));
  self->sorted = FALSE;
}

/**
 * ncm_iset_del:
 * @iset: a #NcmISet
 * @i: index $i$
 *
 * Removes the index $i$ from the set.
 *
 */
void
ncm_iset_del (NcmISet *iset, gint i)
{
  NcmISetPrivate *const self = iset->priv;
  if (self->consistent)
    g_assert (g_queue_remove (self->iq, GINT_TO_POINTER (i)));
  else
    g_queue_remove (self->iq, GINT_TO_POINTER (i));
}

/**
 * ncm_iset_reset:
 * @iset: a #NcmISet
 *
 * Removes all indexes from the set.
 *
 */
void
ncm_iset_reset (NcmISet *iset)
{
  NcmISetPrivate *const self = iset->priv;
  g_queue_clear (self->iq);
  self->sorted = FALSE;
}

typedef struct _NcmISetData
{
  NcmISetPrivate *const self;
  NcmVector *v;
  gdouble max;
  gint max_i;
} NcmISetData;

static void
_ncm_iset_sort (NcmISet *iset)
{
  NcmISetPrivate *const self = iset->priv;
  if (!self->sorted)
  {
    g_queue_sort (self->iq, _ncm_iset_cmp, NULL);
    self->sorted = TRUE;
  }
}

/**
 * ncm_iset_memcpy:
 * @iset: a #NcmISet
 * @target: a #NcmISet
 *
 * Copy the set @iset over @target.
 *
 */
void
ncm_iset_copy (NcmISet *iset, NcmISet *target)
{
  NcmISetPrivate *const self = iset->priv;
  GList *node;

  g_assert_cmpint (ncm_iset_get_max_size (iset), ==, ncm_iset_get_max_size (target));

  _ncm_iset_sort (iset);
  ncm_iset_reset (target);

  node = g_queue_peek_head_link (self->iq);

  while (node != NULL)
  {
    const gint i = GPOINTER_TO_INT (node->data);
    ncm_iset_add (target, i);
    node = node->next;
  }
}

/**
 * ncm_iset_get_len:
 * @iset: a #NcmISet
 *
 * Returns: number of elements in @iset.
 */
guint
ncm_iset_get_len (NcmISet *iset)
{
  NcmISetPrivate *const self = iset->priv;
  return g_queue_get_length (self->iq);
}

/**
 * ncm_iset_get_vector_max:
 * @iset: a #NcmISet
 * @v: a #NcmVector
 * @max_i: (out): Maximum component index
 *
 * Finds the maximum component of the vector @v.
 *
 */
gdouble
ncm_iset_get_vector_max (NcmISet *iset, NcmVector *v, gint *max_i)
{
  NcmISetPrivate *const self = iset->priv;
  gdouble max = GSL_NEGINF;
  GList *node;

  _ncm_iset_sort (iset);
  g_assert_cmpint (ncm_vector_len (v), >=, self->n);

  node = g_queue_peek_head_link (self->iq);
  max_i[0] = -1;

  while (node != NULL)
  {
    gdouble v_i;
    if ((v_i = ncm_vector_get (v, GPOINTER_TO_INT (node->data))) > max)
    {
      max      = v_i;
      max_i[0] = GPOINTER_TO_INT (node->data);
    }

    node = node->next;
  }

  return max;
}

/**
 * ncm_iset_get_subvector:
 * @iset: a #NcmISet
 * @v: a #NcmVector
 * @v_dup: a #NcmVector
 *
 * Construct a continuous vector $s$ using the values from @v
 * and the indexes in @iset. If @v_dup is not null use
 * this vector to build the subvector, otherwise, allocates
 * a new vector.
 *
 * Returns: (transfer full): the vector $s$.
 */
NcmVector *
ncm_iset_get_subvector (NcmISet *iset, NcmVector *v, NcmVector *v_dup)
{
  NcmISetPrivate *const self = iset->priv;
  const guint nsub = g_queue_get_length (self->iq);
  NcmVector *sub;
  GList *node;
  gint i;

  g_assert_cmpuint (self->n, ==, ncm_vector_len (v));

  if (v_dup != NULL)
  {
    g_assert_cmpuint (self->n, ==, ncm_vector_len (v_dup));
    sub = ncm_vector_get_subvector (v_dup, 0, nsub);
  }
  else
    sub = ncm_vector_new (nsub);

  _ncm_iset_sort (iset);
  node = g_queue_peek_head_link (self->iq);
  i    = 0;

  while (node != NULL)
  {
    const gint j = GPOINTER_TO_INT (node->data);

    ncm_vector_set (sub, i, ncm_vector_get (v, j));

    i++;
    node = node->next;
  }

  return sub;
}

/**
 * ncm_iset_get_submatrix:
 * @iset: a #NcmISet
 * @M: a #NcmMatrix
 * @M_dup: a #NcmMatrix
 *
 * Construct a continuous matrix square $S$ using the values
 * from the square matrix @M and the indexes in @iset. If
 * @M_dup is not null use this matrix to build the submatrix,
 * otherwise, allocates a new matrix.
 *
 * Returns: (transfer full): the matrix $S$.
 */
NcmMatrix *
ncm_iset_get_submatrix (NcmISet *iset, NcmMatrix *M, NcmMatrix *M_dup)
{
  NcmISetPrivate *const self = iset->priv;
  const guint nsub = g_queue_get_length (self->iq);
  NcmMatrix *sub;
  GList *node0, *node;
  gint i, j;

  g_assert_cmpuint (self->n, ==, ncm_matrix_nrows (M));
  g_assert_cmpuint (self->n, ==, ncm_matrix_ncols (M));

  if (M_dup != NULL)
  {
    g_assert_cmpuint (self->n, ==, ncm_matrix_nrows (M_dup));
    g_assert_cmpuint (self->n, ==, ncm_matrix_ncols (M_dup));

    sub = ncm_matrix_get_submatrix (M_dup, 0, 0, nsub, nsub);
  }
  else
    sub = ncm_matrix_new (nsub, nsub);

  _ncm_iset_sort (iset);
  node0 = node = g_queue_peek_head_link (self->iq);
  i     = 0;

  while (node != NULL)
  {
    const gint k = GPOINTER_TO_INT (node->data);
    GList *nodei = node0;

    j = 0;
    while (nodei != NULL)
    {
      const gint l = GPOINTER_TO_INT (nodei->data);

      ncm_matrix_set (sub, i, j, ncm_matrix_get (M, k, l));

      nodei = nodei->next;
      j++;
    }

    i++;
    node = node->next;
  }

  return sub;
}

/**
 * ncm_iset_get_submatrix_cols:
 * @iset: a #NcmISet
 * @M: a #NcmMatrix
 * @M_dup: a #NcmMatrix
 *
 * Construct a continuous matrix retangular $S$ using the columns
 * from the rectangular matrix @M and the indexes in @iset. If
 * @M_dup is not null use this matrix to build the submatrix,
 * otherwise, allocates a new matrix.
 *
 * Returns: (transfer full): the matrix $S$.
 */
NcmMatrix *
ncm_iset_get_submatrix_cols (NcmISet *iset, NcmMatrix *M, NcmMatrix *M_dup)
{
  NcmISetPrivate *const self = iset->priv;
  const guint nsub  = g_queue_get_length (self->iq);
  const guint nrows = ncm_matrix_nrows (M);
  const guint ncols = ncm_matrix_ncols (M);
  NcmMatrix *sub;
  GList *node;
  gint i;

  g_assert_cmpuint (self->n, ==, ncols);

  if (M_dup != NULL)
  {
    g_assert_cmpuint (ncols, ==, ncm_matrix_ncols (M_dup));
    g_assert_cmpuint (nrows, ==, ncm_matrix_nrows (M_dup));

    sub = ncm_matrix_get_submatrix (M_dup, 0, 0, nrows, nsub);
  }
  else
    sub = ncm_matrix_new (nrows, nsub);

  _ncm_iset_sort (iset);
  node = g_queue_peek_head_link (self->iq);
  i    = 0;

  while (node != NULL)
  {
    const gint k = GPOINTER_TO_INT (node->data);
    NcmVector *col_k = ncm_matrix_get_col (M, k);
    NcmVector *col_i = ncm_matrix_get_col (sub, i);

    ncm_vector_memcpy (col_i, col_k);

    i++;
    node = node->next;
  }

  return sub;
}


/**
 * ncm_iset_get_sym_submatrix:
 * @iset: a #NcmISet
 * @UL: char indicating 'U'pper or 'L'ower matrix
 * @M: a #NcmMatrix
 * @M_dup: a #NcmMatrix
 *
 * Construct a continuous symmetric matrix $S$ using the values
 * from @M and the indexes in @iset. If @M_dup is not null use
 * this matrix to build the submatrix, otherwise, allocates
 * a new matrix. If @UL == 'U'/'L' only the Upper/Lower triangle
 * will be copied.
 *
 * Returns: (transfer full): the matrix $S$.
 */
NcmMatrix *
ncm_iset_get_sym_submatrix (NcmISet *iset, gchar UL, NcmMatrix *M, NcmMatrix *M_dup)
{
  NcmISetPrivate *const self = iset->priv;
  const guint nsub = g_queue_get_length (self->iq);
  NcmMatrix *sub;
  GList *node;
  gint i, j;

  g_assert_cmpuint (self->n, ==, ncm_matrix_nrows (M));
  g_assert_cmpuint (self->n, ==, ncm_matrix_ncols (M));

  if (M_dup != NULL)
  {
    g_assert_cmpuint (self->n, ==, ncm_matrix_nrows (M_dup));
    g_assert_cmpuint (self->n, ==, ncm_matrix_ncols (M_dup));

    sub = ncm_matrix_get_submatrix (M_dup, 0, 0, nsub, nsub);
  }
  else
    sub = ncm_matrix_new (nsub, nsub);

  _ncm_iset_sort (iset);
  node = g_queue_peek_head_link (self->iq);
  i    = 0;

  switch (UL)
  {
    case 'U':
      while (node != NULL)
      {
        const gint k = GPOINTER_TO_INT (node->data);
        GList *nodei = node;

        j = i;
        while (nodei != NULL)
        {
          const gint l = GPOINTER_TO_INT (nodei->data);

          ncm_matrix_set (sub, i, j, ncm_matrix_get (M, k, l));

          nodei = nodei->next;
          j++;
        }

        i++;
        node = node->next;
      }
      break;
    case 'L':
      while (node != NULL)
      {
        const gint k = GPOINTER_TO_INT (node->data);
        GList *nodei = node;

        j = i;
        while (nodei != NULL)
        {
          const gint l = GPOINTER_TO_INT (nodei->data);

          ncm_matrix_set (sub, j, i, ncm_matrix_get (M, l, k));

          nodei = nodei->next;
          j++;
        }

        i++;
        node = node->next;
      }
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  return sub;
}

/**
 * ncm_iset_get_subset_vec_lt:
 * @iset: a #NcmISet
 * @out: a #NcmISet
 * @v: a #NcmVector
 * @tol: a double $t$
 *
 * Gets the subset (@out) of @iset where $v_i < t$.
 *
 */
void
ncm_iset_get_subset_vec_lt (NcmISet *iset, NcmISet *out, NcmVector *v, const gdouble tol)
{
  NcmISetPrivate *const self = iset->priv;
  GList *node;

  ncm_iset_reset (out);
  _ncm_iset_sort (iset);

  g_assert_cmpint (ncm_vector_len (v), >=, self->n);
  g_assert_cmpint (ncm_iset_get_max_size (out), ==, self->n);

  node = g_queue_peek_head_link (self->iq);

  while (node != NULL)
  {
    const gint i = GPOINTER_TO_INT (node->data);

    if (ncm_vector_get (v, i) < tol)
      ncm_iset_add (out, i);

    node = node->next;
  }
}

/**
 * ncm_iset_remove_subset:
 * @iset: a #NcmISet
 * @target: a #NcmISet
 *
 * Removes indexes of @iset from @target.
 *
 */
void
ncm_iset_remove_subset (NcmISet *iset, NcmISet *target)
{
  NcmISetPrivate *const self = iset->priv;
  GList *node;

  _ncm_iset_sort (iset);

  g_assert_cmpint (ncm_iset_get_max_size (iset), ==, ncm_iset_get_max_size (target));

  node = g_queue_peek_head_link (self->iq);

  while (node != NULL)
  {
    const gint i = GPOINTER_TO_INT (node->data);
    ncm_iset_del (target, i);
    node = node->next;
  }
}

/**
 * ncm_iset_remove_smallest_subset:
 * @iset: a #NcmISet
 * @target: a #NcmISet
 * @v: a #NcmVector
 * @max_remove: maximum number of indexes to be removed
 *
 * Removes indexes of @iset from @target based on the values on @v.
 * The first @max_remove indexes from @target matching the smallest
 * components of @v are removed.
 *
 * Returns: number of indexes removed
 */
guint
ncm_iset_remove_smallest_subset (NcmISet *iset, NcmISet *target, NcmVector *v, guint max_remove)
{
  NcmISetPrivate *const self = iset->priv;
  const gint rsize = ncm_iset_get_len (iset);
  GList *node;

  g_assert_cmpuint (max_remove, >, 0);
  g_assert_cmpint (ncm_vector_len (v), >=, ncm_iset_get_max_size (iset));

  if (max_remove >= rsize)
  {
    ncm_iset_remove_subset (iset, target);
    return rsize;
  }
  else
  {
    NcmVector *invalid_v = ncm_vector_get_subvector (self->tmp, 0, rsize);
    gint j;

    g_assert_cmpint (ncm_iset_get_max_size (iset), ==, ncm_iset_get_max_size (target));

    _ncm_iset_sort (iset);

    node = g_queue_peek_head_link (self->iq);
    j    = 0;

    g_array_set_size (self->atmp, rsize);
    g_array_set_size (self->ptmp, max_remove);

    while (node != NULL)
    {
      const gint i = GPOINTER_TO_INT (node->data);

      ncm_vector_set (invalid_v, j, ncm_vector_get (v, i));
      g_array_index (self->atmp, gint, j) = i;

      j++;

      if (j >= rsize)
        break;

      node = node->next;
    }

    gsl_sort_vector_smallest_index (&g_array_index (self->ptmp, size_t, 0), max_remove, ncm_vector_gsl (invalid_v));

    for (j = 0; j < max_remove; j++)
    {
      const gint vi = g_array_index (self->ptmp, size_t, j);
      const gint ti = g_array_index (self->atmp, gint, vi);

      ncm_iset_del (target, ti);
    }

    ncm_vector_free (invalid_v);

    return max_remove;
  }
}

/**
 * ncm_iset_add_largest_subset:
 * @iset: a #NcmISet
 * @v: a #NcmVector
 * @min: a double $\mu$
 * @add_frac: fraction of indexes to be added
 *
 * Adds indexes to @iset using the largest values of @v
 * satisfying $v_i > \mu$ where $i \in $ complement of @iset.
 *
 * Returns: number of indexes added.
 */
guint
ncm_iset_add_largest_subset (NcmISet *iset, NcmVector *v, const gdouble min, const gdouble add_frac)
{
  NcmISetPrivate *const self = iset->priv;
  const gint max_size = ncm_iset_get_max_size (iset);
  const gint rsize    = ncm_iset_get_len (iset);
  const gint csize    = max_size - rsize;
  NcmVector *v_cmplm;
  GList *node;
  guint adds;
  gint j, k;

  if (csize == 0)
    return 0;

  g_assert_cmpfloat (add_frac, >, 0.0);
  g_assert_cmpfloat (add_frac, <=, 1.0);
  g_assert_cmpuint (csize, >, 0);
  g_assert_cmpint (ncm_vector_len (v), ==, ncm_iset_get_max_size (iset));

  _ncm_iset_sort (iset);

  node = g_queue_peek_head_link (self->iq);
  j    = 0;
  k    = 0;

  v_cmplm = ncm_vector_get_subvector (self->tmp, 0, csize);
  ncm_vector_set_zero (v_cmplm);
  g_array_set_size (self->atmp, csize);
  g_array_set_size (self->ptmp, csize);

  while (node != NULL)
  {
    const gint i = GPOINTER_TO_INT (node->data);
    for (; j < i; j++)
    {
      const gdouble v_j = ncm_vector_get (v, j);
      if (v_j > min)
      {
        /*ncm_message ("Adding to cmplm %d % .2e\n", j, v_j);*/
        ncm_vector_set (v_cmplm, k, v_j);
        g_array_index (self->atmp, gint, k) = j;
        k++;
      }
    }
    j = i + 1;
    node = node->next;
  }

  for (; j < max_size; j++)
  {
    const gdouble v_j = ncm_vector_get (v, j);
    if (v_j > min)
    {
      /*ncm_message ("Adding to cmplm %d % .2e\n", j, v_j);*/
      ncm_vector_set (v_cmplm, k, v_j);
      g_array_index (self->atmp, gint, k) = j;
      k++;
    }
  }

  adds = MIN (k, MAX (k * add_frac, 1));

  /*ncm_vector_log_vals (v_cmplm, "v_cmplm: ", "% .2e", TRUE);*/

  if (adds > 0)
  {
    gsl_sort_largest_index (&g_array_index (self->ptmp, size_t, 0), adds, ncm_vector_data (v_cmplm), ncm_vector_stride (v_cmplm), k);

    for (j = 0; j < adds; j++)
    {
      const gint vi = g_array_index (self->ptmp, size_t, j);
      const gint ti = g_array_index (self->atmp, gint, vi);

      ncm_iset_add (iset, ti);
    }
  }

  ncm_vector_free (v_cmplm);

  return adds;
}



/**
 * ncm_iset_set_complement:
 * @iset: a #NcmISet
 * @cmplm: a #NcmISet
 *
 * Sets @cmplm as the complement of @iset.
 *
 */
void
ncm_iset_set_complement (NcmISet *iset, NcmISet *cmplm)
{
  NcmISetPrivate *const self = iset->priv;
  const guint n = ncm_iset_get_max_size (iset);
  GList *node;

  g_assert_cmpint (n, ==, ncm_iset_get_max_size (cmplm));

  ncm_iset_reset (cmplm);
  if (ncm_iset_get_len (iset) == n)
    return;

  _ncm_iset_sort (iset);
  ncm_iset_add_range (cmplm, 0, n);

  node = g_queue_peek_head_link (self->iq);

  while (node != NULL)
  {
    const gint i = GPOINTER_TO_INT (node->data);
    ncm_iset_del (cmplm, i);
    node = node->next;
  }
}

/**
 * ncm_iset_get_vector_inv_cmp:
 * @iset: a #NcmISet
 * @u: a #NcmVector
 * @v: a #NcmVector
 * @v_dup: a #NcmVector
 *
 * Computes the inverse of the relative difference between
 * vectors @u and @v, namely: $$\left(\frac{u_i - v_i}{u_i}\right)^{-1},$$
 * for indexes $i \in $ @iset.
 *
 * Returns: (transfer full): the subset
 */
NcmVector *
ncm_iset_get_vector_inv_cmp (NcmISet *iset, NcmVector *u, NcmVector *v, NcmVector *v_dup)
{
  NcmISetPrivate *const self = iset->priv;
  const guint nsub = g_queue_get_length (self->iq);
  NcmVector *cmp;
  GList *node;
  gint i;

  g_assert_cmpuint (self->n, ==, ncm_vector_len (v));
  g_assert_cmpuint (self->n, ==, ncm_vector_len (u));

  if (v_dup != NULL)
  {
    g_assert_cmpuint (self->n, ==, ncm_vector_len (v_dup));
    cmp = ncm_vector_get_subvector (v_dup, 0, nsub);
  }
  else
    cmp = ncm_vector_new (nsub);

  _ncm_iset_sort (iset);
  node = g_queue_peek_head_link (self->iq);
  i    = 0;

  while (node != NULL)
  {
    const gint j = GPOINTER_TO_INT (node->data);
    const gdouble u_i = ncm_vector_get (u, j);
    const gdouble v_i = ncm_vector_get (v, j);

    ncm_vector_set (cmp, i, u_i / (u_i - v_i));

    i++;
    node = node->next;
  }

  return cmp;
}

/**
 * ncm_iset_set_subvector:
 * @iset: a #NcmISet
 * @v: a #NcmVector
 * @sub: a #NcmVector
 *
 * Copies the components from @sub to the indexes @iset
 * in @v.
 */
void
ncm_iset_set_subvector (NcmISet *iset, NcmVector *v, NcmVector *sub)
{
  NcmISetPrivate *const self = iset->priv;
  GList *node;
  gint j;

  _ncm_iset_sort (iset);
  g_assert_cmpint (ncm_vector_len (v), >=, self->n);
  g_assert_cmpint (ncm_vector_len (v), >=, ncm_vector_len (sub));

  node = g_queue_peek_head_link (self->iq);

  j = 0;
  while (node != NULL)
  {
    const gint i = GPOINTER_TO_INT (node->data);

    ncm_vector_set (v, i, ncm_vector_get (sub, j));

    j++;
    node = node->next;
  }
}

/**
 * ncm_iset_log_vals:
 * @iset: a #NcmISet
 * @prefix: a string
 *
 * Logs the indexes on @iset with prefix @prefix.
 *
 */
void
ncm_iset_log_vals (NcmISet *iset, const gchar *prefix)
{
  NcmISetPrivate *const self = iset->priv;
  GList *node;

  _ncm_iset_sort (iset);

  ncm_message ("%s:", prefix);

  node = g_queue_peek_head_link (self->iq);

  while (node != NULL)
  {
    const gint i = GPOINTER_TO_INT (node->data);
    ncm_message (" %d", i);
    node = node->next;
  }
  ncm_message ("\n");
}
