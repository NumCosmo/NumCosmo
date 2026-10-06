/***************************************************************************
 *            ncm_vector.c
 *
 *  Tue Jul  8 15:05:41 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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
 * NcmVector:
 *
 * Reference-counted vector of doubles, a view over a #gsl_vector.
 *
 * The data can be allocated by the vector or come from GSL, a #GArray, a #GVariant, another
 * #NcmVector (a subvector) or a user array, with the ownership stated by each constructor.
 * Components are separated by a stride; the "fast" accessors assume a stride of one.
 * Functions taking two vectors require the same length unless stated otherwise.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_cblas.h" /* This must be included before any gsl header */
#include "ncm/algebra/ncm_vector.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_VALS,
};

G_DEFINE_TYPE (NcmVector, ncm_vector, G_TYPE_OBJECT)

static void
ncm_vector_init (NcmVector *cv)
{
  cv->pdata = NULL;
  cv->pfree = NULL;
  cv->type  = 0;
  memset (&cv->vv, 0, sizeof (gsl_vector_view));
}

static void
_ncm_vector_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmVector *cv = NCM_VECTOR (object);

  g_return_if_fail (NCM_IS_VECTOR (object));

  switch (prop_id)
  {
    case PROP_VALS:
    {
      GVariant *var = ncm_vector_get_variant (cv);

      g_value_take_variant (value, var);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_vector_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmVector *cv = NCM_VECTOR (object);

  g_return_if_fail (NCM_IS_VECTOR (object));

  switch (prop_id)
  {
    case PROP_VALS:
    {
      GVariant *var = g_value_get_variant (value);

      ncm_vector_set_from_variant (cv, var);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_vector_dispose (GObject *object)
{
  NcmVector *cv = NCM_VECTOR (object);

  if (cv->pdata != NULL)
  {
    g_assert (cv->pfree != NULL);
    cv->pfree (cv->pdata);
    cv->pdata = NULL;
    cv->pfree = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_vector_parent_class)->dispose (object);
}

static void
_ncm_vector_finalize (GObject *object)
{
  NcmVector *cv = NCM_VECTOR (object);

  switch (cv->type)
  {
    case NCM_VECTOR_SLICE:
      g_slice_free1 (sizeof (gdouble) * ncm_vector_len (cv) * ncm_vector_stride (cv), ncm_vector_data (cv));
      break;
    case NCM_VECTOR_ARRAY:
    case NCM_VECTOR_MALLOC:
    case NCM_VECTOR_GSL_VECTOR:
    case NCM_VECTOR_DERIVED:
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  cv->vv.vector.data = NULL;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_vector_parent_class)->finalize (object);
}

static void
ncm_vector_class_init (NcmVectorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_vector_set_property;
  object_class->get_property = &_ncm_vector_get_property;
  object_class->dispose      = &_ncm_vector_dispose;
  object_class->finalize     = &_ncm_vector_finalize;

  /**
   * NcmVector:values:
   *
   * The components as a #GVariant of type `ad`, used for serialization.
   */
  g_object_class_install_property (object_class, PROP_VALS,
                                   g_param_spec_variant ("values", NULL, "values",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_vector_new:
 * @n: number of components
 *
 * Allocates a vector of @n components, not initialized.
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new (gsize n)
{
  gdouble *d = g_slice_alloc (sizeof (gdouble) * n);

  return ncm_vector_new_data_slice (d, n, 1);
}

/**
 * ncm_vector_new_full:
 * @d: (array) (element-type double): the data
 * @size: number of components
 * @stride: distance between consecutive components, in doubles
 * @pdata: (allow-none): data owned by the vector
 * @pfree: (scope notified) (allow-none): function releasing @pdata
 *
 * Creates a vector over @d; @pfree is called on @pdata when the vector is finalized.
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_full (gdouble *d, gsize size, gsize stride, gpointer pdata, GDestroyNotify pfree)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);

  g_assert (d != NULL);
  g_assert_cmpuint (size,   >, 0);
  g_assert_cmpuint (stride, >, 0);

  if (stride != 1)
    cv->vv = gsl_vector_view_array_with_stride (d, stride, size);
  else
    cv->vv = gsl_vector_view_array (d, size);

  cv->type = NCM_VECTOR_DERIVED;
  g_assert ((pdata == NULL) || (pdata != NULL && pfree != NULL));
  cv->pdata = pdata;
  cv->pfree = pfree;

  return cv;
}

/**
 * ncm_vector_new_fftw:
 * @size: number of components
 *
 * Allocates a vector of @size components with fftw_alloc_real(), aligned for FFTW.
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_fftw (guint size)
{
  gdouble *d    = fftw_alloc_real (size);
  NcmVector *cv = ncm_vector_new_full (d, size, 1, d, (GDestroyNotify) fftw_free);

  cv->type = NCM_VECTOR_MALLOC;

  return cv;
}

/**
 * ncm_vector_new_gsl: (skip)
 * @gv: a #gsl_vector
 *
 * Creates a vector over @gv, which it takes ownership of and frees with gsl_vector_free().
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_gsl (gsl_vector *gv)
{
  NcmVector *cv = ncm_vector_new_full (gv->data, gv->size, gv->stride, gv, (GDestroyNotify) gsl_vector_free);

  cv->type = NCM_VECTOR_GSL_VECTOR;

  return cv;
}

/**
 * ncm_vector_new_gsl_static: (skip)
 * @gv: a #gsl_vector
 *
 * Creates a vector over @gv, which must outlive it.
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_gsl_static (gsl_vector *gv)
{
  NcmVector *cv = ncm_vector_new_full (gv->data, gv->size, gv->stride, NULL, NULL);

  cv->type = NCM_VECTOR_GSL_VECTOR;

  return cv;
}

/**
 * ncm_vector_new_array:
 * @a: (array) (element-type double): a #GArray of doubles
 *
 * Creates a vector over the elements of @a, holding a reference to it. @a must not be empty.
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_array (GArray *a)
{
  g_assert_cmpint (a->len, >, 0);
  {
    NcmVector *cv = ncm_vector_new_full (&g_array_index (a, gdouble, 0), a->len, 1,
                                         g_array_ref (a), (GDestroyNotify) & g_array_unref);

    cv->type = NCM_VECTOR_ARRAY;

    return cv;
  }
}

/**
 * ncm_vector_new_data_slice:
 * @d: (array) (element-type double): the data
 * @size: number of components
 * @stride: distance between consecutive components, in doubles
 *
 * Creates a vector over @d, allocated with g_slice_alloc(), which it takes ownership of.
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_slice (gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = ncm_vector_new_full (d, size, stride,
                                       NULL, NULL);

  cv->type = NCM_VECTOR_SLICE;

  return cv;
}

/**
 * ncm_vector_new_data_malloc:
 * @d: (array) (element-type double): the data
 * @size: number of components
 * @stride: distance between consecutive components, in doubles
 *
 * Creates a vector over @d, allocated with g_malloc() or malloc(), which it takes
 * ownership of and frees with g_free().
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_malloc (gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = ncm_vector_new_full (d, size, stride,
                                       d, &g_free);

  cv->type = NCM_VECTOR_MALLOC;

  return cv;
}

/**
 * ncm_vector_new_data_static:
 * @d: (array) (element-type double): the data
 * @size: number of components
 * @stride: distance between consecutive components, in doubles
 *
 * Creates a vector over @d, which must outlive it.
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_static (gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = ncm_vector_new_full (d, size, stride,
                                       NULL, NULL);

  cv->type = NCM_VECTOR_DERIVED;

  return cv;
}

/**
 * ncm_vector_new_data_dup:
 * @d: (array) (element-type double): the data
 * @size: number of components
 * @stride: distance between consecutive components, in doubles
 *
 * Creates a vector with a copy of @d.
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_dup (gdouble *d, const gsize size, const gsize stride)
{
  NcmVector *s   = ncm_vector_new_data_static (d, size, stride);
  NcmVector *dup = ncm_vector_dup (s);

  ncm_vector_free (s);

  return dup;
}

/**
 * ncm_vector_new_variant:
 * @var: a #GVariant of type `ad`
 *
 * Creates a vector with a copy of the elements of @var.
 *
 * Returns: (transfer full): a new #NcmVector.
 */
NcmVector *
ncm_vector_new_variant (GVariant *var)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR,
                                "values", var,
                                NULL);

  return cv;
}

/**
 * ncm_vector_const_new_data:
 * @d: (array) (element-type double): the data
 * @size: number of components
 * @stride: distance between consecutive components, in doubles
 *
 * Creates a constant vector over @d, which must outlive it.
 *
 * Returns: (transfer full): a new constant #NcmVector.
 */
const NcmVector *
ncm_vector_const_new_data (const gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);

  if (stride != 1)
    cv->vv = gsl_vector_view_array_with_stride ((gdouble *) d, stride, size);
  else
    cv->vv = gsl_vector_view_array ((gdouble *) d, size);

  cv->pdata = NULL;
  cv->pfree = NULL;

  cv->type = NCM_VECTOR_DERIVED;

  return cv;
}

/**
 * ncm_vector_ref:
 * @cv: a #NcmVector
 *
 * Increases the reference count of @cv by one.
 *
 * Returns: (transfer full): @cv.
 */
NcmVector *
ncm_vector_ref (NcmVector *cv)
{
  return g_object_ref (cv);
}

/**
 * ncm_vector_const_ref:
 * @cv: a constant #NcmVector
 *
 * Increases the reference count of @cv by one.
 *
 * Returns: (transfer full): @cv.
 */
const NcmVector *
ncm_vector_const_ref (const NcmVector *cv)
{
  return g_object_ref (NCM_VECTOR ((NcmVector *) cv));
}

/**
 * ncm_vector_const_new_variant:
 * @var: a #GVariant of type `ad`
 *
 * Creates a constant vector over the data of @var, holding a reference to it.
 *
 * Returns: (transfer full): a new constant #NcmVector.
 */
const NcmVector *
ncm_vector_const_new_variant (GVariant *var)
{
  gsize n            = g_variant_n_children (var);
  gconstpointer data = g_variant_get_data (var);
  NcmVector *cv      = (NcmVector *) ncm_vector_const_new_data (data, n, 1);

  NCM_VECTOR (cv)->pdata = g_variant_ref_sink (var);
  NCM_VECTOR (cv)->pfree = (GDestroyNotify) & g_variant_unref;

  return cv;
}

/**
 * ncm_vector_free:
 * @cv: a #NcmVector
 *
 * Decreases the reference count of @cv by one.
 */
void
ncm_vector_free (NcmVector *cv)
{
  g_object_unref (cv);
}

/**
 * ncm_vector_const_free:
 * @cv: a constant #NcmVector
 *
 * Decreases the reference count of @cv by one.
 */
void
ncm_vector_const_free (const NcmVector *cv)
{
  ncm_vector_free (NCM_VECTOR ((NcmVector *) cv));
}

/**
 * ncm_vector_clear:
 * @cv: a #NcmVector
 *
 * If *@cv is not %NULL, decreases its reference count by one and sets *@cv to %NULL.
 */
void
ncm_vector_clear (NcmVector **cv)
{
  g_clear_object (cv);
}

/**
 * ncm_vector_dup:
 * @cv: a constant #NcmVector
 *
 * Returns: (transfer full): a newly allocated copy of @cv, with stride one.
 */
NcmVector *
ncm_vector_dup (const NcmVector *cv)
{
  NcmVector *cv_cp = ncm_vector_new (ncm_vector_len (cv));

  gsl_vector_memcpy (ncm_vector_gsl (cv_cp), ncm_vector_const_gsl (cv));

  return cv_cp;
}

/**
 * ncm_vector_substitute:
 * @cv1: a #NcmVector
 * @cv2: (nullable): a #NcmVector
 * @check_size: whether to require the same length
 *
 * Replaces *@cv1 by a new reference to @cv2, releasing the previous one. With @check_size,
 * aborts if both are set and their lengths differ.
 */
void
ncm_vector_substitute (NcmVector **cv1, NcmVector *cv2, gboolean check_size)
{
  if (*cv1 == cv2)
    return;

  if (*cv1 != NULL)
  {
    if ((cv2 != NULL) && check_size)
      g_assert_cmpuint (ncm_vector_len (*cv1), ==, ncm_vector_len (cv2));

    ncm_vector_clear (cv1);
  }

  if (cv2 != NULL)
    *cv1 = ncm_vector_ref (cv2);
}

/**
 * ncm_vector_get_subvector:
 * @cv: a #NcmVector
 * @k: index of the first component
 * @size: number of components
 *
 * Creates a view of the components @k to @k + @size - 1 of @cv, holding a reference to it.
 *
 * Returns: (transfer full): the subvector.
 */
NcmVector *
ncm_vector_get_subvector (NcmVector *cv, const gsize k, const gsize size)
{
  NcmVector *scv = g_object_new (NCM_TYPE_VECTOR, NULL);

  g_assert_cmpuint (size, >, 0);
  g_assert_cmpuint (size + k, <=, ncm_vector_len (cv));

  scv->vv    = gsl_vector_subvector (ncm_vector_gsl (cv), k, size);
  scv->type  = NCM_VECTOR_DERIVED;
  scv->pdata = ncm_vector_ref (cv);
  scv->pfree = (GDestroyNotify) & ncm_vector_free;

  return scv;
}

/**
 * ncm_vector_get_subvector2:
 * @sub_cv: a #NcmVector
 * @cv: a #NcmVector
 * @k: component index of the original vector
 * @size: number of components of the subvector
 *
 * This function sets @sub_cv to be a subvector of the vector @cv.
 * The start of the new vector is the component @k from the original vector @cv.
 * The new vector has @size elements.
 *
 * It is assumed that @sub_cv is a static vector allocated with
 * ncm_vector_new_data_static(). If a different type of vector is passed
 * then the function will lead to a memory leak.
 *
 */
void
ncm_vector_get_subvector2 (NcmVector *sub_cv, NcmVector *cv, const gsize k, const gsize size)
{
  g_assert_cmpuint (size, >, 0);
  g_assert_cmpuint (size + k, <=, ncm_vector_len (cv));

  sub_cv->vv = gsl_vector_subvector (ncm_vector_gsl (cv), k, size);
}

/**
 * ncm_vector_get_subvector_stride:
 * @cv: a #NcmVector
 * @k: index of the first component
 * @size: number of components
 * @stride: distance between consecutive components, in doubles
 *
 * Creates a view of @size components of @cv starting at @k, @stride doubles apart in
 * memory, holding a reference to @cv.
 *
 * Returns: (transfer full): the subvector.
 */
NcmVector *
ncm_vector_get_subvector_stride (NcmVector *cv, const gsize k, const gsize size, const gsize stride)
{
  g_assert_cmpuint (size,   >, 0);
  g_assert_cmpuint (stride, >, 0);
  {
    NcmVector *scv  = ncm_vector_new_data_static (ncm_vector_ptr (cv, k), size, stride);
    const gsize len = ncm_vector_len (cv);

    g_assert_cmpuint ((size - 1) * stride + k, <, len);

    scv->type  = NCM_VECTOR_DERIVED;
    scv->pdata = ncm_vector_ref (cv);
    scv->pfree = (GDestroyNotify) & ncm_vector_free;

    return scv;
  }
}

/**
 * ncm_vector_get_variant:
 * @cv: a constant #NcmVector
 *
 * Returns: (transfer full): a #GVariant of type `ad` with a copy of the components.
 */
GVariant *
ncm_vector_get_variant (const NcmVector *cv)
{
  guint n = ncm_vector_len (cv);
  GVariantBuilder builder;
  GVariant *var;
  guint i;

  g_variant_builder_init (&builder, G_VARIANT_TYPE ("ad"));

  for (i = 0; i < n; i++)
    g_variant_builder_add (&builder, "d", ncm_vector_get (cv, i));

  var = g_variant_builder_end (&builder);
  g_variant_ref_sink (var);

  return var;
}

/**
 * ncm_vector_peek_variant:
 * @cv: a constant #NcmVector
 *
 * Creates a #GVariant of type `ad` sharing the data of @cv, which must not change while
 * the variant exists; with a stride other than one it copies, as ncm_vector_get_variant().
 *
 * Returns: (transfer full): the #GVariant.
 */
GVariant *
ncm_vector_peek_variant (const NcmVector *cv)
{
  if (ncm_vector_stride (cv) != 1)
  {
    return ncm_vector_get_variant (cv);
  }
  else
  {
    guint n            = ncm_vector_len (cv);
    gconstpointer data = ncm_vector_const_ptr (cv, 0);
    GVariant *vvar     = g_variant_new_from_data (G_VARIANT_TYPE ("ad"),
                                                  data,
                                                  sizeof (gdouble) * n,
                                                  TRUE,
                                                  (GDestroyNotify) & ncm_vector_const_free,
                                                  NCM_VECTOR ((NcmVector *) ncm_vector_const_ref (cv)));

    return g_variant_ref_sink (vvar);
  }
}

/**
 * ncm_vector_log_vals:
 * @cv: a constant #NcmVector
 * @prestr: prefix
 * @format: printf format of one component
 * @cr: whether to end with a newline
 *
 * Logs the components of @cv.
 */
void
ncm_vector_log_vals (const NcmVector *cv, const gchar *prestr, const gchar *format, gboolean cr)
{
  guint i         = 0;
  const guint len = ncm_vector_len (cv);

  g_message ("%s", prestr);

  g_message (format, ncm_vector_get (cv, i));

  for (i = 1; i < len; i++)
  {
    g_message (" ");
    g_message (format, ncm_vector_get (cv, i));
  }

  if (cr)
    g_message ("\n");
}

/**
 * ncm_vector_log_vals_avpb:
 * @cv: a constant #NcmVector
 * @prestr: prefix
 * @format: printf format of one component
 * @a: factor $a$
 * @b: offset $b$
 *
 * Logs $a v_i + b$ for each component $v_i$ of @cv, ending with a newline.
 */
void
ncm_vector_log_vals_avpb (const NcmVector *cv, const gchar *prestr, const gchar *format, const gdouble a, const gdouble b)
{
  guint i         = 0;
  const guint len = ncm_vector_len (cv);

  g_message ("%s", prestr);

  g_message (format, a * ncm_vector_get (cv, i) + b);

  for (i = 1; i < len; i++)
  {
    g_message (" ");
    g_message (format, a * ncm_vector_get (cv, i) + b);
  }

  g_message ("\n");
}

/**
 * ncm_vector_log_vals_func:
 * @cv: a constant #NcmVector
 * @prestr: prefix
 * @format: printf format of one component
 * @f: (scope notified): a #NcmVectorCompFunc
 * @user_data: user data
 *
 * Logs @f of each component of @cv, called with @user_data, ending with a newline.
 */
void
ncm_vector_log_vals_func (const NcmVector *cv, const gchar *prestr, const gchar *format, NcmVectorCompFunc f, gpointer user_data)
{
  guint i         = 0;
  const guint len = ncm_vector_len (cv);

  g_message ("%s", prestr);

  g_message (format, f (ncm_vector_get (cv, i), i, user_data));

  for (i = 1; i < len; i++)
  {
    g_message (" ");
    g_message (format, f (ncm_vector_get (cv, i), i, user_data));
  }

  g_message ("\n");
}

/**
 * ncm_vector_const_new_gsl: (skip)
 * @gv: a constant #gsl_vector
 *
 * Creates a constant vector over @gv, which must outlive it.
 *
 * Returns: (transfer full): a new constant #NcmVector.
 */

/**
 * ncm_vector_get:
 * @cv: a constant #NcmVector
 * @i: component index
 *
 * Returns: the component @i.
 */

/**
 * ncm_vector_fast_get:
 * @cv: a constant #NcmVector
 * @i: component index
 *
 * Same as ncm_vector_get(), for a vector of stride one.
 *
 * Returns: the component @i.
 */

/**
 * ncm_vector_ptr:
 * @cv: a #NcmVector
 * @i: component index
 *
 * Returns: a pointer to the component @i.
 */

/**
 * ncm_vector_fast_ptr:
 * @cv: a #NcmVector
 * @i: component index
 *
 * Same as ncm_vector_ptr(), for a vector of stride one.
 *
 * Returns: a pointer to the component @i.
 */

/**
 * ncm_vector_const_ptr:
 * @cv: a constant #NcmVector
 * @i: component index
 *
 * Returns: a constant pointer to the component @i.
 */

/**
 * ncm_vector_set:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a double
 *
 * Sets the component @i to @val.
 */

/**
 * ncm_vector_fast_set:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a double
 *
 * Same as ncm_vector_set(), for a vector of stride one.
 */

/**
 * ncm_vector_addto:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a double
 *
 * Adds @val to the component @i.
 */

/**
 * ncm_vector_fast_addto:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a double
 *
 * Same as ncm_vector_addto(), for a vector of stride one.
 */

/**
 * ncm_vector_subfrom:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a double
 *
 * Subtracts @val from the component @i.
 */

/**
 * ncm_vector_fast_subfrom:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a double
 *
 * Same as ncm_vector_subfrom(), for a vector of stride one.
 */

/**
 * ncm_vector_mulby:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a double
 *
 * Multiplies the component @i by @val.
 */

/**
 * ncm_vector_fast_mulby:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a double
 *
 * Same as ncm_vector_mulby(), for a vector of stride one.
 */

/**
 * ncm_vector_set_all:
 * @cv: a #NcmVector
 * @val: a double
 *
 * Sets every component to @val.
 */

/**
 * ncm_vector_set_data:
 * @cv: a #NcmVector
 * @array: (array length=size) (element-type double): the values
 * @size: length of @array
 *
 * Sets the components of @cv to @array, whose length must be that of @cv.
 */

/**
 * ncm_vector_set_array:
 * @cv: a #NcmVector
 * @array: (array) (element-type double): a #GArray of doubles
 *
 * Sets the components of @cv to @array, whose length must be that of @cv.
 */

/**
 * ncm_vector_scale:
 * @cv: a #NcmVector
 * @val: a double
 *
 * Multiplies every component by @val.
 */

/**
 * ncm_vector_add_constant:
 * @cv: a #NcmVector
 * @val: a double
 *
 * Adds @val to every component.
 */

/**
 * ncm_vector_mul:
 * @cv1: a #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Multiplies @cv1 by @cv2, component by component.
 */

/**
 * ncm_vector_div:
 * @cv1: a #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Divides @cv1 by @cv2, component by component.
 */

/**
 * ncm_vector_add:
 * @cv1: a #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Adds @cv2 to @cv1.
 */

/**
 * ncm_vector_sub:
 * @cv1: a #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Subtracts @cv2 from @cv1.
 */

/**
 * ncm_vector_sqr_dist:
 * @cv1: a constant #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Returns: $\sum_i (v_{1,i} - v_{2,i})^2$.
 */

/**
 * ncm_vector_set_zero:
 * @cv: a #NcmVector
 *
 * Sets every component to zero.
 */

/**
 * ncm_vector_memcpy:
 * @cv1: a #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Copies @cv2 into @cv1.
 */

/**
 * ncm_vector_memcpy2:
 * @cv1: a #NcmVector
 * @cv2: a constant #NcmVector
 * @cv1_start: component of @cv1
 * @cv2_start: component of @cv2
 * @size: number of components
 *
 * This function copies @size components of the vector @cv2, counting from @cv2_start,
 * to the vector @cv1, starting from the @cv1_start component.
 * It is useful for vectors with different sizes.
 *
 */

/**
 * ncm_vector_get_array:
 * @cv: a #NcmVector
 *
 * @cv must have been created by ncm_vector_new_array().
 *
 * Returns: (transfer full) (element-type double): a new reference to the #GArray of @cv.
 */

/**
 * ncm_vector_dup_array:
 * @cv: a #NcmVector
 *
 * Returns: (transfer full) (element-type double): a new #GArray with a copy of the components.
 */

/**
 * ncm_vector_data:
 * @cv: a #NcmVector
 *
 * Returns: (transfer none): a pointer to the first component.
 */

/**
 * ncm_vector_const_data:
 * @cv: a constant #NcmVector
 *
 * Returns: (transfer none): a constant pointer to the first component.
 */

/**
 * ncm_vector_replace_data: (skip)
 * @cv: a #NcmVector
 * @data: the new data
 *
 * Points @cv to @data, keeping its length and stride. The previous data must not be owned
 * by @cv; nothing is checked.
 */
/**
 * ncm_vector_replace_data_full: (skip)
 * @cv: a #NcmVector
 * @data: the new data
 * @size: number of components
 * @stride: distance between consecutive components, in doubles
 *
 * Points @cv to @data with the given length and stride. The previous data must not be
 * owned by @cv; nothing is checked.
 */
/**
 * ncm_vector_gsl: (skip)
 * @cv: a #NcmVector
 *
 * Returns: the #gsl_vector of @cv.
 */

/**
 * ncm_vector_const_gsl: (skip)
 * @cv: a constant #NcmVector
 *
 * Returns: the constant #gsl_vector of @cv.
 */

/**
 * ncm_vector_dot:
 * @cv1: a constant #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Returns: $\sum_i v_{1,i} v_{2,i}$.
 */
gdouble
ncm_vector_dot (const NcmVector *cv1, const NcmVector *cv2)
{
  return cblas_ddot (ncm_vector_len (cv1), ncm_vector_const_data (cv1), ncm_vector_stride (cv1), ncm_vector_const_data (cv2), ncm_vector_stride (cv2));
}

/**
 * ncm_vector_len:
 * @cv: a constant #NcmVector
 *
 * Returns: the number of components.
 */

/**
 * ncm_vector_stride:
 * @cv: a constant #NcmVector
 *
 * Returns: the distance between consecutive components, in doubles.
 */

/**
 * ncm_vector_get_max:
 * @cv: a constant #NcmVector
 *
 * Returns: the largest component.
 */

/**
 * ncm_vector_get_min:
 * @cv: a constant #NcmVector
 *
 * Returns: the smallest component.
 */

/**
 * ncm_vector_get_max_index:
 * @cv: a constant #NcmVector
 *
 * Returns: the index of the largest component.
 */

/**
 * ncm_vector_get_min_index:
 * @cv: a constant #NcmVector
 *
 * Returns: the index of the smallest component.
 */

/**
 * ncm_vector_get_minmax:
 * @cv: a constant #NcmVector
 * @min: (out): the smallest component
 * @max: (out): the largest component
 *
 * Finds the smallest and largest components.
 */

/**
 * ncm_vector_is_finite:
 * @cv: a constant #NcmVector
 *
 * Returns: whether every component is finite.
 */

/**
 * ncm_vector_lt:
 * @cv1: a constant #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Returns: whether $v_{1,i} < v_{2,i}$ for every $i$.
 */
/**
 * ncm_vector_lteq:
 * @cv1: a constant #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Returns: whether $v_{1,i} \leq v_{2,i}$ for every $i$.
 */
/**
 * ncm_vector_between:
 * @cv: a constant #NcmVector
 * @cv_lb: a constant #NcmVector of lower bounds
 * @cv_ub: a constant #NcmVector of upper bounds
 * @type: 0 or 1
 *
 * Checks every component against its bounds: $l_i \leq v_i < u_i$ for @type 0, and
 * $l_i < v_i \leq u_i$ for @type 1. Aborts for other values of @type.
 *
 * Returns: whether every component is within its bounds.
 */

/**
 * ncm_vector_get_absminmax:
 * @cv: a constant #NcmVector
 * @absmin: (out): the smallest absolute value
 * @absmax: (out): the largest absolute value
 *
 * Finds the smallest and largest absolute values of the components.
 */
void
ncm_vector_get_absminmax (const NcmVector *cv, gdouble *absmin, gdouble *absmax)
{
  guint size = ncm_vector_len (cv);
  guint i;

  *absmin = HUGE_VAL;
  *absmax = 0.0;

  for (i = 0; i < size; i++)
  {
    const gdouble v = fabs (ncm_vector_get (cv, i));

    *absmin = GSL_MIN (*absmin, v);
    *absmax = GSL_MAX (*absmax, v);
  }
}

/**
 * ncm_vector_find_closest_index:
 * @cv: a constant #NcmVector
 * @x: a double $x$
 *
 * Finds by bisection the largest $i < n - 1$ with $v_i \leq x$, for components in
 * increasing order and $v_0 \leq x$, where $n$ is the length of @cv.
 *
 * Returns: the index $i$.
 */
guint
ncm_vector_find_closest_index (const NcmVector *cv, const gdouble x)
{
  gsize ilo = 0;
  gsize ihi = ncm_vector_len (cv) - 1;

  while (ihi > ilo + 1)
  {
    gsize i = (ihi + ilo) / 2;

    if (ncm_vector_get (cv, i) > x)
      ihi = i;
    else
      ilo = i;
  }

  return ilo;
}

/**
 * ncm_vector_set_from_variant:
 * @cv: a #NcmVector
 * @var: a #GVariant of type `ad`
 *
 * Sets the components of @cv to the elements of @var. A vector of length zero is allocated
 * with the length of @var. Aborts if @var has another type or length.
 */
void
ncm_vector_set_from_variant (NcmVector *cv, GVariant *var)
{
  gsize n;
  guint i;

  if (!g_variant_is_of_type (var, G_VARIANT_TYPE ("ad")))
    g_error ("ncm_vector_set_from_variant: Cannot convert `%s' variant to an array of doubles", g_variant_get_type_string (var));

  n = g_variant_n_children (var);

  if (ncm_vector_len (cv) == 0)
  {
    gdouble *d = g_slice_alloc (sizeof (gdouble) * n);

    cv->vv    = gsl_vector_view_array (d, n);
    cv->pdata = NULL;
    cv->pfree = NULL;
    cv->type  = NCM_VECTOR_SLICE;
  }
  else if (n != ncm_vector_len (cv))
  {
    g_error ("set_property: cannot set vector values, variant contains %zu childs but vector dimension is %u", n, ncm_vector_len (cv));
  }

  for (i = 0; i < n; i++)
  {
    gdouble val = 0.0;

    g_variant_get_child (var, i, "d", &val);
    ncm_vector_set (cv, i, val);
  }
}

/**
 * ncm_vector_dnrm2:
 * @cv: a constant #NcmVector
 *
 * Calculates the Euclidean norm of the vector @cv, i.e.,
 * $\vert\text{cv}\vert_2$.
 *
 * Returns: $\vert\text{cv}\vert_2$.
 */
gdouble
ncm_vector_dnrm2 (const NcmVector *cv)
{
  return cblas_dnrm2 (ncm_vector_len (cv),
                      ncm_vector_const_ptr (cv, 0),
                      ncm_vector_stride (cv));
}

/**
 * ncm_vector_axpy:
 * @cv1: a #NcmVector
 * @alpha: $\alpha$
 * @cv2: a constant #NcmVector
 *
 * Sets $v_1 \to v_1 + \alpha v_2$.
 */
void
ncm_vector_axpy (NcmVector *cv1, const gdouble alpha, const NcmVector *cv2)
{
  const guint len = ncm_vector_len (cv1);

  g_assert_cmpuint (len, ==, ncm_vector_len (cv2));

  cblas_daxpy (len, alpha, ncm_vector_const_data (cv2), ncm_vector_stride (cv2), ncm_vector_data (cv1), ncm_vector_stride (cv1));
}

/**
 * ncm_vector_cmp:
 * @cv1: a #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Sets each component of @cv1 to $|v_{1,i} - v_{2,i}| / \min(|v_{1,i}|, |v_{2,i}|)$, or
 * to the absolute value of the other component when one of them is zero.
 */
void
ncm_vector_cmp (NcmVector *cv1, const NcmVector *cv2)
{
  const guint len = ncm_vector_len (cv1);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble x1_i = ncm_vector_get (cv1, i);
    const gdouble x2_i = ncm_vector_get (cv2, i);

    if (G_UNLIKELY (x1_i == 0.0))
    {
      if (G_UNLIKELY (x2_i == 0.0))
        ncm_vector_set (cv1, i, 0.0);
      else
        ncm_vector_set (cv1, i, fabs (x2_i));
    }
    else if (G_UNLIKELY (x2_i == 0.0))
    {
      ncm_vector_set (cv1, i, fabs (x1_i));
    }
    else
    {
      const gdouble abs_x1_i  = fabs (x1_i);
      const gdouble abs_x2_i  = fabs (x2_i);
      const gdouble max_x12_i = GSL_MIN (abs_x1_i, abs_x2_i);

      ncm_vector_set (cv1, i, fabs ((x1_i - x2_i) / max_x12_i));
    }
  }
}

/**
 * ncm_vector_cmp2:
 * @cv1: a constant #NcmVector
 * @cv2: a constant #NcmVector
 * @reltol: the relative tolerance
 * @abstol: the absolute tolerance
 *
 * Performs a comparison, component-wise, of the two vectors and
 * returns the number of components that do not match.
 *
 * Returns: the number of components that do not match.
 */
gint
ncm_vector_cmp2 (const NcmVector *cv1, const NcmVector *cv2, const gdouble reltol, const gdouble abstol)
{
  const guint len1 = ncm_vector_len (cv1);
  const guint len2 = ncm_vector_len (cv2);
  gint n           = 0;
  guint i;

  g_assert_cmpint (len1, ==, len2);

  for (i = 0; i < len1; i++)
  {
    const gdouble x1_i = ncm_vector_get (cv1, i);
    const gdouble x2_i = ncm_vector_get (cv2, i);

    n += abs (ncm_cmp (x1_i, x2_i, reltol, abstol));
  }

  return n;
}

/**
 * ncm_vector_sub_round_off:
 * @cv1: a #NcmVector
 * @cv2: a constant #NcmVector
 *
 * Sets each component of @cv1 to $\epsilon \max(|v_{1,i}|, |v_{2,i}|) / |v_{2,i} - v_{1,i}|$,
 * the relative round-off error of the difference, with $\epsilon$ the double precision;
 * or to one when the components are equal.
 */
void
ncm_vector_sub_round_off (NcmVector *cv1, const NcmVector *cv2)
{
  const guint len = ncm_vector_len (cv1);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble x1_i        = ncm_vector_get (cv1, i);
    const gdouble x2_i        = ncm_vector_get (cv2, i);
    const gdouble abs_x1_i    = fabs (x1_i);
    const gdouble abs_x2_i    = fabs (x2_i);
    const gdouble x2_i_m_x1_i = fabs (x2_i - x1_i);

    if (G_UNLIKELY (x2_i_m_x1_i == 0.0))
    {
      ncm_vector_set (cv1, i, 1.0);
    }
    else
    {
      const gdouble max_x12_i = GSL_MAX (abs_x1_i, abs_x2_i);

      ncm_vector_set (cv1, i, max_x12_i * GSL_DBL_EPSILON / x2_i_m_x1_i);
    }
  }
}

/**
 * ncm_vector_reciprocal:
 * @cv: a #NcmVector
 *
 * Sets each component to its reciprocal.
 */
void
ncm_vector_reciprocal (NcmVector *cv)
{
  const guint len = ncm_vector_len (cv);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble x_i = ncm_vector_get (cv, i);

    ncm_vector_set (cv, i, 1.0 / x_i);
  }
}

/**
 * ncm_vector_square:
 * @cv: a #NcmVector
 *
 * Sets each component to its square.
 */
void
ncm_vector_square (NcmVector *cv)
{
  const guint len = ncm_vector_len (cv);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble x_i = ncm_vector_get (cv, i);

    ncm_vector_set (cv, i, x_i * x_i);
  }
}

/**
 * ncm_vector_sqrt:
 * @cv: a #NcmVector
 *
 * Sets each component to its square root.
 */
void
ncm_vector_sqrt (NcmVector *cv)
{
  const guint len = ncm_vector_len (cv);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble x_i = ncm_vector_get (cv, i);

    ncm_vector_set (cv, i, sqrt (x_i));
  }
}

/**
 * ncm_vector_hypot:
 * @cv1: a #NcmVector
 * @alpha: $\alpha$
 * @cv2: a constant #NcmVector
 *
 * Sets $v_{1,i} \to \sqrt{v_{1,i}^2 + (\alpha v_{2,i})^2}$.
 */
void
ncm_vector_hypot (NcmVector *cv1, const gdouble alpha, const NcmVector *cv2)
{
  const guint len = ncm_vector_len (cv1);
  guint i;

  g_assert_cmpuint (len, ==, ncm_vector_len (cv2));

  for (i = 0; i < len; i++)
  {
    const gdouble x1_i = ncm_vector_get (cv1, i);
    const gdouble x2_i = ncm_vector_get (cv2, i);

    ncm_vector_set (cv1, i, hypot (x1_i, alpha * x2_i));
  }
}

/**
 * ncm_vector_sum_cpts:
 * @cv: a constant #NcmVector
 *
 * Returns: the sum of the components.
 */


/**
 * ncm_vector_mean:
 * @cv: a constant #NcmVector
 *
 * Returns: the mean of the components.
 */

