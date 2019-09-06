/***************************************************************************
 *            ncm_vector.c
 *
 *  Tue Jul  8 15:05:41 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_vector
 * @title: NcmVector
 * @short_description: Vector object representing arrays of doubles.
 * @stability: Stable
 * @include: numcosmo/math/ncm_vector.h
 *
 * This object defines the functions for allocating and accessing vectors.
 * Also includes several vector operations.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_vector.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_VALS,
};

G_DEFINE_TYPE (NcmVector, ncm_vector, G_TYPE_OBJECT);

static void
ncm_vector_init (NcmVector *v)
{
  v->pdata = NULL;
  v->pfree = NULL;
  v->type = 0;
  memset (&v->vv, 0, sizeof (gsl_vector_view));
}

static void
_ncm_vector_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmVector *v = NCM_VECTOR (object);
  g_return_if_fail (NCM_IS_VECTOR (object));

  switch (prop_id)
  {
    case PROP_VALS:
    {
      GVariant *var = ncm_vector_get_variant (v);
      g_value_take_variant (value, var);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_vector_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmVector *v = NCM_VECTOR (object);
  g_return_if_fail (NCM_IS_VECTOR (object));

  switch (prop_id)
  {
    case PROP_VALS:
    {
      GVariant *var = g_value_get_variant (value);
      ncm_vector_set_from_variant (v, var);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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

  g_object_class_install_property (object_class, PROP_VALS,
                                   g_param_spec_variant ("values", NULL, "values",
                                                         G_VARIANT_TYPE_ARRAY, NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}

/**
 * ncm_vector_new:
 * @n: defines the size of the vector
 *
 * This function allocates memory for a new #NcmVector of double
 * with @n components.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new (gsize n)
{
  gdouble *d = g_slice_alloc (sizeof (gdouble) * n);  
  return ncm_vector_new_data_slice (d, n, 1);
}

/**
 * ncm_vector_new_full:
 * @d: (array) (element-type double): pointer to the first double allocated
 * @size: number of doubles allocated
 * @stride: the step-size from one element to the next in physical memory, measured in units of double
 * @pdata: (allow-none): descending data pointer
 * @pfree: (scope notified) (allow-none): free function to be called when destroying the vector
 *
 * This function returns a #NcmVector of the array @d.
 * This function saves @userdata internally and frees it using @free
 * when it is no longer necessary.
 *
 * Returns: A new #NcmVector.
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
 * @size: number of doubles allocated
 *
 * This function allocates memory for a new #NcmVector of double
 * with @n components. It uses fftw_alloc_real in order to be used
 * by fftw* functions.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_fftw (guint size)
{
  gdouble *d    = fftw_alloc_real (size);
  NcmVector *cv = ncm_vector_new_full (d, size, 1, d, (GDestroyNotify) fftw_free);
  cv->type      = NCM_VECTOR_MALLOC;

  return cv;
}

/**
 * ncm_vector_new_gsl: (skip)
 * @gv: vector from GNU Scientific Library (GSL) to be converted into a #NcmVector
 *
 * This function saves @gv internally and frees it when it is no longer necessary.
 * The @gv vector must not be freed.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_gsl (gsl_vector *gv)
{
  NcmVector *cv = ncm_vector_new_full (gv->data, gv->size, gv->stride, gv, (GDestroyNotify)gsl_vector_free);
  cv->type = NCM_VECTOR_GSL_VECTOR;
  return cv;
}

/**
 * ncm_vector_new_gsl_static: (skip)
 * @gv: vector from GNU Scientific Library (GSL) to be converted into a #NcmVector
 *
 * This function saves @gv internally and does not frees.
 * The @gv vector must be valid during the life of the created #NcmVector.
 *
 * Returns: A new #NcmVector.
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
 * @a: (array) (element-type double): array of doubles to be converted into a #NcmVector
 *
 * This function saves @a internally and frees it when it is no longer necessary.
 * The @a array must not be freed.
 *
 * Returns: (transfer full): A new #NcmVector.
 */
NcmVector *
ncm_vector_new_array (GArray *a)
{
	g_assert_cmpint (a->len, >, 0);
	{
		NcmVector *cv = ncm_vector_new_full (&g_array_index (a, gdouble, 0), a->len, 1,
		                                     g_array_ref (a), (GDestroyNotify) &g_array_unref);
		cv->type = NCM_VECTOR_ARRAY;
		return cv;
	}
}

/**
 * ncm_vector_new_data_slice:
 * @d: (array) (element-type double): pointer to the first double allocated
 * @size: number of doubles allocated
 * @stride: the step-size from one element to the next in physical memory, measured in units of double
 *
 * This function returns a #NcmVector of the array @d allocated using g_slice function.
 * This function saves @a internally and frees it when it is no longer necessary.
 * The @a vector must not be freed.
 *
 * Returns: A new #NcmVector.
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
 * @d: (array) (element-type double): pointer to the first double allocated
 * @size: number of doubles allocated
 * @stride: the step-size from one element to the next in physical memory, measured in units of double
 *
 * This function returns a #NcmVector of the array @d allocated using malloc.
 * It saves @d internally and frees it when it is no longer necessary.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_malloc (gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = ncm_vector_new_full (d, size, stride,
                                       d, &g_free);
  cv->type  = NCM_VECTOR_MALLOC;
  return cv;
}

/**
 * ncm_vector_new_data_static:
 * @d: (array) (element-type double): pointer to the first double allocated
 * @size: number of doubles allocated
 * @stride: the step-size from one element to the next in physical memory, measured in units of double
 *
 * This function returns a #NcmVector of the array @d.
 * The memory allocated is kept during all time life of the object and
 * must not be freed during this period.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_static (gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = ncm_vector_new_full (d, size, stride,
                                       NULL, NULL);
  cv->type  = NCM_VECTOR_DERIVED;
  return cv;
}

/**
 * ncm_vector_new_data_dup:
 * @d: (array) (element-type double): pointer to the first double allocated
 * @size: number of doubles allocated
 * @stride: the step-size from one element to the next in physical memory, measured in units of double
 *
 * This function returns a #NcmVector of the array @d.
 * It allocate a new vector and copy the contents of @d into it.
 *
 * Returns: (transfer full): A new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_dup (gdouble *d, const gsize size, const gsize stride)
{
  NcmVector *s = ncm_vector_new_data_static (d, size, stride);
  NcmVector *dup = ncm_vector_dup (s);
  ncm_vector_free (s);
  return dup;
}

/**
 * ncm_vector_new_variant:
 * @var: a #GVariant of the type "ad"
 *
 * This function convert a #GVariant array to a #NcmVector allocating new
 * memory for the vector.
 *
 * Returns: A new #NcmVector.
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
 * @d: (array) (element-type double): pointer to the first double allocated
 * @size: number of doubles allocated
 * @stride: the step-size from one element to the next in physical memory, measured in units of double
 *
 * This function returns a constant #NcmVector of the array @d.
 * The memory allocated is kept during all time life of the object and
 * must not be freed during this period.
 *
 * Returns: A new constant #NcmVector.
 */
const NcmVector *
ncm_vector_const_new_data (const gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);
  if (stride != 1)
    cv->vv = gsl_vector_view_array_with_stride ((gdouble *)d, stride, size);
  else
    cv->vv = gsl_vector_view_array ((gdouble *)d, size);
  cv->pdata = NULL;
  cv->pfree = NULL;

  cv->type = NCM_VECTOR_DERIVED;
  return cv;
}

/**
 * ncm_vector_ref:
 * @cv: a #NcmVector
 *
 * This function increses the reference count of the vector @cv.
 *
 * Returns: (transfer full): @cv
 */
NcmVector *
ncm_vector_ref (NcmVector *cv)
{
  return g_object_ref (cv);
}

/**
 * ncm_vector_const_ref:
 * @cv: a #NcmVector
 *
 * This function increses the reference count of the constant vector @cv.
 *
 * Returns: (transfer full): @cv
 */
const NcmVector *
ncm_vector_const_ref (const NcmVector *cv)
{
  return g_object_ref (NCM_VECTOR (cv));
}

/**
 * ncm_vector_const_new_variant:
 * @var: a #GVariant of the type "ad"
 *
 * This function convert a #GVariant array to a #NcmVector. Since it returns
 * a constant #NcmVector it uses the same memory of @var.
 *
 * Returns: (transfer full): A new #NcmVector.
 */
const NcmVector *
ncm_vector_const_new_variant (GVariant *var)
{
  gsize n = g_variant_n_children (var);
  gconstpointer data = g_variant_get_data (var);
  const NcmVector *v = ncm_vector_const_new_data (data, n, 1);

  NCM_VECTOR (v)->pdata = g_variant_ref_sink (var);
  NCM_VECTOR (v)->pfree = (GDestroyNotify) &g_variant_unref;

  return v;
}


/**
 * ncm_vector_free:
 * @cv: a #NcmVector
 *
 * Atomically decrements the reference count of @cv by one. If the reference count drops to 0,
 * all memory allocated by @cv is released.
 *
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
 * Atomically decrements the reference count of @cv by one. If the reference count drops to 0,
 * all memory allocated by @cv is released.
 *
 */
void
ncm_vector_const_free (const NcmVector *cv)
{
  ncm_vector_free (NCM_VECTOR (cv));
}

/**
 * ncm_vector_clear:
 * @cv: a #NcmVector
 *
 * Atomically decrements the reference count of @cv by one. If the reference count drops to 0,
 * all memory allocated by @cv is released. The pointer is set to NULL.
 *
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
 * This function copies the elements of the constant vector @cv into a new #NcmVector.
 *
 * Returns: (transfer full): A #NcmVector.
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
 * @cv: a #NcmVector
 * @nv: a #NcmVector
 * @check_size: whether to check vector size
 *
 * This function substitute the vector *@cv by @nv, it will unref *@cv first.
 * If @check_size is TRUE then the function asserts that both vectors have the
 * same size.
 *
 */
void
ncm_vector_substitute (NcmVector **cv, NcmVector *nv, gboolean check_size)
{
  if (*cv == nv)
    return;

  if (*cv != NULL)
  {
    if (nv != NULL && check_size)
      g_assert_cmpuint (ncm_vector_len (*cv), ==, ncm_vector_len (nv));
    ncm_vector_clear (cv);
  }

  if (nv != NULL)
    *cv = ncm_vector_ref (nv);
}

/**
 * ncm_vector_get_subvector:
 * @cv: a #NcmVector
 * @k: component index of the original vector
 * @size: number of components of the subvector
 *
 * This function returns a #NcmVector which is a subvector of the vector @cv.
 * The start of the new vector is the component @k from the original vector @cv.
 * The new vector has @size elements.
 *
 * Returns: (transfer full): A #NcmVector.
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
  scv->pfree = (GDestroyNotify) &ncm_vector_free;

  return scv;
}

/**
 * ncm_vector_get_subvector_stride:
 * @cv: a #NcmVector
 * @k: component index of the original vector
 * @size: number of components of the subvector
 * @stride: the step-size from one element to the next in physical memory, measured in units of double
 *
 * This function returns a #NcmVector which is a subvector of the vector @cv.
 * The start of the new vector is the component @k from the original vector @cv.
 * The new vector has @size elements.
 *
 * Returns: (transfer full): A #NcmVector.
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
		scv->pfree = (GDestroyNotify) &ncm_vector_free;

		return scv;
	}
}

/**
 * ncm_vector_get_variant:
 * @v: a #NcmVector
 *
 * Convert @v to a GVariant of the type "ad" without destroying the
 * original vector @v;
 *
 * Returns: (transfer full): A #GVariant of the type "ad".
 */
GVariant *
ncm_vector_get_variant (const NcmVector *v)
{
  guint n = ncm_vector_len (v);
  GVariantBuilder builder;
  GVariant *var;
  guint i;

  g_variant_builder_init (&builder, G_VARIANT_TYPE ("ad"));

  for (i = 0; i < n; i++)
    g_variant_builder_add (&builder, "d", ncm_vector_get (v, i));

  var = g_variant_builder_end (&builder);
  g_variant_ref_sink (var);

  return var;
}

/**
 * ncm_vector_peek_variant:
 * @v: a #NcmVector
 *
 * Convert @v to a GVariant of the type "ad" using the same memory space.
 * The vector @v should not be modified during the variant existance.
 * If the vector has stride != 1 then ncm_vector_get_variant() is called.
 *
 * Returns: (transfer full): A #GVariant of the type "ad".
 */
GVariant *
ncm_vector_peek_variant (const NcmVector *v)
{
  if (ncm_vector_stride (v) != 1)
    return ncm_vector_get_variant (v);
  else
  {
    guint n = ncm_vector_len (v);
    gconstpointer data = ncm_vector_const_ptr (v, 0);
    GVariant *vvar = g_variant_new_from_data (G_VARIANT_TYPE ("ad"),
                                              data,
                                              sizeof (gdouble) * n,
                                              TRUE,
                                              (GDestroyNotify) &ncm_vector_const_free,
                                              NCM_VECTOR (ncm_vector_const_ref (v)));
    return g_variant_ref_sink (vvar);
  }
}

/**
 * ncm_vector_log_vals:
 * @v: a #NcmVector
 * @prestr: initial string
 * @format: float format
 * @cr: whether to include a cariage return
 * 
 * Log the vector values using @prestr and @format.
 *
 */
void 
ncm_vector_log_vals (const NcmVector *v, const gchar *prestr, const gchar *format, gboolean cr)
{
  guint i = 0;
  const guint len = ncm_vector_len (v);
  g_message ("%s", prestr);

  g_message (format, ncm_vector_get (v, i));
  for (i = 1; i < len; i++)
  {
    g_message (" ");
    g_message (format, ncm_vector_get (v, i));
  }
  if (cr)
    g_message ("\n");
}

/**
 * ncm_vector_log_vals_avpb:
 * @v: a #NcmVector
 * @prestr: initial string
 * @format: float format
 * @a: a double
 * @b: a double
 *
 * Log the vector values ($a\vec{v}+b$) using @prestr and @format.
 *
 */
void
ncm_vector_log_vals_avpb (const NcmVector *v, const gchar *prestr, const gchar *format, const gdouble a, const gdouble b)
{
  guint i = 0;
  const guint len = ncm_vector_len (v);
  g_message ("%s", prestr);

  g_message (format, a * ncm_vector_get (v, i) + b);
  for (i = 1; i < len; i++)
  {
    g_message (" ");
    g_message (format, a * ncm_vector_get (v, i) + b);
  }
  g_message ("\n");
}

/**
 * ncm_vector_log_vals_func:
 * @v: a #NcmVector
 * @prestr: initial string
 * @format: float format
 * @f: (scope notified): a #NcmVectorCompFunc
 * @user_data: user data used in @f
 *
 * Log the vector values (f(\vec{v}_i)$) using @prestr and @format.
 *
 */
void
ncm_vector_log_vals_func (const NcmVector *v, const gchar *prestr, const gchar *format, NcmVectorCompFunc f, gpointer user_data)
{
  guint i = 0;
  const guint len = ncm_vector_len (v);
  g_message ("%s", prestr);

  g_message (format, f (ncm_vector_get (v, i), i, user_data));
  for (i = 1; i < len; i++)
  {
    g_message (" ");
    g_message (format, f (ncm_vector_get (v, i), i, user_data));
  }
  g_message ("\n");
}

/**
 * ncm_vector_const_new_gsl: (skip)
 * @v: vector from GNU Scientific Library (GSL)
 *
 * This function converts @v into a constant #NcmVector.
 *
 * Returns: A new constant #NcmVector
 */
/**
 * ncm_vector_get:
 * @cv: a constant #NcmVector
 * @i: component index
 *
 * Returns: The @i-th component of the vector @cv.
 */
/**
 * ncm_vector_fast_get:
 * @cv: a constant #NcmVector
 * @i: component index
 *
 * Returns: The @i-th component of the vector @cv assuming stride == 1.
 */
/**
 * ncm_vector_ptr:
 * @cv: a #NcmVector
 * @i: component index
 *
 * Returns: A pointer to the @i-th component of the vector @cv.
 */
/**
 * ncm_vector_fast_ptr:
 * @cv: a #NcmVector
 * @i: component index
 *
 * Returns: A pointer to the @i-th component of the vector @cv assuming stride == 1.
 */
/**
 * ncm_vector_set:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a constant double
 *
 * This function sets the value of the @i-th component of the vector @cv to @val.
 */
/**
 * ncm_vector_fast_set:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a constant double
 *
 * This function sets the value of the @i-th component of the vector @cv to @val assuming stride == 1.
 */
/**
 * ncm_vector_addto:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a constant double
 *
 * This function adds @val to the value of the @i-th component of @cv.
 *
 */
/**
 * ncm_vector_fast_addto:
 * @cv: a #NcmVector.
 * @i: component index.
 * @val: a constant double.
 *
 * This function adds @val to the value of the @i-th component of @cv assuming stride == 1.
 *
 */
/**
 * ncm_vector_subfrom:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a constant double
 *
 * This function subtracts @val from the value of the @i-th component of @cv.
 *
 */
/**
 * ncm_vector_fast_subfrom:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a constant double
 *
 * This function subtracts @val from the value of the @i-th component of @cv assuming stride == 1.
 *
 */
/**
 * ncm_vector_mulby:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a constant double
 *
 * This function multiplies the @i-th component by @val.
 *
 */
/**
 * ncm_vector_fast_mulby:
 * @cv: a #NcmVector
 * @i: component index
 * @val: a constant double
 *
 * This function multiplies the @i-th component by @val assuming stride == 1.
 * 
 */
/**
 * ncm_vector_set_all:
 * @cv: a #NcmVector
 * @val: a constant double
 *
 * This function sets all the components of the vector @cv to the value @val.
 *
 */
/**
 * ncm_vector_set_data:
 * @cv: a #NcmVector.
 * @array: (array length=size) (element-type double): a pointer to a double array
 * @size: data array size
 * 
 * This function sets all the components of the vector @cv using the data array @array,
 * @size must match the vector size.
 *
 */
/**
 * ncm_vector_set_array:
 * @cv: a #NcmVector.
 * @array: (array) (element-type double): a pointer to a double #GArray
 * 
 * This function sets all the components of the vector @cv using the data array @array,
 * @array->len must match the vector size.
 *
 */
/**
 * ncm_vector_scale:
 * @cv: a #NcmVector
 * @val: a constant double
 *
 * This function multiplies the components of the vector @cv by the constant factor @val.
 *
 */
/**
 * ncm_vector_add_constant:
 * @cv: a #NcmVector.
 * @val: a cosntant double.
 *
 * This function adds the constant @val to all components.
 *
 */
/**
 * ncm_vector_mul:
 * @cv1: a #NcmVector, numerator
 * @cv2: a #NcmVector, denominator
 *
 * This function multiplies the components of the vector @cv1 by the components of the vector @cv2.
 * The two vectors must have the same length.
 *
 */
/**
 * ncm_vector_div:
 * @cv1: a #NcmVector, numerator
 * @cv2: a #NcmVector, denominator
 *
 * This function divides the components of the vector @cv1 by the components of the vector @cv2.
 * The two vectors must have the same length.
 *
 */
/**
 * ncm_vector_add:
 * @cv1: a #NcmVector
 * @cv2: a #NcmVector
 *
 * This function adds the components of the vector @cv2 to the components of the vector @cv1.
 * The two vectors must have the same length.
 *
 */
/**
 * ncm_vector_sub:
 * @cv1: a #NcmVector
 * @cv2: a #NcmVector
 *
 * This function subtracts the components of the vector @cv2 to the components of the vector @cv1.
 * The two vectors must have the same length.
 *
 */
/**
 * ncm_vector_set_zero:
 * @cv: a #NcmVector
 *
 * This function sets all the components of the vector @cv to zero.
 *
 */
/**
 * ncm_vector_memcpy:
 * @cv1: a #NcmVector
 * @cv2: a #NcmVector
 *
 * This function copies the components of the vector @cv2 into the vector @cv1.
 * The two vectors must have the same length.
 *
 */
/**
 * ncm_vector_memcpy2:
 * @cv1: a #NcmVector
 * @cv2: a #NcmVector
 * @cv1_start: component of @cv1
 * @cv2_start: component of @cv2
 * @size: number of components
 *
 * This function copies @size components of @cv2, counting from @cv2_start,
 * to the vector @cv1, starting from the @cv1_start component.
 * It is useful for vectors with different sizes.
 *
 */
/**
 * ncm_vector_get_array:
 * @cv: a #NcmVector
 *
 * FIXME
 *
 * Returns: (transfer container) (element-type double): FIXME
 */
/**
 * ncm_vector_dup_array:
 * @cv: a #NcmVector
 *
 * FIXME
 *
 * Returns: (transfer full) (element-type double): FIXME
 */
/**
 * ncm_vector_data:
 * @cv: a #NcmVector
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
/**
 * ncm_vector_const_data:
 * @cv: a #NcmVector
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
/**
 * ncm_vector_ddot:
 * @a: a #NcmVector
 * @b: a #NcmVector
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_vector_gsl: (skip)
 * @cv: a #NcmVector
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_vector_const_gsl: (skip)
 * @cv: a #NcmVector
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_vector_len:
 * @cv: a #NcmVector
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_vector_stride:
 * @cv: a #NcmVector
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_vector_get_max: 
 * @cv: a @NcmVector.
 * 
 * Gets the maximum value of the vector components.
 * 
 * Returns: the maximum value of the vector components.
 */
/**
 * ncm_vector_get_min: 
 * @cv: a @NcmVector.
 * 
 * Gets the minimum value of the vector components.
 * 
 * Returns: the minimum value of the vector components.
 */
/**
 * ncm_vector_get_max_index: 
 * @cv: a @NcmVector.
 * 
 * Gets the index of the maximal vector component.
 * 
 * Returns: the index of the maximal component.
 */
/**
 * ncm_vector_get_min_index: 
 * @cv: a @NcmVector.
 * 
 * Gets the index of the minimal vector component.
 * 
 * Returns: the index of the minimal component.
 */
/**
 * ncm_vector_get_minmax:
 * @cv: a @NcmVector.
 * @min: (out): minimum component value of @cv
 * @max: (out): maximum component value of @cv
 *
 * Gets the minimum/maximum value of the vector components.
 *
 */
/**
 * ncm_vector_is_finite:
 * @cv: a @NcmVector.
 *
 * Tests all entries, if one or more are not finite return FALSE.
 * Otherwise returns TRUE;
 * 
 * Returns: whether all components of @cv are finite.
 */

/**
 * ncm_vector_get_absminmax:
 * @cv: a @NcmVector.
 * @absmin: (out): minimum component absolute value of @cv
 * @absmax: (out): maximum component absolute value of @cv
 *
 * Gets the minimum/maximum absolute value of the vector components.
 *
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
 * ncm_vector_set_from_variant:
 * @cv: a #NcmVector
 * @var: a #GVariant of type ad
 *
 * Sets the values of @cv using the variant @var. This function fails
 * if @cv and @var differ in size.
 *
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
    cv->vv = gsl_vector_view_array (d, n);
    cv->pdata = NULL;
    cv->pfree = NULL;
    cv->type = NCM_VECTOR_SLICE;
  }
  else if (n != ncm_vector_len (cv))
    g_error ("set_property: cannot set vector values, variant contains %zu childs but vector dimension is %u", n, ncm_vector_len (cv));

  for (i = 0; i < n; i++)
  {
    gdouble val = 0.0;
    g_variant_get_child (var, i, "d", &val);
    ncm_vector_set (cv, i, val);
  }
}

/**
 * ncm_vector_dnrm2:
 * @cv: a @NcmVector
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
 * @cv1: a #NcmVector $y$
 * @alpha: a double $\alpha$
 * @cv2: a #NcmVector $x$
 * 
 * Performs the operation $y = \alpha x + y$.
 *
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
 * @cv1: a #NcmVector $x_1$
 * @cv2: a #NcmVector $x_2$
 * 
 * Performs a comparison, component-wise, of the two vectors and 
 * puts the weighted difference in @cv1.
 *
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
 * ncm_vector_sub_round_off:
 * @cv1: a #NcmVector $x_1$
 * @cv2: a #NcmVector $x_2$
 * 
 * Estimate the round-off error component wise in the
 * subtraction of @cv1 by @cv2.
 *
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
      ncm_vector_set (cv1, i, max_x12_i * GSL_DBL_EPSILON / x2_i_m_x1_i );
    }
  }
}

/**
 * ncm_vector_reciprocal:
 * @cv: a #NcmVector
 * 
 * Calculates the reciprocal of @cv1.
 *
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
 * ncm_vector_sum_cpts:
 * @cv: a @NcmVector
 *
 * Calculates the sum of the components.
 *
 * Returns: the sum of the vector @cv components.
 */
